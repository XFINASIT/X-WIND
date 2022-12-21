C	=====================================================================
C	=====================================================================
      SUBROUTINE MASS_AND_STIFFNESS_BLOCK (MAXA,ALFA,BETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /TIME/ DDT,CTIM,NINC
      DIMENSION MAXA(1)!,DATA_OUTPUT_REAL(1)
      
      ! FREQUENCY DOMAIN INPUT DATA
      
      !CALL MASS_BLOCK (MAXA)
      !CALL STIFFNESS_BLOCK (MAXA)
      
      CALL TRANSFER_FUNCTION (MAXA,NEQ,ALFA,BETA) 
      CALL COMBINE_FREQ_DATA
      
      ! IFFT FREQ > TIME(DISPLACEMENT)
      CALL FAST_FOURIER_BACKWARD_DISPLACEMENT (NINC,NEQ)
      ! PRINT DISPLACEMENT RESULT
      CALL PRINT_OUTPUT_DISPLACEMENT (NINC,NEQ)
 
      RETURN
      END
      
C	=====================================================================
C	===================================================================== 
      SUBROUTINE TRANSFER_FUNCTION (MAXA,NEQ,ALFA,BETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 ASIGN 
      INCLUDE 'MKL_PARDISO.F77'
      CHARACTER*4 TYP1,TYP2  
      DIMENSION RR(1),TT(1) 
      ALLOCATABLE B(:),SW(:,:)
      ALLOCATABLE LPT(:)
      POINTER :: AA(:),CC(:),MAXC_STIFFNESS(:),MAXC_MASS(:),IPERM(:),X(:)
      COMPLEX, POINTER :: CC_COMPLEX(:),X_COMPLEX(:),V_COMPLEX(:),SW_COMPLEX(:,:)
      DIMENSION MAXA(1),IPARM(64),V(NEQ)
      DIMENSION RFV(NEQ),RFC(NEQ)
      DIMENSION BB_OUT(NEQ,NEQ)
      !DIMENSION DD(78)
      !DIMENSION CC_TEST(6),MAXA_TEST(3),MAXC_TEST(6),X_TEST(3),V_TEST(3) 
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /PARDIS/ LPTO(64) 
      COMMON A(7000000),IA(6000000)

      
      PI = 3.141592654
      
      ! ------- STIFFNESS ------
C     IPFAC DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT
      !TYP1 = "STFF"
      TYP2 = "TEMP"
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,0,'FIRS')  

      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0,'FIRS')
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0,'FIRS')
      
	CALL MCALFIL(KSREC,TYP2)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
      
      ALLOCATE(AA(NWK),MAXC_STIFFNESS(NWK)) 
	IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0,'FIRS')
      
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	READ(NFLCH) MAXC_STIFFNESS(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
      ! ------------------------
      
      ! ------- MASS --------
	CALL MCALFIL(NFL,"MASS") 
	REWIND(NFL) 
 
	CALL MCALFIL(NFLCH,'MAXC') 
	REWIND(NFLCH)  
 
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0,'FIRS') 
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0,'FIRS') 
     
      NEQF = 1 
      NEQL = 0 
	DO 800 IBLO = 1,NBLOCK  
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0,'FIRS') 
	CALL MINTFIL('BLOC',NHIG_MASS,3,IBLO,0,'FIRS') 
      NEQL = NEQL+NCOL 
      
	IF(OPER.EQ.'STD') TT(NEQF:NEQL) = 0.0D0  !IF(OPER.EQ.'ADD') ---> NO INITIALIZATION    

	NEQF = NEQF + NCOL 
800   CONTINUE 
      
      ALLOCATE(MAXC_MASS(NWK),B(NWK)) 
C     --------------------------------------
      NEQF = 1 
      NEQL = 0 
      IAA = 1
	ICC = 1
	DO 1000 IBLO = 1,NBLOCK 
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0,'FIRS') 
	CALL MINTFIL('BLOC',NHIG_MASS,3,IBLO,0,'FIRS') 
      NEQL = NEQL+NCOL 
	NPRE = MAXA(NEQF) - 1 
 
      READ(NFL  )    B(IAA:IAA+NHIG_MASS-1)
	READ(NFLCH)    MAXC_MASS(ICC:ICC+NHIG_MASS-1)
      IAA = IAA + NHIG_MASS
	ICC = ICC + NHIG_MASS
      !READ(NFL  )    B(1:NHIG_MASS)
	!READ(NFLCH) MAXC_MASS(1:NHIG_MASS)

1000  CONTINUE
C     --------------------------------------
      
C	------------------------------
C	FORM MASS PROPORTIONAL DAMPING
C	------------------------------  
      NCOMPLEX = 0
      
      IF (ALFA.EQ.0) GOTO 200
      NCOMPLEX = 1
C	-----------------------------------
C	FORM STIFFNESS PROPORTIONAL DAMPING
C	-----------------------------------
200   IF (BETA.EQ.0) GOTO 210
      AA = AA*BETA
      !------------------------
      
      ! DETERMINE TRANSFER FUNCTION
210   IF (OPRES.EQ.1) CALL OFFSHSTEP(TIME,ITIME,NINC,'CALT') !READ TIME DATA FOR OFFSHORE LOAD 
      IF (OPRES.EQ.2) CALL SELECTGRAPH (1,NINC,0,0)
      
      ALLOCATE(SW(NWK,NINC),SW_COMPLEX(NWK,NINC)) 
      ALLOCATE(X(NEQ),X_COMPLEX(NEQ),CC(NWK),CC_COMPLEX(NWK),IPERM(NEQ),LPT(64),V_COMPLEX(NEQ)) 
      SW = 0.
      SW_COMPLEX = 0.
      DO  I =1,NINC/2D0
          
      IF (OPRES.EQ.1.AND.I.EQ.1) THEN
      FREQ =  STEPSTAT   
      WWW  = 2D0*PI*FREQ 
      ELSEIF (OPRES.EQ.1.AND.I.GT.1) THEN
      FREQ = FREQ + STEPINCR
      WWW  = 2D0*PI*FREQ 
      ENDIF
      
      DO J = 1,NWK
        ! remain complex form
        IF (NCOMPLEX.EQ.0) SW(J,I) =  -(WWW**2D0)*B(J)+AA(J)
        IF (NCOMPLEX.EQ.1) THEN
          REAL_VALUE     = -(WWW**2D0)*B(J)+AA(J)
          AIM_VALUE      = WWW*ALFA
          SW_COMPLEX(J,I) = CMPLX(REAL_VALUE,AIM_VALUE)
       ENDIF
        ENDDO
        
        !IF (I.EQ.1)THEN
        !    DO III = 1,NWK
        !    WRITE (252,*) SW(III,I)
        !    ENDDO
        !ENDIF
      

        
      ! ------ FORCE MATRIX -------
	RFV(1:NEQ) = 0.0D0
	RFC(1:NEQ) = 0.0D0
      
      ! GENERAL LOADCASE FOR OFFSHORE (OLD VERSION BEFORCE 04/2018)
      CALL OFFSHFORC(RFV,1,I,NEQ,'VARY','RADD') !VARY OFFSHORE LOAD     -- LOADCASE 1
      CALL OFFSHFORC(RFC,1,I,NEQ,'CONT','RADD') !CONSTANT OFFSHORE LOAD -- LOADCASE 1
      
      REXTERNAL = 0.
      CALL RHSDYNA (A(LLF),A(LLV),A(LLO),RFV,RFC,REXTERNAL,A(LEF),A(LPR),A(LRT1),NEQ,2)
            
      !TESTING FILE
      A(LEF:LEF+NEQ) = 0.0D0
      V_COMPLEX = 0.0
      READ (256,*) A(6817)
      A(6817) = A(6817)*NINC
      
      V(1:NEQ) = A(LEF:LEF+NEQ)/NINC
        
      
      ! ------ PARDISO ------
      IF (NCOMPLEX.EQ.0)THEN ! REAL SOLVER
      X = 0.
C.. Set up PARDISO control parameter
C..
50    NRHS = 1
      MAXFCT = 1
      MNUM = 1
      
!      DO KK = 1,64
!         IPARM(KK) = 0
!      ENDDO
!      IPARM(1)  =  1 ! no solver default
!      IPARM(2)  =  2 ! 2 = fill-in reordering from METIS , 3 = parallel (OpenMP) nested dissection algorithm 
!      IPARM(3)  =  1 !MKL_GET_MAX_THREADS() !1 ! numbers of processors
!      IPARM(4)  =  0 ! no iterative-direct algorithm
!      IPARM(5)  =  0 ! no user fill-in reducing permutation
!      IPARM(6)  =  0 ! =0 solution on the first n compoments of x
!      IPARM(7)  =  0 ! not in use
!      IPARM(8)  =  9 ! numbers of iterative refinement steps
!      IPARM(9)  =  0 ! not in use
!      IPARM(10) =  13 ! perturbe the pivot elements with 1E-13
!      IPARM(11) =  1 ! use nonsymmetric permutation and scaling MPS
!      IPARM(12) =  0 ! not in use
!      IPARM(13) =  0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try IPARM(13) = 1 in case of inappropriate accuracy
!      IPARM(14) =  0 ! Output: number of perturbed pivots
!      IPARM(15) =  0 ! not in use
!      IPARM(16) =  0 ! not in use
!      IPARM(17) =  0 ! not in use
!      IPARM(18) = -1 ! Output: number of nonzeros in the factor LU
!      IPARM(19) = -1 ! Output: Mflops for LU factorization
!      IPARM(20) =  0 ! Output: Numbers of CG Iterations
!      IPARM(27) =  0 ! Output: Numbers of CG Iterations
!      IPARM(60) =  1 ! 
!      IERROR    =  0 ! initialize error flag
!      MSGLVL    =  0 ! print statistical information
!      MTYPE     =  -2 ! symmetric, indefinite
   
      !ALLOCATE( LPT(64) )
!      DO KK = 1, 64
!      LPT(KK) = LPTO(KK) 
!      ENDDO
      
      !CALL CPU_TIME (TIME1)
      IF (NCOMPLEX.EQ.0) CC(1:NWK)         = SW(1:NWK,I) ! REAL 
      IF (NCOMPLEX.EQ.1) CC_COMPLEX(1:NWK) = SW_COMPLEX(1:NWK,I) ! COMPLEX

      
      IPARM(8) = 2 ! max numbers of iterative refinement steps
      IPHASE = 33 ! back substitution and iterative refinement
      
!      CC_TEST(1) = 10D0
!      CC_TEST(2) = 0.0D0
!      CC_TEST(3) = 0.0D0
!      CC_TEST(4) = 20D0
!      CC_TEST(5) = 0.0D00
!      CC_TEST(6) = 30D0
!      
!      MAXA_TEST(1) = 1
!      MAXA_TEST(2) = 4
!      MAXA_TEST(3) = 6
!      
!      MAXC_TEST(1) = 1
!      MAXC_TEST(2) = 2 
!      MAXC_TEST(3) = 3
!      MAXC_TEST(4) = 2
!      MAXC_TEST(5) = 3
!      MAXC_TEST(6) = 3
!      
!      V_TEST(1) = 1D0
!      V_TEST(2) = 1D0
!      V_TEST(3) = 1D0
!      
!      NEQ_TEST = 3
      
       !DO IJK = 1,NWK
       !    READ (256,*) CC(IJK)
       !ENDDO
                   
!      CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ_TEST,CC_TEST,MAXA_TEST,MAXC_TEST,
!     1 IPERM,NRHS,IPARM,MSGLVL,V_TEST,X_TEST,IERROR)
      CALL FULL_MATRIX (CC(1:NWK),BB_OUT,MAXA,MAXC_STIFFNESS,NEQ,NWK)
      
      ! FINAL DISPLACEMENT
      !V = 1D0
      X = MATMUL(BB_OUT,V)
      
!       CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,CC,MAXA,MAXC_STIFFNESS,
!     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
      
      CALL STORE_RESULT_FREQ (X,X_COMPLEX,I,NEQ,NCOMPLEX,"DISP")
     
      !X = 0.D0
      !V = 1D0
      !DD = 0.
      
      ! DEFINE TRANSFER FUNCTION
      ! CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,DD,MAXA,MAXC_STIFFNESS,
      !1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
      
      ! TESTING INV
      !DO III = 1,NWK
      !    AIDEN = CC(III)*X(III)
      !ENDDO
      
      !CALL STORE_RESULT_FREQ (X,X_COMPLEX,I,NEQ,NCOMPLEX,"TRAN")
      
      ELSEIF (NCOMPLEX.EQ.1)THEN ! COMPLEX SOLVER
      ! DISPLACEMENT    
      V_COMPLEX = V
      
      CALL FULL_MATRIX_COMPLEX (CC_COMPLEX(1:NWK),BB_OUT_COMLEX,MAXA,MAXC_STIFFNESS,NEQ,NWK)
      CALL CMATINV(NEQ,NEQ,BB_OUT_COMLEX,V_COMPLEX,DET)
      CALL STORE_RESULT_FREQ (X,V_COMPLEX,I,NEQ,NCOMPLEX,"DISP")    
      ENDIF
      
      ENDDO
      
      
      DEALLOCATE(AA,MAXC_STIFFNESS,MAXC_MASS,SW,SW_COMPLEX,X,X_COMPLEX,CC,CC_COMPLEX,IPERM,LPT,V_COMPLEX)
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE STIFFNESS_BLOCK (MAXA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2  
      DIMENSION MAXA(1)
      ALLOCATABLE AA(:),MAXC(:) 
      
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
C     IPFAC DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT
      TYP2 = "TEMP"
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,0,'FIRS')  

      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0,'FIRS')
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0,'FIRS')
      
	CALL MCALFIL(KSREC,TYP2)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
      
      ALLOCATE(AA(NWK),MAXC(NWK)) 
	IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0,'FIRS')
      
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	READ(NFLCH) MAXC(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
      DEALLOCATE(AA,MAXC)
      
      
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE MASS_BLOCK (MAXA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION MAXA(1),B(1),RR(1),TT(1) 
	ALLOCATABLE MAXC(:) 
C	-------------------------------------------
	CALL MCALFIL(NFL,"MASS") 
	REWIND(NFL) 
 
	CALL MCALFIL(NFLCH,'MAXC') 
	REWIND(NFLCH)  
 
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0,'FIRS') 
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0,'FIRS') 
     
      NEQF = 1 
      NEQL = 0 
	DO 800 IBLO = 1,NBLOCK  
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0,'FIRS') 
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0,'FIRS') 
      NEQL = NEQL+NCOL 
      
	IF(OPER.EQ.'STD') TT(NEQF:NEQL) = 0.0D0  !IF(OPER.EQ.'ADD') ---> NO INITIALIZATION    

	NEQF = NEQF + NCOL 
800	CONTINUE 
C     --------------------------------------
      NEQF = 1 
      NEQL = 0 
	DO 1000 IBLO = 1,NBLOCK 
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0,'FIRS') 
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0,'FIRS') 
      NEQL = NEQL+NCOL 
	NPRE = MAXA(NEQF) - 1 
 
	ALLOCATE(MAXC(NHIG)) 

      CALL READING_VALUE (NFL,NFLCH,NHIG) 

	DEALLOCATE(MAXC)

1000  CONTINUE
C     --------------------------------------
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE READING_VALUE (NFL,NFLCH,NHIG)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION B(NHIG),MAXC(NHIG)
	READ(NFL  )    B(1:NHIG)
	READ(NFLCH) MAXC(1:NHIG)
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE STORE_RESULT_FREQ (RE_SULT,RE_SULT_COMPLEX,KK,NEQ,NCOMPLEX,OPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMPLEX ::RE_SULT_COMPLEX(NEQ)
      DIMENSION  RE_SULT(NEQ)
      CHARACTER*4 OPT

      ! FILE INDEX HAS BEEN DEFINED IN ASSIGN_XSEA_TEXT_DATA.F90
      !OPEN(UNIT=251,FILE='FREQ DOMAIN/Transfer Function.dat'      ,STATUS='UNKNOWN')
      !OPEN(UNIT=252,FILE='FREQ DOMAIN/Displacment.dat'            ,STATUS='UNKNOWN')
      
      ! EXTERNAL STORAGE 
      IF (OPT.EQ."DISP") THEN      ! DISPLACEMENT
      WRITE (252,1) KK             ! HEADER IS INDICATE STEP OF FREQ.
      DO I = 1,NEQ
      IF (NCOMPLEX.EQ.1) WRITE (252,2)  RE_SULT_COMPLEX(I)    ! COMPLEX RESULT 
      IF (NCOMPLEX.EQ.0) WRITE (252,2)  RE_SULT(I)    ! RESULT 
      ENDDO
1     FORMAT (I5)
2     FORMAT (E25.10)
      
      ELSEIF (OPT.EQ."TRAN") THEN ! TRANSFER FUNCTION
      WRITE (251,100) KK        ! HEADER IS INDICATE STEP OF FREQ.
      DO I = 1,NEQ
      IF (NCOMPLEX.EQ.1) WRITE (251,101)  RE_SULT_COMPLEX(I)  !COMPLEX RESULT 
      IF (NCOMPLEX.EQ.0) WRITE (251,101)  RE_SULT(I)  ! RESULT 
      ENDDO
100   FORMAT (I5)
101   FORMAT (E25.10)
      ENDIF
      
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE COMBINE_FREQ_DATA
      USE MKL_DFTI
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      !ALLOCATABLE DISP(:),ATRAN(:),SUM_DISP(:),SUM_TRAN(:)
      ALLOCATABLE DISP(:,:),ATRAN(:,:),DATA_OUTPUT(:)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      !POINTER :: My_Desc1_Handle, My_Desc2_Handle
      !DIMENSION DATA_INPUT(18)
      ! RE-START DATA
      REWIND (251)
      REWIND (252) 
      
      IF (OPRES.EQ.1) CALL OFFSHSTEP(TIME,ITIME,NINC,'CALT') !READ TIME DATA FOR OFFSHORE LOAD 
      IF (OPRES.EQ.2) CALL SELECTGRAPH (1,NINC,0,0)
      
      !ALLOCATE (DISP(NEQ),ATRAN(NEQ),SUM_DISP(NEQ),SUM_TRAN(NEQ))
      ALLOCATE (DISP(NEQ,NINC),ATRAN(NEQ,NINC),DATA_OUTPUT(NEQ))
      ATRAN    = 0.
      DISP     = 0.
      SUM_DISP = 0.
      SUM_TRAN = 0.
      ! 
      
      DO 100 I =1,NINC/2D0
          !READ (251,*) LHEAD_251
          READ (252,*) LHEAD_252
          !IF (I.NE.LHEAD_251) EXIT
          IF (I.NE.LHEAD_252) EXIT
          DO J = 1,NEQ
          !READ (251,*) ATRAN(J,I)
          READ (252,*) DISP(J,I)
          !SUM_TRAN(J) = SUM_TRAN(J) + ATRAN(J)
          !SUM_DISP(J) = SUM_DISP(J) + DISP(J)
          ENDDO
          
100   CONTINUE
      
      !OPEN(UNIT=253,FILE='FREQ DOMAIN/Displacement and Transfer Function.dat'           ,STATUS='UNKNOWN')
      WRITE (253,200) 
      DO J = 1,NEQ
      DO I = 1,NINC/2D0
      IF (I.EQ.1) STEP = STEPSTAT
      IF (I.GT.1) STEP = STEP + STEPINCR
      WRITE (253,201) STEP,J,DISP(J,I)!,ATRAN(J,I)
      ENDDO
      ENDDO
200   FORMAT ("     FREQ.            NEQ      DISPLACMENT ")
201   FORMAT (F17.5,2X,I5,2X,E25.10)
!200   FORMAT ("     FREQ.            NEQ      DISPLACMENT     TRANSFER_FUNCTION ")
!201   FORMAT (F17.5,2X,I5,2X,E25.10,2X,E25.10)
      
      !CALL FAST_FOURIER_BACKWARD (NINC)
      
      !DEALLOCATE (DISP,ATRAN,SUM_DISP,SUM_TRAN)
      DEALLOCATE (DISP,ATRAN,DATA_OUTPUT)
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE PRINT_OUTPUT_DISPLACEMENT (NINC,NEQ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      !DIMENSION DATA_OUTPUT_REAL(NINC,NEQ)
      ! FOR STRESS AND INTERANAL FORCE
      
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE REPLACE_RESULT (NTSTEP,DISP,OPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPT
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
       COMMON /TIME/ DDT,CTIM,NINC
       DIMENSION DATA_OUTPUT_REAL(NINC,NEQ)
       DIMENSION DISP(NEQ)
      
       IF (OPT.EQ."DISP")THEN
          DATA_OUTPUT_REAL = 0.D0
         REWIND (254)
         ! EXTERNAL STORAGE 
         DO I = 1,NEQ
             DO J = 1,NINC
               READ (254,*) DATA_OUTPUT_REAL(J,I)
             ENDDO
         ENDDO
         
         DISP  = 0.
         DO KK = 1, NEQ
         DISP(KK) = DATA_OUTPUT_REAL(NTSTEP,KK)
         ENDDO
      
       ELSEIF (OPT.EQ."FORC")THEN
         DATA_OUTPUT_REAL = 0.0D0
         REWIND (255)
         ! EXTERNAL STORAGE 
         DO I = 1,NEQ
             DO J = 1,NINC
               READ (255,*) DATA_OUTPUT_REAL(J,I)
             ENDDO
         ENDDO
         
         DISP  = 0.
         DO KK = 1, NEQ
         DISP(KK) = DATA_OUTPUT_REAL(NTSTEP,KK)
         ENDDO
           
      ENDIF

      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE IFFT_FORCE 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /LNDEQK/ LDEK,LPR2,LDESTP
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER 
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON A(7000000),IA(6000000)
      
      DIMENSION FORCE_FREQ(NEQ,NINC)
      DIMENSION FORCE_FREQ_INPUT(NINC)
      DIMENSION RFV(NEQ),RFC(NEQ)
      
      DO  INC_FORCE = 1,NINC
      IF(LDEK.EQ.0) THEN
	RFV(1:NEQ) = 0.0D0
	RFC(1:NEQ) = 0.0D0
      
      ! GENERAL LOADCASE FOR OFFSHORE (OLD VERSION BEFORCE 04/2018)
      CALL OFFSHFORC(RFV,1,INC_FORCE,NEQ,'VARY','RADD') !VARY OFFSHORE LOAD     -- LOADCASE 1
      CALL OFFSHFORC(RFC,1,INC_FORCE,NEQ,'CONT','RADD') !CONSTANT OFFSHORE LOAD -- LOADCASE 1
      CALL FASTTOPFORCE_DUMMY (IFAST,KFAST,NODEFAST,NTFASTPARA,NUMCASE,NFASTPARA)
      IF (KFAST.EQ.2) CALL FASTTOPFORCE (IA(LID),NSF,NODEA,NUMCASE,NFASTPARA,NTFASTPARA,RFC,NEQ,INC,'REDD')
      
      IF(KOREV.EQ.1) THEN  ! REALATIVE MOTION OF OFFSHORE 
      CALL RELATIVEMOTION (IA(LLM),A(LCO),IA(LGS),A(LGP),IA(LMS),A(LMP),IA(LGID),MFRMLD,IA(LHG),IA(LHGN),INC,A(LVL),A(LAL),RFV) 
      ENDIF
      
      REXTERNAL = 0.
      CALL RHSDYNA (A(LLF),A(LLV),A(LLO),RFV,RFC,REXTERNAL,A(LEF),A(LPR),A(LRT1),NEQ,2)
      FORCE_FREQ(INC_FORCE,1:NEQ) = A(LEF:LEF+NEQ)
      
      ENDIF
      ENDDO
      
      DO II = 1,NEQ
      FORCE_FREQ_INPUT(1:NINC) = FORCE_FREQ(1:NINC,II)
      CALL FAST_FOURIER_BACKWARD_FORCE (FORCE_FREQ_INPUT,NEQ,NINC)
      ENDDO
      
      
      END
C	=====================================================================
C	===================================================================== 
      SUBROUTINE FULL_MATRIX (CC,BB,MAXA,MAXC,NEQ,NWK)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION CC(NWK),MAXA(NEQ),MAXC(NWK),MAXG(NWK)
      DIMENSION AA(NEQ,NEQ),BB(NEQ,NEQ)
      ! CC = INPUT
      ! CHECKING MAXC
      MAXG        = 0.
      INDEX       = 1
      DO I = 1,NWK
          IF (I.EQ.1) THEN
             MAXG(I) = INDEX
          ELSEIF (I.NE.1)THEN
            FIST_INDEX = MAXC(I-1)
            SECOND_INDEX = MAXC(I)
            IF (SECOND_INDEX.GT.FIST_INDEX) THEN
            MAXG(I) = INDEX    
            ELSEIF (SECOND_INDEX.LE.FIST_INDEX) THEN
            INDEX = INDEX + 1
            MAXG(I) = INDEX    
            ENDIF
          ENDIF
      ENDDO
      AA = 0.
      INDEX_DIA = 1D0
      DO I =1,NWK
          !IF (I.EQ.MAXA(MAXG(I)))THEN
          !   AA(MAXG(I),MAXC(I)) = CC(I)  ! UPPER TRI
          !   !IF (I.NE.MAXA(MAXG(I)))THEN
          !   !IF (MAXG(I).EQ.1) IMAXC  = NEQ
          !   !IF (MAXG(I).NE.1) IMAXC  = NEQ - MAXG(I) + 1.0D0
          !   !AA(IMAXC,1) = CC(I) ! LOWER TRI
          !   !ENDIF
          !IF (I.NE.MAXA(MAXG(I)))THEN
          IF (I.EQ.MAXA(I))THEN
          AA(INDEX_DIA,INDEX_DIA) = CC(I)  ! UPPER TRI
          INDEX_DIA = INDEX_DIA + 1
          ELSEIF (I.NE.MAXA(I))THEN
             AA(MAXG(I),MAXC(I)) = CC(I)  ! UPPER TRI
             
             JINDEX = NEQ-MAXC(I)+1
             IF (MAXG(I).EQ.1.0D0) IMAXC  = NEQ
             IF (MAXG(I).GT.1.0D0) IMAXC  = NEQ - MAXG(I) + 1 
             AA(IMAXC,JINDEX) = CC(I) ! LOWER TRI
          ENDIF
      ENDDO
      ! INPUT  AA
      ! OUTPUT BB
      BB = 0.D0
      CALL inverse (AA,BB,NEQ)
      END
      
C	===================================================================== 
      SUBROUTINE FULL_MATRIX_COMPLEX (CC,BB,MAXA,MAXC,NEQ,NWK)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMPLEX CC(NWK),AA(NEQ,NEQ),BB(NEQ,NEQ)
      DIMENSION MAXA(NEQ),MAXC(NWK),MAXG(NWK)
      !DIMENSION AA(NEQ,NEQ),BB(NEQ,NEQ)
      ! CC = INPUT
      ! CHECKING MAXC
      MAXG        = 0.
      MAXG(1:NEQ) = 1.
      INDEX       = 2
      DO I = NEQ+1,NWK,1
          IF (MAXC(I).EQ.NEQ)THEN
          MAXG(I) = INDEX
          INDEX = INDEX + 1
          ELSEIF (MAXC(I).NE.NEQ)THEN
          MAXG(I) = INDEX 
          ENDIF
      ENDDO
      AA = 0.
      INDEX = 0
      DO I =1,NWK
          INDEX = (I/NEQ) + 1D0
          IF (MAXC(I).EQ.NEQ)THEN
             AA(MAXG(I),MAXC(I)) = CC(I)  ! UPPER TRI
             IF (I.NE.MAXA(MAXG(I)))THEN
             IF (MAXG(I).EQ.1) IMAXC  = NEQ
             IF (MAXG(I).NE.1) IMAXC  = NEQ - MAXG(I) + 1.0D0
             AA(IMAXC,1) = CC(I) ! LOWER TRI
             ENDIF
          ELSEIF (MAXC(I).NE.NEQ)THEN
          AA(MAXG(I),MAXC(I)) = CC(I)  ! UPPER TRI
             IF (I.NE.MAXA(MAXG(I)))THEN
             JINDEX = NEQ-MAXC(I)+1
             IMAXC  = NEQ-INDEX+1
             AA(IMAXC,JINDEX) = CC(I) ! LOWER TRI
             ENDIF
          ENDIF
      ENDDO
      ! INPUT  AA
      ! OUTPUT BB
      BB = 0.D0
      !CALL inverse (AA,BB,NEQ)
      END
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !===========================================================
      subroutine inverse(a,c,n)
      implicit none 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k
      
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0
      
      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do
      
      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do
      
      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse
      
