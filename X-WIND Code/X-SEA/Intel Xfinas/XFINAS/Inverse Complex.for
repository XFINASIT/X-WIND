!***********************************************************************************************
!*            SOLVING A COMPLEX LINEAR MATRIX SYSTEM AX=B BY GAUSS-JORDAN METHOD               *
!* ------------------------------------------------------------------------------------------- *
!*                                                          F90 Version By J-P Moreau, Paris   *
!*                                                                 (www.jpmoreau.fr)           *
!* SAMPLE RUN:                                                                                 *
!* INPUTS  from example file csysmat.dat,                                                      *
!* OUTPUTS to file csysmat.lst.                                                                *
!* ------------------------------------------------------------------------------------------- *
!* SAMPLE RUN:                                                                                 *
!* Input file                                                                                  *
!* 4                                                                                           *
!* 47.0 -15.0  62.0 5.0    0.0 -72.0 61.0  20.0                                                *
!*  6.0  14.0 -17.0 3.0 -102.0  91.0  7.0 -12.0                                                *
!* 13.0 -55.0  32.0 8.0   41.0   7.0 25.0   1.0                                                *
!*111.0  25.0  40.0 0.0   12.0 -82.0 58.0 -30.0                                                *
!* 1                                                                                           *
!* 629.0 988.0                                                                                 *
!*-180.0 825.0                                                                                 *
!* 877.0 441.0                                                                                 *
!* 734.0 -88.0                                                                                 *
!*                                                                                             *
!* Output file                                                                                 *
!* ------------------------------------------------------                                      *
!*   SOLVING THE COMPLEX MATRIX LINEAR SYSTEM AX = B                                           * 
!* ------------------------------------------------------                                      *
!*    (BY GAUSS-JORDAN METHOD WITH FULL PIVOTING)                                              *
!*                                                                                             *
!*  N=  4                                                                                      *
!*                                                                                             *
!*  COMPLEX MATRIX A                                                                           *
!*                                                                                             *
!*(  47.0000, -15.0000) (  62.0000,   5.0000) (   0.0000, -72.0000) (  61.0000,  20.0000)      *
!*(   6.0000,  14.0000) ( -17.0000,   3.0000) (-102.0000,  91.0000) (   7.0000, -12.0000)      *
!*(  13.0000, -55.0000) (  32.0000,   8.0000) (  41.0000,   7.0000) (  25.0000,   1.0000)      *
!*( 111.0000,  25.0000) (  40.0000,   0.0000) (  12.0000, -82.0000) (  58.0000, -30.0000)      *
!*                                                                                             *
!*  M=  1                                                                                      *
!*                                                                                             *
!*  COMPLEX MATRIX B                                                                           *
!*                                                                                             *
!*( 629.0000, 988.0000)                                                                        *
!*(-180.0000, 825.0000)                                                                        *
!*( 877.0000, 441.0000)                                                                        *
!*( 734.0000, -88.0000)                                                                        *
!*                                                                                             *
!*                                                                                             *
!*  INVERSE OF MATRIX A                                                                        *
!*                                                                                             *
!*(   0.0039,  -0.0080) (  -0.0084,  -0.0041) (  -0.0095,   0.0179) (   0.0002,  -0.0017)      *
!*(   0.0368,   0.0218) (   0.0148,  -0.0271) (  -0.0338,  -0.0351) (  -0.0050,  -0.0162)      *
!*(  -0.0042,  -0.0051) (  -0.0077,  -0.0018) (   0.0035,   0.0048) (   0.0007,   0.0036)      *
!*(  -0.0197,  -0.0166) (  -0.0015,   0.0188) (   0.0337,   0.0154) (   0.0054,   0.0172)      *
!*                                                                                             *
!*                                                                                             *
!*  SOLUTION MATRIX X                                                                          *
!*                                                                                             *
!*(  -1.0000,   3.0000)                                                                        *
!*(   2.0000,  10.0000)                                                                        *
!*(   7.0000,  -5.0000)                                                                        *
!*(  17.0000,   6.0000)                                                                        *
!*                                                                                             *
!*  DETERMINANT= ( 0.1741656E+08,-0.1059832E+08)                                               *
!*                                                                                             *
!*                                                                                             *
!*  VERIFICATION A*X=B                                                                         *
!*                                                                                             *
!*( 628.9999, 988.0001)                                                                        *
!*(-180.0000, 824.9999)                                                                        *
!*( 876.9999, 440.9999)                                                                        *
!*( 733.9999, -87.9999)                                                                        *
!* ------------------------------------------------------                                      *
!*                                                                                             *
!***********************************************************************************************
!      SUBROUTINE CSYSMAT
! 
!      COMPLEX, POINTER :: A(:,:),B(:,:),A1(:,:)                                             
!      COMPLEX, POINTER :: C(:,:),D(:,:)
!      COMPLEX, POINTER :: TEMP(:)
!      COMPLEX  DETER
!	  
!	  REAL, POINTER :: X(:), Y(:) 
!
!      CHARACTER*40 INPUT,OUTPUT
!      CHARACTER*36 NOM
!
!      WRITE(*,*) ' '
!      WRITE(*,*) ' SOLVING A LINEAR MATRIX SYSTEM AX=B'
!      WRITE(*,17,advance='no'); READ (*,*) NOM
!      J=0
!      DO I=1,LEN(NOM)
!	    IF(NOM(I:I)<>' ') J=J+1  !J=real length of NOM
!      ENDDO
!      
!      OUTPUT=NOM(1:J)//'.lst'
!      INPUT =NOM(1:J)//'.dat'
!
!      !Open input, output files
!      OPEN(UNIT=2,FILE=OUTPUT,STATUS='UNKNOWN')	                              
!      OPEN(UNIT=1,FILE=INPUT,STATUS='OLD')	                                                                                                                                                                                                                                                      
!                                                                               
!!Read matrix A to be inverted and copy to output file
!!Utility read/write subroutines are defined in optional module BASIS_R
!
!      READ (1,*)  N   !Size of system to be solved
!
!      !dynamic allocation of matrix A and utility vector TEMP
!      allocate(A(1:N,1:N),stat=ialloc)
!      allocate(TEMP(1:N),stat=ialloc)
!      allocate(X(1:N),stat=ialloc)
!      allocate(Y(1:N),stat=ialloc)
!
!      WRITE(2,*) '------------------------------------------------------'
!      WRITE(2,*) '  SOLVING THE COMPLEX MATRIX LINEAR SYSTEM AX = B'      
!      WRITE(2,*) '------------------------------------------------------'
!      WRITE(2,*) '   (BY GAUSS-JORDAN METHOD WITH FULL PIVOTING)' 
!  
!      WRITE(2,12)  N
!
!      DO I=1, N
!	    READ(1,*) (X(J), Y(J), J=1,N)
!		DO J=1,N
!		  A(I,J) = CMPLX(X(J),Y(J))
!        END DO
!      END DO
!
!      WRITE(2,15)
!      WRITE(2,80) ((A(I,J),J=1,N),I=1,N)
!                                                                               
!!Read right hand matrix B if M<>0                                              
!     
!      READ (1,*) M
!   
!      WRITE(2,13)  M                                                                       
!      IF(M.NE.0) THEN     
!!dynamic allocation of matrix B
!   	    allocate(B(1:N,1:M),stat=ialloc)
!
!        DO I=1, N
!	      READ(1,*) (X(J), Y(J), J=1,M)
!		  DO J=1,M
!		    B(I,J) = CMPLX(X(J),Y(J))
!          END DO
!        END DO
!                                  
!        WRITE(2,20)                                                             
!  		WRITE(2,81) ((B(I,J),J=1,M),I=1,N)                                  
!      ENDIF     
!      CLOSE(1)                                                                
!                                                                               
!!Store matrix A in utility matrix A1 for verifications purpose
!! (Original A is destroyed during call to inversion process)                                                     
!                                                            
!      !dynamic allocation of utility matrix A1															                  
!      allocate(A1(1:N,1:N),stat=ialloc)
!      A1=A                                                                                                                       
!                                                                               
!!Call matrix inversion subroutine                                                                                                  
!      CALL CMATINV(N,M,A,B,DETER)                                                    
!                                                                               
!!Print results to output file                                                          
!                                                                               
!      WRITE(2,30)                                        	                                                             
!      WRITE(2,80) ((A(I,J),J=1,N),I=1,N)
!                                                                               
!      IF(M.NE.0) THEN                                                           
!        WRITE(2,40)
!		WRITE(2,81) ((B(I,J),J=1,M),I=1,N)                                    
!      ENDIF                                                                     
!                                                                               
!      WRITE(2,50) DETER                                                         
!                                                                               
!! Verification that new product A*B equals original B matrix               
!! (optional)                                                                
!                                                                               
!      IF(M.NE.0) THEN                                                           
!        WRITE(2,60)                                                             
!!call matrices multiplication subroutine   
!        allocate(C(1:N,1:M),stat=ialloc)                                      
!        CALL MATMUL(A1,B,C,N,M)
!        WRITE(2,81) ((C(I,J),J=1,M),I=1,N)		                                  
!      ELSE                                                                      
!        WRITE(2,70)                                                             
!!call matrices multiplication subroutine
!        allocate(D(1:N,1:N),stat=ialloc)                                         
!        CALL MATMUL(A1,A,D,N,N)
!		WRITE(2,80) ((D(I,J),J=1,N),I=1,N)                                 
!      ENDIF                                                                     
!                                                                         
!      WRITE(2,*) '------------------------------------------------------'																		       
!      CLOSE(2)
!       
!      WRITE(*,*)
!      WRITE(*,*) ' Results in file ',OUTPUT
!	  WRITE(*,*) ' '
!      STOP                                                                      
!                                                                               
!   12 FORMAT(/'  N=',I3)
!   13 FORMAT(/'  M=',I3)               
!   15 FORMAT(/'  COMPLEX MATRIX A '/)
!   17 FORMAT(/'  Input file name (without .dat): ')                                                  
!   20 FORMAT(/'  COMPLEX MATRIX B '/)                                                     
!   30 FORMAT(//'  INVERSE OF MATRIX A'/)                                       
!   40 FORMAT(//'  SOLUTION MATRIX X'/)                                         
!   50 FORMAT(/'  DETERMINANT= (',E14.7,',',E14.7,')'/)                                      
!   60 FORMAT(/'  VERIFICATION A*X=B'/)                                         
!   70 FORMAT(/'  VERIFICATION A1*A=I'/)  
!   80 FORMAT(4('(',F9.4,',',F9.4,') '))    
!   81 FORMAT('(',F9.4,',',F9.4,') ')                                         
!      END                                                                       
                                                                               
!***********************************************
!* SOLVING A COMPLEX LINEAR MATRIX SYSTEM AX=B *
!* with Gauss-Jordan method using full pivoting*
!* at each step. During the process, original  *
!* A and B matrices are destroyed to spare     *
!* storage location.                           *
!* ------------------------------------------- *
!* INPUTS:    A    COMPLEX MATRIX N*N          *
!*            B    COMPLEX MATRIX N*M          *
!* ------------------------------------------- *
!* OUTPUTS:   A    INVERSE OF A N*N            *
!*            DET  COMPLEX DETERMINANT OF A    *
!*            B    SOLUTION MATRIX N*M         *
!* ------------------------------------------- *
!* NOTA - If M=0 inversion of A matrix only.   *
!***********************************************
      SUBROUTINE CMATINV(N,M,AA,BB,DET)
      REAL, PARAMETER :: EPSMACH=2.E-12
      COMPLEX AA(N,N),BB(N,M)
      INTEGER, POINTER :: PC(:), PL(:)
	  COMPLEX, POINTER :: CS(:)
      COMPLEX PV, DET, TT
	  REAL PAV
	                                                    
!Initializations :                       
      allocate(PC(1:N),stat=ialloc)
      allocate(PL(1:N),stat=ialloc)
	  allocate(CS(1:N),stat=ialloc)	                                        
      DET=CMPLX(1.0,0.0)
      DO I=1,N                                                                                                                   
        PC(I)=0
	    PL(I)=0
	    CS(I)=CMPLX(0.0,0.0)
      END DO
!main loop                                                                   
      DO K=1,N                                                                  
!Searching greatest pivot:                                               
        PV=AA(K,K)                                                              
        IK=K                                                                    
        JK=K                                                                    
        PAV=ABS(PV)                                                            
        DO I=K,N                                                                
          DO J=K,N                                                              
            IF (ABS(AA(I,J)).GT.PAV) THEN                                      
              PV=AA(I,J)                                                        
              PAV=ABS(PV)                                                      
              IK=I                                                              
              JK=J                                                              
            ENDIF                                                               
          ENDDO                                                                 
        ENDDO                                                                   
                                                                               
!Search terminated, the pivot is in location I=IK, J=JK.
!Memorizing pivot location: :                                        
        PC(K)=JK                                                                
        PL(K)=IK                                                                
                                                                               
!Determinant DET is actualised
!If DET=0, ERROR MESSAGE and STOP
!Machine dependant EPSMACH equals here 2.E-12                                        
                                                                               
        IF (IK.NE.K) DET=-DET                                                   
        IF (JK.NE.K) DET=-DET                                                   
        DET=DET*PV                                                              
        IF (ABS(DET).LT.EPSMACH) THEN                                          
!Error message and Stop                                                     
          PRINT 10                                                              
          STOP                                                                  
        ENDIF                                                                   
                                                                               
!POSITIONNING PIVOT IN K,K:                                              
        IF(IK.NE.K) THEN                                                        
          DO I=1,N                                                              
!EXCHANGE LINES IK and K:                                            
            TT=AA(IK,I)                                                         
            AA(IK,I)=AA(K,I)                                                    
            AA(K,I)=TT                                                          
          ENDDO                                                                 
        ENDIF                                                                   
      IF (M.NE.0) THEN                                                          
        DO I=1,M                                                               
          TT=BB(IK,I)                                                           
          BB(IK,I)=BB(K,I)                                                      
          BB(K,I)=TT                                                            
        ENDDO                                                                   
      ENDIF                                                                     
!Pivot is at correct line                                                
        IF(JK.NE.K) THEN                                                        
          DO I=1,N                                                              
!Exchange columns JK and K of matrix AA                                         
            TT=AA(I,JK)                                                         
            AA(I,JK)=AA(I,K)                                                    
            AA(I,K)=TT                                                          
          ENDDO                                                                 
        ENDIF                                                                   
!Pivot is at correct column and located in K,K                                              
                                                                               
!Store column K in vector CS then set column K to zero                                             
        DO I=1,N                                                                
          CS(I)=AA(I,K)                                                         
          AA(I,K)=CMPLX(0.0,0.0)                                                          
        ENDDO                                                                   
!                                                                               
        CS(K)=CMPLX(0.0,0.0)                                                                
        AA(K,K)=CMPLX(1.0,0.0)                                                              
!Modify line K :                                            
        IF(ABS(PV).LT.EPSMACH) THEN                                            
          WRITE(*,*) '  PIVOT TOO SMALL - STOP'                               
          STOP                                                                  
        ENDIF                                                                   
        DO I=1,N                                                                
          AA(K,I)=AA(K,I)/PV                                                    
        ENDDO                                                                   
        IF (M.NE.0) THEN                                                        
          DO I=1,M                                                             
            BB(K,I)=BB(K,I)/PV                                                  
          ENDDO                                                                 
        ENDIF                                                                   
!Modify other lines of matrix AA:                                        
        DO J=1,N                                                                
          IF (J.EQ.K) CONTINUE                                                  
          DO I=1,N                                                              
!Modify line J of matrix AA :                                            
            AA(J,I)=AA(J,I)-CS(J)*AA(K,I)                                       
          ENDDO                                                                 
          IF (M.NE.0) THEN                                                      
            DO I=1,M                                                          
              BB(J,I)=BB(J,I)-CS(J)*BB(K,I)                                     
            ENDDO                                                               
          ENDIF                                                                 
        ENDDO                                                                   
!Line K is ready.                                                
      ENDDO                                                                     
!End of K loop                                                              
                                                                               
!The matrix AA is inverted - Rearrange AA                         
                                                                               
!Exchange lines                                                            
      DO I=N,1,-1                                                               
        IK=PC(I)                                                                
        IF (IK.EQ.I) CONTINUE                                                   
!EXCHANGE LINES I AND PC(I) OF AA:                                         
        DO J=1,N                                                                
          TT=AA(I,J)                                                            
          AA(I,J)=AA(IK,J)                                                      
          AA(IK,J)=TT                                                           
        ENDDO                                                                   
        IF (M.NE.0) THEN                                                        
          DO J=1,M                                                             
            TT=BB(I,J)                                                          
            BB(I,J)=BB(IK,J)                                                    
            BB(IK,J)=TT                                                         
          ENDDO                                                                 
        ENDIF                                                                   
!NO MORE EXCHANGE NEEDED                                                      
!GO TO NEXT LINE                                                  
      ENDDO                                                                     
                                                                               
!EXCHANGE COLUMNS                                                          
                                                                               
      DO J=N,1,-1                                                               
        JK=PL(J)                                                                
        IF (JK.EQ.J) CONTINUE                                                   
!EXCHANGE COLUMNS J AND PL(J) OF AA :                                       
        DO I=1,N                                                                
          TT=AA(I,J)                                                            
          AA(I,J)=AA(I,JK)                                                      
          AA(I,JK)=TT                                                           
        ENDDO                                                                   
!NO MORE EXCHANGE NEEDED                                                      
!GO TO NEXT COLUMN   
      ENDDO                                                                     
!REARRANGEMENT TERMINATED.                                                        
      RETURN                                                                    
   10 FORMAT(///'  DETERMINANT EQUALS ZERO, NO SOLUTION!')                    
      END                                                                       

!*******************************************
!* MULTIPLICATION OF TWO COMPLEX MATRICES  *
!* --------------------------------------- *
!* INPUTS:    A1  COMPLEX MATRIX N*N       *
!*            B   COMPLEX MATRIX N*M       *
!*            N,M INTEGER                  *
!* --------------------------------------- *
!* OUTPUTS:   C   COMPLEX MATRIX N*M,      *
!*                PRODUCT A*B              *
!*                                         *
!*******************************************
      SUBROUTINE MATMUL(A,B,C,N,M)
      COMPLEX A(N,N),B(N,M),C(N,M),SUM                                           
      DO I=1,N                                                                  
        DO J=1,M                                                                
          SUM=CMPLX(0.0,0.0)                                                                
          DO K=1,N                                                              
            SUM=SUM+A(I,K)*B(K,J)                                               
          ENDDO                                                                 
          C(I,J)=SUM                                                            
        ENDDO                                                                   
      ENDDO                                                                     
      RETURN                                                                    
      END                                                                       
                                                                               
! End of file csysmat.f90  DIGITAL FORTRAN version, Feb. 2010