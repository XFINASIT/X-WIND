      SUBROUTINE INERTIARELIEF (MAXA,FORCE)
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /DYNA_STEP/ INC
      DIMENSION FORCE(1)
      ALLOCATABLE AK(:,:),AM(:,:),D(:,:),F(:,:)
      DIMENSION MAXA(NEQ)
      ALLOCATABLE Phi(:,:)
      ALLOCATABLE Phi_r(:,:),Phi_rt(:,:)
      ALLOCATABLE Factor_r1(:,:),Factor_r2(:,:),Factor_r(:,:)
      ALLOCATABLE Phi_et(:,:),Phi_e(:,:),Factor_e(:,:)
      ALLOCATABLE Factor_e1(:,:),Factor_e2(:,:)
      ALLOCATABLE diagFactor_r(:,:),diagFactor_e(:,:)
      ALLOCATABLE FE(:,:)
      ALLOCATABLE Omega_e(:,:), Omega_e1 (:,:),Omega_e_INV(:,:) ! OUTPUT RESULT FOR SOLVER
      ALLOCATABLE DISP(:,:),OMEGAINV_FE(:,:),Omega_e_UP(:)
      POINTER ::  AA(:),MAXC(:),IPERM(:),X(:)
      ALLOCATABLE FORCE_TEST(:,:),DISPLACEMENT_STORE(:)
      
      ! --------- DYNAMIC -------
      ALLOCATABLE  c1(:,:),c2(:,:),c3(:,:)
      ALLOCATABLE  U(:,:),U_1(:,:),U_dot(:,:),U_dot_1(:,:),U_2_dot(:,:),U_2_dot_1(:,:)
      ALLOCATABLE  eta_d(:,:),eta_v(:,:),eta_a(:,:)
      ALLOCATABLE  eta_d_1(:,:),eta_v_1(:,:),eta_a_1(:,:)
      ALLOCATABLE  delP(:,:),delu(:,:),delv(:,:),dela(:,:)
      ALLOCATABLE  Iden2(:,:)
      ALLOCATABLE  C(:,:)
      ALLOCATABLE  Kcap(:,:)
      ALLOCATABLE  invKcap(:,:)
      ALLOCATABLE  Fe_1(:,:)
      ! --------- DYNAMIC -------
      
      
      IF (ISOLOP.EQ.1) THEN
      NT = 1    ! LINEAR STATIC ANALYSIS
      ELSEIF (ISOLOP.EQ.5) THEN
      NT = NINC ! LINER DYNAMIC ANALYSIS
      ENDIF 
      
      ! TEST VALUE 
      DDT = 0.05D0
      NT = 200
      
      NR = 6 ! RIGID BODY MODE
      NDOF = NEQ
      NMODE = 20D0 
      NE = NMODE - NR
      
      ALLOCATE (AK(NDOF,NDOF),AM(NDOF,NDOF),D(NMODE,NMODE),F(NDOF,NT))
      ALLOCATE (Phi(NDOF,NMODE))
      ALLOCATE (Phi_r(NDOF,NMODE),Phi_rt(NR,NDOF))
      ALLOCATE (Factor_r1(NR,NDOF),Factor_r2(NR,NR),Factor_r(NR,1))
      ALLOCATE (Phi_et(NE,NDOF),Phi_e(NDOF,NE),Factor_e(NE,1))
      ALLOCATE (Factor_e1(NE,NDOF),Factor_e2(NE,NE))
      ALLOCATE (diagFactor_r(NR,NR),diagFactor_e(NE,NE))
      ALLOCATE (FE(NE,NT))
      ALLOCATE (Omega_e(NE,NE), Omega_e1(NE,NDOF),Omega_e_INV(NE,NE))
      ALLOCATE (DISP(NDOF,NT),OMEGAINV_FE(NE,NE),Omega_e_UP(NWK))
      ALLOCATE (FORCE_TEST(NDOF,200),DISPLACEMENT_STORE(200))
      
      ! --------- DYNAMIC -------
      ALLOCATE (c1(NE,NE),c2(NE,NE),c3(NE,NE))
      ALLOCATE (U(ndof,NT),U_1(ndof,NT),U_dot(ndof,NT),U_dot_1(ndof,NT),U_2_dot(ndof,NT),U_2_dot_1(ndof,NT))
      ALLOCATE (eta_d(NE,NT),eta_v(NE,NT),eta_a(NE,NT))
      ALLOCATE (eta_d_1(NE,1),eta_v_1(NE,1),eta_a_1(NE,1))
      ALLOCATE (delP(NE,1),delu(NE,1),delv(NE,1),dela(NE,1))
      ALLOCATE (Iden2 (NE,NE))
      ALLOCATE (C(NE,NE))
      ALLOCATE (Kcap(NE,NE))
      ALLOCATE (invKcap(NE,NE))
      ALLOCATE (Fe_1(NE,1))
      ! --------- DYNAMIC -------
      
      ! ----- SET INITIAL DATA -----
      AK           = 0.0D0
      AM           = 0.0D0
      D            = 0.0D0
      F            = 0.0D0
      Phi          = 0.0D0
      Phi_r        = 0.0D0
      Phi_rt       = 0.0D0
      Factor_r1    = 0.0D0
      Factor_r2    = 0.0D0
      Factor_r     = 0.0D0
      Phi_et       = 0.0D0
      Phi_e        = 0.0D0
      Factor_e     = 0.0D0
      Factor_e1    = 0.0D0
      Factor_e2    = 0.0D0
      diagFactor_r = 0.0D0
      diagFactor_e = 0.0D0
      Omega_e      = 0.0D0
      Omega_e1     = 0.0D0
      ! --------------------------
      
      open(unit=100,file="Inertia Testing File\FB-Dyn-EMM-Global-Stiffness.csv",status="unknown")
      do i=1, ndof
          read(100,*) (AK(i,j), j=1,ndof)
      end do
      CLOSE (100)

      open(unit=200,file="Inertia Testing File\FB-Dyn-EMM-Global-Mass.csv",status="unknown")
      do i=1,ndof      !define M
          read(200,*) (AM(i,j), j=1,ndof)
      end do
      CLOSE (200)
      
      open(unit=300,file="Inertia Testing File\FB-Dyn-EMM-External Force.csv",status="unknown")
      do i=1,ndof       !define F
          read(300,*) (FORCE_TEST(i,j), j=1,200)
      end do
      CLOSE (300)
      
      open(unit=400,file="Inertia Testing File\FB-Dyn-EMM-Eigenvectors.csv",status="unknown")
      do i=1,ndof       !define Phi
          read(400,*) (Phi(i,j), j=1,nmode)
      end do
      CLOSE (400)

      open(unit=500,file="Inertia Testing File\FB-Dyn-EMM-Eigenvalues.csv",status="unknown")
      do i=1,nmode       !define D
          read(500,*) (D(i,j), j=1,nmode)
      end do
      CLOSE (500)
      
      ! READ IGEN VALUE
      
      ! READ IGEN VECTOR
      
      ! READ STIFFNESS MATRIX
      CALL MDMOVE('STIF','TEMP')
      
      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

      ALLOCATE(AA(NWK),MAXC(NWK))
      AA = 0.
      MAXC = 0.
      
      CALL MCALFIL(KSREC,'TEMP')
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
      
      ! FULL STIFFNESS MATRIX
      CALL READSTIFF(NEQ,NWK,MAXC,MAXA,AA,AK) 
      
      ! READ MASS MATRIX
      CALL MCALFIL(KSREC,'MASS')
	REWIND(KSREC)
      AA   = 0.
      IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	!READ(NFLCH) MAXC(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
      ! FULL MASS MATRIX
      CALL READSTIFF(NEQ,NWK,MAXC,MAXA,AA,AM) 
      
      
      DO IIII = 1,200 
      ! READ FORCE MATIX
      IF (ISOLOP.EQ.1) THEN
      F(1:NDOF,1) = FORCE(1:NDOF)
      ENDIF

      F(1:NDOF,1) = FORCE_TEST(1:NDOF,IIII) 
      
      
      ! RIGID BODY MODE
      CALL Preprocesseigenvectors (NR,NDOF,NMODE,NE,Phi,AM,D,Phi_e,Phi_r,Phi_rt)
      
      CALL Preprocessforinertiareliefoperation (NR,NDOF,NMODE,NE,NT,AM,AK,Phi_e,f,Phi_r,Phi_rt,FE,Omega_e)
      
      ! REAPLACE UPPER TRIANGLE MATRIX AND SOLVE FOR DISPLACMENT
      Omega_e_UP = 0.D0
      INDEX_STA = 1D0
      INDEX_END = 1D0
      DO I = 1,NE
         IF (I.EQ.1) INDEX_END = NE 
         IF (I.GT.1) INDEX_END = INDEX_END + NE - I + 1
         Omega_e_UP(INDEX_STA:INDEX_END) = Omega_e(I,I:NE)
         INDEX_STA = INDEX_END + 1
      ENDDO
      NWK_NEW = INDEX_STA - 1 ! TOTAL NUMBER OF DATA 
      
      !INVERSE MATRIX
      Omega_e_INV = 0.
      OMEGAINV_FE = 0.
      DISP = 0.
      CALL INVMAT(Omega_e,Omega_e_INV,NE)
      OMEGAINV_FE = MATMUL(Omega_e_INV,FE) 
      DISP = MATMUL(Phi_e,OMEGAINV_FE)
      
      ! UPDATE A(LDL) FOR DISPLACEMENT
      IF (ISOLOP.EQ.1)  FORCE(1:NEQ) = DISP(1:NEQ,1)
      
      ! FOR EXAMPLE TEST CASE 1
      DISPLACEMENT_STORE(IIII) = DISP(3,1)
      WRITE (299,100) DISPLACEMENT_STORE(IIII)
100   FORMAT (E12.5)
      
      ! DYNAMIC TIME INTEGRATION
            dt = DDT
            
            !Constants used in Newmark's integration
            gaama=0.5 
            beta=0.25  !Values for unconditional stability of method

            !Define Identity matrix
            do i=1,Ne    
                do j=1,Ne
                    if(i==j) then
                        Iden2(i,j)=1
                    else
                        Iden2(i,j)=0
                    end if
                end do
            end do

            !Define C
            C = 0.

            a1=gaama/(beta*dt)    
            a2=1/(beta*dt*dt) 
            a3=1/(beta*dt)        
            a4=gaama/beta-1 
            a5=1/(2*beta)-1      
            a6=(gaama/(2*beta)-1)*dt 

            c1=a2*Iden2+a1*C 
            c2=a3*Iden2+a4*C 
            c3=a5*Iden2+a6*C 

            !Initialization of generalized coordinates for displacement,velocity and acceleration
            !define eta_d
            eta_d = 0
        
            !define eta_v
            eta_v = 0.

            !define eta_a
            eta_a = 0.

            !Initialization of actual disp,vel and acc
            U = 0.

            U_dot=U

            U_2_dot=U

            !Initial Conditions
            eta_d = 0.
            eta_v = 0.
      CALL TIMEDOMAIN (NR,NDOF,NMODE,NE,NT,Phi_e,Omega_e,FE,eta_v,eta_d)
      
      ENDDO
      
      
      DEALLOCATE (AK,AM,D,F)
      DEALLOCATE (Phi)
      DEALLOCATE (Phi_r,Phi_rt)
      DEALLOCATE (Factor_r1,Factor_r2,Factor_r)
      DEALLOCATE (Phi_et,Phi_e,Factor_e)
      DEALLOCATE (Factor_e1,Factor_e2)
      DEALLOCATE (diagFactor_r,diagFactor_e)
      DEALLOCATE (Omega_e,Omega_e1,Omega_e_INV) 
      DEALLOCATE (DISP,OMEGAINV_FE,Omega_e_UP)
      
      DEALLOCATE (AA,MAXC)
      DEALLOCATE (FE)
      DEALLOCATE (FORCE_TEST,DISPLACEMENT_STORE)
      
      
      
      ! --------------
      DEALLOCATE (c1,c2,c3)
      DEALLOCATE (U,U_1,U_dot,U_dot_1,U_2_dot,U_2_dot_1)
      DEALLOCATE (eta_d,eta_v,eta_a)
      DEALLOCATE (eta_d_1,eta_v_1,eta_a_1)
      DEALLOCATE (delP,delu,delv,dela)
      DEALLOCATE (Iden2)
      DEALLOCATE (C)
      DEALLOCATE (Kcap)
      DEALLOCATE (invKcap)
      DEALLOCATE (Fe_1)
       ! --------------
      
      END
C    ===================================================================================
      SUBROUTINE Preprocesseigenvectors (NR,NDOF,NMODE,NE,Phi,AM,D,Phi_e,Phi_r,Phi_rt)
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N) 
      DIMENSION Phi(NDOF,NMODE),AM(NDOF,NDOF),D(NMODE,NMODE) 
      DIMENSION Phi_r(NDOF,NMODE),Phi_rt(NR,NDOF)
      DIMENSION Factor_r1(NR,NDOF),Factor_r2(NR,NR),Factor_r(NR,1)
      DIMENSION Phi_et(NE,NDOF),Phi_e(NDOF,NE),Factor_e(NE,1)
      DIMENSION Factor_e1(NE,NDOF),Factor_e2(NE,NE)
      DIMENSION diagFactor_r(NR,NR),diagFactor_e(NE,NE)
      ALLOCATABLE W(:,:)
      
      ALLOCATE (W(NMODE,2))
            !preprocess eigen vectors. Steps 3 & 4
            !Separate rigid and elastic modes from each other
           do i=1,ndof  !define Phi_r
                do j=1,Nr
                    Phi_r(i,j)=Phi(i,j)
                end do
            end do

            do i=1,ndof   !define Phi_e
                do j=1,nmode-Nr
                    Phi_e(i,j)=Phi(i,j+Nr)
                end do
            end do

            !Sort eigen vectors
            do i = 1, nmode  !sorting and diagonalizing D
                W(i,1) = D(i,i) 
            end do
        
            do j=1,nmode
                do i=j+1,nmode
                    if (W(j,1) < W(i,1)) then  
                        rmax=i  
                    end if
                end do
                if (rmax/=j) then   
                    temp=W(rmax,1)  
                    W(rmax,1)=W(j,1) 
                    W(j,1)=temp    
                end if
            end do

            do i = 1,nmode
                W(i,2)=i
            end do

            !Normalize the eigen vectors
            !Procedures for calculation of Factor=diag(Phi_rt*M*Phi_r)
            Phi_rt(1:Nr,1:ndof)=transpose(Phi_r(1:ndof,1:Nr))   !Phi_rt is transpose of Phi_r
            Factor_r1(1:Nr,1:ndof)=matmul(Phi_rt(1:Nr,1:ndof),AM(1:ndof,1:ndof))      !multiply Phi_rt by M
            Factor_r2(1:Nr,1:Nr)=matmul(Factor_r1(1:Nr,1:ndof),Phi_r(1:ndof,1:Nr))    !finish multiplication
            do i=1,Nr     !diagonalize Factor_e2
                Factor_r(i,1)=Factor_r2(i,i)
            end do

            !Procedures for calculation of Phi_r=Phi_r*inv(sqrt(diag(Factor)))
            do i=1,Nr       !inv(sqrt(diag(Factor)))
                do j=1,Nr     
                    if (i==j)  then
                        diagFactor_r(i,j)=1/sqrt(Factor_r(i,1))
                    else
                        diagFactor_r(i,j)=0
                    end if
                end do
            end do 
            Phi_r=matmul(Phi_r,diagFactor_r)
            
            Phi_et(1:Ne,1:ndof)=transpose(Phi_e(1:ndof,1:Ne))   !Phi_et is transpose of Phi_e
            Factor_e1(1:Ne,1:ndof)=matmul(Phi_et(1:Ne,1:ndof),AM(1:ndof,1:ndof))      !multiply Phi_et by M
            Factor_e2(1:Ne,1:Ne)=matmul(Factor_e1(1:Ne,1:ndof),Phi_e(1:ndof,1:Ne))    !finish multiplication
            do i=1,Ne     !diagonalize Factor_e2
                Factor_e(i,1)=Factor_e2(i,i)
            end do

            !Procedures for calculation of Phi_e=Phi_e*inv(sqrt(diag(Factor)))
            do i=1,Ne       !inv(sqrt(diag(Factor)))
                do j=1,Ne     
                    if (i==j)  then
                        diagFactor_e(i,j)=1/sqrt(Factor_e(i,1))
                    else
                        diagFactor_e(i,j)=0
                    end if
                end do
            end do 
        
            Phi_e=matmul(Phi_e,diagFactor_e)
      
      DEALLOCATE (W)
      
      end subroutine Preprocesseigenvectors
      
C    ===================================================================================
      SUBROUTINE Preprocessforinertiareliefoperation (NR,NDOF,NMODE,NE,NT,AM,AK,Phi_e,f,Phi_r,Phi_rt,FE,Omega_e)
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N) 
      ALLOCATABLE Iden(:,:)
      ALLOCATABLE A_1(:,:),A_2(:,:)
      ALLOCATABLE A(:,:)
      ALLOCATABLE Phi_et(:,:)
      DIMENSION AK(NDOF,NDOF),AM(NDOF,NDOF),F(NDOF,NT)
      DIMENSION Phi_r(NDOF,NMODE),Phi_rt(NR,NDOF)
      DIMENSION Phi_e(NDOF,NE)
      DIMENSION FE(NE,NT)
      DIMENSION Omega_e(NE,NE),Omega_e1(NE,NDOF)
      
      
      ALLOCATE (Iden(NDOF,NDOF))
      ALLOCATE (A_1(NDOF,NR),A_2(NDOF,NDOF))
      ALLOCATE (A(NDOF,NDOF))
      ALLOCATE (Phi_et(NE,NDOF))
      
      Iden = 0.
      A    = 0.
      A_1  = 0.
      A_2  = 0.
      FE   = 0.
      Phi_et =0.
      
        !Preprocess for inertia relief operation

            !! Elastic calculation
            !Calculation of Inertia Relief Profection Matrix
            !Defining Identity matrix
            do i=1,ndof    
                do j=1,ndof
                    if(i==j) then
                        Iden(i,j)=1
                    else
                        Iden(i,j)=0
                    end if
                end do
            end do

            !Procedures for calculation of A=Iden-M*Phi_r*Phi_rt
            A_1(1:ndof,1:Nr)=matmul(AM(1:ndof,1:ndof),Phi_r(1:ndof,1:Nr))       
            A_2(1:ndof,1:ndof)=matmul(A_1(1:ndof,1:Nr),Phi_rt(1:Nr,1:ndof))    
            A(1:ndof,1:ndof)=Iden(1:ndof,1:ndof)-A_2(1:ndof,1:ndof)

            !Calculation of Elastic modal forces
            Fe(1:Ne,1:nt)=matmul(transpose(Phi_e),F(1:ndof,1:nt))

            !Calculations of elastic modal frequency factors
            Phi_et=transpose(Phi_e(1:ndof,1:Ne))
            Omega_e1=matmul(Phi_et(1:Ne,1:ndof),AK(1:ndof,1:ndof))
            Omega_e=matmul(Omega_e1(1:Ne,1:ndof),Phi_e(1:ndof,1:Ne))
            
            ! CALL SOLVER
            
      DEALLOCATE (Iden)
      DEALLOCATE (A_1,A_2)
      DEALLOCATE (A)
      DEALLOCATE (Phi_et)

      END subroutine Preprocessforinertiareliefoperation
C    ===================================================================================      
      SUBROUTINE TIMEDOMAIN (NR,NDOF,NMODE,NE,NT,Phi_e,Omega_e,FE,eta_v,eta_d)
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)
      !ALLOCATABLE  c1(:,:),c2(:,:),c3(:,:)
      !ALLOCATABLE  U(:,:),U_1(:,:),U_dot(:,:),U_dot_1(:,:),U_2_dot(:,:),U_2_dot_1(:,:)
      !ALLOCATABLE  eta_d(:,:),eta_v(:,:),eta_a(:,:)
      !ALLOCATABLE  eta_d_1(:,:),eta_v_1(:,:),eta_a_1(:,:)
      !ALLOCATABLE  delP(:,:),delu(:,:),delv(:,:),dela(:,:)
      !ALLOCATABLE  Iden2(:,:)
      !ALLOCATABLE  C(:,:)
      !ALLOCATABLE  Kcap(:,:)
      !ALLOCATABLE  invKcap(:,:)
      !ALLOCATABLE  Fe_1(:,:)
      DIMENSION    Omega_e(NE,NE)
      DIMENSION    FE(NE,NT)
      DIMENSION    Phi_e(NDOF,NE)
      
      DIMENSION c1(NE,NE),c2(NE,NE),c3(NE,NE)
      DIMENSION U(ndof,NT),U_1(ndof,NT),U_dot(ndof,NT),U_dot_1(ndof,NT),U_2_dot(ndof,NT),U_2_dot_1(ndof,NT)
      DIMENSION eta_d(NE,NT),eta_v(NE,NT),eta_a(NE,NT)
      DIMENSION eta_d_1(NE,1),eta_v_1(NE,1),eta_a_1(NE,1)
      DIMENSION delP(NE,1),delu(NE,1),delv(NE,1),dela(NE,1)
      DIMENSION Iden2 (NE,NE)
      DIMENSION C(NE,NE)
      DIMENSION Kcap(NE,NE)
      DIMENSION invKcap(NE,NE)
      DIMENSION Fe_1(NE,1)
      
      
      !ALLOCATE (c1(NE,NE),c2(NE,NE),c3(NE,NE))
      !ALLOCATE (U(ndof,NT),U_1(ndof,NT),U_dot(ndof,NT),U_dot_1(ndof,NT),U_2_dot(ndof,NT),U_2_dot_1(ndof,NT))
      !ALLOCATE (eta_d(NE,NT),eta_v(NE,NT),eta_a(NE,NT))
      !ALLOCATE (eta_d_1(NE,1),eta_v_1(NE,1),eta_a_1(NE,1))
      !ALLOCATE (delP(NE,1),delu(NE,1),delv(NE,1),dela(NE,1))
      !ALLOCATE (Iden2 (NE,NE))
      !ALLOCATE (C(NE,NE))
      !ALLOCATE (Kcap(NE,NE))
      !ALLOCATE (invKcap(NE,NE))
      !ALLOCATE (Fe_1(NE,1))
      
            !!!eta_a(:,1) = II\(Fe(:,1)-C*eta_v(:,1)-Omega_e*eta_d(:,1)) ; !!!!!!
            do i=1,Ne
                Fe_1(i,1)=Fe(i,1)
            end do
            
            do i=1,Ne
                eta_v_1(i,1)=eta_v(i,1)
            end do 
            
            do i=1,Ne
                eta_d_1(i,1)=eta_d(i,1)
            end do 
            
            eta_a_1=Fe_1-matmul(C,eta_v_1)-matmul(Omega_e,eta_d_1)
            
            do i=1,Ne
                eta_a(i,1)=eta_a_1(i,1)
            end do
                    
            U_1=matmul(Phi_e(1:ndof,1:Ne),eta_d(1:Ne,1:nt))
            do i=1,ndof
                U(i,1)=U_1(i,1)
            end do

            U_dot_1=matmul(Phi_e(1:ndof,1:Ne),eta_v(1:Ne,1:nt))
            do i=1,ndof
                U_dot(i,1)=U_dot_1(i,1)
            end do

            U_2_dot_1=matmul(Phi_e(1:ndof,1:Ne),eta_a(1:Ne,1:nt))
            do i=1,ndof
                U_2_dot(i,1)=U_2_dot_1(i,1)
            end do

            Kcap=omega_e(1:Ne,1:Ne)+a1*C+a2*Iden2

            !Time step starts
            j=1

30          do i=1,Ne
                    eta_d_1(i,1)=eta_d(i,j)
            end do
             
            do i=1,Ne
                eta_v_1(i,1)=eta_v(i,j)
            end do 
            
            do i=1,Ne
                eta_a_1(i,1)=eta_a(i,j)
            end do 
            
            do i=1,Ne
                Fe_1(i,1)=Fe(i,j+1)
                !Fe_1(i,1)=Fe(i,j)
            end do  

            delP=matmul(c1,eta_d_1)+matmul(c2,eta_v_1)+matmul(c3,eta_a_1)+Fe_1
!
!            call inverse(Kcap,invKcap,Ne)
            CALL INVMAT(Kcap,invKcap,Ne)

            delu=matmul(invKcap,delP)  

            delv=a1*(delu-eta_d_1)-a4*eta_v_1-a6*eta_a_1

            do i=1,Ne
                dela(i,1) = a2*(delu(i,1)-eta_d(i,j))-a3*eta_v(i,j)-a5*eta_a(i,j);
            end do

            do i=1,Ne
                eta_d(i,j+1)=delu(i,1)
            end do

            do i=1,Ne
                eta_v(i,j+1)=delv(i,1)
            end do

            do i=1,Ne
                eta_a(i,j+1)=dela(i,1)
            end do
 
            !Recovery of disp,vel and acc in physical domain
            U_1=matmul(Phi_e(1:ndof,1:Ne),eta_d(1:Ne,1:nt))
            do i=1,ndof
                U(i,j+1)=U_1(i,j+1)
            end do

            U_dot_1=matmul(Phi_e(1:ndof,1:Ne),eta_v(1:Ne,1:nt))
            do i=1,ndof
                U_dot(i,j+1)=U_dot_1(i,j+1)
            end do

            U_2_dot_1=matmul(Phi_e(1:ndof,1:Ne),eta_a(1:Ne,1:nt))
            do i=1,ndof
                U_2_dot(i,j+1)=U_2_dot_1(i,j+1)
            end do

!            j=j+1
!            if (j<nt) then
!                goto 30
!            end if
            
      
      !DEALLOCATE (c1,c2,c3)
      !DEALLOCATE (U,U_1,U_dot,U_dot_1,U_2_dot,U_2_dot_1)
      !DEALLOCATE (eta_d,eta_v,eta_a)
      !DEALLOCATE (eta_d_1,eta_v_1,eta_a_1)
      !DEALLOCATE (delP,delu,delv,dela)
      !DEALLOCATE (Iden2)
      !DEALLOCATE (C)
      !DEALLOCATE (Kcap)
      !DEALLOCATE (invKcap)
      !DEALLOCATE (Fe_1)
      END