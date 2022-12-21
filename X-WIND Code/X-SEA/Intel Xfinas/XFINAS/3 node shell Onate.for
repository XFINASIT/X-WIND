
C     =====================================================================
C     ============ 3 NODE SHELL WITH NO ROTATION DOF BY ONATE =============
C     =====================================================================
      SUBROUTINE ONATESHEL3(PROPM,PROPG,WA,S,COORD,EDIS,RE,MWG,IFFIX,ISIDE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
     
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C     ----------------------------------------------------------------

      DIMENSION PROPM(1),PROPG(1),WA(MWG,1),COORD(3,1)
      DIMENSION S(18,18),EDIS(18),RE(18)
      
      DIMENSION DR(6,6)
      DIMENSION EPS(6),EPSQ(6),SIGR(6)
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION COVR(3),COVS(3)
      
      DIMENSION H(6),HD(2,6),COORDI(3,6)
	DIMENSION BA(6,36),REDIS(36)

      DIMENSION IFFIX(3),ISIDE(3)
      
      DIMENSION T(3,9),T3N(3),A(3,4),B(3,4),AREA2I(4)
      
      DIMENSION STRN(6),STRS(6),BB(6,18),BBV(6,18),BBG(10,18)
      DIMENSION BMEM(3,6,3),BMEMNL(3,6,4),STRNM(3)
      DIMENSION BBEN(3,6,3),BBENNL(3,6,6),STRNB(3)
      DIMENSION DM(10,10)
C     CO-ROTATIONAL DISPLACEMENT      
      DIMENSION CODIS(18)
      
      
C     SHELL THICKNESS      
      TH = PROPG(2)
      
C     INITIAL COORDS (COORDI) 
      CALL COINI(COORDI,COORD,EDIS,NLOPT)

C     -------------------------------------------------
C     REMOVE RIGID BODY TRANSLATIONS AND ROTATIONS FROM TOTAL DISPLACEMENT VECTOR (NLOPT>1)
      REDIS(1:18) = EDIS(1:18)
      IF (NLOPT.LE.1) GOTO 200
      CALL S3MDSH(COORD,COORDI,REDIS)
C     MAKE ZERO ON DISPLACEMENT OF NON-EXIST ADJACENT NODES      
      DO I = 1,3
          K = ISIDE(I)
          IF(K.EQ.0) THEN
              J = I + 3
              REDIS(3*J-2) = 0.0D0
              REDIS(3*J-1) = 0.0D0
              REDIS(3*J-0) = 0.0D0
          ENDIF
      ENDDO      
200   CONTINUE
C     ------------------------------------------------- 
      
      
C     LOCAL BASE VECTOR -- AT ORIGINAL CONFIG
	CALL SHAPD3(COORDI,A,B,AREA2I,ISIDE)

C     DEFORM CONFIG
	CALL ARETTR(COORD,T3N,AREA2N)
      AREA = 0.5*AREA2I(1)

C     RATIO BETWEEN ORIGINAL AND DEFORMED AREA	
	SLAM = AREA2I(1)/AREA2N
	

C     ---------------------
C     MASS MATRIX (ITASK=5)
	IF (ITASK.NE.5) GOTO 100

      GOTO 900
C     ---------------------      
100   CONTINUE
C     ---------------------
           
	
C	CARTESIAN DERIVATIVES AND ASSIGN THEM IN TO T MATRIX
	CALL CARTD3(A,B,T,COORD,ISIDE,IFFIX,SLAM)	
	
C	CALCULATE B MATRIX
      BMEM(1:3,1:6,1:3) = 0.0D0
	CALL BMEM3(A,B,T,BMEM) !MEMBRANE STRAIN AND B MATRIX
 
C     BENDING PART
      FACCMP = 1.0D0
      IF(NLOPT.EQ.3) FACCMP = 10.0D0 !SCALE FACTOR FOR CLAMP EFFECT (INCREASE EFFECT FOR LARGE DEFLECTION ANALYSIS)
      BBEN(1:3,1:6,1:3) = 0.0D0
	CALL BBEND3(COORDI,A,B,T,ISIDE,IFFIX,SLAM,FACCMP,BBEN)  
	 
C     STRAIN DISP. MATRIX BB
	CALL SHSDB3(BB,BMEM,BBEN)

C     BENDING PART
      FACCMP = 1.0D0
      BBEN(1:3,1:6,1:3) = 0.0D0
	CALL BBEND3(COORDI,A,B,T,ISIDE,IFFIX,SLAM,FACCMP,BBEN)  
	 
C     STRAIN DISP. MATRIX BB
	CALL SHSDB3(BBV,BMEM,BBEN)
		
C     CALCULATE STRAIN
      CALL CALSTRN3(BBV,REDIS,STRN)
      STRNM(1:3) = STRN(1:3) !MEMBRANE STRAIN
      STRNB(1:3) = STRN(4:6) !BENDING  STRAIN

C	CALCULATE RIGIDITY MATRIX FOR ELASTIC ANALYSIS
	CALL SHDR3(DR,PROPM,TH)

C     STRESS RESULTANT
      STRS = MATMUL(DR,STRN)


C     ------------------------------------------------------
C     FOR BUCKLINGANALYSIS -- GET GEOMETRIX STIFFNESS	
      IF (ISOLOP.EQ.4.AND.IFEIG.EQ.0) GOTO 800
C     ------------------------------------------------------


C     --------------------------
C     INTERNAL FORCES BLOG
C     --------------------------  
500   CONTINUE   
C     GET INTERNAL FORCES
      RE = AREA*MATMUL(TRANSPOSE(BBV),STRS)  
C     --------------------------   


C     ------------------------------------------------------
C     FLAG TO GET STIFFNESS MATRIX (IFREF=0)
      IF (IFREF ) 900,700,900
C     ------------------------------------------------------
          
          
C     --------------------------
C     LINEAR STIFFNESS BLOG
C     --------------------------  
700   CONTINUE
C     GET LINEAR STIFFNESS
      CALL STIFL3(BB,DR,AREA,S)
C     --------------------------   


C     --------------------------
C     GEOMETRIC STIFFNESS BLOG
C     --------------------------
      IF (NLOPT.LE.1) GOTO 900
800   CONTINUE
C     STRESS RIGIDITY
      CALL RIGD10(STRS,DM)
C     NONLINEAR B MATRIX
      CALL BNONL3(A,B,T,BBEN,BMEMNL,BBENNL)
      CALL SHSDG10(BBG,BMEMNL,BBENNL)
C     GET GEOMETRIC STIFFNESS
      CALL STIFG3(BBG,DM,AREA,S) 
C     --------------------------
      
900   CONTINUE

           
      RETURN
      END
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE SHSDG10(BBG,BMEM,BBEN)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4(I-N)
C	ASSEMBLING MEMBRAIN EFFECT AND BENDING EFFECT TO THE B MATRIX
	DIMENSION BMEM(18,4),BBEN(18,6),BBG(10,18)

C     MEMBRANE	
	DO I=1,4
	    DO J=1,18
	        BBG(I,J) = BMEM(J,I)
	    END DO
	END DO

C     BENDING	
	DO I=1,6
	    DO J=1,18
	        BBG(I+4,J) = BBEN(J,I)
	    END DO
	END DO


	RETURN
	END 
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE BNONL3(A,B,T,BBEN,BMEMNL,BBENNL)
	IMPLICIT REAL*8(A-H,O-Z) 
	IMPLICIT INTEGER*4(I-N) 
C	------------------------------------------------------------------
C     	COMPUTES BENDING B MATRIX FOR NONLINEAR
C	------------------------------------------------------------------
C	 ISIDE(3) 	- VALUE FOR SIDE ELEMENT EXISTS
C	 IFFIX(3) 	- 1 FOR CLAMPED SIDE
C	 A(3,4),B(3,4) 	- ELEMENT SIDE PROJECTIONS OVER LOCAL DIRECTIONS PER ELEMENT IN A PATCH
C	 T(3,9)       	 - DERIVATIVES OF CONFIGURATION 
C	 NBBEN3(3,6,3) 	- BENDING B MATRIX
C	------------------------------------------------------------------
	DIMENSION BBEN(3,6,3),BMEMNL(3,6,4),BBENNL(3,6,6),A(3,1),B(3,1),T(3,3)
      
C     Initialize BBENNL Matrix
      BMEMNL = 0.0
      BBENNL = 0.0
      
C     MEMBRANE B MATRIX FOR NONLINEAR	
	DO J=1,3
	    BMEMNL(1:3,J,1) = -B(J,1)*T(1:3,1)
	    BMEMNL(1:3,J,2) =  A(J,1)*T(1:3,1)
	    BMEMNL(1:3,J,3) = -B(J,1)*T(1:3,2)
	    BMEMNL(1:3,J,4) =  A(J,1)*T(1:3,2)
	ENDDO

C     BENDING B MATRIX FOR NONLINEAR           
      BBENNL(1:3,1:6,1) = BBEN(1:3,1:6,1)
      BBENNL(1:3,1:6,2) = BBEN(1:3,1:6,3)
      BBENNL(1:3,1:6,3) = BBEN(1:3,1:6,3)
      BBENNL(1:3,1:6,4) = BBEN(1:3,1:6,2)
      
      DO I = 1,3
          BBENNL(1:3,I,5) = BBENNL(1:3,I,5) - B(I,1)*T(1:3,3) 
          BBENNL(1:3,I,6) = BBENNL(1:3,I,6) + A(I,1)*T(1:3,3)
      ENDDO
      

      RETURN
      END
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE RIGD10(STRS,DM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------------------------------
C     Solve Componens of Nonlinear Rigidity. RM,RN and Solve Rigidity Matrix for Stiffness      
C     -------------------------------------------------------------------------------------	
c     DM 1,2,3 : Rigidity Matrix
C     RM : Generalized Stresses (Rigidity Moment)
c     RN : Generalized Stresses (Rigidity Force)
C     -------------------------------------------------------------------------------------	
	
      DIMENSION  STRS(6),DM1(10,10),DM2(10,10),DM3(10,10),DM(10,10)     
      
C     INTIALIZE    
      DM1(1:10,1:10) = 0.0D0
      DM2(1:10,1:10) = 0.0D0
      DM3(1:10,1:10) = 0.0D0
    
C     NORMAL FORCES      
      XN11 = STRS(1)
      XN22 = STRS(2)
      XN12 = STRS(3)
C     BENDING MOMENT      
      XM11 = STRS(4)
      XM22 = STRS(5)
      XM12 = STRS(6)
      
c     Calculate Stress Rigidity
     
      !==== DM1 ====!
      
      DM1(1,1) =  XN11      
      DM1(3,3) =  XN11
     
      DM1(1,5) = -XM11
      DM1(3,7) = -XM11

      DM1(5,1) = -XM11
      DM1(7,3) = -XM11
      
      DM1(9,9) =  XN11 

      
      !==== DM2 ====!
      
     
      DM2(2,2) =  XN22
      DM2(4,4) =  XN22
      
      DM2(2,6) = -XM22
      DM2(4,8) = -XM22

      DM2(6,2) = -XM22
      DM2(8,4) = -XM22
      
      DM2(10,10) = XN22 
      
      
      !==== DM3 ====!
      
      DM3(1,2) =  XN12      
      DM3(2,1) =  XN12
      DM3(3,4) =  XN12
      DM3(4,3) =  XN12
     
      DM3(1,6) = -XM12
      DM3(2,5) = -XM12
      DM3(3,8) = -XM12
      DM3(4,7) = -XM12

      DM3(5,2) = -XM12
      DM3(6,1) = -XM12
      DM3(7,4) = -XM12
      DM3(8,3) = -XM12
      
      DM3(9,10) = XN12   
      DM3(10,9) = XN12
    
      DM = DM1 + DM2 + DM3
      
      RETURN
      END
C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE STIFL3(BB,DR,AREA,S) 
      IMPLICIT REAL*8 (A-H,O-Z)
     	IMPLICIT INTEGER*4 (I-N)
C    -----------------------------------------------------------------
C     CALCULATES ELASTIC RIGIDITIES FOR MEMBRANE-PLATE OR SHELL ELEM.
C    -----------------------------------------------------------------
C     DR(6,6) = ELASTIC MEMBRANE,BENDING AND SHEAR RIGIDITIES
C    -----------------------------------------------------------------
      DIMENSION BB(6,18),DR(6,6),SS(18,18),S(18,18)
      
      SS = MATMUL(TRANSPOSE(BB),MATMUL(DR,BB))*AREA
      
      S = S + SS

 	RETURN
    	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE STIFG3(BB,DR,AREA,S) 
      IMPLICIT REAL*8 (A-H,O-Z)
     	IMPLICIT INTEGER*4 (I-N)
C    -----------------------------------------------------------------
C     CALCULATES ELASTIC RIGIDITIES FOR MEMBRANE-PLATE OR SHELL ELEM.
C    -----------------------------------------------------------------
C     DR(6,6) = ELASTIC MEMBRANE,BENDING AND SHEAR RIGIDITIES
C    -----------------------------------------------------------------
      DIMENSION BB(10,18),DR(10,10),SS(18,18),S(18,18)
      
      SS = MATMUL(TRANSPOSE(BB),MATMUL(DR,BB))*AREA
      
      S = S + SS

 	RETURN
    	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE SHDR3(DR,PROPM,TH)
      IMPLICIT REAL*8 (A-H,O-Z)
     	IMPLICIT INTEGER*4 (I-N)
C    -----------------------------------------------------------------
C     CALCULATES ELASTIC RIGIDITIES FOR MEMBRANE-PLATE OR SHELL ELEM.
C    -----------------------------------------------------------------
C     DR(6,6) = ELASTIC MEMBRANE,BENDING AND SHEAR RIGIDITIES
C    -----------------------------------------------------------------
      DIMENSION DR(36),PROPM(1)
      
      E = PROPM(1)
      v = PROPM(2)
      
      A1 = E/(1.0-(v*v))
      B1 = E*v/(1.0-(v*v))
      C1 = E/(2.0*(1.0+v))
	    
      CALL CLEARA (DR,36)
      TH3  = TH*TH*TH/12.0
      	DR(1)  = TH*A1
      	DR(8)  = TH*A1
     		DR(2)  = TH*B1
     		DR(7)  = TH*B1
      	DR(15) = TH*C1
      	DR(22) = TH3*A1
      	DR(29) = TH3*A1
     		DR(23) = TH3*B1
      	DR(28) = TH3*B1
      	DR(36) = TH3*C1
     
 	RETURN
    	END
C     =====================================================================
C     =====================================================================
C     =====================================================================

	SUBROUTINE SHSDB3 (BB,BMEM,BBEN)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4(I-N)
C	ASSEMBLING MEMBRAIN EFFECT AND BENDING EFFECT TO THE B MATRIX
	DIMENSION BMEM(18,3),BBEN(18,3),BB(6,18)

	DO I=1,3
	    DO J=1,18
	      BB(I,J)=BMEM(J,I)
	    END DO
	END DO

	DO I=1,3
	    DO J=1,18
	      BB(I+3,J)=BBEN(J,I)
	    END DO
	END DO

	RETURN
	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE CALSTRN3(BB,EDIS,STRN)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4(I-N)
C	CALCULATE MEMBRANE AND BENDING STRAIN
	DIMENSION BB(6,18),STRN(6),EDIS(18)

      DO I = 1,6
          STRN(I) = 0.0D0
          DO J = 1,18
            STRN(I) = STRN(I) + BB(I,J)*EDIS(J)
          ENDDO
      ENDDO

	RETURN
	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE BBEND3(COORD,A,B,T,ISIDE,IFFIX,SLAM,FACCMP,BBEN)
	IMPLICIT REAL*8(A-H,O-Z) 
	IMPLICIT INTEGER*4(I-N) 
C     COMPUTES BENDING B MATRIX FOR ELEMENT 
C     -----------------------------------------------------------------
C	 ISIDE(3) 	    - VALUE FOR SIDE ELEMENT EXISTS
C	 IFFIX(3) 	    - 1 FOR CLAMPED SIDE
C	 A(3,4),B(3,4) 	- ELEMENT SIDE PROJECTIONS OVER LOCAL DIIReCTIONS PER ELEMENT IN A PATCH
C	 T(3,9)           - DERIVATIVES OF CONFIGURATION 
C	 BBEN(3,6,3) 	    - BENDING B MATRIX
C     -----------------------------------------------------------------
	DIMENSION A(3,4),B(3,4),T(3,9),BBEN(3,6,3),SBBEN(3,6,3)
	DIMENSION K(3,3),NSS(2),D1(3),D2(3),XB(3)
	DIMENSION DN(3),BMN(3,3),BNO(2,2),BND(2,2)
	DIMENSION AN(3),TI(3,2),COORD(3,6)
	DIMENSION ISIDE(3),IFFIX(3)

	DATA K / 4,3,2, 5,1,3, 6,2,1 /    !SIDE ELEMENT CONNECTIVITIES

C     INITIALIZE 
      NSS   = 0
      SBBEN = 0.0
      D1    = 0.0
      D2    = 0.0
      XB    = 0.0
      DN    = 0.0
      BMN   = 0.0
      BNO   = 0.0
      BND   = 0.0
      AN    = 0.0
      TI    = 0.0
      
C	CONTRAVARIANT BASE VECTORS OF THE MAIN ELEMENT
      D1(1) = SLAM*(T(2,2)*T(3,3) - T(3,2)*T(2,3))
	D1(2) = SLAM*(T(3,2)*T(1,3) - T(1,2)*T(3,3))
      D1(3) = SLAM*(T(1,2)*T(2,3) - T(2,2)*T(1,3))
	D2(1) = SLAM*(T(2,3)*T(3,1) - T(3,3)*T(2,1))
	D2(2) = SLAM*(T(3,3)*T(1,1) - T(1,3)*T(3,1))
	D2(3) = SLAM*(T(1,3)*T(2,1) - T(2,3)*T(1,1))

C	-------------------------------------------------------------------------- 
C	-------------------------------------------------------------------------- 
	N = 4                             ! POSITION OF SIDE DERIVATIVE 1
	NS = 0                            ! INITIALIZES NUMBER OF S.S. OR FREE SIDES
	DO 100 I = 1,3 ! FOR EACH SIDE
          II = I+1                      ! POSITION OF SIDE IN ARRAYS
          
	    IF (ISIDE(I).GT.0) THEN ! SIDE ELEMENT EXISTS

C	CONTRIBUTIONS FROM NORMAL DISPLACEMENTS OF THE ADJ ELEMENTS
	        DO J=1,3 ! FOR EACH NODE
                  KK = K(J,I)                   ! LOCAL NODE (ADJACENT ELEMENT)
                  BBEN(1:3,KK,1)= BBEN(1:3,KK,1) +  B(I,1)*B(J,II) *T(1:3,3)                    !K11
                  BBEN(1:3,KK,2)= BBEN(1:3,KK,2) +  A(I,1)*A(J,II) *T(1:3,3)                    !K22
                  BBEN(1:3,KK,3)= BBEN(1:3,KK,3) - (A(I,1)*B(J,II)+B(I,1)*A(J,II)) *T(1:3,3)    !K12
	        END DO
    	
C	CONTRIBUTIONS FROM NORMAL DISPLACEMENTS OF THE MAIN ELEMENT
	        N1 = N+1 ! POSITION OF SIDE DERIVATIVE 2
	        C11 = DOT_PRODUCT(T(1:3,N) ,D1)
	        C12 = DOT_PRODUCT(T(1:3,N) ,D2)
	        C21 = DOT_PRODUCT(T(1:3,N1),D1)
	        C22 = DOT_PRODUCT(T(1:3,N1),D2)
	    
	        DO J=1,3                         ! FOR EACH NODE
                  C1 = C11*B(J,1) - C12*A(J,1)
                  C2 = C21*B(J,1) - C22*A(J,1)
                  BBEN(1:3,J,1) = BBEN(1:3,J,1) - B(I,1)*C1*T(1:3,3)            !K11
                  BBEN(1:3,J,2) = BBEN(1:3,J,2) + A(I,1)*C2*T(1:3,3)            !K22
                  BBEN(1:3,J,3) = BBEN(1:3,J,3) +(A(I,1)*C1-B(I,1)*C2)*T(1:3,3)  !K12
	        END DO
    	
	    ELSEIF(IFFIX(I).EQ.1)THEN !CLAMPED SIDE 

C	SIDE NORMAL COMPONENTS N1,N2
	        BML  = 1D0/SQRT(B(I,1)**2+A(I,1)**2)  ! 2 A(0)/L(I)
	        R11 =  B(I,1)*BML                     ! N1
	        R12 = -A(I,1)*BML                     ! N2
	        
C	 PROJECTION OF DERIVATIVES OF INITIAL CONFIGURATION ALONE 
C	 NORMAL VECTOR TO SYMMETIC PLANE
	        TI(1:3,1) = -MATMUL(COORD(1:3,1:3),B(1:3,1))     ! SAI(1)
	        TI(1:3,2) = +MATMUL(COORD(1:3,1:3),A(1:3,1))     ! SAI(2)
    	    
	        AN(1:3) = -2D0*(R11*TI(1:3,1)+R12*TI(1:3,2))	   ! N1*T1+N2*T2
	        XB = -R12*A(1:3,1)+R11*B(1:3,1)    ! S1*A+S2*B  B-LOCAL = (1,R,1-R)/BML
	        C2 = +     R11*B(I,1)              ! N(I),1 * N(I)1
	        S2 = -     R12*A(I,1)              ! N(I),2 * N(I)2
	        CS = - 2D0*R11*A(I,1)              ! N(I),2 * N(I)1 + N(I),1 * N(I)2
	        DN(1:3) = R11*D1(1:3)+R12*D2(1:3)
    		   
C	 N1*D1+N2*D2 DN (CONTRAVARIANT BASE VECTOR ALONG THE NORMAL VECTOR TO SYMMETRIC AXIS )
	        C11 = DOT_PRODUCT(AN,DN)           ! NORMAL VECTOR TO SYMMETRY PLANE
C	CONTRIBUTIONS FROM NORMAL DISPLACEMENTS OF THE MAIN ELEMENT
	        DO J=1,3                           ! FOR EACH NODE
      	        BNN = C11*XB(J)*FACCMP									!KNN
      	        BBEN(1:3,J,1) = BBEN(1:3,J,1) + C2*BNN*T(1:3,3)
      	        BBEN(1:3,J,2) = BBEN(1:3,J,2) + S2*BNN*T(1:3,3)
      	        BBEN(1:3,J,3) = BBEN(1:3,J,3) + CS*BNN*T(1:3,3)
	        END DO
    	
    	
	    ELSE ! SIMPLE SUPPORTED OR FREE
	    
	        NS      = NS+1                    !INCREASE NUMBER OF S.S. OR FREE SIDES
	        NSS(NS) = I                       !KEEP SIDE ORDER
	        
	    END IF
	    
	    N = N+2                               !POSITION OF SIDE DERIVATIVE 1
100   CONTINUE
C	-------------------------------------------------------------------------- 
C	-------------------------------------------------------------------------- 

         
C	--------------------------------------------------------------------------              
C	ACCORDING TO THE NUMBER OF FREE OR SS SIDES
C	-------------------------------------------------------------------------- 
	SELECT CASE(NS) 

C     --------------
	CASE (1)
C     --------------
C	ONE SIMPLE SUPPORTED SIDE
	    I = NSS(1)                           !SIDE ORDER
	    BML  = SQRT(B(I,1)**2+A(I,1)**2)     !SIDE PSEUDO LENGTH
C	NORMALIZES SIDE TO COMPUTE NORMAL
	    C1 =  B(I,1)/BML                     
	    C2 = -A(I,1)/BML

	    R11 = C1*C1                          !AUXILIAR FACTORS
	    R12 = C1*C2
	    R22 = C2*C2

      BMN = RESHAPE ((/ 1D0-R11*R11,    -R22*R11,    -2*R12*R11, 
     1                     -R11*R22, 1D0-R22*R22,    -2*R12*R22,   
     2                     -R11*R12,    -R22*R12,     1D0-2*R11*R22 /), (/3,3/))

	    CALL PROMA2(BBEN,BBEN,BMN,18,3,3)
    	
c	    BBEN = SBBEN - BBEN
	
C     --------------
	CASE (2)
C     --------------

C	TWO SIDES FREE OR SIMPLE SUPPORTED
C	COMPUTE THE NORMAL VECTORS
	    DO J=1,2                                    
	    I = NSS(J)                                !FREE SIDE
	    BML  = SQRT(B(I,1)**2+A(I,1)**2)          !SIDE PSEUDO LENGTH
	    BNO(1:2,J) = (/ B(I,1), -A(I,1) /)/BML    !SIDE NORMAL
	    END DO
	    
C	SECOND COMPUTE THE DUAL BASE
C	DETERMINANT OF THE BASE
	    BML = BNO(1,1)*BNO(2,2) - BNO(1,2)*BNO(2,1)    	
C	INVERSE OF MATRIX NO
	    BND(1,1) = BNO(2,2)/BML                       
	    BND(1,2) = -BNO(1,2)/BML
	    BND(2,1) = -BNO(2,1)/BML
	    BND(2,2) =  BNO(1,1)/BML
	    
C	AUXILIAR FACTORS
	    R11 = BND(1,1)*BND(1,2)
	    R22 = BND(2,1)*BND(2,2)
	    R12 = BND(1,1)*BND(2,2)+BND(2,1)*BND(1,2)
	    C11 = 2D0*BNO(1,1)*BNO(1,2)
	    C22 = 2D0*BNO(2,1)*BNO(2,2)
	    C12 = BNO(1,1)*BNO(2,2)+BNO(2,1)*BNO(1,2)

	BMN = RESHAPE ((/ 1D0-R11*C11,    -R22*C11, -R12*C11,
	1                     -R11*C22, 1D0-R22*C22, -R12*C22, 
	2                     -R11*C12, -R22*C12, 1D0-R12*C12 /), (/3,3/))

	    CALL PROMA2(BBEN,BBEN,BMN,18,3,3)
          
c          BBEN = BBEN - SBBEN
      
C	-------------------------------------------------------------------------- 
	ENDSELECT
C	-------------------------------------------------------------------------- 


	RETURN
	END

C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE PROMA2(A,B,C,N1,N2,N3)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER(I-N)
C	THIS ROUTINE EVALUATES A MATRIX PRODUCT
C                                                T
C           A(I,J) = B(I,K) * C(J,K)    A = B * C

	DIMENSION B(N1,N3),C(N2,N3),A(N1,N2)

	A = MATMUL(B,TRANSPOSE(C))

	RETURN

	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
		
	SUBROUTINE BMEM3(A,B,T,BMEM)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4(I-N)
C	-----------------------------------------------------
C	CALCULATE MEMBRAIN STRAIN
C	MEMBRANE MATRIX  CONSTANT STRAIN TRIANGLE
C	-----------------------------------------------------
	DIMENSION A(3),B(3),T(3,2)
	DIMENSION BMEM(3,6,3)

	DO J = 1,3 !NODE IN MAIN ELEMENT
	    BMEM(1:3,J,1) =  -B(J)*T(1:3,1)
	    BMEM(1:3,J,2) =   A(J)*T(1:3,2)
	    BMEM(1:3,J,3) =   A(J)*T(1:3,1)-B(J)*T(1:3,2)
	END DO

	RETURN
	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE CARTD3(A,B,T,COORD,ISIDE,IFFIX,SLAM)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------
C	   CALCULATE CARTESIAN DERIVATIVES
C	-----------------------------------------------------
C	 COORD(3,6)	= CURRENT NODAL CO ORDINATES OF THE ELEMENT PATCH
C	 T(3,9)		= DERIVATIVES OF THE ELEMENT N CONFIGURATIO OR CARTESIAN DERIVATIVES
C			    = HERE T IS CHANGING AND UPDATING T3 FOR CUirrENT CONFIGURATION
C	 B	        = SIDE PROYECTIONS ON LOCAL X1
C	 A	        = SIDE PROYECTIONS ON LOCAL X1
C	 ikk	        = SIDE ELEMENT CONNECTIVITIES
C	-----------------------------------------------------
	DIMENSION COORD(3,6),B(3,4),A(3,4),T(3,9)
	DIMENSION KK(3,3),ISIDE(3),IFFIX(3)

	DATA  KK /4,3,2, 5,1,3, 6,2,1 /

C	INITIALIZES DERIVATIVES
	T = 0.0D0         

	N = 4 !FIRST POSITION FOR SIDE ELEMENT VECTORS IN ARRAY T
	DO I = 1,4 !FOR EACH TRIANGLE IN THE PATCH
	    IF( I .EQ. 1)THEN !FOR MAIN TRIANGLE
	        T(1:3,1) = -MATMUL(COORD(1:3,1:3),B(1:3,1))       ! SAI(1)
	        T(1:3,2) = +MATMUL(COORD(1:3,1:3),A(1:3,1))       ! SAI(2)
	        T(1,3) = SLAM*(T(2,1)*T(3,2) - T(3,1)*T(2,2))     !NORMAL * RA
	        T(2,3) = SLAM*(T(3,1)*T(1,2) - T(1,1)*T(3,2))
	        T(3,3) = SLAM*(T(1,1)*T(2,2) - T(2,1)*T(1,2))
	    ELSE !FOR ADJACENT OR SIDE ELEMENTS
	        IF(ISIDE(I-1).GT.0)THEN !IF ADJACANT ELEMENT OR SIDE EXISTS
	            DO J = 1,3 !FOR EACH SIDE OF THE ADJACENT ELEMENT
                      K = KK(J,I-1)                                         !ASSOCIATED NODE       	
                      T(1:3,N)   = T(1:3,N)   - B(J,I)*COORD(1:3,K)   		!SAI(1)(I)
                      T(1:3,N+1) = T(1:3,N+1) + A(J,I)*COORD(1:3,K)   	    !SAI(2)(I)
                  END DO
	        END IF
	        N = N+2 !POINTER TO NEXT SIDE VECTORS IN ARRAY T
	    END IF
	END DO
		
	RETURN
	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE COINI(COORDI,COORD,EDIS,NLOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C    --------------------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C    --------------------------------------------------------------------
C      COORD(3,6) = CUirrENT NODALCOORDINATES
C     COORDI(3,6) = INITIAL NODAL COORDINATES
C     EDIS(18)    = CURRENT NODAL DISPLACEMENTS
C    --------------------------------------------------------------------
C
      DIMENSION COORD(3,6),COORDI(3,6),EDIS(1)

      K=0
      DO 30 I=1,6
      DO 30 J=1,3
	K=K+1
	COORDI(J,I) = 0.0D0
	IF(NLOPT.EQ.3) THEN
      COORDI(J,I) = COORD(J,I) - EDIS(K)
      ELSE
      COORDI(J,I) = COORD(J,I)
      ENDIF
30    CONTINUE

	RETURN
	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
     	SUBROUTINE ARETTR(COORD,T3,AREA2)
	IMPLICIT REAL*8 (A-H,O-Z)
     	IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C	CALCULETE UNIT NOIRNAL TO THE SHELL SURFACE AND AREA
C	OF THE UNDEFORMED ELEMENT
C     ------------------------------------------------------------
C     	SETS THE COORDINATES FOR ANY ELEMENT IN THIS PATCH
C     ------------------------------------------------------------
C     COORD(NCO,NNM)      = COORDINATES FOR PERTICULER ELEMENT
C	TTRE(3)	    	      = UNIT NORMALS PER MAIN ELEMENT
C	AREA			      = TWICE AREA OF THE ELEMENT
C     ------------------------------------------------------------
      DIMENSION COORD(3,6),EL1(3),EL2(3),T3(3)
C     ------------------------------------------------------------

C     INITIALIZE
        T3(1:3)  = 0.0
        EL1(1:3) = 0.0
        EL2(1:3) = 0.0

C	DEFINE CONSTANT
	EL1(1:3) = COORD(1:3,3) - COORD(1:3,2)                !SIDE 1
	EL2(1:3) = COORD(1:3,1) - COORD(1:3,3)                !SIDE 2
	
C     EVALUATE THE CROSS PRODUCT =>  NORMAL TO THE PLAN OF SHELL ELEMENT
	T3(1) = EL1(2)*EL2(3) - EL1(3)*EL2(2)                 !NORMAL * AREA2
	T3(2) = EL1(3)*EL2(1) - EL1(1)*EL2(3)
	T3(3) = EL1(1)*EL2(2) - EL1(2)*EL2(1)

	AREA2=SQRT(T3(1)*T3(1)+T3(2)*T3(2)+T3(3)*T3(3))	    !COMPUTES TWICE AREA
	T3(1:3) = T3(1:3)/AREA2
	
	RETURN
	END
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE SHAPD3(COORD,A,B,AREA2,ISIDE)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------------
C	SET UP UNDEFOREMED CORDINATE TO GET SIDE VECTORES AND UNIT NOMALS IN LOCAL
C	DIIReCTIONS, AND THROUGH THEIR PROJECTIONS OVER LOCAL
c	DIIReCTIONS IN UNDEFORMED CONFIGURATION

C     ---------------------------------------------------------------------
C     SETS THE COORDINATES FOR ANY ELEMENT IN THIS PATCH
C     ---------------------------------------------------------------------
C     COORD(NCO,NNM)        = COORDINATES FOR PERTICULER ELEMENT
C    	LM(NNM)               = ELEMENT COONECTIONS (NODAL INCIDENCES)
C     SIMI(3,3)       	    = MAIN ELEMENT COORDINATE PER ELEMENT
C	A(3,6)		        = PROJECTION OVER LOCAL DIIReCTIN T1
C	B(3,6)		        = PROJECTION OVER LOCAL DIIReCTIN T2
C	T1(3)		            = UNIT VECTOR ALONE LOCAL DIIReCTION ONE
C	T2(3)		            = UNIT VECTOR ALONE LOCAL DIIReCTION TWO
C	T3(3)		            = UNIT NORMAL ALONE LOCAL DIIReCTION 3
C	COORDI(3,6)           = INTIAL COORDINATES AT NODES IN THE ELEMENT (X(3,6))
C	T(3,6)		        = T(3,1:3),T1,T2,T3 FOR MAIN ELEMENT
C				              T(3,4:6),T3 FOR ALL THE PATCH ELEMENTS
C     AREA2(4)        	    = 2 TIMES AREA OF MAIN ELEMENT AND ADJACENT ELEMENTS
C	ISIDE(3)		        = IF GT 1 SIDE ELEMENT( ADJACANT) ELEMENT EXIST
C     ---------------------------------------------------------------------
	DIMENSION COORD(3,6),T(3,6),A(3,4),B(3,4),AREA2(4)
	DIMENSION ISIDE(3),KK(3,3),LL(2,3)
	DIMENSION EL1(3),EL2(3),T1(3),T2(3),T3(3)

	DATA	 KK / 4,3,2, 5,1,3, 6,2,1 /
	DATA	 LL /  3,2, 1,3, 2,1 /

C	 FOR EXISTING ELEMENTS ISIDE(IJ)=NODE>0
C	 FOR MISSING ELEMENT ISIDE(IJ)=NODE=<0


c     -----------------------------------------------------------------
c     INITIALIZE
c     -----------------------------------------------------------------
        T(1:3,1:6)  = 0.0
        A(1:3,1:4)  = 0.0
        B(1:3,1:4)  = 0.0
        AREA2(1:4)  = 0.0
        EL1(1:3)    = 0.0
        EL2(1:3)    = 0.0     

C	-------------------------------------------------------------------------------------------
C	 EVALUATE THE FIRST TWO SIDE VECTORS
 
	EL1(1:3) = COORD(1:3,3) - COORD(1:3,2)                  !SIDE 1
	EL2(1:3) = COORD(1:3,1) - COORD(1:3,3)                  !SIDE 2

C	-------------------------------------------------------------------------------------------
C	 EVALUATE THE CROSS PRODUCT =>  NORMAL TO THE PLANE OF SHELL ELEMENT

	T3(1) = EL1(2)*EL2(3) - EL1(3)*EL2(2)	  !NORMAL * AREA2
	T3(2) = EL1(3)*EL2(1) - EL1(1)*EL2(3)
	T3(3) = EL1(1)*EL2(2) - EL1(2)*EL2(1)
        
	
	AREA2(1) = SQRT(T3(1)*T3(1)+T3(2)*T3(2)+T3(3)*T3(3))	
C	COMPUTES TWICE AREA

	IF (AREA2(1) == 0.0D0)THEN
	    WRITE (*,*) 'ERROR: AREA2 == 0 '
	    WRITE(*,"(3E15.5)")COORD(1:3,1:3)
	    STOP 'SHAPD3: WRONG ELEMENT DEFINITION   '
	END IF
	
	T3 = T3/AREA2(1)      !T3 (UNIT NORMAL)
C	-------------------------------------------------------------------------------------------
C	SELECT LOCAL X ,AND ITS DIRECTION T1 IN THE GLOBAL XY PLANE(T3.T1=0)

	T1 = (/ -T3(2), T3(1) , 0D0 /)
C	-------------------------------------------------------------------------------------------
C	SELECT LOCAL Y ,AND ITS DIRECTION T2 = T3 X T1 IN THE GOLOBAL XYZ PLANE
	T2 = (/ -T3(1)*T3(3), -T3(2)*T3(3), (T3(1)*T3(1)+T3(2)*T3(2)) /)

C     ----------------------------
	IF(ABS(T2(3)) < 1.0D-5) THEN
	    T1 = (/  1D0, 0D0, 0D0 /)
	    IF(T3(3) > 0D0 )THEN
	      T2 = (/  0D0, 1D0, 0D0 /)
	    ELSE
	      T2 = (/  0D0,-1D0, 0D0 /)
	    END IF
	ELSE
C	NORMALIZES T1 & T2 ( UNIT NORMALS)
	    T1 = T1/SQRT(DOT_PRODUCT(T1,T1))
	    T2 = T2/SQRT(DOT_PRODUCT(T2,T2))
	END IF
C     ----------------------------

	T(1:3,1) = T1		!LOCAL X1 DIRECTION
	T(1:3,2) = T2		!LOCAL X2 DIRECTION
	T(1:3,3) = T3		!NORMAL DIRECTION
C	-------------------------------------------------------------------------------------------
C	FIND THE LOCAL COORDINATES

	A(1,1) = EL1(1)*T(1,1)+EL1(2)*T(2,1)+EL1(3)*T(3,1)      ! L1 . T1
	A(2,1) = EL2(1)*T(1,1)+EL2(2)*T(2,1)+EL2(3)*T(3,1)      ! L2 . T1
	A(3,1) = -A(1,1)-A(2,1)
	B(1,1) = EL1(1)*T(1,2)+EL1(2)*T(2,2)+EL1(3)*T(3,2)      ! L1 . T2
	B(2,1) = EL2(1)*T(1,2)+EL2(2)*T(2,2)+EL2(3)*T(3,2)      ! L2 . T2
	B(3,1) = -B(1,1)-B(2,1)

C     ---------------------------
C     ---------------------------
	DO 100 II = 1,3
	IF(ISIDE(II).EQ.0) CYCLE  !NO ADJACENT ELEMT.
	              
	JJ = II+1                                           	!NEXT ELEMENT
	I = KK(1,II)                                        	!LOCAL FIRST NODE
	J = KK(2,II)                                        	!LOCAL SECOND NODE
	K = KK(3,II)                                        	!LOCAL THIRD NODE
	EL1 = COORD(1:3,K) - COORD(1:3,J)                      !SIDE J-K  (I)
	EL2 = COORD(1:3,I) - COORD(1:3,K)                      !SIDE K-I  (J)

	ELL1= SQRT(EL1(1)*EL1(1)+EL1(2)*EL1(2)+EL1(3)*EL1(3))  !LENGTH OF SIDE (I)
	ELL2= SQRT(EL2(1)*EL2(1)+EL2(2)*EL2(2)+EL2(3)*EL2(3))  !LENGTH OF SIDE (J)

	T3(1) = EL1(2)*EL2(3) - EL1(3)*EL2(2)                 !NORMAL OF ELEM (I)
	T3(2) = EL1(3)*EL2(1) - EL1(1)*EL2(3)
	T3(3) = EL1(1)*EL2(2) - EL1(2)*EL2(1)
	AREA2(JJ) = SQRT(T3(1)*T3(1)+T3(2)*T3(2)+T3(3)*T3(3))    !AREA OF ELEM (I)
	
	IF (AREA2(JJ) == 0.0D0)THEN
	    WRITE (*,*) 'ERROR: AREA2 == 0 '
	    WRITE(*,"(3E15.5)")COORD(1:3,I),COORD(1:3,J),COORD(1:3,K)
	    STOP 'SHAPD3: WRONG ELEMENT DEFINITION   '
	END IF

	T(1:3,I) = T3/AREA2(JJ)   !NORMAL OF PATCH ELEM (I)
C	-------------------------------------------------------------------------------------------
C	PROJEC. OF SIDE (I) OF ADJACENT ELEMENT(I)

	A(1,JJ) = -A(II,1)                                  
	B(1,JJ) = -B(II,1)
C	-------------------------------------------------------------------------------------------
C	ANGLE OF  SIDE (I) OF ADJACENT ELEMENT(I)
	
	COSA = A(1,JJ)/ELL1                                  
	SINA = B(1,JJ)/ELL1
C	-------------------------------------------------------------------------------------------
C	ANGLE BETWEEN (I-J) OF ADJACENT ELEMENT(I)

	COSB = (EL1(1)*EL2(1)+EL1(2)*EL2(2)+EL1(3)*EL2(3))/ELL1/ELL2
	SINB = AREA2(JJ)/ELL1/ELL2
C	-------------------------------------------------------------------------------------------
C	ANGLE OF SIDE (J) GAMA = ALPA+BETA
	COSG = COSA*COSB - SINA*SINB                
	SING = COSA*SINB + SINA*COSB
C	-------------------------------------------------------------------------------------------
C	PROJEC. OF SIDE (J) OF ADJACENT ELEMENT(I)
	A(2,JJ) = ELL2*COSG                                  
	B(2,JJ) = ELL2*SING
C	-------------------------------------------------------------------------------------------
C	PROJEC. OF SIDE (K) OF ADJACENT ELEMENT(I)
	A(3,JJ) = -A(1,JJ)-A(2,JJ)                          
	B(3,JJ) = -B(1,JJ)-B(2,JJ)
 
100   CONTINUE
C     ---------------------------
C     ---------------------------

	DO 200 I=1,4
	    IF(I.GT.1) THEN
	        IF(ISIDE(I-1).EQ.0) THEN
	            A(1:3,I) = 0D0
	            B(1:3,I) = 0D0
	            GOTO 200
	        ENDIF	
	    ENDIF
	    A(1:3,I) = A(1:3,I)/AREA2(I)
	    B(1:3,I) = B(1:3,I)/AREA2(I)
200   CONTINUE
	
C	WRITE(7,"('COORDINATES',/,(3E15.5))") COORDI(:,1:6)
C	WRITE(7,"('T1 AND T2  ',/,(3E15.5))") T(:,1:2)
C	WRITE(7,"('2 X AREAS  ',/,(4E15.5))") AREA2(1:4)
C	WRITE(7,"('A/2A       ',/,(3E15.5))") A(:,1:4)
C	WRITE(7,"('B/2A       ',/,(3E15.5))") B(:,1:4)

	RETURN
	END

C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE ESUPFACE(IGIDM,MCONT,ISFAC,NEFC,NEF,NELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     READS, ASSIGN SUPPORT FLAG FOR SHELL WITH NO ROTATION DOF ONATE
C	-------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

C     NEFC = NUMBER OF ELEMENT FACE = NEFC  (SEE ELEMIN)
C     ISFAC -- STORE SUPPORT FLAG FOR EACH ELEMENT FACE (SEE SHELL 3 NODE ONATE)
      DIMENSION IGIDM(1),ISFAC(2*NEFC,1),MCONT(NEF,1)

C     INITIALIZE
      ISFAC(1:2*NEFC,1:NELE) = 0
      
      READ(ITI,*)
      READ(ITI,*) NSUPP
      
C     ----------------------      
      DO 200 ISUPP = 1,NSUPP
      
      READ(ITI,*) IGM,IFACE,ISFLAG
      CALL ELEREODER(IGIDM,NELE,IGM) !DETERMINE THE ELEMENT NUMBER CORRESPONDING TO INPUT GID ELEMENT NUMBER IGM
      IF(IGM.EQ.0) GOTO 200
      
      IF(IFACE.GT.0.AND.IFACE.LE.NEFC) THEN
          CALL FLIPFACE3(IFACE,IFACN)
          ISFAC(IFACN,IGM) = ISFLAG !ISFLAG -- 0=FREE&SIMPLY, 1=CLAMP
      ENDIF
      
200   CONTINUE
C     ----------------------   

      DO I = 1,NEFC
        ISFAC(NEFC+I,1:NELE) = MCONT(NEFC+I,1:NELE)
      ENDDO
      
      
      RETURN
      END
C
C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE FLIPFACE3(IFACO,IFACN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      SELECTCASE(IFACO)
      CASE(1)
      IFACN = 3
      CASE(2)
      IFACN = 1
      CASE(3)
      IFACN = 2
      ENDSELECT
      
      RETURN
      END
C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE ELEMLMS(MCONT,NELE,NEF,NEFC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C	-------------------------------------------------------------
      DIMENSION MCONT(NEF,NELE),ICnA(NEFC,NELE)


C     ---------------------------
C     GETTING ICNA
C     ---------------------------
	! Initialize ICnA Matrix
      ICnA(1:NEFC,1:NELE) = 0

	DO ICD = 1,NELE
        NumADE = 0; ! Number of Adjacent Element
        DO ICA = 1,NELE !Read elemetn to seairch Adjacent Element
            IF(ICD.NE.ICA) THEN ! Compare the nodes of Control domain with 
            ! those of all element (except Control domain)
                ICount = 0
                DO I = 1,NEFC
                    ircO = MCONT(I,ICD) ! The nodes of Control domain
                    IF(ircO.LE.0) CYCLE
                    DO J = 1,NEFC !The Nodes of all element (except control domain)
                        irc = MCONT(J,ICA)
                        IF(irc.LE.0) CYCLE
                        IF(ircO.EQ.irc) THEN
                            ICount = ICount + 1;
                            IF(ICount.EQ.2) THEN ! 2 nodes are same so the elements are connect together
                                NumADE = NumADE + 1;
                                ICnA(NumADE,ICD) = ICA ! Control Domain & Adjacent Element Matrix, First Number = Control Domain Number, 2~4 Number = Adjacent Element Number
                            ENDIF 
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
	ENDDO
C     ---------------------------
C     ---------------------------


C     ---------------------------
C     GETTING CONTROL DOMAIN CONNECTIVITY
C     ---------------------------
	DO IELE = 1,NELE
	  CALL LMSCONTR3(MCONT,ICnA,IELE,NEF)
	ENDDO



      RETURN
      END
C
C     =====================================================================
C     =====================================================================
C     =====================================================================
      SUBROUTINE LMSCONTR3(MCONT,ICnA,IELE,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C	-------------------------------------------------------------
      DIMENSION MCONT(NEF,1),ICnA(3,1),KK(2,3)

	DATA KK / 2,3,  3,1,  1,2 / 

	! Initialize LM Matrix
      MCONT(4:6,IELE) = 0

      DO ISIDE = 1,3
          II = KK(1,ISIDE) ; JJ = KK(2,ISIDE)
          NODI = MCONT(II,IELE) ; NODJ = MCONT(JJ,IELE)
          
          DO IADJE = 1,3
              NADJ = ICnA(IADJE,IELE)
              
              IF(NADJ.NE.0) THEN
                    NODA1 = MCONT(1,NADJ) ; NODA2 = MCONT(2,NADJ) ; NODA3 = MCONT(3,NADJ)
                        IF(NODA1.EQ.NODI.AND.NODA2.EQ.NODJ) GOTO 100
                        IF(NODA2.EQ.NODI.AND.NODA1.EQ.NODJ) GOTO 100
                        IF(NODA1.EQ.NODI.AND.NODA3.EQ.NODJ) GOTO 100
                        IF(NODA3.EQ.NODI.AND.NODA1.EQ.NODJ) GOTO 100
                        IF(NODA2.EQ.NODI.AND.NODA3.EQ.NODJ) GOTO 100
                        IF(NODA3.EQ.NODI.AND.NODA2.EQ.NODJ) GOTO 100
                    GOTO 200
100                 CONTINUE 
                        IF(NODA1.NE.NODI.AND.NODA1.NE.NODJ) THEN
                            MCONT(ISIDE+3,IELE) = NODA1 
                            EXIT
                        ENDIF
                        IF(NODA2.NE.NODI.AND.NODA2.NE.NODJ) THEN
                            MCONT(ISIDE+3,IELE) = NODA2  
                            EXIT
                        ENDIF
                        IF(NODA3.NE.NODI.AND.NODA3.NE.NODJ) THEN
                            MCONT(ISIDE+3,IELE) = NODA3  
                            EXIT
                        ENDIF     
200                 CONTINUE       
              ENDIF
          
          ENDDO
      
      ENDDO
      

      RETURN
      END
C
C     =====================================================================
C     =====================================================================
C     =====================================================================

      SUBROUTINE S3MDSH(COORD,COORDI,REDIS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------
C     MODIFIES TOTAL DISPLACEMENT VECTOR BY DEDUCTING
C     RIGID BODY TRANSLATIONS
C     ----------------------------------------------------
      DIMENSION COORD(3,1),COORDI(3,1),REDIS(1)
      DIMENSION XYZ(3),XYZO(3),CD(3,8),CDO(3,8)
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION VRO(3),VSO(3),VTO(3),TM(3,3)
C
C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
      CALL CLEARA (XYZO,3)
	CALL CLEARA (XYZ ,3)
	
      DO 140 I=1,3
      DO 140 J=1,3
	XYZO(J)=XYZO(J)+COORDI(J,I)/3.0D0
  140 XYZ (J)=XYZ (J)+COORD (J,I)/3.0D0
  
      DO 150 I=1,6
      DO 150 J=1,3
	CDO(J,I)=COORDI(J,I)-XYZO(J)
  150 CD (J,I)=COORD (J,I)-XYZ (J)
  
  	CALL S3BASEV(COORDI,VRO,VSO,VTO)
  	CALL S3BASEV(COORD ,VR ,VS ,VT )
      
      DO 160 I=1,3
      DO 160 J=1,3
 160  TM(I,J)=VR(I)*VRO(J)+VS(I)*VSO(J)+VT(I)*VTO(J)
 
      K=1
      DO 210 I=1,6
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)
      REDIS(K)=CD(J,I)-TCDO
  210 K=K+1
  
      RETURN
      END
C
C     =====================================================================
C     =====================================================================
C     =====================================================================
	SUBROUTINE S3BASEV(COORD,T1,T2,T3)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------------
C     ---------------------------------------------------------------------
	DIMENSION COORD(3,1)
	DIMENSION T1(3),T2(3),T3(3)


C	-------------------------------------------------------------------------------------------
C	 EVALUATE THE FIRST TWO SIDE VECTORS
 
	T1(1:3) = COORD(1:3,3) - COORD(1:3,1)                  !SIDE 12
	T2(1:3) = COORD(1:3,2) - COORD(1:3,1)                  !SIDE 13

      CALL VECPRD (T1,T2,T3)
      CALL SCALEN (T3,T3,DUM,3)
      CALL SCALEN (T1,T1,DUM,3)
      CALL VECPRD (T3,T1,T2)
      CALL SCALEN (T2,T2,DUM,3)
      

      RETURN
      END
C
C     =====================================================================
C     =====================================================================
C     =====================================================================





