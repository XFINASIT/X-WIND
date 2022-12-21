C       ============================================================
C       ============================================================
C       =========================================================
C       Purpose: This program computes the Bessel functions  
C                Jn(x) and Yn(x) ( n=0,1 ) and their derivatives 
C                using subroutine JY01A
C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ?0 )
C       Output:  BJ0 --- J0(x)
C                DJ0 --- J0'(x)
C                BJ1 --- J1(x)
C                DJ1 --- J1'(x)
C                BY0 --- Y0(x)
C                DY0 --- Y0'(x)
C                BY1 --- Y1(x)
C                DY1 --- Y1'(x)
C       Example:
C
C        x       J0(x)        J0'(x)       J1(x)        J1'(x)
C       ---------------------------------------------------------
C        1     .76519769   -.44005059    .44005059    .32514710
C        5    -.17759677    .32757914   -.32757914   -.11208094
C       10    -.24593576   -.04347275    .04347275   -.25028304
C       20     .16702466   -.06683312    .06683312    .16368301
C       30    -.08636798    .11875106   -.11875106   -.08240961
C       40     .00736689   -.12603832    .12603832    .00421593
C       50     .05581233    .09751183   -.09751183    .05776256
C
C        x       Y0(x)        Y0'(x)       Y1(x)        Y1'(x)
C      ---------------------------------------------------------
C        1     .08825696    .78121282   -.78121282    .86946979
C        5    -.30851763   -.14786314    .14786314   -.33809025
C       10     .05567117   -.24901542    .24901542    .03076962
C       20     .06264060    .16551161   -.16551161    .07091618
C       30    -.11729573   -.08442557    .08442557   -.12010992
C       40     .12593642    .00579351   -.00579351    .12608125
C       50    -.09806500    .05679567   -.05679567   -.09692908
C       =========================================================
C       ============================================================
C       ============================================================
        SUBROUTINE JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)

C
C       =======================================================
C       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
C                Y1(x), and their derivatives
C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ?0 )
C       Output:  BJ0 --- J0(x)
C                DJ0 --- J0'(x)
C                BJ1 --- J1(x)
C                DJ1 --- J1'(x)
C                BY0 --- Y0(x)
C                DY0 --- Y0'(x)
C                BY1 --- Y1(x)
C                DY1 --- Y1'(x)
C       =======================================================
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        RP2=0.63661977236758D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           DJ0=0.0D0
           DJ1=0.5D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           DY0=1.0D+300
           DY1=1.0D+300
           RETURN
        ENDIF
        IF (X.LE.12.0D0) THEN
           BJ0=1.0D0
           R=1.0D0
           DO 5 K=1,30
              R=-0.25D0*R*X2/(K*K)
              BJ0=BJ0+R
              IF (DABS(R).LT.DABS(BJ0)*1.0D-15) GO TO 10
5          CONTINUE
10         BJ1=1.0D0
           R=1.0D0
           DO 15 K=1,30
              R=-0.25D0*R*X2/(K*(K+1.0D0))
              BJ1=BJ1+R
              IF (DABS(R).LT.DABS(BJ1)*1.0D-15) GO TO 20
15         CONTINUE
20         BJ1=0.5D0*X*BJ1
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           CS0=0.0D0
           W0=0.0D0
           R0=1.0D0
           DO 25 K=1,30
              W0=W0+1.0D0/K
              R0=-0.25D0*R0/(K*K)*X2
              R=R0*W0
              CS0=CS0+R
              IF (DABS(R).LT.DABS(CS0)*1.0D-15) GO TO 30
25         CONTINUE
30         BY0=RP2*(EC*BJ0-CS0)
           CS1=1.0D0
           W1=0.0D0
           R1=1.0D0
           DO 35 K=1,30
              W1=W1+1.0D0/K
              R1=-0.25D0*R1/(K*(K+1))*X2
              R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS1=CS1+R
              IF (DABS(R).LT.DABS(CS1)*1.0D-15) GO TO 40
35         CONTINUE
40         BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (X.GE.35.0) K0=10
           IF (X.GE.50.0) K0=8
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO 45 K=1,K0
              P0=P0+A(K)*X**(-2*K)
45            Q0=Q0+B(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO 50 K=1,K0
              P1=P1+A1(K)*X**(-2*K)
50            Q1=Q1+B1(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN
        END
C       ============================================================
C       ============================================================
        SUBROUTINE Diffraction_Force_Parameter(ak,Aka,Cm,Ch,Delta)
        ! ka = ak
        IMPLICIT REAL*8 (A-H,O-Z)
        CALL JY01A(ak,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
        pi = 3.1416d0
        Aka = 1d0/( ( DJ1**2d0 + DY1**2d0 )**(0.5d0) )
        Cm  = 2d0*Aka/pi/ak/BJ1
        Ch  = 4d0*Aka/pi/ak**2d0
        Delta  = -atan(DY1/DJ1)
        end
C       ============================================================
C       ============================================================        
        SUBROUTINE Diffraction_Force(ak,z,a,d,period,t,
     &Roh,H,g,Fx,Fy,x,wavenumber)
        ! ka = ak
        IMPLICIT REAL*8 (A-H,O-Z)
     
        CALL Diffraction_Force_Parameter(ak,Aka,Cm,Ch,Delta)
        pi = 3.1416d0
        OMEGA = 2d0*pi/period
        cosss = cos(OMEGA*t-wavenumber*x)
        sinnn = sin(OMEGA*t-wavenumber*x)
       Fx=Roh*g*H*a*Aka/ak*cosh(wavenumber*z)/cosh(wavenumber*d)
     &* cosss 
       
       Fy=Roh*g*H*a*Aka/ak*sinh(wavenumber*z)/sinh(wavenumber*d)
     &* sinnn  
        RETURN
        END SUBROUTINE
C       ============================================================
C       ============================================================        
      SUBROUTINE Diffraction_Pressure(m,ak,z,a,d,period,t,
     &Roh,H,g,Fx,x,wavenumber,zeta)
        
        
        IMPLICIT REAL*8 (A-H,O-Z)
        
        complex(8)  Beta 
        complex(8)  Bi
        
       ! complex(8)  Z 
       ! complex(8)  NM
        complex(8)  CHF1
        complex(8)  CHD1
        complex(8)  CHF2
        complex(8)  CHD2
        
        SumP = 0d0
        
        Do j = 1,100
        
        m=j
        
        if (m.EQ.0d0) then
        Beta = (1.0 ,0)
        else
        
        Bi   = (0 , 1.0 )
        Beta = 2d0*(Bi**m)
        endif
        
        !call CH12N(N,Z,NM,CHF1,CHD1,CHF2,CHD2)
        
        ! P = ( Roh * g * H )*( z/d + cosh(ks)/cosh(kd) )
        
        SumP = SumP + Beta*cos(m*zeta)
        
        
        
        enddo
        
        
        
        
        CALL Diffraction_Force_Parameter(ak,Aka,Cm,Ch,Delta)
        pi = 3.1416d0
        OMEGA = 2d0*pi/period
        cosss = cos(OMEGA*t-wavenumber*x)
       Fx=Roh*g*H*a*Aka/ak*cosh(wavenumber*z)/cosh(wavenumber*d)
     &* cosss 
     
       
        end
C       ============================================================
C       ============================================================        
        SUBROUTINE CH12N(N,Z,NM,CHF1,CHD1,CHF2,CHD2)
C
C       ====================================================
C       Purpose: Compute Hankel functions of the first and
C                second kinds and their derivatives for a
C                complex argument
C       Input :  z --- Complex argument
C                n --- Order of Hn(1)(z) and Hn(2)(z)
C       Output:  CHF1(n) --- Hn(1)(z)
C                CHD1(n) --- Hn(1)'(z)
C                CHF2(n) --- Hn(2)(z)
C                CHD2(n) --- Hn(2)'(z)
C                NM --- Highest order computed
C       Routines called:
C             (1) CJYNB for computing Jn(z) and Yn(z)
C             (2) CIKNB for computing In(z) and Kn(z)
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250),
     &            CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)
        DIMENSION CHF1(0:N),CHD1(0:N),CHF2(0:N),CHD2(0:N)
        !CBJ = cmplx(
        CI=(0.0D0,1.0D0)
        PI=3.141592653589793D0
         XXX = DIMAG(Z)
         
        IF (DIMAG(Z).LT.0.0D0) THEN
           
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 10 K=0,NM
              CHF1(K)=CBJ(K)+CI*CBY(K)
10            CHD1(K)=CDJ(K)+CI*CDY(K)
           ZI=CI*Z
           CALL CIKNB(N,ZI,NM,CBI,CDI,CBK,CDK)
           CFAC=-2.0D0/(PI*CI)
           DO 15 K=0,NM
              CHF2(K)=CFAC*CBK(K)
              CHD2(K)=CFAC*CI*CDK(K)
15            CFAC=CFAC*CI
        ELSE IF (DIMAG(Z).GT.0.0D0) THEN
           ZI=-CI*Z
           CALL CIKNB(N,ZI,NM,CBI,CDI,CBK,CDK)
           CF1=-CI
           CFAC=2.0D0/(PI*CI)
           DO 20 K=0,NM
              CHF1(K)=CFAC*CBK(K)
              CHD1(K)=-CFAC*CI*CDK(K)
20            CFAC=CFAC*CF1
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 25 K=0,NM
            !  CHF2(K)=CBJ(K)-CI*CBY(K)
25          a = 1
            !  CHD2(K)=CDJ(K)-CI*CDY(K)
        ELSE
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 30 K=0,NM
              CHF1(K)=CBJ(K)+CI*CBY(K)
              CHD1(K)=CDJ(K)+CI*CDY(K)
              CHF2(K)=CBJ(K)-CI*CBY(K)
30            CHD2(K)=CDJ(K)-CI*CDY(K)
        ENDIF
        RETURN
        END
C       ============================================================
C       ============================================================
        SUBROUTINE CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
C
C       =======================================================
C       Purpose: Compute Bessel functions Jn(z), Yn(z) and
C                their derivatives for a complex argument
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =======================================================
C

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N),
     &            A(4),B(4),A1(4),B1(4)
        
        EL=0.5772156649015329D0
        PI=3.141592653589793D0
        R2P=.63661977236758D0
        Y0=DABS(DIMAG(Z))
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        IF (A0.LE.300.D0.OR.N.GT.80) THEN
           IF (N.EQ.0) NM=1
           M=MSTA1(A0,200)
           IF (M.LT.NM) THEN
              NM=M
           ELSE
              M=MSTA2(A0,NM,15)
           ENDIF
           CBS=(0.0D0,0.0D0)
           CSU=(0.0D0,0.0D0)
           CSV=(0.0D0,0.0D0)
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
                 IF (Y0.LE.1.0D0) THEN
                    CBS=CBS+2.0D0*CF
                 ELSE
                    CBS=CBS+(-1)**(K/2)*2.0D0*CF
                 ENDIF
                 CSU=CSU+(-1)**(K/2)*CF/K
              ELSE IF (K.GT.1) THEN
                 CSV=CSV+(-1)**(K/2)*K/(K*K-1.0D0)*CF
              ENDIF
              CF2=CF1
15            CF1=CF
           IF (Y0.LE.1.0D0) THEN
              CS0=CBS+CF
           ELSE
              CS0=(CBS+CF)/CDCOS(Z)
           ENDIF
           DO 20 K=0,NM
20            CBJ(K)=CBJ(K)/CS0
           CE=CDLOG(Z/2.0D0)+EL
           CBY(0)=R2P*(CE*CBJ(0)-4.0D0*CSU/CS0)
           CBY(1)=R2P*(-CBJ(0)/Z+(CE-1.0D0)*CBJ(1)-4.0D0*CSV/CS0)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02/
           CT1=Z-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 25 K=1,4
25            CP0=CP0+A(K)*Z**(-2*K)
           CQ0=-0.125D0/Z
           DO 30 K=1,4
30            CQ0=CQ0+B(K)*Z**(-2*K-1)
           CU=CDSQRT(R2P/Z)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CBJ(0)=CBJ0
           CBY(0)=CBY0
           CT2=Z-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 35 K=1,4
35            CP1=CP1+A1(K)*Z**(-2*K)
           CQ1=0.375D0/Z
           DO 40 K=1,4
40            CQ1=CQ1+B1(K)*Z**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
           CBJ(1)=CBJ1
           CBY(1)=CBY1
           DO 45 K=2,NM
              CBJK=2.0D0*(K-1.0D0)/Z*CBJ1-CBJ0
              CBJ(K)=CBJK
              CBJ0=CBJ1
45            CBJ1=CBJK
        ENDIF
        CDJ(0)=-CBJ(1)
        DO 50 K=1,NM
50         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        IF (CDABS(CBJ(0)).GT.1.0D0) THEN
           CBY(1)=(CBJ(1)*CBY(0)-2.0D0/(PI*Z))/CBJ(0)
        ENDIF
        DO 55 K=2,NM
           IF (CDABS(CBJ(K-1)).GE.CDABS(CBJ(K-2))) THEN
              CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           ELSE
              CYY=(CBJ(K)*CBY(K-2)-4.0D0*(K-1.0D0)/(PI*Z*Z))/CBJ(K-2)
           ENDIF
           CBY(K)=CYY
55      CONTINUE
        CDY(0)=-CBY(1)
        DO 60 K=1,NM
60         CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END

C       ============================================================
C       ============================================================
        SUBROUTINE CIKNB(N,Z,NM,CBI,CDI,CBK,CDK)
C
C       ============================================================
C       Purpose: Compute modified Bessel functions In(z) and Kn(z),
C                and their derivatives for a complex argument
C       Input:   z --- Complex argument
C                n --- Order of In(z) and Kn(z)
C       Output:  CBI(n) --- In(z)
C                CDI(n) --- In'(z)
C                CBK(n) --- Kn(z)
C                CDK(n) --- Kn'(z)
C                NM --- Highest order computed
C       Routones called:
C                MSTA1 and MSTA2 to compute the starting point for
C                backward recurrence
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:N),CDI(0:N),CBK(0:N),CDK(0:N)
        PI=3.141592653589793D0
        EL=0.57721566490153D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=(0.0D0,0.0D0)
              CBK(K)=(1.0D+300,0.0D0)
              CDI(K)=(0.0D0,0.0D0)
10            CDK(K)=-(1.0D+300,0.0D0)
           CBI(0)=(1.0D0,0.0D0)
           CDI(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        Z1=Z
        CI=(0.0D0,1.0D0)
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (N.EQ.0) NM=1
        M=MSTA1(A0,200)
        IF (M.LT.NM) THEN
           NM=M
        ELSE
           M=MSTA2(A0,NM,15)
        ENDIF
        CBS=0.0D0
        CSK0=0.0D0
        CF0=0.0D0
        CF1=1.0D-100
        DO 15 K=M,0,-1
           CF=2.0D0*(K+1.0D0)*CF1/Z1+CF0
           IF (K.LE.NM) CBI(K)=CF
           IF (K.NE.0.AND.K.EQ.2*INT(K/2)) CSK0=CSK0+4.0D0*CF/K
           CBS=CBS+2.0D0*CF
           CF0=CF1
15         CF1=CF
        CS0=CDEXP(Z1)/(CBS-CF)
        DO 20 K=0,NM
20         CBI(K)=CS0*CBI(K)
        IF (A0.LE.9.0) THEN
           CBK(0)=-(CDLOG(0.5D0*Z1)+EL)*CBI(0)+CS0*CSK0
           CBK(1)=(1.0D0/Z1-CBI(1)*CBK(0))/CBI(0)
        ELSE
           CA0=CDSQRT(PI/(2.0D0*Z1))*CDEXP(-Z1)
           K0=16
           IF (A0.GE.25.0) K0=10
           IF (A0.GE.80.0) K0=8
           IF (A0.GE.200.0) K0=6
           DO 30 L=0,1
              CBKL=1.0D0
              VT=4.0D0*L
              CR=(1.0D0,0.0D0)
              DO 25 K=1,K0
                 CR=0.125D0*CR*(VT-(2.0*K-1.0)**2)/(K*Z1)
25               CBKL=CBKL+CR
              CBK(L)=CA0*CBKL
30         CONTINUE
        ENDIF
        CG0=CBK(0)
        CG1=CBK(1)
        DO 35 K=2,NM
           CG=2.0D0*(K-1.0D0)/Z1*CG1+CG0
           CBK(K)=CG
           CG0=CG1
35         CG1=CG
        IF (REAL(Z).LT.0.0) THEN
           FAC=1.0D0
           DO 45 K=0,NM
              IF (DIMAG(Z).LT.0.0) THEN
                 CBK(K)=FAC*CBK(K)+CI*PI*CBI(K)
              ELSE
                 CBK(K)=FAC*CBK(K)-CI*PI*CBI(K)
              ENDIF
              CBI(K)=FAC*CBI(K)
              FAC=-FAC
45         CONTINUE
        ENDIF
        CDI(0)=CBI(1)
        CDK(0)=-CBK(1)
        DO 50 K=1,NM
           CDI(K)=CBI(K-1)-K/Z*CBI(K)
50         CDK(K)=-CBK(K-1)-K/Z*CBK(K)
        RETURN
        END
C       ============================================================
C       ============================================================
        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10**(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END
C       ============================================================
C       ============================================================
        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END
C       ============================================================
C       ============================================================
        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
C       ============================================================
C       ============================================================
!        SUBROUTINE WAVE_DIFFRECTION_XFINAS_FORCE
        
      SUBROUTINE DIFFRECTION_WAVE_FORCE(OMEGA,RATIO,H,WATHERDEPTH,RK,RHO,CI,CD,D,COORHORIZON_X,COORHORIZON_Z,ELEVATION
     1                      ,TT,IWAVE,ORDER,VR,AVAL,G,
     1                      VTIDE0,VWIND0,FTX,FTY,FTZ,NC,H0,VCURRENTP,PW,VCL,PERIOD,LWCASE,CS,
     1                      WAVEKINEMATICFACTOR,CURRENTBLOCKAGEFACTOR,WAVEKINEMATICFACTORBREAKING,WAVEBREAKINGKINEMATICFACTOR,
     1                      WAVEFORCECOEFFICIENTS,YLMIN,VELO,ACCE,COA)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        DIMENSION VELO(1),ACCE(1)
        
        !===================================================
        TIME = TT 
        PI = 3.141592653589793D0
        DIA = D
        !WAVELENGTH = G * PERIOD**2 / (2D0*PI )
        
        OMEGA = 2.0*PI/PERIOD      
        CALL NEWTON_RAPHSON(WATHERDEPTH,G,OMEGA,RK)
        WAVELENGTH = 2.0*PI/RK 
        ak = PI * DIA /  WAVELENGTH ! PAGE 393 MECHANICS OF WAVE FORCE ON OFFSHORE STRUCTURE BY TURGUT SARPKAYA
        dk = RK * WATHERDEPTH
        D_OVER_L = DIA / WAVELENGTH
        wavenumber = 2D0 * PI / WAVELENGTH
        a = DIA/2D0
        
        !=========================
        
        !DO I = 1,60 
        !  ELEVATION = I
        
        Call Diffraction_Force_XSEA(ak,ELEVATION,a,WATHERDEPTH,PERIOD,TIME,RHO,H,G,FTX,FTY,
     1  COORHORIZON_X,COORHORIZON_Z,wavenumber,DIA,VELO,ACCE,COA)
        
!        WRITE(790,790)  ELEVATION ,',',FTX 
!790     FORMAT(F12.5,A,E12.5)
        
        !ENDDO
        !CCC = 5
        END SUBROUTINE
C       ============================================================
C       ============================================================
      SUBROUTINE Diffraction_Force_XSEA(ak,z,a,d,period,t,Roh,H,g,Fx,Fy,x,Y_HORIZON,wavenumber,DIA,VELO,ACCE,COA)
        ! ka = ak
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION VELO(3),ACCE(3)
        
        COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV
        
        CALL Diffraction_Force_Parameter(ak,Aka,Cm,Ch,Delta)   
        pi = 3.141592653589793D0
        OMEGA = 2d0*pi/period
        
        cosss = cos(OMEGA*t-wavenumber*x)
        sinnn = sin(OMEGA*t-wavenumber*x)
        
        HYPERBOLIC1 = cosh(wavenumber*z)
        HYPERBOLIC2 = cosh(wavenumber*d) 
         
       !Fx=Roh*g*H*a*Aka/ak*cosh(wavenumber*z)/cosh(wavenumber*d) * cosss 
       
       Fx = pi/8D0*Roh*g*H*wavenumber*(DIA**2)*Cm*cosh(wavenumber*z)/cosh(wavenumber*d) * cosss 
       
       Fy = Roh*g*H*a*Aka/ak*sinh(wavenumber*z)/sinh(wavenumber*d) * sinnn
       
       
       ! REFERENCE KOREA BOOK WAVE THEORY PAGE 150
       
       ! call Diffraction_Pressure_XSEA(m,ak,z,a,d,period,t,Roh,H,g,Fx,x,Y_HORIZON,wavenumber,zeta,DIA)
       
        RETURN
        END SUBROUTINE
C===================================================================
C===================================================================
      SUBROUTINE Diffraction_Pressure_XSEA(m,ak,z,a,d,period,t,Roh,H,g,Fx,x,Y_HORIZON,wavenumber,zeta,DIA,PresureDiff)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        
        complex(8)  Beta 
        complex(8)  Bi
        complex(8)  Coordi
        complex(8)  CHF1
        complex(8)  CHD1
        complex(8)  CHF2
        complex(8)  CHD2
        complex(8)  xxco
        
       ! y = in horizontial 

        xxco = 2 
        
        pi = 3.141592653589793D0
        OMEGA = 2d0*pi/period
        cosss = cos(OMEGA*t-wavenumber*x)
        

        Call Diffraction_Constant_Parameter(ak,x,Y_HORIZON,z,a,DIFFRECTION_PARAMETER,Constant,DIA)
        
        CALL Diffraction_Force_Parameter(ak,Aka,Cm,Ch,Delta)   

        P1 =  Roh*g*H*( z/H )
      
        P2 =  Roh*g*H*( cosh(wavenumber*z)/cosh(wavenumber*d)* DIFFRECTION_PARAMETER ) !* SumP2
      
        PresureDiff =  P1 + P2
        
        PresureDiff =  2d0*pi*g*h*a*Aka/ak*cosh(wavenumber*z)/cosh(wavenumber*d)
        
        PresureDiff =  1.0*P1 + PresureDiff
        
!      WRITE(58,1) z ,  PresureDiff
!1     format(F14.3,4x,F14.3)     
              
          RETURN
          END subroutine
C       ============================================================
C       ============================================================
        SUBROUTINE Diffraction_Constant_Parameter(ak,x,y,z,a,SumP,Constant,DIA)
          
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        
        complex(8)  Beta 
        complex(8)  Bi   
        complex(8)  Coordi
        complex(8)  CHF1  
        complex(8)  CHD1  
        complex(8)  CHF2  
        complex(8)  CHD2  
        complex(8)  ya,xa  
        
        COMMON /DIFFRACTIONCOORDINATE / DIFFRACT_X , DIFFRACT_Y , DIFFRACT_GRAVITY 
        
        pi = 3.141592653589793D0
       
        SumP = 0d0
        
        Do j = 1,50 !50 !100 !10000
        
        m = j
        mj = j
       
        if (m.EQ.0d0) then
        Beta = (1.0 ,0)
        else     
        Bi   = (0 , 1.0 )
        Beta = 2d0*(Bi**m)
        endif


        x = DIFFRACT_X - DIA /2d0
        Y = 0 !DIFFRACT_Y 
        zeta = atan(y/x)
        Coordi = cmplx( x  , DIFFRACT_Y )
        
        IF( y . NE. 0 ) THEN

        call CH12N(1,Coordi,m,CHF1,CHD1,CHF2,CHD2)
        
        ELSEIF(Y . EQ. 0D0) THEN
      

         INDEX_Y = 0D0   
         
       ENDIF
        
        ConH = BJ1 * pi * ak /2d0 
        
        RESULT12 = BESSEL_JN(1, 0.1) 
        
        if ( CHD1 . EQ . 0d0 . OR . INDEX_Y == 0D0    ) then
            
          
        Constant = 0d0    
        else  
        Constant = ( Bi * Beta * cos( mj *zeta) ) / ( pi * ak * CHD1 )
        
        IF( Constant >= 100 ) Constant = 0D0
          
        endif

        
        SumP = SumP + Constant
        
!        write(9600,600) j ,',', Constant ,',', SumP2
!600     format(I5,A,E12.4,A,E12.4)

12        enddo
        
        continue

        RETURN 
      END subroutine
C     ============================================================
C     ============================================================
      SUBROUTINE SHELL_Diffraction_Pressure(OMEGA,RATIO,H,HW,RK,RHO,CI,CD,D,X,Y,TM,IWAVE,ORDER,VR,AVAL,G,
     1                         VTIDE0,VWIND0,PTX,PTY,PTZ,VDUM,NC,H0,VCURRENTP,PW,LC,VCL,PERIOD,LWCASE,CS,
     1                         WAVEKINEMATICFACTOR,CURRENTBLOCKAGEFACTOR,WAVEKINEMATICFACTORBREAKING,WAVEBREAKINGKINEMATICFACTOR,
     1                         WAVEFORCECOEFFICIENTS,YLMIN)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /DIFFRACTIONCOORDINATE / DIFFRACT_X , DIFFRACT_Y , DIFFRACT_GRAVITY 
      
      COMMON /XVWAVE/ VWAVE(3),VWIND(3)
      !XCV = 1
      
        !===================================================
        ELEVATION = DIFFRACT_GRAVITY 
        TIME = TT 
        PI = 3.141592653589793D0
        DIA = D
        WATHERDEPTH = HW
        !WAVELENGTH = G * PERIOD**2 / (2D0*PI )
        
        OMEGA = 2.0*PI/PERIOD      
        CALL NEWTON_RAPHSON(WATHERDEPTH,G,OMEGA,RK)
        WAVELENGTH = 2.0*PI/RK 
         
        ak = PI * DIA /  WAVELENGTH ! PAGE 393 MECHANICS OF WAVE FORCE ON OFFSHORE STRUCTURE BY TURGUT SARPKAYA
        dk = RK * WATHERDEPTH
        D_OVER_L = DIA / WAVELENGTH
        
        wavenumber = 2D0 * PI / WAVELENGTH
        a = DIA/2D0
        !=========================
        
      call Diffraction_Pressure_XSEA(m,ak,ELEVATION,a,d,period,t,RHO,H,g,Fx,DIFFRACT_X,Y_HORIZON,wavenumber,zeta,DIA,
     1 Diffrection_Presure)
            
      PTX = Diffrection_Presure * abs( VWAVE(1) ) * WAVEFORCECOEFFICIENTS
      PTY = Diffrection_Presure * abs( VWAVE(2) ) * WAVEFORCECOEFFICIENTS
      PTZ = Diffrection_Presure * abs( VWAVE(3) ) * WAVEFORCECOEFFICIENTS
      

      
      RETURN
      END SUBROUTINE
C       ============================================================
C       ============================================================
C       ============================================================
C       ============================================================