      SUBROUTINE JACOBIAN_ELIPFUNCTION  
            
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      return
      !CALL CLIPTIC
      
      STOP
      
      DO I = 1,9
      
          AMODULOUS = (I-1)*0.1
      
!      CALL FIRSTKINGELLIPFUNCTION(AMODULOUS,AK,AM)
      
      
!      WRITE(*,*) AMODULOUS , AK ,AM 
      
      ENDDO
      
      
      STOP
      
      RETURN
      
      AMODULOUS = 0.10D0
      

      
!      CALL FIRSTKINGELLIPFUNCTION(AMODULOUS,AK,AM)
      
      
      HK = 0.7D0
      
      STARTVALUE = -5.0D0
      ENDVALUE   =  5.0D0
      AINCREATMENT = 0.01D0
      
      ILOOP = INT(( ENDVALUE - STARTVALUE ) / AINCREATMENT) + 1
      
      U = STARTVALUE
      
!      WRITE(5038,2)                                                                       ! TOEY 10/2021
!2     FORMAT('       m            u0         phi         sn u        cn u        dn u ')  ! TOEY 10/2021
      DO I = 1 , ILOOP
          

      CALL MTHEMATIC_JELP(HK,U,EPH,ESN,ECN,EDN)
      
!      WRITE(5038,1)HK,U,EPH,ESN,ECN,EDN     ! TOEY 10/2021
!1     FORMAT(10F12.5)                       ! TOEY 10/2021
      
      U = U + AINCREATMENT
      
      ENDDO
      
      
      chana = 3
      

      
      
      
      END SUBROUTINE
      
C     ==============================================================
      
      SUBROUTINE MTHEMATIC_JELP(HK,U,EPH,ESN,ECN,EDN)
C
C       ============================================================
C       Purpose: This program computes Jacobian elliptic functions 
C                sn u, cn u and dn u using subroutine JELP
C       Input  : u   --- Argument of Jacobian elliptic fuctions
C                Hk  --- Modulus k ( 0 ?k ?1 )
C       Output : ESN --- sn u
C                ECN --- cn u
C                EDN --- dn u
C                EPH --- phi ( in degrees )
C       Example:
C                k = .5, ( K(k) = 1.68575035 ), and u = u0*K
C
C                u0       phi       sn u        cn u        dn u
C              ----------------------------------------------------
C               0.0      .0000    .0000000   1.0000000   1.0000000
C               0.5    47.0586    .7320508    .6812500    .9306049
C               1.0    90.0000   1.0000000    .0000000    .8660254
C               1.5   132.9414    .7320508   -.6812500    .9306049
C               2.0   180.0000    .0000000  -1.0000000   1.0000000
C       ============================================================
C
      
      
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
        
        !WRITE(*,*)'Please enter k and u '
        ! HK = Modulus k
        !  U = Argument of Jacobian elliptic fuctions
        !READ(*,*)HK,U
!        WRITE(*,*)
!        WRITE(*,*)'   k        u          phi        sn u',
!     &            '        cn u        dn u'
!        WRITE(*,*)' -------------------------------------',
!     &            '---------------------------'
        CALL JELP(U,HK,ESN,ECN,EDN,EPH)
!        WRITE(*,10)HK,U,EPH,ESN,ECN,EDN
10      FORMAT(1X,F5.3,F12.5,2X,F12.5,3F12.7)
        END

C       ============================================================

        SUBROUTINE JELP(U,HK,ESN,ECN,EDN,EPH)
C
C       ========================================================
C       Purpose: Compute Jacobian elliptic functions sn u, cn u
C                and dn u
C       Input  : u   --- Argument of Jacobian elliptic fuctions
C                Hk  --- Modulus k ( 0 ?k ?1 )
C       Output : ESN --- sn u
C                ECN --- cn u
C                EDN --- dn u
C                EPH --- phi ( in degrees )
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION R(40)
        PI=3.14159265358979D0
        A0=1.0D0
        B0=DSQRT(1.0D0-HK*HK)
        DO 10 N=1,40
           A=(A0+B0)/2.0D0
           B=DSQRT(A0*B0)
           C=(A0-B0)/2.0D0
           R(N)=C/A
           IF (C.LT.1.0D-7) GO TO 15
           A0=A
10         B0=B
15      DN=2.0D0**N*A*U
        DO 20 J=N,1,-1
           T=R(J)*DSIN(DN)
           SA=DATAN(T/DSQRT(DABS(1.0D0-T*T)))
           D=.5D0*(DN+SA)
20         DN=D
        EPH=D*180.0D0/PI
        ESN=DSIN(D)
        ECN=DCOS(D)
        EDN=DSQRT(1.0D0-HK*HK*ESN*ESN)
        RETURN
       END
C     ========================================================

      SUBROUTINE FIRSTKINGELLIPFUNCTION(AMODULOUS,AKRESULT,AMRESULT)

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      STARTVAL = 0.0D0
      ENDVAL   = 1.0D0
      AINCREMENT = 0.000010D0
      ILOOP = INT( (ENDVAL-STARTVAL)/AINCREMENT  )
      
      SUMAK = 0
      SUMAM = 0
      
      T = STARTVAL
      
      AK = AMODULOUS
      
      AK = 0.30d0
      
      DO I = 1,ILOOP
          
      T = 0.1
          
      AKKK1 = 1.0D0 / sqrt( ( 1.0D0 - T**2 ) * ( 1.0d0 - ( AK**2 ) *( T**2.0d0 ) ) )
      
      AKKK2 = 1.0D0 / sqrt( ( 1.0D0 - (T+AINCREMENT)**2 ) * ( 1.0d0 - ( AK**2 ) *( (T+AINCREMENT)**2.0d0 ) ) )
      
      AMMM1 = ( sqrt( 1.0d0 - ( AK**2 ) *( T**2.0d0 ) ) )/ ( SQRT( 1 - T**2 ) )
      
      AMMM2 = ( sqrt( 1.0d0 - ( AK**2 ) *( (T+AINCREMENT)**2.0d0 ) ) )/ ( SQRT( 1 - (T+AINCREMENT)**2 ) )
      
      CHANA = 3
      
      AKKK = 0.5D0 * AINCREMENT * ( AKKK1 + AKKK2 )
      AMMM = 0.5D0 * AINCREMENT * ( AMMM1 + AMMM2 )
      
      T = T + AINCREMENT
      
      SUMAK = SUMAK + AKKK * AINCREMENT
      SUMAM = SUMAM + AMMM * AINCREMENT

      ENDDO
      
      AKRESULT = SUMAK
      
      AMRESULT = SUMAM
      
      
      CHANA = 3
      
      
      RETURN
      END SUBROUTINE
C     ========================================================