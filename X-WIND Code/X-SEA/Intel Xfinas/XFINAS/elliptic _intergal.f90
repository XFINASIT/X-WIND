!****************************************************
!*    Program to demonstrate Evaluating elliptic    *
!*  integrals of first and second kinds (complete)  *
!* ------------------------------------------------ *
!* Reference: BASIC Scientific Subroutines, Vol. II *
!* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
!*                                                  *
!*             F90 Version by J.-P. Moreau, Paris.  *
!*                    (www.jpmoreau.fr)             *
!* ------------------------------------------------ *
!* SAMPLE RUN:                                      *
!*                                                  *
!*    K         K(K)          E(K)      STEPS       *
!*  -------------------------------------------     *
!*   0.00     1.5707963     1.5707963     1         *
!*   0.05     1.5717795     1.5698141     3         *
!*   0.10     1.5747456     1.5668619     3         *
!*   0.15     1.5797457     1.5619230     3         *
!*   0.20     1.5868678     1.5549685     3         *
!*   0.25     1.5962422     1.5459573     3         *
!*   0.30     1.6080486     1.5348335     3         *
!*   0.35     1.6225281     1.5215252     3         *
!*   0.40     1.6399999     1.5059416     4         *
!*   0.45     1.6608862     1.4879683     4         *
!*   0.50     1.6857504     1.4674622     4         *
!*   0.55     1.7153545     1.4442435     4         *
!*   0.60     1.7507538     1.4180834     4         *
!*   0.65     1.7934541     1.3886864     4         *
!*   0.70     1.8456940     1.3556611     4         *
!*   0.75     1.9109898     1.3184721     4         *
!*   0.80     1.9953028     1.2763499     4         *
!*   0.85     2.1099355     1.2281083     4         *
!*   0.90     2.2805490     1.1716970     4         *
!*   0.95     2.5900112     1.1027216     5         *
!*   1.00     INFINITY      1.0000000     0         *
!*                                                  *
!****************************************************
      
      SUBROUTINE CLIPTIC
      
      real*8  e,e1,e2,xk
      integer i, n
      
        e=1.d-7
        print *,'  K         K(K)          E(K)      STEPS '
        print *,'------------------------------------------'
        xk=0.d0
        do i = 1, 20
          call CElliptic(e,xk,e1,e2,n)
          write(*,50)  xk,e1,e2,n
          xk = xk + 0.05d0
        end do
        print *,'1.00     INFINITY      1.0000000      0'
        stop
50      format(' ',f4.2,'     ',f9.7,'     ',f9.7,'     ',i2)
        
        RETURN
      END SUBROUTINE


!******************************************************
!* Complete elliptic integral of the first and second *
!* kind. The input parameter is xk, which should be   *
!* between 0 and 1. Technique uses Gauss' formula for *
!* the arithmogeometrical mean. e is a measure of the * 
!* convergence accuracy. The returned values are e1,  *
!* the elliptic integral of the first kind, and e2,   *
!* the elliptic integral of the second kind.          *
!* -------------------------------------------------- *
!* Reference: Ball, algorithms for RPN calculators.   *
!******************************************************
      
      Subroutine CElliptic(e,xk,e1,e2,n)  
      ! Label: 100
      
        real*8 e,xk,e1,e2,pi
        real*8  A(0:99), B(0:99)
        integer j,m,n
        pi = 4.d0*datan(1.d0)
        A(0)=1.d0+xk ; B(0)=1.d0-xk
        n=0
        if (xk < 0.d0) return
        if (xk > 1.d0) return
        if (e <= 0.d0) return
100   n = n + 1
        ! Generate improved values
        A(n)=(A(n-1)+B(n-1))/2.d0
        B(n)=dsqrt(A(n-1)*B(n-1))
        if (dabs(A(n)-B(n)) > e) goto 100
        e1=pi/2.d0/A(n)
        e2=2.d0
        m=1
        do j = 1, n 
          e2=e2-m*(A(j)*A(j)-B(j)*B(j))
          m=m*2
        end do
        e2 = e2*e1/2.d0
        
        RETURN
      END SUBROUTINE


! End of file cliptic.f90

!****************************************************      
      
      
      SUBROUTINE ELLIPTIC_INTERGAL(AMmodulous,Am,AK)
        
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      e=1.d-7
      call CElliptic(e,AMmodulous,Am,AK,n)
      
      RETURN
      END SUBROUTINE

!****************************************************  