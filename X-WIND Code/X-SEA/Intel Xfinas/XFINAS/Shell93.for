
C	================================================================
C	================================================================
C	================================================================
	SUBROUTINE SHAP2D9(R,S,H,P,NODEX,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THERE DERIVATIVES
C     AT THE NODAL POINTS OF A 4 TO 9 NODE,ISOPARAMETRIC QUADRILATERAL
C	----------------------------------------------------------------
C                         NODE NUMBERING CONVENTION
C                         -------------------------
C
C                   2                 5                 1
C
C                     0 . . . . . . . 0 . . . . . . . 0
C                     .                               .
C                     .                               .
C                     .               S               .
C                     .               .               .
C                     .               .               .
C                   6 0             9 0 . . R         0 8
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     0 . . . . . . . 0 . . . . . . . 0
C
C                   3                 7                 4
C
C     R,S        = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     H(9)       = INTERPOLATION (SHAPE) FUNCTIONS
C     P(2,9)     = FUNCTION DERIVATIVES WITH RESPECT TO R,S RESP.
C     NODEX(NEX) = POSITION OF MIDSIDE (EXCESSIVE) NODES
C     NNO        = NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ----------------------------------------------------------------
      DIMENSION  H(9),P(2,9),NODEX(1),IPERM(4)
      DATA (IPERM(I), I=1,4)  /2,3,4,1/
C
      NMI = NNO-4
      RP  = 1.0+R
      SP  = 1.0+S
      RM  = 1.0-R
      SM  = 1.0-S
      R2  = 1.0-R*R
      S2  = 1.0-S*S
C     ---------------------------------------------
C     INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     FOR A FOUR NODE ELEMENT
C     ---------------------------------------------
      H(1)   = 0.25*RP*SP
      H(2)   = 0.25*RM*SP
      H(3)   = 0.25*RM*SM
      H(4)   = 0.25*RP*SM
      P(1,1) = 0.25*SP
      P(1,2) = -P(1,1)
      P(1,3) = -0.25*SM
      P(1,4) = -P(1,3)
      P(2,1) = 0.25*RP
      P(2,2) = 0.25*RM
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)
      IF (NNO.EQ.4)  RETURN
C     -------------------------------------
C     ADD DEGREES OF FREEDOM IN EXCESS OF 4
C     -------------------------------------
      DO 100  IMI=1,NMI
      NN = NODEX(IMI)-4
      GOTO (50,60,70,80,90), NN
C      CALL GOTOER!NON-EXIST SUBROUTINE
 50   H(5)   = 0.5*R2*SP
      P(1,5) = -R*SP
      P(2,5) = 0.5*R2
      GOTO 100
 60   H(6)   = 0.5*RM*S2
      P(1,6) = -0.5*S2
      P(2,6) = -RM*S
      GOTO 100
 70   H(7)   = 0.5*R2*SM
      P(1,7) = -R*SM
      P(2,7) = -0.5*R2
      GOTO 100
 80   H(8)   = 0.5*RP*S2
      P(1,8) = 0.5*S2
      P(2,8) = -RP*S
      GOTO 100
 90   H(9)   =    R2*S2
      P(1,9) =-2.0*R*S2
      P(2,9) =-2.0*S*R2
 100  CONTINUE
C     ----------------------------------------------
C     CORRECT FUNCTIONS AND DERIVATIVES IF 5 OR MORE
C     NODES ARE USED TO DESCRIBE THE ELENENT
C     ----------------------------------------------
	FAC = 0.0D0
	IF(NMI.EQ.5) THEN
	NMI = NMI - 1
	FAC = 1.0D0
	ENDIF

      DO 200  IMI=1,NMI
      IN = NODEX(IMI)
      I1 = IN-4
      I2 = IPERM(I1)
      H(I1)      = H(I1)-0.5*H(IN)+0.25*H(9)*FAC
      H(I2)      = H(I2)-0.5*H(IN)
      H(IMI+4)   = H(IN)          -0.50*H(9)*FAC
      DO 200  J=1,2
      P(J,I1)    = P(J,I1)-0.5*P(J,IN)+0.25*P(J,9)*FAC
      P(J,I2)    = P(J,I2)-0.5*P(J,IN)
 200  P(J,IMI+4) = P(J,IN)            -0.50*P(J,9)*FAC
C
      RETURN
      END
C
C	================================================================
C	================================================================
C	================================================================
	SUBROUTINE SHBMAT9(NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,brr,bss,brs,
     +                brt	,bst,rn,sn)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
c      -----------------------------------------------------------------
c      EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
c
c      NNO               = NUMBER OF NODES FOR ELEMENT
C     H(8)              = SHAPE FUNCTIONS
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     XJI(4)            = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)            = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     -----------------------------------------------------------------
C
      DIMENSION H(9),HD(2,9),VR(3),VS(3),VT(3),XJI(4),HR(9),HS(9)
      DIMENSION B(270),brr(6,54),bss(6,54),brs(4,54),brt(6,54),bst(6,54)
      N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      F1=A2*VR(J)
      F2=A1*VS(J)
      B(N)=F2
      B(N+1)=-F1
      B(N+2)=VS(J)*HS(I)-VR(J)*HR(I)
      N=N+5
   10 CONTINUE
      N=N+15
   20 CONTINUE
      a=sqrt(3.)
      r1 =0.25*(1+a*rn)* sn*(sn+1)
	r2 =0.25*(1-a*rn)* sn*(sn+1)
	r3 = 0.5*(1+a*rn)* (1-sn*sn)
	r4 = 0.5*(1-a*rn)* (1-sn*sn)
	r5 =0.25*(1+a*rn)* sn*(sn-1)
	r6 = 0.25*(1-a*rn)*sn*(sn-1)
c
	s1 =0.25*(1+a*sn)* rn*(rn+1)
	s2 =0.25*(1-a*sn)* rn*(rn+1)
	s3 = 0.5*(1+a*sn)* (1-rn*rn)
	s4 = 0.5*(1-a*sn)* (1-rn*rn)
	s5 =0.25*(1+a*sn)* rn*(rn-1)
	s6 = 0.25*(1-a*sn)*rn*(rn-1)
c	
	rs1 = 0.25*(1+a*rn)*(1+a*sn)
	rs2 = 0.25*(1-a*rn)*(1+a*sn)
	rs3 = 0.25*(1+a*rn)*(1-a*sn)
	rs4 = 0.25*(1-a*rn)*(1-a*sn)
c
	M=1
	L=1
	DO 31 N=1,nno
	DO 30 J=1,3
	brrn=r1*brr(1,l) +r2*brr(2,l) + r3*brr(3,l)
     +    +r4*brr(4,l) +r5*brr(5,l) + r6*brr(6,l)           	 
c
     	bssn=s1*bss(1,l) +s2*bss(2,l) + s3*bss(3,l)
     +    +s4*bss(4,l) +s5*bss(5,l) + s6*bss(6,l)           	 
c     
	brsn=rs1*brs(1,l) +rs2*brs(2,l) +rs3*brs(3,l) + rs4*brs(4,l)
c
	B(M) =brrn*xji(1)*xji(1)+bssn*xji(3)*xji(3)+brsn*2.0*xji(1)*xji(3)
	B(M+1)=brrn*xji(2)*xji(2)+bssn*xji(4)*xji(4)+brsn*2.*xji(2)*xji(4)
	B(M+2)=brrn*xji(1)*xji(2)+bssn*xji(3)*xji(4)+brsn*(xji(1)*xji(4)
     +                                                 +xji(2)*xji(3))
      M=M+5 
      l=l+1
  30	continue
      M=M+15
      l=l+3
  31	continue
      M=1
	do 32 L=1,54
	brtn=r1*brt(1,l) +r2*brt(2,l) + r3*brt(3,l)
     +    +r4*brt(4,l) +r5*brt(5,l) + r6*brt(6,l)           	 
c
     	bstn=s1*bst(1,l) +s2*bst(2,l) + s3*bst(3,l)
     +    +s4*bst(4,l) +s5*bst(5,l) + s6*bst(6,l) 
      B(M+3)=xji(1)*brtn+xji(3)*bstn
	B(M+4)=xji(2)*brtn+xji(4)*bstn
      M=M+5
  32  continue

      RETURN
	END
C	================================================================
C	================================================================
C	================================================================
