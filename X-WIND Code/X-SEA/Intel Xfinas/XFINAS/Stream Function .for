	SUBROUTINE STREAMWAVEFUNCTION (WaveHeight,WavePeriod,WaterDepth,Timedomain,
	1            PeakWater,Gravity,ORDER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION X(50),XF1(61,50),XF2(61,50)
      DIMENSION THETA(61),ETA(61)
C      COMMON /StreamWave / X,XF1,XF2,NN,NTHTS,
C     1                QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,
C     2                OMEGA,time,period,G,WATER_ELEVATION,TOTH,U,V,DUDT,DVDT
     
      COMMON /StreamWave1 / X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,OMEGA
      
      COMMON / StreamInOut / time,period,G,WATER_ELEVATION,TOTH
      ! ----------------------------------------------------------------------------
      ! H                = WAVE HEIGHT           (FT)
      ! T                = WAVE PERIOD           (SEC)
      ! TIME             = TIME                  (SEC)
      ! PI               = PINE (3.141592654)
      ! DPT              = WATER DEPTH
      ! US               = SURFACE CURRENT = 0   (FPS)
      ! UB               = BOTTOM CURRENT  = 0   (FPS)
      ! NORDER           = STREAM FUNCTION ORDER
      ! KMAX             = NUMBER OF ITERATION SET = 100 ( AUTOMATIC )
      ! DAMP             = DAMPING FACTOR
      ! NTHTS            = NUMBER OF SURFACE POINT SET = 61
      ! WATER_ ELEVATION = COORDINATE VERTICAL DIRECTION
      ! -----------------------------------------------------------------------------
      ! ------------------------ METHOD ------------------------
      ! CALL FROM 1) SUBROUTINE FRAMFOC ( FRAMEELEMENT.FOR )
      !           2) SUBROUTINE SOLOFFL ( SELFWEIGHTLOAD.FOR ) 
      !           3) SUBROUTINE SHEOFFL ( SHELLELEMLOAD.FOR )
      ! ALL OF THIS CALL FOR GENERATE BERNOULLI CONSTANT, WAVE CREST, AND TOTAL ELEVATION
      ! USING RELATIVE ERROR TO FORM THIS SOLUTION
      ! 
      ! SECOND, CALL VELO FOR GENERATE VELOCITY AND ACCELERATION
      !
      ! --------------------------           
      H=WaveHeight
      T=WavePeriod
      period=T
      time=Timedomain
      PI=3.141592654
      DPT=WaterDepth
      US=0d0
      UB=0d0 
      NORDER=ORDER
      KMAX=200 
      DAMP=0.0
      NTHTS=61
      G=Gravity
      WATER_ELEVATION = ELEVATION
      ! --------------------------
  888    CONTINUE

      IF(MOD(NTHTS,2).NE.1) THEN
C        WRITE(*,599)
C 599    FORMAT(' NO. OF SURFACE POINTS MUST BE ODD!')
      ENDIF       
      OMEGA=(US-UB)/DPT
C      WRITE(*,733)WVLABEL
C  733 FORMAT('1',9X,' SYMMETRIC STREAM FUNCTION ANALYSIS OF WAVE ',
C     1A20/)
C      OPEN(3,FILE='WAVE.DAT',STATUS='NEW')
C      WRITE(3,431) WVLABEL
  431 FORMAT(A20)
C      WRITE(3,432) H,T,DPT,US,UB,NORDER
  432 FORMAT(5F9.4,I3)
      NN=NORDER+1
      NNP1=NN+1
      IF(KMAX .EQ.0) KMAX=20
      IF( DAMP.EQ.0.0) THEN 
         DAMP=.5
         IF(NN.GT.18) DAMP=.3
      ENDIF
      XN=NTHTS-1
      DTHETA=3.1415927/XN
      THETA(1)=0.0
      DO 10 I=2,NTHTS
   10 THETA(I)=( THETA(I-1)+DTHETA ) 
            
      DBT2=DPT/(T*T)
      IF(DBT2.LT.0.075.AND.DAMP.GT.0.3) DAMP=.3
C*	      DETERMINE BREAKING INDEX
      HBT2T=.873*SINH(.8935*DBT2)/COSH(.8935*DBT2)
C  101 FORMAT(///' WAVE HEIGHT = ',F6.2,' FEET, PERIOD = ',F6.2,' SECONDS
C     1'/' WATER DEPTH = ',F6.2,' FEET ')
C      WRITE(*,105) KMAX,NORDER,UB,US
      HBT2=H/(T*T)
      TT=HBT2/HBT2T
      IF(TT.GT.0.75.AND.DAMP.GT.0.3) DAMP=.3
C      IF((HBT2-.01-HBT2T).GT.0.0) GO TO 6
C      IF(ABS(HBT2-HBT2T)-.01) 3,3,5
C    3 WRITE(*,153)
C      WRITE(*,201)
C  153 FORMAT(' *********************************************************
C     1**********************')
C  201 FORMAT(' **** INPUT WAVE PARAMETERS ARE WITHIN THE ACCURACY OF THE
C     1 BREAKING INDICATOR ***','**** IF SOLUTION IS NOT SUCCESSFUL,
C     2CHECK INPUT PARAMETERS ****' )
C      WRITE(*,153)
      GO TO 5
C    6 WRITE(*,153)
      HB=HBT2T*T*T
C      WRITE(*,200) HB
C  200 FORMAT(' **** WAVE HEIGHT EXCEEDS BREAKING WAVE HEIGHT OF APPROXIM
C     1ATELY',F8.2,' FEET ****'/'    CHECK INPUT PARAMETERS.'   )
C      WRITE(*,153)
    5 CONTINUE
   35 SIG2G=((2.*PI/T)**2)/G
      XK=SIG2G/SQRT(TANH(SIG2G*DPT))
      X(1)=2.*PI/XK 
      X(1)=X(1)+.2*X(1)*(HBT2/HBT2T)**2
      X(2)=-X(1)*H/(T*2.0*SINH(2.0*PI*(DPT+H/2.0)/X(1)))
      DO 132 N=3,NN
  132 X(N)=0.0
      X(NNP1)=0.0
C      WRITE(*,553) (X(K),K=1,NN)
C  553 FORMAT(//' INITIAL STREAM FUNCTION COEFFICIENTS:'/(1P5E12.3))
      CALL STREAM(KMAX,DAMP,*35,*888,IAUTOFF,G)
C  105 FORMAT(' TOTAL NUMBER OF ITERATIONS SPECIFIED= ',I3,', WAVE THEORY
C     1 ORDER= ',I3/' BOTTOM CURRENT VELOCITY= ',F6.2,' FPS., SURFACE  CU
C     1RRENT VELOCITY= ',F6.2,' FPS.')
C       CALL VELOC (U,V,DUDT,DVDT)
C      WRITE(3,433)NTHTS
C  433 FORMAT(I4)
C      WRITE(3,434) (X(IJ),IJ=1,NNP1)
C  434 FORMAT(5E12.5)
C      CLOSE (3)
  800 CONTINUE
      
      !ETA = ETA*1.05d0
      ! --- GENERATE PEAK WATER LEVEL ----
      DO I=1,1 !NTHTS,NTHTS-1
      TOTH=ETA(I)
      ENDDO
      ! ----------------------------------
      
      !CALL STREAMWAVESURFACE(TIMECAL)
      
      !PEAKWATER=TOTH
      
      
      TIMECAL = TIME
      
      
      
      TIMEROTIO = (TIMECAL/T * 60.0D0) + 1
      
      
      NTIME = INT(TIMEROTIO)
      
      RATIO = TIMEROTIO - NTIME
      
      CHANA =  3
      
      WS1 = ETA(NTIME)
      WS2 = ETA(NTIME+1)
      
      WAVES = ETA(NTIME) + (ETA(NTIME)-ETA(NTIME+1)) * RATIO
      
      PEAKWATER = WAVES * 1.03d0
      
      CHANA =  3
      
C      STOP
      
      END
C*------------------------------------------------------------
      SUBROUTINE STREAM(KMAX,DAMP,*,*,IAUTOFF,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION X(50),XF1(61,50),XF2(61,50)
      DIMENSION THETA(61),ETA(61)
      !COMMON X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UU,THETA,ETA,DTHETA
      
C      COMMON /StreamWave / X,XF1,XF2,NN,NTHTS,
C     1                QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,
C     2                OMEGA,time,period,G,WATER_ELEVATION,TOTH,U,V,DUDT,DVDT
     
      COMMON /StreamWave1 / X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,OMEGA

      IEND=0
      DO 5 I=1,NTHTS
    5 ETA(I)=0.0
      DO 6 I=1,NTHTS
      DO 6 N=2,NN
      XN=N-1
      XF1(I,N)=COS(THETA(I)*XN)
    6 XF2(I,N)=SIN(THETA(I)*XN)
      KMAXP1=KMAX+1
      DO 20  KI=1,KMAXP1
      L=1
      IF(KI.EQ.KMAXP1) L=2
      CALL FSCALC(L,ERRORH)
      
      ! ----- AUTOMATIC ITERAION ---
      IF (KI.GE.2.0)THEN
      IF (ERRORH.LE.0.1)THEN
      GOTO 1000 ! END FUNCTION
      ENDIF
      ENDIF
      ! ----------------------------
      
      IF(KI.EQ.KMAXP1) IEND=1
      IKI=KI
      CALL CFF(DAMP,IEND,AMSL,ERRORQ,IKI,G)
      TOTERR=ABS(ERRORH)+ABS(AMSL)+ABS(ERRORQ)
C     WRITE(*,70)ERRORH,AMSL,ERRORQ,TOTERR
C  70 FORMAT(1P10E12.3)
      IF (IAUTOFF.EQ.1) GO TO 21
      IF(TOTERR.LE..01) RETURN
   21 IF(KI.EQ.1) E1=TOTERR
      IF(TOTERR.GT.E1) GO TO 25
   20 CONTINUE
      RETURN
   25 DAMP=DAMP/2.
C      IF(DAMP.LE.0.0) GO TO 30
C      WRITE(*,110)
C  110 FORMAT('          SOLUTION DIVERGING---PROGRAM IS RESTARTING WITH
C     1DAMPING')
      RETURN 1
C   30 WRITE(*,120)
C  120 FORMAT(//' *****************************************'/' **********
C     1CONVERGENCE IMPOSSIBLE *****'/' *********************************
C     2*************'/'                      CHECK INPUT VARIABLES   ')
      RETURN 2
1000  END
C*------------------------------------------------------------
      SUBROUTINE FSCALC(L,ERRORH)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION X(50),XF1(61,50),XF2(61,50)
      DIMENSION THETA(61),ETA(61),ETA1(61)
      !COMMON X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UU,THETA,ETA,DTHETA,H,OMEGA
!      COMMON /StreamWave / X,XF1,XF2,NN,NTHTS,
!     1                QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,
!     2                OMEGA,time,period,G,WATER_ELEVATION,TOTH,U,V,DUDT,DVDT
     
      COMMON /StreamWave1 / X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,OMEGA
     
     
      IEND=0
      J=0
      DO 30 I=1,NTHTS
      CALL FSI(I,ETA(I),N,IEND,OMEGA)
      IF(N.GT.J)J=N
   30 CONTINUE
      IF(IEND.EQ.0) GO TO 10
C      WRITE(*,101)
C  101 FORMAT(I6,' 20 ITERATIONS ON ETA FAILED TO CONVERGE')
   10 CONTINUE
      XH=ETA(1)-ETA(NTHTS)
      ERRORH=H-XH
C      WRITE(*,102) XH,ERRORH
C  102 FORMAT('        WAVE HEIGHT =',F8.2,' FT., ERROR =',1PE10.3,' FT.'
C     1)
      IF(L.EQ.1) GO TO 12
C      WRITE(*,105)(THETA(I),ETA(I),I=1,NTHTS)
C  105 FORMAT('        THETAS AND CALCULATED ETAS '/(6(1X,F8.2,1X
C     1,F8.4)))
   12 RETURN
      END
C*------------------------------------------------------------
      SUBROUTINE CFF(DAMP,IEND,AMSL,ERRORQ,KI,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION X(50),XF1(61,50),XF2(61,50)
      DIMENSION ETA(61),DCDX(53),DETADX(61,53),DUDX(53),DVDX(53),
     1DQDX(61,53),Q(61),B(2809),D(53),ERROR(61),THETA(61),DQBAR(53)
      !COMMON X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UO,THETA,ETA,DTHETA,H,OMEGA

C      COMMON /StreamWave / X,XF1,XF2,NN,NTHTS,
C     1                QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,
C     2                OMEGA,time,period,G,WATER_ELEVATION,TOTH,U,V,DUDT,DVDT
     
      COMMON /StreamWave1 / X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,OMEGA

      CON=6.2831853/X(1)
C      G=32.17 !===============
      NTHTM1=NTHTS-1
      QBAR=0.0
      NNP1=NN+1
      NNP3=NN+3
      DI=DTHETA/(3.0*3.1415927)
      XNTHTS=NTHTS
      DCDX(1)=1.0/T
      DO 1 N=2,NNP1
    1 DCDX(N)=0.0
      DO 20 I=1,NTHTS
      U=0.0
      V=0.0
      DUDX(1)=0.0
      DVDX(1)=0.0
      DVDETA=0.0
      DUDETA=0.0
      DETADX(I,1)=ETA(I)/T
      TOTH=DPT+ETA(I)
      UU=UO
      IF(OMEGA.NE.0.0) UU=UO+OMEGA*TOTH
      DO 14 N=2,NN
      XN=N-1
      ZN=CON*XN
      COEF1=XF1(I,N)*X(N)
      COEF2=XF2(I,N)*X(N)
      SH=SINH(ZN*TOTH)
      CH=COSH(ZN*TOTH)
      U=U-ZN*CH*COEF1
      V=V-ZN*COEF2*SH
      DETADX(I,1)=DETADX(I,1)-ZN*TOTH/X(1)*CH*COEF1
      DVDX(N)=-ZN*SH*XF2(I,N)
      DVDETA=DVDETA-ZN*ZN*CH*COEF2
      DETADX(I,N)=SH*XF1(I,N)
      DUDX(1)=DUDX(1)+(ZN*ZN*TOTH/X(1)*SH+ZN/X(1)*CH)*COEF1
      DUDX(N)=-ZN*CH*XF1(I,N)
      DUDETA=DUDETA-ZN*ZN*SH*COEF1
   14 DVDX(1)=DVDX(1)+(ZN*ZN*TOTH/X(1)*CH+ZN/X(1)*SH)*COEF2
      CC=U+UU-X(1)/T
      DQDETA=1.0
      DETADX(I,NNP1)=-1.0
      DUDX(NNP1)=0.0
      DVDX(NNP1)=0.0
      DQDU=CC/G
      DQDC=-DQDU
      DQDV=V/G
      DO 15 N=1,NNP1
   15 DETADX(I,N)=DETADX(I,N)/CC
      DO 16 N=1,NNP1
   16 DQDX(I,N)=DQDETA*DETADX(I,N)+ DQDU*(DUDX(N)+DUDETA*DETADX(I,N))
     1+DQDV*(DVDX(N)+DVDETA*DETADX(I,N))+DQDC*DCDX(N)
   20 Q(I)=ETA(I)+1.0/(2.0*G)*(CC*CC+V*V)
      DO 31 I=2,NTHTM1,2
   31 QBAR=QBAR+Q(I-1)+4.0*Q(I)+Q(I+1)
      QBAR=QBAR*DI
      DO 22 N=1,NNP1
      DQBAR(N)=0.0
      DO 21 I=2,NTHTM1,2
   21 DQBAR(N)=DQBAR(N)+DQDX(I-1,N)+4.0*DQDX(I,N)+DQDX(I+1,N)
   22 DQBAR(N)=DQBAR(N)*DI
      ERRORQ=0.0
      DO 23 I=1,NTHTS
   23 ERROR(I)=Q(I)-QBAR
      DO 25 I=2,NTHTM1,2
   25 ERRORQ=ERRORQ+ERROR(I-1)**2+4.0*ERROR(I)**2+ERROR(I+1)**2
      ERRORQ=SQRT(ERRORQ*DI)
C      WRITE(*,19) QBAR
C   19 FORMAT(' BERNOULLI CONSTANT = ',F7.2,' FT.')
C      WRITE(*,24) ERRORQ
C   24 FORMAT(' RMS ERROR IN DFSBC = ',1PE9.3,' FT.')
C     WRITE(*,124)(ERROR(I),I=1,NTHTS)
C 124 FORMAT('          THE DISTRIBUTED DFSBC ERRORS:'/(1P10E11.3))
      NNP2=NN+2
      K=0
C*	      DEFINE AN NNP3 X NNP1 ARRAY FOR B(K)
      DO 30 N=1,NNP1
      DO 30 M=1,NNP3
      K=K+1
      B(K)=0.0
      IF(M.EQ.NNP2) GO TO 50
      IF(M.EQ.NNP3) GO TO 170
      DO 26 I=2,NTHTM1,2
   26 B(K)=B(K)+(DQDX(I-1,N)-DQBAR(N))*(DQDX(I-1,M)-DQBAR(M))+4.0*(DQDX
     1(I,N)-DQBAR(N))*(DQDX(I,M)-DQBAR(M))+(DQDX(I+1,N)-DQBAR(N))*(DQDX
     2(I+1,M)-DQBAR(M))
      B(K)=B(K)*DTHETA/3.0
      GO TO 30
   50 DO 51 I=2,NTHTM1,2
   51 B(K)=B(K)+DETADX(I-1,N)+4.0*DETADX(I,N)+DETADX(I+1,N)
      B(K)=B(K)*DTHETA/3.0
      GO TO 30
  170 B(K)=DETADX(1,N)-DETADX(NTHTS,N)
   30 CONTINUE
C*	      LAMBDA1 COEFFICIENTS
      DO 55 M=1,NNP3
      D(M)=0.0
      K=K+1
      B(K)=0.0
      IF(M.EQ.NNP2) GO TO 56
      IF(M.EQ.NNP3) GO TO 172
      DO 53 I=2,NTHTM1,2
   53 B(K)=B(K)+DETADX(I-1,M)+4.0*DETADX(I,M)+DETADX(I+1,M)
      B(K)=B(K)*DTHETA/3.0
      DO 57 I=2,NTHTM1,2
   57 D(M)=D(M)+(QBAR-Q(I-1))*(DQDX(I-1,M)-DQBAR(M))+4.0*(QBAR-Q(I))*(DQ
     1DX(I,M)-DQBAR(M))+(QBAR-Q(I+1))*(DQDX(I+1,M)-DQBAR(M))
      D(M)=D(M)*DTHETA/3.0
      GO TO 55
   56 DO 58 I=2,NTHTM1,2
   58 D(M)=D(M)+ETA(I-1)+4.0*ETA(I)+ETA(I+1)
      D(M)=D(M)*DTHETA/3.0
      AMSL=D(M)/3.1415927
C      WRITE(*,96)AMSL
C   96 FORMAT('          AMSL ERROR = ',1PE9.3,' FT.')
      IF(IEND.EQ.1) RETURN
      D(M)=-D(M)
      GO TO 55
  172 D(M)=H+ETA(NTHTS)-ETA(1)
   55 CONTINUE
C*	      LAMBDA2 COEFFICIENTS
      DO 175 M=1,NNP3
      K=K+1
      B(K)=0.0
      IF(M.GE.NNP2) GO TO 175
      B(K)=DETADX(1,M)-DETADX(NTHTS,M)
  175 CONTINUE
      CALL SIMQ(B,D,NNP3,KS)
C      IF(KS.EQ.1) WRITE(*,87)
C   87 FORMAT('          CFF MATRIX IS SINGULAR')
      DO 40 N=1,NNP1
   40 X(N)=X(N)+D(N)*DAMP
C      WRITE(*,92) KI
C   92 FORMAT(///'           ITERATION',I3)
C      WRITE(*,91)(X(N),N=1,NNP1)
C   91 FORMAT('           X(N)=',1P5E12.4)
      RETURN
      END
C*------------------------------------------------------------
      SUBROUTINE FUNC(ET,Y,YP,I,OMEGA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION X(50),XF1(61,50),XF2(61,50),ETA(61),ETA1(61)
      !COMMON X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UU
      
!      COMMON /StreamWave / X,XF1,XF2,NN,NTHTS,
!     1                QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,
!     2                OMEGA,time,period,G,WATER_ELEVATION,TOTH,U,V,DUDT,DVDT
     
      COMMON /StreamWave1 / X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UB,THETA,ETA,DTHETA,H !,OMEGA
     
     
      CON=6.2831853/X(1)
      NNP1=NN+1
      C=X(1)/T-UU
      ELEV=DPT+ET
      Y=X(NNP1)-C*ET+OMEGA*ELEV*ELEV/2.
      YP=-C+OMEGA*ELEV
      DO 14 N=2,NN
      ZN=CON*FLOAT(N-1)
      Y=Y-X(N)*SINH(ZN*ELEV)*XF1(I,N)
   14 YP=YP-X(N)*COSH(ZN*ELEV)*XF1(I,N)*ZN
   
      RETURN
      END
C*------------------------------------------------------------
      SUBROUTINE FSI(I,X,N,IEND,OMEGA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      M=20
      EPSIL=0.001
      XPP=X
      CALL FUNC(XPP,YPP,YDPP,I,OMEGA)
      XP=XPP-YPP/YDPP
      DO 21 N=1,M
      CALL FUNC(XP,YP,YDP,I,OMEGA)
      D=XP-XPP
      A=(2.*YDP+YDPP-3.*(YP-YPP)/D)/D
      U=YP/YDP
      X=XP-U*(1.+U*A/YDP)
      IF(ABS((X-XP)/X)-EPSIL)3,3,4
    4 XPP=XP
      XP=X
      YPP=YP
   21 YDPP=YDP
      IEND=1
    3 RETURN
      END
C*------------------------------------------------------------
      SUBROUTINE SIMQ(A,B,N,KS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION A(2809),B(53)
C*
C*		 FORWARD SOLUTION
C*
      TOL=0.0
      KS=0
      JJ=-N
      DO 65 J=1,N
      JY=J+1
      JJ=JJ+N+1
      BIGA=0
      IT=JJ-J
      DO 30 I=J,N
C*
C*		 SEARCH FOR MAXIMUM COEFFICIENT IN COLUMN
C*
      IJ=IT+I
      IF(ABS(BIGA)-ABS(A(IJ))) 20,30,30
   20 BIGA=A(IJ)
      IMAX=I
   30 CONTINUE
C*
C*		 TEST FOR PIVOT LESS THAN TOLERANCE (SINGULAR MATRIX)
C*
      IF(ABS(BIGA)-TOL) 35,35,40
   35 KS=1
      RETURN
C*
C*		 INTERCHANGE ROWS IF NECESSARY
C*
   40 I1=J+N*(J-2)
      IT=IMAX-J
      DO 50 K=J,N
      I1=I1+N
      I2=I1+IT
      SAVE=A(I1)
      A(I1)=A(I2)
      A(I2)=SAVE
C*
C*		 DIVIDE EQUATION BY LEADING COEFFICIENT
C*
   50 A(I1)=A(I1)/BIGA
      SAVE=B(IMAX)
      B(IMAX)=B(J)
      B(J)=SAVE/BIGA
C*
C*		 ELIMINATE NEXT VARIABLE
C*
      IF(J-N) 55,70,55
   55 IQS=N*(J-1)
      DO 65 IX=JY,N
      IXJ=IQS+IX
      IT=J-IX
      DO 60 JX=JY,N
      IXJX=N*(JX-1)+IX
      JJX=IXJX+IT
   60 A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX))
   65 B(IX)=B(IX)-(B(J)*A(IXJ))
C*
C*		 BACK SOLUTION
C*
   70 NY=N-1
      IT=N*N
      DO 80 J=1,NY
      IA=IT-J
      IB=N-J
      IC=N
      DO 80 K=1,J
      B(IB)=B(IB)-A(IA)*B(IC)
      IA=IA-N
   80 IC=IC-1
      RETURN
      END
C*-----------------------------------------------------------
      SUBROUTINE VELOC (XC,YC,U,V,DUDT,DVDT,TT,RAMDA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION X(50),XF1(61,50),XF2(61,50)
      DIMENSION THETA(61),ETA(61)
      INTEGER AUTOFF
      CHARACTER*20 WVLABEL
      
      COMMON /StreamWave1 / X,XF1,XF2,NN,NTHTS,QBAR,T,DPT,UB,THETA,ETA,DTHETA,H,OMEGA
      COMMON / StreamInOut / time,period,G,WATER_ELEVATION,TOTH
      
      T=TT
      CON=6.2831853/X(1)
      rALAMDA=X(1)
      RAMDA=X(1)
      
      DO 200 I=1,1 !NTHTS,NTHTS-1
      TOTH=DPT+ETA(I)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !write(*,*)'TOTH',TOTH
      !ELEV=TOTH+DY
      !=====================
      !ELEV=10.88
      U=0.0
      V=0.0
      DUDX=0.0
      DUDZ=0.0
      !ELEV=ELEV-DY
      !==================
      ELEV=YC
      UU=UO
      IF(OMEGA.NE.0.0) UU=UO+OMEGA*ELEV
      DO 14 N=2,NN
      
      Xcoordinate=XC
      Ycoordinate=YC
      
      XN=N-1
      ZN=CON*XN
      !COEF1=XF1(I,N)*X(N)
      !COEF2=XF2(I,N)*X(N)
      COEF1=X(N)!*XF1(I,N)
      COEF2=X(N)!*XF2(I,N)
      SH=SINH(ZN*Ycoordinate)
      CH=COSH(ZN*Ycoordinate)
      
      !t=time
      !xxx=cos(XN*2d0*3.141592654*t/period)
      !! OMEAGA=2D0*3.141592654*t/period
      !U=U-ZN*CH*COEF1 *cos(XN*2d0*3.141592654*t/period)
      !V=V-ZN*COEF2*SH *sin(XN*2d0*3.141592654*t/period)
      
      PI=3.141592654
      TIME=TT
      COSfunction=COS(XN*2d0*PI*Xcoordinate/rALAMDA-XN*2D0*PI*TIME/period)
      Sinfunction=SIN(XN*2d0*PI*Xcoordinate/rALAMDA-XN*2D0*PI*TIME/period)
      U=U-ZN*CH*COEF1*COSfunction
      V=V-ZN*COEF2*SH*Sinfunction
      
      DUDX=DUDX+ZN*ZN*COEF2*CH * Sinfunction
   14 DUDZ=DUDZ-ZN*ZN*COEF1*SH * COSfunction
      CC=X(1)/period
      DUDT=(U+UU-CC)*DUDX+V*DUDZ
      DVDT=(U+UU-CC)*DUDZ-V*DUDX
C      WRITE(6,99)ELEV, U+UU, DUDT ,V, DVDT
C   99 FORMAT(5F8.2)
  200 CONTINUE 
      RETURN
      END
C*-----------------------------------------------------------      

