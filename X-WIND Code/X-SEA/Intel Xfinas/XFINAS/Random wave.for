      Subroutine  Pierson_Moskowitz_RANDOM (OPT,ATIME,OFFSELECT,SURFACE,X,Y,UAX,UAY,AUX,AUY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /RANDOM_NUMBER_STORE/ RADNUM(1000,100)
      !COMMON/RANDOM_ELV/ EVL_RAN(20000),ITIME,IOFFL
      COMMON /SEA_PARA/ NSPECTRUM
      CHARACTER*3 OPT
      
      REWIND (4904)
      IF (OPT.EQ."SUR") THEN
      SURFACE = 0.
      WN      = 0.
      AN      = 0.
      FP      = 0.
      SPM     = 0.
      C       = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP

      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      !CALL RANDOM_NUMBER(ANUMBR_RANDOM)
      !ANUMBR_RANDOM = ANUMBR_RANDOM*2D0*PI
      
      ! ----- validation random function from matlab ------
      !READ (4904,*) ANUMBR_RANDOM
      !RADNUM(I,1)     = ANUMBR_RANDOM
      
      JJ = OFFSELECT
      ANUMBR_RANDOM  = RADNUM(I,JJ)
      C   = WN*TP/(2D0*PI) 
      !SPM = 0.04974*(HS**2D0)*TP*(C**-5D0)*EXP(-1.25D0*C**-4D0)
      SPM=0.3125d0*(HS**2)*(FP**4D0)*(WN**-5d0)*exp(-1.25D0*(FP/WN)**4D0)
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 100
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      !RK = WN**2D0/9.806D0 ! VALIDATION DATA
      SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      !SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X)
      !UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      !AUX     = AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
100   ENDDO
           !IF (IOFFL.EQ.1) THEN
           !EVL_RAN(ITIME) = SURFACE
           !ENDIF
      RETURN
      ELSEIF (OPT.EQ."DAI") THEN
      WN  = 0.
      FP  = 0.
      SPM = 0.
      C   = 0.
      UAX = 0.
      UAY = 0.
      AUX = 0.
      AUY = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP
       
      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      !CALL RANDOM_NUMBER(ANUMBR_RANDOM)
      ANUMBR_RANDOM = RADNUM(I,OFFSELECT)
      C   = WN*TP/(2D0*PI) 
      !SPM = 0.04974*(HS**2D0)*TP*(C**-5D0)*EXP(-1.25D0*C**-4D0)
      SPM=0.3125d0*(HS**2)*(FP**4D0)*(WN**-5d0)*exp(-1.25D0*(FP/WN)**4D0)
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 200
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      !OMEGA = 2.0*PI/(1D0/(WN/2D0*PI))
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      WAVE_LENGTH = 2D0*PI/RK
      !SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      UAY     = UAY + AN*WN*SINH(RK*Y)/(SINH(RK*WDEPTH))*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUX     = AUX + AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUY     = AUY - AN*(WN**2)*SINH(RK*Y)/SINH(RK*WDEPTH)*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      ! --- validation mathlab---
      !RK      = WN**2D0/9.806D0
      !AAA     = AN*WN*(COSH(RK*20)/(SINH(RK*20)))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      !UAX     = UAX + AAA
200   ENDDO
      RETURN
      ENDIF  
      
      END
C   ==========================================================================================================          
      Subroutine  Jonswap_RANDOM (OPT,ATIME,OFFSELECT,SURFACE,X,Y,UAX,UAY,AUX,AUY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /RANDOM_NUMBER_STORE/ RADNUM(1000,100)
      CHARACTER*3 OPT
      
      IF (OPT.EQ."SUR") THEN
      SURFACE = 0.
      WN      = 0.
      AN      = 0.
      FP      = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP

      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      ANUMBR_RANDOM  = RADNUM(I,OFFSELECT)
      !  ---- GRAMMA PARAMETER ---
        IF(Tp/sqrt(Hs).LE.3.6d0)THEN
        Gamma=5D0
        ELSEIF(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) THEN
        Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
        ELSEIF(Tp/sqrt(Hs).GT.5.0d0)THEN
        Gamma=1d0
        ENDIF
      ! ---------------------------
      
      ! ---- ST PARAMETER ---
        IF (WN.LE.Fp) THEN
        St=0.07d0
        ELSE
        St=0.09d0
        ENDIF
      ! ---------------------------  
      SPM=0.3125d0*(Hs**2)*(1/Fp)*((WN/Fp)**-5d0)*exp(-1.25*(Fp/WN)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((WN/Fp-1)/st)**2)
      
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 100
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
100   ENDDO
      RETURN
      ELSEIF (OPT.EQ."DAI") THEN
      WN = 0.
      FP  = 0.
      UAX = 0.
      UAY = 0.
      AUX = 0.
      AUY = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP
       
      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      ANUMBR_RANDOM  = RADNUM(I,OFFSELECT)
      !  ---- GRAMMA PARAMETER ---
        IF(Tp/sqrt(Hs).LE.3.6d0)THEN
        Gamma=5D0
        ELSEIF(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) THEN
        Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
        ELSEIF(Tp/sqrt(Hs).GT.5.0d0)THEN
        Gamma=1d0
        ENDIF
      ! ---------------------------
      
      ! ---- ST PARAMETER ---
        IF (WN.LE.Fp) THEN
        St=0.07d0
        ELSE
        St=0.09d0
        ENDIF
      ! ---------------------------  
      SPM=0.3125d0*(Hs**2)*(1/Fp)*((WN/Fp)**-5d0)*exp(-1.25*(Fp/WN)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((WN/Fp-1)/st)**2)
      
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 200
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      !OMEGA = 2.0*PI/(1D0/(WN/2D0*PI))
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      !SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      UAY     = UAY + AN*WN*SINH(RK*Y)/(SINH(RK*WDEPTH))*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUX     = AUX + AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUY     = AUY - AN*(WN**2)*SINH(RK*Y)/SINH(RK*WDEPTH)*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
200   ENDDO
      RETURN
      ENDIF  
      END
      
C   ==========================================================================================================     
      Subroutine Ochi_RANDOM (OPT,ATIME,OFFSELECT,SURFACE,X,Y,UAX,UAY,AUX,AUY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /RANDOM_NUMBER_STORE/ RADNUM(1000,100)
      CHARACTER*3 OPT
      
      IF (OPT.EQ."SUR") THEN
      SURFACE = 0.
      WN      = 0.
      AN      = 0.
      FP      = 0.
      SOH1 = 0.
      SOH2 = 0.
      SPM  = 0.
      CALL SELECTRANDOM_OCHI (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH,Hs1,Hs2,TP1,TP2,AShape1,AShape2)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      Call GAMMA_FUNCTION(AShape1,GAM1)
      Call GAMMA_FUNCTION(AShape2,GAM2)
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP

      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      ANUMBR_RANDOM  = RADNUM(I,OFFSELECT)
        SOH1 = (Hs1**2)/(4d0*GAM1)*(((4d0*AShape1+1d0)/4d0*OMEGAM1**4d0)**AShape1)*(1d0/(WN**(4d0*AShape1+1d0)))
     1       *EXP(((4d0*AShape1+1d0)/-4d0)*(OMEGAM1/WN)**4d0)
        SOH2 = (Hs2**2)/(4d0*GAM2)*(((4d0*AShape2+1d0)/4d0*OMEGAM2**4d0)**AShape2)*(1d0/(WN**(4d0*AShape2+1d0)))
     1       *EXP(((4d0*AShape2+1d0)/-4d0)*(OMEGAM2/WN)**4d0)
        SPM = SOH1+SOH2
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 100
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
100   ENDDO
      RETURN
      ELSEIF (OPT.EQ."DAI") THEN
      WN = 0.
      FP  = 0.
      UAX = 0.
      UAY = 0.
      AUX = 0.
      AUY = 0.
      SOH1 = 0.
      SOH2 = 0.
      SPM  = 0.
      CALL SELECTRANDOM_OCHI (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH,Hs1,Hs2,TP1,TP2,AShape1,AShape2)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      Call GAMMA_FUNCTION(AShape1,GAM1)
      Call GAMMA_FUNCTION(AShape2,GAM2)
      OMEGAM1=2d0*Pi/TP1
      OMEGAM2=2d0*Pi/TP2
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP
       
      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      ANUMBR_RANDOM  = RADNUM(I,OFFSELECT)
        SOH1 = (Hs1**2)/(4d0*GAM1)*(((4d0*AShape1+1d0)/4d0*OMEGAM1**4d0)**AShape1)*(1d0/(WN**(4d0*AShape1+1d0)))
     1       *EXP(((4d0*AShape1+1d0)/-4d0)*(OMEGAM1/WN)**4d0)
        SOH2 = (Hs2**2)/(4d0*GAM2)*(((4d0*AShape2+1d0)/4d0*OMEGAM2**4d0)**AShape2)*(1d0/(WN**(4d0*AShape2+1d0)))
     1       *EXP(((4d0*AShape2+1d0)/-4d0)*(OMEGAM2/WN)**4d0)
        SPM = SOH1+SOH2
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 200
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      !OMEGA = 2.0*PI/(1D0/(WN/2D0*PI))
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      !SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      UAY     = UAY + AN*WN*SINH(RK*Y)/(SINH(RK*WDEPTH))*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUX     = AUX + AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUY     = AUY - AN*(WN**2)*SINH(RK*Y)/SINH(RK*WDEPTH)*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
200   ENDDO
      RETURN
      ENDIF  
      END

C   ==========================================================================================================     
      Subroutine  Bretschneider_RANDOM (OPT,ATIME,OFFSELECT,SURFACE,X,Y,UAX,UAY,AUX,AUY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /RANDOM_NUMBER_STORE/ RADNUM(1000,100)
      CHARACTER*3 OPT
      
      ! SIGNIFICANT WAVE HEIGHT
      IF (OPT.EQ."SUR") THEN
      SURFACE = 0.
      WN      = 0.
      AN      = 0.
      FP      = 0.
      SPM     = 0.
      C       = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP

      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      !CALL RANDOM_NUMBER(ANUMBR_RANDOM)
      !ANUMBR_RANDOM = ANUMBR_RANDOM*2D0*PI
      !READ (4904,*) ANUMBR_RANDOM
      !RADNUM(I)     = ANUMBR_RANDOM
      JJ = OFFSELECT
      ANUMBR_RANDOM  = RADNUM(I,JJ)
      C   = WN*TP/(2D0*PI) 
      !SPM = 0.04974*(HS**2D0)*TP*(C**-5D0)*EXP(-1.25D0*C**-4D0)
      SPM=0.3125d0*(HS**2)*(FP**4D0)*(WN**-5d0)*exp(-1.25D0*(FP/WN)**4D0)
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 100
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      !RK = WN**2D0/9.806D0 ! VALIDATION DATA
      SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      !UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      !AUX     = AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
100   ENDDO

      RETURN
      ELSEIF (OPT.EQ."DAI") THEN
      WN  = 0.
      FP  = 0.
      SPM = 0.
      C   = 0.
      UAX = 0.
      UAY = 0.
      AUX = 0.
      AUY = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP
       
      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      !CALL RANDOM_NUMBER(ANUMBR_RANDOM)
      ANUMBR_RANDOM = RADNUM(I,OFFSELECT)
      C   = WN*TP/(2D0*PI) 
      !SPM = 0.04974*(HS**2D0)*TP*(C**-5D0)*EXP(-1.25D0*C**-4D0)
      SPM=0.3125d0*(HS**2)*(FP**4D0)*(WN**-5d0)*exp(-1.25D0*(FP/WN)**4D0)
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 200
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      !OMEGA = 2.0*PI/(1D0/(WN/2D0*PI))
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK)
      !SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      UAY     = UAY + AN*WN*SINH(RK*Y)/(SINH(RK*WDEPTH))*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUX     = AUX + AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUY     = AUY - AN*(WN**2)*SINH(RK*Y)/SINH(RK*WDEPTH)*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
200   ENDDO
      RETURN
      ENDIF  
      
      END
C   ==========================================================================================================     
      Subroutine  TMA_SPectrum_RANDOM (OPT,ATIME,OFFSELECT,SURFACE,X,Y,UAX,UAY,AUX,AUY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /RANDOM_NUMBER_STORE/ RADNUM(1000,100)
      CHARACTER*3 OPT
      
            ! SIGNIFICANT WAVE HEIGHT
      IF (OPT.EQ."SUR") THEN
      SURFACE = 0.
      WN      = 0.
      AN      = 0.
      FP      = 0.
      SPM     = 0.
      C       = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP

      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      !CALL RANDOM_NUMBER(ANUMBR_RANDOM)
      !ANUMBR_RANDOM = ANUMBR_RANDOM*2D0*PI
      !READ (4904,*) ANUMBR_RANDOM
      !RADNUM(I)     = ANUMBR_RANDOM
      JJ = OFFSELECT
      ANUMBR_RANDOM  = RADNUM(I,JJ)
      
      Alpha=exp(-((WN-Fp)**2)/(2d0*St**2d0*Fp**2))
      IF(Tp/sqrt(Hs).LE.3.6d0)THEN
      Gamma=5
      ELSEIF(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) THEN
      Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
      ELSEIF(Tp/sqrt(Hs).GT.5.0d0)THEN
      Gamma=1d0
      ENDIF
       
      IF (WN.LE.Fp) THEN
      ST=0.07d0
      ELSE
      ST=0.09d0
      ENDIF
       
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK) 
      SJS=0.3125d0*(Hs**2)*(1/Fp)*((WN/Fp)**-5d0)*exp(-1.25*(Fp/WN)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((WN/Fp-1)/ST)**2)
      Phi=((cosh(RK*WDEPTH))**2d0)/((sinh(RK*WDEPTH))**2d0+ ((WN)**2)*WDEPTH/GRAV)
      SPM=SJS*Phi
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 100
      ENDIF
      
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      !RK = WN**2D0/9.806D0 ! VALIDATION DATA
      SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      !UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      !AUX     = AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
100   ENDDO

      RETURN
      ELSEIF (OPT.EQ."DAI") THEN
      WN  = 0.
      FP  = 0.
      SPM = 0.
      C   = 0.
      UAX = 0.
      UAY = 0.
      AUX = 0.
      AUY = 0.
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
      
      PI  = 3.141592654D0
      FP  = 2D0*PI/TP
       
      DO I = 1,FREQ_INTER+1
      IF (I.EQ.1) WN = FREQ_START*2D0*PI
      IF (I.NE.1) WN = WN + FREQ_INTEVAL*2D0*PI
      !CALL RANDOM_NUMBER(ANUMBR_RANDOM)
      ANUMBR_RANDOM = RADNUM(I,OFFSELECT)
      Alpha=exp(-((WN-Fp)**2)/(2d0*St**2d0*Fp**2))
      IF(Tp/sqrt(Hs).LE.3.6d0)THEN
      Gamma=5
      ELSEIF(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) THEN
      Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
      ELSEIF(Tp/sqrt(Hs).GT.5.0d0)THEN
      Gamma=1d0
      ENDIF
       
      IF (WN.LE.Fp) THEN
      ST=0.07d0
      ELSE
      ST=0.09d0
      ENDIF
       
      CALL Airy_Wave_Number(GRAV,WN,WDEPTH,RK) 
      SJS=0.3125d0*(Hs**2)*(1/Fp)*((WN/Fp)**-5d0)*exp(-1.25*(Fp/WN)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((WN/Fp-1)/ST)**2)
      Phi=((cosh(RK*WDEPTH))**2d0)/((sinh(RK*WDEPTH))**2d0+ ((WN)**2)*WDEPTH/GRAV)
      SPM=SJS*Phi
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 100
      ENDIF
      
      IF (WN.EQ.0) THEN
      SPM = 0.0D0
      GOTO 200
      ENDIF
      AN  = SQRT(2D0*SPM*FREQ_INTEVAL*2D0*PI)
      !SURFACE = SURFACE + AN*SIN(WN*ATIME-RK*X+ANUMBR_RANDOM)
      UAX     = UAX + AN*WN*COSH(RK*Y)/(SINH(RK*WDEPTH))*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
      UAY     = UAY + AN*WN*SINH(RK*Y)/(SINH(RK*WDEPTH))*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUX     = AUX + AN*(WN**2)*COSH(RK*Y)/SINH(RK*WDEPTH)*SIN(RK*X - WN*ATIME + ANUMBR_RANDOM)
      AUY     = AUY - AN*(WN**2)*SINH(RK*Y)/SINH(RK*WDEPTH)*COS(RK*X - WN*ATIME + ANUMBR_RANDOM)
200   ENDDO
      RETURN
      ENDIF  
      END
C   ==========================================================================================================
      SUBROUTINE RANDOM_DATA 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /RANDOM_NUMBER_STORE/ RADNUM(1000,100)
      COMMON /SEA_PARA/ NSPECTRUM
      PI  = 3.141592654D0
      DO J = 1,NSPECTRUM    
       OFFSELECT = J   
       !REWIND (4904) ! TESTING VALUE MATLAB
      CALL SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,WDEPTH)
      FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_INTER)
         DO I = 1,FREQ_INTER+1
         CALL RANDOM_NUMBER(ANUMBR_RANDOM)
         ANUMBR_RANDOM = ANUMBR_RANDOM*2D0*PI
         !READ (4904,*) ANUMBR_RANDOM ! TESTING VALUE MATLAB
         RADNUM(I,J)     = ANUMBR_RANDOM
         ENDDO
      ENDDO
      
      
      END
C   ==========================================================================================================      
      SUBROUTINE SELECTRANDOM (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,AWDEPTH)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION VWAVE(3)
      DIMENSION FWIND(3)
      DIMENSION VG(3),V1(3),V2(3),VWIND(3)
      COMMON/ RES_OFFSHOREX / SEABEDRX(1000),WVHIGHTRX(1000),WDEPTHRX(1000),H1POSRX(1000),H2POSRX(1000),
     1                        NRESWAVERX(1000),PERIODRX(1000),GRAVRX(1000),RHOWRX(1000),RHOARX(1000),WKFRX(1000),WFCRX(1000),
     1                        FREQENCYX(1000),DAMPINGX(1000),
     1                        WVHIGHTRX1(1000),WVHIGHTRX2(1000),PERIODRX1(1000),PERIODRX2(1000),SHAPEX1(1000),SHAPEX2(1000),
     1                        TAPX(1000),UCURRENTX(1000),
     1                        NSWINDX(1000),UMEANWINDSPECX(1000),POWERLAWSPECX(1000),ROUGHSPECX(1000),REFERSPECX(1000),
     1                        TURFACTORX(1000),SCALEFACTORX(1000),STDFACTORX(1000),
     1                        POINT(1000,1000),SFREQ(1000,1000),AMPLI(1000,1000),
     1                        VGX(1000,3),V1X(1000,3),V2X(1000,3),VWX(1000,3),
     1                        FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000)
C    1                        NOPSPECX(1000),TIMESTARTX(1000),TIMEINX(1000),TIMEENDX(1000),
      
      AWDEPTH       = WDEPTHRX(OFFSELECT)
      GRAV          = GRAVRX(OFFSELECT)
      FREQ_START    = FSpectrumStartX(OFFSELECT)
      FREQ_END      = FSpectrumEndX(OFFSELECT)
      FREQ_INTER    = NUMBEROFRANDOMWANVEXInterationX(OFFSELECT)
      HS            = WVHIGHTRX(OFFSELECT)
      TP            = PERIODRX(OFFSELECT)
      
      END
C   ==========================================================================================================      
      SUBROUTINE SELECTRANDOM_OCHI (OFFSELECT,FREQ_START,FREQ_END,FREQ_INTER,GRAV,HS,TP,AWDEPTH,WH1,WH2,TP1,TP2,SHAP1,SHAP2)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION VWAVE(3)
      DIMENSION FWIND(3)
      DIMENSION VG(3),V1(3),V2(3),VWIND(3)
      COMMON/ RES_OFFSHOREX / SEABEDRX(1000),WVHIGHTRX(1000),WDEPTHRX(1000),H1POSRX(1000),H2POSRX(1000),
     1                        NRESWAVERX(1000),PERIODRX(1000),GRAVRX(1000),RHOWRX(1000),RHOARX(1000),WKFRX(1000),WFCRX(1000),
     1                        FREQENCYX(1000),DAMPINGX(1000),
     1                        WVHIGHTRX1(1000),WVHIGHTRX2(1000),PERIODRX1(1000),PERIODRX2(1000),SHAPEX1(1000),SHAPEX2(1000),
     1                        TAPX(1000),UCURRENTX(1000),
     1                        NSWINDX(1000),UMEANWINDSPECX(1000),POWERLAWSPECX(1000),ROUGHSPECX(1000),REFERSPECX(1000),
     1                        TURFACTORX(1000),SCALEFACTORX(1000),STDFACTORX(1000),
     1                        POINT(1000,1000),SFREQ(1000,1000),AMPLI(1000,1000),
     1                        VGX(1000,3),V1X(1000,3),V2X(1000,3),VWX(1000,3),
     1                        FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000)
C    1                        NOPSPECX(1000),TIMESTARTX(1000),TIMEINX(1000),TIMEENDX(1000),
      
      AWDEPTH       = WDEPTHRX(OFFSELECT)
      GRAV          = GRAVRX(OFFSELECT)
      FREQ_START    = FSpectrumStartX(OFFSELECT)
      FREQ_END      = FSpectrumEndX(OFFSELECT)
      FREQ_INTER    = NUMBEROFRANDOMWANVEXInterationX(OFFSELECT)
      HS            = WVHIGHTRX(OFFSELECT)
      TP            = PERIODRX(OFFSELECT)
      WH1           = WVHIGHTRX1(OFFSELECT)
      WH2           = WVHIGHTRX2(OFFSELECT)
      TP1           = PERIODRX1(OFFSELECT)
      TP2           = PERIODRX2(OFFSELECT)
      SHAP1         = SHAPEX1(OFFSELECT)
      SHAP2         = SHAPEX2(OFFSELECT)
      END

