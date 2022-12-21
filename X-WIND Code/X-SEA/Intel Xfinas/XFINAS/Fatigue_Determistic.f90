      SUBROUTINE  STATICS_SCF
      
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV
      COMMON / STATIC_JOINT_SCF / STATICSCF
      COMMON / JOINT_FIXITY / JOINTFIXITY
      COMMON /INPUTJOINT/ NNODEJOINT(1000),JOINTCASE(1000),SAFTYFACTOR(1000),DESIGNLIFE(1000),GAP_JOINT(1000),JOINTLOOP
      COMMON /TURBULAR_GAP/ GAP
      
!      KOFFL   = OFFSHORE 
!      KSPEC   = SPECTRAL
!      KDESIGN = STEEL DESIGN
!      KFATM   = FATIGUE MEMBER
!      KFATJ   = FATIGUE JOINT
!      KFATL   = FATIGUE LIFE
!      KFAST   = FAST OPTION
!      KOREV   = RELATIVE FLOW ( MORISON EQ. )
      !========================================
      ! INNITIAL SETUP
      !========================================
      
      ! INNITIAL SETUP
      !KFATJ = 1
      
      JOINTFIXITY = 1
      
      ! INPUT
      STATICSCF = 1
      
      JOINTLOOP = 1
      GAP_JOINT(1) = 1
      
      GAP = 0.1D0
      
      DO I=1,JOINTLOOP
          
          NNODEJOINT(I) = 3
          GAP_JOINT(I)  = -8.28D0*0.01D0
          
      ENDDO
      
!       WRITE(5000,1)
!1      FORMAT('                              |--------------- AXIAL ------------------|---IN PLANE-------|----OUT PLANE------|')       
!       WRITE(5000,2)
!2      FORMAT('  JOINT   NUM_BR   NUM_CH TYPE SADDLE_CH  CROWN_CH SADDLE_BR  CROWN_BR  CROWN_CH  CROWN_BR SADDLE_CH SADDLE_BR')    
      
      WRITE(5000,3)
3     FORMAT('  JOINT   TYPE  NUM_BR   AX-SD-CH  AX-CR-CH  IN-PL-CH  OU-PL-CH 
     1  NUM_CH  AX-SD_BR  AX-CR_BR  IN-PL-BR  OU-PL-BR')
      
      write(5001,4)
4     FORMAT('  JOINT    NELE   MEMBER  NPLANE   AXIAL_XX    SHEAR_TT    SHEAR_SS   MOMENT_TO   MOMENT_TT   MOMENT_SS ')
      
      write(5002,5)
5     FORMAT('  JOINT    NELE  MEMBER  NPLANE   STRESS_XX   STRESS_TT   STRESS_SS   STRESS_TO   STRESS_TT   STRESS_SS ')
      
      CALL FATIGUE_CONNECT (NNODE,ITIME,NFATIGUELOOPOUT)
      

      END
      !==========================================================================================================================
      Subroutine SCF_X_JOINT_AXIAL (Betaparameter,Shearparameter,Gammaparameter
     1                ,Alpha,ZETA,NOPTION,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      !C
      !C1
      !C2
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
      C = 0.7D0  
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 )  
            
      !Chord saddle
      EQ12 = 3.87d0 * Gammaparameter * Shearparameter * Betaparameter *(1.10d0 - Betaparameter **1.8d0)*(sin(zeta))**1.7d0
      
      
      CHANA = (sin(zeta))**1.7d0
      
      SCFCHORDSADDLE = EQ12
      
      IF( Alpha .LT. 12 ) SCFCHORDSADDLE = EQ12 * F3
      
      
      !Chord crown
      EQ13 = ( Gammaparameter**0.2 ) *  Shearparameter * (2.65d0+5d0*(Betaparameter-0.65d0)**2d0)  
     & - 3d0 * Shearparameter * Betaparameter * sin(zeta)
      
      SCFCHORDCROWN = EQ13
      
      !Brace saddle
      EQ14 = 1d0 + 1.9d0 * Gammaparameter * (Shearparameter **0.5d0) * (Betaparameter**0.90d0) *
     & (1.09d0 - Betaparameter**1.7d0) *  (sin(zeta))**2.5d0
      
      SCFBRACESADDLE = EQ14
      
       IF( Alpha .LT. 12 ) SCFBRACESADDLE = EQ14 * F3
      
      !Brace crown
      EQ15 = 3.0d0 +(Gammaparameter**1.2d0)* (0.12d0 * (exp(-4d0*Betaparameter)) + 0.011*(Betaparameter**2d0) -0.045d0 ) 
      
       SCFBRACECROWN = EQ15
      
      ! In joints with short cords ( Alpha < 12 ) the saddle SCF can be reduced by the factor F1 (fixed chord ends) of F2 ( pinned chord ends)
      ! where

      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F2 = 1d0 - (1.43d0 * Betaparameter - 0.97d0 *  Betaparameter**2d0 - 0.03d0 ) *  Gammaparameter**(0.04d0) 
     & *exp( -0.71d0 * ( Gammaparameter**-1.38d0 )  * Alpha**2.5d0 )
      
      RETURN
      End  subroutine
      !==========================================================================================================================
      
      Subroutine SCF_X_JJOINT_IN_PLANE (Betaparameter,Shearparameter,Gammaparameter
     1                ,Alpha,ZETA,NOPTION,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      !C
      !C1
      !C2
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
      C = 0.7D0  
      
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 ) 
      
      !IF (NOPTION.EQ.1) THEN ! IN PLANE BENDING
      !Chord crown
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0 
      
      SCFCHORDCROWN = EQ8
      
      !Brace crown
      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
      
      SCFBRACECROWN = EQ9
      
      !ELSEIF (NOPTION.EQ.2)THEN ! OUT PLANE BENDING
      !Chord saddle
      EQ16 =  Gammaparameter *  Shearparameter  * Betaparameter* 
     & (1.56d0 - 1.34d0*Betaparameter**4.0d0) *  (sin(zeta))**1.6d0
      
      SCFCHORDSADDLE = EQ16
      
      IF( Alpha .LT. 12 )  SCFCHORDSADDLE = EQ16 * F3
      
      !Brace saddle
      
      EQ17 = Shearparameter**(-0.54d0)*( Gammaparameter **(-0.05) )*( 0.99d0 - 0.47 * Betaparameter + 0.08d0* Betaparameter**4d0)
     & * EQ16
      
      SCFBRACESADDLE = EQ17
      
      IF( Alpha .LT. 12 ) SCFBRACESADDLE = EQ17 * F3
      
      ! In joints with short chords ( Alpha < 12 ) Equations (16) and (17) can be reduced by the factor F3 where:
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 )  
      
      ! ENDIF
      Return
      End  Subroutine
      !==========================================================================================================================
      !==========================================================================================================================
      
      Subroutine SCF_X_JJOINT_OUT_PLANE (Betaparameter,Shearparameter,Gammaparameter
     1                ,Alpha,ZETA,NOPTION,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      !C
      !C1
      !C2
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
      C = 0.7D0  
      
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 ) 
      
      !IF (NOPTION.EQ.1) THEN ! IN PLANE BENDING
      !Chord crown
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0 
      
      SCFCHORDCROWN = EQ8
      
      !Brace crown
      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
      
      SCFBRACECROWN = EQ9
      
      !ELSEIF (NOPTION.EQ.2)THEN ! OUT PLANE BENDING
      !Chord saddle
      EQ16 =  Gammaparameter *  Shearparameter  * Betaparameter* 
     & (1.56d0 - 1.34d0*Betaparameter**4.0d0) *  (sin(zeta))**1.6d0
      
      SCFCHORDSADDLE = EQ16
      
      IF( Alpha .LT. 12 )  SCFCHORDSADDLE = EQ16 * F3
      
      !Brace saddle
      
      EQ17 = Shearparameter**(-0.54d0)*( Gammaparameter **(-0.05) )*( 0.99d0 - 0.47 * Betaparameter + 0.08d0* Betaparameter**4d0)
     & * EQ16
      
      SCFBRACESADDLE = EQ17
      
      IF( Alpha .LT. 12 ) SCFBRACESADDLE = EQ17 * F3
      
      ! In joints with short chords ( Alpha < 12 ) Equations (16) and (17) can be reduced by the factor F3 where:
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 )  
      
      ! ENDIF
      Return
      End  Subroutine
      !==========================================================================================================================
      !==========================================================================================================================
      Subroutine SCF_KT_JOINT_AXIAL (Betaparameter,Shearparameter,Gammaparameter,ZINE
     1  ,DUMMYGammaparameter,dummyShearparameter,zetaA 
     1 ,Alpha,ZETA,DUMMYZETA,IMAXA,DUMMYBETA,IMAXB,NOPTION,SCFCHORD,SCFBRACE,
     1   SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWNM,IOPTION,KT_rankX,Zine_II)
      !===============================================
      ! According to DNV
      !SCF_Joint_fomular_fixity_Condigion_B5part1
      !===============================================
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION DUMMYZETA(IMAXA),DUMMYBETA(IMAXB),DUMMYGammaparameter(IMAXA),dummyShearparameter(IMAXA)
      DIMENSION ANGLE(IMAXA)
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      ! zine
      !C
      !C1
      !C2
      ! GET ZETA MAX
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
      C = 0.7D0 
      
      
      
      zetamax            = MAXVAL(DUMMYZETA)
      zetamin            = MINVAL(DUMMYZETA)
      AmaxBetaparameter  = MAXVAL(DUMMYBETA)
      AminBetaparameter  = MINVAL(DUMMYBETA)
      
      gap = 0.1d0
      
      zi_AB = gap
      
      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F2 = 1d0 - (1.43d0 * Betaparameter - 0.97d0 *  Betaparameter**2d0 - 0.03d0 ) *  Gammaparameter**(0.04d0) 
     & *exp( -0.71d0 * ( Gammaparameter**-1.38d0 )  * Alpha**2.5d0 )
      
      !=======================================================
      ! Balanced axial load
      !=======================================================
      !IF (NOPTION.EQ.1)THEN
      !===================    
      !Chord:
      ! K1 : API Page 211
      !===================    
          
          
       
      if ( KT_rankX .EQ. 2 )    Then   ! Center Brace
                
      EQ20 = ( Shearparameter ** 0.9d0 ) * (Gammaparameter**0.5d0) * ( 0.67d0 - Betaparameter**2d0 + 1.16 * Betaparameter )
     &  * sin(zeta) *( ( sin(zetamax)/sin(zetamin))**0.30d0 ) *( (AmaxBetaparameter/AminBetaparameter)**0.3d0 )
     &  *(1.64d0 + 0.29d0 * (Betaparameter**(-0.38d0)) * atan(8d0*( Zine_II )))
      
      chana = (1.64d0 + 2.9d0 * (Betaparameter**-0.38d0) * atan(8d0*( 0.949 )))
      
      SCFCHORD = EQ20 !! Eq: K1 : API Page 211 
      
      if ( Alpha .LT. 12 ) SCFCHORD = EQ20 * F1
      
      
      SCFBRACE = 1d0 + SCFCHORD *(1.97d0-1.57d0*(Betaparameter**0.25d0) ) * ( Shearparameter ** -0.14d0 ) *( (sin(zeta))**0.30d0 )
      
      chana =3
      

       
      
      elseif ( KT_rankX .NE. 2 )    then
           
      EQ20 = ( Shearparameter ** 0.9d0 ) * (Gammaparameter**0.5d0) * ( 0.67d0 - Betaparameter**2d0 + 1.16 * Betaparameter )
     &  * sin(zeta) *( ( sin(zetamax)/sin(zetamin))**0.30d0 ) *( (AmaxBetaparameter/AminBetaparameter)**0.3d0 )
     &  *(1.64d0 + 0.29d0 * (Betaparameter**(-0.38d0)) * atan(8d0*( Zine_II )))
      
      
      XXXXchana = (1.64d0 + 2.9d0  )
      
      SCFCHORD = EQ20 !! Eq: K1 : API Page 211   
      
      !if ( Alpha .LT. 12 ) SCFCHORD = EQ20 * F1
      
      SCFBRACE = 1d0 + SCFCHORD *(1.97d0-1.57d0*(Betaparameter**0.25d0) ) * ( Shearparameter ** -0.14d0 ) *( (sin(zeta))**0.70d0 )
      
      
      chana =3
      
      endif
      
      return
      End subroutine
      !==========================================================================================================================
      !==========================================================================================================================
      Subroutine SCF_KT_JOINT_IN_PLANE (Betaparameter,Shearparameter,Gammaparameter,ZINE
     1  ,DUMMYGammaparameter,dummyShearparameter,zetaA 
     1 ,Alpha,ZETA,DUMMYZETA,IMAXA,DUMMYBETA,IMAXB,NOPTION,SCFCHORD,SCFBRACE,
     1   SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWNM,IOPTION,KT_rankX,Zine_A,Zine_B,Zine_C)
      !===============================================
      ! According to DNV
      !SCF_Joint_fomular_fixity_Condigion_B5part1
      !===============================================
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION DUMMYZETA(IMAXA),DUMMYBETA(IMAXB),DUMMYGammaparameter(IMAXA),dummyShearparameter(IMAXA)
      DIMENSION ANGLE(IMAXA)
      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      ! zine
      !C
      !C1
      !C2
      ! GET ZETA MAX
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
      C = 0.7D0 
      
      
      
      zetamax            = MAXVAL(DUMMYZETA)
      zetamax_degree     = zetamax * 180/3.1416d0
      zetamin            = MINVAL(DUMMYZETA)
      zetamin_degree     = zetamin * 180/3.1416d0
      AmaxBetaparameter  = MAXVAL(DUMMYBETA)
      AminBetaparameter  = MINVAL(DUMMYBETA)
      
      gap = 0.1d0
      
      zi_AB = gap
      
      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F2 = 1d0 - (1.43d0 * Betaparameter - 0.97d0 *  Betaparameter**2d0 - 0.03d0 ) *  Gammaparameter**(0.04d0) 
     & *exp( -0.71d0 * ( Gammaparameter**-1.38d0 )  * Alpha**2.5d0 )
      
      !=======================================================
      ! Balanced axial load
      !=======================================================
      !IF (NOPTION.EQ.1)THEN
      !===================    
      !Chord:
      ! K1 : API Page 211
      !===================    
          
          IF ( Zine_A .LE. Zine_C )  Zine_MIN_CAL =  Zine_A
          IF ( Zine_A .GT. Zine_C )  Zine_MIN_CAL =  Zine_C
       
      if ( KT_rankX .EQ. 2 )    Then   ! Center Brace
                
       ! IN plane

      
      EQ8 = 1.45d0*Betaparameter*(Shearparameter**(0.85d0) )*( Gammaparameter**(1.0d0 - 0.68d0*Betaparameter))*(sin(zeta))**0.7d0
!     & * (1.0d0 * 0.46d0 * (Betaparameter**12.0d0) * exp( -3.0d0 * ( Zine_II ) ))
      
      
      xchana = (1.0d0 * 0.46d0 * (Betaparameter**12.0d0) * exp( -3.0d0 * ( Zine_MIN_CAL ) ))
      
      SCFCHORDCROWN = EQ8
      
      !Brace crown:
       EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
       
       SCFBRACECROWNM = EQ9
      
      elseif ( KT_rankX .NE. 2 )    then
           
      EQ8 = 1.45d0*Betaparameter*(Shearparameter**(0.85d0) )*( Gammaparameter**(1.0d0 - 0.68d0*Betaparameter))*(sin(zeta))**0.7d0
     & * (1.0d0 + 0.46d0 * (Betaparameter**12.0d0) * exp( -3.0d0 * ( Zine_B ) ))
      
      
      xchana = (1.0d0 + 0.46d0 * (Betaparameter**12.0d0) * exp( -3.0d0 * ( Zine_B ) ))
      
      SCFCHORDCROWN = EQ8
      
      !Brace crown:
       EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
       
      SCFBRACECROWNM = EQ9
      
      chana = 3

      Endif
      
      Return
      End subroutine
      !==========================================================================================================================
      !==========================================================================================================================
      !==========================================================================================================================
      Subroutine SCF_KT_JOINT_OUT_PLANE (Betaparameter,Shearparameter,Gammaparameter,ZINE
     1  ,DUMMYGammaparameter,dummyShearparameter,zetaA 
     1 ,Alpha,ZETA,DUMMYZETA,IMAXA,DUMMYBETA,IMAXB,NOPTION,SCFCHORD,SCFBRACE,
     1   SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWNM,IOPTION,KT_rankX,Zine_A,Zine_B,Zine_C,XX_AB, XX_BC,XX_AC)
      !===============================================
      ! According to DNV
      !SCF_Joint_fomular_fixity_Condigion_B5part1
      !===============================================
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION DUMMYZETA(IMAXA),DUMMYBETA(IMAXB),DUMMYGammaparameter(IMAXA),dummyShearparameter(IMAXA)
      DIMENSION ANGLE(IMAXA)
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      ! zine
      !C
      !C1
      !C2
      ! GET ZETA MAX
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
      C = 0.7D0 
      
      
      
      zetamax            = MAXVAL(DUMMYZETA)
      zetamin            = MINVAL(DUMMYZETA)
      AmaxBetaparameter  = MAXVAL(DUMMYBETA)
      AminBetaparameter  = MINVAL(DUMMYBETA)
      
      gap = 0.1d0
      
      zi_AB = gap
      
      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F2 = 1d0 - (1.43d0 * Betaparameter - 0.97d0 *  Betaparameter**2d0 - 0.03d0 ) *  Gammaparameter**(0.04d0) 
     & *exp( -0.71d0 * ( Gammaparameter**-1.38d0 )  * Alpha**2.5d0 )
      
      !=======================================================
      ! Balanced axial load
      !=======================================================
      !IF (NOPTION.EQ.1)THEN
      !===================    
      !Chord:
      ! K1 : API Page 211
      !===================    
          
          
       
      If ( KT_rankX .EQ. 2 )    Then   ! Center Brace
                
      SCF_Chord_center = ( Gammaparameter*Shearparameter*Betaparameter*(1.7-1.05*(Betaparameter**3.0d0)* ((sin(ZETA) )**1.6)) )*
     1 (1-0.08* ( ( Betaparameter * Gammaparameter )**0.5d0 ) * exp(-0.8d0* XX_BC ))
      
      
       SCFCHORDSADDLE  = SCF_Chord_center  !! Eq: K1 : API Page 211 
       
       SCF_Brace_center = (Shearparameter**-0.54d0) * (Gammaparameter**-0.05) *(0.99d0-0.47d0*Betaparameter+0.08*(Betaparameter**4)) 
     1  * SCF_Chord_center 
       
       SCFBRACESADDLE = SCF_Brace_center

      chana = 3
       
      Elseif ( KT_rankX .NE. 2 )    then
          
      SCF_Chord_center = ( Gammaparameter*Shearparameter*Betaparameter*(1.7-1.05*(Betaparameter**3.0d0)) * (( sin(ZETA) )**1.6d0) )
     1 *(1-0.08* ( ( Betaparameter * Gammaparameter )**0.5d0 ) * exp(-0.8d0* 3.762 )) 
      
       SCFCHORDSADDLE  = SCF_Chord_center  !! Eq: K1 : API Page 211  
       
       SCF_Brace_center = (Shearparameter**-0.54d0) * (Gammaparameter**-0.05) *(0.99d0-0.47d0*Betaparameter+0.08*(Betaparameter**4)) 
     1  * SCF_Chord_center  
       
       SCFBRACESADDLE = SCF_Brace_center

      chana = 3
      
      Endif
      
      return
      End subroutine
      !==========================================================================================================================
          subroutine SCF_BALANCE_K_JOINT_AXIAL (SCFChord,SCFBrace,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,ZINE
     1                                            ,SCFBRACECROWN,Betaparameter,Shearparameter,dummyShearparameter,Gammaparameter
     1                                                  ,DUMMYGammaparameter,Alpha,ZETA,DUMMYZETA,IMAXA,DUMMYBETA,IMAXB,NOPTION)
       
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION DUMMYZETA(IMAXA),DUMMYBETA(IMAXB),DUMMYGammaparameter(IMAXA),dummyShearparameter(IMAXA)
      DIMENSION ANGLE(IMAXA)

      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      ! zine
      ! C
      ! C1
      ! C2
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
      C = 0.7D0   
      
      zetamax            = MAXVAL(DUMMYZETA)
      zetamin            = MINVAL(DUMMYZETA)
      AmaxBetaparameter  = MAXVAL(DUMMYBETA)
      AminBetaparameter  = MINVAL(DUMMYBETA)
      
      !=============================================
      IF (NOPTION.EQ.1)THEN ! BALANCED AXIAL CHORD
      !=============================================
          
      !zine = -0.078
         
      C = 0.7
          
      ! Chord
      EQ20 = ( Shearparameter ** 0.9d0 ) * (Gammaparameter**0.5d0) * ( 0.67d0 - Betaparameter**2d0 + 1.16 * Betaparameter )
     &  * sin(zeta) *( ( sin(zetamax)/sin(zetamin))**0.30d0 ) *( (AmaxBetaparameter/AminBetaparameter)**0.3d0 )
     &  *(1.64d0 + 0.29d0 * (Betaparameter**-0.38d0) * atan(8d0*( zine )))
            
      SCFChord = EQ20
      
      
      
      ! Brace:
      EQ21 = 1d0 + (1.97d0 - 1.57d0 *Betaparameter**0.25d0 )*(Shearparameter**-0.14)*(sin(zeta)**0.7d0)*EQ20 +
     & (sin(zetamax + zetamin)**1.8d0)*(0.131d0-0.084d0*atan(14d0*zine + 4.2d0 * Betaparameter)) 
     &  * C * (Betaparameter**1.5d0) * (Gammaparameter**0.5d0) * (Shearparameter**-1.22d0)
        
      SCFBrace = EQ21

      
      ENDIF
      

      RETURN
      
      END SUBROUTINE

      !==========================================================================================================================
      !==========================================================================================================================
      !==========================================================================================================================
      
      subroutine SCF_BALANCE_K_JOINT_INPANE (SCFChord,SCFBrace,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,ZINE
     1                                            ,SCFBRACECROWN,Betaparameter,Shearparameter,dummyShearparameter,Gammaparameter
     1                                                  ,DUMMYGammaparameter,Alpha,ZETA,DUMMYZETA,IMAXA,DUMMYBETA,IMAXB,NOPTION)
                
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION DUMMYZETA(IMAXA),DUMMYBETA(IMAXB),DUMMYGammaparameter(IMAXA),dummyShearparameter(IMAXA)
      DIMENSION ANGLE(IMAXA)

      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      ! zine
      !C
      !C1
      !C2
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
            
      C = 0.7D0      
      
      
      zetamax            = MAXVAL(DUMMYZETA)
      zetamin            = MINVAL(DUMMYZETA)
      AmaxBetaparameter  = MAXVAL(DUMMYBETA)
      AminBetaparameter  = MINVAL(DUMMYBETA)
      
      F4 = 1.0D0 - 1.07D0*(Betaparameter**1.88)*EXP(0.16D0*(Gammaparameter**(-1.06D0))*(Alpha**(2.4D0)))
      
      !=============================================
      !IF (NOPTION.EQ.1)THEN ! BALANCED AXIAL CHORD
      !=============================================
 
            
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0 
     & *(1D0+0.46D0*(Betaparameter**1.2D0)*EXP(-3.0D0*zine))
      
      SCFCHORDCROWN = EQ8
      
      IF(Alpha.LT.12) SCFCHORDCROWN = ABS(EQ8*F4)
      
      
      
      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
      
      SCFBRACECROWN = EQ9
      
      !==================================================
      !ELSEIF (NOPTION.EQ.2)THEN ! UNBALANCED AXIAL CHORD
      !==================================================
          

      !ENDIF
      
      
      RETURN
      END SUBROUTINE
      
      !==========================================================================================================================     
      !==========================================================================================================================
      
      subroutine SCF_BALANCE_K_JOINT_OUT_OF_PANE (SCFChord,SCFBrace,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,ZINE
     1                                            ,SCFBRACECROWN,Betaparameter,Shearparameter,dummyShearparameter,Gammaparameter
     1                                                  ,DUMMYGammaparameter,Alpha,ZETA,DUMMYZETA,IMAXA,DUMMYBETA,IMAXB,NOPTION)
                
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION DUMMYZETA(IMAXA),DUMMYBETA(IMAXB),DUMMYGammaparameter(IMAXA),dummyShearparameter(IMAXA)
      DIMENSION ANGLE(IMAXA)

      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      ! zine
      !C
      !C1
      !C2
      
      SCFChord         = 5.0
      SCFBrace         = 5.0
      SCFCHORDSADDLE   = 5.0
      SCFCHORDCROWN    = 5.0
      SCFBRACESADDLE   = 5.0
      SCFBRACECROWN    = 5.0
            
      C                = 0.7D0      
      
      zetamax            = MAXVAL(DUMMYZETA)
      zetamin            = MINVAL(DUMMYZETA)
      AmaxBetaparameter  = MAXVAL(DUMMYBETA)
      AminBetaparameter  = MINVAL(DUMMYBETA)
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 )  
      
      F4 = 1.0D0 - 1.07D0*(Betaparameter**1.88)*EXP(0.16D0*(Gammaparameter**(-1.06D0))*(Alpha**(2.4D0)))
      
      !=============================================
      IF (NOPTION.EQ.1)THEN ! BALANCED AXIAL CHORD
      !=============================================
 
      X = 1.0D0 + ZINE * ( sin(zeta) ) / Betaparameter
            
      SCF_CS = (Gammaparameter * Shearparameter * Betaparameter * (1.7D0-1.05D0*(Betaparameter**3))* (sin(zeta))**1.6d0 )
     & *(1D0- 0.08D0*((Betaparameter*Gammaparameter)**0.5D0)*EXP(-0.8D0*X))
      
      SCFCHORDSADDLE = SCF_CS 
      
      !IF(Alpha.LT.12) SCFCHORDCROWN = ABS( SCF_CS * F3 )
      
      EQ10 = Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0
      
      SHORTCORD = EQ10 * F3

      EBRACE = EQ10*(1D0- 0.08D0*((Betaparameter*Gammaparameter)**0.5D0)*EXP(-0.8D0*X))

      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
     
      
      SCF_BS  =(Shearparameter **(-0.54))* (Gammaparameter**( -0.05d0 ) )* (0.99 - 0.47*Betaparameter + 0.08*Betaparameter**4d0 )
     & * SCF_CS 

      SCFBRACESADDLE = SCF_BS

      !==================================================
      ELSEIF (NOPTION.EQ.2)THEN ! UNBALANCED AXIAL CHORD
      !==================================================

      SCF_CS = (Gammaparameter * Shearparameter * Betaparameter * (1.7D0-1.05D0*(Betaparameter**3))* (sin(zeta))**1.6d0 )
     & *(1D0- 0.08D0*((Betaparameter*Gammaparameter)**0.5D0)*EXP(-0.8D0*X)) 




      ENDIF
      
      
      RETURN
      END SUBROUTINE
      
      !==========================================================================================================================     
      subroutine SCF_Joint_fomular_AXIAL_T_JOINT (Betaparameter,Shearparameter,Gammaparameter,Alpha,zeta
     1                                                            ,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      COMMON / JOINT_FIXITY / JOINTFIXITY
      
      !==============================================================================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      !==============================================================================
      
      Pi = 3.14166d0
      

      zetadegree = zeta * 180.0d0 / Pi
      !Alpha = 9.37207d0
      
      C = 0.7d0  ! chord and fixity parameter
      C1 = 2.0 * ( C - 0.5d0 ) 
      C2 = C / 2.0d0
      C3 = C / 5.0d0
      

      !==========================
      ! SHORT CHORD CORRECTION
      !==========================
      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F2 = 1d0 - (1.43d0 * Betaparameter - 0.97d0 *  Betaparameter**2d0 - 0.03d0 ) *  Gammaparameter**(0.04d0) 
     & *exp( -0.71d0 * ( Gammaparameter**-1.38d0 )  * Alpha**2.5d0 )
            
            
       !=============================================================================
      
      !xxxxx = sin(30d0*Pi/180d0)
  
      
      !Chord saddle 
      EQ1 = Gammaparameter*Shearparameter**(1.1d0)*( 1.11d0 - 3d0 * ( Betaparameter -0.52d0 )**(2.0d0)) *(sin(zeta))**1.6d0
      
      SCFCHORDSADDLE = EQ1
      
      IF (Alpha .LT. 12.0D0)   SCFCHORDSADDLE = EQ1 * F1   ! MULTIPLILE SHORT CHORD CORRECTION IF alpha   more than 12
      

      
      !Chord crown  
      EQ2 = Gammaparameter**(0.2d0) * Shearparameter * ( 2.65 + 5 * ( Betaparameter -0.65d0 )**(2.0d0) ) 
     & + Shearparameter * Betaparameter * (0.25d0 * Alpha -3.0d0)*sin(zeta)
      
      SCFCHORDCROWN = EQ2
      
      !Brace saddle 
      
      EQ3 = 1.3d0 + Gammaparameter * ( Shearparameter **(0.52d0)) * (Alpha**0.1d0) * 
     & (0.187d0 - 1.25d0 *( Betaparameter**(1.1d0)) * ( Betaparameter - 0.96d0 ) )*( sin(zeta ) )**(2.7d0 - 0.01d0 *Alpha)
      
      
      SCFBRACESADDLE = EQ3
      
       IF (Alpha .LT. 12.0D0)   SCFBRACESADDLE = EQ3 * F1   ! MULTIPLILE SHORT CHORD CORRECTION IF alpha   more than 12
      
      !Brace crown 
      EQ4 = 3d0 +(Gammaparameter**(1.2d0)) * (0.12d0 * exp(-4d0 * Betaparameter) +0.011 *(Betaparameter**2d0)- 0.045d0 )
     & + Betaparameter * Shearparameter * ( 0.1d0 * Alpha -1.2d0 ) 
      
      SCFBRACECROWN = EQ4
      
      !==========================================================================================================================
      ! Axial load-General fixity conditions
      !==========================================================================================================================
      
      IF (JOINTFIXITY .EQ. 1 ) Then ! Axial load-General fixity conditions
       
          !Chord saddle 
      EQ1 = EQ1 + C1 *(0.8d0*Alpha - 6.0d0)*Shearparameter*(Betaparameter**2d0)*(1d0-((Betaparameter**2d0))**0.5d0)*
     1   ( (sin(2d0*zeta))**2d0 ) 
      
      SCFCHORDSADDLE = EQ1
      
       if (Alpha .LT. 12.0D0)  SCFCHORDSADDLE = EQ1 * F2
       
          !Chord crown  
      EQ2 = Gammaparameter**(0.2d0) * Shearparameter * ( 2.65 + 5 * ( Betaparameter -0.65d0 )**(2.0d0) ) 
     & + Shearparameter * Betaparameter * (C2 * Alpha -3.0d0)*sin(zeta)
      
      SCFCHORDCROWN = EQ2
      
       !Brace saddle 
      if (Alpha .LT. 12.0D0) SCFBRACESADDLE = EQ3*F2
      
       !Brace crown 
      EQ4 = 3d0 +(Gammaparameter**(1.2d0)) * (0.12d0 * exp(-4d0 * Betaparameter) +0.011 *(Betaparameter**2d0)- 0.045d0 )
     & + Betaparameter * Shearparameter * ( C3 * Alpha -1.2d0 ) 
      
      SCFBRACECROWN = EQ4

      Endif
      !==========================================================================================================================

      RETURN
      End  Subroutine
      !==========================================================================================================================
      !==========================================================================================================================
      Subroutine SCF_Joint_fomular_IN_PLAND_BENDING_T_JOINT(Betaparameter,Shearparameter,Gammaparameter
     1                                    ,Alpha,zeta,SCFCHORDCROWN,SCFBRACECROWN)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      !C
      !C1
      !C2
      !=======================
      ! OUT-PLANE BENDING
      !=======================
      
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0 
      
      SCFCHORDCROWN = EQ8
      
      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
      
      SCFBRACECROWN = EQ9
      
      
      RETURN
      End  Subroutine
      !==========================================================================================================================
      !==========================================================================================================================
      Subroutine SCF_Joint_fomular_OUT_PLAND_BENDING_T_JOINT(Betaparameter,Shearparameter,Gammaparameter
     1                                    ,Alpha,zeta,SCFCHORDSADDLE,SCFBRACESADDLE)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      !C
      !C1
      !C2
      !=======================
      ! OUT-PLANE BENDING
      !=======================
      
      EQ10 = Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0
      
      SCFCHORDSADDLE = EQ10
      
      EQ11 =(Shearparameter **(-0.54))* (Gammaparameter**( -0.05d0 ) )* (0.99 - 0.47*Betaparameter + 0.08*Betaparameter**4d0 )*EQ10 
     
      SCFBRACESADDLE = EQ11
      
            
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 )  
      
      IF (Alpha .LT. 12.0D0) THEN ! SHORT CHORD CORRECTION
          
           SCFCHORDSADDLE = EQ10 * F3 ! MULTIPLILE SHORT CHORD CORRECTION IF alpha   more than 12
           SCFBRACESADDLE = EQ11 * F3 ! MULTIPLILE SHORT CHORD CORRECTION IF alpha   more than 12
          
      ENDIF
      
      RETURN
      End  Subroutine
      !==========================================================================================================================
      !==========================================================================================================================
      !==========================================================================================================================