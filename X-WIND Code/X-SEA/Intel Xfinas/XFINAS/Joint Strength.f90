!=====================================================================================================================================
!======= JOINT STRENGTH  =============================================================================================================
!=====================================================================================================================================      

      SUBROUTINE TABLEVALUEAPIQU (BETA,SHEAR,GRAMMA,ALPHA,ZETA,GAP,NNODE,CHORDDIAMETERI,INDEX
     1                            ,NOUTCORD,NNCORD,CHORDDIMETER,CHORDTHICKNESSI,NOPTION,NANALYSIS,NJOINT,IPLANE)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 NJOINT
      CHARACTER*5 NANALYSIS
      CHARACTER*5 NOPTION
      CHARACTER*3 IPLANE
      DIMENSION CHORDDIMETER(NNCORD),NOUTCORD(NNCORD)
      COMMON / LOCALFORCE / DESIGNAXIAL(500000),DESIGNSHEART(500000),DESIGNSHEARS(500000)
     1                     ,DESIGNTORSION(500000),DESIGNMOMENTT(500000),DESIGNMOMENTS(500000)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /STOREN/ IELEMENT(2000000),N1(2000000),N2(2000000)
      COMMON /STOREGEO/  IGEO(1000000),NM(1000000)
      COMMON /DESIGN/ AX(999999),AY(999999),AZ(999999),NN,NND(999999),NCOM
      COMMON /MGRAV/ NGRAV
      
      SELECTCASE (NANALYSIS)
          
      CASE ("BRACE")
      
      IF (NOPTION.EQ."AXIAL")THEN
           IF (DESIGNAXIAL(INDEX).GT.0.0D0)THEN
           ! TENSION
           NDESINGAXIAL = 1
           ELSEIF (DESIGNAXIAL(INDEX).LE.0.0D0)THEN
           ! COMPRESSION
           NDESIGNAXIAL = 2    
           ENDIF
      ELSE
      NDESIGNAXIAL = 0      
      ENDIF         
                  
      
C     GEOMETRIC FACTOR DEFINED BY
      IF (BETA.GT.0.6D0) QB  = 0.3D0/(BETA*(1-0.833*BETA))
      IF (BETA.LT.0.6D0) QB  = 1.0D0
      
C     GAP FACTOR DEFINED BY
      GD  = GAP/CHORDDIAMETERI
      IF (GD.GE.0.05D0)THEN
      QG  = 1D0+0.2D0*(1.0D0-(2.8D0*GD))**(3D0)    
            IF (QG.LT.1.0) QG = 1.0D0  
      ELSEIF (GD.LE.-0.05D0)THEN
      FREE  = THICKNESS*FYB/(T*FY)
      QG    = 0.13D0*0.65D0*FREE*(GRAMMA**(0.5D0))
      ELSE
      ! SEE IN COMMENTARY 4.3.3
      ENDIF
      
      SELECTCASE(NJOINT)
      CASE ('K')
          
          IF (NOPTION.EQ."AXIAL")THEN
          QU  = (16D0+1.2*GRAMMA)*(BETA**1.2D0)*QG
          A   = 40*(BETA**1.2D0)*QG
                IF (QU.GT.A) QU = A
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3   
          C1  = 0.2D0
          C2  = 0.2D0
          C3  = 0.3D0
                
          ELSEIF (NOPTION.EQ."INMOMENT")THEN
          QU  = (5D0+0.7D0*GRAMMA)*(BETA**1.2D0) 
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ELSEIF (NOPTION.EQ."OUTMOMENT")THEN
          QU  =  2.5D0+(4.5D0+0.2D0*GRAMMA)*(BETA**2.6D0)  
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ENDIF
                
      CASE ('T','Y')
          
          IF (NOPTION.EQ."AXIAL")THEN
          IF (NDESIGNAXIAL.EQ.1)THEN
          QU  = 30*BETA   
          ELSEIF (NDESIGNAXIAL.EQ.2)THEN
          QU  = 2.8D0+(20D0+0.8*GRAMMA)*(BETA**1.6D0)  
          A   = 2.8D0+(3.6D0*(BEATA**1.6D0))
                IF (QU.GT.A) QU = A
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
          C1  = 0.3D0
          C2  = 0D0
          C3  = 0.8D0
          
          ENDIF
          ELSEIF (NOPTION.EQ."INMOMENT")THEN
          QU  = (5D0+0.7D0*GRAMMA)*(BETA**1.2D0)
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ELSEIF (NOPTION.EQ."OUTMOMENT")THEN
          QU  =  2.5D0+(4.5D0+0.2D0*GRAMMA)*(BETA**2.6D0)      
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3      
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ENDIF
          
      CASE ('X')
          
          IF (NOPTION.EQ."AXIAL")THEN
          IF (NDESIGNAXIAL.EQ.1)THEN
             IF (BETA.LT.0.9D0)THEN
             QU  = 23D0*BETA
             ELSEIF (BETA.GT.0.9D0)THEN
             QU  = 20.7D0+(BETA-0.9D0)*(17D0*GRAMMA-220D0)
             ENDIF
C            TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             IF (BETA.LE.0.9D0)THEN
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.5D0
             ELSEIF (BETA.EQ.1.0D0)THEN
             C1  = -0.2D0
             C2  = 0D0
             C3  = 0.2D0
             ELSEIF (BETA.LT.0.9D0.AND.BETA.LT.1.0)THEN
C            LINEAR INTERPOLATION VALUE FOR X JOINT UNDER AXIAL LOAD
             BETAINTER1 = 0.4D0
             BETAINTER3 = 0.3D0
             C1 = 0.2D0-BETAINTER1*(BETA-0.9D0)
             C3 = 0.5D0-BETAINTER3*(BETA-0.9D0)
             ENDIF
          
          ELSEIF (NDESIGNAXIAL.EQ.2)THEN
          QU  = (2.8D0+(12D0+0.1*GRAMMA)*BETA)*QB   
          ENDIF
      ELSEIF (NOPTION.EQ."INMOMENT")THEN
          QU  = (5D0+0.7D0*GRAMMA)*(BETA**1.2D0)    
C            TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ELSEIF (NOPTION.EQ."OUTMOMENT")THEN
          QU  =  2.5D0+(4.5D0+0.2D0*GRAMMA)*(BETA**2.6D0)  
C            TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ENDIF
          
      ENDSELECT
      
      
      CALL LOADFACTORQF (INDEX,C1,C2,C3,CHORDDIMETER,CHORDTHICKNESSI,ZETA,QU,NJOINT,NNODE,IPLANE)
      
      
      CASE ("CHORD")
          
      DO I = 1, NNCORD
      INDEX = NOUTCORD(I)
      
      IF (DESIGNAXIAL(INDEX).GT.0.0D0)THEN
      ! TENSION
      NDESINGAXIAL = 1
      ELSEIF (DESIGNAXIAL(INDEX).LE.0.0D0)THEN
      ! COMPRESSION
      NDESIGNAXIAL = 2    
      ENDIF

      
C     GEOMETRIC FACTOR
      IF (BETA.GT.0.6D0) QB  = 0.3D0/(BETA*(1-0.833*BETA))
      IF (BETA.LT.0.6D0) QB  = 1.0D0
      
C     GAP FACTOR
      GD  = GAP/CHORDDIAMETERI
      IF (GD.GE.0.05D0)THEN
      QG  = 1D0+0.2D0*(1.0D0-(2.8D0*GD))**(3D0)    
            IF (QG.LT.1.0) QG = 1.0D0  
      ELSEIF (GD.LE.-0.05D0)THEN
      FREE  = THICKNESS*FYB/(T*FY)
      QG    = 0.13D0*0.65D0*FREE*(GRAMMA**(0.5D0))
      ELSE
      ! SEE IN COMMENTARY 4.3.3
      ENDIF
      
      SELECTCASE(NJOINT)
      CASE ('K')
          
          IF (NOPTION.EQ."AXIAL")THEN
          QU  = (16D0+1.2*GRAMMA)*(BETA**1.2D0)*QG
          A   = 40*(BETA**1.2D0)*QG
                IF (QU.GT.A) QU = A
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3   
          C1  = 0.2D0
          C2  = 0.2D0
          C3  = 0.3D0
                
          ELSEIF (NOPTION.EQ."INMOMENT")THEN
          QU  = (5D0+0.7D0*GRAMMA)*(BETA**1.2D0) 
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ELSEIF (NOPTION.EQ."OUTMOMENT")THEN
          QU  =  2.5D0+(4.5D0+0.2D0*GRAMMA)*(BETA**2.6D0)  
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ENDIF
                
      CASE ('T','Y')
          
          IF (NOPTION.EQ."AXIAL")THEN
          IF (NDESIGNAXIAL.EQ.1)THEN
          QU  = 30*BETA   
          ELSEIF (NDESIGNAXIAL.EQ.2)THEN
          QU  = 2.8D0+(20D0+0.8*GRAMMA)*(BETA**1.6D0) 
          A   = 2.8D0+(3.6D0*(BEATA**1.6D0))
                IF (QU.GT.A) QU = A
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
          C1  = 0.3D0
          C2  = 0D0
          C3  = 0.8D0
          
          ENDIF
          ELSEIF (NOPTION.EQ."INMOMENT")THEN
          QU  = (5D0+0.7D0*GRAMMA)*(BETA**1.2D0)
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ELSEIF (NOPTION.EQ."OUTMOMENT")THEN
          QU  =  2.5D0+(4.5D0+0.2D0*GRAMMA)*(BETA**2.6D0)      
C         TABLE 4.3-2 - VALUES FOR C1,C2,C3      
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ENDIF
          
      CASE ('X')
          
          IF (NOPTION.EQ."AXIAL")THEN
          IF (NDESIGNAXIAL.EQ.1)THEN
             IF (BETA.LT.0.9D0)THEN
             QU  = 23D0*BETA
             ELSEIF (BETA.GT.0.9D0)THEN
             QU  = 20.7D0+(BETA-0.9D0)*(17D0*GRAMMA-220D0)
             ENDIF
C            TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             IF (BETA.LE.0.9D0)THEN
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.5D0
             ELSEIF (BETA.EQ.1.0D0)THEN
             C1  = -0.2D0
             C2  = 0D0
             C3  = 0.2D0
             ELSEIF (BETA.LT.0.9D0.AND.BETA.LT.1.0)THEN
C            LINEAR INTERPOLATION VALUE FOR X JOINT UNDER AXIAL LOAD
             BETAINTER1 = 0.4D0
             BETAINTER3 = 0.3D0
             C1 = 0.2D0-BETAINTER1*(BETA-0.9D0)
             C3 = 0.5D0-BETAINTER3*(BETA-0.9D0)
             ENDIF
          
          ELSEIF (NDESIGNAXIAL.EQ.2)THEN
          QU  = (2.8D0+(12D0+0.1*GRAMMA)*BETA)*QB  
          ENDIF
          ELSEIF (NOPTION.EQ."INMOMENT")THEN
          QU  = (5D0+0.7D0*GRAMMA)*(BETA**1.2D0)    
C            TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ELSEIF (NOPTION.EQ."OUTMOMENT")THEN
          QU  =  2.5D0+(4.5D0+0.2D0*GRAMMA)*(BETA**2.6D0)  
C            TABLE 4.3-2 - VALUES FOR C1,C2,C3  
             C1  = 0.2D0
             C2  = 0D0
             C3  = 0.4D0
             
          ENDIF
          
      ENDSELECT
      
      ENDDO
      
          
          
      ENDSELECT
      
      END
!=====================================================================================================================================

!=====================================================================================================================================
      SUBROUTINE LOADFACTORQF (NELEMENT,C1,C2,C3,CHORDDIMETER,CHORDTHICKNESS,ZETA,QU,NJOINT,NNODE,IPLANE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)  
      CHARACTER*3 NJOINT,IPLANE
      DIMENSION AXIAL(5)
      DIMENSION FA(2)
      COMMON / LOCALFORCE / DESIGNAXIAL(500000),DESIGNSHEART(500000),DESIGNSHEARS(500000)
     1                     ,DESIGNTORSION(500000),DESIGNMOMENTT(500000),DESIGNMOMENTS(500000)
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE            
     1                        ,DROUND                 
     1                        ,BREC,HREC          
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT     
     1                        ,PLASTICT,PLASTICS,AMODULUS                                  
      COMMON /OPTION/ NDESIGN
      
C     PA    = ALLOWABLE CAPACITY FOR BRACE AXIAL LOAD
C     AMA   = ALLOWABLE FOR BRACE BENDING MOMENT
C     FS    = SAFETY FACTOR 1.6
C     PY    = THE YEILD AXIAL CAPACITY OF THE CHORD
C     AMP   = THE PLASTIC MOMENT CAPACITY OF THE CHORD
C     C1,C2,C3 = COEFFICENT DEPENDING ON JOINT AND LOAD TYPE
C     PC    = NOMINAL AXIAL LOAD
C     AMC   = BENDING MC**2 = MIPB**2 + MOPB**2
C     AMP   = PLASTIC MOMENT
C     Z     = PLASTIC SETION MODULUS
C     FA(1) = ALLOWABLE AXIAL CAPACITY T-T
C     FA(2) = ALLOWABLE AXIAL CAPACITY S-S
      NDESIGN = 2
      CALL SELECTPROP (NELEMENT,AMODULUS,POSSION,DENSITY,NSTANDARD,FU,FY)
      
      FY = 500D0*(10**6)
      FS = 1.6D0

      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      
      CALL CALLENGTH (NELEMENT,VREW,0,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                ,ATMIDX,ATMIDY,ATMIDZ)
      
         AKT        = 1.0D0
         AKS        = 1.0D0
         AUNBRACEDT = ALENGTH
         AUNBRACEDS = ALENGTH
         AC         = 0.3  ! RECOMMENDED FOR USED
         ARGS       = SQRT(AIS/AREA)
         ARGT       = SQRT(AIT/AREA)
         

      IF (DESIGNAXIAL(NELEMENT).LT.0.0D0)THEN
           ! DESIGN AXIAL VALUE
           AXIAL(1:5) = -DESIGNAXIAL(NELEMENT)
           ! ------ COMPRESSION -------
         IF (NSECTION.EQ.21)THEN
           CALL LOCALBUCKLINGAPI (AC,AKS,AKT,ARGS,ARGT,DPIPE,TPIPE,AMODULUS,ALENGTH,FY,AXIAL,AUNBRACEDS,AUNBRACEDT,FXIE)
           RATIO = DPIPE/TPIPE
           IF (RATIO.GT.60.AND.RATIO.LT.300.AND.THICKNESS.GT.VALUET) FY = FXIE
         ENDIF
           CALL COMPRESSTIONAPI (AKS,AKT,ALENGTH,FY,DPIPE,TPIPE,AC,AXIAL,FA
     1                           ,NSECTION,NCOMPACT,AUNBRACEDT,AUNBRACEDS)
          ELSEIF (DESIGNAXIAL(NELEMENT).GT.0.0D0)THEN
           ! DESIGN AXIAL VALUE
           AXIAL(1:5) = DESIGNAXIAL(NELEMENT)
           ! ------ TENSION -------
           CALL TENSIONAPI (FY,AXIAL,AREA,NSECTION,FA)
           ENDIF
      
      PY  = FA(1)                             ! ALLOWABLE AXIAL CAPACITY
      
      AMP = FY*PLASTICT                       ! PLASTIC MOMENT OF THE CHORD
      
      AMOPB = ABS(DESIGNMOMENTT(NELEMENT))
      AMIPB = ABS(DESIGNMOMENTS(NELEMENT))
      
      AMC = SQRT((AMIPB**2D0) + (AMOPB**2D0)) ! IN THE CHORD
      
C     A PARAMETER
      A  = (((FS*ABS(DESIGNAXIAL(NELEMENT))/PY)**2D0)+((FS*AMC/AMP)**2D0))**0.5D0 
      
      QFI = (1D0+C1*(FS*ABS(DESIGNAXIAL(NELEMENT))/PY)-C2*(FS*AMIPB/AMP)-C3*(A**2)) ! IN-PLANE BENDING
      QFO = (1D0+C1*(FS*ABS(DESIGNAXIAL(NELEMENT))/PY)-C2*(FS*AMIPB/AMP)-C3*(A**2)) ! OUT-PLANE BENDING
      
C     BASIC CAPACITY
      PAI  = QU*QFI*FY*(CHORDTHICKNESS**2)/(FS*SIN(ZETA))               ! IN-PLANE BENDING
      AMAI = QU*QFI*FY*(CHORDTHICKNESS**2)*CHORDDIMETER/(FS*SIN(ZETA))  ! OUT-PLANE BENDING
      
      PAO  = QU*QFO*FY*(CHORDTHICKNESS**2)/(FS*SIN(ZETA))               ! IN-PLANE BENDING
      AMAO = QU*QFO*FY*(CHORDTHICKNESS**2)*CHORDDIMETER/(FS*SIN(ZETA))  ! OUT-PLANE BENDING
C     STRENGTH CHECK 
      AIR1 = ABS(DESIGNAXIAL(NELEMENT))/(PAI)
      AIR2 = ABS(DESIGNAXIAL(NELEMENT))/(PAO)
      AIR3 = (DESIGNMOMENTT(NELEMENT)/(AMAI))**2D0
      AIR4 = ABS(DESIGNMOMENTS(NELEMENT))/(AMAO)
      
      ! MAXIMUM ALLOWABLE VALUE
      IF (AIR1.GT.AIR2) AIR5 = AIR1
      IF (AIR2.GT.AIR1) AIR5 = AIR2
      TOTALAIR = AIR5+AIR3+AIR4

      IF (TOTALAIR.GT.1.0D0)THEN     ! CONDITION FAIL > 1
      WRITE (970,10) NNODE,IPLANE,NELEMENT,NJOINT,FS,QU,QFI,TOTALAIR
10    FORMAT(I5,7X,A3,7X,I5,7X,A5,7X,F5.2,4X,F10.3,4X,F10.3,4X,F10.3,4X,"FAIL")
      ELSEIF (TOTALAIR.LE.1.0D0)THEN ! CONDITION PASS <= 1
      WRITE (970,20) NNODE,IPLANE,NELEMENT,NJOINT,FS,QU,QFI,TOTALAIR
20    FORMAT(I5,7X,A3,7X,I5,7X,A5,7X,F5.2,4X,F10.3,4X,F10.3,4X,F10.3)
      ENDIF
      
      RETURN
      END SUBROUTINE
      
!=====================================================================================================================================

!=====================================================================================================================================

!=================================================================================
!=================================================================================
!   Stress Concenstration Factor
!=================================================================================
!=================================================================================
      
      subroutine SCF_Joint_fomular_AxialLoadChordEndFixed_B1part1 (Betaparameter,Shearparameter,Gammaparameter,Alpha,zeta
     1                                                            ,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      !==============================================================================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D    (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T) (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      !==============================================================================
      
      !Chord saddle 
      EQ1 = Gammaparameter * Shearparameter**(1.1d0)  * ( 1.11d0 - 3d0 * ( Betaparameter -0.52d0 )**(2.0d0)  )
      
      SCFCHORDSADDLE = EQ1
      
      !Chord crown  
      EQ2 = Gammaparameter**(0.2d0) * Shearparameter * ( 2.65 + 5 * ( Betaparameter -0.52d0 )**(2.0d0) ) 
     & + Shearparameter * Betaparameter
      
      SCFCHORDCROWN = EQ2
      
      !Brace saddle 
      EQ3 = 1.3d0 + Gammaparameter * Shearparameter **(0.52d0) * (0.187-1.25*Betaparameter**(1.1d0) * 
     & ( Betaparameter - 0.96 ) ) * ( sin(zeta ) )**(2.7-0.01*Alpha)
      
      SCFBRACESADDLE = EQ3
      
      !Brace crown 
      EQ4 = 3d0 + Gammaparamete**(1.2d0) * (0.12 * exp(-4 * Betaparameter) +0.011 * Betaparameter**2d0
     &-0.045d0 ) + Betaparameter * Shearparameter * ( 0.1 * Alpha -1.2 ) 
      SCFBRACECROWN = EQ4
      
      end ! subroutine
!================================================================================== 
      
      subroutine SCF_Joint_fomular_fixity_Condigion_B1part2
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
      
      !Chord saddle
      
      EQ1 = Gammaparameter * Shearparameter**(1.1d0)  * ( 1.11d0 - 3d0 * ( Betaparameter -0.52d0 )**(2.0d0)  )
      
      EQ5 =  EQ1 + C1 *( 0.8d0 * Alpha - 6d0 ) * Shearparameter * Betaparameter **(2d0) *( 1-Betaparameter)**(0.5)
     &  * ( sin(zeta ) )**(2.7-0.01*Alpha) 
      
      !Chord crown 
      EQ6 = Gammaparameter**(0.2d0) * Shearparameter * ( 2.65d0 + 5d0 * ( Betaparameter - 0.65*d0 )**2 )
     & + Shearparameter * Betaparameter * ( C2 * Alpha - 3d0 ) * sin(zeta)
 
      !Brace sabble
      EQ3 = 1.3d0 + Gammaparameter * Shearparameter **(0.52d0) * (0.187-1.25*Betaparameter**(1.1d0) * 
     & ( Betaparameter - 0.96 ) ) * ( sin(zeta ) )**(2.7-0.01*Alpha)
       
      !Brace crown
      EQ7 = 3d0 + Gammaparameter**(1.2d0) * ( 0.12d0 * exp(-4 * Alpha) + 0.011d0 * Betaparameter**2d0 - 0.045 )
     &  + Betaparameter * Shearparameter * ( C3 * Alpha - 1.20d0 )
      
      
      end ! subroutine
!================================================================================== 
      
      subroutine SCF_Joint_fomular_fixity_Condigion_B1part3 (Betaparameter,Shearparameter,Gammaparameter
     1                                    ,Alpha,zeta,NOPTION,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,SCFBRACECROWN)
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
      
      !Chord saddle
       
      IF (NOPTION.EQ.1)THEN ! IN-PLANE BENDING
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0 
      
      SCFCHORDCROWN = EQ8
      
      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
      
      SCFBRACECROWN = EQ9
      
      ELSEIF (NOPTION.EQ.2)THEN ! OUT-PLANE BENDING

      EQ10 = Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0
      
      SCFCHORDSADDLE = EQ10
      
      EQ11 =(Shearparameter **(-0.54))* (Gammaparameter**( -0.05d0 ) )* (0.99 - 0.47*Betaparameter + 0.08*Betaparameter**4d0 )*EQ10 
     
      SCFBRACESADDLE = EQ11
      
      ENDIF
      
      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F2 = 1d0 - (1.43d0 * Betaparameter - 0.97d0 *  Betaparameter**2d0 - 0.03d0 ) *  Gammaparameter**(0.04d0) 
     & *exp( -0.71d0 * ( Gammaparameter**-1.38d0 )  * Alpha**2.5d0 )
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 )  
      
      end ! subroutine
!==================================================================================
      
      subroutine SCF_X_Joint_B2part1 (Betaparameter,Shearparameter,Gammaparameter
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
      EQ12 = 3.87d0 * Gammaparameter * Shearparameter * Betaparameter * (1.10d0 - Betaparameter **1.8d0)
     & *  (sin(zeta))**1.7d0
      
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
      
      
      end ! subroutine
!================================================================================= 
      
      subroutine SCF_X_Joint_B2part2 (Betaparameter,Shearparameter,Gammaparameter
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
      
      end ! subroutine
!================================================================================= 
!================================================================================= 
      
      subroutine SCF_Joint_fomular_fixity_Condigion_B2part3
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
      
       EQ1 = Gammaparameter * Shearparameter**(1.1d0)  * ( 1.11d0 - 3d0 * ( Betaparameter -0.52d0 )**(2.0d0)  )
      
      EQ5 =  EQ1 + C1 *( 0.8d0 * Alpha - 6d0 ) * Shearparameter * Betaparameter **(2d0) *( 1-Betaparameter)**(0.5)
     &  * ( sin(zeta ) )**(2.7-0.01*Alpha) 
      
      
 
      !Brace sabble
      EQ3 = 1.3d0 + Gammaparameter * Shearparameter **(0.52d0) * (0.187-1.25*Betaparameter**(1.1d0) * 
     & ( Betaparameter - 0.96 ) ) * ( sin(zeta ) )**(2.7-0.01*Alpha)
       
     
      !Chord saddle
      EQ18 = (1d0 - 0.26*Betaparameter**3d0)* EQ5
      
      !Chord crown 
      EQ6 = Gammaparameter**(0.2d0) * Shearparameter * ( 2.65d0 + 5d0 * ( Betaparameter - 0.65*d0 )**2 )
     & + Shearparameter * Betaparameter * ( C2 * Alpha - 3d0 ) * sin(zeta)
     
      !Brace saddle
      EQ19 = (1d0 - 0.26*Betaparameter**3d0)* EQ3
      
      !Brace crown
      EQ7 = 3d0 + Gammaparameter**(1.2d0) * ( 0.12d0 * exp(-4 * Alpha) + 0.011d0 * Betaparameter**2d0 - 0.045 )
     &  + Betaparameter * Shearparameter * ( C3 * Alpha - 1.20d0 )
      
      ! in joints with short chord ( Alpha < 12 ) the sabble SCFs can be reduced by the Factor 
      ! F1(fixed chord ends) or F2(pinned chord ends) where:
      
      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F2 = 1d0 - (1.43d0 * Betaparameter - 0.97d0 *  Betaparameter**2d0 - 0.03d0 ) *  Gammaparameter**(0.04d0) 
     & *exp( -0.71d0 * ( Gammaparameter**-1.38d0 )  * Alpha**2.5d0 )
      
      
      end ! subroutine
!================================================================================= 
      
      subroutine SCF_K_JOINT_B3part3 (SCFChord,SCFBrace,SCFCHORDSADDLE,SCFCHORDCROWN,SCFBRACESADDLE,ZINE
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
      zetamax            = MAXVAL(DUMMYZETA)
      zetamin            = MINVAL(DUMMYZETA)
      AmaxBetaparameter  = MAXVAL(DUMMYBETA)
      AminBetaparameter  = MINVAL(DUMMYBETA)
      
      !=============================================
      IF (NOPTION.EQ.1)THEN ! BALANCED AXIAL CHORD
      !=============================================
          
          zine = -0.078
          
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
      
      
      !=======================================================
      ELSEIF (NOPTION.EQ.2)THEN ! UNBALANCED IN PLANE BENDING
      !=======================================================   
      
      ! Chord crown
      ! for overlap exceeding 30% of contact length use 1.2 ( Eqn.8)
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0
      
      SCFCHORDCROWN = EQ8
      
      ! Gap joint brace crown:(Eqn(9))
      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
      
      ASCFBRACECROWN = EQ9
      
      ! Overlap joint brace crown:
      EQ22 = EQ9 * (0.9d0 + 0.4d0*Betaparameter)
      
      
      !=================================================================
      ELSEIF (NOPTION.EQ.3)THEN ! UNBALANCED OUT OF PLANE BENDING
      !=================================================================    
                                    
      ! Chord saddle SCF adjacent to brace A:
      EQ10= Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0

      SCFCHORDSADDLE = EQ10   
      
      EQ10A = EQ10
      EQ10B = EQ10

      ! where
      Betaparameter_A  = DUMMYBETA(1)
      Betaparameter_B  = DUMMYBETA(2)
      zetaA            = DUMMYZETA(1)
      ZetaB            = DUMMYZETA(2)
      Gammaparameter_A = DUMMYGammaparameter(1)
      Gammaparameter_B = DUMMYGammaparameter(2)
      Shearparameter_A = dummyShearparameter(1)
      Shearparameter_B = dummyShearparameter(2)
      
      EQ10A=Gammaparameter_A * Shearparameter_A * Betaparameter_A * ( 1.7d0 - 1.05d0 * Betaparameter_A**3d0 ) * (sin(zetaA))**1.6d0
      EQ10B=Gammaparameter_B * Shearparameter_B * Betaparameter_B * ( 1.7d0 - 1.05d0 * Betaparameter_B**3d0 ) * (sin(zetaB))**1.6d0
      
      
      XX=1d0 + ( ZINE * sin(zetaA) / Betaparameter_A )
      
      EQ23 =  (EQ10A)*(1d0 - 0.08d0 *((Betaparameter_B*Gammaparameter)**0.5d0 )* exp(-0.8* XX))
     & + (EQ10B) *(1d0 - 0.08d0 * (Betaparameter_A*Gammaparameter)**0.5d0* exp(-0.8* XX))
     &  *(2.05d0*(AmaxBetaparameter**0.05d0) * exp(-1.3*XX))
      
      SCFCHORDSADDLE = EQ23
      
      
      ! Brace A saddle SCF
      EQ24 = ( Shearparameter**-0.54d0 ) * (Gammaparameter**-0.05) * ( 0.99d0 - 0.47d0 * Betaparameter + 0.08d0*Betaparameter**4d0)
     &  * EQ23
      ENDIF
      
      SCFBRACESADDLE = EQ24
      
      F4 = 1d0 - 1.07d0 * (Betaparameter**1.88d0) * exp(-0.16d0 * (Gammaparameter**-1.06d0) * (Alpha**2.4d0) )
      
         
      
      end ! subroutine
      
!=================================================================================
      
      subroutine SCF_Joint_fomular_B3part1
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
      ! zine
      !C
      !C1
      !C2
      
      ! Chord
      EQ20 = ( Shearparameter ** 0.9d0 ) * (Alpha**0.5d0) * ( 0.67d0 - Betaparameter**2d0 + 1.16 * Betaparameter )
     &  * sin(zeta) *( ( sin(zetamax/zetamin))**0.30d0 ) *( (maxBetaparameter/minBetaparameter)**0.3d0 )
     &  *(1.64d0 + 2.9d0 * Betaparameter**-0.38d0 * atan(8d0*( zine )))
      
      ! Brace:
      EQ21 = 1d0 + (1.97d0 - 1.57d0 *Betaparameter**0.25d0 )*(Shearparameter**-0.14)*(sin(zeta)**0.7d0)*EQ20 +
     & (sin(zetamax + zetamin)**1.8d0)*(0.131d0-0.084d0*atan(14d0*zine + 4.2d0 * Betaparameter)) 
     &  * C * (Betaparameter**1.5d0) * (Gammaparameter**0.5d0) * (Shearparameter**-1.22d0)
      
      
      ! Chord crown
      ! for overlap exceeding 30% of contact length use 1.2 ( Eqn.8)
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0
      
      ! Gap joint brace crown:(Eqn(9))
      EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
      
      ! Overlap joint brace crown:
      EQ22 = EQ9 * (0.9d0 + 0.4d0*Betaparameter)
      
      ! Chord saddle SCF adjacent to brace A:
      EQ10 = Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0
      
      ! where
      
      XX=1d0 + ( zine * sin(zetaA) / Betaparameter_A )
      
      EQ23 =  EQ10*A*(1d0 - 0.08d0 *((Betaparameter*Gammaparameter)**0.5d0 )* exp(-0.8* XX))
     & + EQ10 * B *(1d0 - 0.08d0 * (Betaparameter*Gammaparameter)**0.5d0* exp(-0.8* XX))
     &  *(2.05d0*(MAX_Betaparameter**0.05d0) * exp(-1.3*XX))
      
      ! Brace A saddle SCF
      EQ24 = ( Shearparameter**-0.54d0 ) * (Gammaparameter**-0.05) * ( 0.99d0 - 0.47d0 * Betaparameter + 0.08d0*Betaparameter**4d0)
     &  * EQ23
      
      
      F4 = 1d0 - 1.07d0 * (Betaparameter**1.88d0) * exp(-0.16d0 * (Gammaparameter**-1.06d0) * (Alpha**2.4d0) )
      
            
      end  subroutine
!=====================================================================================================================================
!=====================================================================================================================================
!=====================================================================================================================================  
      
      subroutine SCF_Joint_fomular_fixity_Condigion_B4part1
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
      ! zine
      !C
      !C1
      !C2
      
      !=======================================================
      ! Axial Load in one Load Only
      !=======================================================
      
      ! Chord Saddle: 
      EQ5 =  EQ1 + C1 *( 0.8d0 * Alpha - 6d0 ) * Shearparameter * Betaparameter **(2d0) *( 1-Betaparameter)**(0.5)
     &  * ( sin(zeta ) )**(2.7-0.01*Alpha) 
      
      ! Chord Crown: 
      EQ6 = Gammaparameter**(0.2d0) * Shearparameter * ( 2.65d0 + 5d0 * ( Betaparameter - 0.65*d0 )**2 )
     & + Shearparameter * Betaparameter * ( C2 * Alpha - 3d0 ) * sin(zeta)
      
      !Brace saddle:
      EQ3 = 1.3d0 + Gammaparameter * Shearparameter **(0.52d0) * (0.187-1.25*Betaparameter**(1.1d0) * 
     & ( Betaparameter - 0.96 ) ) * ( sin(zeta ) )**(2.7-0.01*Alpha)
      
      !Brace Crown:
      EQ7 = 3d0 + Gammaparameter**(1.2d0) * ( 0.12d0 * exp(-4 * Alpha) + 0.011d0 * Betaparameter**2d0 - 0.045 )
     &  + Betaparameter * Shearparameter * ( C3 * Alpha - 1.20d0 )
           
      !Note that all ageometric paramenters and the resulting SCFs relate to the loaded brace.
      
      !=======================================================
      ! In-plane-bending on one brace only
      !=======================================================
      
      ! Chord saddle:
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0
      
      !Brace crown:
       EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
       
      ! Note that all geometric parameters and the resulting SCFs relate to the loaded brace.
       
      !=======================================================
      ! Out-of- plane bending on one brace only 
      !=======================================================
       
       ! where
       
       XX= 1d0 + (zine* sin(zetaA)/Betaparameter_A)
       
       EQ10 = Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0
       
       EQ25 = EQ10*A * (1d0-0.08d0*sqrt(Betaparameter * Gammaparameter )*exp(-0.8d0* XX))
       
       !Brace saddle:
       
       EQ26=(Shearparameter**-0.54d0)*(Gammaparameter**-0.05d0)*(0.99d0-0.47d0*Betaparameter+0.08*Betaparameter**4d0)*EQ25
       
      !=======================================================
      ! Short chord correction factors:
      !=======================================================
       
      F1 = 1d0 - (0.83d0 * Betaparameter - 0.56d0 *  Betaparameter**2d0 - 0.02d0 ) *  Gammaparameter**(0.23d0) 
     & *exp( -0.21d0 * ( Gammaparameter**-1.16d0 )  * Alpha**2.5d0 )
      
      F3 = 1 - 0.55d0 * (Betaparameter**1.8d0) * (Gammaparameter**(0.16d0)) 
     &*exp( -0.49d0 * ( Gammaparameter**-0.89d0 )  * Alpha**1.8d0 )  
            
      RETURN
      End  Subroutine
!=====================================================================================================================================
!=====================================================================================================================================
!=====================================================================================================================================   
      subroutine SCF_KT_JOINT (Betaparameter,Shearparameter,Gammaparameter,ZINE
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
      IF (NOPTION.EQ.1)THEN
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
      
      if ( Alpha .LT. 12 ) SCFCHORD = EQ20 * F1
      
      SCFBRACE = 1d0 + SCFCHORD *(1.97d0-1.57d0*(Betaparameter**0.25d0) ) * ( Shearparameter ** -0.14d0 ) *( (sin(zeta))**0.30d0 )
      
      
      chana =3
           
      Endif
      
      IF (IOPTION.EQ.1) RETURN
      ! Brace:
      ! K2 : API Page 211
      EQ21 = 1d0 + (1.97d0 - 1.57d0 *Betaparameter**0.25d0 )*(Shearparameter**-0.14)*(sin(zeta)**0.7d0)*EQ20 +
     & (sin(zetamax + zetamin)**1.8d0)*(0.131d0-0.084d0*atan(14d0*zine + 4.2d0 * Betaparameter)) 
     &  * C * (Betaparameter**1.5d0) * (Gammaparameter**0.5d0) * (Shearparameter**-1.22d0)
      
      SCFBRACE = EQ21!! Eq: K2 : API Page 211 
      
      !For the diagonal braces A&C use zine = zineAB + zineBC +zineBB
      !For the central brace, B , use zine = maximum of zineAB , zineBC
      
      !=======================================================
      ! In-plane bending 
      !=======================================================
      ELSEIF (NOPTION.EQ.2)THEN
      ! Chord saddle:
      ! T8 1.2
      EQ8 = 1.45d0 * Betaparameter * Shearparameter**(0.85d0)  * Gammaparameter**(1.0d0 - 0.68d0*Betaparameter)* (sin(zeta))**0.7d0
      
      SCFCHORDSADDLE = EQ8
      
      !Brace crown:
       EQ9 = 1.0d0 + 0.65d0 * Betaparameter * Shearparameter**(0.4d0)  
     &* Gammaparameter**(1.09d0 - 0.77d0*Betaparameter)* (sin(zeta))**(0.06d0*Gammaparameter - 1.16d0) 
       
       SCFBRACECROWN = EQ9
      
      !=======================================================
      ! Unbalanced out-of-plane bending
      !=======================================================
      ELSEIF (NOPTION.EQ.3)THEN
          
      BetaparameterA  = DUMMYBETA(1)
      BetaparameterB  = DUMMYBETA(2)
      BetaparameterC  = DUMMYBETA(3)
      Betaparameter_A = DUMMYBETA(1)
      Betaparameter_B = DUMMYBETA(2)
      Betaparameter_C = DUMMYBETA(3)
      zetaA            = DUMMYZETA(1)
      zetaB            = DUMMYZETA(2)
      zetaC            = DUMMYZETA(3)
      Gammaparameter_A = DUMMYGammaparameter(1)
      Gammaparameter_B = DUMMYGammaparameter(2)
      Gammaparameter_C = DUMMYGammaparameter(3)
      Shearparameter_A = dummyShearparameter(1)
      Shearparameter_B = dummyShearparameter(2)
      Shearparameter_C = dummyShearparameter(3)
      
      B_BRACE = MAXVAL(DUMMYZETA)
      
      
       ! Chord saddle SCF adjacent to diagonal brace A:
       EQ10 = Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0
       
       SCFCHORDSADDLE = EQ10
       
     
      EQ10A = Gammaparameter_A * Shearparameter_A * Betaparameter_A *( 1.7d0 - 1.05d0 * Betaparameter_A**3d0 )*(sin(zetaA))**1.6d0
      
      EQ10B = Gammaparameter_B * Shearparameter_B * Betaparameter_B *( 1.7d0 - 1.05d0 * Betaparameter_B**3d0 )*(sin(zetaB))**1.6d0
      
      EQ10C = Gammaparameter_C * Shearparameter_C * Betaparameter_C *( 1.7d0 - 1.05d0 * Betaparameter_C**3d0 )*(sin(zetaC))**1.6d0
       

       ! Where
       zineAB = ZINE
       
       zineB  = ZINE
       
       zineA  = ZINE
       
       XAB = 1d0 + ( zineAB * sin(zetaA) ) / BetaparameterA
       
       XAC = 1d0 + ( (zineAB + zineBC +zineB) * sin(zetaA) ) / BetaparameterA
       

       EQ27 = EQ10A *(1d0-0.08d0*sqrt(BetaparameterB*Gammaparameter)*exp(-0.8*XAB))
     & * ( 1d0 - 0.08d0*sqrt(BetaparameterB*Gammaparameter)*exp(-0.8d0*XAB) ) + 
     &         EQ10B *(1d0-0.08d0*sqrt(BetaparameterB*Gammaparameter)*exp(-0.8*XAB))
     & * ( 1d0 - 0.08d0*sqrt(BetaparameterB*Gammaparameter)*exp(-0.8*XAB))  +
     &         EQ10C *(1d0-0.08d0*sqrt(BetaparameterB*Gammaparameter)*exp(-0.8*XAB))
     & * ( 1d0 - 0.08d0*sqrt(BetaparameterB*Gammaparameter)*exp(-0.8*XAB))
       
        !SCFCHORDSADDLE_BRACE_A = EQ27
        
         SCFCHORDSADDLE = EQ27
         SCFCHORDCROWN  = EQ27
         
         BMAXBRACESCF = EQ27
        
      ! Where
         
       zineAB = ZINE
       
       zineBC = ZINE
       
       XAB = 1d0 + ( zineAB * sin(zetaB) ) / BetaparameterB
       
       XAC = 1d0 + ( (zineBC) * sin(zetaB) ) / BetaparameterB
       
       P1 = ( BetaparameterA / BetaparameterB )**2d0
       
       P2 = ( BetaparameterC / BetaparameterB )**2d0
       
       IF (B_BRACE.EQ.ZETA) THEN
       
       ! Chord saddle SCF adjacent to central brace: B:  
       EQ28 = EQ10B*((1d0-0.08d0*sqrt( BetaparameterB * Gammaparameter )*exp(-0.8d0*XAB))**P1 )
     & *( ( 1d0 - 0.08d0 * sqrt( BetaparameterC * Gammaparameter )  * exp(-0.8d0*XBC) )** P2 ) + 
     &  EQ10A*((1d0-0.08d0*sqrt( BetaparameterB * Gammaparameter )*exp(-0.8d0*XAB)) )
     & *( ( 2.05d0 * sqrt( AmaxBetaparameter )  * exp(-1.3d0*XAB) )) + 
     &  EQ10C*((1d0-0.08d0*sqrt( BetaparameterB * Gammaparameter )*exp(-0.8d0*XBC)) )
     & *( ( 2.05d0 * sqrt( AmaxBetaparameter )  * exp(-1.3d0*XBC) ))   
       
       SCFCHORDSADDLE_BRACE_B  = EQ28
       
        BMAXBRACESCF = EQ28
        
       ENDIF
       
       !===========================================================================================
       ! Out of plane bending bracd SCFs are obtained directly from the adjacent chord SCFs using:
       !===========================================================================================
       
      
       ! where SCF chord = ( EQ27 for EQ28 )
       
       EQ29 = (Shearparameter**-0.54d0) * (Gammaparameter**-0.05d0) * (0.99d0 - 0.47d0 * Betaparameter + 0.08d0*Betaparameter**4d0)
     &  * BMAXBRACESCF
       
       SCFBRACESADDLE   = EQ29
       SCFBRACECROWNM   = EQ29
       
       ENDIF
      end ! subroutine
!=====================================================================================================================================
!=====================================================================================================================================
!=====================================================================================================================================   
!=====================================================================================================================================   
      subroutine SCF_Joint_fomular_fixity_Condigion_B5part2
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      !=======================
      ! Parameter 
      !=======================
      ! Diameter ratio    => Betaparameter  = d/D     (brace diameter/ chord diameter)
      ! Chord Slenderness => Gammaparameter = D/(2T)  (Diameter/2/Theckness)
      ! Wall Thick ness ratio => Shearparameter = t/T (brace theckness/ chord theckness)
      ! zeta  =>  the anble between brace and chord
      ! Alpha = 2* L /D 
      ! zine
      !C
      !C1
      !C2
      
      !=======================================================
      ! Axial load on one brace only
      !=======================================================
      
      ! Chord Saddle: 
      EQ5 =  EQ1 + C1 *( 0.8d0 * Alpha - 6d0 ) * Shearparameter * Betaparameter **(2d0) *( 1-Betaparameter)**(0.5)
     &  * ( sin(zeta ) )**(2.7-0.01*Alpha) 
      
      ! Chord Crown: 
      EQ6 = Gammaparameter**(0.2d0) * Shearparameter * ( 2.65d0 + 5d0 * ( Betaparameter - 0.65*d0 )**2 )
     & + Shearparameter * Betaparameter * ( C2 * Alpha - 3d0 ) * sin(zeta)
      
      !Brace saddle:
      EQ3 = 1.3d0 + Gammaparameter * Shearparameter **(0.52d0) * (0.187-1.25*Betaparameter**(1.1d0) * 
     & ( Betaparameter - 0.96 ) ) * ( sin(zeta ) )**(2.7-0.01*Alpha)
      
      !Brace Crown:
      EQ7 = 3d0 + Gammaparameter**(1.2d0) * ( 0.12d0 * exp(-4 * Alpha) + 0.011d0 * Betaparameter**2d0 - 0.045 )
     &  + Betaparameter * Shearparameter * ( C3 * Alpha - 1.20d0 )
      
      !=======================================================
      ! Out of plane bending on one blrace only
      !=======================================================
      
      ! Chord saddle SCF adjacent to diagonal brace A:
       EQ10 = Gammaparameter * Shearparameter * Betaparameter * ( 1.7d0 - 1.05d0 * Betaparameter**3d0 ) * (sin(zeta))**1.6d0
       
       ! Where
       XAB = 1d0 + ( zineAB * sin(zetaA) ) / BetaparameterA
       
       XAC = 1d0 + ( (zineAB + zineBC +zineB) * sin(zetaA) ) / BetaparameterA
       
       EQ30 = EQ10 * (1d0 - 0.08d0 *sqrt(BetaparameterA * Gammaparameter)*exp(-0.8d0*XAB))
     & * ( ( 1d0 - 0.08d0 * sqrt( BetaparameterC * Gammaparameter )  * exp(-0.8d0*XAC) ) )
       
      !===================================================================
      ! Chord SCF adjacent to central braced B:
      !===================================================================
      
       ! Where
       XAB = 1d0 + ( zineAB * sin(zetaB) ) / BetaparameterB
       
       XAC = 1d0 + ( ( zineBC ) * sin(zetaB) ) / BetaparameterB
      
       P1 = ( BetaparameterA / BetaparameterB )**2d0
       
       P2 = ( BetaparameterC / BetaparameterB )**2d0
       
       EQ31 = EQ10 * ( (1d0 - 0.08d0 *sqrt(BetaparameterA * Gammaparameter)*exp(-0.8d0*XAB))**P1 )
     &   *( ( 1d0 - 0.08d0 * sqrt( BetaparameterC * Gammaparameter )  * exp(-0.8d0*XBC) )** P2 )
       
      !===================================================================
      ! Out of plane brace SCFs
      !===================================================================
      ! Out of plane brace SCFs are obtained directly from the adjacent chord SCFs using:
       
      ! where SCF chord = ( EQ27 for EQ28 )
       
       EQ32 = (Shearparameter**-0.54d0) * (Gammaparameter**-0.05d0) * (0.99d0 - 0.47d0 * Betaparameter + 0.08d0*Betaparameter**4d0)
     &  * SCFchord      
              
      end ! subroutine
!=====================================================================================================================================



!=====================================================================================================================================