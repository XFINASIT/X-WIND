      SUBROUTINE WAVEVECTORCOEFFICIENT(CCX,CCY,CCZ)
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)   
      
      COMMON /WAVEFRAMEELEMENTVECT/ ELEMENTFXYZ(6)
      COMMON /MGRAV/ NGRAV   
      COMMON / WRITEDATATPZANGLE / WRITEDATATPZ
      
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER 
            
      PI  = 3.141592653589793
      
      DiffXELE = ELEMENTFXYZ(4) - ELEMENTFXYZ(1)
      DiffYELE = ELEMENTFXYZ(5) - ELEMENTFXYZ(2)
      DiffZELE = ELEMENTFXYZ(6) - ELEMENTFXYZ(3)
      
      
      
      !===================================================
      IF(NGRAV .EQ. 2 ) THEN ! GRAVITY IN Z-DIRECTION
      !===================================================
          
      PHI  = ATAN( SQRT( DiffXELE**2.0D0 + DiffZELE**2.0D0 ) / DiffYELE )
      If(DiffYELE.eq.0)then
          PHI =0.0d0
      Endif
      
      ZETA = ATAN( DiffZELE / DiffXELE )
      if(DiffXELE.eq.0)then
          ZETA =0.0d0
      Endif
      
      !=============================================================
      ! GRAVITY IN Y-DIRECTION
      !
      !                       
      !                      
      !                      
      !                       
      !                       Q7
      !                       ||
      !                       ||
      !             Q3        ||          Q4
      !                       ||
      !                       ||
      !    Q6  ===============Q9==================  Q8 -------->  FX
      !                       ||
      !                       ||
      !             Q2        ||           Q1
      !                       ||
      !                       Q5
      !                       ||
      !                       ||
      !                       \/
      !                       FZ
      !
      !=============================================================
      
      if( abs(DiffXELE) .LT. 0.000000001 ) then
          
          DiffXELE = 0.0d0
          
      Endif   
      
      if( abs(DiffYELE) .LT. 0.000000001 ) then
          
          DiffYELE = 0.0d0
          
            Endif   
            
      if( abs(DiffZELE) .LT. 0.000000001 ) then
          
          DiffZELE = 0.0d0
          
      Endif   
      
      
      if( DiffXELE .GT. 0 .AND. DiffZELE .LT. 0 )THEN      ! Q1
          
          ZETA = ATAN( ABS(DiffZELE) / ABS(DiffXELE) )
      
      ELSEIF( DiffXELE .LT. 0 .AND. DiffZELE .LT. 0 )THEN  ! Q2
      
          ZETA = ATAN( ABS(DiffZELE) / ABS(DiffXELE) ) + PI*0.5D0
      
      ELSEIF( DiffXELE .LT. 0 .AND. DiffZELE .GT. 0 )THEN  ! Q3
      
          ZETA = ATAN( ABS(DiffZELE) / ABS(DiffXELE) ) + PI
      
      ELSEIF( DiffXELE .GT. 0 .AND. DiffZELE .LT. 0 )THEN  ! Q4
          
          ZETA = ATAN( ABS(DiffZELE) / ABS(DiffXELE) ) + 1.50D0*PI
      
      ELSEIF( DiffXELE .EQ. 0 .AND. DiffZELE .LT. 0 )THEN  ! Q5
          
          ZETA = PI * 0.5D0
          
      ELSEIF( DiffXELE .LT. 0 .AND. DiffZELE .EQ. 0 )THEN  ! Q6
          
          ZETA = PI
          
      ELSEIF( DiffXELE .EQ. 0 .AND. DiffZELE .GT. 0 )THEN  ! Q7
          
          ZETA = 1.50D0 * PI
          
      ELSEIF( DiffXELE .GT. 0 .AND. DiffZELE .EQ. 0 )THEN  ! Q8
          
          ZETA = 0.0D0
      
      ELSEIF( DiffXELE .EQ. 0 .AND. DiffZELE .EQ. 0 )then ! Q9
          
          ZETA =0.0D0
          
      ENDIF
     
      degree_ZETA = ZETA *180.0D0 /PI
      
      !=============================================================
      

      CCX = SIN(PHI)*COS(ZETA)
      CCY = COS(PHI)
      CCZ = SIN(PHI)*SIN(ZETA)
      
      !===================================================
      ELSEIF(NGRAV .EQ. 3 ) THEN ! GRAVITY IN Z-DIRECTION
      !===================================================
      
      PHI  = ATAN( SQRT( DiffXELE**2.0D0 + DiffYELE**2.0D0 ) / DiffZELE )
      
      degree_PHI = PHI *180.0D0 /PI
      !if(DiffZELE.eq.0)then
      !    PHI =0.0d0
      !endif
      
      !=============================================================
      ! GRAVITY IN Z-DIRECTION
      !                       FY
      !                        ^
      !                        !
      !                        !
      !                       Q5
      !                       !!
      !                       !!
      !             Q2        !!          Q1
      !                       !!
      !                       !!
      !    Q6  ===============Q9==================  Q8 -------->  FX
      !                       !!
      !                       !!
      !             Q3        !!           Q5
      !                       !!
      !                       Q7
      !
      !=============================================================
      
      if( abs(DiffXELE) .LT. 0.000000001 ) then
          
          DiffXELE = 0.0d0
          
      Endif   
      
      if( abs(DiffYELE) .LT. 0.000000001 ) then
          
          DiffYELE = 0.0d0
          
            Endif   
            
      if( abs(DiffZELE) .LT. 0.000000001 ) then
          
          DiffZELE = 0.0d0
          
      Endif   
      
      
      if( DiffXELE .GT. 0 .AND. DiffYELE .GT. 0 )THEN      ! Q1
          
          ZETA = ATAN( ABS(DiffYELE) / ABS(DiffXELE) )
      
      ELSEIF( DiffXELE .LT. 0 .AND. DiffYELE .GT. 0 )THEN  ! Q2
      
          ZETA = ATAN( ABS(DiffYELE) / ABS(DiffXELE) ) + PI*0.5D0
      
      ELSEIF( DiffXELE .LT. 0 .AND. DiffYELE .LT. 0 )THEN  ! Q3
      
          ZETA = ATAN( ABS(DiffYELE) / ABS(DiffXELE) ) + PI
      
      ELSEIF( DiffXELE .GT. 0 .AND. DiffYELE .LT. 0 )THEN  ! Q4
          
          ZETA = ATAN( ABS(DiffYELE) / ABS(DiffXELE) ) + 1.50D0*PI
      
      ELSEIF( DiffXELE .EQ. 0 .AND. DiffYELE .GT. 0 )THEN  ! Q5
          
          ZETA = PI * 0.5D0
          
      ELSEIF( DiffXELE .LT. 0 .AND. DiffYELE .EQ. 0 )THEN  ! Q6
          
          ZETA = PI
          
      ELSEIF( DiffXELE .EQ. 0 .AND. DiffYELE .LT. 0 )THEN  ! Q7
          
          ZETA = 1.50D0 * PI
          
      ELSEIF( DiffXELE .GT. 0 .AND. DiffYELE .EQ. 0 )THEN  ! Q8
          
          ZETA = 0.0D0
      
      ELSEIF( DiffXELE .EQ. 0 .AND. DiffYELE .EQ. 0 )then ! Q9
          
          ZETA =0.0D0
          
      ENDIF
     
      degree_ZETA = ZETA *180.0D0 /PI
      
      !=============================================================
      
      
      CCX = SIN(PHI)*COS(ZETA)
      CCZ = COS(PHI)
      CCY = SIN(PHI)*SIN(ZETA)

      
      ENDIF
      
      IF ( KFTTD .NE. 1 ) THEN
!      IF (WRITEDATATPZ.EQ.1) WRITE(*,1) degree_PHI , degree_ZETA , CCX , CCY , CCZ 
!1     format('PHI,ZETA,CX,CY,CZ', 2x,f10.3,2x,f10.3, 2x,f10.3,2x,f10.3,2x,f10.3 )
      ENDIF
      WRITEDATATPZ = 0
      
      CHANA = 3
      
      RETURN
      END SUBROUTINE
      
!     ================================================================================
!     ================================================================================
      
      SUBROUTINE WAVEVENORMALVECTOR(UX,UY,VNORMAL,UNX,UNY,UNZ)
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)   
      
      COMMON /MGRAV/ NGRAV 
      
      PI  = 3.141592653589793
      
      NGRAV = 3
      
      CALL WAVEVECTORCOEFFICIENT(CCX,CCY,CCZ)
      
      IF ( NGRAV .EQ. 3  ) THEN
          
      AAA =  ( CCX * UX + CCZ * UY )
      
      BBB =  ( UX**2.0D0 + UY**2.0D0 - ( CCX * UX + CCZ * UY ) )
       
      VNORMAL = SQRT( ABS( UX**2.0D0 + UY**2.0D0 - ( CCX * UX + CCZ * UY )**2.0d0 ) )
      
      UNX = UX - CCX * ( CCX * UX + CCZ * UY ) 
      
      UNY = -CCY*( CCX * UX + CCZ * UY ) 
      
      UNZ = UY - CCZ * ( CCX * UX + CCZ * UY ) 
      
      chana = 3
      
      ELSEIF ( NGRAV .EQ. 2 ) THEN
          
      VNORMAL = SQRT( ABS( UX**2.0D0 + UY**2.0D0 - ( CCX * UX + CCY * UY )**2.0d0 ) )
      
      UNX = UX - CCX * ( CCX * UX + CCY * UY ) 
      
      UNY = UY - CCY * ( CCX * UX + CCY * UY ) 
      
      UNZ = -CCZ*( CCX * UX + CCY * UY ) 
      
      ENDIF
      
      
      RETURN
      END SUBROUTINE
      
!     ================================================================================
!     ================================================================================
      
      SUBROUTINE WAVE_ANGLE_VECTOR( WFX , WFY , WFZ , VGH )
            
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)  
      
      DIMENSION VGH(3)
      
      COMMON /MGRAV/ NGRAV 
      
      PI  = 3.141592653589793
      
      CHANA = 3
      
!      SELECTCASE(NGRAV)
!      CASE(2)   
!      WAVE(1)=XWAVE
!      WAVE(2)=YWAVE
!      WAVE(3)=ZWAVE
!      CASE(3)
!      WAVE(1)=XWAVE
!      WAVE(3)=YWAVE
!      WAVE(2)=ZWAVE
!      ENDSELECT
      
      SELECTCASE(NGRAV)
      CASE(2)    
       
      IF(VGH(3).GE.0.0D0) THEN
          
      ANGLE1 = ACOS(VGH(1))
      ANGLE2 = ASIN(VGH(3))
      
      ELSE
      
      ANGLE1 = 2.0D0*PI - ACOS(VGH(1))
      ANGLE2 = 2.0D0*PI - ASIN(VGH(3))
      
      ENDIF
           
      DEG_ANGLE1 = ANGLE1*180.0D0/PI
      DEG_ANGLE2 = ANGLE2*180.0D0/PI
      
      WFX = WFX * COS(ANGLE1) + WFZ * SIN(ANGLE1)
      WFZ = WFZ * COS(ANGLE1) - WFX * SIN(ANGLE1)
          
          
      CASE(3)
          
      IF(VGH(2).GE.0.0D0) THEN
      
      ANGLE1 = ACOS(VGH(1))
      ANGLE2 = ASIN(VGH(2))
      
      ELSE
      
      ANGLE1 = 2.0D0*PI - ACOS(VGH(1))
      ANGLE2 = 2.0D0*PI - ASIN(VGH(2))
      
      ENDIF
           
      DEG_ANGLE1 = ANGLE1*180.0D0/PI
      DEG_ANGLE2 = ANGLE2*180.0D0/PI
      
      WFX = WFX * COS(ANGLE1) + WFY * SIN(ANGLE1)
      WFY = WFY * COS(ANGLE1) - WFX * SIN(ANGLE1)
      
      ENDSELECT
      
      CHANA = 3
      
      RETURN
      END SUBROUTINE
      
!     ================================================================================
!     ================================================================================
      
      SUBROUTINE INCLINEFACTOR(XYZ,MLE,OPTION,FACTOR)
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)  
      
      CHARACTER*1 OPTION
      
      DIMENSION XYZ(NCO*NNM,NELE)
      
      COMMON /MGRAV/ NGRAV 
            
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
      
      PI  = 3.141592653589793
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      
      IF( OPTION .EQ. 'X' ) THEN
          
          SELECTCASE(NGRAV)
              
          CASE(3)

          FACTOR2 = SQRT(DiffZELE**2.0D0 + DiffXELE**2.0D0)/ABS(DiffZELE)
          !IF(DiffZELE.EQ.0.0D0
          
          !WRITE(*,*) 'FACTOR', FACTOR,FACTOR2,DEGREE_AINCLINEFACTORX
          

          ANGLE = ATAN( abs(DiffZELE) / abs(DiffXELE) )
          
          !ANGLE = 0.2*PI
          
          DEGREE_ANGLE = ANGLE*180.0D0 / PI
          
          IF( DEGREE_ANGLE .GE. 45 ) THEN
              
          FACTOR = 1.0D0 / SIN(ANGLE)
          
          !FACTOR2 = SQRT( DiffZELE**2.0D0 + DiffXELE**2.0D0 ) / ABS( DiffZELE )
          
          !FACTOR3 = 1.0D0 / SIN(ANGLE)
          
          !CHANA = 3
          
          ELSE
              
          ANGLE = 0.50D0 * PI - ATAN( abs(DiffZELE) / abs(DiffXELE) )
              
          FACTOR = 1.0D0 / SIN(ANGLE)    
              
          
          ENDIF
          
          
          
          IF( DiffZELE .EQ. 1 ) THEN
              
              
              
          ENDIF
              
              
          ENDSELECT
          
      ENDIF
      

      
      
      
      RETURN 
      END SUBROUTINE

!     ================================================================================
      
      SUBROUTINE WAVEVENORMALVECTOR_Matrix(UX,UY,UZ,VNORMAL,UNX,UNY,UNZ)
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)   
      
      !COMMON /MGRAV/ NGRAV 
      
      DIMENSION A(3),UUU(3),CC(3)  , UUN(3) , VVN(3) , VDUM(3)
      
      COMMON /WAVEFRAMEELEMENTVECT/ ELEMENTFXYZ(6)
      COMMON /MGRAV/ NGRAV   
      COMMON / WRITEDATATPZANGLE / WRITEDATATPZ
      
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER 
            
      PI  = 3.141592653589793
      
      DiffXELE = ELEMENTFXYZ(4) - ELEMENTFXYZ(1)
      DiffYELE = ELEMENTFXYZ(5) - ELEMENTFXYZ(2)
      DiffZELE = ELEMENTFXYZ(6) - ELEMENTFXYZ(3)
      
      TOTAL_LENGTH = SQRT( DiffXELE**2.0D0 + DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      
      Cx = DiffXELE / TOTAL_LENGTH
      
      Cy = DiffYELE / TOTAL_LENGTH
      
      Cz = DiffZELE / TOTAL_LENGTH
      
      
      VNORMAL = SQRT( UX**2.0D0 + UY**2.0D0 + UZ**2.0D0 - ( Cx*Cx + Cy*Cy + Cz*Cz )**2.0d0 )
      
      
      UNX = (1-Cx**2)*UX -        CxCy* UY          -Cx*Cz*UZ
      UNY =  -(Cx*Cy)*UX +   (1-Cy**2)* UY          -Cy*Cz*UZ
      UNZ =  -(Cx*Cz)*UX -        CyCz* UY      +(1-Cz**2)*UZ
      
      IF( UX .EQ. 0 .AND.  UY.EQ. 0  .AND.  UZ .EQ. 0 ) then
           UNX = 0.0D0
           UNY = 0.0D0
           UNZ = 0.0D0 
           VNORMAL = 0.0D0
      ENDIF
            

      
      RETURN
      END SUBROUTINE
!     ================================================================================
!C	=======================================================================	
      SUBROUTINE VECTORCROSSWave(A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
!C	-----------------------------------------------------------------------  
      DIMENSION A(1,3),B(1,3),C(1,3)     
      
!      A x B = C

      A1 = A(1,1)
      A2 = A(1,2)
      A3 = A(1,3)

      B1 = B(1,1)
      B2 = B(1,2)
      B3 = B(1,3)

      C1 = A2*B3 - B2*A3
      C2 = A3*B1 - B3*A1
      C3 = A1*B2 - B1*A2

      C(1,1) = C1
      C(1,2) = C2
      C(1,3) = C3
      

!C    -----------------------------------------------------------------------------------   
      RETURN
      END  
!C	===================================================================================
      SUBROUTINE WAVEVENORMALVECTOR_Math_Matrix(UX,UY,UZ,VNORMAL,UNX,UNY,UNZ)
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)   
      
      !COMMON /MGRAV/ NGRAV 
      
      DIMENSION A(3),UUU(3),CC(3)  , UUN(3) , VVN(3) , VDUM(3)
      
      COMMON /WAVEFRAMEELEMENTVECT/ ELEMENTFXYZ(6), NXWAVEFRAME
      COMMON /MGRAV/ NGRAV   
      COMMON / WRITEDATATPZANGLE / WRITEDATATPZ
      
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER 
            
      PI  = 3.141592653589793d0
      
      DiffXELE = ELEMENTFXYZ(4) - ELEMENTFXYZ(1)
      DiffYELE = ELEMENTFXYZ(5) - ELEMENTFXYZ(2)
      DiffZELE = ELEMENTFXYZ(6) - ELEMENTFXYZ(3)
      
      TOTAL_LENGTH = SQRT( DiffXELE**2.0D0 + DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      
      Cx = DiffXELE / TOTAL_LENGTH
      
      Cy = DiffYELE / TOTAL_LENGTH
      
      Cz = DiffZELE / TOTAL_LENGTH
      
      ! Unit Vector of  the Structural Member
      CC(1) = Cx
      CC(2) = Cy
      CC(3) = Cz
      
      !if (NGRAV.eq.3)
      UUU(1) = UX
      UUU(2) = UY
      UUU(3) = UZ
      
      CALL VECTORCROSSWave(UUU,CC,UUN)
      CALL VECTORCROSSWave(CC,UUN,VVN)
      CALL VECTORUNIT(VVN,VDUM,VNORMAL) 
      
      UNX = VVN(1)
      UNY = VVN(2)
      UNZ = VVN(3)
      
      
      ! IF THE INPUT VECTOR IS EQUAL TO ZERO, IT WILL GENERATE THE ENTRY OF NORMAL VECTOR ENTRY EQUAL TO "NAN"
      
      IF( UX .EQ. 0 .AND.  UY.EQ. 0  .AND.  UZ .EQ. 0 ) then
           UNX = 0.0D0
           UNY = 0.0D0
           UNZ = 0.0D0 
           VNORMAL = 0.0D0
      ENDIF
            
      
      RETURN
      END SUBROUTINE
!     ================================================================================
      SUBROUTINE WAVEVERAMEDIRECTION(UX,UY,UNX,UNY,UNZ)
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)   
      
      COMMON / FRAME_WAVE_DIRECTION / BARING(3),VERTICALDIR(3)
            
      UNX = BARING(1)*UX + VERTICALDIR(1)*UY
      UNY = BARING(2)*UX + VERTICALDIR(2)*UY
      UNZ = BARING(3)*UX + VERTICALDIR(3)*UY
      
      RETURN
      END SUBROUTINE
!     ================================================================================