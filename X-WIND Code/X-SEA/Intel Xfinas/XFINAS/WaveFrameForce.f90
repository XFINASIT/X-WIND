C	=======================================================================
C                         WAVE FORCE BY PARMAIN      
C	=======================================================================
C	=======================================================================
      SUBROUTINE WAVE_FORCE(OMEGA,RATIO,H,HW,RK,RHO,CI,CD,D,X,Y,TT,IWAVE,ORDER,VR,AVAL,G,
     1                      VTIDE0,VWIND0,FTX,FTY,FTZ,NC,H0,VCURRENTP,PW,VCL,PERIOD,LWCASE,CS,
     1                      WAVEKINEMATICFACTOR,CURRENTBLOCKAGEFACTOR,WAVEKINEMATICFACTORBREAKING,WAVEBREAKINGKINEMATICFACTOR,
     1                      WAVEFORCECOEFFICIENTS,YLMIN,
     1                      VELO,ACCE,COA,IORRE,IRWAVE)
     

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)    
C	---------------------------------------------------
      DIMENSION VX(1,3),VY(1,3),VZ(1,3),VTOL(1,3),VNOL(1,3),VR(3)
      DIMENSION AX(1,3),AY(1,3),ATOL(1,3),ANOL(1,3),TRANN(14,14)
      DIMENSION VWS(3),VWZ(3),VDUM(1,3),VCL(5)
      DIMENSION DX(32),AFY(32)
      DIMENSION LWCASE(7)
      DIMENSION VELO(3),ACCE(3)
      
      
      COMMON / FRAME_WAVE_DIRECTION / BARING(3),VERTICALDIR(3)
      COMMON /MGRAV/ NGRAV
      
      
      
C      DIMENSION FUX(32),FUY(32),AUX(32),AUY(32)
      COMMON/ OFFSHOREFORCE / WAVEVELOCITYX,WAVEVELOCITYY,WAVEACCX,WAVEACCY,WAVEFORCEX,WAVEFORCEY
     1                        ,CURRENTVELOCITY,CURRENTFORCE,WAVEPLUSCURRENTX,WAVEPLUSCURRENTY
     1                        ,WAVEPLUSCURRENTZ,COORDINATEX,COORDINATEY,FDX
      
      
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV
      COMMON /WAVEFRAMEELEMENTVECT/ ELEMENTFXYZ(6), NXWAVEFRAME
      
      
      PI  = 3.141592653589793d0
      
C	---------------------------------------------------
C                       VELOCITY FORCE
C	---------------------------------------------------      
C	1) DRAG VELOCITY 
      IF (LWCASE(1).EQ.1)THEN
      CALL DRAG_VELOCITY(OMEGA,RATIO,H,HW,RK,X,Y,TT,IWAVE,ORDER,AVAL,UDX,UDY,         
     1                   AXSTREAM,AYSTREAM,G,PERIOD,IORRE,IRWAVE,AXRAN,AYRAN)
          
      ! ----- ACCELERATION FROM STREM FUNCTION ------ 
        AAX=AXSTREAM
        AAY=AYSTREAM     
      ! ---------------------------------------------  
      ELSE
          UDX = 0.0D0
          UDY = 0.0D0
          AAX = 0.0D0
          AAY = 0.0D0
      ENDIF
      
      ! ----- ACCELERATION FROM RANDOM WAVE ------ 
      IF (IORRE.EQ.2) THEN
        AAX=AXRAN
        AAY=AYRAN     
      ENDIF
        
C     ---------------------------------------------------      
C                GENERATE CURRENT PROFILES 
C     ---------------------------------------------------
C	2) CURRENT VELOCITY
      IF (NC.EQ.1D0)THEN ! DNV
       UT    = VTIDE0*((HW+Y)/HW)**(1.0/7.0)
       UW    = VWIND0*((H0+Y)/H0)
       UC    = UT+UW
      ELSEIF (NC.EQ.3D0)THEN !IEC
       UT    = VTIDE0*((HW+Y)/HW)**(1.0/7.0)
       UW    = VWIND0*((1+Y)/20D0)
       UC    = UT+UW  
      ELSEIF (NC.EQ.4)THEN ! POWER LAW
       UC    = VCURRENTP*(((HW+Y)/HW)**(PW))
      ELSEIF (NC.EQ.5.OR.NC.EQ.2)THEN ! API AND LINEAR CURRENT
       AM1   =  (VCL(1)-VCL(2))/(((HW*4D0)/4D0)-((HW*3D0)/4D0))
       AM2   =  (VCL(2)-VCL(3))/(((HW*3D0)/4D0)-((HW*2D0)/4D0))
       AM3   =  (VCL(3)-VCL(4))/(((HW*2D0)/4D0)-((HW*1D0)/4D0))
       AM4   =  (VCL(4)-VCL(5))/((HW*1D0)/4D0)
       A2    =   HW*3.0D0/4.0D0
       A3    =   HW*2.0D0/4.0D0
       A4    =   HW/4.0D0
         IF     (Y.LE.HW.AND.Y.GE.HW*3D0/4D0)         THEN
         UC  =  VCL(2)+(AM1*(Y-(HW*3D0/4D0)))
         ELSEIF (Y.LE.HW*3D0/4D0.AND.Y.GE.HW*2D0/4D0) THEN 
         UC  =  VCL(3)+(AM2*(Y-(HW*2D0/4D0)))
         ELSEIF (Y.LE.HW*2D0/4D0.AND.Y.GE.HW*1D0/4D0) THEN 
         UC  =  VCL(4)+(AM3*(Y-(HW*1D0/4D0)))
         ELSEIF (Y.LE.HW*1D0/4D0.AND.Y.GE.0.0D0)     THEN
         UC  =  VCL(5)+(AM4*(Y))
         ENDIF 
      ENDIF
          
C     --------------------------------------------------
      
C     ---------------------------------------------------         
C	      TOTAL VELOCITY SELECT BY LOAD CASE     
C     ---------------------------------------------------
C     1) COMBINATION WAVE VELOCITY + CURRENT VELOCITY  
      IF (LWCASE(1).EQ.1.AND.LWCASE(2).EQ.1.0)THEN
      UUX = UDX*WAVEKINEMATICFACTOR + UC*CURRENTBLOCKAGEFACTOR
C     2) PURE WAVE VELOCITY 
      ELSEIF (LWCASE(1).EQ.1)THEN
      UUX = UDX*WAVEKINEMATICFACTOR
C     3) PURE CURRENT VELOCITY
      ELSEIF (LWCASE(2).EQ.1.0)THEN
      UUX = UC*CURRENTBLOCKAGEFACTOR
C     4) PURE WAVE BREAKLING ( PREPARE FOR CALCULATE WAVE VELOCITY BELOW MEAN SEA LEVEL )  
      ELSEIF (LWCASE(5).EQ.1.OR.LWCASE(6).EQ.1)THEN ! BREAKING WAVE
      UUX = UDX*WAVEKINEMATICFACTORBREAKING
      ENDIF
      
      !==================================
      ! Incline Member Option
      !NXWAVEFRAME = 1  !Drawson Method
      !NXWAVEFRAME = 2  !Angle Method
      !NXWAVEFRAME = 3  !Matrix Method
      !NXWAVEFRAME = 4  !Parmni Method
      !==================================
      
      NXWAVEFRAME = 3
      
      IF( NXWAVEFRAME .EQ. 1 .or. NXWAVEFRAME .EQ. 2 .or. NXWAVEFRAME .EQ. 3  ) THEN
      
      ! Incline Vector by chana
          
      Call WAVEVERAMEDIRECTION(UUX,UDY,UdirX,UdirY,UdirZ) ! Wave Direction

      IF( NXWAVEFRAME .EQ. 1 ) THEN
      CALL WAVEVENORMALVECTOR(UUX,UDY,VNORV,UNX,UNY,UNZ)
      ELSEIF ( NXWAVEFRAME .EQ. 2 ) THEN
      CALL WAVEVENORMALVECTOR_Matrix(UdirX,UdirY,UdirZ,VNORV,UNX,UNY,UNZ)
      ELSEIF ( NXWAVEFRAME .EQ. 3 ) THEN
      CALL WAVEVENORMALVECTOR_Math_Matrix(UdirX,UdirY,UdirZ,VNORV,UNX,UNY,UNZ)
      ! PREVENT ERROR 2018
        IF (VNORV.EQ.0) CALL WAVEVENORMALVECTOR(UUX,UDY,VNORV,UNX,UNY,UNZ)
      ENDIF
      
      
      ELSE
      
      UUY = UDY   
      UUZ = 0.0D0    

C	NORMAL VELOCITY
      VX(1,1) = 1.0D0
      VX(1,2) = 0.0D0
      VX(1,3) = 0.0D0      
      VY(1,1) = 0.0D0
      VY(1,2) = 0.0D0
      VY(1,3) = 1.0D0     
      DO IV = 1,3      
      VX(1,IV)   = UUX*VX(1,IV)
      VY(1,IV)   = UUY*VY(1,IV)      
      VTOL(1,IV) = VX(1,IV) + VY(1,IV)
      ENDDO
      CALL VECTORCROSS(VTOL,VR,VZ)
      CALL VECTORCROSS(VR,VZ,VNOL)
      CALL VECTORUNIT(VNOL,VDUM,VV)             
      UN = VNOL(1,1)  
      VN = VNOL(1,2)  
      WN = VNOL(1,3)   
      
      
      ENDIF
              
C	DRAG FORCE 
      IF (KOREV.EQ.0)THEN
      STD = 0.50D0*RHO*CD*D            
            
      IF( NXWAVEFRAME .EQ. 4) THEN
      FDX = STD*VV*UN  
      FDY = STD*VV*VN 
      
      CHANA = 3
      
      ELSEIF( NXWAVEFRAME .EQ. 1 .or. NXWAVEFRAME .EQ. 2 .or. NXWAVEFRAME .EQ. 3  ) THEN
          
      STD = 0.50D0*RHO*CD*D  
      
      FDX = STD*VNORV*UNX
      FDY = STD*VNORV*UNY  
      FDZ = STD*VNORV*UNZ  
      
      
      CHANA = 3
      
      ENDIF
      
      FDXtest = STD*ABS(UUX)*UUX ! CHANA
      FDYtest = STD*ABS(UUY)*UUY ! CHANA
      

      ELSEIF (KOREV.EQ.1)THEN
      STD = 0.50D0*RHO*CD*D            
      FDX = STD*ABS(UN-VELO(1))*(UN-VELO(1)) 
      FDY = STD*ABS(UN-VELO(2))*(UN-VELO(2))
      FDZ = STD*ABS(WN-VELO(3))*(WN-VELO(3))
      ENDIF
      
!      write(*,*) FDX, FDZ
      
C     -------------------------------------------------------
C                  PLUNGING AND SURGING       
C     ------------------------------------------------------- 
      IF (LWCASE(5).EQ.1.OR.LWCASE(6).EQ.(1))THEN
        IF (LWCASE(1).EQ.0)THEN ! DISABLE WAVE ANALYSIS
           IF (Y.LE.YLMIN)THEN  
              IF (LWCASE(5).EQ.1)THEN
C	         FOR PLUNGING
               FDXP = 0.5D0*RHO*CS*D*((UDX*WAVEBREAKINGKINEMATICFACTOR)**2)
               FDYP = 0.0D0
               FDZP = 0.0D0
              ELSEIF (LWCASE(6).EQ.1)THEN
C              FOR SURGING
               FDXS = 0.5D0*RHO*CS*D*((UDX*WAVEBREAKINGKINEMATICFACTOR)**2)
               FDYS = 0.0D0
               FDZS = 0.0D0
              ENDIF
           ENDIF 
C              PDX  = WAVE FORCE + PLUNGING + SURGING               
               FDX  = FDXP + FDXS + FDX
               FDY  = FDYP + FDYS + FDY
               FDZ  = FDZP + FDZS + FDZ
        ENDIF
      ENDIF 
      
       
C	---------------------------------------------------
C                    ACCELERATION FORCE
C	---------------------------------------------------  
C	INERTIA FORCE        
      IF (IORRE.EQ.1) CALL DRAG_ACCELERATION(OMEGA,RATIO,H,HW,RK,X,Y,TT,IWAVE,ORDER,AVAL,G,AAX,AAY)  
      
      ! CURRENT VELOCITY
      IF (LWCASE(2).EQ.1.0D0.AND.LWCASE(1).EQ.0.0)THEN
      AAX = 0.0D0
      AAY = 0.0D0
      ENDIF       
     
      IF( NXWAVEFRAME .EQ. 1 .or. NXWAVEFRAME .EQ. 2 .or. NXWAVEFRAME .EQ. 3  ) THEN

      Call WAVEVERAMEDIRECTION(AAX,AAY,AAdirX,AAdirY,AAdirZ) ! Wave Direction
     
      IF( NXWAVEFRAME .EQ. 1 ) THEN
      CALL WAVEVENORMALVECTOR(AAX,AAY,AVNORA,ANNX,ANNY,ANNZ)
      ELSEIF ( NXWAVEFRAME .EQ. 2 ) THEN
      CALL WAVEVENORMALVECTOR_Matrix(AAdirX,AAdirY,AAdirZ,VNORV,ANNX,ANNY,ANNZ)
      ELSEIF ( NXWAVEFRAME .EQ. 3 ) THEN
      CALL WAVEVENORMALVECTOR_Math_Matrix(AAdirX,AAdirY,AAdirZ,VNORV,ANNX,ANNY,ANNZ)
      ENDIF

      ELSE

C	NORMAL ACCELERATION
      AX(1,1) = 1.0D0
      AX(1,2) = 0.0D0
      AX(1,3) = 0.0D0
      AY(1,1) = 0.0D0
      AY(1,2) = 1.0D0
      AY(1,3) = 0.0D0
      DO IV = 1,3            
      AX(1,IV)   = AAX*AX(1,IV)
      AY(1,IV)   = AAY*AY(1,IV)      
      ATOL(1,IV) = AX(1,IV) + AY(1,IV)
      ENDDO
      CALL VECTORCROSS(ATOL,VR,AZ)
      CALL VECTORCROSS(VR,AZ,ANOL)                
      ANX = ANOL(1,1)  
      ANY = ANOL(1,2)  
      ANZ = ANOL(1,3)
      
      ENDIF
      
      !CALL WAVEVENORMALVECTOR(AAX,AAY,AVNORA,ANNX,ANNY,ANNZ)
      !CALL WAVEVENORMALVECTOR_Matrix(AAX,AAY,VNORV,UNX,UNY,UNZ)
      

 
      IF(KOREV.EQ.0)THEN 
      !STI = RHO*CI*PI*(D**2.0)/4.0D0  
      IF( NXWAVEFRAME .EQ. 1 .or. NXWAVEFRAME .EQ. 2 .or. NXWAVEFRAME .EQ. 3  ) THEN   
      
      STI = RHO*CI*PI*(D**2.0)/4.0D0  
          
      FIX = STI*ANNX
      FIY = STI*ANNY
      FIZ = STI*ANNZ 
      
      chana = 3
      
      ELSE

      STI = RHO*CI*PI*(D**2.0)/4.0D0  
          
      FIX = STI*ANX
      FIY = STI*ANY
      FIZ = STI*ANZ      
          
      ENDIF
      
      ELSEIF(KOREV.EQ.1)THEN 
      STI = RHO*PI*(D**2.0)/4.0D0  
      FIX = STI*((CI*ANX)-(COA*ACCE(1)))
      FIY = STI*((CI*ANY)-(COA*ACCE(2)))
      FIZ = 0.0 !STI*((CI*ANZ)-(COA*ACCE(3)))
      ENDIF
      
!C     =======================================================================
!C                  FOR PRINTING RESULT MODIFY BY TOEY
!C     =======================================================================
!      COORDINATEX      = X
!      COORDINATEY      = Y
!      WAVEVELOCITYX    = UDX
!      WAVEVELOCITYY    = UDY
!      WAVEACCX         = AAX
!      WAVEACCY         = AAY
!      WAVEFORCEX       = 0.50D0*RHO*CD*D*ABS(UDX)*UDX
!      WAVEFORCEY       = 0.50D0*RHO*CD*D*ABS(UDY)*UDY
!      CURRENTVELOCITY  = UC
!      CURRENTFORCE     = 0.50D0*RHO*CD*D*ABS(UC)*UC
!      WAVEPLUSCURRENTX = FIX
!      WAVEPLUSCURRENTY = FIY
!      WAVEPLUSCURRENTZ = FIZ
!      
!      
!C     =======================================================================
!          
!C	-----------------------------------------------------------------------
!C                    TOTAL WAVE FORCE
!C	-----------------------------------------------------------------------  
          
      FTX = (FDX + FIX)*WAVEFORCECOEFFICIENTS
      FTY = (FDY + FIY)*WAVEFORCECOEFFICIENTS
      FTZ = (FDZ + FIZ)*WAVEFORCECOEFFICIENTS
      
!      WRITE (299,501) COORDINATEX,COORDINATEY,WAVEVELOCITYX,WAVEVELOCITYY,FDX,FIX,FTX
!501   FORMAT (F12.5,2X,F12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)
      
      
C	-----------------------------------------------------------------------  
      RETURN
	END   
C	=======================================================================