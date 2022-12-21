C	=======================================================================
      SUBROUTINE STOKE_COEFFICIENT(DL,A11,A13,A15,A22,A24,A33,A35,A44,A55,B22,B24,B33,B35,B44,B55,
     1                             C1,C2,C3,C4)
            
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      !H=  18.6667
      !D=  30.0
      !T=  7.72
      !G=  32.18
      
      
      pi=3.141592654
      Lo=g*(T**2)/(2*pi)
      Lox=Lo
      
      
      
      C=cosh(2d0*pi*DL)

      S=sinh(2d0*pi*DL)
!!      
      !%==========================================================================
      !%   Wave Velocity Parameters depend on kd. 
      !%==========================================================================
      
      
      
      A11 =  1/S
      
      A13 =  -(C**2)*(5.0D0*(C**2)+ 1.0D0 )/(8.0D0*(S**5))
      
      A15 = -(1184.0D0*(C**10)-1440*C**8-1992.0D0*C**6 + 2641.0*C**4 - 249.0*C**2 + 18.0D0)/(1536.0D0*S**11)
      
      A22 = 3.0D0/(8.0D0*S**4)
      
      A24 = (192.0D0*C**8 - 424.0D0*C**6 - 312.0D0*C**4 + 480.0D0*C**2 - 17.0D0 )/(768.0D0 * S**10)
      
      A33 = (13.0D0-4.0D0*C**2)/(64.0D0*S**7)
      
      A35 = (512.0D0*C**12 + 4224.0d0*C**10 - 6800*C**8 -12808.0d0*C**6  + 16704.0d0*C**4 - 3154.0d0*C**2 + 107.0d0) / 
     1  ( (4096.0d0 * S**13 )*( 6.0d0*C**2 -1))
      
      A44 = ( 80.0d0 *C**6 -816.0d0*C**4 + 1338.0d0*C**2 -197.0d0 )/((1536.0d0*S**10.0 )*(6.0d0*C**2-1.0d0))
      
      A55 = -(2880.0d0*C**10 - 72480.0d0*C**8 + 324000.0d0*C**6 - 432000.0d0*C**4 + 163470.0d0*C**2 
     1     - 16245.0d0) / ( (61440.0d0*S**11 ) *( 6.0d0*C**2 -1.0d0 )*(8.0d0*C**4.0d0 - 11.0d0*C**2 + 3.0d0))
      
      !%==========================================================================
      !%    Parameters depend on kd.
      !%==========================================================================
      
      
      B22 = (2.0d0*c**2+1)*c / ( 4.0d0 * s**3 )
      
      B24 = c*(272.0d0*c**8-504.0d0*c**6-192.0d0*c**4+322.0d0*c**2+21.0d0)/(384.0d0*s**9)
      
      B33=3d0*((8d0*c**6d0)+1d0)/64d0/s**6d0
     
      bc=88128d0*c**14d0-208224d0*c**12d0+70848d0*c**10d0

      B35=(bc+54000d0*c**8d0-21816d0*c**6d0+6264d0*c**4d0-54d0*c**2d0-81d0)/(12288d0*s**12d0*(6d0*c**2d0-1d0))
      
      B44 = c*(768.0d0*c**10-448.0d0*c**8-48.0d0*c**6+48.0d0*c**4+106.0d0*c**2-21.0d0)/(384.0d0*s**9*(6.0d0*c**2-1.0d0))
     
      bc5=192000d0*c**16d0-262720d0*c**14d0+83680d0*c**12d0+20160d0*c**10d0
      !---------------------------------------------------------------------
      A1=(bc5-7280d0*c**8d0+7160d0*c**6d0-1800d0*c**4d0-1050d0*c**2d0+225d0)
      A2=(12288d0*s**10d0*(6d0*c**2d0-1d0)*(8d0*c**4d0-11d0*c**2d0+3d0))
      !---------------------------------------------------------------------
      B55=A1/A2
      
      !%==========================================================================
      !%    Parameters depend on kd.
      !%==========================================================================
      C1=(8d0*c**4d0-8d0*c**2d0+9d0)/8d0/s**4d0
      cc=3840d0*c**12d0-4096d0*c**10d0
      C2=(cc+2592d0*c**8d0-1008d0*c**6d0+5944d0*c**4d0-1830d0*c**2d0+147d0)/(512d0*s**10d0*(6d0*c**2d0-1d0))
      C3=-1.0d0/(4.0d0*S*C)
      C4= (12.0d0*C**8 + 36.0d0*C**6 - 162.0d0*C**4 +141.0d0*C**2 -27.0d0)/(192.0d0*C*S**9)
      
      
      chana = 3

      
      RETURN
      END SUBROUTINE
!     %=============================================================================
C	=======================================================================	
C	MODIFINE BY CHANA 22-08-2012
C	=======================================================================
      SUBROUTINE STOKES_WAVELENGTH (GRAV,WDEPTH,PERIOD,WVHIGHT,ACOEFFICIENT_STOKE,RAMDA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION DLL(50001),aLamdaa1(50001),aLamdaa2(50001),aLamdaa1_2(50001),AR(50001)
 !    Varidate Value     
 !     !  H = Wave Height
 !     !  d = Water depth
 !     !  T = Time Period
 !     H=  10.66
 !     D=  22.86
 !     T=  9.27
 !     G=  9.81
 !     
 !     H=  18.6667
 !     D=  30.0
 !     T=  7.72
 !     G=  32.18
      
 !     pi=3.141592654
 !     Lo=g*(T**2)/(2*pi)
 !     Lox=Lo
 !     n=10
      
      H = WVHIGHT
      D = WDEPTH
      T = PERIOD
      G = GRAV
      
      
      pi=3.141592653589793D0
      Lo=g*(T**2)/(2*pi)
      Lox=Lo
      
      OMEGA = 2.0d0*pi/T
      
      
      n=10
     
      alamPoint1=0.1d0
      aLamPoint2=0.16d0 
      
      startDL= (d*(2*pi))/(G*T**2)
      
      CALL Airy_Wave_Number(GRAV,OMEGA,D,RK)
      lamin = 2.0d0*pi/RK
      startDL = WDEPTH/lamin
      
      
      
      write(5026,*)'Trial Stoke Fifth Order Wave Length' 
      
      write(5026,*)'DL   trial_a1   trial_a2'
      
     
      DO 10000 i=1,5000000
     
      ! DL=0.05d0+i*0.0001d0
      
      DL=startDL-i*0.00001 
     
      C=cosh(2d0*pi*DL)

      S=sinh(2d0*pi*DL)
      
      B33=3d0*((8d0*c**6d0)+1d0)/64d0/s**6d0
     
      bc=88128d0*c**14d0-208224d0*c**12d0+70848d0*c**10d0
      B35=(bc+54000d0*c**8d0-21816d0*c**6d0+6264d0*c**4d0-54d0*c**2d0-81d0)/(12288d0*s**12d0*(6d0*c**2d0-1d0))
     
      bc5=192000d0*c**16d0-262720d0*c**14d0+83680d0*c**12d0+20160d0*c**10d0
      !---------------------------------------------------------------------
      A1=(bc5-7280d0*c**8d0+7160d0*c**6d0-1800d0*c**4d0-1050d0*c**2d0+225d0)
      A2=(12288d0*s**10d0*(6d0*c**2d0-1d0)*(8d0*c**4d0-11d0*c**2d0+3d0))
      !---------------------------------------------------------------------
      B55=A1/A2
      C1=(8d0*c**4d0-8d0*c**2d0+9d0)/8d0/s**4d0
      cc=3840d0*c**12d0-4096d0*c**10d0
      C2=(cc+2592d0*c**8d0-1008d0*c**6d0+5944d0*c**4d0-1830d0*c**2d0+147d0)/(512d0*s**10d0*(6d0*c**2d0-1d0))
     
      alam11=alamPoint1
      aLam12=aLamPoint2
      
      ERROR1 = 1.0D0
      
      ERROR1YY1 = 1
      ERROR1YY2 = 1
      
      Do While(ERROR1YY1.GT.0.00001.AND.ERROR1YY2.GT.0.00001)
     
         YY1=-pi*H/d+1/DL*(aLam11+aLam11**3*B33+aLam11**5*(B35+B55))
         YY2=-pi*H/d+1/DL*(aLam12+aLam12**3*B33+aLam12**5*(B35+B55))
         
       !  write(*,*)YY1,YY2
         
         aLam13=aLam11
         aLam14=aLam12
         
         !ERROR1=abs(abs(aLam11)-abs(aLam12))
      
        aLam100=aLam13+(aLam14-aLam13)/(YY2-YY1)*(0d0-YY1)
        
        aLam11=aLam14
        aLam12=abs(aLam100)
        
        
        !a=abs(YY2)
        
        ERROR1YY1 = ABS(YY1)
        ERROR1YY2 = ABS(YY2)
        
        Resul_aLamda1=aLam12
        
        If (ERROR1YY1.LE.0.00001.AND.ERROR1YY2.LE.0.00001)then
        Resul_aLamda1=aLam12
        GOTO 2001
        endif
        
      ENDDO   
 2001 CONTINUE  
      
       
       aLamZ1=alamPoint1
       aLamZ2=alamPoint2
       
       ERROR2ZZ1 = 1
       ERROR2ZZ2 = 1
      Do While (ERROR2ZZ1.GT.0.00001.AND.ERROR2ZZ2.GT.0.00001)
     
         ZZ1=-d/Lo + DL*tanh(2d0*pi*DL)*(1+aLamZ1**2*C1+aLamZ1**4d0*C2)
         ZZ2=-d/Lo + DL*tanh(2d0*pi*DL)*(1+aLamZ2**2*C1+aLamZ2**4d0*C2)
         
         aLamZ3=aLamZ1
         aLamZ4=aLamZ2
         
         ERROR2=abs(ABS(aLamZ1)-(aLamZ2))  

         aLamZ100=aLamZ3+(aLamZ4-aLamZ3)/(ZZ2-ZZ1)*(0d0-ZZ1)
         
         aLamZ1=aLamZ4
         aLamZ2=aLamZ100
         
         ERROR2ZZ1 = ABS(ZZ1)
         ERROR2ZZ2 = ABS(ZZ2)
        
         Resul_aLamda2 = aLamZ2
            
        If (ERROR2ZZ1.LE.0.00001.AND.ERROR2ZZ2.LE.0.00001)then
            
            Resul_aLamda2 = aLamZ2
            
            CHANA = 3
            
         GOTO 4001
         endif
   
      ENDDO 
 4001 CONTINUE
      
      
      !write(5026,1)DL,Resul_aLamda1,Resul_aLamda2,(D/DL)    ! TOEY 10/2021
      
      !write(*,1)DL,Resul_aLamda1,Resul_aLamda2,(D/DL)
      
!1      format(4F14.5)      
      
      ERRORTOTAL = ABS(ABS(Resul_aLamda1)-ABS(Resul_aLamda2))
      
      if( ERRORTOTAL.LT.0.0001)then
      DLx=DLL(i-1)
      GOTO 10002
      endif
            
      
10000 CONTINUE 
10002 CONTINUE             

      ACOF1 = aLam12
      
      ACOF2 = aLamZ2
      
      ACOEFFICIENT_STOKE = aLamZ2
      
      RAMDA=D/DL    !WDEPTH/DLX
      
      chana = 3
            
      RETURN
      END SUBROUTINE
C	=======================================================================	
      
      SUBROUTINE STOKES_WAVELENGTH_SURFACE_PLOT
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
      
      COMMON /OFFSHOREx/ SEABEDx(1000),WVHIGHTx(1000),WDEPTHx(1000),THIGHTx(1000),H1POSx(1000),H2POSx(1000),
	1                  RWAVEx(1000),ORDERx(1000),PERIODx(1000),GRAVx(1000),RHOWx(1000),RHOAx(1000),
	1                  VWIND0x(1000),VTIDEx(1000),H0x(1000),
	1                  APx(1000),SPx(1000),CSx(1000),HMX(1000),HWX(1000),HCX(1000),VBREAKINGX1(1000),VBREAKINGX2(1000),
     1                  VBREAKINGX3(1000),IBREAKING(1000),
	1                  RHIGHx(1000),UHx(1000),ALPHAx(1000),Z0x(1000),
     1                  VWAVEx1(1000),VWAVEx2(1000),VWAVEx3(1000), VWINDx1(1000),VWINDx2(1000),VWINDx3(1000),
     1                  NCURRENTX(1000),POWERLAWX(1000),VCURRENTLX(1000,5),VCURRENTAPIX(1000),UHAPIX(1000),NWINDX(1000),
     1                  FACTORX(1000),VCURRENTPX(1000),AVERAGEX(1000),UHDX(1000),
     1                  WKFX(1000),CBFX(1000),WFCX(1000),
     1                  Offshoreparameterx(1000),
     1                  FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000),
     1                  TIMEX(1000)

      
      COMMON /stroke_k_and_Ramda/RKx(1000),RAMDAx(1000),ACOEFFICIENT(1000)
      
      pi = 3.141592653589793D0
      
      Interval = 5 !degree
      
      DO J = 1,NUM_OF_OFFSHORE_PARAMETER
          
      IF (RWAVEx(J) == 2) THEN
          
      CALL STOKES_WAVELENGTH(abs(GRAVx(J)),WDEPTHx(J),PERIODx(J),WVHIGHTx(J),aLamda2,RAMDA)
      
      RATIO  =   WDEPTHx(J)/RAMDA
      
      RK    = 2*pi/RAMDA
      
      CALL PARAMETER_F(RATIO,F22,F24,F33,F35,F44,F55) 
      
          AVAL = aLamda2
      
          F1 = AVAL
          F2 = (AVAL**2.0)*F22 + (AVAL**4.0)*F24
          F3 = (AVAL**3.0)*F33 + (AVAL**5.0)*F35
          F4 = (AVAL**4.0)*F44
          F5 = (AVAL**5.0)*F55
          
!      write(5026,1)                                             ! TOEY 10/2021
!1     format(' Phase Angle(degree)     Surface Elevation')      ! TOEY 10/2021
      
      DO I=1,(360/Interval + 1 )
          
      zeta = Interval*(i-1.0d0)*pi/(180.0d0)
          

      WNUDawaon = (1.0D0/RK)*(F1*COS(1.0D0*(zeta)) + 
     1                  F2*COS(2.0D0*(zeta)) + 
     1                  F3*COS(3.0D0*(zeta)) + 
     1                  F4*COS(4.0D0*(zeta)) + 
     1                  F5*COS(5.0D0*(zeta)))
      
      
      CALL STOKE_COEFFICIENT(RATIO,A11,A13,A15,A22,A24,A33,A35,A44,A55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
      
          F1 = AVAL
          F2 = (AVAL**2.0)*F22 + (AVAL**4.0)*F24
          F3 = (AVAL**3.0)*F33 + (AVAL**5.0)*F35
          F4 = (AVAL**4.0)*F44
          F5 = (AVAL**5.0)*F55
          
      WNUoriginal = (1.0D0/RK)*(F1*COS(1.0D0*(zeta)) + 
     1                  F2*COS(2.0D0*(zeta)) + 
     1                  F3*COS(3.0D0*(zeta)) + 
     1                  F4*COS(4.0D0*(zeta)) + 
     1                  F5*COS(5.0D0*(zeta)))
      
      if(i.eq.1) WAVECEST = WNUDawaon
                  
      !write(5026,2) zeta*180.0d0/pi , WNUDawaon , WNUoriginal 
!      write(5026,2) zeta*180.0d0/pi , WNUoriginal                ! TOEY 10/2021
      
2      format(4x,F10.3,10X,F12.5,4X,F12.5)
      
      CHANA = 3
      
      
      ENDDO
      
      CALL Velocity_Profile_Plot(WAVECEST)
      
      ENDIF
      
      ENDDO
      
      

      
      END SUBROUTINE
C	=======================================================================	
      SUBROUTINE Velocity_Profile_Plot(WAVECEST)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)

      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
      
      COMMON /OFFSHOREx/ SEABEDx(1000),WVHIGHTx(1000),WDEPTHx(1000),THIGHTx(1000),H1POSx(1000),H2POSx(1000),
	1                  RWAVEx(1000),ORDERx(1000),PERIODx(1000),GRAVx(1000),RHOWx(1000),RHOAx(1000),
	1                  VWIND0x(1000),VTIDEx(1000),H0x(1000),
	1                  APx(1000),SPx(1000),CSx(1000),HMX(1000),HWX(1000),HCX(1000),VBREAKINGX1(1000),VBREAKINGX2(1000),
     1                  VBREAKINGX3(1000),IBREAKING(1000),
	1                  RHIGHx(1000),UHx(1000),ALPHAx(1000),Z0x(1000),
     1                 VWAVEx1(1000),VWAVEx2(1000),VWAVEx3(1000), VWINDx1(1000),VWINDx2(1000),VWINDx3(1000),
     1                  NCURRENTX(1000),POWERLAWX(1000),VCURRENTLX(1000,5),VCURRENTAPIX(1000),UHAPIX(1000),NWINDX(1000),
     1                  FACTORX(1000),VCURRENTPX(1000),AVERAGEX(1000),UHDX(1000),
     1                  WKFX(1000),CBFX(1000),WFCX(1000),
     1                  Offshoreparameterx(1000),
     1                  FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000),
     1                  TIMEX(1000)
      
      COMMON /stroke_k_and_Ramda/RKx(1000),RAMDAx(1000),ACOEFFICIENT(1000)
      
      
      pi = 3.141592653589793D0
      
      Interval = 5 !degree
      
      write(5026,*)
!      write(5026,1)
!1     format('      ELEVATION            VELOCITY(U)     VELOCITY(W)  ACCELERATION(AX) ACCELERATION(AY)')  ! TOEY 10/2021
      
      DO J = 1,NUM_OF_OFFSHORE_PARAMETER
          
      CALL STOKES_WAVELENGTH(abs(GRAVx(J)),WDEPTHx(J),PERIODx(J),WVHIGHTx(J),aLamda2,RAMDA)
      
      HW = WDEPTHx(J)
      
      RATIO  =   WDEPTHx(J)/RAMDA
      
      RK    = 2*pi/RAMDA
      
      AVAL = aLamda2
      
      CALL PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55)       
C	-----------------------------------
      G1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
      G2 = 2.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
      G3 = 3.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
      G4 = 4.0D0*(AVAL**4.0)*G44
      G5 = 5.0D0*(AVAL**5.0)*G55
      
      
      X = H1POSx(j)
      OMEGA = 2.0D0*pi /PERIODx(J)
      
      !CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') 
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') 
      
      TT = TIME
      
      NN = (WDEPTHx(J) + WAVECEST)
!      
      DO I=1,(NN + 1 )
!          
      Y = I
      
      IF(Y.GT.(WDEPTHx(J) + WAVECEST)) Y = (WDEPTHx(J) + WAVECEST)
      
C	-----------------------------------
C	DRAG FORCE        
      UDXOLD =(G1*COSH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW)*COS(1.0D0*(RK*X
     1   - OMEGA*TT))
     1   + G2*COSH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW)*COS(2.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G3*COSH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW)*COS(3.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G4*COSH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW)*COS(4.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G5*COSH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW)*COS(5.0D0*(RK*X 
     1   - OMEGA*TT)))
     1    *(OMEGA/RK)
      UDYOLD =(G1*SINH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW)*SIN(1.0D0*(RK*X 
     1   - OMEGA*TT)) 
     1   + G2*SINH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW)*SIN(2.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G3*SINH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW)*SIN(3.0D0*(RK*X 
     1   - OMEGA*TT)) 
     1   + G4*SINH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW)*SIN(4.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G5*SINH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW)*SIN(5.0D0*(RK*X 
     1   - OMEGA*TT)))
     1    *(OMEGA/RK)  
      
      !==================================
      !    MODIFY BY CHANA 23 FEB 2016
      !==================================
      
      !WVHIGHT = H
      !WDEPTH  = HW
      !GRAV = G
!      RK=(2*(AVAL+((AVAL**3)*F33)+((AVAL**5)*(F35+F55))))/WVHIGHT 
!      
!      CWAVE = SQRT(GRAV/RK*(1.0D0 + C1*(AVAL**2.0) + C2*AVAL**4.0)*TANH(RK*WDEPTH))
!      
!      chana = (OMEGA/RK)
!      
!      FG1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
!      FG2 = 1.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
!      FG3 = 1.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
!      FG4 = 1.0D0*(AVAL**4.0)*G44
!      FG5 = 1.0D0*(AVAL**5.0)*G55
!      
!      UDX = CWAVE*( 1.0D0*FG1*COSH(1.0D0*RK*Y)*COS(1.0D0*(RK*X - OMEGA*TT))
!     1              + 2.0D0*FG2*COSH(2.0D0*RK*Y)*COS(2.0D0*(RK*X - OMEGA*TT))
!     1              + 3.0D0*FG3*COSH(3.0D0*RK*Y)*COS(3.0D0*(RK*X - OMEGA*TT))
!     1              + 4.0D0*FG4*COSH(4.0D0*RK*Y)*COS(4.0D0*(RK*X - OMEGA*TT))
!     1              + 5.0D0*FG5*COSH(5.0D0*RK*Y)*COS(5.0D0*(RK*X - OMEGA*TT)))
!      
!      UDY = CWAVE*( 1.0D0*FG1*SINH(1.0D0*RK*Y)*SIN(1.0D0*(RK*X - OMEGA*TT))
!     1              + 2.0D0*FG2*SINH(2.0D0*RK*Y)*SIN(2.0D0*(RK*X - OMEGA*TT))
!     1              + 3.0D0*FG3*SINH(3.0D0*RK*Y)*SIN(3.0D0*(RK*X - OMEGA*TT))
!     1              + 4.0D0*FG4*SINH(4.0D0*RK*Y)*SIN(4.0D0*(RK*X - OMEGA*TT))
!     1              + 5.0D0*FG5*SINH(5.0D0*RK*Y)*SIN(5.0D0*(RK*X - OMEGA*TT)))
      
      
        
C	--------------------------------------------
      U1 = G1*(COSH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW))
      U2 = G2*(COSH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW))
      U3 = G3*(COSH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW))
      U4 = G4*(COSH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW))
      U5 = G5*(COSH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW))
C	--------------------------------------------
      V1 = G1*(SINH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW))
      V2 = G2*(SINH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW))
      V3 = G3*(SINH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW))
      V4 = G4*(SINH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW))
      V5 = G5*(SINH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW))
C	---------------------------------------------
      R1 =  2.0D0*U1 - 1.0D0*U1*U2 - 1.0D0*V1*V2 - 1.0D0*U2*U3 
     1      - 1.0D0*V2*V3
      R2 =  4.0D0*U2 - 1.0D0*U1*U1 + 1.0D0*V1*V1 - 2.0D0*U1*U3 
     1      - 2.0D0*V1*V3
      R3 =  6.0D0*U3 - 3.0D0*U1*U2 + 3.0D0*V1*V2 - 3.0D0*U1*U4 
     1      - 3.0D0*V1*V4
      R4 =  8.0D0*U4 - 2.0D0*U2*U2 + 2.0D0*V2*V2 - 4.0D0*U1*U3 
     1      + 4.0D0*V1*V3
      R5 = 10.0D0*U5 - 5.0D0*U1*U4 - 5.0D0*U2*U3 + 5.0D0*V1*V4 
     1      + 5.0D0*V2*V3
C	-----------------------------------
      S0 = -2.0D0*U1*V1 + 4.0D0*U2*V2
      S1 =  2.0D0*V1    - 3.0D0*U1*V2  - 3.0D0*U2*V1  - 5.0D0*U2*V3 
     1     -5.0D0*U3*V2
      S2 =  4.0D0*V2    - 4.0D0*U1*V3  - 4.0D0*U3*V1
      S3 =  6.0D0*V3    - 1.0D0*U1*V2  + 1.0D0*U2*V1  - 5.0D0*U1*V4 
     1     -5.0D0*U4*V1
      S4 =  8.0D0*V4    - 2.0D0*U1*V3  + 2.0D0*U3*V1  + 4.0D0*U2*V2
      S5 = 10.0D0*V5    - 3.0D0*U1*V4  + 3.0D0*U4*V1  - 1.0D0*U2*V3 
     1     +1.0D0*U3*V2
      C = SQRT(G/RK*(1.0D0 + (AVAL**2.0)*C1 + (AVAL**4.0)*C2)*TANH(RK*HW))
      AAX =  (0.50D0*RK*C**2.0)*(R1*SIN(1.0D0*(RK*X - OMEGA*TT)) +
     1                           R2*SIN(2.0D0*(RK*X - OMEGA*TT)) +
     1                           R3*SIN(3.0D0*(RK*X - OMEGA*TT)) +
     1                           R4*SIN(4.0D0*(RK*X - OMEGA*TT)) +
     1                           R5*SIN(5.0D0*(RK*X - OMEGA*TT)))       
      AAY = -(0.50D0*RK*C**2.0)*(S1*COS(1.0D0*(RK*X - OMEGA*TT)) +
     1                           S2*COS(2.0D0*(RK*X - OMEGA*TT)) +
     1                           S3*COS(3.0D0*(RK*X - OMEGA*TT)) +
     1                           S4*COS(4.0D0*(RK*X - OMEGA*TT)) +
     1                           S5*COS(5.0D0*(RK*X - OMEGA*TT)))
           
                  
!       write(5026,2) Y , UDX , UDY , AAX , AAY                   ! TOEY 10/2021
! 2      format(4x,F10.3,10X,F12.5,4X,F12.5,4X,F12.5,4X,F12.5)    ! TOEY 10/2021
      
      CHANA = 3
      
      
      ENDDO
     
      ENDDO

      
      
      END SUBROUTINE
C	=======================================================================
C	=======================================================================
      SUBROUTINE Airy_Wave_Number(Gravity,OMEGA,WATERDEPTH,WAVENUMBER)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      ! MAKE BY CHANA   24.FEB.2016.
      
      ERROR = 1      ! INITIAL SETUP
      
      AK1 = 0.100D0  ! INITIAL SETUP
      
      DO WHILE ( ERROR .GT. 0.00001D0 )
                    
          AK2 = (OMEGA**2.0D0) / ( ABS(Gravity) * TANH( AK1 * WATERDEPTH))
          
          ERROR = ABS(AK2-AK1)
          
          AK1 = AK2
          
          WAVENUMBER = AK2
          
      ENDDO
      
      RETURN
      END SUBROUTINE
C	=======================================================================
C	=======================================================================