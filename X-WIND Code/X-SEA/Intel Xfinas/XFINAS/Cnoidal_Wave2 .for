        SUBROUTINE Cnoidal_wave_Length (H,waterdepth,g,T,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
       
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
       
        Dimension AKK(13,1),AMM(13,1),AEE(13,1),AHL(13,1),ALAM(13,1),OMEGA(13,1)
        Dimension AXX(13,1),EERROR(13,1),WaveK(13,1)
        Dimension AM00(6,1),AM02(6,1),AM04(6,1),AM06(6,1),AM08(6,1),AM10(6,1)
        Dimension AM12(6,1),AM14(6,1),AM16(6,1),AM18(6,1),AM20(6,1)

        pi=3.141592654d0
       
        d=waterdepth
        
        NOPTIONJACOBIAN = 3   ! Optin = 1: Dawson , 2: table approximation , 3: calculation
        
        IF ( NOPTIONJACOBIAN .EQ. 1 .or.  NOPTIONJACOBIAN .EQ. 2 ) THEN
       
        !*****************************************************************************
        ! AKK is K denotes a parameter (know as the complete elliptical integral)
        ! Reference to Offshore Structural Engineering , Thomas H. Dawson page 109
        !*****************************************************************************
       
        IF ( NOPTIONJACOBIAN .EQ. 1 ) THEN
        
        AKK(1,1)=1.571d0
        AKK(2,1)=1.612d0
        AKK(3,1)=1.660d0
        AKK(4,1)=1.714d0
        AKK(5,1)=1.778d0
        AKK(6,1)=1.854d0
        AKK(7,1)=1.950d0
        AKK(8,1)=2.075d0
        AKK(9,1)=2.257d0
        AKK(10,1)=2.578d0
        AKK(11,1)=2.908d0
        AKK(12,1)=3.696d0
        
        ELSEIF ( NOPTIONJACOBIAN .EQ. 2 ) THEN
            
        AKK(1,1) = 1.5707963d0
        AKK(2,1) = 1.5747456d0
        AKK(3,1) = 1.5868678d0
        AKK(4,1) = 1.6080486d0
        AKK(5,1) = 1.6399999d0
        AKK(6,1) = 1.6857504d0
        AKK(7,1) = 1.7507538d0
        AKK(8,1) = 1.8456940d0
        AKK(9,1) = 1.9953028d0
        AKK(10,1)= 2.2805490d0
        AKK(11,1)= 2.5900112d0
        AKK(12,1)= 3.696d0
            
        ENDIF
       
        !******************************************************************************
        ! E is a parameter know as the complete elliptic integral of the second kind
        !******************************************************************************
        
        IF ( NOPTIONJACOBIAN .EQ. 1 ) THEN
       
        AEE(1,1)=1.571d0
        AEE(2,1)=1.531d0
        AEE(3,1)=1.489d0
        AEE(4,1)=1.445d0
        AEE(5,1)=1.399d0
        AEE(6,1)=1.351d0
        AEE(7,1)=1.298d0
        AEE(8,1)=1.242d0
        AEE(9,1)=1.178d0
        AEE(10,1)=1.105d0
        AEE(11,1)=1.060d0
        AEE(12,1)=1.016d0
        
        ELSEIF ( NOPTIONJACOBIAN .EQ. 2 ) THEN
        
        
        AEE(1,1)   =    1.5707963                           
        AEE(2,1)   =    1.5668619                          
        AEE(3,1)   =    1.5549685                         
        AEE(4,1)   =    1.5348335                          
        AEE(5,1)   =    1.5059416                          
        AEE(6,1)   =    1.4674622                         
        AEE(7,1)   =    1.4180834                         
        AEE(8,1)   =    1.3556611                       
        AEE(9,1)   =    1.2763499                         
        AEE(10,1)  =    1.1716970                            
        AEE(11,1)  =    1.1027216          
        AEE(12,1)  =    1.0000000
            
          ENDIF
        
        !**************************************************************************
        ! m = Modoulus of Jacobian ellipction fucntion
        !**************************************************************************
       
        AMM(1,1)=0.0d0
        AMM(2,1)=0.1d0
        AMM(3,1)=0.2d0
        AMM(4,1)=0.3d0
        AMM(5,1)=0.4d0
        AMM(6,1)=0.5d0
        AMM(7,1)=0.6d0
        AMM(8,1)=0.7d0
        AMM(9,1)=0.8d0
        AMM(10,1)=0.9d0
        AMM(11,1)=0.95d0
        AMM(12,1)=0.99d0   
        

       
        !******************************************************************************
        
         IF ( NOPTIONJACOBIAN .EQ. 1 ) THEN
       
        AHL(1,1)=0.0d0
        AHL(2,1)=1.38d0
        AHL(3,1)=2.94d0
        AHL(4,1)=4.71d0
        AHL(5,1)=6.74d0
        AHL(6,1)=9.16d0
        AHL(7,1)=12.17d0
        AHL(8,1)=16.09d0
        AHL(9,1)=21.74d0
        AHL(10,1)=31.90d0
        AHL(11,1)=42.85d0
        AHL(12,1)=72.13d0
        
         ELSEIF ( NOPTIONJACOBIAN .EQ. 2 ) THEN
        
        AHL(1,1)   =      16.0D0/3.0D0*AMM(1,1)*( AKK(1,1)**2.0D0 )
        AHL(2,1)   =      16.0D0/3.0D0*AMM(2,1)*( AKK(2,1)**2.0D0 )
        AHL(3,1)   =      16.0D0/3.0D0*AMM(3,1)*( AKK(3,1)**2.0D0 )
        AHL(4,1)   =      16.0D0/3.0D0*AMM(4,1)*( AKK(4,1)**2.0D0 )
        AHL(5,1)   =      16.0D0/3.0D0*AMM(5,1)*( AKK(5,1)**2.0D0 )
        AHL(6,1)   =      16.0D0/3.0D0*AMM(6,1)*( AKK(6,1)**2.0D0 )
        AHL(7,1)   =      16.0D0/3.0D0*AMM(7,1)*( AKK(7,1)**2.0D0 )
        AHL(8,1)   =      16.0D0/3.0D0*AMM(8,1)*( AKK(8,1)**2.0D0 )
        AHL(9,1)   =      16.0D0/3.0D0*AMM(9,1)*( AKK(9,1)**2.0D0 )
        AHL(10,1)  =     16.0D0/3.0D0*AMM(10,1)*( AKK(10,1)**2.0D0 )
        AHL(11,1)  =     16.0D0/3.0D0*AMM(11,1)*( AKK(11,1)**2.0D0 )
        AHL(12,1)  =     16.0D0/3.0D0*AMM(12,1)*( AKK(12,1)**2.0D0 )
       
        ENDIF
           
        !*******************************************************************************
       
        ALAM(1,1)=SQRT(AHL(1,1)*(d**3)/H)
        ALAM(2,1)=SQRT(AHL(2,1)*(d**3)/H)
        ALAM(3,1)=SQRT(AHL(3,1)*(d**3)/H)
        ALAM(4,1)=SQRT(AHL(4,1)*(d**3)/H)
        ALAM(5,1)=SQRT(AHL(5,1)*(d**3)/H)
        ALAM(6,1)=SQRT(AHL(6,1)*(d**3)/H)
        ALAM(7,1)=SQRT(AHL(7,1)*(d**3)/H)
        ALAM(8,1)=SQRT(AHL(8,1)*(d**3)/H)
        ALAM(9,1)=SQRT(AHL(9,1)*(d**3)/H)
        ALAM(10,1)=SQRT(AHL(10,1)*(d**3)/H)
        ALAM(11,1)=SQRT(AHL(11,1)*(d**3)/H)
        ALAM(12,1)=SQRT(AHL(12,1)*(d**3)/H)
       

        
        !***************************************************************************
        ! OMEGA = 2* pi / T
        !***************************************************************************
       
        OMEGA(1,1)=2d0*AKK(1,1)/T
        OMEGA(2,1)=2d0*AKK(2,1)/T
        OMEGA(3,1)=2d0*AKK(3,1)/T
        OMEGA(4,1)=2d0*AKK(4,1)/T
        OMEGA(5,1)=2d0*AKK(5,1)/T
        OMEGA(6,1)=2d0*AKK(6,1)/T
        OMEGA(7,1)=2d0*AKK(7,1)/T
        OMEGA(8,1)=2d0*AKK(8,1)/T
        OMEGA(9,1)=2d0*AKK(9,1)/T
        OMEGA(10,1)=2d0*AKK(10,1)/T
        OMEGA(11,1)=2d0*AKK(11,1)/T
        OMEGA(12,1)=2d0*AKK(12,1)/T
       
       !***********************************************************************************
       
        Do 3000 j=1,12
       
        WaveK(j,1)=2d0*AKK(j,1)/ALAM(j,1)
       
            
        AXX(j,1)=g*d*(WaveK(j,1)**2d0)*((1+((0.5d0- AEE(j,1)/AKK(j,1))*H/d/AMM(j,1)))**2)
    
        EERROR(j,1)= -OMEGA(j,1)**2+AXX(j,1)
       
3000    Continue
       
        Do 4000 i=3,12
       
        CHECK=EERROR(i,1)/EERROR(i+1,1)
       
        if (CHECK.LE.0.1d0)THEN
        YYY1=EERROR(i,1)
        YYY2=EERROR(i+1,1)
       
       ALAMDA=ALAM(i,1)+(ALAM(i+1,1)-ALAM(i,1))*abs(EERROR(i,1))/(abs(EERROR(i,1))+abs(EERROR(i+1,1)))
       AKKK=AKK(i,1)+(AKK(i+1,1)-AKK(i,1))*abs(EERROR(i,1))/(abs(EERROR(i,1))+abs(EERROR(i+1,1)))
       Ammm=AMM(i,1)+(AMM(i+1,1)-AMM(i,1))*abs(EERROR(i,1))/(abs(EERROR(i,1))+abs(EERROR(i+1,1)))
       AEEE=AEE(i,1)+(AEE(i+1,1)-AEE(i,1))*abs(EERROR(i,1))/(abs(EERROR(i,1))+abs(EERROR(i+1,1)))
       
       wavenumber=2d0*AKKK/ALAMDA
       
       Omega=2d0*AKKK/T
       
       GOTO 5000
       endif
       
4000   CONTINUE
       
5000   Continue
       
     
      if (ALAMDA.EQ.0.0d0)THEN
      Write (*,*)
      Write (*,*)'*******************************************************************'
      Write (*,*)'* Environmental Wave Load is not suitable for Cnoidal Wave Theory *'
      Write (*,*)'*******************************************************************'
      Write (*,*)
      stop
      endif
      
        ELSEIF ( NOPTIONJACOBIAN .EQ. 3 ) THEN
            
      INTERVAL = 10000
      
      AINTERVAL = INTERVAL
            
      DO I=1,INTERVAL-1
            
      AMmodulous = (I-1)/AINTERVAL ! STEP INCREASE
      
      AMmodulous = 1.0D0 - (I)/AINTERVAL ! STEP DECREASE
            
            
      Call    ELLIPTIC_INTERGAL(AMmodulous,AK,AE)     
      
      
      ALAMDA = SQRT(16.0D0/3.0D0*(waterdepth**3.0D0)/H*AMmodulous*(AK**2.0D0))
      
      WAVENUM = 2.0D0*AK/ALAMDA
      
      OMEGAII = 2.0D0*AK/T
      
      FACTOR1 = g*d*(WAVENUM**2d0)*((1+((0.5d0- AE/AK)*H/d/AMmodulous))**2)
     
      FACTOR2 = OMEGAII**2.0D0
      
      
      ERROR = FACTOR1 - FACTOR2
      
      IF(ERROR.GE.0) GOTO 6000
      
!      WRITE(5039,1)AMmodulous,AK,AE,ALAMDA,FACTOR1,FACTOR2,ERROR         ! TOEY 10/2021
!1     FORMAT(10f12.5)                                                    ! TOEY 10/2021
      
      ENDDO
            
6000  CONTINUE     
      
      ALAMDA = ALAMDA
      Ammm = AMmodulous
      AKKK = AK
      AEEE = AE
            
      ENDIF      
     
        
      RETURN   
      End SUBROUTINE
    
!****************************************************************************************************************
!****************************************************************************************************************
    
      SUBROUTINE Jacobian_elliptic_function (ZZETAA,Ammm,CN,CN2,SSN,ddn)
            
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
             
      Dimension AKK(13,1),AMM(13,1),AEE(13,1),AHL(13,1),ALAM(13,1),OMEGA(13,1)
      Dimension AXX(13,1),EERROR(13,1),WaveK(13,1)
      Dimension AM00(6,1),AM02(6,1),AM04(6,1),AM06(6,1),AM08(6,1),AM10(6,1)
      Dimension AM12(6,1),AM14(6,1),AM16(6,1),AM18(6,1),AM20(6,1)
             
      !-------------------------------------------------------------------------------------------------------
            
            AM00(1,1)= 1d0
            AM00(2,1)= 1d0
            AM00(3,1)= 1d0
            AM00(4,1)= 1d0
            AM00(5,1)= 1d0
            AM00(6,1)= 1d0
            
            !---------------------------------------------------------------
            
            AM02(1,1)= 0.960d0
            AM02(2,1)= 0.960d0
            AM02(3,1)= 0.960d0
            AM02(4,1)= 0.960d0
            AM02(5,1)= 0.960d0
            AM02(6,1)= 0.960d0
            
            !--------------------------------------------------------------
            
            AM04(1,1)= 0.848d0
            AM04(2,1)= 0.850d0
            AM04(3,1)= 0.852d0
            AM04(4,1)= 0.852d0
            AM04(5,1)= 0.854d0
            AM04(6,1)= 0.856d0
            
            !--------------------------------------------------------------
            
            AM06(1,1)= 0.681d0
            AM06(2,1)= 0.687d0
            AM06(3,1)= 0.694d0
            AM06(4,1)= 0.699d0
            AM06(5,1)= 0.706d0
            AM06(6,1)= 0.712d0
            
            !--------------------------------------------------------------
            
            AM08(1,1)= 0.487d0
            AM08(2,1)= 0.500d0
            AM08(3,1)= 0.516d0
            AM08(4,1)= 0.530d0
            AM08(5,1)= 0.545d0
            AM08(6,1)= 0.560d0
            
            !--------------------------------------------------------------
            
            AM10(1,1)= 0.292d0
            AM10(2,1)= 0.317d0
            AM10(3,1)= 0.342d0
            AM10(4,1)= 0.368d0
            AM10(5,1)= 0.394d0
            AM10(6,1)= 0.420d0
            
            !--------------------------------------------------------------
            
            AM12(1,1)= 0.131d0
            AM12(2,1)= 0.162d0
            AM12(3,1)= 0.194d0
            AM12(4,1)= 0.229d0
            AM12(5,1)= 0.266d0
            AM12(6,1)= 0.305d0
            
            !--------------------------------------------------------------
            
            AM14(1,1)= 0.029d0
            AM14(2,1)= 0.053d0
            AM14(3,1)= 0.085d0
            AM14(4,1)= 0.123d0
            AM14(5,1)= 0.166d0
            AM14(6,1)= 0.216d0
            
            !--------------------------------------------------------------
            
            AM16(1,1)= 0.001d0
            AM16(2,1)= 0.003d0
            AM16(3,1)= 0.019d0
            AM16(4,1)= 0.049d0
            AM16(5,1)= 0.094d0
            AM16(6,1)= 0.151d0
            
            !--------------------------------------------------------------
            
            AM18(1,1)= 0.052d0
            AM18(2,1)= 0.016d0
            AM18(3,1)= 0.000d0
            AM18(4,1)= 0.009d0
            AM18(5,1)= 0.044d0
            AM18(6,1)= 0.104d0
            
            !--------------------------------------------------------------
            
            AM20(1,1)= 0.175d0
            AM20(2,1)= 0.062d0
            AM20(3,1)= 0.028d0
            AM20(4,1)= 0.001d0
            AM20(5,1)= 0.013d0
            AM20(6,1)= 0.071d0
            
            !--------------------------------------------------------------
            ! Inter pulation of CN() in Jacobian Elip Function
            !--------------------------------------------------------------
            
           ZZETAA=abs(ZZETAA)
   
           Do i=1,1000000
             if (ZZETAA.GE.0d0.AND.ZZETAA.LT.2.0d0)then
                Goto 991
             else if (ZZETAA.GE.2.0d0.AND.ZZETAA.LT.4.0d0)then
               AZZETAA=4d0-ZZETAA
               ZZETAA=AZZETAA
                Goto 991
  
           end if
               ZZETAA = ZZETAA-4.0d0
    
           end do
991         CONTINUE

   
           if (ZZETAA.GE.0d0.AND.ZZETAA.LT.0.2d0)Then 
           
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM00(1,1)+(Ammm-0.0d0)*(AM00(2,1)-AM00(1,1))/0.2d0
              ANNNT2=AM02(1,1)+(Ammm-0.0d0)*(AM02(2,1)-AM02(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1=AM00(2,1)+(Ammm-0.2d0)*(AM00(3,1)-AM00(2,1))/0.2d0
              ANNNT2=AM02(2,1)+(Ammm-0.2d0)*(AM02(3,1)-AM02(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM00(3,1)+(Ammm-0.4d0)*(AM00(4,1)-AM00(3,1))/0.2d0
              ANNNT2=AM02(3,1)+(Ammm-0.4d0)*(AM02(4,1)-AM02(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM00(4,1)+(Ammm-0.6d0)*(AM00(5,1)-AM00(4,1))/0.2d0
              ANNNT2=AM02(4,1)+(Ammm-0.6d0)*(AM02(5,1)-AM02(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM00(5,1)+(Ammm-0.8d0)*(AM00(6,1)-AM00(5,1))/0.2d0
              ANNNT2=AM02(5,1)+(Ammm-0.8d0)*(AM02(6,1)-AM02(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              end if
              
          else if (ZZETAA.GE.0.2d0.AND.ZZETAA.LT.0.4d0)Then 
          
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM02(1,1)+(ZZETAA-0.0d0)*(AM02(2,1)-AM02(1,1))/0.2d0
              ANNNT2=AM04(1,1)+(ZZETAA-0.0d0)*(AM04(2,1)-AM04(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.2d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1=AM02(2,1)+(ZZETAA-0.2d0)*(AM02(3,1)-AM02(2,1))/0.2d0
              ANNNT2=AM04(2,1)+(ZZETAA-0.2d0)*(AM04(3,1)-AM04(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.2d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM02(3,1)+(ZZETAA-0.4d0)*(AM02(4,1)-AM02(3,1))/0.2d0
              ANNNT2=AM04(3,1)+(ZZETAA-0.4d0)*(AM04(4,1)-AM02(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.2d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM02(4,1)+(Ammm-0.6d0)*(AM02(5,1)-AM02(4,1))/0.2d0
              ANNNT2=AM04(4,1)+(Ammm-0.6d0)*(AM04(5,1)-AM04(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.2d0)/0.2d0*(ANNNT2-ANNNT1)
            
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM02(5,1)+(Ammm-0.8d0)*(AM02(6,1)-AM02(5,1))/0.2d0
              ANNNT2=AM04(5,1)+(Ammm-0.8d0)*(AM04(6,1)-AM04(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.2d0)/0.2d0*(ANNNT2-ANNNT1)
              
              end if
          
          else if (ZZETAA.GE.0.4d0.AND.ZZETAA.LT.0.6d0)Then 
          
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM04(1,1)+(Ammm-0.0d0)*(AM04(2,1)-AM04(1,1))/0.2d0
              ANNNT2=AM06(1,1)+(Ammm-0.0d0)*(AM06(2,1)-AM06(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.4d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1 = AM04(2,1)+(Ammm-0.2d0)*(AM04(3,1)-AM04(2,1))/0.2d0
              ANNNT2 = AM06(2,1)+(Ammm-0.2d0)*(AM06(3,1)-AM06(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.4d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM04(3,1)+(Ammm-0.4d0)*(AM04(4,1)-AM04(3,1))/0.2d0
              ANNNT2=AM06(3,1)+(Ammm-0.4d0)*(AM06(4,1)-AM06(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.4d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM04(4,1)+(Ammm-0.6d0)*(AM04(5,1)-AM04(4,1))/0.2d0
              ANNNT2=AM06(4,1)+(Ammm-0.6d0)*(AM06(5,1)-AM06(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.4d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM04(5,1)+(Ammm-0.8d0)*(AM04(6,1)-AM04(5,1))/0.2d0
              ANNNT2=AM06(5,1)+(Ammm-0.8d0)*(AM06(6,1)-AM06(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.4d0)/0.2d0*(ANNNT2-ANNNT1)
              end if    
              
           else if (ZZETAA.GE.0.6d0.AND.ZZETAA.LT.0.8d0)Then 
          
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM06(1,1)+(Ammm-0.0d0)*(AM06(2,1)-AM06(1,1))/0.2d0
              ANNNT2=AM08(1,1)+(Ammm-0.0d0)*(AM08(2,1)-AM08(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.6d0)/0.2d0*(ANNNT2-ANNNT1)
                
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1=AM06(2,1)+(Ammm-0.2d0)*(AM06(3,1)-AM06(2,1))/0.2d0
              ANNNT2=AM08(2,1)+(Ammm-0.2d0)*(AM08(3,1)-AM08(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.6d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM06(3,1)+(Ammm-0.4d0)*(AM06(4,1)-AM06(3,1))/0.2d0
              ANNNT2=AM08(3,1)+(Ammm-0.4d0)*(AM08(4,1)-AM08(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.6d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM06(4,1)+(Ammm-0.6d0)*(AM06(5,1)-AM06(4,1))/0.2d0
              ANNNT2=AM08(4,1)+(Ammm-0.6d0)*(AM08(5,1)-AM08(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.6d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM06(5,1)+(Ammm-0.8d0)*(AM06(6,1)-AM06(5,1))/0.2d0
              ANNNT2=AM08(5,1)+(Ammm-0.8d0)*(AM08(6,1)-AM08(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.6d0)/0.2d0*(ANNNT2-ANNNT1)
              end if    
              
           else if (ZZETAA.GE.0.8d0.AND.ZZETAA.LT.1.0d0)Then 
          
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM08(1,1)+(Ammm-0.0d0)*(AM08(2,1)-AM08(1,1))/0.2d0
              ANNNT2=AM10(1,1)+(Ammm-0.0d0)*(AM10(2,1)-AM10(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.8d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1=AM08(2,1)+(Ammm-0.2d0)*(AM08(3,1)-AM08(2,1))/0.2d0
              ANNNT2=AM10(2,1)+(Ammm-0.2d0)*(AM10(3,1)-AM10(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.8d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM08(3,1)+(Ammm-0.4d0)*(AM08(4,1)-AM08(3,1))/0.2d0
              ANNNT2=AM10(3,1)+(Ammm-0.4d0)*(AM10(4,1)-AM10(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.8d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM08(4,1)+(Ammm-0.6d0)*(AM08(5,1)-AM08(4,1))/0.2d0
              ANNNT2=AM10(4,1)+(Ammm-0.6d0)*(AM10(5,1)-AM10(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.8d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM08(5,1)+(Ammm-0.8d0)*(AM08(6,1)-AM08(5,1))/0.2d0
              ANNNT2=AM10(5,1)+(Ammm-0.8d0)*(AM10(6,1)-AM10(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-0.8d0)/0.2d0*(ANNNT2-ANNNT1)
              end if
          
           else if (ZZETAA.GE.1.0d0.AND.ZZETAA.LT.1.2d0)Then 
          
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM10(1,1)+(Ammm-0.0d0)*(AM10(2,1)-AM10(1,1))/0.2d0
              ANNNT2=AM12(1,1)+(Ammm-0.0d0)*(AM12(2,1)-AM12(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.0d0)/0.2d0*(ANNNT2-ANNNT1)
                  
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1=AM10(2,1)+(Ammm-0.2d0)*(AM10(3,1)-AM10(2,1))/0.2d0
              ANNNT2=AM12(2,1)+(Ammm-0.2d0)*(AM12(3,1)-AM12(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.0d0)/0.2d0*(ANNNT2-ANNNT1)
                  
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM10(3,1)+(Ammm-0.4d0)*(AM10(4,1)-AM10(3,1))/0.2d0
              ANNNT2=AM12(3,1)+(Ammm-0.4d0)*(AM12(4,1)-AM12(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM10(4,1)+(Ammm-0.6d0)*(AM10(5,1)-AM10(4,1))/0.2d0
              ANNNT2=AM12(4,1)+(Ammm-0.6d0)*(AM12(5,1)-AM12(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM10(5,1)+(Ammm-0.8d0)*(AM10(6,1)-AM10(5,1))/0.2d0
              ANNNT2=AM12(5,1)+(Ammm-0.8d0)*(AM12(6,1)-AM12(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.0d0)/0.2d0*(ANNNT2-ANNNT1)
              
              end if
          
           else if (ZZETAA.GE.1.2d0.AND.ZZETAA.LT.1.4d0)Then 
          
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM12(1,1)+(Ammm-0.0d0)*(AM12(2,1)-AM12(1,1))/0.2d0
              ANNNT2=AM14(1,1)+(Ammm-0.0d0)*(AM14(2,1)-AM14(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.2d0)/0.2d0*(ANNNT2-ANNNT1)
                  
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1=AM12(2,1)+(Ammm-0.2d0)*(AM12(3,1)-AM12(2,1))/0.2d0
              ANNNT2=AM14(2,1)+(Ammm-0.2d0)*(AM14(3,1)-AM14(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.2d0)/0.2d0*(ANNNT2-ANNNT1)
                  
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM12(3,1)+(Ammm-0.4d0)*(AM12(4,1)-AM12(3,1))/0.2d0
              ANNNT2=AM14(3,1)+(Ammm-0.4d0)*(AM14(4,1)-AM14(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.2d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM12(4,1)+(Ammm-0.6d0)*(AM12(5,1)-AM12(4,1))/0.2d0
              ANNNT2=AM14(4,1)+(Ammm-0.6d0)*(AM14(5,1)-AM14(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.2d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM12(5,1)+(Ammm-0.8d0)*(AM12(6,1)-AM12(5,1))/0.2d0
              ANNNT2=AM14(5,1)+(Ammm-0.8d0)*(AM14(6,1)-AM14(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.2d0)/0.2d0*(ANNNT2-ANNNT1)
              
              end if
              
           else if (ZZETAA.GE.1.4d0.AND.ZZETAA.LT.1.6d0)Then 
          
              if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
              ANNNT1=AM14(1,1)+(Ammm-0.0d0)*(AM14(2,1)-AM14(1,1))/0.2d0
              ANNNT2=AM16(1,1)+(Ammm-0.0d0)*(AM16(2,1)-AM16(1,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.4d0)/0.2d0*(ANNNT2-ANNNT1)
                  
              else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
              ANNNT1=AM14(2,1)+(Ammm-0.2d0)*(AM14(3,1)-AM14(2,1))/0.2d0
              ANNNT2=AM16(2,1)+(Ammm-0.2d0)*(AM16(3,1)-AM16(2,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.4d0)/0.2d0*(ANNNT2-ANNNT1)
                  
              else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
              ANNNT1=AM14(3,1)+(Ammm-0.4d0)*(AM14(4,1)-AM14(3,1))/0.2d0
              ANNNT2=AM16(3,1)+(Ammm-0.4d0)*(AM16(4,1)-AM16(3,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.4d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
              ANNNT1=AM14(4,1)+(Ammm-0.6d0)*(AM14(5,1)-AM14(4,1))/0.2d0
              ANNNT2=AM16(4,1)+(Ammm-0.6d0)*(AM16(5,1)-AM16(4,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.4d0)/0.2d0*(ANNNT2-ANNNT1)
              
              else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
              ANNNT1=AM14(5,1)+(Ammm-0.8d0)*(AM14(6,1)-AM14(5,1))/0.2d0
              ANNNT2=AM16(5,1)+(Ammm-0.8d0)*(AM16(6,1)-AM16(5,1))/0.2d0
              ANNNT = ANNNT1 + (ZZETAA-1.4d0)/0.2d0*(ANNNT2-ANNNT1)
              end if
              
       else if (ZZETAA.GE.1.6d0.AND.ZZETAA.LT.1.8d0)Then 

      if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
      ANNNT1=AM16(1,1)+(Ammm-0.0d0)*(AM16(2,1)-AM16(1,1))/0.2d0
      ANNNT2=AM18(1,1)+(Ammm-0.0d0)*(AM18(2,1)-AM18(1,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.6d0)/0.2d0*(ANNNT2-ANNNT1)
        
      else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
      ANNNT1=AM16(2,1)+(Ammm-0.2d0)*(AM16(3,1)-AM16(2,1))/0.2d0
      ANNNT2=AM18(2,1)+(Ammm-0.2d0)*(AM18(3,1)-AM18(2,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.6d0)/0.2d0*(ANNNT2-ANNNT1)
        
      else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
      ANNNT1=AM16(3,1)+(Ammm-0.4d0)*(AM16(4,1)-AM16(3,1))/0.2d0
      ANNNT2=AM18(3,1)+(Ammm-0.4d0)*(AM18(4,1)-AM18(3,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.6d0)/0.2d0*(ANNNT2-ANNNT1)
    
      else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
      ANNNT1=AM16(4,1)+(Ammm-0.6d0)*(AM16(5,1)-AM16(4,1))/0.2d0
      ANNNT2=AM18(4,1)+(Ammm-0.6d0)*(AM18(5,1)-AM18(4,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.6d0)/0.2d0*(ANNNT2-ANNNT1)
    
      else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
      ANNNT1=AM16(5,1)+(Ammm-0.8d0)*(AM16(6,1)-AM16(5,1))/0.2d0
      ANNNT2=AM18(5,1)+(Ammm-0.8d0)*(AM18(6,1)-AM18(5,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.6d0)/0.2d0*(ANNNT2-ANNNT1)
      end if    
    
      else if (ZZETAA.GE.1.8d0.AND.ZZETAA.LE.2.0d0)Then 

      if (Ammm.GE.0d0.AND.Ammm.LT.0.2d0)then 
      ANNNT1=AM18(1,1)+(Ammm-0.0d0)*(AM18(2,1)-AM18(1,1))/0.2d0
      ANNNT2=AM20(1,1)+(Ammm-0.0d0)*(AM20(2,1)-AM20(1,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.8d0)/0.2d0*(ANNNT2-ANNNT1)
        
      else if ((Ammm.GE.0.2d0).AND.(Ammm.LT.0.4d0))then 
      ANNNT1=AM18(2,1)+(Ammm-0.2d0)*(AM18(3,1)-AM18(2,1))/0.2d0
      ANNNT2=AM20(2,1)+(Ammm-0.2d0)*(AM20(3,1)-AM20(2,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.8d0)/0.2d0*(ANNNT2-ANNNT1)
    
      else if ((Ammm.GE.0.4d0).AND.(Ammm.LT.0.6d0))then
      ANNNT1=AM18(3,1)+(Ammm-0.4d0)*(AM18(4,1)-AM18(3,1))/0.2d0
      ANNNT2=AM20(3,1)+(Ammm-0.4d0)*(AM20(4,1)-AM20(3,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.8d0)/0.2d0*(ANNNT2-ANNNT1)
    
      else if ((Ammm.GE.0.6d0).AND.(Ammm.LT.0.8d0))then 
      ANNNT1=AM18(4,1)+(Ammm-0.6d0)*(AM18(5,1)-AM18(4,1))/0.2d0
      ANNNT2=AM20(4,1)+(Ammm-0.6d0)*(AM20(5,1)-AM20(4,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.8d0)/0.2d0*(ANNNT2-ANNNT1)
    
      else if ((Ammm.GE.0.8d0).AND.(Ammm.LE.1.0d0))then 
      ANNNT1=AM18(5,1)+(Ammm-0.8d0)*(AM18(6,1)-AM18(5,1))/0.2d0
      ANNNT2=AM20(5,1)+(Ammm-0.8d0)*(AM20(6,1)-AM20(5,1))/0.2d0
      ANNNT = ANNNT1 + (ZZETAA-1.8d0)/0.2d0*(ANNNT2-ANNNT1)
      end if    

      end if  
      write (10000,*)'m=', Ammm
      write (10000,*)'NNT=', ANNNT
  
    
      CN2=ANNNT
      CN=sqrt(ANNNT)
   
      SSN2=1d0-ANNNT
      SSN=sqrt(SSN2)
    !
      ddn2=1d0-Ammm*(1d0-ANNNT)
      ddn=sqrt(ddn2)
    
      write (10000,*)'SSN**2=', SSN2
      write (10000,*)'SSN  =', SSN
      write (10000,*)'ddn**2=', ddn2
      write (10000,*)'ddn  =', ddn
                                             
      RETURN
      End SUBROUTINE
    
!     ***********************************************************************************
!     Peak calculation                                                                   
!     ***********************************************************************************
    
      SUBROUTINE Cnoidal_Peak  (H,waterdepth,g,T,peak,ht)
    
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
    
      Dimension AKK(13,1),AMM(13,1),AEE(13,1),AHL(13,1),ALAM(13,1),OMEGA(13,1)
      Dimension AXX(13,1),EERROR(13,1),WaveK(13,1)
      Dimension AM00(6,1),AM02(6,1),AM04(6,1),AM06(6,1),AM08(6,1),AM10(6,1)
      Dimension AM12(6,1),AM14(6,1),AM16(6,1),AM18(6,1),AM20(6,1)
      
      
      
      ZZETAA=0d0
    
      CALL Cnoidal_wave_Length (H,waterdepth,g,T,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
    
      WAVE_trough=H*(AKKK*(1d0-Ammm)-AEEE)/(Ammm*AKKK)
      CALL Jacobian_elliptic_function (ZZETAA,Ammm,CN,CN2,SSN,ddn)
    
      peak = H*(CN**2d0)+waterdepth-abs(WAVE_trough)
      ht   = peak-H

      RETURN
      End SUBROUTINE

!*********************************************************************************************************************
!*********************************************************************************************************************
!*********************************************************************************************************************    

      SUBROUTINE Cnoidal_Wave_Kinematic  ( H,waterdepth,g,T,xx,ZZ,wavenumber,AKKK,ALAMDA,Ammm,AEEE,CN,SSN,ddn,ZZETAA,UUUXX,
     1VVVXX,AAAXX,AAAZZ )
     
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      Dimension AKK(13,1),AMM(13,1),AEE(13,1),AHL(13,1),ALAM(13,1),OMEGA(13,1)
      Dimension AXX(13,1),EERROR(13,1),WaveK(13,1)
      Dimension AM00(6,1),AM02(6,1),AM04(6,1),AM06(6,1),AM08(6,1),AM10(6,1)
      Dimension AM12(6,1),AM14(6,1),AM16(6,1),AM18(6,1),AM20(6,1) 
      
      
      d=waterdepth
 
      WAVE_trough=-H*(AKKK*(1-Ammm)-AEEE)/Ammm/AKKK 
 
      peak = H+waterdepth-abs(WAVE_trough)
      ht   = peak-H
    
    
        AMINNT=ht
        
        Uphase = 0.0d0
        
        Call MTHEMATIC_JELP(Ammm,U,EPH,ESN,ECN,EDN)

      chana = 3
      
       FactorU=(-1.25d0+1.5d0*AMINNT/d-0.25d0*AMINNT**2d0/d**2d0
     1+(1.5d0*H/d-0.5d0*H*AMINNT/d**2d0)*CN**2d0-0.25d0*H**2d0/d**2d0*CN**4d0
     1-8d0*H*AKKK**2d0/ALAMDA**2d0*(d/3d0-0.5d0*(abs(ZZ))**2d0/d)*(-Ammm*(SSN**2d0)*(CN**2d0)+(CN**2d0)*(ddn**2d0)
     1-(SSN**2d0)*ddn**2d0))
       
       AFACTOR1 = sqrt(g*d)
       
       AFACTOR2 = -1.25d0+1.5d0*AMINNT/d -0.25d0*AMINNT**2d0/d**2d0
     1+(1.5d0*H/d-0.5d0*H*AMINNT/d**2d0)*CN**2d0
       
       AFACTOR3 = -8d0*H*AKKK**2d0/ALAMDA**2d0*(d/3d0-0.5d0*(abs(ZZ))**2d0/d)

      UUUXX=(sqrt(g*d))*FactorU
 
      FactorV=ZZ*2d0*H*AKKK/ALAMDA/d*(1d0+ht/d+H/d*(CN**2d0)+32d0/3d0*(AKKK**2d0)/(ALAMDA**2)*(d**2d0-0.5*d0*(ZZ**2))*(Ammm*
     1(SSN**2d0)-Ammm*(CN**2d0)-(ddn**2d0)))*SSN *CN*ddn
    
      VVVXX=   (sqrt(g*d))*FactorV  
    
    
      FactorAX= (4d0*H*AKKK/(T*d))*((1.5d0-(ht/(2d0*d)))-(0.5d0*H*(CN**2d0)/d)+(16.0d0*(AKKK**2d0)/(ALAMDA**2d0))*
     1(((d**2.0d0)/3.0d0)-(zz**2.0d0))*(Ammm*(SSN**2d0)-Ammm*(CN**2d0)-ddn**2d0))*SSN*CN*ddn
      AAAXX =  (sqrt(g*d))*FactorAX  
    

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    
        FactorAZ= 4d0*H*(AKKK**2d0)/ALAMDA/T/d*((1d0+ht/d)*((SSN**2d0)*(ddn**2d0)-(CN**2d0)*(ddn**2d0)+Ammm*(SSN**2d0)*(ddn**2d0))
     1+H/d *(3d0*(SSN**2d0)*(ddn**2d0)-(CN**2d0)*(ddn**2d0)+Ammm*(SSN**2d0))*(CN**2d0)-32d0/3d0*(AKKK**2d0)/(ALAMDA**2d0)*
     1(d**2d0-0.5d0*ZZ) 
     1*(9d0*Ammm*(SSN**2d0)*(CN**2d0)*(ddn**2d0)-Ammm*(SSN**4d0)*(Ammm*(CN**2d0)+(ddn**2d0)) 
     1+Ammm*(CN**4d0)*(Ammm*(SSN**2d0)+(ddn**2d0))+(ddn**4d0)*((SSN**2d0)-(CN**2d0)) ))
  
      AAAZZ =  (sqrt(g*d))*FactorAZ 
      
      RETURN
      End SUBROUTINE
      
 !***************************************************************************************************************************************************************
    
      SUBROUTINE Cnoidal_Wave_Surface  (H,waterdepth,g,T,xx,ZZ,times,wavenumber,AKKK,ALAMDA,Ammm,AEEE,CN,Suface_heignt)
        
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
    
      Dimension AKK(13,1),AMM(13,1),AEE(13,1),AHL(13,1),ALAM(13,1),OMEGA(13,1)
      Dimension AXX(13,1),EERROR(13,1),WaveK(13,1)
      Dimension AM00(6,1),AM02(6,1),AM04(6,1),AM06(6,1),AM08(6,1),AM10(6,1)
      Dimension AM12(6,1),AM14(6,1),AM16(6,1),AM18(6,1),AM20(6,1) 
    
      d=waterdepth
    
!    CALL Cnoidal_wave_Length (H,waterdepth,g,T,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
    
      ZZETAA=(wavenumber*xx-2d0*AKKK/T*times)
    
    
      WAVE_trough=-H*(AKKK*(1-Ammm)-AEEE)/Ammm/AKKK 
    
    !H*(AKKK*(1d0-Ammm)-AEEE)/(Ammm*AKKK)
      CALL Jacobian_elliptic_function (ZZETAA,Ammm,CN,CN2,SSN,ddn)
    
!    write(*,7)ZZETAA,Ammm,CN,SSN,ddn
7     Format(F12.4,F12.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)

      peak = H+waterdepth-abs(WAVE_trough)
      ht   = peak-H  
    
    
      Suface_heignt = ht+H*(CN**2d0)
    
      RETURN
      END SUBROUTINE
      
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
      
      SUBROUTINE Cnoidal_Wave_Crest  (H,waterdepth,g,T,TTime,Wave_Crest)
        
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
    
      Dimension AKK(13,1),AMM(13,1),AEE(13,1),AHL(13,1),ALAM(13,1),OMEGA(13,1)
      Dimension AXX(13,1),EERROR(13,1),WaveK(13,1)
      Dimension AM00(6,1),AM02(6,1),AM04(6,1),AM06(6,1),AM08(6,1),AM10(6,1)
      Dimension AM12(6,1),AM14(6,1),AM16(6,1),AM18(6,1),AM20(6,1) 
      
      call Cnoidal_wave_Length (H,waterdepth,g,T,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)

      ZZETAA=2d0*AKKK/T*TTime
      call Jacobian_elliptic_function (ZZETAA,Ammm,CN,CN2,SSN,ddn)
      WAVE_trough=-H*(AKKK*(1-Ammm)-AEEE)/Ammm/AKKK 
      Wave_Crest = (H-abs(WAVE_trough))  *1.11
      
      TEST = Wave_Crest + waterdepth

      RETURN
      END SUBROUTINE
    
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************