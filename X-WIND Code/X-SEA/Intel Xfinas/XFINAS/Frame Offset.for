      SUBROUTINE Reduce_Offset_member (XYZ,XR,YR,ZR,MLE,WFX,WFY,WFZ
     1           ,IELETPZ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      
      
      !------------------------------ sorn mark wave force change    --------------------------
!          WAVEFORCEXELE is the final wave force (begin + end)
!         1. check location on offset or not
!         2. if on location set force = 0
         
       COMMON A(9000000),IA(9000000)
       COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN
	  
	  !sorn added 
       COMMON /sornOffset/ NELO,Nmemoffset(1000),NmemOffsetSelect(1000) ,Nmemtrueelement(1000),GloLoSw(1000) !total number of off set node
       
       COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,                      
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,                  
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT 
       
       
       DIMENSION XYZ(NCO*NNM,NELE)
       Dimension sor(12)
      
       
       !if local axis must tranform to global first
       
                       
        ! offset value
           
        numoffset = NELO
        sorn = 3
        NumLoop = 6*NELO
        
   !     Nmemoffset
   !     NmemOffsetSelect
        
        do 333 iNumOffcase = 1,NELO
           if(MLE == Nmemoffset(iNumOffcase)) then
               iSelOffcase = NmemOffsetSelect(iNumOffcase)
               exit
           endif
            
333   continue
 !
        
      
        
        ! offset mark i
         Rmarkix =  XYZ(1,MLE) +       A(LOP+(iSelOffcase-1)*6)
         Rmarkiy =  XYZ(2,MLE) +       A(LOP+(iSelOffcase-1)*6 + 1)
         Rmarkiz =  XYZ(3,MLE) +       A(LOP+(iSelOffcase-1)*6 + 2)
        
        ! offset mark j
        
         Rmarkjx =  XYZ(4,MLE) +        A(LOP+(iSelOffcase-1)*6 + 3)
         Rmarkjy =  XYZ(5,MLE) +        A(LOP+(iSelOffcase-1)*6 + 4)
         Rmarkjz =  XYZ(6,MLE) +        A(LOP+(iSelOffcase-1)*6 + 5)
        
        Soffvalue1  =  A(LOP+(iSelOffcase-1)*6 + 0)
        Soffvalue2  =  A(LOP+(iSelOffcase-1)*6 + 1)
        Soffvalue3  =  A(LOP+(iSelOffcase-1)*6 + 2)
        Soffvalue4  =  A(LOP+(iSelOffcase-1)*6 + 3)
        Soffvalue5  =  A(LOP+(iSelOffcase-1)*6 + 4)
        Soffvalue6  =  A(LOP+(iSelOffcase-1)*6 + 5)
      
      
       ! rearrang 
       
      if (XYZ(1,MLE).GT.Rmarkix) then
      
      RMark2Xi = XYZ(1,MLE)
      RMark1Xi = Rmarkix
      
      else 
      
      RMark1Xi = XYZ(1,MLE)
      RMark2Xi = Rmarkix   
          
      endif
      
      if (XYZ(2,MLE).GT.Rmarkiy) then
      
      RMark2Yi = XYZ(2,MLE)
      RMark1Yi = Rmarkiy
      
      else 
      
      RMark1Yi = XYZ(2,MLE)
      RMark2Yi = Rmarkiy   
          
      endif
      
      if (XYZ(3,MLE).GT.Rmarkiz) then
      
      RMark2Zi = XYZ(3,MLE)
      RMark1Zi = Rmarkiz
      
      else 
      
      RMark1Zi = XYZ(3,MLE)
      RMark2Zi = Rmarkiz   
          
      endif
  ! ---------------------------------------------------------------!    
       if (XYZ(4,MLE).GT.Rmarkjx) then
      
      RMark2Xj = XYZ(4,MLE)
      RMark1Xj = Rmarkjx
      
      else 
      
      RMark1Xj = XYZ(4,MLE)
      RMark2Xj = Rmarkjx   
          
      endif
      
      if (XYZ(5,MLE).GT.Rmarkjy) then
      
      RMark2Yj = XYZ(5,MLE)
      RMark1Yj = Rmarkjy
      
      else 
      
      RMark1Yj = XYZ(5,MLE)
      RMark2Yj = Rmarkjy   
          
      endif
      
      if (XYZ(6,MLE).GT.Rmarkjz) then
      
      RMark2Zj = XYZ(6,MLE)
      RMark1Zj = Rmarkjz
      
      else 
      
      RMark1Zj = XYZ(6,MLE)
      RMark2Zj = Rmarkjz   
          
      endif
 

         
        SornOriX1 = XYZ(1,MLE)
        SornOriX2 = XYZ(4,MLE)
        SornOriY1 = XYZ(2,MLE) 
        SornOriY2 = XYZ(5,MLE)
        SornOriZ1 = XYZ(3,MLE)
        SornOriZ2 = XYZ(6,MLE)
         
           
       ! key
        ijkey = 0
        ioffsetkey = 0
        Soffsetkey = 0
        
       if (XR == Rmarkix) then
           if  (YR == Rmarkiy) then 
               if (ZR == Rmarkiz) then
                   ijkey = 1
                   endif
               endif
       elseif (XR == Rmarkjx) then
           if (YR == Rmarkjy) then
               if (ZR == Rmarkjz) then 
                   ijkey = 1
                   sorn =4
               endif
           endif
      endif
       
       if ( (XR >= RMark1Xi) .and. (XR <= RMark2Xi)  .and. (YR >= RMark1Yi) .and. (YR <= RMark2Yi).and. (ZR >=   RMark1Zi)                         ! offset 1
     1          .and. (ZR <=  RMark2Zi ) .and. (ijkey /= 1) )
     1  then 
       Soffsetkey = 1 
        
      else if ((XR >= RMark1Xj).and.(XR <= RMark2Xj).and.(YR >= RMark1Yj).and.(YR <=  RMark2Yj).and.(ZR >=  RMark1Zj)                                 ! offset 1
     1         .and.(ZR <=  RMark2Zj).and. (ijkey /= 1))
     1 then 
      Soffsetkey = 1 
      endif
        
      if (Soffsetkey == 1) then

           WFX = 0
           WFY = 0 
           WFZ = 0   
           
      endif
           
1414      ss = 0
      
      RETURN
      END
  