      SUBROUTINE CNOIEAL_CALL_MODULOUS(WH,WDEP,GRAVM,AMODOUS)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      Call    ELLIPTIC_INTERGAL(AMmodulous,Am,AK)

      RETURN
      END SUBROUTINE
C	=======================================================================
C	======================================================================= 
      SUBROUTINE CNOIDALSURFACEELEVPLOT(H,d,PEROID,Wave_Crest,AModolous)
      
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      WNT= Wave_Crest-H
      
      ai=1.0D0
      
      EPH = 0.0
      
      DO WHILE ( EPH .LE. 360.0D0 )
          
          UUU =  (ai-1.0D0)/100.0D0
          
      CALL MTHEMATIC_JELP(AModolous,UUU,EPH,ESN,ECN,EDN)
          
      WSURFACE = WNT +  H*(ECN**2) + d
      !WAIRY = (H/2)*COS(EPH*2*(22/7)/180.0D0)+d
      
      EPH = EPH*2.0D0
      
      WRITE(5040,1)UUU,EPH,WSURFACE!,WAIRY

1     FORMAT(10F12.5)
      
      
      
      ai = ai + 1.0D0
      
          
      ENDDO    
      
      
      
      
      
      
      END SUBROUTINE
C	=======================================================================
C	======================================================================= 