      Subroutine Mlife_FAST_Fatigue_Analysis
      !USE IFPORT
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
    
      write(*,*) 'Hello Mlife'
      
      
      !MAKEDIR = MAKEDIRQQ('DataMLife') 
      !OPEN(UNIT=189654  ,FILE='DataMLife\MLife_input.mlif'       ,STATUS='UNKNOWN'    )
      
      !Call system ('"C:\Users\r751j\Documents\MATLAB\Fortran Mlife\MalifeFortran\MalifeFortran\MLife1.exe" <IN.DAT')
      !CALL EXECUTE_COMMAND_LINE(COMMAND [, WAIT, EXITSTAT, CMDSTAT, CMDMSG ]) 
      
      CALL WriteMLifeInput
      
      CALL WriteMLifeLoadResopnseInput

      
      !CALL system ("MLife1.exe")
      Call Execute_Command_Line ("MLife1.exe")
      
      
      chana = 3
      
 !     CALL READMLifeOUTPUT!(TOTAL_DAMAGE,Life_Time)
      
      write(*,*) 'End of MLife Analysis'
      
      !STOP
      
      RETURN
      END SUBROUTINE
!     ===============================================================================================================================================================================