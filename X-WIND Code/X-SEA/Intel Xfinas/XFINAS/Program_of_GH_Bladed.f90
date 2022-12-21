      SUBROUTINE GH_BLADED_FATIGUE_ANALYSIS_PROGRAM
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      ! THE GH_BLADED.EXE PROGRAM WAS USED FOR FATIGUE ANALYSIS
      ! PLEASE SEE THE DETAIL OF CODE IN GH_BLADED PROGRAM FILES
      Call Execute_Command_Line ("GH_bladed_XSEA_program.exe")

      ! -------------------------------GH BLADED-----------------------------
      ! TO OPERATED GH BLADED   , USER SHOULD INPUT XOPT = 4 IN "XFASTOPT.DAT" 
      ! THE PROGRAM IS REQUIRED SACS.IN AND SOME FILES, WHICH IS REQUIRE IN SACS PROGRAM
      ! ----------------------------------------------------------------------
      
      ! -------------------------------X-SEA FATIGUE -------------------------
      ! TO OPERATED X-SEA   , USER SHOULD INPUT XOPT = 4 IN "XFASTOPT.DAT" 
      ! THE REQUIREMENT OF FILE HAS A SERVERL FILES, PLEASE USING GID FOR RUNNING 
      ! ----------------------------------------------------------------------
      
      
      RETURN
      END SUBROUTINE
      !========================================================================================
      
      SUBROUTINE FAST_FLS
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      ! THE GH_BLADED.EXE PROGRAM WAS USED FOR FATIGUE ANALYSIS
      ! PLEASE SEE THE DETAIL OF CODE IN GH_BLADED PROGRAM FILES
      Call Execute_Command_Line ("Tranfer_FAST_Force.exe")
      Call Execute_Command_Line ("GH_bladed_XSEA_program.exe")

      ! -------------------------------GH BLADED-----------------------------
      ! TO OPERATED GH BLADED   , USER SHOULD INPUT XOPT = 4 IN "XFASTOPT.DAT" 
      ! THE PROGRAM IS REQUIRED SACS.IN AND SOME FILES, WHICH IS REQUIRE IN SACS PROGRAM
      ! ----------------------------------------------------------------------
      
      ! -------------------------------X-SEA FATIGUE -------------------------
      ! TO OPERATED X-SEA   , USER SHOULD INPUT XOPT = 4 IN "XFASTOPT.DAT" 
      ! THE REQUIREMENT OF FILE HAS A SERVERL FILES, PLEASE USING GID FOR RUNNING 
      ! ----------------------------------------------------------------------
      
      
      RETURN
      END SUBROUTINE
      !========================================================================================