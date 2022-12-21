      SUBROUTINE READMLifeOUTPUT!(TOTAL_DAMAGE,Life_Time)
      
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      !====================================================
      ! FatigueMlife_Lifetime_Damage.txt'
      ALLOCATABLE ALifetime_Damage(:)
      ALLOCATABLE Time_Until_Failure(:)
      ALLOCATABLE ALifetime_Damage_without_Goodman(:)
      ALLOCATABLE Time_Until_Failure_without_Goodman(:)
      !====================================================
      ! FatigueMlife_Lifetime_DELs.txt
      ALLOCATABLE ALifetime_DEL(:)
      ALLOCATABLE ALifetime_DELs_ZeroMean(:)
      ALLOCATABLE ALifetime_DELs_Without_Goodman(:)
      !====================================================
      
      

      
      
      NUMBER_OF_EXPONENT = 2
      
      
      !===========================================================================================
      ! Data_Statistics
      !===========================================================================================
!      CALL READMLifeOUTPUT_STATIC(DataMinimum,DataMean,DataMeanMaximum,DataMeanStdDev,
!     1 DataMeanSkewness,DataMeanKurtosis,DataMeanRange)
      
      
      !===========================================================================================
      ! FatigueMlife_Lifetime_Damage
      !===========================================================================================
      ALLOCATE(ALifetime_Damage(NUMBER_OF_EXPONENT))
      ALLOCATE(Time_Until_Failure(NUMBER_OF_EXPONENT))
      ALLOCATE(ALifetime_Damage_without_Goodman(NUMBER_OF_EXPONENT))
      ALLOCATE(Time_Until_Failure_without_Goodman(NUMBER_OF_EXPONENT))
      
      
      CALL READMLifeOUTPUT_Lifetime_Damage(NUMBER_OF_EXPONENT,ALifetime_Damage,
     1 Time_Until_Failure,ALifetime_Damage_without_Goodman,Time_Until_Failure_without_Goodman)
      
      !===========================================================================================
      ! FatigueMlife_Lifetime_DELs
      !===========================================================================================
      ALLOCATE (ALifetime_DEL(NUMBER_OF_EXPONENT)                  )
      ALLOCATE (ALifetime_DELs_ZeroMean(NUMBER_OF_EXPONENT)        )
      ALLOCATE (ALifetime_DELs_Without_Goodman(NUMBER_OF_EXPONENT) )
      
      
      
      
      CALL READMLifeOUTPUT_Lifetime_DEL(NUMBER_OF_EXPONENT,ALifetime_DEL,
     1  ALifetime_DELs_ZeroMean,ALifetime_DELs_Without_Goodman)
      
      !===========================================================================================
      ! FatigueMlife_Short-term_Damage_Rate
      !===========================================================================================
      
      
      
      
      
      
      
      

      
      
      DEALLOCATE(ALifetime_Damage,Time_Until_Failure,ALifetime_Damage_without_Goodman,Time_Until_Failure_without_Goodman)
      DEALLOCATE(ALifetime_DEL,ALifetime_DELs_ZeroMean,ALifetime_DELs_Without_Goodman)
      
      RETURN
      END SUBROUTINE
!     ====================================================================================================================================      
      SUBROUTINE READMLifeOUTPUT_STATIC(DataMinimum,DataMean,DataMeanMaximum,DataMeanStdDev,
     1 DataMeanSkewness,DataMeanKurtosis,DataMeanRange)
      
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*200   NAME
      
      
      
       WRITE(*,*)"READ MLife ANALYSIS RESULTS"
       
       chana = 3
       
      REWIND(5046)
      
      READ(5046,*)!1
      READ(5046,*)!2
      READ(5046,*)!3
      READ(5046,*)!4
      READ(5046,*)!5
      READ(5046,*)!6
      
      !READ
      !DataMinimum,DataMean,DataMeanMaximum,DataMeanStdDev,DataMeanSkewness,DataMeanKurtosis,DataMeanRange
      
      READ(5046,*)  NAME ,DataMinimum,DataMean,DataMeanMaximum,DataMeanStdDev,DataMeanSkewness,DataMeanKurtosis,DataMeanRange
      
      
      chana = 3
      
      
      
      
      RETURN
      END SUBROUTINE
!     ====================================================================================================================================       
      SUBROUTINE READMLifeOUTPUT_Lifetime_Damage(NUMBER_OF_EXPONENT,ALifetime_Damage,
     1 Time_Until_Failure,ALifetime_Damage_without_Goodman,Time_Until_Failure_without_Goodman)
      
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*200   NAME
      
      DIMENSION ALifetime_Damage(NUMBER_OF_EXPONENT)
      DIMENSION Wohler_EXPONENT(NUMBER_OF_EXPONENT)
      DIMENSION Time_Until_Failure(NUMBER_OF_EXPONENT)
      DIMENSION ALifetime_Damage_without_Goodman(NUMBER_OF_EXPONENT)
      DIMENSION Time_Until_Failure_without_Goodman(NUMBER_OF_EXPONENT)
      
      
      REWIND(5047)
      
      NHEADER = 20
      
      DO I = 1,NHEADER
      READ(5047,*)! HEADDER LINE  
      ENDDO
      
      
      READ(5047,*)NAME,LUlT
      DO I =1,NUMBER_OF_EXPONENT      
      READ(5047,*)Wohler_EXPONENT(I),ALifetime_Damage(I)
      ENDDO
      
      DO I = 1,5
      READ(5047,*)! HEADDER LINE  
      ENDDO
      
      READ(5047,*)NAME,LUlT
      DO I =1,NUMBER_OF_EXPONENT      
      READ(5047,*)Wohler_EXPONENT(I),Time_Until_Failure(I)
      ENDDO
      
      DO I = 1,5
      READ(5047,*)! HEADDER LINE  
      ENDDO
      
      READ(5047,*)NAME,LUlT
      DO I =1,NUMBER_OF_EXPONENT      
      READ(5047,*)Wohler_EXPONENT(I),ALifetime_Damage_without_Goodman(I)
      ENDDO
      
      DO I = 1,5
      READ(5047,*)! HEADDER LINE  
      ENDDO
      
      READ(5047,*)NAME,LUlT
      DO I =1,NUMBER_OF_EXPONENT      
      READ(5047,*)Wohler_EXPONENT(I),Time_Until_Failure_without_Goodman(I)
      ENDDO

      RETURN
      END SUBROUTINE
!     ====================================================================================================================================  
!     ====================================================================================================================================       
      SUBROUTINE READMLifeOUTPUT_Lifetime_DEL(NUMBER_OF_EXPONENT,ALifetime_DEL,
     1 ALifetime_DELs_ZeroMean,ALifetime_DELs_Without_Goodman)
      
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*200   NAME
      
      DIMENSION ALifetime_DEL(NUMBER_OF_EXPONENT)
      DIMENSION Wohler_EXPONENT(NUMBER_OF_EXPONENT)
      DIMENSION ALifetime_DELs_ZeroMean(NUMBER_OF_EXPONENT)
      DIMENSION ALifetime_DELs_Without_Goodman(NUMBER_OF_EXPONENT)
      !DIMENSION Time_Until_Failure_without_Goodman(NUMBER_OF_EXPONENT)
      
      
      ! FatigueMlife_Lifetime_DELs.txt
      
      REWIND(5048)
      
      NHEADER = 28
      
      DO I = 1,NHEADER
      READ(5048,*)! HEADDER LINE  
      ENDDO
      
      !'FATIGUE'  Lifetime DELs at Fixed Mean for various S/N Curves
      
      READ(5048,*)NAME,LUlT
      READ(5048,*)NAME,AL_MF_FixMean
      DO I =1,NUMBER_OF_EXPONENT      
      READ(5048,*)Wohler_EXPONENT(I),ALifetime_DEL(I)
      ENDDO
      
      NHEADER = 5
      DO I = 1,NHEADER
      READ(5048,*)! HEADDER LINE  
      ENDDO
      
      !'FATIGUE'  Lifetime DELs at Zero Mean for various S/N Curves
      
      READ(5048,*)NAME,LUlT
      READ(5048,*)NAME,AL_MF_zeroMean
      DO I =1,NUMBER_OF_EXPONENT      
      READ(5048,*)Wohler_EXPONENT(I),ALifetime_DELs_ZeroMean(I)
      ENDDO
      
      
      NHEADER = 5
      DO I = 1,NHEADER
      READ(5048,*)! HEADDER LINE  
      ENDDO
      
      !'FATIGUE'  Lifetime DELs without Goodman Correction for various S/N Curves
      
      READ(5048,*)NAME,LUlT
      READ(5048,*)NAME,AL_MF_Without_Goodman
      DO I =1,NUMBER_OF_EXPONENT      
      READ(5048,*)Wohler_EXPONENT(I),ALifetime_DELs_Without_Goodman(I)
      ENDDO

      RETURN
      END SUBROUTINE
!     ====================================================================================================================================
!     ====================================================================================================================================       
!      SUBROUTINE READMLifeOUTPUT_Short-term_Damage_Rate(NUMBER_OF_EXPONENT,ALifetime_DEL,
!     1 ALifetime_DELs_ZeroMean,ALifetime_DELs_Without_Goodman)
!      
!      IMPLICIT INTEGER*4 (I-N)
!      IMPLICIT REAL*8 (A-H,O-Z)
!      CHARACTER*200   NAME
!      
!      DIMENSION ALifetime_DEL(NUMBER_OF_EXPONENT)
!      DIMENSION Wohler_EXPONENT(NUMBER_OF_EXPONENT)
!      DIMENSION ALifetime_DELs_ZeroMean(NUMBER_OF_EXPONENT)
!      DIMENSION ALifetime_DELs_Without_Goodman(NUMBER_OF_EXPONENT)
!      !DIMENSION Time_Until_Failure_without_Goodman(NUMBER_OF_EXPONENT)
!      
!      
!      ! FatigueMlife_Lifetime_DELs.txt
!      
!      REWIND(5049)
!      
!      NHEADER = 6
!      
!      DO I = 1,NHEADER
!      READ(5049,*)! HEADDER LINE  
!      ENDDO
!      
!      !'FATIGUE'  Short-term Damage-rate (-/s) for various S/N Curves
!      
!      READ(5048,*)NAME,LUlT
!      READ(5048,*)NAME,AL_MF_FixMean
!      DO I =1,NUMBER_OF_EXPONENT      
!      READ(5048,*)Wohler_EXPONENT(I),ALifetime_DEL(I)
!      ENDDO
!      
!      NHEADER = 5
!      DO I = 1,NHEADER
!      READ(5048,*)! HEADDER LINE  
!      ENDDO
!      
!      !'FATIGUE'  Lifetime DELs at Zero Mean for various S/N Curves
!      
!      READ(5048,*)NAME,LUlT
!      READ(5048,*)NAME,AL_MF_zeroMean
!      DO I =1,NUMBER_OF_EXPONENT      
!      READ(5048,*)Wohler_EXPONENT(I),ALifetime_DELs_ZeroMean(I)
!      ENDDO
!      
!      
!      NHEADER = 5
!      DO I = 1,NHEADER
!      READ(5048,*)! HEADDER LINE  
!      ENDDO
!      
!      !'FATIGUE'  Lifetime DELs without Goodman Correction for various S/N Curves
!      
!      READ(5048,*)NAME,LUlT
!      READ(5048,*)NAME,AL_MF_Without_Goodman
!      DO I =1,NUMBER_OF_EXPONENT      
!      READ(5048,*)Wohler_EXPONENT(I),ALifetime_DELs_Without_Goodman(I)
!      ENDDO
!
!      RETURN
!      END SUBROUTINE
!     ==================================================================================================================================== 