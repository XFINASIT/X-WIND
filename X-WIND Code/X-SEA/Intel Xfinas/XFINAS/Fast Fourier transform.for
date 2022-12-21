      SUBROUTINE FAST_FOURIER_FORWARD (DATA_INPUT,DATA_OUTPUT,INDEX)
      USE MKL_DFTI
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
     
      !Complex :: X_in(32)
      !Complex :: X_out(32)
      REAL :: DATA_INPUT_REAL(INDEX)
      REAL :: DATA_OUTPUT_REAL(INDEX)
      COMPLEX :: DATA_INPUT_COMPLEX(INDEX)
      COMPLEX :: DATA_OUTPUT_COMPLEX(INDEX)
      type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
      Integer :: Status
      
      DATA_INPUT_COMPLEX = DATA_INPUT_REAL
      
      
      ! ------------------- ONE-DIMENSION -------------------
      !...put input data into X_in(1),...,X_in(32); Y_in(1),...,Y_in(32)
      ! Perform a complex to complex transform
      !Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_SINGLE,DFTI_COMPLEX , 1, 32 )
      !Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      !Status = DftiCommitDescriptor( My_Desc1_Handle)
      !Status = DftiComputeForward( My_Desc1_Handle, X_in, X_out )
      !Status = DftiFreeDescriptor(My_Desc1_Handle)
      ! result is given by {X_out(1),X_out(2),...,X_out(32)}
      ! Perform a real to complex conjugate-even transform
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INDEX)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeForward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeForward
      !Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)
      
      ! SELECTED REALVALUE FORM COMPLEX
      DATA_OUTPUT_REAL = DATA_OUTPUT_COMPLEX
      ! ------------------- ONE-DIMENSION -------------------
      END
! ======================================================================================
      SUBROUTINE FAST_FOURIER_FORWARD_TEST (INDEX,SPM_TEST)
      USE MKL_DFTI
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
     
      !Complex :: X_in(32)
      !Complex :: X_out(32)
      REAL :: DATA_INPUT_REAL(INDEX)
      REAL :: DATA_OUTPUT_REAL(INDEX)
      COMPLEX :: DATA_INPUT_COMPLEX(INDEX)
      COMPLEX :: DATA_OUTPUT_COMPLEX(INDEX)
      type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
      Integer :: Status
      DIMENSION SPM_TEST(INDEX)
      
      !OPEN(UNIT=11111,FILE='Read_Data.dat'      ,STATUS='UNKNOWN')
      !
      !DO I = 1,INDEX
      !READ (11111,*) TIME,DATA_INPUT_REAL(I)
      !ENDDO
      
      !DATA_INPUT_COMPLEX = DATA_INPUT_REAL
      DATA_INPUT_COMPLEX = SPM_TEST
      
      
      ! ------------------- ONE-DIMENSION -------------------
      !...put input data into X_in(1),...,X_in(32); Y_in(1),...,Y_in(32)
      ! Perform a complex to complex transform
      !Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_SINGLE,DFTI_COMPLEX , 1, 32 )
      !Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      !Status = DftiCommitDescriptor( My_Desc1_Handle)
      !Status = DftiComputeForward( My_Desc1_Handle, X_in, X_out )
      !Status = DftiFreeDescriptor(My_Desc1_Handle)
      ! result is given by {X_out(1),X_out(2),...,X_out(32)}
      ! Perform a real to complex conjugate-even transform
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INDEX)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      !Status = DftiComputeForward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeForward
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      !Status = DftiFreeDescriptor(My_Desc2_Handle)
      
      ! SELECTED REALVALUE FORM COMPLEX
      DATA_OUTPUT_REAL = DATA_OUTPUT_COMPLEX
      ! ------------------- ONE-DIMENSION -------------------
      END
! ======================================================================================
      SUBROUTINE FAST_FOURIER_BACKWARD_DISPLACEMENT (INDEX,NEQ)!(DATA_INPUT_REAL,DATA_OUTPUT_REAL,INDEX)
      USE MKL_DFTI
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 ASIGN 
      
      !Complex :: X_in(32)
      !Complex :: X_out(32)
      REAL :: DATA_INPUT_REAL(INDEX)
      REAL :: DATA_OUTPUT_REAL(INDEX,NEQ)
      COMPLEX :: DATA_INPUT_COMPLEX(INDEX)
      COMPLEX :: DATA_COMPLEX(INDEX)
      COMPLEX :: DATA_OUTPUT_COMPLEX(INDEX)
      type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
      Integer :: Status
      COMMON /TIME/ DDT,CTIM,NINC
      
      DATA_INPUT_COMPLEX = 0.
      DATA_OUTPUT_REAL   = 0.
      REWIND (253)  
      READ (253,*)
      DO J = 1,NEQ
      DO I = 1,NINC/2D0
      READ (253,*) DUM,DUM,DATA_INPUT_COMPLEX(I)
      ENDDO
      DATA_INPUT_COMPLEX = DATA_INPUT_COMPLEX
      ! TESTING FILES
      !REWIND (256)
      !DO JJ = 1,NINC
      !READ (256,*) DATA_REAL,DATA_IMG,ASIGN
      !IF (TRIM(ASIGN).EQ."+")THEN
      !DATA_INPUT_COMPLEX(JJ) = CMPLX(DATA_REAL,+DATA_IMG)/NINC
      !ELSEIF (TRIM(ASIGN).EQ."-")THEN
      !DATA_INPUT_COMPLEX(JJ) = CMPLX(DATA_REAL,-DATA_IMG)/NINC
      !ENDIF
      !ENDDO
      
      ! WITHOUT COMPLEXS
      !REWIND (256)
      !DO JJ = 1,NINC
      !    READ (256,*) DATA_INPUT_REAL(JJ)
      !    DATA_INPUT_COMPLEX(JJ) = DATA_INPUT_REAL(JJ)!/NINC
      !ENDDO
      
      ! FILTER RESULT
      DO KK = 1,NINC/2D0  
      IF (KK.GT.1)THEN
          ! FOR REAL VALUE
          DATA_INPUT_COMPLEX(NINC+2-KK) = CONJG(DATA_INPUT_COMPLEX(KK))
          ! FOR COMLEX VALUE
          !DATA_INPUT_COMPLEX(NINC+2-KK)  = CONJG(DATA_INPUT_COMPLEX(KK))
      ENDIF
         
      ENDDO

      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INDEX)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)

      ! SELECTED REAL VALUE FORM COMPLEX
       DATA_OUTPUT_REAL(1:NINC,J) = DATA_OUTPUT_COMPLEX
       DO K = 1,NINC
       WRITE (254,100) DATA_OUTPUT_REAL(K,J)
       ENDDO
100    FORMAT (E25.10)
       ENDDO
      
      !ENDDO
    
      ! ------------------- ONE-DIMENSION -------------------
      END
C	===================================================================== 
      SUBROUTINE FAST_FOURIER_BACKWARD_FORCE (FORCE_FREQ,INDEX,NINC)
      USE MKL_DFTI
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      REAL :: DATA_INPUT_REALfREPLACE_RESULT_FORCE(NINC)
      REAL :: DATA_OUTPUT_REAL(NINC)
      COMPLEX :: DATA_INPUT_COMPLEX(NINC)
      COMPLEX :: DATA_OUTPUT_COMPLEX(NINC)
      type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
      Integer :: Status
      DIMENSION FORCE_FREQ(NINC)
      DATA_OUTPUT_COMPLEX = 0.
      
      DO II = 1,NINC
      DATA_INPUT_COMPLEX(II) = FORCE_FREQ(II)
      ENDDO
      
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, NINC)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) 
      Status = DftiFreeDescriptor(My_Desc2_Handle)
      
       ! SELECTED REALVALUE FORM COMPLEX
       DO II = 1,NINC
       DATA_OUTPUT_REAL(II) = DATA_OUTPUT_COMPLEX(II)
       ENDDO
       
       DO K = 1,NINC
       WRITE (255,100) DATA_OUTPUT_REAL(K)
       ENDDO
100    FORMAT (E25.10)
      
      END
! ======================================================================================
      SUBROUTINE TRANFORM (NELE)
      ! THIS SUBROUTINE IS USED FOR TRANFORM THE VALUE OF THE FREQUENCY TO TIME DOMAIN
      ! CAN BE TRANFORM ONLY INTERNAL FORCE FOR FULLFILL FATIGUE ANALYSIS
      ! IN ADDITION, REACTION FORCE AND MOMENT WERE EXECUTED. 
      USE MKL_DFTI ! OPEN MKL FOR FFT, ***REQUIRE FILE***
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*200 NAME
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /TIME/ DDT,CTIM,NINC
      DIMENSION FX_OUT(NINC),FY_OUT(NINC),FZ_OUT(NINC),AMX_OUT(NINC),AMY_OUT(NINC),AMZ_OUT(NINC)
      
      COMMON /LHEADER/ LHDSTD(500),LHDCOM(500)
      
      COMMON /BUNDARY_CONDITION/ INDEXBUNDARY(9999,6),NODEABUNDARY(9999)
      
      DIMENSION FX(NINC),FY(NINC),FZ(NINC),AMX(NINC),AMY(NINC),AMZ(NINC)
      DIMENSION FX_RESULT(5),FY_RESULT(5),FZ_RESULT(5) 
      DIMENSION AMX_RESULT(5),AMY_RESULT(5),AMZ_RESULT(5)
      
      DIMENSION FORCEX(NINC),FORCEY(NINC),FORCEZ(NINC)
      DIMENSION FORCE_X_OUT(NINC),FORCE_Y_OUT(NINC),FORCE_Z_OUT(NINC)
      
      DIMENSION AMOMENT_X(NINC),AMOMENT_Y(NINC),AMOMENT_Z(NINC)
      DIMENSION AMOMENT_X_OUT(NINC),AMOMENT_Y_OUT(NINC),AMOMENT_Z_OUT(NINC)
      
      REAL :: DATA_INPUT_REAL(NINC)
      
      REAL :: DATA_OUTPUT_REAL(NINC)
      COMPLEX :: DATA_INPUT_COMPLEX(NINC)
      COMPLEX :: DATA_OUTPUT_COMPLEX(NINC)
      type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
      Integer :: Status
      
      ! INTIAL SETTING 
      FX   = 0.
      FY   = 0.
      FZ   = 0.
      AMX  = 0.
      AMY  = 0.
      AMZ  = 0.
      
      ! TRANFORM DATA FOR INPUT IN FFT FUNCTION (INTEL MKL)
      WRITE(5201,90)
90    FORMAT ("NO.ELE     POINT       TIME_STEP     FX      FY      FZ      MX      MY      MZ")
      
      DO J = 1,NELE ! TOTAL NUMBER OF ELEMENT
      DO K = 1,5
      INDEX = 1
      REWIND (5200)
      DO I = 1,NELE*NINC*5
         READ (5200,*) INCC,IFRAMESTORE,IG
         IF (IFRAMESTORE.EQ.J.AND.IG.EQ.K)THEN
         READ (5200,*) FX(INDEX),FY(INDEX),FZ(INDEX),AMX(INDEX),AMY(INDEX),AMZ(INDEX)
         INDEX = INDEX +1
         ELSEIF (IFRAMESTORE.NE.J.OR.IG.NE.K)THEN
         READ (5200,*)
         ENDIF
      ENDDO
      
      !  FFT FREQ. > TIME.
      !  ============= FX =============
      DO I =1,NINC
      DATA_INPUT_COMPLEX(I) = FX(I)
      ENDDO
      
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)

      DO I =1,NINC
      FX_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
      ENDDO
      !  ============= FX =============
      
      !  ============= FY =============
      DO I =1,NINC
      DATA_INPUT_COMPLEX(I) = FY(I)
      ENDDO
      
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)

      DO I =1,NINC
      FY_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
      ENDDO
      !  ============= FY =============
      
      !  ============= FZ =============
      DO I =1,NINC
      DATA_INPUT_COMPLEX(I) = FZ(I)
      ENDDO
      
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)

      DO I =1,NINC
      FZ_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
      ENDDO
      !  ============= FZ =============
      
      !  ============= AMX =============
      DO I =1,NINC
      DATA_INPUT_COMPLEX(I) = AMX(I)
      ENDDO
      
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)

      DO I =1,NINC
      AMX_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
      ENDDO
      !  ============= AMX =============
      
      !  ============= AMY =============
      DO I =1,NINC
      DATA_INPUT_COMPLEX(I) = AMY(I)
      ENDDO
      
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)

      DO I =1,NINC
      AMY_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
      ENDDO
      !  ============= AMY =============
      
      
      !  ============= AMZ =============
      DO I =1,NINC
      DATA_INPUT_COMPLEX(I) = AMZ(I)
      ENDDO
      
      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
      Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
      Status = DftiFreeDescriptor(My_Desc2_Handle)

      DO I =1,NINC
      AMZ_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
      ENDDO
      !  ============= AMZ =============
      
      
      ! ---------------------------------------------------
      ! PRINTING DATA
      ! J = NUMBER OF ELEMENT (NELE)
      ! K = NUMBER OF POINT (5)
      DO I = 1 ,NINC
      WRITE(5201,100) J,K,I,FX_OUT(I),FY_OUT(I),FZ_OUT(I),AMX_OUT(I),AMY_OUT(I),AMZ_OUT(I)
100   FORMAT (I5,2X,I5,2X,I5,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7)
      ENDDO
      ! ---------------------------------------------------
      
      ENDDO
      ENDDO
      ! END OF TRANFORM 
      
      
      ! TRANFORM DATA TO GID OUTPUT

      DO I = 1,NINC 
      ! PRINTING HEAD
               CALL PRNLCHD(NAME,NAML,1,0,1)
                WRITE (5202,150) NAME(1:NAML),I
                WRITE (5202,151)
                WRITE (5202,152)
150             FORMAT ('Result "XFrame-Local-Force (FREQ>TIME)"',1X,'"',A,'"',I5,1X,'Matrix  OnGaussPoints  "XFrame"')
151             FORMAT ('ComponentNames "Axial-Force" "Shear-S" "Shear-T" "Torsion" "Moment-S" "Moment-T"')
152             FORMAT ("Values") 
         DO J = 1,NELE
             INDEX = 1
             REWIND (5201)
             READ (5201,*)  ! READ HEAD
             FX_RESULT   = 0.
             FY_RESULT   = 0.
             FZ_RESULT   = 0.
             AMX_RESULT  = 0.
             AMY_RESULT  = 0.
             AMZ_RESULT  = 0.
             DO 1000 K =1,5*NINC*NELE
                READ (5201,*)  NUMBER_ELE,NPOINT,NTIME_STEP,ATEM_FX,ATEM_FY,ATEM_FZ,ATEM_AMX,ATEM_AMY,ATEM_AMZ
                ! STORAGE
                IF (I.EQ.NTIME_STEP.AND.INDEX.EQ.NPOINT.AND.J.EQ.NUMBER_ELE) THEN
                    FX_RESULT(INDEX)  = ATEM_FX 
                    FY_RESULT(INDEX)  = ATEM_FY 
                    FZ_RESULT(INDEX)  = ATEM_FZ 
                    AMX_RESULT(INDEX) = ATEM_AMX
                    AMY_RESULT(INDEX) = ATEM_AMY
                    AMZ_RESULT(INDEX) = ATEM_AMZ
                    INDEX      = INDEX + 1
                    IF (INDEX.GT.5) EXIT
                  
                ENDIF
1000          CONTINUE
             
              DO IJ = 1,5
                IF (IJ.EQ.1) THEN               
                WRITE (5202,200) J,FX_RESULT(IJ),FY_RESULT(IJ),FZ_RESULT(IJ),AMX_RESULT(IJ),AMY_RESULT(IJ),AMZ_RESULT(IJ)
200             FORMAT (I5,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7)
                ELSEIF (IJ.NE.1)THEN
                WRITE (5202,201) FX_RESULT(IJ),FY_RESULT(IJ),FZ_RESULT(IJ),AMX_RESULT(IJ),AMY_RESULT(IJ),AMZ_RESULT(IJ)
201             FORMAT (7X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7)  
                ENDIF
              ENDDO
         ENDDO
               WRITE (5202,153)
153            FORMAT ("End Values") 
               WRITE (5202,154)
154            FORMAT ("")
      ENDDO
      ! FUTURE NEED TO CHANGE THE DATA INTO GIDOUT.
      
      ! INTIAL SETTING
      FORCEX = 0.
      FORCEY = 0.
      FORCEZ = 0.
      ! REACTION FORCE AND MOMENT
      DO KK = 1,NBS
      NBOUNDARY = NODEABUNDARY(KK)
      REWIND (5203)
      DO I  = 1,NINC
      READ (5203,*) NDUM
      IF (NBOUNDARY.EQ.NDUM) THEN
      READ (5203,*) FORCEX(I),FORCEY(I),FORCEZ(I)
      ELSEIF (NBOUNDARY.EQ.NDUM)THEN
      READ (5203,*)
      ENDIF
      ENDDO
      
         !  ============= REACTION FORCE =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = FORCEX(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         FORCE_X_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= REACTION FORCE =============
      
         !  ============= REACTION FORCE =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = FORCEY(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         FORCE_Y_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= REACTION FORCE =============
      
         !  ============= REACTION FORCE =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = FORCEZ(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         FORCE_Z_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= REACTION FORCE =============
      
      ! WRITE RESULT OF REACTION FORCE
      DO I =1,NINC
      WRITE (5205,500) NBOUNDARY,FORCE_X_OUT(I),FORCE_Y_OUT(I),FORCE_Y_OUT(I)
      ENDDO
      
      ENDDO
      
      AMOMENT_X  = 0.
      AMOMENT_Y  = 0.
      AMOMENT_Z  = 0.
      DO KK = 1,NBS
      NBOUNDARY = NODEABUNDARY(KK)
      REWIND (5204)
      DO I  = 1,NINC
      READ (5204,*) NDUM
      IF (NBOUNDARY.EQ.NDUM)THEN
      READ (5204,*) AMOMENT_X(I),AMOMENT_Y(I),AMOMENT_Z(I)
      ELSEIF (NBOUNDARY.NE.NDUM)THEN
      READ (5204,*) 
      ENDIF
      ENDDO
      
         !  ============= REACTION MOMENT =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = AMOMENT_X(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         AMOMENT_X_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= REACTION MOMENT =============
      
         !  ============= REACTION MOMENT =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = AMOMENT_Y(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         AMOMENT_Y_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= REACTION MOMENT =============
      
         !  ============= REACTION MOMENT =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = AMOMENT_Z(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         AMOMENT_Z_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
      ENDDO
         !  ============= REACTION MOMENT =============
      
      ! WRITE RESULT OF REACTION FORCE
      DO I =1,NINC
      WRITE (5206,500) NBOUNDARY,AMOMENT_X_OUT(I),AMOMENT_Y_OUT(I),AMOMENT_Z_OUT(I)
      ENDDO
      
      ENDDO
500   FORMAT (I5,2X,F20.7,2X,F20.7,2X,F20.7)
      
      ! REACTION FORCE AND MOMENT TRANFORM SUBROUTINE 
      CALL REACTION_TRANFORM
      ! DISPLACEMENT TRANFORM SUBROUTINE
      !CALL DISPLACEMENT_TRANFORM
      
      END
      ! ======================================================================================
      SUBROUTINE REACTION_TRANFORM
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*200 NAME
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /BUNDARY_CONDITION/ INDEXBUNDARY(9999,6),NODEABUNDARY(9999)
      DIMENSION FX(NSN),FY(NSN),FZ(NSN)
      
      ! REACTION FORCE
      DO K =1,NINC
      REWIND (5205)
      CALL PRNLCHD(NAME,NAML,1,0,1)
      WRITE (5207,100) NAME(1:NAML),K
      WRITE (5207,101) 
      WRITE (5207,102) 
      FX  = 0.
      FY  = 0.
      FZ  = 0.
        DO I =1,NBS
          NBOUNDARY_READ = NODEABUNDARY(I)
             DO  J  = 1,NINC 
             READ (5205,*) NBOUNDARY,READ_FX,READ_FY,READ_FZ
             IF (NBOUNDARY.EQ.NBOUNDARY_READ.AND.K.EQ.J)THEN
             FX(NBOUNDARY) = READ_FX
             FY(NBOUNDARY) = READ_FY
             FZ(NBOUNDARY) = READ_FZ   
             ENDIF
             ENDDO
        ENDDO
        DO IJK = 1,NSN
        WRITE (5207,150) IJK,FX(IJK),FY(IJK),FZ(IJK)
        ENDDO
        WRITE (5207,103)
        WRITE (5207,104)
      ENDDO
      ! FUTURE NEED TO WRITTING THE DATA INTO GIDOUT FILE.
      ! DELETE THE FILE NUMBER 5203,5205,5207
      ! CLOSE(UNIT=5203,STATUS='DELETE')
      ! CLOSE(UNIT=5205,STATUS='DELETE')
      ! CLOSE(UNIT=5207,STATUS='DELETE')
      
      
      ! FOR REACTION MOMENT 
      DO K =1,NINC
      REWIND (5206)
      CALL PRNLCHD(NAME,NAML,1,0,1)
      WRITE (5208,120) NAME(1:NAML),K
      WRITE (5208,101) 
      WRITE (5208,102) 
      FX  = 0.
      FY  = 0.
      FZ  = 0.
        DO I =1,NBS
          NBOUNDARY_READ = NODEABUNDARY(I)
             DO  J  = 1,NINC 
             READ (5206,*) NBOUNDARY,READ_FX,READ_FY,READ_FZ
             IF (NBOUNDARY.EQ.NBOUNDARY_READ.AND.K.EQ.J)THEN
             FX(NBOUNDARY) = READ_FX
             FY(NBOUNDARY) = READ_FY
             FZ(NBOUNDARY) = READ_FZ   
             ENDIF
             ENDDO
        ENDDO
        DO IJK = 1,NSN
        WRITE (5208,150) IJK,FX(IJK),FY(IJK),FZ(IJK)
        ENDDO
        WRITE (5208,103)
        WRITE (5208,104)
      ENDDO
      
      ! FUTURE NEED TO WRITTING THE DATA INTO GIDOUT FILE.
      ! DELETE THE FILE NUMBER 5204,5206,5208
      ! CLOSE(UNIT=5204,STATUS='DELETE')
      ! CLOSE(UNIT=5206,STATUS='DELETE')
      ! CLOSE(UNIT=5208,STATUS='DELETE')
      
100   FORMAT ('Result "Global Reaction-Force (FREQ>TIME)"',1X,'"',A,'"',I5,1X,'Vector  OnNodes')
101   FORMAT ('ComponentNames  "F-X" "F-Y" "F-Z"')
102   FORMAT ('Values')
103   FORMAT ('End Values')
104   FORMAT ('')
      
120   FORMAT ('Result "Global Reaction-Moment (FREQ>TIME)"',1X,'"',A,'"',I5,1X,'Vector  OnNodes')
          
150   FORMAT (I5,2X,E15.7,2X,E15.7,2X,E15.7)
      END
      ! ======================================================================================
      SUBROUTINE TRANFORM_FATIGUE (NELE,INDEXSCAT,NWAVESCAT)
      USE MKL_DFTI
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*200 NAME
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /TIME/ DDT,CTIM,NINC
      DIMENSION FX_OUT(NINC),FY_OUT(NINC),FZ_OUT(NINC),AMX_OUT(NINC),AMY_OUT(NINC),AMZ_OUT(NINC)
      DIMENSION FX_OUT_MAX(NINC,5),FY_OUT_MAX(NINC,5),FZ_OUT_MAX(NINC,5),AMX_OUT_MAX(NINC,5),AMY_OUT_MAX(NINC,5),AMZ_OUT_MAX(NINC,5)
      
      COMMON /LHEADER/ LHDSTD(500),LHDCOM(500)
      
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND,FREQ_TIME_ANA
      COMMON /FREQ_FATIGUE/ NINC_FATIGUE
      
      DIMENSION FX(NINC),FY(NINC),FZ(NINC),AMX(NINC),AMY(NINC),AMZ(NINC)
      DIMENSION FX_RESULT(5),FY_RESULT(5),FZ_RESULT(5) 
      DIMENSION AMX_RESULT(5),AMY_RESULT(5),AMZ_RESULT(5)
      
      
      REAL :: DATA_INPUT_REAL(NINC)
      
      REAL :: DATA_OUTPUT_REAL(NINC)
      COMPLEX :: DATA_INPUT_COMPLEX(NINC)
      COMPLEX :: DATA_OUTPUT_COMPLEX(NINC)
      type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
      Integer :: Status
      
      ! TRANFORM DATA FOR INPUT IN FFT FUNCTION (INTEL MKL)
      WRITE(5201,90)
90    FORMAT ("NO.ELE     POINT       TIME_STEP     FX      FY      FZ      MX      MY      MZ")
      
      WRITE(*,2000)
2000  FORMAT ("START FFT")
      FX   = 0.
      FY   = 0.
      FZ   = 0.
      AMX  = 0.
      AMY  = 0.
      AMZ  = 0.
      
      !DO M = 1, NWAVESCAT
      
      DO J = 1,NELE
      DO K = 1,5
      INDEX = 1
      REWIND (5200)
        DO I = 1,NELE*NINC*5
           READ (5200,*) INCC,IFRAMESTORE,IG
           IF (IFRAMESTORE.EQ.J.AND.IG.EQ.K)THEN
           READ (5200,*) FX(INDEX),FY(INDEX),FZ(INDEX),AMX(INDEX),AMY(INDEX),AMZ(INDEX)
           INDEX = INDEX +1
           ELSEIF (IFRAMESTORE.NE.J.OR.IG.NE.K)THEN
           READ (5200,*)
           ENDIF
        ENDDO
      
         !  FFT FREQ. > TIME.
         !  ============= FX =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = FX(I)
         ENDDO
      
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)

         DO I =1,NINC
         FX_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= FX =============
         
         !  ============= FY =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = FY(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         FY_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= FY =============
         
         !  ============= FZ =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = FZ(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         FZ_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= FZ =============
         
         !  ============= AMX =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = AMX(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         AMX_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= AMX =============
         
         !  ============= AMY =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = AMY(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         AMY_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= AMY =============
         
         
         !  ============= AMZ =============
         DO I =1,NINC
         DATA_INPUT_COMPLEX(I) = AMZ(I)
         ENDDO
         
         Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_SINGLE,DFTI_COMPLEX, 1, INCC)
         Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
         Status = DftiCommitDescriptor(My_Desc2_Handle)
         Status = DftiComputeBackward(My_Desc2_Handle, DATA_INPUT_COMPLEX, DATA_OUTPUT_COMPLEX) !DftiComputeBackward
         Status = DftiFreeDescriptor(My_Desc2_Handle)
         
         DO I =1,NINC
         AMZ_OUT(I) = DATA_OUTPUT_COMPLEX(I)/NINC
         ENDDO
         !  ============= AMZ =============
      
         ! ---------------------------------------------------
         ! PRINTING DATA
         ! J = NUMBER OF ELEMENT (NELE)
         ! K = NUMBER OF POINT (5)
         DO I = 1 ,NINC
         WRITE(5201,100) J,K,I,FX_OUT(I),FY_OUT(I),FZ_OUT(I),AMX_OUT(I),AMY_OUT(I),AMZ_OUT(I)
         FX_OUT_MAX(I,K)   =  FX_OUT(I)
         FY_OUT_MAX(I,K)   =  FY_OUT(I)
         FZ_OUT_MAX(I,K)   =  FZ_OUT(I)
         AMX_OUT_MAX(I,K)  =  AMX_OUT(I)
         AMY_OUT_MAX(I,K)  =  AMY_OUT(I)
         AMZ_OUT_MAX(I,K)  =  AMZ_OUT(I)
100      FORMAT (I5,2X,I5,2X,I5,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7)
         
         ENDDO
      ! ---------------------------------------------------
      ENDDO
      
      CALL CALLSCATDATA (INDEXSCAT,WVHIGHT,PERIOD,MENITUDE,1,SEABED,WDEPTH,H1POS,H2POS,IRWAVE,GRAV,RHOW
     1                  ,RHOA,WKF,WFC,FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                   PER1,PER2,SHAPE1,SHAPE2)
      
      WRITE (5101,500) J,INDEXSCAT
500   FORMAT ("NUMBER OF ELEMENT",I5,2X,"WAVE SCAT",I5)
501   FORMAT (F10.5,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7)
      STEP          = PERIOD/NINC
      ATIME         = 0.
      NINC_FATIGUE  = FREQ_TIME_ANA/PERIOD
      DO II = 1,NINC_FATIGUE
      DO I = 1,NINC
      IF (MAXVAL(ABS(FX_OUT_MAX(I,1:5))).GT.MINVAL(ABS(FX_OUT_MAX(I,1:5)))) FX_OUT_FATIGUE = MAXVAL(FX_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(FX_OUT_MAX(I,1:5))).LE.MINVAL(ABS(FX_OUT_MAX(I,1:5)))) FX_OUT_FATIGUE = MINVAL(FX_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(FY_OUT_MAX(I,1:5))).GT.MINVAL(ABS(FY_OUT_MAX(I,1:5)))) FY_OUT_FATIGUE = MAXVAL(FY_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(FY_OUT_MAX(I,1:5))).LE.MINVAL(ABS(FY_OUT_MAX(I,1:5)))) FY_OUT_FATIGUE = MINVAL(FY_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(FZ_OUT_MAX(I,1:5))).GT.MINVAL(ABS(FZ_OUT_MAX(I,1:5)))) FZ_OUT_FATIGUE = MAXVAL(FZ_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(FZ_OUT_MAX(I,1:5))).LE.MINVAL(ABS(FZ_OUT_MAX(I,1:5)))) FZ_OUT_FATIGUE = MINVAL(FZ_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(AMX_OUT_MAX(I,1:5))).GT.MINVAL(ABS(AMX_OUT_MAX(I,1:5)))) AMX_OUT_FATIGUE = MAXVAL(AMX_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(AMX_OUT_MAX(I,1:5))).LE.MINVAL(ABS(AMX_OUT_MAX(I,1:5)))) AMX_OUT_FATIGUE = MINVAL(AMX_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(AMY_OUT_MAX(I,1:5))).GT.MINVAL(ABS(AMY_OUT_MAX(I,1:5)))) AMY_OUT_FATIGUE = MAXVAL(AMY_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(AMY_OUT_MAX(I,1:5))).LE.MINVAL(ABS(AMY_OUT_MAX(I,1:5)))) AMY_OUT_FATIGUE = MINVAL(AMY_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(AMZ_OUT_MAX(I,1:5))).GT.MINVAL(ABS(AMZ_OUT_MAX(I,1:5)))) AMZ_OUT_FATIGUE = MAXVAL(AMZ_OUT_MAX(I,1:5))
      IF (MAXVAL(ABS(AMZ_OUT_MAX(I,1:5))).LE.MINVAL(ABS(AMZ_OUT_MAX(I,1:5)))) AMZ_OUT_FATIGUE = MINVAL(AMZ_OUT_MAX(I,1:5))
      ATIME = ATIME + STEP
      WRITE(5101,501) ATIME,FX_OUT_FATIGUE,FY_OUT_FATIGUE,FZ_OUT_FATIGUE,AMX_OUT_FATIGUE,AMY_OUT_FATIGUE,AMZ_OUT_FATIGUE
      ENDDO
      ENDDO
      
      ENDDO
      ! END OF TRANFORM 
      
      
      ! TRANFORM DATA TO GID OUTPUT

      DO I = 1,NINC 
      ! PRINTING HEAD
               CALL PRNLCHD(NAME,NAML,1,0,1)
                WRITE (5202,150) NAME(1:NAML),I
                WRITE (5202,151)
                WRITE (5202,152)
150             FORMAT ('Result "XFrame-Local-Force (FREQ>TIME)"',1X,'"',A,'"',I5,1X,'Matrix  OnGaussPoints  "XFrame"')
151             FORMAT ('ComponentNames "Axial-Force" "Shear-S" "Shear-T" "Torsion" "Moment-S" "Moment-T"')
152             FORMAT ("Values") 
         DO J = 1,NELE
             INDEX = 1
             REWIND (5201)
             READ (5201,*)  ! READ HEAD
             FX_RESULT   = 0.
             FY_RESULT   = 0.
             FZ_RESULT   = 0.
             AMX_RESULT  = 0.
             AMY_RESULT  = 0.
             AMZ_RESULT  = 0.
             DO 1000 K =1,5*NINC*NELE
                READ (5201,*)  NUMBER_ELE,NPOINT,NTIME_STEP,ATEM_FX,ATEM_FY,ATEM_FZ,ATEM_AMX,ATEM_AMY,ATEM_AMZ
                ! STORAGE
                IF (I.EQ.NTIME_STEP.AND.INDEX.EQ.NPOINT.AND.J.EQ.NUMBER_ELE) THEN
                    FX_RESULT(INDEX)  = ATEM_FX 
                    FY_RESULT(INDEX)  = ATEM_FY 
                    FZ_RESULT(INDEX)  = ATEM_FZ 
                    AMX_RESULT(INDEX) = ATEM_AMX
                    AMY_RESULT(INDEX) = ATEM_AMY
                    AMZ_RESULT(INDEX) = ATEM_AMZ
                    INDEX      = INDEX + 1
                    IF (INDEX.GT.5) EXIT
                  
                ENDIF
1000          CONTINUE
             
              DO IJ = 1,5
                IF (IJ.EQ.1) THEN               
                WRITE (5202,200) J,FX_RESULT(IJ),FY_RESULT(IJ),FZ_RESULT(IJ),AMX_RESULT(IJ),AMY_RESULT(IJ),AMZ_RESULT(IJ)
200             FORMAT (I5,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7)
                ELSEIF (IJ.NE.1)THEN
                WRITE (5202,201) FX_RESULT(IJ),FY_RESULT(IJ),FZ_RESULT(IJ),AMX_RESULT(IJ),AMY_RESULT(IJ),AMZ_RESULT(IJ)
201             FORMAT (7X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7,2X,E15.7)  
                ENDIF
              ENDDO
         ENDDO
               WRITE (5202,153)
153            FORMAT ("End Values") 
               WRITE (5202,154)
154            FORMAT ("")
      ENDDO
      ! FUTURE NEED TO CHANGE THE DATA INTO GIDOUT.
      
      
      
      
      
      END
      
      