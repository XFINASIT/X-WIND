      SUBROUTINE READ_EXTERNALFORCE_TIMESERIES (ID)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*200 FILENUMBER,DIRECTORY,ROOTFST
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON/EXTERNALFORCE/ NEXTERNAL_GEN,NEXTERNAL_OFF
      DIMENSION ID(NSF,1)
      DIMENSION R(NEQ)
      DIMENSION NODE(1000),ATIME(100000),AFX(100000),AFY(100000),AFZ(100000),AMX(100000),AMY(100000),AMZ(100000)
      !ALLOCATABLE NODE(:),ATIME(:),AFX(:),AFY(:),AFZ(:),AMX(:),AMY(:),AMZ(:)
      
      !ALLOCATE (NODE(1000),ATIME(100000),AFX(100000),AFY(100000),AFZ(100000),AMX(100000),AMY(100000),AMZ(100000))
      
      ! INITIALIZE VALUE
      ROOTFST    = ".TXT"
      ACTIM      = 0.
      NODE       = 0.
      ATIME      = 0.
      AFX        = 0.
      AFY        = 0.
      AFZ        = 0.
      AMX        = 0.
      AMY        = 0.
      AMZ        = 0.

      OPEN (6002,STATUS='SCRATCH',ACCESS='DIRECT', FORM='UNFORMATTED',RECL=NEQ*2)
      OPEN (6003,STATUS='SCRATCH',ACCESS='DIRECT', FORM='UNFORMATTED',RECL=NEQ*2)
      ! EXTERNAL FORCE IN TIME SERIES FOR OFFSHORE LOAD
      DO JJ = 1,NEXTERNAL_OFF
           IF (JJ.LT.1000) WRITE( FILENUMBER,'(I3)') JJ
           IF (JJ.LT.100) WRITE( FILENUMBER,'(I2)') JJ
           IF (JJ.LT.10) WRITE( FILENUMBER,'(I1)') JJ
           DIRECTORY =  "Time series force/Timeseries_Offshore_"//TRIM(FILENUMBER)//TRIM(ROOTFST) 
           OPEN(UNIT=6001,FILE=DIRECTORY  ,STATUS='UNKNOWN')
           READ (6001,*)
           READ (6001,*) NOPTION
           READ (6001,*)
           READ (6001,*) NNODE
           READ (6001,*)
           
           DO I = 1,NNODE
           READ (6001,*) NODE(I)
           ENDDO
           
           READ (6001,*)
           READ (6001,*) NTOTAL
           READ (6001,*)
           
           DO II = 1,NTOTAL
           READ (6001,*) ATIME(II),AFX(II),AFY(II),AFZ(II),AMX(II),AMY(II),AMZ(II)
           ENDDO
           ATIME(NTOTAL + 1) = ATIME(NTOTAL) + 0.0001D0
           
           DO III = 1,NINC
           R     = 0.
           ACTIM = ACTIM + DDT
           IF (III.EQ.200)THEN
           ATEST = 1
           ENDIF
           ! STARTING INTERPOLATION (LINEAR INTERPOLATION)
           IF (NOPTION.EQ.1) THEN 
             DO 100 IIII = 1,NTOTAL+1
                 N2     = 0.
                 N3     = 0.
                IF (ACTIM.GT.MAXVAL(ATIME(1:NTOTAL+1))) EXIT
             IF (ACTIM.EQ.ATIME(IIII))THEN
                 ATIME3 = ATIME(IIII)
                 N3     = IIII
                 GOTO 110
             ELSEIF (ACTIM.GT.ATIME(IIII))THEN
                 ATIME1 = ATIME(IIII)
                 N1     = IIII
             ELSEIF (ACTIM.LT.ATIME(IIII))THEN
                 ATIME2 = ATIME(IIII)
                 N2     = IIII
                 EXIT
             ENDIF
             
100          CONTINUE  
             
           !IF (III.GT.NTOTAL) GOTO 120
           IF (ACTIM.GT.MAXVAL(ATIME(1:NTOTAL+1))) GOTO 120
           
           FX  = AFX(N1)+((ACTIM-ATIME1)*(AFX(N2)-AFX(N1))/(ATIME2-ATIME1))
           FY  = AFY(N1)+((ACTIM-ATIME1)*(AFY(N2)-AFY(N1))/(ATIME2-ATIME1))
           FZ  = AFZ(N1)+((ACTIM-ATIME1)*(AFZ(N2)-AFZ(N1))/(ATIME2-ATIME1))
           CMX = AMX(N1)+((ACTIM-ATIME1)*(AMX(N2)-AMX(N1))/(ATIME2-ATIME1))
           CMY = AMY(N1)+((ACTIM-ATIME1)*(AMY(N2)-AMY(N1))/(ATIME2-ATIME1))
           CMZ = AMZ(N1)+((ACTIM-ATIME1)*(AMZ(N2)-AMZ(N1))/(ATIME2-ATIME1))
           GOTO 120
           
110        FX   =  AFX(N3)
           FY   =  AFY(N3)
           FZ   =  AFZ(N3)
           CMX  =  AMX(N3)
           CMY  =  AMY(N3)
           CMZ  =  AMZ(N3)
           GOTO 120
           
120        IF (N2.EQ.0.AND.N3.EQ.0)THEN
           FX   =  0.
           FY   =  0.
           FZ   =  0.
           CMX  =  0.
           CMY  =  0.
           CMZ  =  0.
           ENDIF

           ! ADD INTO DOF
           DO J = 1,NNODE
           INO = NODE(J)
               DO K = 1,6
               ISF = IDOFCALL(IDOF,K)
               IEQ = ID(K,INO)
               IF (K.EQ.1) R(IEQ) = R(IEQ) + FX
               IF (K.EQ.2) R(IEQ) = R(IEQ) + FY
               IF (K.EQ.3) R(IEQ) = R(IEQ) + FZ
               IF (K.EQ.4) R(IEQ) = R(IEQ) + CMX
               IF (K.EQ.5) R(IEQ) = R(IEQ) + CMY
               IF (K.EQ.6) R(IEQ) = R(IEQ) + CMZ
               ENDDO
           ENDDO
           
           ELSEIF (NOPTION.EQ.0) THEN
               
           DO 200 IIII = 1,NTOTAL
               NEXIT = 0
               IF (ACTIM.EQ.ATIME(IIII)) THEN
               ATIME1 = IIII
               ATIME2 = IIII
               NEXIT = 1    
               EXIT
               ENDIF
200        CONTINUE  
           
           IF (NEXIT.EQ.1)THEN
               FX  = AFX(ATIME1)
               FY  = AFY(ATIME1)
               FZ  = AFZ(ATIME1)
               CMX = AMX(ATIME1)
               CMY = AMY(ATIME1)
               CMZ = AMZ(ATIME1)
           ELSEIF (NEXIT.EQ.0)THEN
               FX  = 0.
               FY  = 0.
               FZ  = 0.
               CMX = 0.
               CMY = 0.
               CMZ = 0.
           ENDIF
           
           ! ADD INTO DOF
           DO J = 1,NNODE
           INO = NODE(J)
               DO K = 1,6
               ISF = IDOFCALL(IDOF,K)
               IEQ = ID(K,INO)
               IF (K.EQ.1) R(IEQ) = R(IEQ) + FX
               IF (K.EQ.2) R(IEQ) = R(IEQ) + FY
               IF (K.EQ.3) R(IEQ) = R(IEQ) + FZ
               IF (K.EQ.4) R(IEQ) = R(IEQ) + CMX
               IF (K.EQ.5) R(IEQ) = R(IEQ) + CMY
               IF (K.EQ.6) R(IEQ) = R(IEQ) + CMZ
               ENDDO
           ENDDO
           
           
           ENDIF
           ! STORAGE FORCE VIA TIME 
           WRITE (6002,REC=III) R(1:NEQ)
           ENDDO
           
      ENDDO
      
      ! INITIALIZE INTERNAL STORAGE
      IF (NEXTERNAL_OFF.EQ.0)THEN
          DO II = 1,NINC
             R = 0.
             WRITE (6002,REC=II) R(1:NEQ) 
          ENDDO
      ENDIF
      CLOSE (6001)
      
      
      ! EXTERNAL FORCE IN TIME SERIES FOR GENERAL LOAD
      ACTIM = 0.
      DO JJ = 1,NEXTERNAL_GEN
           IF (JJ.LT.1000) WRITE(FILENUMBER,'(I3)') JJ
           IF (JJ.LT.100)  WRITE(FILENUMBER,'(I2)') JJ
           IF (JJ.LT.10)   WRITE(FILENUMBER,'(I1)') JJ
           DIRECTORY =  "Time series force/Timeseries_General_"//TRIM(FILENUMBER)//TRIM(ROOTFST) 
           OPEN(UNIT=6001,FILE=DIRECTORY  ,STATUS='UNKNOWN')
           READ (6001,*)
           READ (6001,*) NOPTION
           READ (6001,*)
           READ (6001,*) NNODE
           READ (6001,*)
           
           DO I = 1,NNODE
           READ (6001,*) NODE(I)
           ENDDO
           
           READ (6001,*)
           READ (6001,*) NTOTAL
           READ (6001,*)
           
           DO II = 1,NTOTAL
           READ (6001,*) ATIME(II),AFX(II),AFY(II),AFZ(II),AMX(II),AMY(II),AMZ(II)
           ENDDO
           
           DO III = 1,NINC
           R     = 0.
           ACTIM = ACTIM + DDT
           ! STARTING INTERPOLATION (LINEAR INTERPOLATION)
           IF (NOPTION.EQ.1) THEN 
             DO 300 IIII = 1,NTOTAL+1
                 N2     = 0.
                 N3     = 0.
                 IF (ACTIM.GT.MAXVAL(ATIME(1:NTOTAL+1))) EXIT
             IF (ACTIM.EQ.ATIME(IIII))THEN
                 ATIME3 = ATIME(IIII)
                 N3     = IIII
                 GOTO 310
             ELSEIF (ACTIM.GT.ATIME(IIII))THEN
                 ATIME1 = ATIME(IIII)
                 N1     = IIII
             ELSEIF (ACTIM.LT.ATIME(IIII))THEN
                 ATIME2 = ATIME(IIII)
                 N2     = IIII
                 EXIT
             ENDIF
300          CONTINUE   
             
           IF (III.GT.NTOTAL) GOTO 320
           IF (ACTIM.GT.MAXVAL(ATIME(1:NTOTAL+1))) GOTO 320
           
           FX  = AFX(N1)+((ACTIM-ATIME1)*(AFX(N2)-AFX(N1))/(ATIME2-ATIME1))
           FY  = AFY(N1)+((ACTIM-ATIME1)*(AFY(N2)-AFY(N1))/(ATIME2-ATIME1))
           FZ  = AFZ(N1)+((ACTIM-ATIME1)*(AFZ(N2)-AFZ(N1))/(ATIME2-ATIME1))
           CMX = AMX(N1)+((ACTIM-ATIME1)*(AMX(N2)-AMX(N1))/(ATIME2-ATIME1))
           CMY = AMY(N1)+((ACTIM-ATIME1)*(AMY(N2)-AMY(N1))/(ATIME2-ATIME1))
           CMZ = AMZ(N1)+((ACTIM-ATIME1)*(AMZ(N2)-AMZ(N1))/(ATIME2-ATIME1))
           GOTO 320
           
310        FX   =  AFX(N3)
           FY   =  AFY(N3)
           FZ   =  AFZ(N3)
           CMX  =  AMX(N3)
           CMY  =  AMY(N3)
           CMZ  =  AMZ(N3)
           GOTO 320
           
320        IF (N2.EQ.0.AND.N3.EQ.0)THEN
           FX   =  0.
           FY   =  0.
           FZ   =  0.
           CMX  =  0.
           CMY  =  0.
           CMZ  =  0.
           ENDIF
           
           ! ADD INTO DOF
           DO J = 1,NNODE
           INO = NODE(J)
               DO K = 1,6
               ISF = IDOFCALL(IDOF,K)
               IEQ = ID(K,INO)
               IF (K.EQ.1) R(IEQ) = R(IEQ) + FX
               IF (K.EQ.2) R(IEQ) = R(IEQ) + FY
               IF (K.EQ.3) R(IEQ) = R(IEQ) + FZ
               IF (K.EQ.4) R(IEQ) = R(IEQ) + CMX
               IF (K.EQ.5) R(IEQ) = R(IEQ) + CMY
               IF (K.EQ.6) R(IEQ) = R(IEQ) + CMZ
               ENDDO
           ENDDO
           
           ELSEIF (NOPTION.EQ.0) THEN
               
           DO 400 IIII = 1,NTOTAL
               NEXIT = 0
               IF (ACTIM.EQ.ATIME(IIII)) THEN
               ATIME1 = IIII
               ATIME2 = IIII
               NEXIT = 1    
               EXIT
               ENDIF
400        CONTINUE  
           
           IF (NEXIT.EQ.1)THEN
               FX  = AFX(ATIME1)
               FY  = AFY(ATIME1)
               FZ  = AFZ(ATIME1)
               CMX = AMX(ATIME1)
               CMY = AMY(ATIME1)
               CMZ = AMZ(ATIME1)
           ELSEIF (NEXIT.EQ.0)THEN
               FX  = 0.
               FY  = 0.
               FZ  = 0.
               CMX = 0.
               CMY = 0.
               CMZ = 0.
           ENDIF
           
           ! ADD INTO DOF
           DO J = 1,NNODE
           INO = NODE(J)
               DO K = 1,6
               ISF = IDOFCALL(IDOF,K)
               IEQ = ID(K,INO)
               IF (K.EQ.1) R(IEQ) = R(IEQ) + FX
               IF (K.EQ.2) R(IEQ) = R(IEQ) + FY
               IF (K.EQ.3) R(IEQ) = R(IEQ) + FZ
               IF (K.EQ.4) R(IEQ) = R(IEQ) + CMX
               IF (K.EQ.5) R(IEQ) = R(IEQ) + CMY
               IF (K.EQ.6) R(IEQ) = R(IEQ) + CMZ
               ENDDO
           ENDDO
           ENDIF
           ! STORAGE FORCE VIA TIME 
           WRITE (6003,REC=III) R(1:NEQ)
           ENDDO
      ENDDO
      
      ! INITIALIZE INTERNAL STORAGE
      IF (NEXTERNAL_GEN.EQ.0)THEN
          DO II = 1,NINC
             R = 0.
             WRITE (6003,REC=II) R(1:NEQ) 
          ENDDO
      ENDIF
      
      !DEALLOCATE (NODE,ATIME,AFX,AFY,AFZ,AMX,AMY,AMZ)
      CLOSE (6001)
      END
      
C	--------------------------------------------------------
      SUBROUTINE EXTERNAL_FORCE (INC,R)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      DIMENSION R(NEQ),RGEN(NEQ),ROFF(NEQ)
      
      READ (6003,REC=INC,IOSTAT=ios) RGEN(1:NEQ)
      IF (IOSTAT.NE.0) RGEN(1:NEQ) = 0.
      READ (6002,REC=INC,IOSTAT=ios) ROFF(1:NEQ)
      IF (IOSTAT.NE.0) ROFF(1:NEQ) = 0.
      
      R(1:NEQ) = RGEN(1:NEQ) + ROFF(1:NEQ)
      
      END