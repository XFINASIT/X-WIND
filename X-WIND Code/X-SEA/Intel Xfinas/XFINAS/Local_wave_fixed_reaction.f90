!================================================================================================
!================================================================================================
      
      SUBROUTINE KEEP_WAVE_FORCE_REACTION(FRMFOC,MLE,ILCAS)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N) 
      
      DIMENSION FRMFOC(7,2)
      DIMENSION FRMFOC_WAVE(14)
      
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
      !COMMON / WAVEREACFIXPARAMETER /  WAVEREACFIX(9999,7,2)
      
      INDEX = 1
      
      DO J = 1,7
	DO K = 1,2  
      !WAVEREACFIX(MLE,J,K) = FRMFOC(J,K)
      FRMFOC_WAVE(INDEX) = FRMFOC(J,K)
      INDEX = INDEX + 1
	ENDDO
      ENDDO

!      WRITE(60,1) ILCAS , MLE ,  FRMFOC_WAVE(1:14)
!1     FORMAT(I5,2X,I5,14E20.8)  
      
      IF (ILCAS.EQ.1) THEN
      NEC_RE = MLE
      ELSEIF (ILCAS.GT.1D0)THEN
      NEC_RE = ILCAS*NELE + MLE
      ENDIF
      
      WRITE(810,REC=NEC_RE) ILCAS,MLE,FRMFOC_WAVE(1:14)   
      
      ! WRITE(64,1) ILCAS , MLE ,  FRMFOC_WAVE(1:14) ! TOEY 10/21
      
      RETURN
      END SUBROUTINE
      !================================================================================================
      !================================================================================================
      !================================================================================================
      SUBROUTINE UNVECTER_ELEMENT_OF_WAVE_FORCE(IEL,VR,VS,VT)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC  
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /STCAS/ ILC   
      DIMENSION VR(3)
      DIMENSION VS(3)
      DIMENSION VT(3)
      
!      TOEY 10/2021      
!      REWIND(61)
!     
!      DO WHILE(1.NE.0) 
!          
!          READ(61,1,END=100) IELE  , VR(1:3) , VS(4:6)  ,VT(7:9)   !KEEPLOCALVECTOR(1:9)
!1         FORMAT(I6,2X,10F20.7)            
!          
!          IF( IELE .EQ. IEL  ) EXIT
!          
!      END DO
      
      IF (ISOLOP.NE.0) LCS_IN = 1D0
      IF (ISOLOP.EQ.5) LCS_IN = 1D0
      IF (ISOLOP.EQ.1) THEN
      CALL OFFSHORE_LOAD_CASE_NUMBER (ILC,LCS_IN,OPT_OFF,"CALL")
      ENDIF
      
      IF (LCS_IN.EQ.1) THEN
      NEC_VEC = IEL
      ELSEIF (LCS_IN.GT.1D0)THEN
      NEC_VEC = LCS_IN*NELE + IEL
      ENDIF
      READ (811,REC=NEC_VEC) IELE , VR(1:3) , VS(4:6)  ,VT(7:9)
      
      RETURN
      
100   CONTINUE
      
      !WRITE(*,*) "  ERROR AT SUBROUTINE UNVECTER_ELEMENT_OF_WAVE_FORCE"
      !STOP
      
      RETURN
      END SUBROUTINE
      !================================================================================================
      SUBROUTINE WAVEFORCE_FIXED_REACTION(IEL,LCS_IN,WAVEFIXEDFORCEMEM) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      DIMENSION WAVEFIXEDFORCEMEM(7,2)
      DIMENSION WAVEFIXEDFORCEMEMREAD(14)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC  
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT     
      !REWIND(60)
            
      
!      DO WHILE(1.NE.0) 
!          READ(60,*,END=100)  LCSDATA , IELE , WAVEFIXEDFORCEMEMREAD(1:14)         
!          
!          IF( IELE .EQ. IEL .AND. LCSDATA .EQ. LCS ) EXIT
!          
!      END DO
      IF (ISOLOP.NE.0) LCS_INPUT = 1
      IF (ISOLOP.EQ.5) LCS_INPUT = 1
      
      IF (LCS_INPUT.EQ.1) THEN
      NEC_RE = IEL
      ELSEIF (LCS_INPUT.GT.1D0)THEN
      NEC_RE = LCS_INPUT*NELE + IEL
      ENDIF
      
      READ (810,REC=NEC_RE,IOSTAT=IReason) LCSDATA , IELE , WAVEFIXEDFORCEMEMREAD(1:14)     
      IF (IReason.NE.0.0D0) THEN
          WAVEFIXEDFORCEMEMREAD(1:14)  = 0.
          IELE = 0.
          LCSDATA = 0.
      ENDIF
      
      INDEX = 1
      DO J = 1,7
	DO K = 1,2
          
       WAVEFIXEDFORCEMEM(J,K) = WAVEFIXEDFORCEMEMREAD(INDEX)
       !FRMFOC_WAVE(INDEX) = FRMFOC(J,K)
       INDEX = INDEX + 1
	ENDDO
      ENDDO      

      RETURN
      
100   CONTINUE
           
      WAVEFIXEDFORCEMEM=0.0d0
      !WRITE(*,*) " ADDING DATA ELEMENT", IEL , "LOAD CASE", LCS
      !WRITE(*,*) "  ERROR AT SUBROUTINE  WAVEFORCE_FIXED_REACTION"
      !STOP
      
      RETURN
      END SUBROUTINE
      !================================================================================================
      
!================================================================================================
      
      SUBROUTINE TRAPIZOIDAL_TRANSFORMATION(VR,VS,VT,OPTION,WFLOCAL,LCS_IN)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPTION
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      !COMMON /TRAPIZOIDALLCS/ LCSANA(100), NLCSDATA(100)
      !COMMON /TRAPIZOIDALLCSS/ LCSTOTAL

C	==================================================================
C	TRAPIZOIDAL RULE   
      !COMMON / WAVEREACFIXPARAMETER /  WAVEREACFIX(9999,7,2)
C	==================================================================      
      !COMMON / ELEMENT_LOCALWAVEREAC /  LOCALELEMENTWAVEVECTOR(9999,9)
C     SAVE ELEMENT NUMBER      
      COMMON /NELEM/ IEL  
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC  
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /STCAS/ ILC
      
      DIMENSION TT(14,14),FIXLRWAVEFRAME(14),TTT(14,14)
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION WFLOCAL(14)
      
      DIMENSION WAVEFIXEDFORCEMEM(7,2)
      
      DIMENSION AKEEPLOCALVECTOR(9)
      
      WFLOCAL = 0.0D0
      
          
      IF(OPTION.EQ.'RED'.AND.ITASK.NE.5)THEN
          
         ! VR(1:3)   = LOCALELEMENTWAVEVECTOR(IEL,1:3)
         ! VS(1:3)   = LOCALELEMENTWAVEVECTOR(IEL,4:6)
         ! VT(1:3)   = LOCALELEMENTWAVEVECTOR(IEL,7:9)
          
      CALL UNVECTER_ELEMENT_OF_WAVE_FORCE(IEL,VR,VS,VT)

      CALL TT1A (VR,VS,VT,TT)
 
      CALL WAVEFORCE_FIXED_REACTION(IEL,LCS_IN,WAVEFIXEDFORCEMEM) 
      
      FIXLRWAVEFRAME(1)	=    WAVEFIXEDFORCEMEM(1,1)           !WAVEREACFIX(IEL,1,1)   
      FIXLRWAVEFRAME(2)	=    WAVEFIXEDFORCEMEM(2,1)           !WAVEREACFIX(IEL,2,1)   
      FIXLRWAVEFRAME(3)	=    WAVEFIXEDFORCEMEM(3,1)           !WAVEREACFIX(IEL,3,1)   
      FIXLRWAVEFRAME(4)	=    WAVEFIXEDFORCEMEM(4,1)           !WAVEREACFIX(IEL,4,1)   
      FIXLRWAVEFRAME(5)	=    WAVEFIXEDFORCEMEM(5,1)           !WAVEREACFIX(IEL,5,1)   
      FIXLRWAVEFRAME(6)	=    WAVEFIXEDFORCEMEM(6,1)           !WAVEREACFIX(IEL,6,1)   
      FIXLRWAVEFRAME(7)	=    WAVEFIXEDFORCEMEM(7,1)           !WAVEREACFIX(IEL,7,1)   
      FIXLRWAVEFRAME(8)	=    WAVEFIXEDFORCEMEM(1,2)           !WAVEREACFIX(IEL,1,2)   
      FIXLRWAVEFRAME(9)	=    WAVEFIXEDFORCEMEM(2,2)           !WAVEREACFIX(IEL,2,2)   
      FIXLRWAVEFRAME(10)	=    WAVEFIXEDFORCEMEM(3,2)           !WAVEREACFIX(IEL,3,2)   
      FIXLRWAVEFRAME(11)	=    WAVEFIXEDFORCEMEM(4,2)           !WAVEREACFIX(IEL,4,2)   
      FIXLRWAVEFRAME(12)	=    WAVEFIXEDFORCEMEM(5,2)           !WAVEREACFIX(IEL,5,2)   
      FIXLRWAVEFRAME(13)	=    WAVEFIXEDFORCEMEM(6,2)           !WAVEREACFIX(IEL,6,2)   
      FIXLRWAVEFRAME(14)	=    WAVEFIXEDFORCEMEM(7,2)           !WAVEREACFIX(IEL,7,2) 
      
!      WRITE(65,11) IEL,LCS_IN,FIXLRWAVEFRAME(1:14)
!11    FORMAT(I5,2X,I5,14E20.8)       

      !TTT = TRANSPOSE(TT)
      NEF = 14    
      DO IEF = 1,NEF
	!FIXLR(IEF) = 0.0
	!FIXLO(IEF) = 0.0
      WFLOCAL(IEF) = 0.0 
	DO JEF = 1,NEF
	WFLOCAL(IEF) = WFLOCAL(IEF) + TT(JEF,IEF)*FIXLRWAVEFRAME(JEF)  !VARY FIXEND
      !FIXLR(IEF) = FIXLR(IEF) + TT(IEF,JEF)*FIXEN(JEF)  !VARY FIXEND
	!FIXLO(IEF) = FIXLO(IEF) + TT(IEF,JEF)*FIXEO(JEF)  !CONT FIXEND
	ENDDO
      ENDDO

      RETURN
      ELSEIF(OPTION.EQ.'WRT')THEN
          
      !LOCALELEMENTWAVEVECTOR(IEL,1:3) = VR(1:3)
      !LOCALELEMENTWAVEVECTOR(IEL,4:6) = VS(1:3)    
      !LOCALELEMENTWAVEVECTOR(IEL,7:9) = VT(1:3)
      
      AKEEPLOCALVECTOR(1:3) = VR(1:3)
      AKEEPLOCALVECTOR(4:6) = VS(1:3) 
      AKEEPLOCALVECTOR(7:9) = VT(1:3)
      
!      WRITE(61,1) IEL , AKEEPLOCALVECTOR(1:9)   ! TOEY 10/2021
!1     FORMAT(I6,2X,10F20.7) 
      ! TOEY 10/2021
      IF (ISOLOP.NE.0) LCS_IN = 1D0
      IF (ISOLOP.EQ.5) LCS_IN = 1D0
      IF (ISOLOP.EQ.1) THEN
      CALL OFFSHORE_LOAD_CASE_NUMBER (ILC,LCS_IN,OPT_OFF,"CALL")
      ENDIF
      
      IF (LCS_IN.EQ.1) THEN
      NEC_VEC = IEL
      ELSEIF (LCS_IN.GT.1D0)THEN
      NEC_VEC = LCS_IN*NELE + IEL
      ENDIF
      WRITE (811,REC=NEC_VEC) IEL , AKEEPLOCALVECTOR(1:9)
      
      !write(63,1) IEL, AKEEPLOCALVECTOR(1:9)    ! TOEY 10/2021
      
      ENDIF
      
      RETURN
      END SUBROUTINE
      
!================================================================================================  