      SUBROUTINE ADDSTIFFNESS (DD,NWKM,MAXAA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
       COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
       COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
       COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /OPT_ADDED_STIFFNESS_MATRIX/ NMSTIFF
      DIMENSION DD(1),MAXAA(1)
      COMMON A(9000000),IA(9000000)

      IF (NMSTIFF.NE.0D0)THEN
      CALL ADDED_STIFFESS (DD,NWKM,NSF,IA(LID),MAXAA)
      ENDIF

      END
C     ===============================================================================
      SUBROUTINE ADDED_STIFFESS (DD,NWKM,NSF,ID,MAXAA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
       COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
       COMMON /NUMB2/ IDOFB(9)
      COMMON /ADD_STIFFNESS/ NTOTOAL(100),NODE_STIFF(100,100),STIFFNESSV(81,100)
      COMMON /OPT_ADDED_STIFFNESS_MATRIX/ NMSTIFF
      DIMENSION DD(1),MAXAA(1)
      DIMENSION ID(NSF,1)
      
      NDOF = 0
      DO I = 1,9
          IF (IDOFB(I).EQ.0) NDOF = NDOF + 1
      ENDDO
      
      DO I = 1,NMSTIFF
      NGROUP = NTOTOAL(I) !NODE_STIFF(I,1)
          DO J = 1,NGROUP
          ! ADDING NODE POSITION
          NNODE = NODE_STIFF(J,I) 
               INDEXSTART = 0.0D0
               INDEXEND   = 0.0D0
               DO K = 1,NDOF
               IEQ   = ID(K,NNODE)
               INDEXSTART = INDEXEND + (I-1)*(NDOF*NDOF) + 1D0
               INDEXEND   = INDEXSTART + NDOF - 1
               DD(MAXAA(IEQ):MAXAA(IEQ)+NDOF-1) = DD(MAXAA(IEQ):MAXAA(IEQ)+NDOF-1) + STIFFNESSV(INDEXSTART:INDEXEND,I)
               ENDDO
          ENDDO
      ENDDO
      
      
      END
C     ===============================================================================
      SUBROUTINE READADDEDSTIFFNESS
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /NUMB2/ IDOFB(9)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /OPT_ADDED_STIFFNESS_MATRIX/ NMSTIFF
      COMMON /ADD_STIFFNESS/ NTOTOAL(100),NODE_STIFF(100,100),STIFFNESSV(81,100)
      DIMENSION STIFF_VALUE(9)
      
      
      NDOF = 0
      DO I = 1,9
          IF (IDOFB(I).EQ.0) NDOF = NDOF + 1
      ENDDO
      
      DO I = 1,NMSTIFF
      READ (ITI,*) NTOTOAL(I),NODE_STIFF(1:NTOTOAL(I),I)
          DO J = 1,NDOF
          READ (ITI,*) STIFF_VALUE(1:NDOF)
               DO K =1,NDOF
               STIFFNESSV(INDEX,I) = STIFF_VALUE(K)  
               ENDDO
          ENDDO
      ENDDO
      
      ! WARNING
      IF (NMSTIFF.GT.100D0) WRITE (*,100) 
100   FORMAT ("********* STORAGE FULL OF ADDED STIFFNESS  **********")      
      END
      