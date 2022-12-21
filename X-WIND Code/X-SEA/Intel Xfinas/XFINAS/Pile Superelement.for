      SUBROUTINE PILE_SUPERELEMENT (MAXA,ILC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2  
      !ALLOCATABLE AK_MST(:,:),AK_SLASS(:,:),AK_SLAMS(:,:),AK_SLASM(:,:)
      !ALLOCATABLE AK_SUM(:,:)
      !ALLOCATABLE FORCE(:),FORCE_MST(:),FORCE_SLA(:),DISP(:),DISP_MST(:),DISP_SLA(:)
      ! 
      !DIMENSION MST_MAT(100)
      ALLOCATABLE AKS(:,:) 
C     ===============================================================================
      POINTER :: AA(:),CC(:),MAXC_STIFFNESS(:),MAXC_MASS(:),IPERM(:),X(:)
      DIMENSION MAXA(1)
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /PARDIS/ LPTO(64) 
      COMMON A(7000000),IA(6000000)
      
      DIMENSION MST_MAT(100)
      
      ! ------- STIFFNESS FOR TOTAL STFFNESS FILE------
C     IPFAC DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT
      !TYP1 = "STFF"
      TYP2 = "TEMP"
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,0,'FIRS')  

      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0,'FIRS')
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0,'FIRS')
      
	CALL MCALFIL(KSREC,TYP2)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
      
      ALLOCATE(AA(NWK),MAXC_STIFFNESS(NWK)) 
	IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0,'FIRS')
      
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	READ(NFLCH) MAXC_STIFFNESS(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
      
          !REWIND(13000)
          !WRITE(13000,*) NEQ,NWK
          !WRITE(13000,*) MAXA(1:NEQ+1)
          !WRITE(13000,*) MAXC_STIFFNESS(1:NWK)
          !WRITE(13000,*) 'TOTAL STIFFNESS'
          !WRITE(13000,*) AA(1:NWK)         
      ! ------------------------
      
      !! OPEN FILE FOR PILE-SUPER ELEMENT BY BJ ---------        
      !OPEN(UNIT=13000,FILE='Substructuring/Total_Stiffness.dat'             ,STATUS='OLD')  
      !OPEN(UNIT=13001,FILE='Substructuring/FULL_Stiffness.dat'              ,STATUS='UNKNOWN')  
      !OPEN(UNIT=13002,FILE='Substructuring/Substructuring_data.dat'         ,STATUS='UNKNOWN')  
      !OPEN(UNIT=13003,FILE='Substructuring/Condensed_Stiffness.dat'         ,STATUS='UNKNOWN')  
      !OPEN(UNIT=13004,FILE='Substructuring/Condensed_Diagonal_Stiffness.dat'         ,STATUS='UNKNOWN') 
      
      !READ TOTAL STIFFNESS(NWK) AND EXPAND TO FULL STIFFNESS MATRIX (NEQ,NEQ)
      ALLOCATE(AKS(NEQ,NEQ)) ! FULL STIFFNESS MATRIX (NEQ,NEQ)  
      CALL READSTIFF(NEQ,NWK,MAXC_STIFFNESS,MAXA,AA,AKS) 
      DO I =1,NEQ
          DO J = 1,NEQ
      WRITE (300,100) I,J,AKS(I,J)
      ENDDO
      ENDDO
100   FORMAT (I5,2X,I5,2X,E20.10)
      !READ SUBSTRUCTURING DATA
      CALL READSUB(NMST,MST_MAT)      
      
      CALL STIFFPILE(AKS,NEQ,NWK,NMST,MST_MAT,ILC)

      
      DEALLOCATE(AA,MAXC_STIFFNESS)
      DEALLOCATE(AKS)
      
      END
C	=============================================================
C	=============================================================
C	=============================================================    
      SUBROUTINE READSTIFF(NEQ,NWK,MAXC,MAXA,AA,AKS_OUT) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)     
      DIMENSION MAXA(NEQ+1),MAXC(NWK),AA(NWK)
      DIMENSION AKS_OUT(NEQ,NEQ)
      !DIMENSION AKS(NEQ,NEQ) ! FULL STIFFNESS MATRIX (NEQ,NEQ)   
      ALLOCATABLE AKS(:,:) ! FULL STIFFNESS MATRIX (NEQ,NEQ)
      ALLOCATE(AKS(NEQ,NEQ)) ! FULL STIFFNESS MATRIX (NEQ,NEQ)  

      !EXPAND UPPER-TRIANGULAR STIFFNESS MATRIX(NWK) TO FULL STIFFNESS MATRIX (NEQ,NEQ)
          !INITIATE THE MATRIX
          AKS = 0.
          
          INEQ = 0      
          NPOS = 1
          DO IEQ = 1,NEQ
              IMA = MAXA(IEQ+1)
              DO INWK = 1,NWK
                  IMC = MAXC(NPOS)   
                  AKS(IEQ,IMC) = AA(NPOS)
                  IF (AA(NPOS).LT.0.000000001D0.AND.AA(NPOS).GT.-0.000000001D0) AKS(IEQ,IMC) = 0.0D0
                  NPOS = NPOS + 1 
                  IF(NPOS.EQ.IMA) EXIT                                   
              ENDDO    
          ENDDO
          
          DO II = 1,NEQ
              DO JJ = 1,NEQ
                 AKS(JJ,II) = AKS(II,JJ)
              ENDDO    
          ENDDO   
          
          AKS_OUT = AKS
          
      DEALLOCATE(AKS)
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================   
      SUBROUTINE READSUB(NMST_OUT,MST_MAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)  
      COMMON /SOIL_LOCATE/NOUT_NODE_SUPER(1000,500),INDEX_OUT_G(500),NGROUP
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /MGRAV/ NGRAV  
      DIMENSION MST_MAT(100)
      ALLOCATABLE AX(:),AY(:),AZ(:)
      
      !READ(13002,*) 
      !READ(13002,*) 
      !READ(13002,*) NMST, MST_MAT(1:NMST)
      ALLOCATE (AX(NSN),AY(NSN),AZ(NSN))
      AX = 0.
      AY = 0.
      AZ = 0.
      IF (NGRAV.EQ.1)THEN
      DO I = 1 ,INDEX_OUT_G(1)-1
      CALL RELFILL('@XYZ',AX(I),1,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY(I),2,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ(I),3,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
      ENDDO
      
      NCOOR_MAX = MINLOC(AX(1:NSN),DIM=1)
      NCOOR_MIN = MAXLOC(AX(1:NSN),DIM=1)
      
      IF (AX(NCOOR_MAX).GT.AX(NCOOR_MIN)) MST_MAT(1) = NCOOR_MAX
      IF (AX(NCOOR_MAX).LE.AX(NCOOR_MIN)) MST_MAT(1) = NCOOR_MIN
      
      ELSEIF (NGRAV.EQ.2)THEN
      DO I = 1 ,INDEX_OUT_G(1)-1
      CALL RELFILL('@XYZ',AX(I),1,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY(I),2,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ(I),3,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
      ENDDO
      
      NCOOR_MAX = MINLOC(AY(1:NSN),DIM=1)
      NCOOR_MIN = MAXLOC(AY(1:NSN),DIM=1)
      
      IF (AY(NCOOR_MAX).GT.AY(NCOOR_MIN)) MST_MAT(1) = NCOOR_MAX
      IF (AY(NCOOR_MAX).LE.AY(NCOOR_MIN)) MST_MAT(1) = NCOOR_MIN
      ELSEIF (NGRAV.EQ.3)THEN
      DO I = 1 ,INDEX_OUT_G(1)-1
      CALL RELFILL('@XYZ',AX(I),1,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY(I),2,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ(I),3,NOUT_NODE_SUPER(I,1),0)  !GETTING NODAL COORDINATE
      ENDDO
      
      NCOOR_MAX = MINLOC(AZ(1:NSN),DIM=1)
      NCOOR_MIN = MAXLOC(AZ(1:NSN),DIM=1)
      
      IF (AZ(NCOOR_MAX).GT.AZ(NCOOR_MIN)) MST_MAT(1) = NCOOR_MAX
      IF (AZ(NCOOR_MAX).LE.AZ(NCOOR_MIN)) MST_MAT(1) = NCOOR_MIN
          
      ENDIF
      NMST_OUT = 1D0
      
      DEALLOCATE (AX,AY,AZ)
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================      
      SUBROUTINE STIFFPILE(AKS,NEQ,NWK,NMST,MST_MAT,ILC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)   
      CHARACTER*200 LHDSTD,LHDCOM,NAMEL
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON/OFFSHORE_CASE/ LGEN,IOFFL,IORI_OFFSHORE
      ALLOCATABLE AK_COND(:,:)
      ALLOCATABLE AK_MST(:,:),AK_SLAMS(:,:),AK_SLASS(:,:),AK_SLASM(:,:)
      ALLOCATABLE AK_COND1(:,:),AK_COND2(:,:),AK_SLASS_INV(:,:) 
      ALLOCATABLE ISLA_MAT(:)
      ALLOCATABLE AK_COND_INV(:,:),FORCE_MST(:),FORCE_SLA(:),FORCE(:),DISP_MST(:)
      ALLOCATABLE NSTART(:),NEND(:)
      COMMON /LHEADER/ LHDSTD(500),LHDCOM(500)
      COMMON /BOUNDARY/ NBOUNDARY(9,500)
      COMMON /NUMB2/ IDOFB(9)
      
      !AK_COND : CONDENSED STIFFNESS OF ONE SUPER ELEMENT
      !AK_MST : MASTER STIFFNESS 
      !AK_SLAMS : MASTER-SLAVE STIFFNESS
      !AK_SLASS : SLAVE-SLAVE STIFFNESS
      !AK_SLASM : SLAVE-MASTER STIFFNESS
      !AK_COND1 : AK_SLASS * AK_SLASM
      !AK_COND2 : AK_SLAMS * AK_COND1
      !NMST : THE NUMBER OF MASTER NODES
      !MST_MAT : LIST OF MASTER NODES
      
      DIMENSION AKS(NEQ,NEQ),MST_MAT(100)
      DIMENSION AK_COND_DIAG(6)
      
      NDOF_STRUCTURE = 0
      DO I = 1,9
          IF (IDOFB(I).EQ.0D0) NDOF_STRUCTURE = NDOF_STRUCTURE + 1
      ENDDO
      
      !REWIND(13001)
      !READ(13001,*) 
      !READ(13001,*) AKS(1:NEQ,1:NEQ)
      
      !DEFINE SIZES OF MASTER MATRIX(NEQ_MST) AND SLAVE MATRIX(NEQ_SLA)
      NEQ_MST = NMST*NDOF_STRUCTURE !6 = THE NUMBER OF D.O.F OF MASTER NODES
      NEQ_SLA = NEQ-NEQ_MST
      NUM_TOT = NEQ/NDOF_STRUCTURE ! TOTAL NUMBER OF NODES
      !NMST =  THE NUMBER OF MATER NODES
      NSLA = NEQ_SLA/NDOF_STRUCTURE ! THE NUMBER OF SALVE NODES
      
      !DEFINE THE LIST OF SLAVE NODES
      ALLOCATE(ISLA_MAT(NSLA))
      NUM = 0 
      NUM_SLA_MAT = 0
      DO INUM = 1,NUM_TOT
          NUM = NUM + 1
          NUM_SLA = 0
          DO IMST = 1,NMST
              NODE_MST = MST_MAT(IMST)
              IF(NUM.NE.NODE_MST) THEN
                  NUM_SLA = NUM_SLA + 1                  
              ENDIF    
              IF(NUM_SLA.EQ.NMST) THEN
                  NUM_SLA_MAT = NUM_SLA_MAT + 1
                  ISLA_MAT(NUM_SLA_MAT) = NUM
              ENDIF
          ENDDO                 
      ENDDO
      
      
      ! REMOVE SLAVE NODE FROM BOUNDARY CONDITION
      
      ALLOCATE(AK_COND(NEQ_MST,NEQ_MST))          
      ALLOCATE(AK_MST(NEQ_MST,NEQ_MST),AK_SLAMS(NEQ_MST,NEQ_SLA),AK_SLASS(NEQ_SLA,NEQ_SLA),AK_SLASM(NEQ_SLA,NEQ_MST))
      ALLOCATE(AK_COND1(NEQ_SLA,NEQ_MST),AK_COND2(NEQ_MST,NEQ_MST),AK_SLASS_INV(NEQ_SLA,NEQ_SLA))  
      
      ALLOCATE(AK_COND_INV(NEQ_MST,NEQ_MST),FORCE_MST(NEQ_MST),FORCE_SLA(NEQ_SLA),FORCE(NEQ),DISP_MST(NEQ_MST))
      ALLOCATE (NSTART(NEQ),NEND(NEQ))
      
      ! CHECK DOFS ARE USED.
      IF (NBS.NE.0)THEN
      NBOUNDARY_NEQ = 0D0
      NEND = 0.
      NSTART = 0.
      DO I = 1,NSN
          DO J = 1,9
          IF (NBOUNDARY(J,I).EQ.0.AND.IDOFB(J).EQ.0) THEN
              NBOUNDARY_NEQ = NBOUNDARY_NEQ + 1
              NEND(I)   = NBOUNDARY_NEQ
              NSTART(I) = NEND(I)  - NDOF_STRUCTURE + 1
          ENDIF
          ENDDO
      ENDDO
      ENDIF
      
      !INITIATE MATRIX
      AK_COND = 0.0D0
      AK_MST = 0.0D0
      AK_SLAMS = 0.0D0
      AK_SLASS = 0.0D0
      AK_SLASM = 0.0D0
      AK_COND1 = 0.0D0      
      AK_COND2 = 0.0D0
      AK_SLASS_INV = 0.0D0
      
      !SAVE STIFFNESS OF MASTER NODES (AK_MST) ---------------------------------------
      DO IMST = 1,NMST
          IF (NBS.EQ.0)THEN
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          MST_NODE = MST_MAT(IMST)
          IDOF1 = (MST_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF2 = IDOF1 + 5 
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF2 = MDOF1 + 5 
         ELSEIF (NBS.NE.0)THEN
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          MST_NODE = MST_MAT(IMST)
          IDOF1 = NSTART(MST_NODE)
          IDOF2 = NEND(MST_NODE)
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF2 = MDOF1 + 5 
         ENDIF
          
         AK_MST(MDOF1:MDOF2,MDOF1:MDOF2) = AKS(IDOF1:IDOF2,IDOF1:IDOF2)
       ENDDO      
         
      !SAVE STIFFNESS OF SLAVE NODES (AK_SLASS) ---------------------------------------
         ! DO ISLA = 1,NSLA
         ! !DEFINE OF POSITION OF SLAVE NODES IN AKS MATRIX
         !  ISLA_NODE = ISLA_MAT(ISLA)
         !  IDOF1 = (ISLA_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
         !  IDOF2 = IDOF1 + 5 
         ! !DEFINE OF POSITION IN SLAVE MATRIX
         !  MDOF1 = (ISLA-1)*6 + 1 !6 = THE NUMBER OF D.O.F
         !  MDOF2 = MDOF1 + 5 
         !  
         ! ! AK_SLASS(MDOF1:MDOF2,MDOF1:MDOF2) = AKS(IDOF1:IDOF2,IDOF1:IDOF2)
         ! ENDDO  
         
      DO ISLA = 1,NSLA
         ! START AND END POSITION OF SLAVE NODE IN AKS MATRIX
          ISLA_NODE = ISLA_MAT(ISLA)
          IDOF1     = (ISLA_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF2     = IDOF1 + 5 
      
          MDOF1 = (ISLA-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF2 = MDOF1 + 5 
          
          DO J = 1,NSLA
          ISLA_NODE = ISLA_MAT(J)
          IDOF3     = (ISLA_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF4     = IDOF3 + 5 
      
          MDOF3 = (J-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF4 = MDOF3 + 5
              
          AK_SLASS(MDOF1:MDOF2,MDOF3:MDOF4) = AKS(IDOF1:IDOF2,IDOF3:IDOF4)
          ENDDO
         ENDDO  
         
!      DO I =1,NEQ_SLA
!          DO J = 1,NEQ_SLA
!      WRITE (299,999) I,J,AK_SLASS(I,J)
!      ENDDO
!      ENDDO
!999   FORMAT (I5,2X,I5,2X,E20.10)      
        !AK_SLASS(1:18,1:18) = AKS(1:18,1:18)
      
      ! UPPER  
      !SAVE STIFFNESS OF MASTER-SLAVE NODES (AK_SLAMS(NEQ_MST,NEQ_SLA)) ---------------------------------------
      DO IMST = 1,NMST
          IF (NBS.EQ.0)THEN
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          MST_NODE = MST_MAT(IMST)
          IDOF_MST1 = (MST_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF_MST2 = IDOF_MST1 + 5 
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF_MST1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_MST2 = MDOF_MST1 + 5 
         ELSEIF (NBS.NE.0)THEN
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          MST_NODE = MST_MAT(IMST)
          IDOF_MST1 = NSTART(MST_NODE)
          IDOF_MST2 = NEND(MST_NODE)
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF_MST1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_MST2 = MDOF_MST1 + 5 
         ENDIF
          
          DO ISLA = 1,NSLA 
          !DEFINE OF POSITION OF SLAVE NODES IN AKS MATRIX
          ISLA_NODE = ISLA_MAT(ISLA)
          IDOF_SLA1 = (ISLA_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF_SLA2 = IDOF_SLA1 + 5 
         !DEFINE OF POSITION IN SLAVE MATRIX
          MDOF_SLA1 = (ISLA-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_SLA2 = MDOF_SLA1 + 5     
          
          AK_SLAMS(MDOF_MST1:MDOF_MST2,MDOF_SLA1:MDOF_SLA2) = AKS(IDOF_MST1:IDOF_MST2,IDOF_SLA1:IDOF_SLA2)
         ENDDO
          ENDDO 
      !AK_SLAMS(1:6,1:126) = AKS(121:126,1:126)
          
      ! LOWER   
      !SAVE STIFFNESS OF SLAVE-MASTER NODES (AK_SLASM(NEQ_SLA,NEQ_MST)) ---------------------------------------
      DO ISLA = 1,NSLA 
          !DEFINE OF POSITION OF SLAVE NODES IN AKS MATRIX
          ISLA_NODE = ISLA_MAT(ISLA)
          IDOF_SLA1 = (ISLA_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF_SLA2 = IDOF_SLA1 + 5 
         !DEFINE OF POSITION IN SLAVE MATRIX
          MDOF_SLA1 = (ISLA-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_SLA2 = MDOF_SLA1 + 5         
          
          DO IMST = 1,NMST
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          IF (NBS.EQ.0)THEN
          MST_NODE = MST_MAT(IMST)
          IDOF_MST1 = (MST_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF_MST2 = IDOF_MST1 + 5 
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF_MST1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_MST2 = MDOF_MST1 + 5   
          ELSEIF (NBS.NE.0)THEN
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          MST_NODE =  MST_MAT(IMST)
          IDOF_MST1 = NSTART(MST_NODE)
          IDOF_MST2 = NEND(MST_NODE)
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF_MST1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_MST2 = MDOF_MST1 + 5 
         ENDIF
            
          !AK_SLASM(MDOF_SLA1:MDOF_SLA2,MDOF_MST1:MDOF_MST2) = AKS(IDOF_SLA1:IDOF_SLA2,IDOF_MST1:IDOF_MST2)
          AK_SLASM(MDOF_SLA1:MDOF_SLA2,MDOF_MST1:MDOF_MST2) = AKS(IDOF_SLA1:IDOF_SLA2,IDOF_MST1:IDOF_MST2)
         ENDDO
      ENDDO 
      !AK_SLASM(1:120,1:6) = AKS(1:120,121:126)    
          
      !CONDENSATION OF STIFFNESS MATRIX      
      !(CONDENSED STIFFNESS) = (MASTER STIFFNESS) - (MASTER-SLAVE STIFFNESS)*Inverse(SLAVE-SLAVE STIFFNESS)*(SLAVE-MASTER STIFFNESS)
          !TRANSFORMATION MATRIX T 
          ! 
          ! T = [                              I
          !       -Inverse(SLAVE-SLAVE STIFFNESS)*(SLAVE-MASTER STIFFNESS) ]   => SIZE OF THE I MATRIX IS SAME WITH MASTER STIFFNESS      
          !=========================================================      
          ! MULTIPLYING MATRICES USING 'DGEMM'
          ! M, N, KIntegers indicating the size of the matrices:
          ! C = A*B
          ! A: M rows by K columns
          ! B: K rows by N columns
          ! C: M rows by N columns
          ! CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
          !=========================================================    
          ALPHA = 1.0
          BETA = 1.0
          AK_COND1 = 0.0D0
          !!AK_COND2 = 0.0D0
          
          !INVERSE OF SLAVE-SLAVE STIFFNESS
          CALL INVMAT(AK_SLASS,AK_SLASS_INV,NEQ_SLA)
          
           
          AK_COND1 = MATMUL(AK_SLASS_INV,AK_SLASM)        
          AK_COND2 = MATMUL(AK_SLAMS,AK_COND1)    
          
          !CALL DGEMM('N','N',NEQ_SLA,NEQ_MST,NEQ_SLA,ALPHA,AK_SLASS_INV,NEQ_SLA,AK_SLASM,NEQ_SLA,BETA,AK_COND1,NEQ_SLA)
          
          !CALL DGEMM('N','N',NEQ_MST,NEQ_MST,NEQ_SLA,ALPHA,AK_SLAMS,NEQ_MST,AK_COND1,NEQ_SLA,BETA,AK_COND2,NEQ_MST)
          
          AK_COND = AK_MST - AK_COND2  !(CONDENSED STIFFNESS)
          
          !IF (AK_COND(6,6).EQ.0.0) = AK_COND(6,6) = AK_COND2(6,6)
          
      !CONDENSED STIFFNESS FOR EACH PILES          
      !WRITE(13003,*) 'The Number of Piles'   
      !WRITE(13003,*) NMST    
      !DO IMST = 1,NMST
      !    MDOF1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
      !    MDOF2 = MDOF1 + 5 
      !    WRITE(13003,100) IMST
      !    WRITE(13003,200) AK_COND(MDOF1:MDOF2,MDOF1:MDOF2)
      !ENDDO    
                  
      !CONDENSED DIAGONAL STIFFNESS FOR EACH PILES
      LENGTH = LEN_TRIM(LHDSTD(ILC))    
      NAMEL  = LHDSTD(ILC)
      NLOCASE = LGEN + IOFFL
      IF (ILC.EQ.1) WRITE(13004,*) '--- PILE-SUPERELEMENT STIFFNESS ---'
      IF (ILC.EQ.1) WRITE(13004,400) NLOCASE
400   FORMAT (I5,2X,"TOTAL LOAD CASE")
      WRITE(13004,'')
      DO IMST = 1,NMST
          MDOF1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF2 = MDOF1 + 5 
          WRITE(13004,100) NAMEL(1:LENGTH)
          WRITE(13004,150) 
          NN = 0
          DO IONE = MDOF1,MDOF2
              DO ITWO = MDOF1,MDOF2 
                  IF(IONE.EQ.ITWO) THEN
                      NN = NN + 1
                      AK_COND_DIAG(NN) = AK_COND(IONE,ITWO)
                  ENDIF
              ENDDO
          ENDDO
          !WRITE(13004,300) IMST,AK_COND_DIAG(1:6)
      ENDDO      
      
      DO I = 1,NDOF_STRUCTURE
          WRITE(13004,300) I,AK_COND(I,1:6)
      ENDDO
      
      FORCE(1:NEQ) = 0.0D0
      FORCE(1) = -344263
      FORCE(2) = 17599
      FORCE(3) = -7603990
      
      FORCE_MST(1:NEQ_MST) = 0.0D0
      FORCE_SLA(1:NEQ_SLA) = 0.0D0
      FORCE_MST(1) =   -297225 
      FORCE_MST(2) =   78850.6 
      FORCE_MST(3) =  -2860530 
      !AK_COND(6,6) = 1d0
      
      !CALCULATE CONDENSED DISPLACEMENT 
          !INVERSE OF CONDENSED STIFFNESS
          CALL INVMAT_BJ(AK_COND,AK_COND_INV,NEQ_MST)
          
          FORCE_MST = FORCE_MST - MATMUL(TRANSPOSE(AK_COND1),FORCE_SLA)
          DISP_MST = MATMUL(AK_COND_INV,FORCE_MST)
                  
          
      DEALLOCATE(AK_COND)
      DEALLOCATE(AK_MST,AK_SLAMS,AK_SLASS,AK_SLASM)
      DEALLOCATE(AK_COND1,AK_COND2,AK_SLASS_INV)    
      DEALLOCATE(ISLA_MAT)
      
      DEALLOCATE(AK_COND_INV,FORCE_MST,FORCE,FORCE_SLA,DISP_MST)
      DEALLOCATE (NSTART,NEND)

100   FORMAT ("NAME OF LOAD CASE:",2X,A)  
150   FORMAT ("DATA NO.    SX               SY               SZ              SRX              SRY              SRZ")
200   FORMAT (E20.10,2X,E20.10,2X,E20.10,2X,E20.10,2X,E20.10,2X,E20.10,2X)      
300   FORMAT (I3,1X,E20.10,2X,E20.10,2X,E20.10,2X,E20.10,2X,E20.10,2X,E20.10,2X)   

      WRITE (13004,"")
      WRITE (13004,350)
350   FORMAT ("*************")
      
      END
C	=============================================================
C	=============================================================
C	=============================================================     
      
C=====================================================================
	SUBROUTINE INVMAT_BJ(A,B,N)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------
C	PURPOSE:	MATRIX INVERSION BY ELIMINATION WITH PARTIAL PIVOTING
C					AND DETERMINANT OF MATRIX A
C
C	INPUT VARIABLES
C     A	= INPUT SQUARE MATRIX
C     N	= SIZE OF A AND B
C
C	LOCAL VARIABLES
C     EPS	= CONTROL VARIABLE
C     DET	= DETERMINANT OF A
C
C	OUTPUT VARIABLE
C     B	= INVERSE MATRIX OF A
C	-----------------------------------------------------------------
      DIMENSION A(N,N),B(N,N),C(N,N)

      EPS=1.0e-20

C	----------------------------------
C     CONSTRUCT IDENTITY MATRIX B(I,J)=I
C	----------------------------------
	DO 88 I=1,N
	 DO 88 J=1,N
	  C(I,J)=A(I,J)
88	CONTINUE

	DO 6 I=1,N
	 DO 5 J=1,N
	  IF(I-J) 4,3,4
3	   B(I,J)=1.0
	   GOTO 5
4	   B(I,J)=0.0
5	 CONTINUE
6	CONTINUE
C	---------------------------------------------------------
C     LOCATE MAXIMUM MAGNITUDE A(I,K) ON OR BELOW MAIN DIAGONAL
C	---------------------------------------------------------
      DET=1.0

	DO 45 K=1,N
	 IF (K-N) 12,30,30
12	 IMAX=K
       AMAX=ABS(C(K,K))
	 KP1=K+1
	 DO 20 I=KP1,N
	   IF (AMAX-ABS(C(I,K))) 15,20,20
15	   IMAX=I
         AMAX=ABS(C(I,K))
20	 CONTINUE
C	 --------------------------------------------------
C      INTERCHANGE ROWS IMAX AND K IF IMAX NOT EQUAL TO K
C	 --------------------------------------------------
       IF (IMAX-K) 25,30,25
25	 DO 29 J=1,N
         ATMP=C(IMAX,J)
	   C(IMAX,J)=C(K,J)
	   C(K,J)=ATMP
	   BTMP=B(IMAX,J)
	   B(IMAX,J)=B(K,J)
29	 B(K,J)=BTMP
       DET=-DET
30	 CONTINUE
C	 ------------------------
C      TEST FOR SINGULAR MATRIX
C	 ------------------------
       IF (ABS(C(K,K))-EPS) 33,33,35
33	 STOP
35	 DET=C(K,K)*DET
C	 ---------------------------------------------
C      DIVIDE PIVOT ROW BY ITS MAIN DIAGONAL ELEMENT
C	 ---------------------------------------------
       DIV=C(K,K)
	 DO 38 J=1,N
	   C(K,J)=C(K,J)/DIV
38	 B(K,J)=B(K,J)/DIV
C	 -----------------------------------------------------
C      REPLACE EACH ROW BY LINEAR COMBINATION WITH PIVOT ROW
C	 -----------------------------------------------------
       DO 43 I=1,N
	   AMULT=C(I,K)
	   IF (I-K) 39,43,39
39	   DO 42 J=1,N
	      C(I,J)=C(I,J)-AMULT*C(K,J)
42	   B(I,J)=B(I,J)-AMULT*B(K,J)
43	 CONTINUE
45	CONTINUE
C
      RETURN
	END
C
C=====================================================================      