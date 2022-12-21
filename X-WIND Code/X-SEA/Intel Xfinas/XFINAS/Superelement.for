      SUBROUTINE SUPERELEMENT (MAXA,ILC)
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
      
      !READ SUBSTRUCTURING DATA
      CALL READSUB(NMST,MST_MAT)      
            
      CALL STIFFSUB(AKS,NEQ,NWK,NMST,MST_MAT)

      
      DEALLOCATE(AA,MAXC_STIFFNESS)
      DEALLOCATE(AKS)
      
      END
C	=============================================================
C	=============================================================
C	============================================================= 
      
      SUBROUTINE STIFFSUB(AKS,NEQ,NWK,NMST,MST_MAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)   
      CHARACTER*200 LHDSTD,LHDCOM,NAMEL
      ALLOCATABLE AK_COND(:,:)
      ALLOCATABLE AK_MST(:,:),AK_SLAMS(:,:),AK_SLASS(:,:),AK_SLASM(:,:)
      ALLOCATABLE AK_COND1(:,:),AK_COND2(:,:),AK_SLASS_INV(:,:) 
      ALLOCATABLE ISLA_MAT(:)
      ALLOCATABLE AK_COND_INV(:,:),FORCE_MST(:),FORCE_SLA(:),FORCE(:),DISP_MST(:)
      COMMON /LHEADER/ LHDSTD(500),LHDCOM(500)
      
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
      
      !REWIND(13001)
      !READ(13001,*) 
      !READ(13001,*) AKS(1:NEQ,1:NEQ)
      
      !DEFINE SIZES OF MASTER MATRIX(NEQ_MST) AND SLAVE MATRIX(NEQ_SLA)
      NEQ_MST = NMST*6 !6 = THE NUMBER OF D.O.F OF MASTER NODES
      NEQ_SLA = NEQ-NEQ_MST
      NUM_TOT = NEQ/6 ! TOTAL NUMBER OF NODES
      !NMST =  THE NUMBER OF MATER NODES
      NSLA = NEQ_SLA/6 ! THE NUMBER OF SALVE NODES
      
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
      
      ALLOCATE(AK_COND(NEQ_MST,NEQ_MST))          
      ALLOCATE(AK_MST(NEQ_MST,NEQ_MST),AK_SLAMS(NEQ_MST,NEQ_SLA),AK_SLASS(NEQ_SLA,NEQ_SLA),AK_SLASM(NEQ_SLA,NEQ_MST))
      ALLOCATE(AK_COND1(NEQ_SLA,NEQ_MST),AK_COND2(NEQ_MST,NEQ_MST),AK_SLASS_INV(NEQ_SLA,NEQ_SLA))  
      
      ALLOCATE(AK_COND_INV(NEQ_MST,NEQ_MST),FORCE_MST(NEQ_MST),FORCE_SLA(NEQ_SLA),FORCE(NEQ),DISP_MST(NEQ_MST))
      
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
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          MST_NODE = MST_MAT(IMST)
          IDOF1 = (MST_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF2 = IDOF1 + 5 
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF2 = MDOF1 + 5 
          
         AK_MST(MDOF1:MDOF2,MDOF1:MDOF2) = AKS(IDOF1:IDOF2,IDOF1:IDOF2)
       ENDDO      
         
      !SAVE STIFFNESS OF SLAVE NODES (AK_SLASS) ---------------------------------------
      DO ISLA = 1,NSLA
         !DEFINE OF POSITION OF SLAVE NODES IN AKS MATRIX
          ISLA_NODE = ISLA_MAT(ISLA)
          IDOF1 = (ISLA_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF2 = IDOF1 + 5 
         !DEFINE OF POSITION IN SLAVE MATRIX
          MDOF1 = (ISLA-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF2 = MDOF1 + 5 
          
         AK_SLASS(MDOF1:MDOF2,MDOF1:MDOF2) = AKS(IDOF1:IDOF2,IDOF1:IDOF2)
      ENDDO  
         
      !SAVE STIFFNESS OF MASTER-SLAVE NODES (AK_SLAMS(NEQ_MST,NEQ_SLA)) ---------------------------------------
      DO IMST = 1,NMST
         !DEFINE OF POSITION OF MASTER NODES IN AKS MATRIX
          MST_NODE = MST_MAT(IMST)
          IDOF_MST1 = (MST_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF_MST2 = IDOF_MST1 + 5 
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF_MST1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_MST2 = MDOF_MST1 + 5 
          
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
          MST_NODE = MST_MAT(IMST)
          IDOF_MST1 = (MST_NODE-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          IDOF_MST2 = IDOF_MST1 + 5 
         !DEFINE OF POSITION IN MASTER MATRIX
          MDOF_MST1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF_MST2 = MDOF_MST1 + 5           
            
          AK_SLASM(MDOF_SLA1:MDOF_SLA2,MDOF_MST1:MDOF_MST2) = AKS(IDOF_SLA1:IDOF_SLA2,IDOF_MST1:IDOF_MST2)
         ENDDO
      ENDDO 
          
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
!      LENGTH = LEN_TRIM(LHDSTD(ILC))    
!      NAMEL  = LHDSTD(ILC)
      IF (ILC.EQ.1) WRITE(13005,*) '--- SUPERELEMENT STIFFNESS ---'
      WRITE(13005,'')
      DO IMST = 1,NMST
          MDOF1 = (IMST-1)*6 + 1 !6 = THE NUMBER OF D.O.F
          MDOF2 = MDOF1 + 5 
!          WRITE(13005,100) NAMEL(1:LENGTH)
          WRITE(13005,150) 
          NN = 0
          DO IONE = MDOF1,MDOF2
              DO ITWO = MDOF1,MDOF2 
                  IF(IONE.EQ.ITWO) THEN
                      NN = NN + 1
                      AK_COND_DIAG(NN) = AK_COND(IONE,ITWO)
                  ENDIF
              ENDDO
          ENDDO
          WRITE(13005,300) IMST,AK_COND_DIAG(1:6)
      ENDDO      
      
      FORCE(1:NEQ) = 0.0D0
      FORCE(1) = -1e+006
      FORCE(2) = 20000
      FORCE(3) = -3e+005 
      
      FORCE_MST(1:NEQ_MST) = 0.0D0
      FORCE_SLA(1:NEQ_SLA) = 0.0D0
      FORCE_MST(1) = -1e+006
      FORCE_MST(2) = 20000
      FORCE_MST(3) = -3e+005     
      
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

100   FORMAT ("NAME OF LOAD CASE:",2X,A)  
150   FORMAT ("DATA NO.    SFX               SY               SZ              SRX              SRY              SRZ")
200   FORMAT (E15.5,2X,E15.5,2X,E15.5,2X,E15.5,2X,E15.5,2X,E15.5,2X)      
300   FORMAT (I3,1X,E15.5,2X,E15.5,2X,E15.5,2X,E15.5,2X,E15.5,2X,E15.5,2X)   

      WRITE (13005,"")
      WRITE (13005,350)
350   FORMAT ("*************")
      
      END
C	=============================================================
C	=============================================================
C	=============================================================   