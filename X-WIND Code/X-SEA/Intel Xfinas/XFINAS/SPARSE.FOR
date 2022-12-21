C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE MCSPARS(MAXA,NEQ,KSC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------------------
C     REFORM THE SIZE OF STIFFNESS BY REMOVE THE ZERO TERM UNDER THE SKYLINE
C	THE NEW MATRIC MAXC (COLUMN INDEX OF EACH COEF.) IS STORE IN FILE 'MAXC'
C	MAXA VECTOR WILL BE MODIFIED HERE
C     ------------------------------------------------------------------------
      DIMENSION MAXA(1)

	ALLOCATABLE LM(:),S(:),MAMC(:),MAXC(:),MCOH(:) 

      ALLOCATABLE LMEMB(:,:)

C	-------------------------------------------------------------------------------------------
C	IMOPT = 0   !SKYLINE + COLSOL
C	IMOPT = 1   !SPARSE  + ITERATION
C	ICPRO = 0   !STIFNESS SIZE REDUCTION PROCESS IS NOT YET DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1)
C	ICPRO = 1   !STIFNESS SIZE REDUCTION PROCESS IS ALREADY DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1)
	CALL MINTFIL('SOLV',IMOPT ,1,1 ,0)  !STIFFNESS PROFILE AND SOLVER OPTION
	CALL MINTFIL('BLOK',ICPRO ,1,4 ,0)  !INDEX FOR STIFFNESS SIZE REDUCTION PROCESS
	IF(IMOPT.EQ.0) RETURN
	IF(ICPRO.EQ.1) RETURN
	ICPRO = 1
	CALL MINTFIL('BLOK',ICPRO ,1,4 ,1)  !UPDATE HERE INDEX FOR STIFFNESS SIZE REDUCTION PROCESS
C	-------------------------------------------------------------------------------------------

	CALL CPU_TIME (TIME1)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	CALL MINTFIL('BLOK',MLINE ,1,3, 0)  !NUMBER OF LINE
      CALL MINTFIL('BLOK',NEFM  ,1,7 ,0)  !MAX ELEMENT NDOF (MAX OF OVERALL GROUP)
      
      NWKM = (NEFM*NEFM+NEFM)/2
	ALLOCATE(LM(NEFM),S(NWKM))

	ALLOCATE(MAMC(MSTOR),MCOH(NEQ))
      
      ALLOCATE(LMEMB(MLINE,NEFM+2))

	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)

C     --------------
C     INITIALISATION
C     --------------
      NEQF = 1
      NEQL = 0
C	NEW STORAGE SIZE
	NSTOR= 0
C	UPDATE VECTOR OF COLUMN HIGH OF EACH EQUATION
	MCOH(1:NEQ) = 0
      
C     FLAG TO CHECK WHETHER ELEMENT HAS BEEN CALL AT LEAST ONCE OR NOT      
      LMEMB = 0
      
C     --------------------------
C     LOOP OVER NUMBER OF BLOCKS
C     --------------------------
      DO 2000  IBLO=1,NBLOCK
      CALL CLEARI (MAMC,MSTOR)

	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
      NEQL   = NEQL+NCOL
      NPRE   = MAXA(NEQF)-1	

	CALL MCALFIL(KEREC,'ETIF')
      REWIND(KEREC)
C     -----------------------------------------------
C     LOOP OVER NUMBER OF LINE RECORDS
C     SCAN ELEMENT EQUATION NUMBERS FOR CONTRIBUTIONS
C     -----------------------------------------------
	DO 500 MELM = 1,MLINE
       
      ICOUNT = LMEMB(MELM,1)
      
      IF(ICOUNT.EQ.0) THEN  
	    CALL MESTIF(LM,S,NWK,NEF,'REL')
          
          CALL HEAPSORI(NEF,LM)
          NEFL = 0
          DO IEF = 1,NEF
              IF(LM(IEF).GT.0) THEN
                  NEFL = NEF-IEF+1
                  LM(1:NEFL) = LM(IEF:NEF)
                  EXIT
              ENDIF
          ENDDO
          NEF = NEFL

          LMEMB(MELM,1) = 1
          LMEMB(MELM,2) = NEFL
          LMEMB(MELM,2+1:2+NEFL) = LM(1:NEFL)

      ELSE
          NEF = LMEMB(MELM,2)
          LM(1:NEF) = LMEMB(MELM,2+1:2+NEF)
      ENDIF
      
C	------------------------------
	DO 120 IEF = 1,NEF
	IEQ = LM(IEF)
	IF(IEQ.LE.0) GOTO 120

      IF(IEQ.GE.NEQF .AND. IEQ.LE.NEQL)  GOTO 50
	GOTO 120
50    CONTINUE
      
	MC = MCOH(IEQ)
C	----------------------
	DO 110 JEF = IEF,NEF
	JEQ = LM(JEF)
	IF(JEQ.LE.0  ) GOTO 110
	IF(JEQ.LT.IEQ) GOTO 110
      
      IF (MC.EQ.0) GOTO 100

      II = MC + MAXA(IEQ) - NPRE - 1
      IF(JEQ.EQ.MAMC(II)) GOTO 110
      IF(JEQ.GT.MAMC(II)) GOTO 100
      
	DO IC = 1,MC
	    II = IC + MAXA(IEQ) - NPRE - 1
	    IF(JEQ.EQ.MAMC(II)) GOTO 110
          IF(JEQ.LT.MAMC(II)) THEN
              JJ = MC + MAXA(IEQ) - NPRE - 1
              MAMC(II+1:JJ+1) = MAMC(II:JJ)
              MAMC(II) = JEQ
              MC = MC + 1
	        MCOH(IEQ) = MC
              GOTO 110
          ENDIF
	ENDDO

100	MC = MC + 1
	II = MC + MAXA(IEQ) - NPRE - 1
	MAMC(II ) = JEQ
	MCOH(IEQ) = MC

110	CONTINUE
C	----------------------

120	CONTINUE
C	------------------------------

C     -----------------------------------------------
500	CONTINUE
C     -----------------------------------------------


C	-------------------------------
C     SORT & UPDATE MAX BLOCK STORAGE SIZE    
      NHIG = SUM(MCOH(NEQF:NEQL))   
	IF(NHIG.GT.NSTOR) NSTOR = NHIG
C	-------------------------------
      
      ALLOCATE(MAXC(NHIG))
C	-------------------------------
C	SORT COLUMN INDEX FROM IEQ TO UPPER COLUMN
C	-------------------------------
      NHIG = 0
	DO IEQ = NEQF,NEQL
          MC = MCOH(IEQ)
	    II = MAXA(IEQ) - NPRE - 1
          MAXC(NHIG+1:NHIG+MC) = MAMC(II+1:II+MC)
      
          NHIG = NHIG+MC
	ENDDO
C	-------------------------------
	WRITE(NFLCH) MAXC(1:NHIG)
C	-------------------------------
	DEALLOCATE(MAXC)

	CALL MINTFIL('BLOC',NHIG,3,IBLO,1)  !UPDATE STORAGE SIZE OF CURRENT BLOCK

	NEQF = NEQF + NCOL
C     -----------------------------------------------
2000  CONTINUE
C     -----------------------------------------------

C	-----------------------------------------------
C	UPDATE MAXA VECTOR
C	-----------------------------------------------
	MAXA(1) = 1
	DO IEQ = 1,NEQ
	MAXA(IEQ+1) = MAXA(IEQ) + MCOH(IEQ)
	ENDDO
C	-----------------------------------------------

C	-----------------------------------------------
	MSTOR = NSTOR
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,1)   !UPDATE STORAGE SIZE
C	-----------------------------------------------
	NWK = MAXA(NEQ+1)-1
	NWM = MAXA(NEQ+1)-1
C	-----------------------------------------------


	DEALLOCATE(LM,S)
	DEALLOCATE(MAMC)
	DEALLOCATE(MCOH)

	DEALLOCATE(LMEMB)

      IF (KSC.NE.1)THEN
	WRITE(*,'(/X,A34)')     'MEMORY SIZE AFTER SPARSE REDUCTION'
	WRITE(*,'( X,A33,I10)') 'NUMBER OF DOF . . . . . . . . . .',NEQ
	WRITE(*,'( X,A33,I10)') 'SIZE OF STIFFNESS MATRIX. . . . .',NWK
	WRITE(*,'( X,A33,I10)') 'SIZE OF MASS MATRIX . . . . . . .',NWM
	WRITE(*,'( X,A33,I10/)') 'SIZE OF IN-CORE STORAGE BLOCK . .',MSTOR
	ENDIF

c	WRITE(110,'(/X,A34)')     'MEMORY SIZE AFTER SPARSE REDUCTION'
c	WRITE(110,'( X,A33,I10)') 'NUMBER OF DOF . . . . . . . . . .',NEQ
c	WRITE(110,'( X,A33,I10)') 'SIZE OF STIFFNESS MATRIX. . . . .',NWK
c	WRITE(110,'( X,A33,I10)') 'SIZE OF MASS MATRIX . . . . . . .',NWM
c	WRITE(110,'( X,A33,I10/)') 'SIZE OF IN-CORE STORAGE BLOCK . .',MSTOR

      RETURN
      END

C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE MCSPARS_WORK_WELL_2(MAXA,NEQ,KSC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------------------
C     REFORM THE SIZE OF STIFFNESS BY REMOVE THE ZERO TERM UNDER THE SKYLINE
C	THE NEW MATRIC MAXC (COLUMN INDEX OF EACH COEF.) IS STORE IN FILE 'MAXC'
C	MAXA VECTOR WILL BE MODIFIED HERE
C     ------------------------------------------------------------------------
      DIMENSION MAXA(1)

	ALLOCATABLE LM(:),S(:),MAMC(:),MAXC(:),MCOH(:) 
      ALLOCATABLE DUMMY(:)
      ALLOCATABLE LMEMB(:,:)

C	-------------------------------------------------------------------------------------------
C	IMOPT = 0   !SKYLINE + COLSOL
C	IMOPT = 1   !SPARSE  + ITERATION
C	ICPRO = 0   !STIFNESS SIZE REDUCTION PROCESS IS NOT YET DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1)
C	ICPRO = 1   !STIFNESS SIZE REDUCTION PROCESS IS ALREADY DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1)
	CALL MINTFIL('SOLV',IMOPT ,1,1 ,0)  !STIFFNESS PROFILE AND SOLVER OPTION
	CALL MINTFIL('BLOK',ICPRO ,1,4 ,0)  !INDEX FOR STIFFNESS SIZE REDUCTION PROCESS
	IF(IMOPT.EQ.0) RETURN
	IF(ICPRO.EQ.1) RETURN
	ICPRO = 1
	CALL MINTFIL('BLOK',ICPRO ,1,4 ,1)  !UPDATE HERE INDEX FOR STIFFNESS SIZE REDUCTION PROCESS
C	-------------------------------------------------------------------------------------------

	ALLOCATE(LM(200),S(21000))

	CALL CPU_TIME (TIME1)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	CALL MINTFIL('BLOK',MLINE ,1,3, 0)  !NUMBER OF LINE
      CALL MINTFIL('BLOK',NEFM  ,1,7 ,0)  !MAX ELEMENT NDOF (MAX OF OVERALL GROUP)

	ALLOCATE(MAMC(MSTOR),MAXC(MSTOR),MCOH(NEQ))
      
      ALLOCATE(DUMMY(NEQ))
      
      ALLOCATE(LMEMB(MLINE,NEFM+2))

	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)

C     --------------
C     INITIALISATION
C     --------------
      NEQF = 1
      NEQL = 0
C	NEW STORAGE SIZE
	NSTOR= 0
C	UPDATE VECTOR OF COLUMN HIGH OF EACH EQUATION
	MCOH(1:NEQ) = 0
      
C     FLAG TO CHECK WHETHER ELEMENT HAS BEEN CALL AT LEAST ONCE OR NOT      
      LMEMB = 0
      
C     --------------------------
C     LOOP OVER NUMBER OF BLOCKS
C     --------------------------
      DO 2000  IBLO=1,NBLOCK
      CALL CLEARI (MAMC,MSTOR)

	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
      NEQL   = NEQL+NCOL
      NPRE   = MAXA(NEQF)-1	

	CALL MCALFIL(KEREC,'ETIF')
	REWIND(KEREC)
C     -----------------------------------------------
C     LOOP OVER NUMBER OF LINE RECORDS
C     SCAN ELEMENT EQUATION NUMBERS FOR CONTRIBUTIONS
C     -----------------------------------------------
	DO 500 MELM = 1,MLINE
       
      ICOUNT = LMEMB(MELM,1)
      
      IF(ICOUNT.EQ.0) THEN  
	    CALL MESTIF(LM,S,NWK,NEF,'REL')
          
          IF(NEF.LE.NEFM) THEN
              LMEMB(MELM,1) = 1
              LMEMB(MELM,2) = NEF
          
              LMEMB(MELM,2+1:2+NEF) = LM(1:NEF)
          ENDIF
      ELSE
          NEF = LMEMB(MELM,2)
          LM(1:NEF) = LMEMB(MELM,2+1:2+NEF)
      ENDIF
      

      
C	------------------------------
	DO 120 IEF = 1,NEF
	IEQ = LM(IEF)
	IF(IEQ.LE.0) GOTO 120

      IF(IEQ.GE.NEQF .AND. IEQ.LE.NEQL)  GOTO 50
	GOTO 120
50    CONTINUE
      
	MC = MCOH(IEQ)
C	----------------------
	DO 110 JEF = 1,NEF
	JEQ = LM(JEF)
	IF(JEQ.LE.0  ) GOTO 110
	IF(JEQ.LT.IEQ) GOTO 110

	DO IC = 1,MC
	II = IC + MAXA(IEQ) - NPRE - 1
	IF(JEQ.EQ.MAMC(II)) GOTO 110
	ENDDO

	MC = MC + 1
	II = MC + MAXA(IEQ) - NPRE - 1
	MAMC(II ) = JEQ
	MCOH(IEQ) = MC

110	CONTINUE
C	----------------------

120	CONTINUE
C	------------------------------

C     -----------------------------------------------
500	CONTINUE
C     -----------------------------------------------


C	-------------------------------
C	SORT COLUMN INDEX FROM IEQ TO UPPER COLUMN
C	-------------------------------
      NHIG = 0
      
	DO IEQ = NEQF,NEQL

      MC = MCOH(IEQ)
	II = MAXA(IEQ) - NPRE - 1
      CALL HEAPSORI(MC,MAMC(II+1))
      
      II = MAXA(IEQ) - NPRE - 1
      MAXC(NHIG+1:NHIG+MC) = MAMC(II+1:II+MC)
      
      NHIG = NHIG+MC
      
	ENDDO
C	-------------------------------
	
	IF(NHIG.GT.NSTOR) NSTOR = NHIG

	WRITE(NFLCH) MAXC(1:NHIG)

	CALL MINTFIL('BLOC',NHIG,3,IBLO,1)  !UPDATE STORAGE SIZE OF CURRENT BLOCK

	NEQF = NEQF + NCOL
C     -----------------------------------------------
2000	CONTINUE
C     -----------------------------------------------

C	-----------------------------------------------
C	UPDATE MAXA VECTOR
C	-----------------------------------------------
	MAXA(1) = 1
	DO IEQ = 1,NEQ
	MAXA(IEQ+1) = MAXA(IEQ) + MCOH(IEQ)
	ENDDO
C	-----------------------------------------------

C	-----------------------------------------------
	MSTOR = NSTOR
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,1)   !UPDATE STORAGE SIZE
C	-----------------------------------------------
	NWK = MAXA(NEQ+1)-1
	NWM = MAXA(NEQ+1)-1




	DEALLOCATE(LM,S)
	DEALLOCATE(MAMC)
	DEALLOCATE(MAXC)
	DEALLOCATE(MCOH)

	DEALLOCATE(DUMMY)

	DEALLOCATE(LMEMB)

      IF (KSC.NE.1)THEN
	WRITE(*,'(/X,A34)')     'MEMORY SIZE AFTER SPARSE REDUCTION'
	WRITE(*,'( X,A33,I10)') 'NUMBER OF DOF . . . . . . . . . .',NEQ
	WRITE(*,'( X,A33,I10)') 'SIZE OF STIFFNESS MATRIX. . . . .',NWK
	WRITE(*,'( X,A33,I10)') 'SIZE OF MASS MATRIX . . . . . . .',NWM
	WRITE(*,'( X,A33,I10/)') 'SIZE OF IN-CORE STORAGE BLOCK . .',MSTOR
	ENDIF

c	WRITE(110,'(/X,A34)')     'MEMORY SIZE AFTER SPARSE REDUCTION'
c	WRITE(110,'( X,A33,I10)') 'NUMBER OF DOF . . . . . . . . . .',NEQ
c	WRITE(110,'( X,A33,I10)') 'SIZE OF STIFFNESS MATRIX. . . . .',NWK
c	WRITE(110,'( X,A33,I10)') 'SIZE OF MASS MATRIX . . . . . . .',NWM
c	WRITE(110,'( X,A33,I10/)') 'SIZE OF IN-CORE STORAGE BLOCK . .',MSTOR

      RETURN
      END

C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE MCSPARS_BACKUP_WORK_WELL (MAXA,NEQ,KSC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------------------
C     REFORM THE SIZE OF STIFFNESS BY REMOVE THE ZERO TERM UNDER THE SKYLINE
C	THE NEW MATRIC MAXC (COLUMN INDEX OF EACH COEF.) IS STORE IN FILE 'MAXC'
C	MAXA VECTOR WILL BE MODIFIED HERE
C     ------------------------------------------------------------------------
      DIMENSION MAXA(1)

	ALLOCATABLE LM(:),S(:),MAMC(:),MAXC(:),MCOH(:) 

C	-------------------------------------------------------------------------------------------
C	IMOPT = 0   !SKYLINE + COLSOL
C	IMOPT = 1   !SPARSE  + ITERATION
C	ICPRO = 0   !STIFNESS SIZE REDUCTION PROCESS IS NOT YET DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1)
C	ICPRO = 1   !STIFNESS SIZE REDUCTION PROCESS IS ALREADY DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1)
	CALL MINTFIL('SOLV',IMOPT ,1,1 ,0)  !STIFFNESS PROFILE AND SOLVER OPTION
	CALL MINTFIL('BLOK',ICPRO ,1,4 ,0)  !INDEX FOR STIFFNESS SIZE REDUCTION PROCESS
	IF(IMOPT.EQ.0) RETURN
	IF(ICPRO.EQ.1) RETURN
	ICPRO = 1
	CALL MINTFIL('BLOK',ICPRO ,1,4 ,1)  !UPDATE HERE INDEX FOR STIFFNESS SIZE REDUCTION PROCESS
C	-------------------------------------------------------------------------------------------

	ALLOCATE(LM(200),S(21000))

	CALL CPU_TIME (TIME1)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	CALL MINTFIL('BLOK',MLINE ,1,3, 0)  !NUMBER OF LINE

	ALLOCATE(MAMC(MSTOR),MAXC(MSTOR),MCOH(NEQ))

	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)

C     --------------
C     INITIALISATION
C     --------------
      NEQF = 1
      NEQL = 0
C	NEW STORAGE SIZE
	NSTOR= 0
C	UPDATE VECTOR OF COLUMN HIGH OF EACH EQUATION
	MCOH(1:NEQ) = 0
C     --------------------------
C     LOOP OVER NUMBER OF BLOCKS
C     --------------------------
      DO 2000  IBLO=1,NBLOCK
      CALL CLEARI (MAMC,MSTOR)

	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
      NEQL   = NEQL+NCOL
      NPRE   = MAXA(NEQF)-1	

	CALL MCALFIL(KEREC,'ETIF')
	REWIND(KEREC)
C     -----------------------------------------------
C     LOOP OVER NUMBER OF LINE RECORDS
C     SCAN ELEMENT EQUATION NUMBERS FOR CONTRIBUTIONS
C     -----------------------------------------------
	DO 500 MELM = 1,MLINE
	CALL MESTIF(LM,S,NWK,NEF,'REL')

C	------------------------------
	DO 120 IEF = 1  ,NEF
	IEQ = LM(IEF)
	IF(IEQ.LE.0) GOTO 120

      IF(IEQ.GE.NEQF .AND. IEQ.LE.NEQL)  GOTO 50
	GOTO 120
50	CONTINUE

	MC = MCOH(IEQ)
C	----------------------
	DO 110 JEF = 1,NEF
	JEQ = LM(JEF)
	IF(JEQ.LE.0  ) GOTO 110
	IF(JEQ.LT.IEQ) GOTO 110

	DO IC = 1,MC
	II = IC + MAXA(IEQ) - NPRE - 1
	IF(JEQ.EQ.MAMC(II)) GOTO 110
	ENDDO

	MC = MC + 1
	II = MC + MAXA(IEQ) - NPRE - 1
	MAMC(II ) = JEQ
	MCOH(IEQ) = MC

110	CONTINUE
C	----------------------

120	CONTINUE
C	------------------------------

C     -----------------------------------------------
500	CONTINUE
C     -----------------------------------------------


C	-------------------------------
C	SORT COLUMN INDEX FROM IEQ TO UPPER COLUMN
C	-------------------------------

	DO IEQ = NEQF,NEQL

      GOTO 510
      
	MC = MCOH(IEQ)
	DO IC = 1,MC
	II = IC + MAXA(IEQ) - NPRE - 1
	NCI = MAMC(II)
	DO JC = IC,MC
	JJ = JC + MAXA(IEQ) - NPRE - 1
	NCJ = MAMC(JJ)
	IF(NCJ.LT.NCI) THEN
	MAMC(JJ) = NCI
	NCI      = NCJ 
	ENDIF
	ENDDO
	MAMC(II) = NCI
	ENDDO

510	MC = MCOH(IEQ)
	II = MAXA(IEQ) - NPRE
      CALL HEAPSORI(MC,MAMC(II))

C      DUMMY(1:MC) = REAL(MAMC(II+1:II+MC))
C      CALL DLASRT('I',MC,DUMMY,INFO)
C      MAMC(II+1:II+MC) = INT(DUMMY(1:MC))      
      
	ENDDO
C	-------------------------------


	NHIG = 0
	DO IEQ = NEQF,NEQL
	MC = MCOH(IEQ)
	DO IC = 1,MC
	NHIG = NHIG + 1
	II = IC + MAXA(IEQ) - NPRE - 1
	MAXC(NHIG) = MAMC(II)
	ENDDO
	ENDDO	
	IF(NHIG.GT.NSTOR) NSTOR = NHIG

	WRITE(NFLCH) MAXC(1:NHIG)

	CALL MINTFIL('BLOC',NHIG,3,IBLO,1)  !UPDATE STORAGE SIZE OF CURRENT BLOCK

	NEQF = NEQF + NCOL
C     -----------------------------------------------
2000	CONTINUE
C     -----------------------------------------------

C	-----------------------------------------------
C	UPDATE MAXA VECTOR
C	-----------------------------------------------
	MAXA(1) = 1
	DO IEQ = 1,NEQ
	MAXA(IEQ+1) = MAXA(IEQ) + MCOH(IEQ)
	ENDDO
C	-----------------------------------------------

C	-----------------------------------------------
	MSTOR = NSTOR
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,1)   !UPDATE STORAGE SIZE
C	-----------------------------------------------
	NWK = MAXA(NEQ+1)-1
	NWM = MAXA(NEQ+1)-1




	DEALLOCATE(LM,S)
	DEALLOCATE(MAMC)
	DEALLOCATE(MAXC)
	DEALLOCATE(MCOH)



      IF (KSC.NE.1)THEN
	WRITE(*,'(/X,A34)')     'MEMORY SIZE AFTER SPARSE REDUCTION'
	WRITE(*,'( X,A33,I10)') 'NUMBER OF DOF . . . . . . . . . .',NEQ
	WRITE(*,'( X,A33,I10)') 'SIZE OF STIFFNESS MATRIX. . . . .',NWK
	WRITE(*,'( X,A33,I10)') 'SIZE OF MASS MATRIX . . . . . . .',NWM
	WRITE(*,'( X,A33,I10/)') 'SIZE OF IN-CORE STORAGE BLOCK . .',MSTOR
	ENDIF

c	WRITE(110,'(/X,A34)')     'MEMORY SIZE AFTER SPARSE REDUCTION'
c	WRITE(110,'( X,A33,I10)') 'NUMBER OF DOF . . . . . . . . . .',NEQ
c	WRITE(110,'( X,A33,I10)') 'SIZE OF STIFFNESS MATRIX. . . . .',NWK
c	WRITE(110,'( X,A33,I10)') 'SIZE OF MASS MATRIX . . . . . . .',NWM
c	WRITE(110,'( X,A33,I10/)') 'SIZE OF IN-CORE STORAGE BLOCK . .',MSTOR

      RETURN
      END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE BLOCKM(MAXA,MAXC,A,B,C,IEQF,IEQL,NEQ,NRB,IND,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 IND,OPER
C	=============================================================
	DIMENSION MAXA(1),MAXC(1),A(1),B(1),C(1)  !B(NEQ,NRB),C(NEQ,NRB)
C	=============================================================
	IOPT = -1
	IF(IND.EQ.'CAB') IOPT = 0
	IF(IND.EQ.'CBA') IOPT = 1
	IF(IOPT.EQ.-1) RETURN

C	--------------------------------------------------------------
	SELECTCASE(IOPT)
C	--------------------------------------------------------------

	CASE(0)
C	C(NEQ,NRB) = A(NEQ,NCB)*B(NCB,NRB)  ! [C] = [A]*[B]
C	NCB = IEQF-IEQL

	NPRE = MAXA(IEQF) - 1
	
	DO 1000 IRB = 1,NRB

	DO IEQ = IEQF,IEQL
	NCOH = MAXA(IEQ+1) - MAXA(IEQ)
	DO ICOH = 1,NCOH
	INDX = MAXA(IEQ) + ICOH - 1 - NPRE
	ICOL = MAXC(INDX)
	VALV = A(INDX)

	C(IEQ+NEQ*(IRB-1)) = C(IEQ+NEQ*(IRB-1)) + VALV*B(ICOL+NEQ*(IRB-1))
 	IF(ICOH.GT.1) C(ICOL+NEQ*(IRB-1)) = C(ICOL+NEQ*(IRB-1)) + 
	1									VALV*B(IEQ+NEQ*(IRB-1))
C	C(IEQ,IRB) = C(IEQ,IRB) + VALV*B(ICOL,IRB)
C	IF(ICOH.GT.1) C(ICOL,IRB) = C(ICOL,IRB) + VALV*B(IEQ,IRB)

	ENDDO
	ENDDO

1000	CONTINUE

C	--------------------------------------------------------------
C	--------------------------------------------------------------

	CASE(1)
C	C(NRB,NEQ) = B(NCB,NRB)*A(NCB,NEQ)  ! [C] = [B]T*[A]
C	NCB = IEQF-IEQL

	NPRE = MAXA(IEQF) - 1

	DO 2000 IRB = 1,NRB

	DO IEQ = IEQF,IEQL
	NCOH = MAXA(IEQ+1) - MAXA(IEQ)
	DO ICOH = 1,NCOH
	INDX = MAXA(IEQ) + ICOH - 1 - NPRE
	ICOL = MAXC(INDX)
	VALV = A(INDX)

	C(IRB+NRB*(IEQ-1)) = C(IRB+NRB*(IEQ-1)) + VALV*B(ICOL+NEQ*(IRB-1))
	IF(ICOH.GT.1) C(IRB+NRB*(ICOL-1)) = C(IRB+NRB*(ICOL-1)) + 
     1									VALV*B(IEQ+NEQ*(IRB-1))
C	C(IRB,IEQ) = C(IRB,IEQ) + VALV*B(ICOL,IRB)
C	IF(ICOH.GT.1) C(IRB,ICOL) = C(IRB,ICOL) + VALV*B(IEQ,IRB)

	ENDDO
	ENDDO

2000	CONTINUE

C	--------------------------------------------------------------
	ENDSELECT
C	--------------------------------------------------------------


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE MASULT (MAXA,B,RR,TT,TYP,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP
      CHARACTER*3 OPER

C     ----------------------------------------------------------------
C     EVALUATES MATRIX PRODUCT OF {TT} =        [B]*{RR} IF OPER='STD' 
C     EVALUATES MATRIX PRODUCT OF {TT} = {TT} + [B]*{RR} IF OPER='ADD' 
C	----------------------------------------------------------------
      DIMENSION MAXA(1),B(1),RR(1),TT(1)
	ALLOCATABLE MAXC(:)
C	-------------------------------------------
	CALL MCALFIL(NFL,TYP)
	REWIND(NFL)

	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH) 

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)


      NEQF = 1
      NEQL = 0
	DO 800 IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
      NEQL = NEQL+NCOL
      
	IF(OPER.EQ.'STD') TT(NEQF:NEQL) = 0.0D0  !IF(OPER.EQ.'ADD') ---> NO INITIALIZATION   

	NEQF = NEQF + NCOL
800	CONTINUE
C     --------------------------------------



      NEQF = 1
      NEQL = 0
	DO 1000 IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
      NEQL = NEQL+NCOL
	NPRE = MAXA(NEQF) - 1

	ALLOCATE(MAXC(NHIG))


	READ(NFL  )    B(1:NHIG)
	READ(NFLCH) MAXC(1:NHIG)
C     --------------------------------------
C     LOOP OVER NUMBER OF COLUMNS IN BLOCK
C     MULTIPLY COLUMN COEFFICIENTS WITH (RR)
C     --------------------------------------
	NRB = 1
	NEQ = 1  !NEQ IS NOT EFFECTED IF NRB = 1

	CALL BLOCKM(MAXA,MAXC,B,RR,TT,NEQF,NEQL,NEQ,NRB,'CAB',OPER)  !SPARSE MULTIPLICATION

	NEQF = NEQF + NCOL
	DEALLOCATE(MAXC)
1000	CONTINUE
C     --------------------------------------

C
      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE ASSEMS (LM,MAXA,A,S,NEQF,NEQL,NPRE,NEF,MAXC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C     ----------------------------------------------------------------
C     ASSEMBLES UPPER TRIANGULAR ELEMENT STIFFNESS INTO COMPACTED
C     GLOBAL STIFFNESS BLOCK (CONTRIBUTIONS BETWEEN NEQF AND NEQL)
C
C     LM(NEF)       = EQUATION NUMBERS FOR ELEMENT DEGREE OF FREEDOMS
C     MAXA(NEQ1)    = ADDRESSES OF DIAGONAL ELEMENTS IN A
C     MAXC(ISTOR)   = COLUMN INDEX NUMBER
C     A(ISTOR)      = GLOBAL COMPACTED STIFFNESS BLOCK
C     S(NWS)        = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     NEQF,NEQL     = FIRST AND LAST EQUATION CONTAINED IN BLOCK
C     NPRE          = NUMBER OF PREVIOUS ELEMENTS IN A
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ----------------------------------------------------------------
C
      DIMENSION LM(1),MAXA(1),A(1),S(1),MAXC(1)
C


      NDI = 0
      DO 200  IEF=1,NEF
      II = LM(IEF)
      IF (II.LT.NEQF .OR. II.GT.NEQL)  GOTO 200
      MI = MAXA(II) - NPRE
      IEF1 = IEF-1
      NIE = 1 + IEF1*NEF - (IEF1*IEF1-IEF1)/2

      DO 220  JEF=1,NEF
      JJ = LM(JEF)
      IF (JJ)  220,220,110
 110  IJ = JJ-II
      IF (IJ)  220,205,205

C	------------------------------------
C	SEARCH FOR COLUMN NUMBER
C	------------------------------------
 205	NCOL = MAXA(II+1) - MAXA(II)
	IJ = 0
	DO ICOL = 1,NCOL
	IF(JJ.EQ.MAXC(MI+ICOL-1)) THEN 	
	IJ = ICOL-1
	GOTO 210
	ENDIF
	ENDDO
	GOTO 220
C	------------------------------------

 210  KK = MI+IJ
      KSS = NIE+JEF-IEF      
      IF (JEF.LT.IEF) THEN
      JEF1 = JEF-1
      NJE = 1 + JEF1*NEF - (JEF1*JEF1-JEF1)/2
      KSS = NJE+IEF-JEF
      ENDIF
      A(KK) = A(KK)+S(KSS)
 220  CONTINUE


 200  CONTINUE
C
C
      RETURN

      END
C	=============================================================
C	=============================================================
C	=============================================================	
	SUBROUTINE CONJUGS(MAXA,A,B,NEQ,IND,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C	=============================================================
C	CONJUGATE GRADIENT ITERATION WITHOUT PRECONDITIONER
	DIMENSION MAXA(1),MAXC(1),A(1),B(1)
C	=============================================================
	ALLOCATABLE X(:),P(:),R(:),V(:)


	IF(IND.EQ.1) THEN
	CALL MDMOVE(TYP1,TYP2)
	RETURN
	ENDIF

	ALLOCATE(X(NEQ),P(NEQ),R(NEQ),V(NEQ))


	CALL CPU_TIME (TIME1)

C	------------------------------------
	DO IEQ = 1,NEQ
	IF(B(IEQ).NE.0.0D0) GOTO 100
	ENDDO
	DEALLOCATE(X,P,R,V)
	RETURN
C	------------------------------------

 100	TOL = 1.0E-2
	X(1:NEQ) = 0.0D0
	P(1:NEQ) = 0.0D0
	R(1:NEQ) = 0.0D0

C	------------------------------------
C	{R} = {B} - [A]*{X}
	CALL MASULT(MAXA,A,X,P,TYP2,'STD')
	R(1:NEQ) = B(1:NEQ) - P(1:NEQ)
C	{P} = {R}
	P(1:NEQ) = R(1:NEQ)
C	------------------------------------



	ITR = 0
C	----------------------------------------
2000	CONTINUE
	ITR = ITR + 1

C	{V} = [A]*{P}
	CALL MASULT(MAXA,A,P,V,TYP2,'STD')
	
	ALPHA1 = 0.0D0
	ALPHA2 = 0.0D0
	DO IEQ = 1,NEQ
C	{ALP1} = {R}T*{R}
	ALPHA1 = ALPHA1 + R(IEQ)*R(IEQ)
C	{ALP2} = {P}T*{R}
	ALPHA2 = ALPHA2 + P(IEQ)*V(IEQ)
	ENDDO

	ALPHA3 = 0.0D0
	DO IEQ = 1,NEQ
	X(IEQ) = X(IEQ) + ALPHA1*P(IEQ)/ALPHA2 
	R(IEQ) = R(IEQ) - ALPHA1*V(IEQ)/ALPHA2
C	{ALP3} = {R i+1}T*{R i+1}
	ALPHA3 = ALPHA3 + R(IEQ)*R(IEQ)
	ENDDO
	BETA = ALPHA3/ALPHA1


	TEST = SQRT(ALPHA1)

	IF(TEST.LT.TOL) GOTO 3500
	IF(ITR.GT.5000) THEN
	WRITE(*,7000) 
	STOP
	ENDIF

	DO IEQ = 1,NEQ
	P(IEQ) = R(IEQ) + BETA*P(IEQ)
	ENDDO
	
	GOTO 2000

C	----------------------------------------
 3500	CONTINUE

	B(1:NEQ) = X(1:NEQ)

	CALL CPU_TIME (TIME2)
C	WRITE(*,6000) NEQ,ITR,TIME2-TIME1


	DEALLOCATE(X,P,R,V)


 6000	FORMAT (X,'SOLVING.',I6,' EQNS. WITH NO. ITR.',I5,
	1		'  TIME TO PERFORM. .',E15.6)
 7000	FORMAT (X,'**** ITERATION PROCESS NOT CONVERGE ****')


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================	
	SUBROUTINE CONJUGJ(MAXA,A,B,NEQ,IND,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C	=============================================================
C	CONJUGATE GRADIENT ITERATION WITH JACOBI PRECONDITIONER
	DIMENSION MAXA(1),MAXC(1),A(1),B(1)
C	=============================================================
	ALLOCATABLE X(:),P(:),R(:),V(:)
	ALLOCATABLE D(:),Z(:)


	IF(IND.EQ.1) THEN
	CALL MDMOVE(TYP1,TYP2)
	RETURN
	ENDIF

	ALLOCATE(X(NEQ),P(NEQ),R(NEQ),V(NEQ),D(NEQ),Z(NEQ))

	CALL CPU_TIME (TIME1)

C	------------------------------------
	DO IEQ = 1,NEQ
	IF(B(IEQ).NE.0.0D0) GOTO 100
	ENDDO
	DEALLOCATE(X,P,R,V,D,Z)
	RETURN
C	------------------------------------

 100	TOL = 1.0E-2
	X(1:NEQ) = 0.0D0
	P(1:NEQ) = 0.0D0
	R(1:NEQ) = 0.0D0

C	------------------------------------
C	CALL DIAGONAL TERM OF STIFFNESS [C]
	CALL MDCALL_DIA(TYP2,D,'RED',MAXA)


	DO IEQ = 1,NEQ

	IF(D(IEQ).EQ.0.0D0) D(IEQ) = 1.0E-6
	IF(ABS(D(IEQ)).LT.1.0E-6) D(IEQ) = 1.0E-6*D(IEQ)/ABS(D(IEQ))

	D(IEQ) = 1.0D0/D(IEQ)
	X(IEQ) = D(IEQ)*B(IEQ)

	ENDDO
C	------------------------------------
C	{R} = {B} - [A]*{X}
	CALL MASULT(MAXA,A,X,P,TYP2,'STD')
	R(1:NEQ) = B(1:NEQ) - P(1:NEQ)
C	{P} = {Z} = [C-1]*{R}
	DO IEQ = 1,NEQ
	Z(IEQ) = R(IEQ)*D(IEQ)
	P(IEQ) = Z(IEQ)
	ENDDO
C	------------------------------------


	ITR = 0
C	----------------------------------------
2000	CONTINUE
	ITR = ITR + 1

C	{V} = [A]*{P}
	CALL MASULT(MAXA,A,P,V,TYP2,'STD')
	
	ALPHA1 = 0.0D0
	ALPHA2 = 0.0D0
	DO IEQ = 1,NEQ
C	{ALP1} = {R}T*{R}
	ALPHA1 = ALPHA1 + Z(IEQ)*R(IEQ)
C	{ALP2} = {P}T*{R}
	ALPHA2 = ALPHA2 + P(IEQ)*V(IEQ)
	ENDDO


	TEST   = 0.0D0
	ALPHA3 = 0.0D0
	DO IEQ = 1,NEQ
	X(IEQ) = X(IEQ) + ALPHA1*P(IEQ)/ALPHA2 
	R(IEQ) = R(IEQ) - ALPHA1*V(IEQ)/ALPHA2
C	{Z}    = [C-1]*{R i+1}
	Z(IEQ) = R(IEQ)*D(IEQ)
C	{ALP3} = {Z i+1}T*{R i+1}
	ALPHA3 = ALPHA3 + Z(IEQ)*R(IEQ)
	TEST   = TEST   + R(IEQ)*R(IEQ)
	ENDDO
	BETA = ALPHA3/ALPHA1


	TEST = SQRT(TEST)

	IF(TEST.LT.TOL) GOTO 3500
	IF(ITR.GT.4000) THEN
	WRITE(*,7000) 
	STOP
	ENDIF

	DO IEQ = 1,NEQ
	P(IEQ) = Z(IEQ) + BETA*P(IEQ)
	ENDDO
	
	GOTO 2000

C	----------------------------------------
 3500	CONTINUE

	B(1:NEQ) = X(1:NEQ)

	CALL CPU_TIME (TIME2)
C	WRITE(*,6000) NEQ,ITR,TIME2-TIME1


	DEALLOCATE(X,P,R,V,D,Z)


 6000	FORMAT (X,'SOLVING.',I6,' EQNS. WITH NO. ITR.',I5,
	1		'  TIME TO PERFORM. .',E15.6)
 7000	FORMAT (X,'**** ITERATION PROCESS NOT CONVERGE ****')

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================	
	SUBROUTINE CONJUGI(MAXA,A,B,NEQ,IND,TYP1,TYP2,DD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C	=============================================================
C	CONJUGATE GRADIENT ITERATION WITH IMCOMPLETE CHOL PRECONDITIONER
	DIMENSION MAXA(1),MAXC(1),A(1),B(1),DD(1),DIA(NEQ)
C	=============================================================
	ALLOCATABLE X(:),P(:),R(:),V(:)
	ALLOCATABLE D(:),Z(:)


	INDPD = 2

	IF(IND.EQ.1) THEN
	CALL MDMOVE(TYP1,TYP2)
	CALL IMC_LDL (MAXA,A,DD,B,NEQ,IND,INDPD,TYP1,'ICCG')
	RETURN
	ENDIF

	ALLOCATE(X(NEQ),P(NEQ),R(NEQ),V(NEQ),D(NEQ),Z(NEQ))

	CALL CPU_TIME (TIME1)

C	------------------------------------
	DO IEQ = 1,NEQ
	IF(B(IEQ).NE.0.0D0) GOTO 100
	ENDDO
	DEALLOCATE(X,P,R,V,D,Z)
	RETURN
C	------------------------------------

 100	TOL = 1.0E-2
	X(1:NEQ) = 0.0D0
	P(1:NEQ) = 0.0D0
	R(1:NEQ) = 0.0D0

C	------------------------------------
C	ESTIMATE OF X  --- X = [C-1]*B
	X(1:NEQ) = B(1:NEQ)
	CALL IMC_LDL (MAXA,A,DD,X,NEQ,IND,INDPD,'ICCG','ICCG')
C	------------------------------------
C	{R} = {B} - [A]*{X}
	CALL MASULT(MAXA,A,X,P,TYP2,'STD')
	R(1:NEQ) = B(1:NEQ) - P(1:NEQ)
C	{P} = {Z} = [C-1]*{R}
	Z(1:NEQ) = R(1:NEQ)
	CALL IMC_LDL (MAXA,A,DD,Z,NEQ,IND,INDPD,'ICCG','ICCG')
	P(1:NEQ) = Z(1:NEQ)
C	------------------------------------


	ITR = 0
C	----------------------------------------
2000	CONTINUE
	ITR = ITR + 1

C	{V} = [A]*{P}
	CALL MASULT(MAXA,A,P,V,TYP2,'STD')

	ALPHA1 = 0.0D0
	ALPHA2 = 0.0D0
	DO IEQ = 1,NEQ
C	{ALP1} = {R}T*{R}
	ALPHA1 = ALPHA1 + Z(IEQ)*R(IEQ)
C	{ALP2} = {P}T*{R}
	ALPHA2 = ALPHA2 + P(IEQ)*V(IEQ)
	ENDDO

	DO IEQ = 1,NEQ
	X(IEQ) = X(IEQ) + ALPHA1*P(IEQ)/ALPHA2 
	R(IEQ) = R(IEQ) - ALPHA1*V(IEQ)/ALPHA2
	ENDDO

C	{Z} = [C-1]*{R i+1}
	Z(1:NEQ) = R(1:NEQ)
	CALL IMC_LDL (MAXA,A,DD,Z,NEQ,IND,INDPD,'ICCG','ICCG')

	TEST   = 0.0D0
	ALPHA3 = 0.0D0
	DO IEQ = 1,NEQ
C	{ALP3} = {Z i+1}T*{R i+1}
	ALPHA3 = ALPHA3 + Z(IEQ)*R(IEQ)
	TEST   = TEST   + R(IEQ)*R(IEQ)
	ENDDO
	BETA = ALPHA3/ALPHA1

	TEST = SQRT(TEST)

	IF(TEST.LT.TOL) GOTO 3500
	IF(ITR.GT.4000.OR.ITR.GT.NEQ) THEN
	WRITE(*,7000) 
	STOP
	ENDIF

	DO IEQ = 1,NEQ
	P(IEQ) = Z(IEQ) + BETA*P(IEQ)
	ENDDO
	
	GOTO 2000

C	----------------------------------------
3500	CONTINUE


	B(1:NEQ) = X(1:NEQ)

	CALL CPU_TIME (TIME2)
C	WRITE(*,6000) NEQ,ITR,TIME2-TIME1


	DEALLOCATE(X,P,R,V,D,Z)

 6000	FORMAT (X,'SOLVING.',I6,' EQNS. WITH NO. ITR.',I5,
	1		'  TIME TO PERFORM. .',E15.6)
 7000	FORMAT (X,'**** ITERATION PROCESS NOT CONVERGE ****')

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE IMC_LDL (MAXA,A,D,V,NEQ,IND,INDPD,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C     -------------------------------------------------------------------------
C	IND = 1
C	DOING IMCOMPLETE CHOLESKI FACTORIZATION OF SPARSE MATRIX
C     -------------------------------------------------------------------------
C	IND = 2
C	V = [C-1]*V WHERE C = IMCOMPLETE CHOLESKI FACTORIZATION OF SPARSE MATRIX
C     -------------------------------------------------------------------------
C
	DIMENSION MAXA(1),A(1),D(1),V(1)
	ALLOCATABLE NCLM(:),ICOP(:),IHIG(:),B(:),MAXC(:),MAXB(:)

	STOL  = 1.0E-8

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)


	ALLOCATE(NCLM(NBLOCK),ICOP(NBLOCK),IHIG(NBLOCK))

	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCLM(IBLO),1,IBLO,0)
	CALL MINTFIL('BLOC',ICOP(IBLO),2,IBLO,0)
	CALL MINTFIL('BLOC',IHIG(IBLO),3,IBLO,0)
	ENDDO

	STOL = 1.0E-6

      IF (IND-2) 100,800,800

100	FACD = 1.0D0
110	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)
	CALL MCALFIL(KTREC,TYP2)
	REWIND(KTREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)

	CALL MDMOVI('MAXC','TEM2')
	CALL MCALFIL(NBLCH,'TEM2')
	REWIND(NBLCH)


	ALLOCATE(B(MSTOR))
	
	CALL CPU_TIME (TIME1)
C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	NCOLB = 0

	DO 1000 IBLO=1,NBLOCK

	NCOL = NCLM(IBLO)
      NPRE = MAXA(NCOLB+1)-1
	
	NHIG = IHIG(IBLO)
	READ(KSREC) A(1:NHIG)

	ALLOCATE(MAXC(NHIG))
	READ(NFLCH) MAXC(1:NHIG)
C     --------------------------------------------------------
C     FACTORIZE STIFFNESS MATRIX A (A = L*D*L'T DECOMPOSISION)
C     --------------------------------------------------------
	IF (IBLO.EQ.ICOP(IBLO))  GOTO 400

C     ------------------------------------------------
C     FIND SUM (MCOLB) OF PRECEEDING UNCOUPLED COLUMNS
C     FIRST (JBF) AND LAST (JBL) COUPLED BLOCK
C     ------------------------------------------------
      JBF = ICOP(IBLO)-1
      MCOLB = 0
      IF (JBF)  400,150,120
 120  DO 130  JBLO=1,JBF
 130  MCOLB = MCOLB+NCLM(JBLO)
 150  MCOL  = NCOLB-MCOLB
      JBF = JBF+1
      JBL = IBLO-1
C
C     ----------------------------------------------
C     REDUCE BLOCK BY THE PRECEEDING COUPLING BLOCKS
C     LOOP OVER NUMBER OF COUPLING BLOCKS
C     READ BLOCK FROM RANDOM ACESS TAPE
C     ----------------------------------------------
	REWIND(KTREC)
	REWIND(NBLCH)
	DO II = 1,JBF-1
	NHIG = IHIG(II)
	READ(KTREC) !B(1:NHIG)
	READ(NBLCH) !MAXB(1:NHIG)
	ENDDO

      DO 390  JBLO=JBF,JBL
	NHIG = IHIG(JBLO)
	READ(KTREC) B(1:NHIG)
      MCOL = MCOL-NCLM(JBLO)
      MPRE = MAXA(MCOLB+1)-1

	ALLOCATE(MAXB(NHIG))
	READ(NBLCH) MAXB(1:NHIG)
C
C     ------------------------------------------
C     LOOP OVER NUMBER OF COLUMNS OF BLOCK TO BE
C     REDUCED,ADDRESSES OF COLUMN ELEMENTS
C     ------------------------------------------
      DO 295  ICOL=1,NCOL
      KDIA = MAXA(NCOLB+ICOL)-NPRE
      KLOW = KDIA+1
      KUPP = MAXA(NCOLB+ICOL+1)-1-NPRE
      KHEI = KUPP-KLOW-ICOL+1
      KMOD = KHEI-MCOL

	IF (KMOD.LE.0) GOTO 295
C     -------------------------------------------------
C     LOOP OVER NUMBER OF COUPLING COLUMNS (JCF TO JCL)
C     ADDRESSES OF ELEMENTS IN COUPLING COLUMN
C     -------------------------------------------------
      LCOP = 0
      JCF  = NCLM(JBLO)-KMOD+1
      JCL  = NCLM(JBLO)
      IF (JCF.GT.0)  GOTO 210
      LCOP = KMOD-NCLM(JBLO)
      JCF  = 1
 210  KTOP = KUPP-LCOP
C
      DO 290  JCOL=JCF,JCL
      LCOP = LCOP+1
      KTOP = KTOP-1
      LDIA = MAXA(MCOLB+JCOL)-MPRE
      LHEI = MAXA(MCOLB+JCOL+1)-1-LDIA-MPRE
      IF (LHEI)  290,290,230
 230  NL = MIN0(LCOP,LHEI)
      C = 0.0D0

      DO 260  IL=1,NL
	ICN = MAXB(LDIA+IL)
      DO 250  JL=1,NL
	JCN = MAXC(KTOP+JL)
	IF(JCN.EQ.ICN) C = C + B(LDIA+IL)*A(KTOP+JL)
	IF(JCN.EQ.ICN) GOTO 260
 250	CONTINUE
 260	CONTINUE


      A(KTOP) = A(KTOP)-C
 290  CONTINUE
 295  CONTINUE
C
	DEALLOCATE(MAXB)
 390  MCOLB = MCOLB+NCLM(JBLO)

C     -----------------------------------------
C     REDUCE BLOCK, LOOP OVER NUMBER OF COLUMNS
C     -----------------------------------------
 400  DO 690  ICOL=1,NCOL
      KDIA = MAXA(NCOLB+ICOL)-NPRE
	KLOW = KDIA+1
	KUPP = MAXA(NCOLB+ICOL+1)-1-NPRE
	KCOL = ICOL + NCOLB
	KHEI = KUPP-KLOW
	KMOD = MIN0(KHEI,ICOL-1)

C	SCALE DIAGONAL TERM TO PREVENT NEGATIVE PIVOT
	A(KDIA) = A(KDIA)*FACD

	IF (KMOD) 600,500,410
C     ---------------------------------------------
C     LOOP OVER NUMBER OF COUPLING COLUMNS JCF, JCL
C     ---------------------------------------------
 410  LCOP = 0
      JCF  = ICOL-KMOD
      JCL  = ICOL-1
      KTOP = KLOW+KMOD
      IF (ICOL-1.LT.KHEI) LCOP = KHEI-ICOL+1	
      DO 490 JCOL=JCF,JCL
      LCOP = LCOP+1
      KTOP = KTOP-1
      LDIA = MAXA(NCOLB+JCOL)-NPRE
      LHEI = MAXA(NCOLB+JCOL+1)-1-LDIA-NPRE
      IF (LHEI) 490,490,430
 430  NL = MIN0(LCOP,LHEI)
      C = 0.0D0

      DO 460  IL=1,NL
	ICN = MAXC(LDIA+IL)
      DO 450  JL=1,NL
	JCN = MAXC(KTOP+JL)
	IF(JCN.EQ.ICN) C = C + A(LDIA+IL)*A(KTOP+JL)
	IF(JCN.EQ.ICN) GOTO 460
 450	CONTINUE
 460	CONTINUE

      A(KTOP) = A(KTOP) - C
 490  CONTINUE

C     -------------------------------------------
C     FINAL COLUMN TERMS LIJ = GIJ/DII  AND
C     DIAGONAL TERMS     DJJ = KJJ - SUM(IRJ*GRJ)
C     -------------------------------------------
 500  IEQ = KCOL
      SUM = 0.0
      DO 590  KA=KLOW,KUPP
      IEQ = IEQ-1
	JEQ = MAXC(KA) !!!!
	C = A(KA)/D(JEQ)
      SUM = SUM + C*A(KA)
 590  A(KA) = C
      A(KDIA) = A(KDIA) - SUM
C     ----------------------------------------------------------
C     SET DIAGONAL TERMS IN D AND TEST WHETHER POSITIVE DEFINITE
C     ----------------------------------------------------------
 600  D(KCOL) = A(KDIA)


C     ----------------------------------------------------------
C	IF DIAGONAL IS NEGATIVE LET INCREASE THE FACTOR OF DIAGONAL MULTIPLIER AND RE-FACTORIZE AGAIN
	IF(D(KCOL).LE.0.0D0) THEN

	IF(FACD.GT.10.0D0) INDPD = 0		!!!
	IF(FACD.GT.10.0D0) GOTO 605			!!!

	DEALLOCATE(B,MAXC)	
	FACD = FACD + 0.15
	GOTO 110
	ENDIF
C     ----------------------------------------------------------


 605  PIV = A(KDIA)
      IF (PIV-STOL)       610,610,690
 610  IF (DABS(PIV)-STOL) 650,650,620
 620  IF (INDPD-1)        630,690,690
 630  CALL ERRORS (12,KCOL,PIV,' SOLUTION ')
 650  IF (INDPD-1)        630,630,660
 660  PIV = STOL
      IF (PIV.EQ.0.0) PIV = -1.0D-16
      D(KCOL) = PIV
      A(KDIA) = PIV
 690  CONTINUE
C     ------------------------------------------------
C     END OF LOOP
C     CALCULATE DETERMINANT OF STIFFNESS MATRIX (DETK)
C     ------------------------------------------------

	NHIG = IHIG(IBLO)
	WRITE(KTREC) A(1:NHIG)

	DEALLOCATE(MAXC)
C     ------------------------------------------------
1000  NCOLB = NCOLB + NCOL  !END BLOCK LOOP
C     ------------------------------------------------

	CALL CPU_TIME (TIME2)
C	WRITE(*,6000) TIME2-TIME1


	DEALLOCATE(NCLM,ICOP,IHIG)
	DEALLOCATE(B)

C     ------------------------------------------------
 795  RETURN
C     ------------------------------------------------



 800	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	NCOLB = 0

	DO 2000 IBLO=1,NBLOCK

	NCOL = NCLM(IBLO)
      NPRE = MAXA(NCOLB+1)-1

	NHIG = IHIG(IBLO)
	READ(KSREC) A(1:NHIG)

	ALLOCATE(MAXC(NHIG))
	READ(NFLCH) MAXC(1:NHIG)
C     --------------------------------
C     FORWARD REDUCTION OF LOAD VECTOR
C     --------------------------------

C     ------------------------------------------------
C     LOOP OVER NUMBER OF COLUMNS, VI = RI-SUM(IRI*VR)
C     ------------------------------------------------
	DO 880  ICOL=1,NCOL
      KLOW = MAXA(NCOLB+ICOL)+1-NPRE
      KUPP = MAXA(NCOLB+ICOL+1)-1-NPRE
      IF (KUPP-KLOW) 880,830,830
 830  IEQ = NCOLB+ICOL
      KV  = IEQ
      C   = 0.0D0
      DO 850  KA=KLOW,KUPP
      KV = KV-1
	KV = MAXC(KA)  !!!
 850  C = C + A(KA)*V(KV)
      V(IEQ) = V(IEQ) - C
 880  CONTINUE

	DEALLOCATE(MAXC)
C     ------------------------------------------------
2000  NCOLB = NCOLB + NCOL  !END BLOCK LOOP
C     ------------------------------------------------

C     -----------------------------
C     BACK SUBSTITUTION , V=[D]-1*V
C     -----------------------------
 890	DO 900  IEQ=1,NEQ
 900  V(IEQ) = V(IEQ)/D(IEQ)


 905	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	KBLO = NBLOCK

	DO 3000 IBLO=1,NBLOCK
	REWIND(KSREC)
	REWIND(NFLCH)

	NCOL  =  NCLM(KBLO)
	NCOLB = NCOLB - NCOL
      NPRE  = MAXA(NCOLB+1)-1

	DO II = 1,KBLO-1
	NHIG = IHIG(II)
	READ(KSREC) !A(1:NHIG)
	READ(NFLCH) !MAXC(1:NHIG)
	ENDDO
	NHIG = IHIG(KBLO)
	READ(KSREC) A(1:NHIG)

	ALLOCATE(MAXC(NHIG))
	READ(NFLCH) MAXC(1:NHIG)
C     ---------------------------
C     LOOP OVER NUMBER OF COLUMNS
C     ---------------------------
 910	KCOL = NCOL
      DO 980  ICOL=1,NCOL
      KLOW = MAXA(NCOLB+KCOL)+1-NPRE
      KUPP = MAXA(NCOLB+KCOL+1)-1-NPRE
      IF (KUPP-KLOW) 980,930,930
 930  IEQ = NCOLB+KCOL
      KV = IEQ
      DO 950  KA=KLOW,KUPP
      KV = KV-1
	KV = MAXC(KA)  !!!
 950  V(KV) = V(KV) - A(KA)*V(IEQ)
 980  KCOL = KCOL-1

	DEALLOCATE(MAXC)
	KBLO = KBLO-1
C     ------------------------------------------------
3000  CONTINUE  !END BLOCK LOOP
C     ------------------------------------------------

	DEALLOCATE(NCLM,ICOP,IHIG)


 6000	FORMAT (X,'IM. CHOL. FACTORIZ. TIME TO PERFORM. .',E15.6)


      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================

      SUBROUTINE HEAPSORR(N,RA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION RA(N)
	
	IF(N.LE.1) RETURN
	
	L=N/2+1
	IR = N

   10 CONTINUE 
      IF(L.GT.1)THEN
	   L=L-1
	   RRA = RA(L)
	ELSE 
	   RRA = RA(IR)
	   RA(IR)=RA(1)
	   IR = IR-1
	   IF(IR.EQ.1)THEN
	    RA(1)=RRA
	    RETURN
	   ENDIF
	ENDIF
	
      I=L
	J=L+L
   20	IF(J.LE.IR)THEN
	    IF(J.LT.IR)THEN
	      IF(RA(J).LT.RA(J+1))J=J+1
	    ENDIF
	    IF(RRA.LT.RA(J))THEN
	      RA(I)=RA(J)
	      I=J
	      J=J+J
	    ELSE
	      J=IR+1
	    ENDIF
	    GO TO 20
	ENDIF
	RA(I)=RRA
	GO TO 10
	
	END

C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE HEAPSORI(N,IA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION IA(N)
	
	IF(N.LE.1) RETURN
	
	L=N/2+1
	IR = N

   10 CONTINUE 
   
      IF(L.GT.1)THEN
	   L=L-1
	   IIA = IA(L)
	ELSE 
	   IIA = IA(IR)
	   IA(IR)=IA(1)
	   IR = IR-1
	   IF(IR.EQ.1)THEN
	    IA(1)=IIA
	    RETURN
	   ENDIF
	ENDIF
	
      I=L
	J=L+L
   20	IF(J.LE.IR)THEN
	    IF(J.LT.IR)THEN
	      IF(IA(J).LT.IA(J+1)) J=J+1
	    ENDIF
	    IF(IIA.LT.IA(J)) THEN
	      IA(I)=IA(J)
	      I=J
	      J=J+J
	    ELSE
	      J=IR+1
	    ENDIF
	    GO TO 20
	ENDIF
	IA(I)=IIA
	GO TO 10
	
	END

C	=============================================================
C	=============================================================
C	=============================================================

