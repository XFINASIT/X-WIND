C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ELLWOK(IEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	ELEMENT WORKING ARRAY INITIALIZATION


	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !


	DO 5000 IEL = 1,NELE

	CALL MDUWOK(IEG,IEL)

5000	CONTINUE



	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE MDUWOK(IEG,IEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	ELEMENT WORKING ARRAY INITIALIZATION


	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !

	SELECTCASE(ITYPE)

	CASE(5)
	CALL FMEWOK(IEG,IEL)


	ENDSELECT



	RETURN
	END




C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FMEWOK(IEG,IEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	ELEMENT WORKING ARRAY INITIALIZATION


	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !

	SELECTCASE(ISTYP)

	CASE(5)
	CALL FRMWOK(IEG,IEL)


	ENDSELECT



	RETURN
	END




C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FRMWOK(IEG,IEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	INITIALIZE WORKING ARRAY
	DIMENSION NPM(10),NPI(10)

	ALLOCATABLE IF3(:),AF4(:),AF5(:)
	ALLOCATABLE IAX(:),AAX(:)
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) ! MAXP
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) ! NPT
	CALL INTFILL('%GRP',MFIB ,6 ,IEG,0) !
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	CALL INTFILL('%GRP',NNM  ,8 ,IEG,0) !
	CALL INTFILL('%GRP',NNF  ,9 ,IEG,0) !
	CALL INTFILL('%GRP',NGR  ,10,IEG,0) !
	CALL INTFILL('%GRP',NGS  ,11,IEG,0) !
	CALL INTFILL('%GRP',NGT  ,12,IEG,0) !
C	GROUP FILE
	CALL INTFILL('%GRF',N3   ,3 ,IEG,0) !
	CALL INTFILL('%GRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('%GRF',N5   ,5 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	
	CALL INTFILL('%GRF',NFLI ,21,IEG,0) !
	CALL INTFILL('%GRF',NFLR ,22,IEG,0) !
	
	CALL INTFILL('%GRF',LASI ,23,IEG,0) !
	CALL INTFILL('%GRF',LASR ,24,IEG,0) !
	
C	------------------------------------
	ALLOCATE (IF3(N3),AF4(N4),AF5(N5))
	ALLOCATE (IAX(LASI),AAX(LASR))

      READ(NFLI,REC=IEL) IAX
      READ(NFLR,REC=IEL) AAX
      
C	------------------------------------  
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')
	
	ISEC = IF3(2)     !SECTION SET
	LSTY = IF3(5)     !SEC TYPE FLAG
	NPRT = IF3(6)     !NUMBER OF PART
C	------------------------------------ 

C     ----------------------------------------------------------
C     LOOP OVER GAUSS 
C     ----------------------------------------------------------
      DO 400 IGR = 1,NGR


	DO 300 IPRT = 1,NPRT

	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
	NPIB    = INT(AF4(2))   !NUM FIBER FOR THIS PART


	SELECTCASE(LSTY)

	CASE(0) !READY PROP
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   1,IGR ,IPRT]
	CALL REDWOK(ISEC,NPM,NPI,AAX,NDX5)

	CASE(1) !FIBER PROP
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   0,IGR ,IPRT]   !NFIB TO BE MODIFIED IN NEXT SUBRUTINE
	CALL FIBWOK(ISEC,NPM,NPI,NPIB,AAX,NDX5)

	ENDSELECT

300	CONTINUE


C     ----------------------------------------------------------
400	CONTINUE
C     ----------------------------------------------------------

      WRITE(NFLI,REC=IEL) IAX
      WRITE(NFLR,REC=IEL) AAX
      
	DEALLOCATE(IAX,AAX)		 
	DEALLOCATE(IF3,AF4,AF5)


	RETURN

	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE REDWOK(ISEC,NPM,NPI,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION NPM(1),NPI(1),ASEC(15),NSM(5),NSI(5),AAX(1)
	ALLOCATABLE AF5(:)
C	LTMOPT(7)  INITIALIZE WORKING ARRAY
C	A				1
C	OFFSET SC		2
C	OFFSET TC		3
C	QS            4
C	QT            5
C	IS		    6
C	IT		    7
C	IST		    8
C	J				9
C	IW		    10
C	MATERIAL FLAG	11	
C	CONC or STEEL	12	
C	SHEAR COEF S	13
C	SHEAR COEF T	14

	N5   = NPM(1)
	IGR  = NPI(3)
	IPRT = NPI(4)
	ALLOCATE(AF5(N5))


	CALL INTFILL('$FSC',NSEC ,1,1 ,0)
	CALL INTFILL('$FSC',MAXP ,1,4 ,0)
	CALL INTFILL('$FSC',MX   ,1,5 ,0)
	CALL INTFILL('$FSC',NGP  ,1,6 ,0)
	CALL INTFILL('$FSC',NFSCL,1,8 ,0)  !FILE ID FOR SECTION DATA
	CALL INTFILL('$FSC',MRX  ,1,9 ,0)  !MAX DATA PER RECORD
	CALL INTFILL('$FSC',MAXF ,1,10,0)  !MAX FIBER

	NSM(1:5) = [MRX,MAXF,NGP,MAXP,NSEC]
	NSI(1:5) = [MRX,1   ,IGR,IPRT,ISEC]
	CALL PUTREAL(NFSCL,ASEC,NSM,NSI,5,0)  

	AREA = ASEC(1)
	MATN = INT(ASEC(11))
	MATY = INT(ASEC(12))

	RHO  = 0.0D0
	IF(MATY.EQ.0) CALL RELFILL('%FNC',RHO,4,MATN,0)
	IF(MATY.EQ.1) CALL RELFILL('%MSL',RHO,9,MATN,0)
C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')     !CALLING DATA OF THIS PART

	AF5(1:14) = ASEC(1:14)				  !TRANSFER DATA HERE
	AF5(20  ) = RHO   				      !MASS DEN

	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')     !BACKUP  DATA OF THIS PART
C	----------------------------------------
	

	DEALLOCATE(AF5)

	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FIBWOK(ISEC,NPM,NPI,NPIB,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(7)  INITIALIZE WORKING ARRAY
	DIMENSION NPM(1),NPI(1),ASEC(15),NSM(5),NSI(5),AAX(1)
	ALLOCATABLE AF5(:)

	N5   = NPM(1)
	IGR  = NPI(3)
	IPRT = NPI(4)
	ALLOCATE(AF5(N5))

	CALL INTFILL('$FSC',NSEC ,1,1 ,0)
	CALL INTFILL('$FSC',MAXP ,1,4 ,0)
	CALL INTFILL('$FSC',MX   ,1,5 ,0)
	CALL INTFILL('$FSC',NGP  ,1,6 ,0)
	CALL INTFILL('$FSC',NFSCL,1,8 ,0)  !FILE ID FOR SECTION DATA
	CALL INTFILL('$FSC',MRX  ,1,9 ,0)  !MAX DATA PER RECORD
	CALL INTFILL('$FSC',MAXF ,1,10,0)  !MAX FIBER

C	----------------------------------------	
	DO 1000 IFIB = 1,NPIB
C	----------------------------------------

C	CALLING DATA FROM WORKING ARRAY HERE
	NPI(2) = IFIB
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART

C	----------------------------------------
C	CALLING DATA FROM THIS SECTION PART
	NSM(1:5) = [MRX,MAXF,NGP,MAXP,NSEC]
	NSI(1:5) = [MRX,IFIB,IGR,IPRT,ISEC]
	CALL PUTREAL(NFSCL,ASEC,NSM,NSI,5,0)  

C	----------------------------------------
C	UPDATE DATA TO WORKING ARRAY HERE
C	FIBER AREA		1
C	FIBER DISTANCE S	2
C	FIBER DISTANCE T	3
C	MATERIAL NO.		4
C	MATERIAL TYPE		5
C	WARPING COEF		6	EACH FIBER
C	TORSION COEF		7	EACH FIBER
C	SHEAR COEF S		8
C	SHEAR COEF T		9
C	SHRINKAGE FACTOR	27
	AF5(1:5) = ASEC(1:5)
	AF5(27)  = ASEC(6)
	AF5(6:7) = ASEC(7:8)
	AF5(8:9) = ASEC(10:11)
C	----------------------------------------
	RHO  = 0.0D0
	AREA = ASEC(1)
	MATN = INT(ASEC(4))
	MATY = INT(ASEC(5))
	IF(MATY.EQ.0) CALL RELFILL('%FNC',RHO,4,MATN,0)
	IF(MATY.EQ.1) CALL RELFILL('%MSL',RHO,9,MATN,0)

	AF5(20)  = RHO   				      !MASS DEN

C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')   


C	----------------------------------------
1000	CONTINUE
C	----------------------------------------
	

	DEALLOCATE(AF5)


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
