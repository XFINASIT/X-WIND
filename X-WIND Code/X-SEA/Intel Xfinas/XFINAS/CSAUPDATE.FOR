C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ELLUPD(IEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
C     FOR SOLID 24EAS ONLY -----------------       
      COMMON /CSAMM/ MM,MM1,MM2,NDIMC
C     --------------------------------------      
      
C	STRESS AND RIGIDITY & MATERIAL STATE
	DIMENSION NPM(10),NPI(10)
	ALLOCATABLE IF1(:),IF2(:),IF3(:),AF6(:)
	ALLOCATABLE COORD(:,:),EDIS(:)
	ALLOCATABLE IAX(:),AAX(:)


	CALL INTFILL('%NUB',NSN,1,1,0)
	CALL INTFILL('%NUB',NSF,1,2,0)
	CALL INTFILL('%NUB',NSC,1,4,0)

      CALL INTFILL('%SOL',LOPT1,1,1,0)
      CALL INTFILL('%SOL',LOPT2,1,2,0)

C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) !
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
	CALL INTFILL('%GRP',MFIB ,6 ,IEG,0) !
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	CALL INTFILL('%GRP',NNM  ,8 ,IEG,0) !
	CALL INTFILL('%GRP',NNF  ,9 ,IEG,0) !
	CALL INTFILL('%GRP',NGR  ,10,IEG,0) !
	CALL INTFILL('%GRP',NGS  ,11,IEG,0) !
	CALL INTFILL('%GRP',NGT  ,12,IEG,0) !
C	GROUP FILE
	CALL INTFILL('%GRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('%GRF',N2   ,2 ,IEG,0) !
	CALL INTFILL('%GRF',N3   ,3 ,IEG,0) !
	CALL INTFILL('%GRF',N6   ,6 ,IEG,0) !

	CALL INTFILL('%GRF',NDX1 ,11,IEG,0) !
	CALL INTFILL('%GRF',NDX2 ,12,IEG,0) !
	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX6 ,16,IEG,0) !
	
	CALL INTFILL('%GRF',NFLI ,21,IEG,0) !
	CALL INTFILL('%GRF',NFLR ,22,IEG,0) !
	
	CALL INTFILL('%GRF',LASI ,23,IEG,0) !
	CALL INTFILL('%GRF',LASR ,24,IEG,0) !

	ALLOCATE(IF1(N1),IF2(N2),IF3(N3),AF6(N6))
	ALLOCATE(COORD(3,NNM),EDIS(NEF))
	ALLOCATE(IAX(LASI),AAX(LASR))

C     FOR SOLID 24EAS ONLY -----------------	
      IF (ITYPE == 10) THEN
      IF (ISTYP ==  1) THEN	  !SOLICO
      MM    = 11
      MM1   = 24
      MM2   = 0
      NDIMC = 0 
      ELSEIF (ISTYP == 6) THEN 	  !SOLID EAS 24
      MM    = 24
      MM1   = 24
      MM2   = 0
      NDIMC = 0    
      ENDIF
      ENDIF 	
C     --------------------------------------  

C	------------------------------------
C	LOOP OVER ELEMENT
C	------------------------------------
	DO 5000 IEL = 1,NELE
	
      READ(NFLI,REC=NELE+IEL,ERR=5) IXIST !READ EXISTING FLAG TO LAST RECORD OF ELEMENT INTEGER DATA (INCREASE THE SPEED PURPOSE)
      GOTO 6
5     IXIST = 0
6     IF(IXIST.LE.0) GOTO 5000
      
      READ(NFLI,REC=IEL) IAX
      READ(NFLR,REC=IEL) AAX
	
C	------------------------------------
	IF(LOPT2.GT.0) THEN !UPDAT NON-MECHANICAL BEFOR MECHANICAL
C	--------------------------
	CALL MDUOPT2(IEG,IEL,IAX,AAX)
	ENDIF
C	------------------------------------

C	------------------------------------
	IF(LOPT1.GT.0) THEN
C	--------------------------
C	EXTRACT NODAL COORDINATE
	NPM(1) = N1
	NPI(1) = N1
	CALL CALLINTG(IF1,IAX,NDX1,NPM,NPI,1,'READ')
	DO INM = 1,NNM
	NOD = IF1(INM)
	COORD(1:3,INM) = 0.0D0
	DO ISC = 1,NSC
	CALL RELFILL('$XYZ',COORD(ISC,INM),ISC,NOD,0)  !GETTING HERE NODAL COORDINATE
	ENDDO
	ENDDO
C	--------------------------
C	--------------------------
C	UPDATE THE ELEMENT DISPLACEMENT DUE TO STRESSING CABLE
	IF(LOPT1.EQ.1) THEN	   !WHEN UPDATE STRESS
	CALL INTFILL('%COV',ICONV,1,10,0)                     !ITERATION CONVERGENCE FLAG
	IF(ITYPE.EQ.2) THEN
	IF(ISTYP.EQ.4.OR.ISTYP.EQ.7) THEN		                !CABLE PARA 
	NPM(1) = N2
	NPI(1) = N2
	CALL CALLINTG(IF2,IAX,NDX2,NPM,NPI,1,'READ')          !IF2(1:NEF) = ORIGINAL EQAUATION
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')          !CALLING TENDON DATA
	
	IF(ICONV.GT.0) THEN
	IXIST = IF3(4) !EXISTING FLAG  UPDATE WHEN CONVERGENCE ACHIEVE
	IF(IXIST.EQ.1) THEN 
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')          !INITIAL DISPLACEMENT                   
	IF3(4) = 2						                    !CHANGE IXIST TO 2 (MEAN THAT WE UPDATE THE EFFECT DURING JACKING ALREADY)
	EDIS(1:NEF) = 0.0D0
	DO IEF = 1,NEF
	IEQ = IF2(IEF)
	CALL RELFILL('$XDS',EDIS(IEF),IEQ,1,0)                !GETTING HERE NODAL DISPLACEMENT
	AF6(IEF) = EDIS(IEF)			                        !AF6(1:NEF) = INITAIL ERECTED DISP.
	ENDDO
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')          !UPDATE INITIAL DISPLACEMENT DUE TO JACKING
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'WRIT')          !CALLING TENDON DATA
	ENDIF
	ENDIF
	
	ENDIF
	ENDIF
	ENDIF
C	--------------------------
C	EXTRACT NODAL DISPLACEMENT	
	NPM(1) = N2
	NPI(1) = N2
	CALL CALLINTG(IF2,IAX,NDX2,NPM,NPI,1,'READ')        !IF2(1:NEF) = ORIGINAL EQAUATION
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')        !INITIAL ERECTION DISPLACEMENT
	EDIS(1:NEF) = 0.0D0
	DO IEF = 1,NEF
	IEQ = IF2(IEF)
	CALL RELFILL('$XDS',EDIS(IEF),IEQ,1,0)    !GETTING HERE NODAL DISPLACEMENT
	IF(ITYPE.NE.5) EDIS(IEF) = EDIS(IEF) - AF6(IEF)          !AF6(1:NEF) = INITAIL ERECTED DISP.  !FRAME ELEMENT WILL CORRECT THE ERECTED DISPLACEMENT IN THEIR SUBROUTINE (BECAUSE IT HAVE MORE THAN ONE SECTION PART)
	ENDDO
C	--------------------------
	CALL MDUOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)
	ENDIF
C	------------------------------------

      WRITE(NFLI,REC=IEL) IAX
      WRITE(NFLR,REC=IEL) AAX
C	------------------------------------
5000	CONTINUE
C	------------------------------------

	DEALLOCATE(IAX,AAX)
	DEALLOCATE(IF1,IF2,IF3,AF6)
	DEALLOCATE(COORD,EDIS)


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE MDUOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION COORD(3,1),EDIS(1),IAX(1),AAX(1)

	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) 

C	------------------------------------
	SELECTCASE(ITYPE)
C	------------------------------------

	CASE(2)  !TRUSS
	CALL LINOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)

	CASE(5)  !FRAME
	CALL FMEOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)

	CASE(9)  !SHELL
	CALL SHEOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)

	CASE(10) !SOLID
	CALL SOEOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)

C	------------------------------------
	ENDSELECT
C	------------------------------------


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE MDUOPT2(IEG,IEL,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION IAX(1),AAX(1)
      
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) 

C	------------------------------------
	SELECTCASE(ITYPE)
C	------------------------------------

	CASE(2)  !TRUSS
	CALL LINOPT2(IEG,IEL,IAX,AAX)
		
	CASE(5)  !FRAME
	CALL FMEOPT2(IEG,IEL,IAX,AAX)

	CASE(9)  !SHELL
	CALL SHEOPT2(IEG,IEL,IAX,AAX)

	CASE(10) !SOLID
	CALL SOEOPT2(IEG,IEL,IAX,AAX)


C	------------------------------------
	ENDSELECT
C	------------------------------------


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FMEOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION COORD(3,1),EDIS(1),IAX(1),AAX(1)


	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) 

	SELECTCASE(ISTYP)

	CASE(5)
	CALL FRMOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)


	ENDSELECT



	RETURN
	END




C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FMEOPT2(IEG,IEL,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION IAX(1),AAX(1)

	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) 


	SELECTCASE(ISTYP)

	CASE(5)
	CALL FRMOPT2(IEG,IEL,IAX,AAX)


	ENDSELECT



	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
      SUBROUTINE SOEOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION COORD(3,1),EDIS(1),IAX(1),AAX(1)


	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) 

	SELECTCASE(ISTYP)
            
	CASE(1) !SOLICO
	CALL SOCOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)
      
      CASE(6) !SOLID-EAS
	CALL SODOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)

	ENDSELECT

	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
      SUBROUTINE SOEOPT2(IEG,IEL,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION IAX(1),AAX(1)

	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) 

	SELECTCASE(ISTYP)
      
      CASE(1) !SOLICO
	CALL SOCOPT2(IEG,IEL,IAX,AAX)
            
	CASE(6) !SOLID-EAS
	CALL SODOPT2(IEG,IEL,IAX,AAX)

	ENDSELECT

	RETURN
	END






C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FRMOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(1)  WORKING ON MECHANICAL DATA
C	 0 = NO UPDATE STRESS & MATERIAL STATE
C	 1 = UPDATE STRESS & MATERIAL STATE
C	 2 = STIFFNESS ONLY

	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)

	DIMENSION COORD(3,1),EDIS(1),IAX(1),AAX(1)

	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION TT(14,14),TRANH(14,14),LREAS(14),TRANO(14,14),PROPO(6)
	DIMENSION BPG(10),BWG(10),BBX(14)
	DIMENSION EPS(7),NPM(10),NPI(10)

	DIMENSION BMATX(7,14)

	ALLOCATABLE IF3(:),AF4(:),AF5(:),AF6(:),AF7(:),AF10(:)
	ALLOCATABLE COORO(:,:),EDPRT(:,:)
	ALLOCATABLE REDIS(:,:),EDISL(:,:),EDISLP(:)
	ALLOCATABLE GPL(:),GPW(:)
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) ! MAXP
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
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
	CALL INTFILL('%GRF',N6   ,6 ,IEG,0) !
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !
	CALL INTFILL('%GRF',N10  ,10,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX6 ,16,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
	CALL INTFILL('%GRF',NDX10,20,IEG,0) !
C	----------------------------------------
	ALLOCATE(IF3(N3),AF4(N4),AF5(N5),AF6(N6),AF7(N7),AF10(N10))
	ALLOCATE(COORO(3,NNM),EDPRT(NEF,MPRT))
	ALLOCATE(REDIS(NEF,MPRT),EDISL(NEF,MPRT),EDISLP(14))


C	---------------------------------------------------------
	ALLOCATE(GPL(NGR),GPW(NGR))

	DO IGR = 1,NGR
	IF(IGR.EQ.1  ) GPL(IGR) = -1.0D0
	IF(IGR.EQ.NGR) GPL(IGR) =  1.0D0 
	IF(IGR.NE.1.AND.IGR.NE.NGR) GPL(IGR) =  GAUSP(IGR-1,NGR-2)
	IF(IGR.EQ.1  ) GPW(IGR) =  0.0D0
	IF(IGR.EQ.NGR) GPW(IGR) =  0.0D0 
	IF(IGR.NE.1.AND.IGR.NE.NGR) GPW(IGR) =  GAUSW(IGR-1,NGR-2)
	ENDDO
C	---------------------------------------------------------


C	CALLING LOPT1
      CALL INTFILL('%SOL',LOPT1,1,1,0)

C	---------------------------------------------------------	
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')
C	GET OFFSET VALUE
	DO I = 1,6
	NUM = I + NEF + NEF + NEF + NEF
	PROPO(I) = AF6(NUM)
	ENDDO

C	---------------------------------------------------------

	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ') !ELEMENT DATA

	ISEC = IF3(2)     !SECTION SET
	LSTY = IF3(5)     !SEC TYPE FLAG
	NPRT = IF3(6)     !NUMBER OF PART
		
C	GET RELEASE FLAG
	DO I = 1,14
	LREAS(I) = IF3(7+I)
	ENDDO

C	----------------------------------------------------------
C	CORRECTION OF ERECTED DISPLACEMENT FOR ALL PART		
	DO IPRT = 1,NPRT
	NPM(1:2)  = [N10,MPRT]
	NPI(1:2)  = [N10,IPRT]
	CALL CALLREAL(AF10,AAX,NDX10,NPM,NPI,2,'READ') 
	DO IEF = 1,NEF
	EDPRT(IEF,IPRT) = EDIS(IEF) - AF10(IEF)
	ENDDO
	ENDDO
C	---------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	CALL TRNOOF(TRANO,PROPO)
	CALL TRNOFT1(COORD,COORO,EDIS,PROPO,TRANO)

C	CALCULATE ORIGINAL S BASE VECTOR
	DO I = 1,3
	VR(I) = COORO(I,2) - COORO(I,1)
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)            !GET LENGTH OF ELEMENT HERE
	CALL FMVEVR(VR,VS,VT)
	CALL RELFILL('$RNG',RANG,1,ISEC,0)  !SECTION ROTATION
	CALL ROMBAC(VR,VS,VT,RANG)



C	---------------------------------------------------------
C	TRANSFER DIRECTLY ELEMENT TOTAL DISPLACEMENT (EDIS) TO ELEMENT PURE DEFORMATION (REDIS)
C	DO I = 1,14
C	REDIS(I)=EDIS(I)
C	ENDDO

C	DISPLACEMENT FOR ALL PART		
	DO IPRT = 1,NPRT
	DO IEF = 1,NEF
	REDIS(IEF,IPRT) = EDPRT(IEF,IPRT)
	ENDDO
	ENDDO
C	---------------------------------------------------------


C	OBTAIN END RELEASE TRANSFORMATION MATRIX
	CALL TRNHIG(TRANH,ELN,LREAS)

C	DETERMINE GAUSS POINT POSITION AND WEIGHT
	DO IGR = 1,NGR
	RI = GPL(IGR)  !GAUSP(IGR,NGR)
	RW = GPW(IGR)  !GAUSW(IGR,NGR)
	BPG(IGR) = 0.5*ELN*(1.0 + RI)
	BWG(IGR) = 0.5*ELN*RW
	ENDDO

C	INITILIZE ARRAY
	EPS(1:7) = 0.0				!STRAIN

C	---------------------------------------------------------
C	TRANSFORM DISP. TO LOCAL COORDINATE
C	CALL TT1A (VR,VS,VT,TT)
C	EDISL = MATMUL(TRANSPOSE(TT),REDIS)

C	DISPLACEMENT FOR ALL PART	
	CALL TT1A (VR,VS,VT,TT)	
	DO IPRT = 1,NPRT
	DO IEF = 1,NEF
	EDISL(IEF,IPRT) = 0.0D0
	DO JEF = 1,NEF
	EDISL(IEF,IPRT) = EDISL(IEF,IPRT) + TT(JEF,IEF)*REDIS(JEF,IPRT)
	ENDDO
	ENDDO
	ENDDO

C	---------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION {Ur} = [T]{U}
C	CALL TRNMUL(TRANH,EDISL,1)

C	DISPLACEMENT FOR ALL PART		
	DO IPRT = 1,NPRT
	CALL TRNMUL(TRANH,EDISL(1,IPRT),1)
	ENDDO
C	---------------------------------------------------------

C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ----------------------------------------------------------
      DO 400 IGR = 1,NGR

	NPM(1:2) = [N7,NGR]
	NPI(1:2) = [N7,IGR]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ') 
	N7_R = 1           !POINTER OF RIGIDITY
	N7_S = 1 + 20       !POINTER OF MECHANICAL STRESS
	N7_N = 1 + 20 + 8   !POINTER OF NON-MECHANICAL STRESS
	

	AF7(N7_R:N7_R+20-1) = 0.0D0	!INITIALIZE OF RIGIDITY
	AF7(N7_S:N7_S+8 -1) = 0.0D0	!INITIALIZE OF MECHANICAL STRESS


C     GAUSS LOCATION ALONG ELEMENT AXIS
	BXD = BPG(IGR)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
C	CALL BBXFRM(BBX,ELN,BXD)
C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
	CALL BBXFRMW(BMATX,ELN,BXD)

C     ------------------------------------------------------------
C     STRAIN TERMS MEMBRANE - BENDING - TORSION
C     ------------------------------------------------------------
C	EPS(1)= BBX(1)*EDISL(1)+BBX(8)*EDISL(8)
C	EPS(4)= BBX(3)*EDISL(3)+BBX(4)*EDISL(5)+BBX(10)*EDISL(10)+
C	+		BBX(11)*EDISL(12)
C	EPS(5)= BBX(2)*EDISL(2)+BBX(5)*EDISL(6)+BBX(9)*EDISL(9)+
C	+		BBX(12)*EDISL(13)
C	EPS(6)= BBX(6)*EDISL(4)+BBX(13)*EDISL(11)

C	EPS = MATMUL(BMATX,EDISL)
	
	DO 300 IPRT = 1,NPRT

C	DISPLACEMENT FOr EACH PART JAN09
	EDISLP(1:14) = EDISL(1:14,IPRT)
	EPS = MATMUL(BMATX,EDISLP)

	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
	IXIST   = INT(AF4(1))   !EXISTING FLAG
	NPIB    = INT(AF4(2))   !NUM FIBER FOR THIS PART

	IF (IXIST.EQ. 0) GOTO 300 !PART HAS NEVER EXISTED & CURRENTLY REMOVED

	SELECTCASE(LSTY)

	CASE(0) !READY PROP
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,1   ,IGR ,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
	CALL REDOPT1(AF7(N7_R),AF7(N7_S),AF5,EPS)

	CASE(1) !FIBER PROP
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   0,IGR ,IPRT]   !NFIB TO BE MODIFIED IN NEXT SUBRUTINE
	CALL FIBOPT1(NPM,NPI,NPIB,AF7(N7_R),
	1			 AF7(N7_S),EPS,LOPT1,AAX,NDX5)

	ENDSELECT


300	CONTINUE

	NPM(1:2) = [N7,NGR]
	NPI(1:2) = [N7,IGR]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT')   !BACK UP HERE RIGIDITY AND MCHA STRESS
C     ----------------------------------------------------------
400	CONTINUE
C     ----------------------------------------------------------


500	DEALLOCATE(IF3,AF4,AF5,AF6,AF7,AF10)
	DEALLOCATE(COORO,EDPRT)
	DEALLOCATE(REDIS,EDISL,EDISLP)
	DEALLOCATE(GPL,GPW)



	RETURN

	END





C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SODOPT1(IEG,IEL,COORD,EDIS,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(1)  WORKING ON MECHANICAL DATA
C	 0 = NO UPDATE STRESS & MATERIAL STATE
C	 1 = UPDATE STRESS & MATERIAL STATE
C	 2 = STIFFNESS ONLY


	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)		
C     SOLID ELEMENT -------------------------------------------------		 
C      COMMON /MMENH/ MM,MM1,MM2,NDIMC 
      COMMON /CSAMM/ MM,MM1,MM2,NDIMC 
          
C	BASIC MATRIX FOR COMPUTATION	
	DIMENSION COORD(3,1),EDIS(1),IAX(1),AAX(1)
	DIMENSION EPS(6),NPM(10),NPI(10)  
	
	DIMENSION DP(6,6),H(21),HD(3,21),XJ(3,3),XJI(3,3),B(63),DISD(9)
	DIMENSION BM(6,24),DE(6,6),DEP(6,6),STRAIN(6)   	
C     EAS METHOD -----------------------
	DIMENSION GE(6,MM)
	DIMENSION T0(6,6),TT(6,6),TT0(6,6)   
	DIMENSION XJO(3,3),TM(6,MM)
	DIMENSION EAS(6),ALPHA(MM),RE1(24)     
C    -------------------------------------------------------------	
	ALLOCATABLE IF3(:),AF4(:),AF7(:)
	ALLOCATABLE REDIS(:),NODEX(:)
C     ---------------------------------	
	  	
	
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) !MAXP
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
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
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
C	----------------------------------------
	ALLOCATE(IF3(N3),AF4(N4),AF7(N7))
	ALLOCATE(REDIS(NEF))
	NDEX = NNM - 8 + 1
	ALLOCATE(NODEX(NDEX))

C	CALLING SOLUTION CONTROL (LOPT1 = 1, UPDATE MATERIAL STATE)
      CALL INTFILL('%SOL',LOPT1,1,1,0)
      
C	IF3-R(ELEMENT DATA)---------------------------------------
	NPM(1) = N3   
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ') !ELEMENT DATA
	ISEC = IF3(2)     !SECTION SET
	NPRT = IF3(6)     !NUMBER OF PART
C	----------------------------------------------------------

C     24 DOFS  FOR EACH ELEMENT
      DO I = 1,24
	REDIS(I) = EDIS(I)
	ENDDO

C	INITILIZE ARRAY
      DO I = 1,6
	EPS(I)=0.0D0                  !STRAIN
	ENDDO
	
C     ----------------------------------------------------------    
C     EAS METHOD
C     ----------------------------------------------------------   
C     COMPUTE JACOBIAN MATRIX AT THE NATURAL CENTER (R=0,S=0,T=0)
      CALL SHAP3D_S (0.0D0,0.0D0,0.0D0,H,HD,NODEX,NNM)
	CALL JACO3D_S (COORD,HD,XJO,XJI,DETO,IEL,NNM)	
	
C     COMPUTE T MATRIX (FOR EAS METHOD)	
      CALL MATRIXT_S (XJO,T0)
      CALL INVMATRIX (T0,TT,6)      
      TT0 = TRANSPOSE(TT)
      
C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ----------------------------------------------------------
C     SETTING INDEX FOR GAUSS POINT
	IPT = 0
      DO 400 IGR = 1,NGR
      DO 400 IGS = 1,NGS
      DO 400 IGT = 1,NGT
	RI = GAUSP(IGR,NGR)
	SI = GAUSP(IGS,NGS)
	TI = GAUSP(IGT,NGT)
	WT = GAUSW(IGR,NGR)*GAUSW(IGS,NGS)*GAUSW(IGT,NGT)
	IPT = IPT + 1
	
C     POINTER FOR RIGIDITY AND MECHANICAL STRESS
C	AF7 ------------------------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  
C	----------------------------------------------------------
	N7_R = 1            !POINTER OF RIGIDITY
	N7_S = 1 + 12       !POINTER OF MECHANICAL STRESS (6 STRESSES)
	N7_N = 1 + 12 + 6   !POINTER OF NON-MECHANICAL STRESS (6 STRESSES)
		
	AF7(N7_R:N7_R+12-1) = 0.0D0	!INITIALIZE OF RIGIDITY
	AF7(N7_S:N7_S+6-1)  = 0.0D0	!INITIALIZE OF MECHANICAL STRESS
	
C     -----------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C     -----------------------------------------------------      
      CALL SHAP3D_S (RI,SI,TI,H,HD,NODEX,NNM)
      CALL JACO3D_S (COORD,HD,XJ,XJI,DET,IEL,NNM)
      DVOL = WT*DET
C     STRAIN-DISPLACEMENT MATRIX      
      CALL BMATSLD_S (HD,XJI,B,BM,NNM)
C     MATERIAL PROPERTIES (LINEAR)       
	!CALL DMATSLD_S (DE)
	
C     -----------------------------------------------------     
C     FOR ENHANCED STRAIN METHOD AND COMPUTE MATRIX [GE]
C     -----------------------------------------------------
	!CALL SOLIDGE_S (RI,SI,TI,DETO,DET,TT0,MM,GE)	
	CALL MTEAS7(RI,SI,TI,TT0,DETO,DET,GE,MM) 
	
C     -----------------------------------------------------
C     COMPUTE DISPLACEMENT GRADIENT (DISD)
C     -----------------------------------------------------
	CALL CLEARA (DISD,9)
      
      DO IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDIS(IEF)   ! H,x *u
      DISD(2) = DISD(2) + B(JEF)*EDIS(JEF)   ! H,y *v
      DISD(3) = DISD(3) + B(KEF)*EDIS(KEF)   ! H,z *w
      DISD(4) = DISD(4) + B(JEF)*EDIS(IEF)   ! H,y *u
      DISD(5) = DISD(5) + B(KEF)*EDIS(IEF)   ! H,z *u
      DISD(6) = DISD(6) + B(IEF)*EDIS(JEF)   ! H,x *v 
      DISD(7) = DISD(7) + B(KEF)*EDIS(JEF)   ! H,z *v
      DISD(8) = DISD(8) + B(IEF)*EDIS(KEF)   ! H,x *w
      DISD(9) = DISD(9) + B(JEF)*EDIS(KEF)   ! H,y *w
	END DO

C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      EPS(1) = DISD(1)
      EPS(2) = DISD(2)
      EPS(3) = DISD(3)
      EPS(4) = DISD(4) + DISD(6)
      EPS(5) = DISD(5) + DISD(8)
      EPS(6) = DISD(7) + DISD(9)      
      
C     -------------------------------
C     FOR EAS TERM  EAS=[GE]*{ALPHA}
C     -------------------------------
	!EAS = MATMUL(GE,ALPHA)

C	-------------------------------
C	COMPUTE TOTAL COMPATIBLE STRAIN
C	------------------------------- 
	!DO I=1,6
      !EPS(I)=EPS(I) - EAS(I)      
	!END DO      
      
	DO 300 IPRT = 1,NPRT   
C	AF4 ------------------------------------------------------
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
C     ----------------------------------------------------------	
	IXIST   = INT(AF4(1))   !EXISTING FLAG
	
C	AF5&AF7 --------------------------------------------------
	IF (IXIST.EQ. 0) GOTO 300 !PART HAS NEVER EXISTED & CURRENTLY REMOVED
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]   !MFIB = 1
	NPI(1:4) = [N5,   0,IPT ,IPRT]   !NFIB TO BE MODIFIED IN NEXT SUBRUTINE
C     ----------------------------------------------------------	
	CALL GUSOPT1(NPM,NPI,ISEC,AF7(N7_R),
	1			 AF7(N7_S),EPS,LOPT1,AAX,NDX5)
300	CONTINUE



C	--------------------------------
C	STORE HERE JAN09 LOADCOMB
C	--------------------------------
	IC = 0
	CALL INTFILL('@COM',IC,1,1,0)
	IF(IC.GT.0) THEN

C	PREVIOUS
	MS1_P = N7 - 6*11 + 6*(IC-1) + 1
	MS2_P = N7 - 6*11 + 6*(IC-1) + 6
C	CURRENT
	MS1_C = N7 - 6*11 + 6*(11-1) + 1
	MS2_C = N7 - 6*11 + 6*(11-1) + 6
      
      IF(IC.EQ.8) THEN !FOR CREEP&SHRINKAGE PRIMARY  STORE DIRECTLY TOTAL NON-MECHANICAL STRESS RESULTANT
	AF7(MS1_P:MS2_P) = AF7(N7_S:N7_S+6 -1)
      ELSE
	AF7(MS1_P:MS2_P) = AF7(MS1_P:MS2_P) + AF7(N7_S:N7_S+6 -1) -
	1				   AF7(MS1_C:MS2_C)			!BACK UP THE STRESS RESULTANT HERE FOR LOAD COMBINATION
      ENDIF
      
      
	ELSEIF(IC.EQ.-1) THEN

C	CURRENT
	MS1_C = N7 - 6*11 + 6*(11-1) + 1
	MS2_C = N7 - 6*11 + 6*(11-1) + 6
	AF7(MS1_C:MS2_C) = AF7(N7_S:N7_S+6 -1)

      
	ENDIF
C	--------------------------------

C     BACK UP HERE RIGIDITY AND MCHA STRESS
C	AF7 ------------------------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT')   
C     ----------------------------------------------------------
400	CONTINUE
C     ----------------------------------------------------------

	DEALLOCATE(IF3,AF4,AF7)
	DEALLOCATE(REDIS)
	DEALLOCATE(NODEX)


	RETURN

	END
	
	
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SOCOPT1(IEG,IEL,COORDI,EDIS,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(1)  WORKING ON MECHANICAL DATA
C	 0 = NO UPDATE STRESS & MATERIAL STATE
C	 1 = UPDATE STRESS & MATERIAL STATE
C	 2 = STIFFNESS ONLY


	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)		
C     SOLID ELEMENT -------------------------------------------------		 
C      COMMON /MMENH/ MM,MM1,MM2,NDIMC 
      COMMON /CSAMM/ MM,MM1,MM2,NDIMC 
          
C	BASIC MATRIX FOR COMPUTATION	
	DIMENSION COORD(3,8),EDIS(1),IAX(1),AAX(1)
	DIMENSION EPS(6),NPM(10),NPI(10) 
	
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION VRN(3),VSN(3),VTN(3)	 
	DIMENSION COORDI(3,8)	
	DIMENSION XY(3,8),BB(6,24)   
	DIMENSION EDISO(24),EDISI(24),RE(24)
	DIMENSION EDISLO(24)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
C	DIMENSION FIN(NEF)
	DIMENSION REP(24),SEP1(24,24),SEP2(24,24),EDISP(24),MPNOD(8),MPDIS(24),MPGUS(8)
C     REMOVE RIGID BODY MOTION
      DIMENSION COORDO(3,8),REDIS(24)
      		
	
	DIMENSION DP(6,6),H(21),HD(3,21),XJ(3,3),XJI(3,3),B(63),DISD(9)
	DIMENSION BM(6,24),DE(6,6),DEP(6,6),STRAIN(6)   	
C     EAS METHOD -----------------------
	DIMENSION GE(6,MM)
	DIMENSION T0(6,6),TT(6,6),TT0(6,6)   
	DIMENSION XJO(3,3),TM(6,MM)
	DIMENSION EAS(6),ALPHA(MM),RE1(24)     
C    -------------------------------------------------------------	
	ALLOCATABLE IF3(:),AF4(:),AF7(:)
	ALLOCATABLE NODEX(:)
C     ---------------------------------	
	  	
	
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) !MAXP
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
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
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
C	----------------------------------------
	ALLOCATE(IF3(N3),AF4(N4),AF7(N7))
	NDEX = NNM - 8 + 1
	ALLOCATE(NODEX(NDEX))

C	CALLING SOLUTION CONTROL (LOPT1 = 1, UPDATE MATERIAL STATE)
      CALL INTFILL('%SOL',LOPT1,1,1,0)
      
C	IF3-R(ELEMENT DATA)---------------------------------------
	NPM(1) = N3   
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ') !ELEMENT DATA
	ISEC = IF3(2)     !SECTION SET
	NPRT = IF3(6)     !NUMBER OF PART
C	----------------------------------------------------------

C	DETECTED THE THIN DIRECTION
C	----------------------------------------------------------
	XYZ1 = SQRT( (COORDI(1,2)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,2)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,2)-COORDI(3,1))**2.0 )
	XYZ2 = SQRT( (COORDI(1,4)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,4)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,4)-COORDI(3,1))**2.0 )
	XYZ3 = SQRT( (COORDI(1,5)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,5)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,5)-COORDI(3,1))**2.0 )

	NMX = 0
	IF(XYZ1.LT.XYZ2.AND.XYZ1.LT.XYZ3) NMX = 1
	IF(XYZ2.LT.XYZ1.AND.XYZ2.LT.XYZ3) NMX = 2
	IF(XYZ3.LT.XYZ1.AND.XYZ3.LT.XYZ2) NMX = 3

	IF(NMX.EQ.0) THEN
	IF(XYZ1.EQ.XYZ2) NMX = 1
	IF(XYZ1.EQ.XYZ3) NMX = 1
	IF(XYZ2.EQ.XYZ3) NMX = 2
	ENDIF

	CALL SLCOFACE(NMX,MPNOD,MPDIS,MPGUS)
	
C     FLIPPING COORDINATE
      DO I = 1,NNM
      NOD = MPNOD(I)
      COORD(1:3,I) = COORDI(1:3,NOD)
      ENDDO
      
C     FLIPPING DISPLACEMENT
      EDISP(1:24) = EDIS(1:24)
      DO I = 1,24
      IEF = MPDIS(I)
      EDISO(IEF) = EDISP(I)
      ENDDO
C     ------------------------------------------------------------
C	FIND THE LOCAL VECTORS,LOCAL COORDINATES AND DISPLACEMENTS
C     -----------------------------------------------------------   
C     COMPUTE JACOBIAN MATRIX AT THE NATURAL CENTER (R=0,S=0,T=0)   
      CALL SHAP3D (0.0D0,0.0D0,0.0D0,H,HD,NODEX,NNM)
C     CALL SHAP3D_S (0.0D0,0.0D0,0.0D0,H,HD,NODEX,NNM)      
  	CALL SOLIRST(COORD,HD,VR,VS,VT)   ! LOCAL VECTORS
  	
C     REMOVE RIGID BODY MOTION
  	REDIS(1:24) = EDISO(1:24)
      IF(NLOPT.EQ.3) THEN
          K = 0
          DO I = 1,NNM
              DO J = 1,3
                  K = K + 1
                  COORDO(J,I) = COORD(J,I) - EDISO(K)
              ENDDO
          ENDDO
          CALL SOMDSP(COORD,COORDO,EDISO,REDIS,H,HD,VR,VS,VT,NNM)
      ENDIF
C     REDIS = DISPLACEMENT REMOVING RIGID BOBY      
      CALL SCOOR(COORD,REDIS,XY,EDISLO,VR,VS,VT) !REDIS INSTEAD OF EDIS HERE 
C     ----------------------------------------------------------        
C	INITILIZE ARRAY
      DO I = 1,6
	EPS(I)=0.0D0                  !STRAIN
	ENDDO      
      
C     --------------------------------------------
C     EAS TRANSFORMATION MATRIX
C     --------------------------------------------	
	CALL JACO3D (XY,HD,XJO,XJI,DETO,IEL,NNM)	
C     CALL JACO3D_S (COORD,HD,XJO,XJI,DETO,IEL,NNM)		
	CALL TNEAS6(XJO,TT0)                       

C     24 DOFS  FOR EACH ELEMENT
      !DO I = 1,24
	!REDIS(I) = EDIS(I)
	!ENDDO	
C     ----------------------------------------------------------    
C     COMPUTE T MATRIX (FOR EAS METHOD)	
      !CALL MATRIXT_S (XJO,T0)
      !CALL INVMATRIX (T0,TT,6)      
      !TT0 = TRANSPOSE(TT)
      
C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ----------------------------------------------------------
C     SETTING INDEX FOR GAUSS POINT
	IPT = 0
      DO 400 IGR = 1,NGR
      DO 400 IGS = 1,NGS
      DO 400 IGT = 1,NGT
	RI = GAUSP(IGR,NGR)
	SI = GAUSP(IGS,NGS)
	TI = GAUSP(IGT,NGT)
	WT = GAUSW(IGR,NGR)*GAUSW(IGS,NGS)*GAUSW(IGT,NGT)
	IPT = IPT + 1
	
C     POINTER FOR RIGIDITY AND MECHANICAL STRESS
C	AF7 ------------------------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  
C	----------------------------------------------------------
	N7_R = 1            !POINTER OF RIGIDITY
	N7_S = 1 + 12       !POINTER OF MECHANICAL STRESS (6 STRESSES)
	N7_N = 1 + 12 + 6   !POINTER OF NON-MECHANICAL STRESS (6 STRESSES)
	
	N7_V = 1 + 12 + 6 + 6  !POINTER OF LOCAL VECTORS (6 STRESSES)
		
	AF7(N7_R:N7_R+12-1) = 0.0D0	!INITIALIZE OF RIGIDITY
	AF7(N7_S:N7_S+6-1)  = 0.0D0	!INITIALIZE OF MECHANICAL STRESS
	
C     -----------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C     -----------------------------------------------------      
      !CALL SHAP3D_S (RI,SI,TI,H,HD,NODEX,NNM)
      CALL SHAP3D (RI,SI,TI,H,HD,NODEX,NNM)
      !CALL JACO3D_S (COORD,HD,XJ,XJI,DET,IEL,NNM)
      CALL JACO3D (XY,HD,XJ,XJI,DET,IEL,NNM)
      DVOL = WT*DET
C     STRAIN-DISPLACEMENT MATRIX      
      !CALL BMATSLD_S (HD,XJI,B,BM,NNM)
50    CALL SOBMAT (HD,XJI,B,NNM)      
C     MATERIAL PROPERTIES (LINEAR)       
	!CALL DMATSLD_S (DE)
C     ASSUMED STRAIN-DISPLACEMENT MATRIX (BB)
      CALL SOBMATANS(XY,BM,B,RI,SI,TI,NNM,IEL,NODEX)   ! NEW ANS METHOD		
	
C     -----------------------------------------------------     
C     FOR ENHANCED STRAIN METHOD AND COMPUTE MATRIX [GE]
C     -----------------------------------------------------
	!CALL SOLIDGE_S (RI,SI,TI,DETO,DET,TT0,MM,GE)	
	CALL MTEAS7(RI,SI,TI,TT0,DETO,DET,GE,MM)
	
C     -----------------------------------------------------
C     COMPUTE DISPLACEMENT GRADIENT (DISD)
C     -----------------------------------------------------
	CALL CLEARA (DISD,9)
      
      DO IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDISLO(IEF)   ! H,x *u
      DISD(2) = DISD(2) + B(JEF)*EDISLO(JEF)   ! H,y *v
      DISD(3) = DISD(3) + B(KEF)*EDISLO(KEF)   ! H,z *w
      DISD(4) = DISD(4) + B(JEF)*EDISLO(IEF)   ! H,y *u
      DISD(5) = DISD(5) + B(KEF)*EDISLO(IEF)   ! H,z *u
      DISD(6) = DISD(6) + B(IEF)*EDISLO(JEF)   ! H,x *v 
      DISD(7) = DISD(7) + B(KEF)*EDISLO(JEF)   ! H,z *v
      DISD(8) = DISD(8) + B(IEF)*EDISLO(KEF)   ! H,x *w
      DISD(9) = DISD(9) + B(JEF)*EDISLO(KEF)   ! H,y *w
	END DO

C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      EPS(1) = DISD(1)
      EPS(2) = DISD(2)
      EPS(3) = DISD(3)
      EPS(4) = DISD(4) + DISD(6)
      !EPS(5) = DISD(5) + DISD(8)
      !EPS(6) = DISD(7) + DISD(9) 
      EPS(5)=0.0D0
	EPS(6)=0.0D0
	DO I=1,24     
      EPS(5)=EPS(5)+BM(5,I)*EDISLO(I)
	EPS(6)=EPS(6)+BM(6,I)*EDISLO(I)
      ENDDO     
      
C     -------------------------------
C     FOR EAS TERM  EAS=[GE]*{ALPHA}
C     -------------------------------
	!EAS = MATMUL(GE,ALPHA)

C	-------------------------------
C	COMPUTE TOTAL COMPATIBLE STRAIN
C	------------------------------- 
	!DO I=1,6
      !EPS(I)=EPS(I) - EAS(I)      
	!END DO      
      
	DO 300 IPRT = 1,NPRT   
C	AF4 ------------------------------------------------------
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
C     ----------------------------------------------------------	
	IXIST   = INT(AF4(1))   !EXISTING FLAG
	
C	AF5&AF7 --------------------------------------------------
	IF (IXIST.EQ. 0) GOTO 300 !PART HAS NEVER EXISTED & CURRENTLY REMOVED
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]   !MFIB = 1
	NPI(1:4) = [N5,   0,IPT ,IPRT]   !NFIB TO BE MODIFIED IN NEXT SUBRUTINE
C     ----------------------------------------------------------	
	CALL GUCOPT1(NPM,NPI,ISEC,AF7(N7_R),
	1			 AF7(N7_S),EPS,LOPT1,AAX,NDX5)
300	CONTINUE

C     ----------------------------------------------------------     
C     BACK UP LOCAL VECTORS
	AF7(N7_V+0:N7_V+2) = VR(1:3) 
	AF7(N7_V+3:N7_V+5) = VS(1:3) 
	AF7(N7_V+6:N7_V+8) = VT(1:3) 	



C	--------------------------------
C	STORE HERE JAN09 LOADCOMB
C	--------------------------------
	IC = 0
	CALL INTFILL('@COM',IC,1,1,0)
	IF(IC.GT.0) THEN

C	PREVIOUS
	MS1_P = N7 - 6*11 + 6*(IC-1) + 1
	MS2_P = N7 - 6*11 + 6*(IC-1) + 6
C	CURRENT
	MS1_C = N7 - 6*11 + 6*(11-1) + 1
	MS2_C = N7 - 6*11 + 6*(11-1) + 6
      
      IF(IC.EQ.8) THEN !FOR CREEP&SHRINKAGE PRIMARY  STORE DIRECTLY TOTAL NON-MECHANICAL STRESS RESULTANT
	AF7(MS1_P:MS2_P) = AF7(N7_S:N7_S+6 -1)
      ELSE
	AF7(MS1_P:MS2_P) = AF7(MS1_P:MS2_P) + AF7(N7_S:N7_S+6 -1) -
	1				   AF7(MS1_C:MS2_C)			!BACK UP THE STRESS RESULTANT HERE FOR LOAD COMBINATION
      ENDIF
      
      
	ELSEIF(IC.EQ.-1) THEN

C	CURRENT
	MS1_C = N7 - 6*11 + 6*(11-1) + 1
	MS2_C = N7 - 6*11 + 6*(11-1) + 6
	AF7(MS1_C:MS2_C) = AF7(N7_S:N7_S+6 -1)

      
	ENDIF
C	--------------------------------


C     BACK UP HERE RIGIDITY AND MCHA STRESS
C	AF7 ------------------------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT')   
C     ----------------------------------------------------------
400	CONTINUE
C     ----------------------------------------------------------

	DEALLOCATE(IF3,AF4,AF7)	
	DEALLOCATE(NODEX)


	RETURN

	END	





C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE REDOPT1(PA,SA,AF5,EPS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION PA(1),AF5(1),SA(1),EPS(8)

C	LTMOPT(1)  WORKING ON MECHANICAL DATA
C	 0 = NO UPDATE STRESS & MATERIAL STATE
C	 1 = UPDATE STRESS & MATERIAL STATE

C	----------------------------------------	
C	N = 1  A     
C	N = 2  SC     
C	N = 3  TC     
C	N = 4  QS     
C	N = 5  QT     
C	N = 6  IS     
C	N = 7  IT     
C	N = 8  IST     
C	N = 9  J     
C	N = 10 IW     
C	N = 11 MATERIAL FLAG     
C	N = 12 CONC or STEEL   

C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE


	Ao    = AF5(1 )
	SCo   = AF5(2 )
	TCo   = AF5(3 )
	QSo   = AF5(4 )
	QTo   = AF5(5 )
	SISo  = AF5(6 )
	SITo  = AF5(7 )
	SISTo = AF5(8 )
	SJo   = AF5(9 )
	SJW   = AF5(10)
	MATN  = INT(AF5(11))
	MATY  = INT(AF5(12))
	FSHS  = AF5(13)
	FSHT  = AF5(14)

	EMOD  = AF5(17)  
	POISN = AF5(18)
	GG    = AF5(19)
	RHO   = AF5(20)


	A    = Ao 
	QS   = QSo   + TCo*Ao
	QT   = QTo   + SCo*Ao
	SIT  = SITo  + 2.0*SCo*QTo + SCo*SCo*Ao
	SIS  = SISo  + 2.0*TCo*QSo + TCo*TCo*Ao
	SIST = SISTo + SCo*QSo + TCo*QTo + SCo*TCo*Ao
	SJ   = SIS + SIT


	PA(1) = PA(1) + A*EMOD
	PA(2) = PA(2) + QS*EMOD
	PA(3) = PA(3) + QT*EMOD
	PA(4) = PA(4) + SIS*EMOD
	PA(5) = PA(5) + SIT*EMOD
	PA(6) = PA(6) + SIST*EMOD
	PA(7) = PA(7) + SJ*GG
	PA(8) = PA(8) + A*RHO

	PA(9)  = SJW*EMOD	!EIw
	PA(10) = 0.0D0		!EIsw
	PA(11) = 0.0D0		!EItw
	PA(12) = 0.0D0		!EIw

	PA(13) = Ao*FSHS*GG
	PA(14) = Ao*FSHT*GG

	AE   = PA(1)
	QSE  = PA(2)
	QTE  = PA(3)
	TIE  = PA(4)
	SIE  = PA(5)
	SITE = PA(6)
	PJL  = PA(7)
	SJWE = PA(9)

	SA(1) = SA(1) + AE*EPS(1)  + QSE*EPS(5)  - QTE*EPS(6)   !AXIAL FORCE
	SA(5) = SA(5) + QSE*EPS(1) + SIE*EPS(5)  - SITE*EPS(6)  !MOMENT ABOUT S-AXIS (MY FOR XY PLANE)
	SA(6) = SA(6) - QTE*EPS(1) - SITE*EPS(5) + TIE*EPS(6)   !MOMENT ABOUT T-AXIS (MZ FOR XY PLANE)
	SA(4) = SA(4) + PJL*EPS(4)                              !TORSIONAL MOMENT
	SA(7) = SA(7) + SJWE*EPS(7)                             !WARPING MOMENT

	
	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FIBOPT1(NPM,NPI,NPIB,PA,SA,EPS,LOPT1,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(1)  WORKING ON MECHANICAL DATA
C	 0 = NO UPDATE STRESS & MATERIAL STATE
C	 1 = UPDATE STRESS & MATERIAL STATE

	DIMENSION EPS(7),PA(1),SA(1),AAX(1)
	DIMENSION NPM(1),NPI(1)
	ALLOCATABLE AF5(:)


	N5 = NPM(1)
	ALLOCATE(AF5(N5))

      EPS1 = EPS(1)
      EPS4 = EPS(4)
      EPS5 = EPS(5)
      EPS6 = EPS(6)
      EPS7 = EPS(7)
C      IF(ABS(EPS1).LT.1.0E-10) EPS1 = 0.0D0
C      IF(ABS(EPS4).LT.1.0E-10) EPS4 = 0.0D0
C      IF(ABS(EPS5).LT.1.0E-10) EPS5 = 0.0D0
C      IF(ABS(EPS6).LT.1.0E-10) EPS6 = 0.0D0
C      IF(ABS(EPS7).LT.1.0E-10) EPS7 = 0.0D0
      
	DO 1000 IFIB = 1,NPIB

	NPI(2) = IFIB
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART

C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE

	DA      = AF5(1)  !FIBER AREA
	SC      = AF5(2)  !FIBER DISTANCE S
	TC      = AF5(3)  !FIBER DISTANCE T
	MATN    = INT(AF5(4))
	MATY    = INT(AF5(5))
	WP      = AF5(6)  !WARPING COEF
	TR      = AF5(7)  !TORSION COEF
	FSHS    = AF5(8)
	FSHT    = AF5(9)


	MATFL   = INT(AF5(16))
	EMOD    = AF5(17)
	POISN   = AF5(18)
	GG      = AF5(19)
	RHO     = AF5(20)
	FPC     = AF5(21)
	FPT     = AF5(22)
	EMODP   = AF5(26)     !YOUNG MODULUS RESPECTED TO CODE

	EMP		= AF5(30)     !TOTAL MECHANICAL STRAIN
	ENP		= AF5(31)     !TOTAL NON-MECHANICAL STRAIN
	SCP		= AF5(32)     !OLD STRESS
	EPSK	= AF5(33)     !TOTAL SHRNK STRAIN
	EPCP	= AF5(34)     !TOTAL CREEP STRAIN
	SINC    = AF5(35)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION
C	----------------------------------------	


C	DEFINE TOTAL CONCRETE FIBER STRAIN
	EC  = EPS1  + TC*EPS5 - SC*EPS6 + WP*EPS7

	EMC = EC  - ENP  !TOTAL CURRENT MECHANICAL STRAIN
	DMC = EMC - EMP  !INCR  CURRENT MECHANICAL STRAIN
	EMP = EMC        !UPDATE  TOTAL MECHANICAL STRAIN

C	--------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MATY)

C	-------------
	CASE(0) !CONC	
C	TRACKING STRESS FROM CORRESPONDING MATERIAL STATE (TO BE MODIFIED)
	MATFL = 1
	EMOD  = EMODP !*********
	SCR   = SCP + EMOD*DMC
	EMOD  = EMOD
	GG    = EMOD/(2.0+2.0*POISN)


C	WRITE(*,*) EC,ENP,EMC 

	CALL RELFILL('%FNC',C1,11,MATN,0) !CALL CONC. PROPERTIES
	CALL RELFILL('%FNC',C2,12,MATN,0)
	CALL RELFILL('%FNC',R1,13,MATN,0)

C	STRESS CHANGE VALUE FOR CREEP CALCULATION
	SN1 = SCR
      IF(SN1.LT.0.D0) THEN
         IF(ABS(SN1) .GT. R1*FPC) SN1 = C1*SN1 + C2*FPC
      ENDIF
	SN2 = SCP
      IF(SN2.LT.0.D0) THEN
         IF(ABS(SN2) .GT. R1*FPC) SN2 = C1*SN2 + C2*FPC
      ENDIF
      SINC = SN1 - SN2

C	-------------	
	CASE(1) !STEEL

C	TRACKING STRESS FROM CORRESPONDING MATERIAL STATE (TO BE MODIFIED)
	MATFL = 1
	CALL RELFILL('%MSL',EMOD,1,MATN,0) !***********
	SCR   = SCP + EMOD*DMC
	EMOD  = EMOD
	GG    = EMOD/(2.0+2.0*POISN)
	SINC  = 0.0D0

	ENDSELECT
C	--------------------------------------------

	SCP = SCR  !UPDATE STRESS HERE

C	SECTION RIGIDITY
	A    = DA
	QS   = TC*DA
	QT   = SC*DA
	SIT  = SC*SC*DA
	SIS  = TC*TC*DA
	SIST = SC*TC*DA
	SJ   = SIS + SIT

	
	PA(1) = PA(1) + A*EMOD
	PA(2) = PA(2) + QS*EMOD
	PA(3) = PA(3) + QT*EMOD
	PA(4) = PA(4) + SIS*EMOD
	PA(5) = PA(5) + SIT*EMOD
	PA(6) = PA(6) + SIST*EMOD
	PA(7) = PA(7) + SJ*GG - TR*GG*DA
	PA(8) = PA(8) + A*RHO                           !MASS DENSITY

	PA(9)  = PA(9)  + WP*WP*EMOD*DA		!EIw
	PA(10) = PA(10) + WP*SC*EMOD*DA		!EIsw
	PA(11) = PA(11) + WP*TC*EMOD*DA		!EItw
	PA(12) = PA(12) + WP*EMOD*DA		!EIw

	PA(13) = PA(13) + FSHS*GG*DA
	PA(14) = PA(14) + FSHT*GG*DA

C	MECHANICAL STRESS CONTRIBUTION
	SA(1) = SA(1) + SCR*DA							!AXIAL FORCE
	SA(5) = SA(5) + TC*SCR*DA						!MOMENT ABOUT S-AXIS (MY FOR XY PLANE)
	SA(6) = SA(6) - SC*SCR*DA						!MOMENT ABOUT T-AXIS (MZ FOR XY PLANE)
	SA(4) = SA(4) + GG*(SC*SC + TC*TC - TR)*EPS4*DA	!TORSIONAL MOMENT
	SA(7) = SA(7) + WP*SCR*DA							!WARPING MOMENT

C	UPDATE WORKING ARRAY
C	----------------------------------------
	AF5(16) = FLOAT(MATFL) !MATERIAL STATE
	AF5(17) = EMOD  !UPDATE YOUNG MODULUS
	AF5(18) = POISN !UPDATE POISN
	AF5(19) = GG    !UPDATE SHEAR MODULUS
	AF5(30) = EMP   !TOTAL MECHANICAL STRAIN
	AF5(31) = ENP   !TOTAL NON-MECHANICAL STRAIN
	AF5(32) = SCP   !OLD STRESS
	AF5(33) = EPSK  !TOTAL SHRNK STRAIN
	AF5(34) = EPCP  !TOTAL CREEP STRAIN
	AF5(35) = SINC  !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION

C	WRITE(*,*) IFIB,SCP,SINC
C	----------------------------------------

C	STORE WORKING ARRAY
	IF(LOPT1.EQ.1) CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    

1000	CONTINUE

	DEALLOCATE(AF5)

	
	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FRMOPT2(IEG,IEL,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

	DIMENSION NPM(10),NPI(10),DATF(2),IAX(1),AAX(1)
	ALLOCATABLE IF3(:),AF4(:),AF5(:),AF7(:)
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) ! MAXP
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
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
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !

      CALL INTFILL('%SOL',LOPT2,1,2,0)

C	----------------------------------------
	ALLOCATE(IF3(N3),AF4(N4),AF5(N5),AF7(N7))


	CALL RELFILL('-SOV',TODAY,1,1,0) !CALL TODAY
C	---------------------------------------------------------
C	CALLING SECTION SET NUMBER
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')

	ISEC = IF3(2)     !SECTION SET
C	ICRP = IF3(3)     !CREEP INTEGRATION OPTION  !ALREADY CALL IN MATERIAL
	LSTY = IF3(5)     !SEC TYPE FLAG
	NPRT = IF3(6)     !NUMBER OF PART
	
C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ----------------------------------------------------------
      DO 400 IGR = 1,NGR

	NPM(1:2) = [N7,NGR]
	NPI(1:2) = [N7,IGR]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  
	N7_R = 1           !POINTER OF RIGIDITY
	N7_S = 1 + 20       !POINTER OF MECHANICAL STRESS
	N7_N = 1 + 20 + 8   !POINTER OF NON-MECHANICAL STRESS
	AF7(N7_N:N7_N+8-1) = 0.0D0	!INITIALIZE OF NON-MECHANICAL STRESS


	DO 300 IPRT = 1,NPRT


	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
	IXIST   = INT(AF4(1))   !EXISTING FLAG
	NPIB    = INT(AF4(2))   !NUM FIBER FOR THIS PART
	DATF(1) = AF4(3)        !CASTING DAY
	DATF(2) = AF4(4)        !INSTALLATION DAY


	IF (IXIST.EQ. 0) GOTO 400 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED
      AGEC  = 0.0001
	IF(TODAY-DATF(1).GT.AGEC) AGEC = TODAY-DATF(1)
      TINC  = TODAY - DATF(2)

	IF(IGR.EQ.NGR) AF4(4) = TODAY  !BACKUP TODAY TO THIS PART (AT LAST GAUSS)
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')

	SELECTCASE(LSTY)

	CASE(0) !READY PROP
	IF(LOPT2.EQ.1) THEN
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,1   ,IGR ,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
	CALL REDOPT2(AGEC,TINC,AF7(N7_N),AF5)
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    !BACKUP  DATA OF THIS PART
	ENDIF

	CASE(1) !FIBER PROP
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   0,IGR ,IPRT]   !NFIB TO BE MODIFIED IN NEXT SUBRUTINE
	IF(LOPT2.EQ.1) CALL FIBOPT21(AGEC,TINC,NPM,NPI,NPIB,
	1							 AF7(N7_N),AAX,NDX5)
	IF(LOPT2.EQ.2) CALL FIBOPT22(NPM,NPI,NPIB,AAX,NDX5)
	IF(LOPT2.EQ.3) CALL FIBOPT23(NPM,NPI,NPIB,AF7(N7_N),AAX,NDX5)

	ENDSELECT

300	CONTINUE


	NPM(1:2) = [N7,NGR]
	NPI(1:2) = [N7,IGR]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT')   !BACKUP HERE NON-MECHANICAL STRESS
C     ----------------------------------------------------------
400	CONTINUE
C     ----------------------------------------------------------

	DEALLOCATE(IF3,AF4,AF5,AF7)


	RETURN

	END




C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE GUSOPT1(NPM,NPI,ISEC,PA,SA,EPS,LOPT1,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(1)  WORKING ON MECHANICAL DATA
C	 0 = NO UPDATE STRESS & MATERIAL STATE
C	 1 = UPDATE STRESS & MATERIAL STATE


	DIMENSION EPS(6),PA(1),SA(1)
	DIMENSION NPM(1),NPI(1),AAX(1)

	DIMENSION EMP(6),ENP(6),SCP(6),EPSK(6),EPCP(6),SINC(6)
	DIMENSION EMC(6),DMC(6),SCR(6),DMATX(6,6),DP(6,6),EC(6)

C	FOR FLOW VECTOR
	DIMENSION AVECT(6)

	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
	CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN
	
	MATN = INT(FMAT)  !MAT NO
	MTYP = INT(FTYP)	!MAT TYPE  0=CONCRETE 1=STEEL
	MATS = INT(FATS)  !MAT No. FOR STEEL
	
C	=========================================================	
C     CONCRETE LAYER
C	=========================================================	
	NPI(2) = 1   ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE
      NSTAT	= INT(AF5(1)) !MATERIAL STATE NUMBER
	EMOD    = AF5(2)
	POISN   = AF5(3)
	FPC     = AF5(4)
	FPT     = AF5(5)
	EMODP   = AF5(6)     !YOUNG MODULUS RESPECTED TO CODE
C	FUNLO   = AF5(7)     !UNLOADING FLAG --> 0.0 FOR UNLOAD, 1.0 FOR LOADING 
      ! 6 STRAINS IN SOLID ELEMENT
	EMP(1:6) = AF5(11:16)     !TOTAL MECHANICAL STRAIN
	ENP(1:6) = AF5(17:22)     !TOTAL NON-MECHANICAL STRAIN
	SCP(1:6) = AF5(23:28)     !OLD STRESS
	EPSK(1:6)= AF5(29:34)     !TOTAL SHRNK STRAIN
	EPCP(1:6)= AF5(35:40)     !TOTAL CREEP STRAIN
	SINC(1:6)= AF5(41:46)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION	
	 
C	--------------------------------------------
C	           MECHANICAL STRAINS
C	--------------------------------------------
C     MECHANICAL STRAIN (Em = E - En)
	DO 11 I = 1,6  
11	EMC(I)  = EPS(I)  - ENP(I)

C     INCREMENTAL MECHANICAL STRAIN (DEm = Em - Em0)
	DO 12 I = 1,6  
12	DMC(I)  = EMC(I) - EMP(I)
C     UPDATE  TOTAL MECHANICAL STRAIN (Em0 = Em)
	EMP(1:6) = EMC(1:6) 

C	--------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	============================================
	CASE(0) !CONC	
C	============================================	
C	TRACKING STRESS FROM CORRESPONDING MATERIAL STATE (TO BE MODIFIED)
	NSTAT = 1
	EMOD  = EMODP !*********

	CALL SOMODUC(DMATX,EMOD,POISN)

	CALL RELFILL('%FNC',C1,11,MATN,0) !CALL CONC. PROPERTIES
	CALL RELFILL('%FNC',C2,12,MATN,0)
	CALL RELFILL('%FNC',R1,13,MATN,0)
	
C	--------------------------------------------
C                 MECHANICAL STRESS
C	--------------------------------------------
	DO I = 1,6
C	SCR(I) = SCP(I) !INCREMENTAL FORM
	SCR(I) = 0.0D0  !TOTAL FORM
	DO J = 1,6
C	SCR(I) = SCR(I) + DMATX(I,J)*DMC(J) !INCREMENTAL FORM
	SCR(I) = SCR(I) + DMATX(I,J)*EMC(J) !TOTAL FORM
	ENDDO
	ENDDO
	SCR = MATMUL(DMATX,EMC)

C	DO I = 1,2
C	IF(ABS(SCR(I)).GT.0.6*FPC.AND.SCR(I).NE.0.0D0) THEN
C	SCR(I) = 0.6*FPC*SCR(I)/ABS(SCR(I))
C	ENDIF
C	ENDDO

C	--------------------------------------------
C            INCREMENTAL MECHANICAL STRESS
C	--------------------------------------------
C     NONLINEAR CREEP EFFECTS AT HIGH STRESS LEVELS
	DO I = 1,6
	SN1 = SCR(I)
      IF(SN1.LT.0.D0) THEN
         IF(ABS(SN1) .GT. R1*FPC) SN1 = C1*SN1 + C2*FPC 
      ENDIF
	SN2 = SCP(I)
      IF(SN2.LT.0.D0) THEN
         IF(ABS(SN2) .GT. R1*FPC) SN2 = C1*SN2 + C2*FPC
      ENDIF
	SINC(I) = SN1 - SN2
	ENDDO
C	----------------------------------------------------------	
C     DISPLACEMENTS 
      !WRITE(IWRIT,*) 'SINC1 = ',SINC(1)         
C	----------------------------------------------------------		
		
C	--------------------------------------------	
C	DRUCKER PRAGER YIELD SURFACE
C	--------------------------------------------	
	ALPI1 = 0.0D0
	CALL SOYSURF(SCR,SN1,AVECT,ALPI1)  !GET SN1 HERE
	CALL SOYSURF(SCP,SN2,AVECT,ALPI1)  !GET SN2  HERE
	FUNLO = 0.0D0
	IF(SN1.GE.SN2) FUNLO = 1.0D0

C	============================================	
	CASE(1) !STEEL
C	============================================	

C	TRACKING STRESS FROM CORRESPONDING MATERIAL STATE (TO BE MODIFIED)
	NSTAT = 1
	CALL RELFILL('%MSL',EMOD,1,MATN,0) !***********
	CALL RELFILL('%MSL',POISN,8,MATN,0) !***********
	CALL SOMODUC(DMATX,EMOD,POISN)

	DO I = 1,6
C	SCR(I) = SCP(I) !INCREMENTAL FORM
	SCR(I) = 0.0D0  !TOTAL FORM
	DO J = 1,6
C	SCR(I) = SCR(I) + DMATX(I,J)*DMC(J) !INCREMENTAL FORM
	SCR(I) = SCR(I) + DMATX(I,J)*EMC(J) !TOTAL FORM
	ENDDO
	ENDDO
	SCR = MATMUL(DMATX,EMC)
	SINC(1:6) = 0.0D0
	
C	============================================
	ENDSELECT
	
C     UPDATE STRESS HERE		
	SCP(1:6) = SCR(1:6)  

C	------------------------------------------------------------------
C	INITIALIZE ELASTO-PLASTIC RIGIDITY MATRIX	
C	------------------------------------------------------------------
	DP(1:6,1:6) = 0.0D0
	DO I = 1,3
	DO J = 1,3	
	DP(I,J) = DMATX(I,J)
	ENDDO
	ENDDO
	DO I = 4,6	
	DP(I,I) = DMATX(I,I)	
	ENDDO
	
C	------------------------------------------------------------------
C	RIGIDITY ARRAY
C	------------------------------------------------------------------
	K = 0
	DO I = 1,3
	DO J = 1,3
	K = K + 1
	PA(K) = DP(I,J)
	ENDDO
	ENDDO
	DO I = 4,6	
	K = K + 1
	PA(K) = DP(I,I)	
	ENDDO
	
C	------------------------------------------------------------------
C	MECHANICAL STRESS CONTRIBUTION
C	------------------------------------------------------------------
	DO I = 1,6	
	SA(I) = SA(I) + SCR(I)
      ENDDO	

C	UPDATE WORKING ARRAY
C	----------------------------------------	
      AF5(1)	= FLOAT(NSTAT)	!MATERIAL STATE NUMBER
	AF5(2)	= EMOD  !UPDATE YOUNG MODULUS    
	AF5(3)	= POISN !UPDATE POISN  
	AF5(7)    = FUNLO !UNLOADING FLAG --> 0.0 FOR UNLOAD, 1.0 FOR LOADING 
	AF5(11:16) = EMP(1:6)     !TOTAL MECHANICAL STRAIN
	AF5(17:22) = ENP(1:6)     !TOTAL NON-MECHANICAL STRAIN
	AF5(23:28) = SCP(1:6)     !OLD STRESS
	AF5(29:34) = EPSK(1:6)    !TOTAL SHRNK STRAIN
	AF5(35:40) = EPCP(1:6)    !TOTAL CREEP STRAIN
	AF5(41:46) = SINC(1:6)    !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION
C	----------------------------------------	
C	STORE WORKING ARRAY
	IF(LOPT1.EQ.1) THEN
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    !UPDATE DATA OF THIS LAYER	   
	ELSEIF(LOPT1.EQ.2) THEN	!NO BACK UP ...ONLY CALLING STIFFNESS  (BACKUP ONLY STRESS) JAN09 LOADCOMB	
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !UPDATE DATA OF THIS LAYER (ONLY STRESS) 
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    !UPDATE DATA OF THIS LAYER (ONLY STRESS) 
	ENDIF
C	----------------------------------------	

C	=================================================== CONCRETE 	

C	--------------------------------------------------------------
C	START FOR STEEL LAYER (WAITING FOR MODIFIED)
C	--------------------------------------------------------------
	IF(SRTO.LE.0.0D0) THEN
	DEALLOCATE(AF5)
	RETURN
	ENDIF
	
C	--------------------------------------------------------------	

	DEALLOCATE(AF5)

	RETURN
	END
	
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE GUCOPT1(NPM,NPI,ISEC,PA,SA,EPS,LOPT1,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(1)  WORKING ON MECHANICAL DATA
C	 0 = NO UPDATE STRESS & MATERIAL STATE
C	 1 = UPDATE STRESS & MATERIAL STATE


	DIMENSION EPS(6),PA(1),SA(1)
	DIMENSION NPM(1),NPI(1),AAX(1)

	DIMENSION EMP(6),ENP(6),SCP(6),EPSK(6),EPCP(6),SINC(6)
	DIMENSION EMC(6),DMC(6),SCR(6),DMATX(6,6),DP(6,6),EC(6)

C	FOR FLOW VECTOR
	DIMENSION AVECT(6)

	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
	CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN
	
	MATN = INT(FMAT)  !MAT NO
	MTYP = INT(FTYP)	!MAT TYPE  0=CONCRETE 1=STEEL
	MATS = INT(FATS)  !MAT No. FOR STEEL
	
C	=========================================================	
C     CONCRETE LAYER
C	=========================================================	
	NPI(2) = 1   ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE
      NSTAT	= INT(AF5(1)) !MATERIAL STATE NUMBER
	EMOD    = AF5(2)
	POISN   = AF5(3)
	FPC     = AF5(4)
	FPT     = AF5(5)
	EMODP   = AF5(6)     !YOUNG MODULUS RESPECTED TO CODE
C	FUNLO   = AF5(7)     !UNLOADING FLAG --> 0.0 FOR UNLOAD, 1.0 FOR LOADING 
      ! 6 STRAINS IN SOLID ELEMENT
	EMP(1:6) = AF5(11:16)     !TOTAL MECHANICAL STRAIN
	ENP(1:6) = AF5(17:22)     !TOTAL NON-MECHANICAL STRAIN
	SCP(1:6) = AF5(23:28)     !OLD STRESS
	EPSK(1:6)= AF5(29:34)     !TOTAL SHRNK STRAIN
	EPCP(1:6)= AF5(35:40)     !TOTAL CREEP STRAIN
	SINC(1:6)= AF5(41:46)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION	
	 
C	--------------------------------------------
C	           MECHANICAL STRAINS
C	--------------------------------------------
C     MECHANICAL STRAIN (Em = E - En)
	DO 11 I = 1,6  
11	EMC(I)  = EPS(I)  - ENP(I)

C     INCREMENTAL MECHANICAL STRAIN (DEm = Em - Em0)
	DO 12 I = 1,6  
12	DMC(I)  = EMC(I) - EMP(I)
C     UPDATE  TOTAL MECHANICAL STRAIN (Em0 = Em)
	EMP(1:6) = EMC(1:6) 

C	--------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	============================================
	CASE(0) !CONC	
C	============================================	
C	TRACKING STRESS FROM CORRESPONDING MATERIAL STATE (TO BE MODIFIED)
	NSTAT = 1
	EMOD  = EMODP !*********

	CALL SCMODUC(DMATX,EMOD,POISN)

	CALL RELFILL('%FNC',C1,11,MATN,0) !CALL CONC. PROPERTIES
	CALL RELFILL('%FNC',C2,12,MATN,0)
	CALL RELFILL('%FNC',R1,13,MATN,0)
	
C	--------------------------------------------
C                 MECHANICAL STRESS
C	--------------------------------------------
	DO I = 1,6
C	SCR(I) = SCP(I) !INCREMENTAL FORM
	SCR(I) = 0.0D0  !TOTAL FORM
	DO J = 1,6
C	SCR(I) = SCR(I) + DMATX(I,J)*DMC(J) !INCREMENTAL FORM
	SCR(I) = SCR(I) + DMATX(I,J)*EMC(J) !TOTAL FORM
	ENDDO
	ENDDO
	SCR = MATMUL(DMATX,EMC)

C	DO I = 1,2
C	IF(ABS(SCR(I)).GT.0.6*FPC.AND.SCR(I).NE.0.0D0) THEN
C	SCR(I) = 0.6*FPC*SCR(I)/ABS(SCR(I))
C	ENDIF
C	ENDDO

C	--------------------------------------------
C            INCREMENTAL MECHANICAL STRESS
C	--------------------------------------------
C     NONLINEAR CREEP EFFECTS AT HIGH STRESS LEVELS
	DO I = 1,6
	SN1 = SCR(I)
      IF(SN1.LT.0.D0) THEN
         IF(ABS(SN1) .GT. R1*FPC) SN1 = C1*SN1 + C2*FPC 
      ENDIF
	SN2 = SCP(I)
      IF(SN2.LT.0.D0) THEN
         IF(ABS(SN2) .GT. R1*FPC) SN2 = C1*SN2 + C2*FPC
      ENDIF
	SINC(I) = SN1 - SN2
	ENDDO
C	----------------------------------------------------------	
C     DISPLACEMENTS 
      !WRITE(IWRIT,*) 'SINC1 = ',SINC(1)         
C	----------------------------------------------------------		
		
C	--------------------------------------------	
C	DRUCKER PRAGER YIELD SURFACE
C	--------------------------------------------	
	ALPI1 = 0.0D0
	CALL SOYSURF(SCR,SN1,AVECT,ALPI1)  !GET SN1 HERE
	CALL SOYSURF(SCP,SN2,AVECT,ALPI1)  !GET SN2  HERE
	FUNLO = 0.0D0
	IF(SN1.GE.SN2) FUNLO = 1.0D0

C	============================================	
	CASE(1) !STEEL
C	============================================	

C	TRACKING STRESS FROM CORRESPONDING MATERIAL STATE (TO BE MODIFIED)
	NSTAT = 1
	CALL RELFILL('%MSL',EMOD,1,MATN,0) !***********
	CALL RELFILL('%MSL',POISN,8,MATN,0) !***********
	CALL SCMODUC(DMATX,EMOD,POISN)

	DO I = 1,6
C	SCR(I) = SCP(I) !INCREMENTAL FORM
	SCR(I) = 0.0D0  !TOTAL FORM
	DO J = 1,6
C	SCR(I) = SCR(I) + DMATX(I,J)*DMC(J) !INCREMENTAL FORM
	SCR(I) = SCR(I) + DMATX(I,J)*EMC(J) !TOTAL FORM
	ENDDO
	ENDDO
	SCR = MATMUL(DMATX,EMC)
	SINC(1:6) = 0.0D0
	
C	============================================
	ENDSELECT
	
C     UPDATE STRESS HERE		
	SCP(1:6) = SCR(1:6)  

C	------------------------------------------------------------------
C	INITIALIZE ELASTO-PLASTIC RIGIDITY MATRIX	
C	------------------------------------------------------------------
	DP(1:6,1:6) = 0.0D0
	DO I = 1,3
	DO J = 1,3	
	DP(I,J) = DMATX(I,J)
	ENDDO
	ENDDO
	DO I = 4,6	
	DP(I,I) = DMATX(I,I)	
	ENDDO
	
C	------------------------------------------------------------------
C	RIGIDITY ARRAY
C	------------------------------------------------------------------
	K = 0
	DO I = 1,3
	DO J = 1,3
	K = K + 1
	PA(K) = DP(I,J)
	ENDDO
	ENDDO
	DO I = 4,6	
	K = K + 1
	PA(K) = DP(I,I)	
	ENDDO
	
C	------------------------------------------------------------------
C	MECHANICAL STRESS CONTRIBUTION
C	------------------------------------------------------------------
	DO I = 1,6	
	SA(I) = SA(I) + SCR(I)
	!SA(I) = SCR(I)
      ENDDO	

C	UPDATE WORKING ARRAY
C	----------------------------------------	
      AF5(1)	= FLOAT(NSTAT)	!MATERIAL STATE NUMBER
	AF5(2)	= EMOD  !UPDATE YOUNG MODULUS    
	AF5(3)	= POISN !UPDATE POISN  
	AF5(7)    = FUNLO !UNLOADING FLAG --> 0.0 FOR UNLOAD, 1.0 FOR LOADING 
	AF5(11:16) = EMP(1:6)     !TOTAL MECHANICAL STRAIN
	AF5(17:22) = ENP(1:6)     !TOTAL NON-MECHANICAL STRAIN
	AF5(23:28) = SCP(1:6)     !OLD STRESS
	AF5(29:34) = EPSK(1:6)    !TOTAL SHRNK STRAIN
	AF5(35:40) = EPCP(1:6)    !TOTAL CREEP STRAIN
	AF5(41:46) = SINC(1:6)    !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION
C	----------------------------------------	
C	STORE WORKING ARRAY
	IF(LOPT1.EQ.1) THEN
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    !UPDATE DATA OF THIS LAYER	   
	ELSEIF(LOPT1.EQ.2) THEN	!NO BACK UP ...ONLY CALLING STIFFNESS  (BACKUP ONLY STRESS) JAN09 LOADCOMB	
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !UPDATE DATA OF THIS LAYER (ONLY STRESS) 
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    !UPDATE DATA OF THIS LAYER (ONLY STRESS) 
		
	ENDIF
C	----------------------------------------	

C	=================================================== CONCRETE 	

C	--------------------------------------------------------------
C	START FOR STEEL LAYER (WAITING FOR MODIFIED)
C	--------------------------------------------------------------
	IF(SRTO.LE.0.0D0) THEN
	DEALLOCATE(AF5)
	RETURN
	ENDIF
	
C	--------------------------------------------------------------	

	DEALLOCATE(AF5)

	RETURN
	END	
	


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SODOPT2(IEG,IEL,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

	DIMENSION NPM(10),NPI(10),DATF(2),IAX(1),AAX(1),STn(6)
	ALLOCATABLE IF3(:),AF4(:),AF7(:)
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) ! MAXP
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
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
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !

C	SOLUTION CONTROL (NON-MECHANICAL STRAIN)
      CALL INTFILL('%SOL',LOPT2,1,2,0)

C	----------------------------------------
C     IF3 = GENERAL SECTION AND CREEP DATA 
C     AF4  = CASTING DAY, INSTALLATION DAY
C     AF5  = WORKING ARRAY FOR FIBER/LAYER
C     AF7  = WORKING ARRAY FOR GAUSS POINT
C	----------------------------------------
	ALLOCATE(IF3(N3),AF4(N4),AF7(N7))

	CALL RELFILL('-SOV',TODAY,1,1,0) !CALL TODAY

C	CALLING SECTION SET NUMBER
C     IF3-R(SECTION) ------------------------------------------- 
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')
	ISEC = IF3(2)     !SECTION SET (ISEC = 1)
C	ICRP = IF3(3)     !CREEP INTEGRATION OPTION  !ALREADY CALL IN MATERIAL
	NPRT = IF3(6)     !NUMBER OF PART (NPRT = 1)
C	----------------------------------------------------------	
	
		
C     ==========================================================
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ==========================================================

      DO 400 IPT = 1,MGAS
      
C     AF7-R(NON-MECH STRAIN)------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  
	N7_R = 1            !POINTER OF RIGIDITY
	N7_S = 1 + 12       !POINTER OF MECHANICAL STRESS
	N7_N = 1 + 12 + 6   !POINTER OF NON-MECHANICAL STRESS
	AF7(N7_N:N7_N+6-1) = 0.0D0	!INITIALIZE OF NON-MECHANICAL STRESS
C	----------------------------------------------------------
      
      DO 300 IPRT = 1,NPRT
      
C     AF4-R(DAYS)-----------------------------------------------        
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
	IXIST   = INT(AF4(1))   !EXISTING FLAG
	DATF(1) = AF4(3)        !CASTING DAY
	DATF(2) = AF4(4)        !INSTALLATION DAY
C	----------------------------------------------------------
	IF (IXIST.EQ. 0) GOTO 400 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED
      AGEC  = 0.0001
	IF(TODAY-DATF(1).GT.AGEC) AGEC = TODAY-DATF(1)
      TINC  = TODAY - DATF(2)
      
C     AF4-W(DAYS)-----------------------------------------------
	IF(IPT.EQ.MGAS) AF4(4) = TODAY  !BACKUP TODAY TO THIS PART (AT LAST GAUSS)
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')
C	----------------------------------------------------------	

C     AF5,AF7-R(MECH STRESS & NON MECH STRAIN)------------------
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   0,IPT ,IPRT]   !NFIB TO BE MODIFIED IN NEXT SUBRUTINE
	IF(LOPT2.EQ.1) CALL GUSOPT21(AGEC,TINC,NPM,NPI,ISEC,
	1							   AF7(N7_N),AAX,NDX5)
	IF(LOPT2.EQ.2) CALL GUSOPT22(NPM,NPI,ISEC,AAX,NDX5)
	IF(LOPT2.EQ.3) CALL GUSOPT23(NPM,NPI,ISEC,AF7(N7_N),AAX,NDX5)
C	----------------------------------------------------------

300	CONTINUE
C     AF7-W(NON-MECH STRAIN)------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT')   !BACKUP HERE NON-MECHANICAL STRESS
C     ----------------------------------------------------------
      IF(LOPT2.EQ.1) THEN
      DO II = 1,6
      STn(II) = AF7(N7_N+II-1)
      ENDDO
      ENDIF
      
400	CONTINUE
C     ==========================================================


	DEALLOCATE(IF3,AF4,AF7)

	RETURN
	END
	
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SOCOPT2(IEG,IEL,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

	DIMENSION NPM(10),NPI(10),DATF(2),IAX(1),AAX(1),STn(6)
	ALLOCATABLE IF3(:),AF4(:),AF7(:)
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) ! MAXP
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
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
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !

C	SOLUTION CONTROL (NON-MECHANICAL STRAIN)
      CALL INTFILL('%SOL',LOPT2,1,2,0)

C	----------------------------------------
C     IF3 = GENERAL SECTION AND CREEP DATA 
C     AF4  = CASTING DAY, INSTALLATION DAY
C     AF5  = WORKING ARRAY FOR FIBER/LAYER
C     AF7  = WORKING ARRAY FOR GAUSS POINT
C	----------------------------------------
	ALLOCATE(IF3(N3),AF4(N4),AF7(N7))

	CALL RELFILL('-SOV',TODAY,1,1,0) !CALL TODAY

C	CALLING SECTION SET NUMBER
C     IF3-R(SECTION) ------------------------------------------- 
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')
	ISEC = IF3(2)     !SECTION SET (ISEC = 1)
C	ICRP = IF3(3)     !CREEP INTEGRATION OPTION  !ALREADY CALL IN MATERIAL
	NPRT = IF3(6)     !NUMBER OF PART (NPRT = 1)
C	----------------------------------------------------------	
	
		
C     ==========================================================
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ==========================================================

      DO 400 IPT = 1,MGAS
      
C     AF7-R(NON-MECH STRAIN)------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  
	N7_R = 1            !POINTER OF RIGIDITY
	N7_S = 1 + 12       !POINTER OF MECHANICAL STRESS
	N7_N = 1 + 12 + 6   !POINTER OF NON-MECHANICAL STRESS
	AF7(N7_N:N7_N+6-1) = 0.0D0	!INITIALIZE OF NON-MECHANICAL STRESS
C	----------------------------------------------------------
      
      DO 300 IPRT = 1,NPRT
      
C     AF4-R(DAYS)-----------------------------------------------        
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
	IXIST   = INT(AF4(1))   !EXISTING FLAG
	DATF(1) = AF4(3)        !CASTING DAY
	DATF(2) = AF4(4)        !INSTALLATION DAY
C	----------------------------------------------------------
	IF (IXIST.EQ. 0) GOTO 400 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED
      AGEC  = 0.0001
	IF(TODAY-DATF(1).GT.AGEC) AGEC = TODAY-DATF(1)
      TINC  = TODAY - DATF(2)
      
C     AF4-W(DAYS)-----------------------------------------------
	IF(IPT.EQ.MGAS) AF4(4) = TODAY  !BACKUP TODAY TO THIS PART (AT LAST GAUSS)
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')
C	----------------------------------------------------------	

C     AF5,AF7-R(MECH STRESS & NON MECH STRAIN)------------------
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   0,IPT ,IPRT]   !NFIB TO BE MODIFIED IN NEXT SUBRUTINE
	IF(LOPT2.EQ.1) CALL GUCOPT21(AGEC,TINC,NPM,NPI,ISEC,
	1							   AF7(N7_N),AAX,NDX5)
	IF(LOPT2.EQ.2) CALL GUCOPT22(NPM,NPI,ISEC,AAX,NDX5)
	IF(LOPT2.EQ.3) CALL GUCOPT23(NPM,NPI,ISEC,AF7(N7_N),AAX,NDX5)
C	----------------------------------------------------------

300	CONTINUE
C     AF7-W(NON-MECH STRAIN)------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT')   !BACKUP HERE NON-MECHANICAL STRESS
C     ----------------------------------------------------------
      IF(LOPT2.EQ.1) THEN
      DO II = 1,6
      STn(II) = AF7(N7_N+II-1)
      ENDDO
      ENDIF
      
400	CONTINUE
C     ==========================================================


	DEALLOCATE(IF3,AF4,AF7)

	RETURN
	END	
	



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE REDOPT2(TAU,TINC,SA,AF5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION AF5(1),SA(1),CCOF(1000)

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	NO NON- MECHANICAL STRESS FOR READY PROPERTY SECTION

	SA(1) = 0.0D0   !AXIAL FORCE
	SA(4) = 0.0D0   !MOMENT ABOUT S-AXIS (MY FOR XY PLANE)
	SA(5) = 0.0D0   !MOMENT ABOUT T-AXIS (MZ FOR XY PLANE)
	SA(6) = 0.0D0   !TORSIONAL MOMENT
	SA(7) = 0.0D0   !WARPING MOMENT

C	UPDATE MATERIAL STATE
	MATN = INT(AF5(11))
	MATY = INT(AF5(12))
	IF(MATY.EQ.0) THEN
	TAUO = TAU-TINC
	CALL TMPCNC(TAUO,MATN,EMODP,SHRNO,FPC,FPT,RHO,POISN,GG,CCOF) !ONLY FOR INITIAL SHIRNKAGE STRAIN SHRNO AT THE TIME OF ERECTION UPDATE EMODP HERE
C	CALL TMPCNC(TAU ,MATN,EMODC,SHRN ,FPC,FPT,RHO,POISN,GG,CCOF) 
	ENDIF
	IF(MATY.EQ.1) THEN
	CALL RELFILL('%MSL',EMODC,1,MATN,0)
	CALL RELFILL('%MSL',RHO  ,6,MATN,0)
	CALL RELFILL('%MSL',POISN,8,MATN,0)
	GG  = EMODC/(2.0+2.0*POISN)
	FPC = 0.0D0
	FPT = 0.0D0
	ENDIF

	EMOD    = EMODP
	AF5(17) = EMOD        
	AF5(18) = POISN 
	AF5(19) = GG   

	AF5(20) = RHO     
	AF5(21) = FPC     
	AF5(22) = FPT    
	AF5(26) = EMODP       !YOUNG MODULUS RESPECTED TO CODE

	
	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FIBOPT21(TAU,TINC,NPM,NPI,NPIB,SA,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL AT START OF TIME STEP
C	IM = CREEP INTEGRATION OPTION  

	DIMENSION SA(1),NPM(1),NPI(1),AAX(1)
	DIMENSION RLXK(2)
	ALLOCATABLE CCOP(:,:),CCOF(:,:),GTJ(:),APHI(:,:)
	ALLOCATABLE ELIT(:),ALFA(:,:)
	ALLOCATABLE AF5(:)


	N5 = NPM(1)
	ALLOCATE(AF5(N5))


	DO 1000 IFIB = 1,NPIB



	NPI(2) = IFIB
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART

C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE

	DA      = AF5(1)  !FIBER AREA
	SC      = AF5(2)  !FIBER DISTANCE S
	TC      = AF5(3)  !FIBER DISTANCE T
	MATN    = INT(AF5(4))
	MATY    = INT(AF5(5))
	WP      = AF5(6)  !WARPING COEF
	TR      = AF5(7)  !TORSION COEF

	MATFL   = INT(AF5(16))
	EMOD    = AF5(17)
	POISN   = AF5(18)
	GG      = AF5(19)
	RHO     = AF5(20)
	FPC     = AF5(21)
	FPT     = AF5(22)
	SHRN    = AF5(23) 
	RLXK(1) = AF5(24) 
	RLXK(2) = AF5(25) 
	EMODP   = AF5(26)     !YOUNG MODULUS RESPECTED TO CODE
	FSRK    = AF5(27)     !SHRINKAGE FACTOR
	
	EMP		= AF5(30)     !TOTAL MECHANICAL STRAIN
	ENP		= AF5(31)     !TOTAL NON-MECHANICAL STRAIN
	SCP		= AF5(32)     !OLD STRESS
	EPSK	= AF5(33)     !TOTAL SHRNK STRAIN
	EPCP	= AF5(34)     !TOTAL CREEP STRAIN
	SINC    = AF5(35)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION
C	----------------------------------------	

C	--------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MATY)

C	-------------
	CASE(0) !CONC	

	CALL INTFILL('%CPT',KM ,1,MATN,0)
	CALL INTFILL('%CPT',NAT,2,MATN,0)
	CALL INTFILL('%CPT',NTT,3,MATN,0)
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)

      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG


	ALLOCATE(CCOP(NRT,KL),CCOF(NRT,KL),GTJ(NRT),APHI(NRT,KL))
	ALLOCATE(ELIT(NRT),ALFA(NRT,KL))

	TAUO = TAU-TINC
	CALL TMPCNC(TAUO,MATN,EMODP,SHRNO,FPC,FPT,RHO,POISN,GG,CCOF) !ONLY FOR INITIAL SHIRNKAGE STRAIN SHRNO AT THE TIME OF ERECTION UPDATE EMODP HERE
	CALL TMPCNC(TAU ,MATN,EMODC,SHRN ,FPC,FPT,RHO,POISN,GG,CCOF) 

C	FIND SHRINKAGE STRAIN INCREMENT
	DESK  = 0.0D0
	ECSK  = SHRN*FSRK
      IF(TINC.NE.0.D0) THEN
C	EPSK = ECSK         !CHANGE TO NEXT TO ACCOUNT THE EFFECT OF INITIAL SHIRNKAGE STRAIN SHRNO AT THE TIME OF ERECTION
	EPSK = SHRNO*FSRK  
	DESK = ECSK - EPSK
	ENDIF

C	FIND AGING STRAIN INCREMENT
	DEAG = 0.0D0
      IF(TINC.NE.0.D0) THEN
	DEAG = SCP/EMODP - SCP/EMODC
      ENDIF

C	---------------------------------------
C	FIND CREEP STRAIN INCREMENT
C	---------------------------------------
C	----------------------------
C	FIND CREEP RELATED VARIABLES
	CALL RELFILL('%FNC',GAMA,6,MATN,0)
      GAMA = GAMA * 1.0D-1


	K = 50
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	CCOP(I,J) = AF5(K)   !CREEP COEF  (51-80) ALLOCATE 30 VALUE  
	ENDDO
	ENDDO
	K = 80
	DO I = 1,NRT
	K = K + 1
	GTJ(I)    = AF5(K)   !CREEP COEF  (81-95) ALLOCATE 15 VALUE  
	ENDDO	
	K = 100
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	ALFA(I,J) = AF5(K)   !CREEP COEF  (101-130) ALLOCATE 30 VALUE  
	ENDDO
	ENDDO
	K = 130
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	APHI(I,J) = AF5(K)   !CREEP COEF  (131-160) ALLOCATE 30 VALUE 
	ENDDO
	ENDDO

	RLXK(1:2)   = 0.0D0
C
      DO 13 J=1,NRT
      GAMA = GAMA * 1.0D+1
      YIOT = TINC / GAMA
      ELIM = DEXP(-YIOT)
      ELIT(J) = 1.0D+0 - ELIM

	LZE =0
	IF(ABS(YIOT).LT.1.0E-15) THEN
	YIOT = 0.0D0
	LZE = 1
	ENDIF
C
	PHIM = 0.0D0
      IF(LZE.EQ.0) 
	1 PHIM = 1.0D+0 - ( ELIT(J) / YIOT )         !(1-Phi)  Eq-3.3.14b
      IF(IM.EQ.3) BBIT = (YIOT / 2.0D+0) - PHIM   !(dt/2G)-(1-Phi)  Eq-3.319b
C
      DO 12 K=1,KL
      RLXF = CCOF(J,K) * PHIM                     !Ai(1-Phi) Eq-3.3.14b        ***Used CCOF = CCOP
	APHI(J,K) = CCOF(J,K)                       !A(tj-1)   Eq 3.3.11d  C-C   ***Used CCOF = CCOP
	IF(IM.LT.2) GOTO 12	
	APHI(J,K) = APHI(J,K)*(1.0-PHIM)            !A(tj-1)*phi Eq 3.3.14d C-L
	IF(IM.NE.3) GOTO 12							!L-L CREEP INTEGRATION METHOD
	APHI(J,K) = APHI(J,K) + ALFA(J,K)*PHIM      !A(tj-1)*phi + B(tj-1)*(1-phi) Eq 3.3.19e L-L
	ALFA(J,K) = 0.0D0
      IF(LZE.EQ.0) 
	1 ALFA(J,K) = (CCOF(J,K) - CCOP(J,K)) / YIOT !Beta=(dA)*G/dt   Eq-3.3.19d
      RLXF = RLXF + ALFA(J,K) * BBIT              !Add term   Beta*((dt/2G)-(1-Phi))  Eq-3.3.19b
12	RLXK(K) = RLXK(K) + RLXF                    !Sum over NRT  Eq-3.3.14b and  3.3.19b
C
13	CONTINUE
C	----------------------------
	DECR = 0.0D0
      DO 45 J=1,NRT
	DECR = DECR + GTJ(J) * ELIT(J)              !The 2nd term in Eq-3.3.11a  3.3.14a  3.3.19a  ,GTJ = g(tj-1)  
45	GTJ(J) = GTJ(J)*(1.0-ELIT(J))               !UPDATE g(tj) 
      EPCP = EPCP + DECR                          !UPDATE CREEP STRAIN
C	---------------------------------------

C	SUM OF INCREMENTAL NON-MECHANICAL STRAIN
	DENN = DESK + DEAG + DECR

C	UPDATE TOTAL NON-MECHANICAL STRAIN
	ENP  = ENP  + DENN

C-----ASSEMBLE EQUIVALENT ELEMENT LOADS DUE TO NON-MECHANICAL STRAINS---
      DENN = DENN !- DEAG   !CORRECTION OF AGING STRAIN FROM NON-MECHANICAL STRAIN

C	INCREMENTAL STRESS DUE TO NON-MECHA. STRAIN
	DSIG = DENN*EMODP


C	WRITE(*,*) EMP,DENN
C	WRITE(*,*) IEL,SINC,PHIM
C	-----------------------
C	UPDATE WORKING ARRAY
C	-----------------------
	AF5(31) = ENP   !TOTAL NON-MECHANICAL STRAIN
	AF5(33) = EPSK  !TOTAL SHRNK STRAIN
	AF5(34) = EPCP  !TOTAL CREEP STRAIN

	AF5(18) = POISN
	AF5(20) = RHO     
	AF5(21) = FPC     
	AF5(22) = FPT    
	AF5(23) = SHRN    
	AF5(24) = RLXK(1) 
	AF5(25) = RLXK(2)  
	AF5(26) = EMODP       !YOUNG MODULUS RESPECTED TO CODE        !EMODP OR EMODC !!!!!!!!!! (EMODC WILL MAKE A TROUBLE IN CREEP&SHRINKAGE SECONDARY)

	K = 50
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = CCOF(I,J)		!CREEP COEF  (51-80) ALLOCATE 30 VALUE  
	ENDDO
	ENDDO
	K = 80
	DO I = 1,NRT
	K = K + 1
	AF5(K) =  GTJ(I)		!CREEP COEF  (81-95) ALLOCATE 15 VALUE  
	ENDDO	
	K = 100
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) =  ALFA(I,J)		!CREEP COEF  (101-130) ALLOCATE 30 VALUE  
	ENDDO
	ENDDO
	K = 130
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) =  APHI(I,J)		!CREEP COEF  (131-160) ALLOCATE 30 VALUE 
	ENDDO
	ENDDO


	DEALLOCATE(CCOP,CCOF,GTJ,APHI)
	DEALLOCATE(ELIT,ALFA)
C	-----------------------

C	-------------	
	CASE(1) !STEEL

	DSIG = 0.0D0 !NOTHING FOR STEEL

	ENDSELECT
C	--------------------------------------------


C	STRESS CONTRIBUTION
	SA(1)=SA(1) + DSIG*DA             !AXIAL FORCE
	SA(5)=SA(5) + TC*DSIG*DA          !MOMENT ABOUT S-AXIS (MY FOR XY PLANE)
	SA(6)=SA(6) - SC*DSIG*DA          !MOMENT ABOUT T-AXIS (MZ FOR XY PLANE)
	SA(7)=SA(7) + WP*DSIG*DA          !MOMENT WARPING

C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    

1000	CONTINUE

	

	DEALLOCATE(AF5)


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE FIBOPT22(NPM,NPI,NPIB,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL DATA AT THE END OF TIME STEP
C	IM = CREEP INTEGRATION OPTION  

	DIMENSION NPM(1),NPI(1),AAX(1)
	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))


	DO 1000 IFIB = 1,NPIB

	NPI(2) = IFIB
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART

C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE
	MATN    = INT(AF5(4))
	MATY    = INT(AF5(5))
C	----------------------------------------	

C	--------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MATY)

C	-------------
	CASE(0) !CONC	


	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)

      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG


	AF5(36) = 0.0D0  !DICP OLD ITERATIVE CREEP STRAIN

C	KL = 2 FOR OLD CEB (NORMAL AND RECOVERY PART)

	K = 23
	DO J = 1,KL
	K = K + 1
	AF5(K) =  0.0D0 !RLXK(J)		 
	ENDDO
	K = 100
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) =  0.0D0 !ALFA(I,J)		!CREEP COEF  (101-130) ALLOCATE 30 VALUE  
	ENDDO
	ENDDO
	K = 130
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) =  0.0D0 !APHI(I,J)		!CREEP COEF  (131-160) ALLOCATE 30 VALUE 
	ENDDO
	ENDDO
C	-----------------------

C	-------------	
	CASE(1) !STEEL

	DKAK= 0.0D0 !NOTHING FOR STEEL

	ENDSELECT
C	--------------------------------------------


C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    

1000	CONTINUE

	
	DEALLOCATE(AF5)


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FIBOPT23(NPM,NPI,NPIB,SA,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL DATA DURING ITERATION
C	IM = CREEP INTEGRATION OPTION  

	DIMENSION SA(1),NPM(1),NPI(1),AAX(1)
	DIMENSION RLXK(2)
	ALLOCATABLE GTJ(:),APHI(:,:)
	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))

	DO 1000 IFIB = 1,NPIB

	NPI(2) = IFIB
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART

C	----------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE

	DA      = AF5(1)  !FIBER AREA
	SC      = AF5(2)  !FIBER DISTANCE S
	TC      = AF5(3)  !FIBER DISTANCE T
	MATN    = INT(AF5(4))
	MATY    = INT(AF5(5))
	WP      = AF5(6)  !WARPING COEF
	TR      = AF5(7)  !TORSION COEF

	MATFL   = INT(AF5(16))
	EMOD    = AF5(17)
	RLXK(1) = AF5(24) 
	RLXK(2) = AF5(25) 

	ENP		= AF5(31)     !TOTAL NON-MECHANICAL STRAIN
	EPCP	= AF5(34)     !TOTAL CREEP STRAIN
	SINC    = AF5(35)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION

	DICP    = AF5(36)     !OLD ITERATIVE CREEP STRAIN
C	----------------------------------------	


C	--------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MATY)

C	-------------
	CASE(0) !CONC	

	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG
	ALLOCATE(GTJ(NRT),APHI(NRT,KL))


	K = 80
	DO I = 1,NRT
	K = K + 1
	GTJ(I)    = AF5(K)    !CREEP COEF  (81-95) ALLOCATE 15 VALUE  
	ENDDO	
	K = 130
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	APHI(I,J) = AF5(K)   !CREEP COEF  (131-160) ALLOCATE 30 VALUE 
	ENDDO
	ENDDO

	DSIG = 0.0D0
C	HERE FOR C-C
      LL = 1
      IF(KL.EQ.2) THEN
	IF(SINC*EPCP .LT. 0.D0) LL = 2
      ENDIF

	IF(IM.LT.2) GOTO 100

      LL = 1
      IF(KL.EQ.2) THEN
	CSTS = SINC * RLXK(1)
	IF(CSTS*EPCP .LT. 0.D0) LL = 2
      ENDIF
      DECP = RLXK(LL) * SINC  !INCREMENTAL CREEP STRAIN

C	COREECTION DUE TO CREEP CONVERGENCE ACCELERATION
	CR = 0.71
	DECP = (DECP + CR*DICP)/(1.0D0+CR)

C	RRR = MOD(DECP,1.0E-6)
C	DECP = DECP - RRR

      EPCP = EPCP + DECP      !UPDATE TOTAL CREEP STRAIN

C	UPDATE TOTAL NON-MECHANICAL STRAIN-----
	ENP  = ENP + DECP

C	INCREMENTAL STRESS DUE TO NON-MECHA. STRAIN
	DSIG = DECP*EMOD


100	CONTINUE


	DO I = 1,NRT
	GTJ(I) = GTJ(I) + APHI(I,LL)*SINC         !UPDATE g(tj)
	ENDDO
C	-----------------------
C	UPDATE WORKING ARRAY
C	-----------------------
	AF5(31) = ENP   !TOTAL NON-MECHANICAL STRAIN
	AF5(34) = EPCP  !TOTAL CREEP STRAIN

	K = 80
	DO I = 1,NRT
	K = K + 1
	AF5(K) =  GTJ(I)		!CREEP COEF  (81-95) ALLOCATE 15 VALUE  
	ENDDO

	AF5(36) = DECP !UPDATE ITERATIVE CREEP STRAIN

C	WRITE(*,222) IFIB,SINC
C	WRITE(110,222) IFIB,SINC,APHI(1,1),GTJ(1)
C222	FORMAT(I5,20E15.6)

	DEALLOCATE(GTJ,APHI)
C	-----------------------


C	-------------	
	CASE(1) !STEEL

	DSIG = 0.0D0 !NOTHING FOR STEEL

	ENDSELECT
C	--------------------------------------------


C	STRESS CONTRIBUTION
	SA(1)=SA(1) + DSIG*DA             !AXIAL FORCE
	SA(5)=SA(5) + TC*DSIG*DA          !MOMENT ABOUT S-AXIS (MY FOR XY PLANE)
	SA(6)=SA(6) - SC*DSIG*DA          !MOMENT ABOUT T-AXIS (MZ FOR XY PLANE)
	SA(7)=SA(7) + WP*DSIG*DA          !MOMENT WARPING


C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    

1000	CONTINUE

	
	DEALLOCATE(AF5)


	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
      SUBROUTINE GUSOPT21(TAU,TINC,NPM,NPI,ISEC,SA,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL AT START OF TIME STEP
C	IM = CREEP INTEGRATION OPTION  


	DIMENSION SA(1),NPM(1),NPI(1),AAX(1)
	DIMENSION RLXK(2)

	DIMENSION EMP(6),ENP(6),SCP(6),EPSK(6),EPCP(6)
	DIMENSION DESK(6),ECSK(6),DEAG(6)
	DIMENSION DENN(6),DECR(6),DSIG(6),SINC(6),DECR0(6)
	DIMENSION DMATX(6,6)

C	FOR FLOW VECTOR
	DIMENSION AVECT(6)

	ALLOCATABLE CCOP(:,:),CCOF(:,:),GTJ(:,:),APHI(:,:)
	ALLOCATABLE ELIT(:),ALFA(:,:)
	ALLOCATABLE AF5(:)

	
C     ----------------------------------------------------	
C     AF5  = WORKING ARRAY FOR FIBER/LAYER
C     ----------------------------------------------------	
	N5 = NPM(1)	
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
	CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN

	MATN = INT(FMAT)
	MTYP = INT(FTYP)	
	MATS = INT(FATS)
	ISKP = INT(SPAT)
	
C	SHRINKAGE FACTOR -----------------------------------
      !CALL RELFILL('SKP1',FSPK,1,ISKP,0)     !CALLING NUMBER OF SHRINKAGE POINT
	!NSPK = INT(FSPK)
	!NSPK = 1            
C     ----------------------------------------------------	
	
	!DO 1000 IGUS = 1,NLAY 
	
C	SHRINKAGE FACTOR ----------------------------------	
      FSRK = 1.0D0  !uniform shrinkage distribution
	IF(MTYP.EQ.0) THEN	
	FSRK = 1.0D0
	ENDIF
C     AF5-R(MECH & NON-MECH STRAI) ----------------------
	NPI(2) = 1    ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C	---------------------------------------------------
	NSTAT	= INT(AF5(1)) !MATERIAL STATE NUMBER
	EMOD    = AF5(2)
	POISN   = AF5(3)
	FPC     = AF5(4)
	FPT     = AF5(5)
	EMODP   = AF5(6)     !YOUNG MODULUS RESPECTED TO CODE
	FUNLO   = AF5(7)     !UNLOADING FLAG --> 0.0 FOR UNLOAD, 1.0 FOR LOADING 
      ! 6 STRAINS IN SOLID ELEMENT
	EMP(1:6) = AF5(11:16)     !TOTAL MECHANICAL STRAIN
	ENP(1:6) = AF5(17:22)     !TOTAL NON-MECHANICAL STRAIN
	SCP(1:6) = AF5(23:28)     !OLD STRESS
	EPSK(1:6)= AF5(29:34)     !TOTAL SHRNK STRAIN
	EPCP(1:6)= AF5(35:40)     !TOTAL CREEP STRAIN
	SINC(1:6)= AF5(41:46)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION

	RLXK(1:2)= AF5(47:48)     !R Factor for creep calculation 3.3.14b 3.3.19b	
C	---------------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	======================================================
	CASE(0) !CONC	
C	======================================================

	CALL INTFILL('%CPT',KM ,1,MATN,0)
	CALL INTFILL('%CPT',NAT,2,MATN,0)
	CALL INTFILL('%CPT',NTT,3,MATN,0)
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG
      
      
	ALLOCATE(CCOP(NRT,KL),CCOF(NRT,KL),GTJ(6,NRT),APHI(NRT,KL))
	ALLOCATE(ELIT(NRT),ALFA(NRT,KL))

	TAUO = TAU-TINC
	!ONLY FOR INITIAL SHIRNKAGE STRAIN SHRNO AT THE TIME OF ERECTION UPDATE EMODP HERE
	CALL TMPCNC(TAUO,MATN,EMODP,SHRNO,FPC,FPT,RHO,POISN,GG,CCOF) 
	CALL TMPCNC(TAU ,MATN,EMODC,SHRN ,FPC,FPT,RHO,POISN,GG,CCOF) 
	
C     ------------------------------------------------------------
C                     SHRINKAGE INCREMENTAL STRAIN
C     ------------------------------------------------------------
C	FIND SHRINKAGE STRAIN INCREMENT (r, s and t diraction)
	DESK(1:6)  = 0.0D0
	ECSK(1:3)  = SHRN*FSRK*(1.0D0/3.0D0)
      IF(TINC.NE.0.D0) THEN
C	EPSK = ECSK         !CHANGE TO NEXT TO ACCOUNT THE EFFECT OF INITIAL SHIRNKAGE STRAIN SHRNO AT THE TIME OF ERECTION
	EPSK(1:3) = SHRNO*FSRK*(1.0D0/3.0D0) 
	DESK(1:3) = ECSK(1:3) - EPSK(1:3)
	ENDIF	
	
C     ------------------------------------------------------------
C                       AGING INCREMENTAL STRAIN
C     ------------------------------------------------------------
C	FIND AGING STRAIN INCREMENT (r, s and t diraction)
	DEAG(1:6) = 0.0D0
      IF(TINC.NE.0.D0) THEN
	DEAG(1:3) = SCP(1:3)/EMODP - SCP(1:3)/EMODC
      ENDIF

C     ------------------------------------------------------------
C                      CREEP INCREMENTAL STRAIN
C     ------------------------------------------------------------
C	FIND CREEP RELATED VARIABLES
	CALL RELFILL('%FNC',GAMA,6,MATN,0)	
	GAMA = GAMA * 1.0D-1

      K = 55
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	CCOP(I,J) = AF5(K)   !CREEP COEF  (56-85) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	GTJ(I,J)  = AF5(K)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO

	K = 185
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	ALFA(I,J) = AF5(K)   !CREEP COEF  (186-215) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	APHI(I,J) = AF5(K)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO

	RLXK(1:2)   = 0.0D0
	
C	----------------------------
C	FIND CREEP STRAIN INCREMENT
C	----------------------------
      DO 13 J=1,NRT
      GAMA = GAMA * 1.0D+1
      YIOT = TINC / GAMA
      ELIM = DEXP(-YIOT)
      ELIT(J) = 1.0D+0 - ELIM

	LZE =0
	IF(ABS(YIOT).LT.1.0E-15) THEN
	YIOT = 0.0D0
	LZE = 1
	ENDIF
C
	PHIM = 0.0D0
      IF(LZE.EQ.0) 
	1 PHIM = 1.0D+0 - ( ELIT(J) / YIOT )         !(1-Phi)  Eq-3.3.14b
      IF(IM.EQ.3) BBIT = (YIOT / 2.0D+0) - PHIM   !(dt/2G)-(1-Phi)  Eq-3.3.19b
C
      DO 12 K=1,KL
      RLXF = CCOF(J,K) * PHIM                     !Ai(1-Phi) Eq-3.3.14b        ***Used CCOF = CCOP
	APHI(J,K) = CCOF(J,K)                       !A(tj-1)   Eq 3.3.11d  C-C   ***Used CCOF = CCOP
	IF(IM.LT.2) GOTO 12	
	APHI(J,K) = APHI(J,K)*(1.0-PHIM)            !A(tj-1)*phi Eq 3.3.14d C-L
	IF(IM.NE.3) GOTO 12							!L-L CREEP INTEGRATION METHOD
	APHI(J,K) = APHI(J,K) + ALFA(J,K)*PHIM      !A(tj-1)*phi + B(tj-1)*(1-phi) Eq 3.3.19e L-L
	ALFA(J,K) = 0.0D0
      IF(LZE.EQ.0) 
	1 ALFA(J,K) = (CCOF(J,K) - CCOP(J,K)) / YIOT !Beta=(dA)*G/dt   Eq-3.3.19d
      RLXF = RLXF + ALFA(J,K) * BBIT              !Add term   Beta*((dt/2G)-(1-Phi))  Eq-3.3.19b
12	RLXK(K) = RLXK(K) + RLXF                    !Sum over NRT  Eq-3.3.14b and  3.3.19b
C
13	CONTINUE
C	----------------------------
	DECR0(1:6) = 0.0D0
      DO 45 I=1,6
      DO 45 J=1,NRT
	DECR0(I) = DECR0(I) + GTJ(I,J) * ELIT(J)            !The 2nd term in Eq-3.3.11a  3.3.14a  3.3.19a  ,GTJ = g(tj-1)  
45	GTJ(I,J) = GTJ(I,J)*(1.0-ELIT(J))					!UPDATE g(tj) 


C     ------------------------------------------------------------
C	DRUCKER PRAGER YIELD SURFACE
C     ------------------------------------------------------------
	ALPI1 = 0.0D0
	CALL SOYSURF(SCP,SN2,AVECT,ALPI1)  !GET SN1 HERE

	SUMW  = SCP(1)*DECR0(1) + SCP(2)*DECR0(2) + SCP(3)*DECR0(3) + 
	1		  SCP(4)*DECR0(4) + SCP(5)*DECR0(5) + SCP(6)*DECR0(6)
	SUME  = SCP(1)*AVECT(1) + SCP(2)*AVECT(2) + SCP(3)*AVECT(3) + 
	1		  SCP(4)*AVECT(4) + SCP(5)*AVECT(5) + SCP(6)*AVECT(6)

	FLAMD = 0.0D0
	IF(SUME .NE.0.0D0) FLAMD = SUMW/SUME
	
	
C     ------------------------------------------------------------
C                   INCREMENTAL NON-MECHANICAL STRAIN
C     ------------------------------------------------------------
C	FOR UNLOADING
C	IF(FUNLO.EQ.0.0D0) THEN
	FLAMD = 1.0D0
	AVECT(1:6) = DECR0(1:6)
C	ENDIF

C     INCREMENTAL CREEP STRAIN (W/O R-CONTRIBUTION)
	DECR(1:6) = FLAMD*AVECT(1:6)
		
C     UPDATE CREEP STRAIN
      EPCP(1:6) = EPCP(1:6) + DECR(1:6)                   
C	---------------------------------------
C	SUM OF INCREMENTAL NON-MECHANICAL STRAIN
	DENN(1:6) = DESK(1:6) + DEAG(1:6) + DECR(1:6)

C	UPDATE TOTAL NON-MECHANICAL STRAIN
	ENP(1:6)  = ENP(1:6)  + DENN(1:6)
	
C-----ASSEMBLE EQUIVALENT ELEMENT LOADS DUE TO NON-MECHANICAL STRAINS---
      DENN(1:6) = DENN(1:6) !- DEAG(1:5)   !CORRECTION OF AGING STRAIN FROM NON-MECHANICAL STRAIN

C     ------------------------------------------------------------
C                   INCREMENTAL NON-MECHANICAL STRESS
C     ------------------------------------------------------------
C	INCREMENTAL STRESS DUE TO NON-MECHA. STRAIN	  
	CALL SOMODUC(DMATX,EMODP,POISN)	  
	DO I = 1,6
	DSIG(I) = 0.0D0
	DO J = 1,6
	DSIG(I) = DSIG(I) + DMATX(I,J)*DENN(J)
	ENDDO
	ENDDO
	DSIG = MATMUL(DMATX,DENN)
	
      
C	WRITE(*,*) EMP,DENN
C	WRITE(*,*) IEL,SINC,PHIM
C     ------------------------------------------------------------
C	UPDATE WORKING ARRAY
C     ------------------------------------------------------------
C     AF5-W(MATERIALS, NON-MECHANICAL STRAIN)----------------
	AF5(3) = POISN   
	AF5(4) = FPC     
	AF5(5) = FPT       
	AF5(6) = EMODP       !YOUNG MODULUS RESPECTED TO CODE       
      !EMODP OR EMODC !!!!!!!!!! (EMODC WILL MAKE A TROUBLE IN CREEP&SHRINKAGE SECONDARY)
	AF5(17:22) = ENP(1:6)   !TOTAL NON-MECHANICAL STRAIN
	AF5(29:34) = EPSK(1:6)  !TOTAL SHRNK STRAIN
	AF5(35:40) = EPCP(1:6)  !TOTAL CREEP STRAIN
	AF5(47:48) = RLXK(1:2)  !R Factor for creep calculation 3.3.14b 3.3.19b

	AF5(246:251) = AF5(246:251) + DSIG(1:6)  !STORE TOTAL NON-MECHANICAL STRESS
C	-------------------------------------------------------
	
      K = 55
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = CCOP(I,J)   !CREEP COEF  (56-85) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	AF5(K) = GTJ(I,J)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO

	K = 185
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = ALFA(I,J)   !CREEP COEF  (186-215) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = APHI(I,J)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO

	DEALLOCATE(CCOP,CCOF,GTJ,APHI)
	DEALLOCATE(ELIT,ALFA)
C	-------------------------------------------------------

C	======================================================
	CASE(1) !STEEL
C	======================================================
	DSIG(1:6) = 0.0D0 !NOTHING FOR STEEL

	ENDSELECT
C	-------------------------------------------------------

C	STRESS CONTRIBUTION
	DO I = 1,6	
	SA(I) = DSIG(I)
      ENDDO

C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT') 	

C1000	CONTINUE
	

	DEALLOCATE(AF5)
	


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
	SUBROUTINE GUSOPT22(NPM,NPI,ISEC,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL DATA AT THE END OF TIME STEP
C	IM = CREEP INTEGRATION OPTION  

	DIMENSION NPM(1),NPI(1),AAX(1)
	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
      CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN

	MATN = INT(FMAT)
	MTYP = INT(FTYP)	
	MATS = INT(FATS)
	ISKP = INT(SPAT)

	!DO 1000 IGUS = 1,NLAY
C	AF5-R(NON-MECH STRAIN, RIGIDITY)------------------------
	NPI(2) = 1  ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C	--------------------------------------------------------

C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	========================================================
	CASE(0) !CONC	
C	========================================================
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG

	AF5(50:55) = 0.0D0  !DICP OLD ITERATIVE CREEP STRAIN FOR CONVERGENCE ACCEL OF CREEP

C	KL = 2 FOR OLD CEB (NORMAL AND RECOVERY PART)

	AF5(47:48) =  0.0D0 !RLXK(J)		 

	K = 185
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = 0.0D0  !ALFA(I,J)   !CREEP COEF  (186-215) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = 0.0D0  !APHI(I,J)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO


C	========================================================	
	CASE(1) !STEEL
C	========================================================
	DKAK= 0.0D0 !NOTHING FOR STEEL

	ENDSELECT

C	STORE WORKING ARRAY
C	AF5-R(NON-MECH STRAIN, RIGIDITY)------------------------
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    
C	--------------------------------------------------------

C1000	CONTINUE
	
	DEALLOCATE(AF5)
	
	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
      SUBROUTINE GUSOPT23(NPM,NPI,ISEC,SA,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL DATA DURING ITERATION
C	IM = CREEP INTEGRATION OPTION  



	DIMENSION SA(6),NPM(1),NPI(1),AAX(1)
	DIMENSION RLXK(2)
	DIMENSION EMP(6),ENP(6),EPCP(6),SCP(6),SINC(6)
	DIMENSION DSIG(6),DMATX(6,6)
	DIMENSION DECP(6),DICP(6)

C	FOR FLOW VECTOR
	DIMENSION AVECT(6)

	ALLOCATABLE GTJ(:,:),APHI(:,:)
	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
	CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN	

	MATN = INT(FMAT)
	MTYP = INT(FTYP)	
	MATS = INT(FATS)
	ISKP = INT(SPAT)
	
	!DO 1000 ILAY = 1,NLAY
	
C     AF5-R(SECTION DATA)----------------------------------------
	NPI(2) = 1  ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C     -----------------------------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE
	EMOD     = AF5(2)
	POISN    = AF5(3)
	FUNLO    = AF5(7)         !UNLOADING FLAG  0.0 FOR UNLOAD  1.0 FOR LOAD	
	
	RLXK(1:2)= AF5(47:48)     !R Factor for creep calculation 3.3.14b 3.3.19b
	
	EMP(1:6) = AF5(11:16)     !TOTAL MECHANICAL STRAIN
	ENP(1:6) = AF5(17:22)     !TOTAL NON-MECHANICAL STRAIN
	SCP(1:6) = AF5(23:28)     !OLD STRESS
	EPCP(1:6)= AF5(35:40)     !TOTAL CREEP STRAIN
	SINC(1:6)= AF5(41:46)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION

	DICP(1:6)= AF5(50:55)     !OLD ITERATIVE CREEP STRAIN
C     -----------------------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	===========================================================
	CASE(0) !CONC	
C	===========================================================
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG
C     -----------------------------------------------------------      
	ALLOCATE(GTJ(6,NRT),APHI(NRT,KL))	
	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	GTJ(I,J)  = AF5(K)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO
	
	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	APHI(I,J) = AF5(K)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO
	
C     -----------------------------------------------------------
	DO 100 ISIG = 1,6

	DECP(ISIG) = 0.0D0
C	HERE FOR C-C
      LL = 1
      IF(KL.EQ.2) THEN
	IF(SINC(ISIG).LT. 0.D0) LL = 2
      ENDIF
	
C     -----------------------------------------------------------
C                       INCREMENTAL CREEP STRAIN
C     -----------------------------------------------------------
      IF(IM.LT.2) GOTO 90
      LL = 1
      IF(KL.EQ.2) THEN
	CSTS = SINC(ISIG) * RLXK(1)
	IF(CSTS .LT. 0.D0) LL = 2
      ENDIF
      
      DECP(ISIG) = RLXK(LL) * SINC(ISIG) !INCREMENTAL CREEP STRAIN

C	RRR = MOD(DECP,1.0E-6)
C	DECP = DECP - RRR
90	CONTINUE

	DO J = 1,NRT
	GTJ(ISIG,J) = GTJ(ISIG,J) + APHI(J,LL)*SINC(ISIG)    !UPDATE g(tj)
	ENDDO
	

100	CONTINUE

C	----------------------------------------------------------	
C     DISPLACEMENTS 
      !WRITE(IWRIT,*) 'SINC23 = ',SINC(1)         
C	----------------------------------------------------------	

C     -----------------------------------------------------------	
C	DRUCKER PRAGER YIELD SURFACE
C     -----------------------------------------------------------
	ALPI1 = 0.0D0
	CALL SOYSURF(SCP,SN2,AVECT,ALPI1)  !GET SN1 HERE

	SUMW  = SCP(1)*DECP(1)  + SCP(2)*DECP(2)  + SCP(3)*DECP(3) + 
	1		  SCP(4)*DECP(4)  + SCP(5)*DECP(5)  + SCP(6)*DECP(6)
	SUME  = SCP(1)*AVECT(1) + SCP(2)*AVECT(2) + SCP(3)*AVECT(3) + 
	1		  SCP(4)*AVECT(4) + SCP(5)*AVECT(5) + SCP(6)*AVECT(6)
	FLAMD = 0.0D0
	IF(SUME .NE.0.0D0) FLAMD = SUMW/SUME

C     -----------------------------------------------------------
C                   INCREMENTAL NON-MECHANICAL STRAIN
C     -----------------------------------------------------------
C	FOR UNLOADING
C	IF(FUNLO.EQ.0.0D0) THEN
	FLAMD = 1.0D0
	AVECT(1:6) = DECP(1:6)
C	ENDIF

	DECP(1:6) = FLAMD*AVECT(1:6)

	IF(IM.GE.2) THEN
C	COREECTION DUE TO CREEP CONVERGENCE ACCELERATION
	CR = 0.71 !0.71
	DECP(1:6) = (DECP(1:6) + CR*DICP(1:6))/(1.0D0+CR)
	ENDIF

C	UPDATE TOTAL CREEP STRAIN
      EPCP(1:6) = EPCP(1:6) + DECP(1:6)  
	    
C	UPDATE TOTAL NON-MECHANICAL STRAIN
	ENP(1:6)  = ENP(1:6)  + DECP(1:6)

C     -----------------------------------------------------------
C                 INCREMENTAL NON-MECHANICAL STRESS
C     -----------------------------------------------------------
C	INCREMENTAL STRESS DUE TO NON-MECHA. STRAIN
	CALL SOMODUC(DMATX,EMOD,POISN)
	DO I = 1,6
	DSIG(I) = 0.0D0
	DO J = 1,6
	DSIG(I) = DSIG(I) + DMATX(I,J)*DECP(J)
	ENDDO
	ENDDO
	DSIG = MATMUL(DMATX,DECP)


C     -----------------------------------------------------------
C	UPDATE WORKING ARRAY
C     -----------------------------------------------------------
	AF5(17:22) = ENP(1:6)   !TOTAL NON-MECHANICAL STRAIN
	AF5(35:40) = EPCP(1:6)  !TOTAL CREEP STRAIN
	
	AF5(246:251) = AF5(246:251) + DSIG(1:6) !STORE TOTAL NON-MECHANICAL STRESS
	
	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	AF5(K) = GTJ(I,J)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO

	AF5(50:55) = DECP(1:6)		    !UPDATE ITERATIVE CREEP STRAIN

C	WRITE(*,222) IFIB,SINC
C	WRITE(110,222) IFIB,SINC,APHI(1,1),GTJ(1)
C222	FORMAT(I5,20E15.6)

	DEALLOCATE(GTJ,APHI)


C	===========================================================
	CASE(1) !STEEL
C	===========================================================	
	DSIG = 0.0D0 !NOTHING FOR STEEL
	ENDSELECT
C	--------------------------------------------

C	STRESS CONTRIBUTION
	DO I = 1,6
	!SA(I) = SA(I) + DSIG(I)
	SA(I) = DSIG(I)
	ENDDO

C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')  
	
C1000	CONTINUE
	DEALLOCATE(AF5)

	
	RETURN
	END	
	
C	=======================================================================
C	=======================================================================
C	=======================================================================	



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
      SUBROUTINE GUCOPT21(TAU,TINC,NPM,NPI,ISEC,SA,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL AT START OF TIME STEP
C	IM = CREEP INTEGRATION OPTION  


	DIMENSION SA(1),NPM(1),NPI(1),AAX(1)
	DIMENSION RLXK(2)

	DIMENSION EMP(6),ENP(6),SCP(6),EPSK(6),EPCP(6)
	DIMENSION DESK(6),ECSK(6),DEAG(6)
	DIMENSION DENN(6),DECR(6),DSIG(6),SINC(6),DECR0(6)
	DIMENSION DMATX(6,6)

C	FOR FLOW VECTOR
	DIMENSION AVECT(6)

	ALLOCATABLE CCOP(:,:),CCOF(:,:),GTJ(:,:),APHI(:,:)
	ALLOCATABLE ELIT(:),ALFA(:,:)
	ALLOCATABLE AF5(:)

	
C     ----------------------------------------------------	
C     AF5  = WORKING ARRAY FOR FIBER/LAYER
C     ----------------------------------------------------	
	N5 = NPM(1)	
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
	CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN

	MATN = INT(FMAT)
	MTYP = INT(FTYP)	
	MATS = INT(FATS)
	ISKP = INT(SPAT)
	
C	SHRINKAGE FACTOR -----------------------------------
      !CALL RELFILL('SKP1',FSPK,1,ISKP,0)     !CALLING NUMBER OF SHRINKAGE POINT
	!NSPK = INT(FSPK)
	!NSPK = 1            
C     ----------------------------------------------------	
	
	!DO 1000 IGUS = 1,NLAY 
	
C	SHRINKAGE FACTOR ----------------------------------	
      FSRK = 1.0D0
	IF(MTYP.EQ.0) THEN	
	FSRK = 1.0D0
	ENDIF
C     AF5-R(MECH & NON-MECH STRAI) ----------------------
	NPI(2) = 1    ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C	---------------------------------------------------
	NSTAT	= INT(AF5(1)) !MATERIAL STATE NUMBER
	EMOD    = AF5(2)
	POISN   = AF5(3)
	FPC     = AF5(4)
	FPT     = AF5(5)
	EMODP   = AF5(6)     !YOUNG MODULUS RESPECTED TO CODE
	FUNLO   = AF5(7)     !UNLOADING FLAG --> 0.0 FOR UNLOAD, 1.0 FOR LOADING 
      ! 6 STRAINS IN SOLID ELEMENT
	EMP(1:6) = AF5(11:16)     !TOTAL MECHANICAL STRAIN
	ENP(1:6) = AF5(17:22)     !TOTAL NON-MECHANICAL STRAIN
	SCP(1:6) = AF5(23:28)     !OLD STRESS
	EPSK(1:6)= AF5(29:34)     !TOTAL SHRNK STRAIN
	EPCP(1:6)= AF5(35:40)     !TOTAL CREEP STRAIN
	SINC(1:6)= AF5(41:46)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION

	RLXK(1:2)= AF5(47:48)     !R Factor for creep calculation 3.3.14b 3.3.19b	
C	---------------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	======================================================
	CASE(0) !CONC	
C	======================================================

	CALL INTFILL('%CPT',KM ,1,MATN,0)
	CALL INTFILL('%CPT',NAT,2,MATN,0)
	CALL INTFILL('%CPT',NTT,3,MATN,0)
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG
      
      
	ALLOCATE(CCOP(NRT,KL),CCOF(NRT,KL),GTJ(6,NRT),APHI(NRT,KL))
	ALLOCATE(ELIT(NRT),ALFA(NRT,KL))

	TAUO = TAU-TINC
	!ONLY FOR INITIAL SHIRNKAGE STRAIN SHRNO AT THE TIME OF ERECTION UPDATE EMODP HERE
	CALL TMPCNC(TAUO,MATN,EMODP,SHRNO,FPC,FPT,RHO,POISN,GG,CCOF) 
	CALL TMPCNC(TAU ,MATN,EMODC,SHRN ,FPC,FPT,RHO,POISN,GG,CCOF) 
	
C     ------------------------------------------------------------
C                     SHRINKAGE INCREMENTAL STRAIN
C     ------------------------------------------------------------
C	FIND SHRINKAGE STRAIN INCREMENT (r, s and t diraction)
	DESK(1:6)  = 0.0D0
	ECSK(1:3)  = SHRN*FSRK*(1.0D0/3.0D0)
      IF(TINC.NE.0.D0) THEN
C	EPSK = ECSK         !CHANGE TO NEXT TO ACCOUNT THE EFFECT OF INITIAL SHIRNKAGE STRAIN SHRNO AT THE TIME OF ERECTION
	EPSK(1:3) = SHRNO*FSRK*(1.0D0/3.0D0) 
	DESK(1:3) = ECSK(1:3) - EPSK(1:3)
	ENDIF	
	
C     ------------------------------------------------------------
C                       AGING INCREMENTAL STRAIN
C     ------------------------------------------------------------
C	FIND AGING STRAIN INCREMENT (r, s and t diraction)
	DEAG(1:6) = 0.0D0
      IF(TINC.NE.0.D0) THEN
	DEAG(1:3) = SCP(1:3)/EMODP - SCP(1:3)/EMODC
      ENDIF

C     ------------------------------------------------------------
C                      CREEP INCREMENTAL STRAIN
C     ------------------------------------------------------------
C	FIND CREEP RELATED VARIABLES
	CALL RELFILL('%FNC',GAMA,6,MATN,0)	
	GAMA = GAMA * 1.0D-1

      K = 55
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	CCOP(I,J) = AF5(K)   !CREEP COEF  (56-85) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	GTJ(I,J)  = AF5(K)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO

	K = 185
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	ALFA(I,J) = AF5(K)   !CREEP COEF  (186-215) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	APHI(I,J) = AF5(K)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO

	RLXK(1:2)   = 0.0D0
	
C	----------------------------
C	FIND CREEP STRAIN INCREMENT
C	----------------------------
      DO 13 J=1,NRT
      GAMA = GAMA * 1.0D+1
      YIOT = TINC / GAMA
      ELIM = DEXP(-YIOT)
      ELIT(J) = 1.0D+0 - ELIM

	LZE =0
	IF(ABS(YIOT).LT.1.0E-15) THEN
	YIOT = 0.0D0
	LZE = 1
	ENDIF
C
	PHIM = 0.0D0
      IF(LZE.EQ.0) 
	1 PHIM = 1.0D+0 - ( ELIT(J) / YIOT )         !(1-Phi)  Eq-3.3.14b
      IF(IM.EQ.3) BBIT = (YIOT / 2.0D+0) - PHIM   !(dt/2G)-(1-Phi)  Eq-3.3.19b
C
      DO 12 K=1,KL
      RLXF = CCOF(J,K) * PHIM                     !Ai(1-Phi) Eq-3.3.14b        ***Used CCOF = CCOP
	APHI(J,K) = CCOF(J,K)                       !A(tj-1)   Eq 3.3.11d  C-C   ***Used CCOF = CCOP
	IF(IM.LT.2) GOTO 12	
	APHI(J,K) = APHI(J,K)*(1.0-PHIM)            !A(tj-1)*phi Eq 3.3.14d C-L
	IF(IM.NE.3) GOTO 12							!L-L CREEP INTEGRATION METHOD
	APHI(J,K) = APHI(J,K) + ALFA(J,K)*PHIM      !A(tj-1)*phi + B(tj-1)*(1-phi) Eq 3.3.19e L-L
	ALFA(J,K) = 0.0D0
      IF(LZE.EQ.0) 
	1 ALFA(J,K) = (CCOF(J,K) - CCOP(J,K)) / YIOT !Beta=(dA)*G/dt   Eq-3.3.19d
      RLXF = RLXF + ALFA(J,K) * BBIT              !Add term   Beta*((dt/2G)-(1-Phi))  Eq-3.3.19b
12	RLXK(K) = RLXK(K) + RLXF                    !Sum over NRT  Eq-3.3.14b and  3.3.19b
C
13	CONTINUE
C	----------------------------
	DECR0(1:6) = 0.0D0
      DO 45 I=1,6
      DO 45 J=1,NRT
	DECR0(I) = DECR0(I) + GTJ(I,J) * ELIT(J)            !The 2nd term in Eq-3.3.11a  3.3.14a  3.3.19a  ,GTJ = g(tj-1)  
45	GTJ(I,J) = GTJ(I,J)*(1.0-ELIT(J))					!UPDATE g(tj) 


C     ------------------------------------------------------------
C	DRUCKER PRAGER YIELD SURFACE
C     ------------------------------------------------------------
	ALPI1 = 0.0D0
	CALL SOYSURF(SCP,SN2,AVECT,ALPI1)  !GET SN1 HERE

	SUMW  = SCP(1)*DECR0(1) + SCP(2)*DECR0(2) + SCP(3)*DECR0(3) + 
	1		  SCP(4)*DECR0(4) + SCP(5)*DECR0(5) + SCP(6)*DECR0(6)
	SUME  = SCP(1)*AVECT(1) + SCP(2)*AVECT(2) + SCP(3)*AVECT(3) + 
	1		  SCP(4)*AVECT(4) + SCP(5)*AVECT(5) + SCP(6)*AVECT(6)

	FLAMD = 0.0D0
	IF(SUME .NE.0.0D0) FLAMD = SUMW/SUME
	
	
C     ------------------------------------------------------------
C                   INCREMENTAL NON-MECHANICAL STRAIN
C     ------------------------------------------------------------
C	FOR UNLOADING
C	IF(FUNLO.EQ.0.0D0) THEN
	FLAMD = 1.0D0
	AVECT(1:6) = DECR0(1:6)
C	ENDIF

C     INCREMENTAL CREEP STRAIN (W/O R-CONTRIBUTION)
	DECR(1:6) = FLAMD*AVECT(1:6)
		
C     UPDATE CREEP STRAIN
      EPCP(1:6) = EPCP(1:6) + DECR(1:6)                   
C	---------------------------------------
C	SUM OF INCREMENTAL NON-MECHANICAL STRAIN
	DENN(1:6) = DESK(1:6) + DEAG(1:6) + DECR(1:6)

C	UPDATE TOTAL NON-MECHANICAL STRAIN
	ENP(1:6)  = ENP(1:6)  + DENN(1:6)
	
C-----ASSEMBLE EQUIVALENT ELEMENT LOADS DUE TO NON-MECHANICAL STRAINS---
      DENN(1:6) = DENN(1:6) !- DEAG(1:5)   !CORRECTION OF AGING STRAIN FROM NON-MECHANICAL STRAIN

C     ------------------------------------------------------------
C                   INCREMENTAL NON-MECHANICAL STRESS
C     ------------------------------------------------------------
C	INCREMENTAL STRESS DUE TO NON-MECHA. STRAIN	
	CALL SCMODUC(DMATX,EMODP,POISN)  
	DO I = 1,6
	DSIG(I) = 0.0D0
	DO J = 1,6
	DSIG(I) = DSIG(I) + DMATX(I,J)*DENN(J)
	ENDDO
	ENDDO
	DSIG = MATMUL(DMATX,DENN)
	
        
	
C	WRITE(*,*) IEL,SINC,PHIM
C     ------------------------------------------------------------
C	UPDATE WORKING ARRAY
C     ------------------------------------------------------------
C     AF5-W(MATERIALS, NON-MECHANICAL STRAIN)----------------
	AF5(3) = POISN   
	AF5(4) = FPC     
	AF5(5) = FPT       
	AF5(6) = EMODP       !YOUNG MODULUS RESPECTED TO CODE       
      !EMODP OR EMODC !!!!!!!!!! (EMODC WILL MAKE A TROUBLE IN CREEP&SHRINKAGE SECONDARY)
	AF5(17:22) = ENP(1:6)   !TOTAL NON-MECHANICAL STRAIN
	AF5(29:34) = EPSK(1:6)  !TOTAL SHRNK STRAIN
	AF5(35:40) = EPCP(1:6)  !TOTAL CREEP STRAIN
	AF5(47:48) = RLXK(1:2)  !R Factor for creep calculation 3.3.14b 3.3.19b

	AF5(246:251) = AF5(246:251) + DSIG(1:6)  !STORE TOTAL NON-MECHANICAL STRESS
C	-------------------------------------------------------
	
      K = 55
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = CCOP(I,J)   !CREEP COEF  (56-85) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	AF5(K) = GTJ(I,J)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO

	K = 185
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = ALFA(I,J)   !CREEP COEF  (186-215) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = APHI(I,J)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO

	DEALLOCATE(CCOP,CCOF,GTJ,APHI)
	DEALLOCATE(ELIT,ALFA)
C	-------------------------------------------------------

C	======================================================
	CASE(1) !STEEL
C	======================================================
	DSIG(1:6) = 0.0D0 !NOTHING FOR STEEL

	ENDSELECT
C	-------------------------------------------------------

C	STRESS CONTRIBUTION
	DO I = 1,6	
	SA(I) = DSIG(I)
      ENDDO

C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT') 	

C1000	CONTINUE
	

	DEALLOCATE(AF5)
	


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
	SUBROUTINE GUCOPT22(NPM,NPI,ISEC,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL DATA AT THE END OF TIME STEP
C	IM = CREEP INTEGRATION OPTION  

	DIMENSION NPM(1),NPI(1),AAX(1)
	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
      CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN

	MATN = INT(FMAT)
	MTYP = INT(FTYP)	
	MATS = INT(FATS)
	ISKP = INT(SPAT)

	!DO 1000 IGUS = 1,NLAY
C	AF5-R(NON-MECH STRAIN, RIGIDITY)------------------------
	NPI(2) = 1  ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C	--------------------------------------------------------

C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	========================================================
	CASE(0) !CONC	
C	========================================================
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG

	AF5(50:55) = 0.0D0  !DICP OLD ITERATIVE CREEP STRAIN FOR CONVERGENCE ACCEL OF CREEP

C	KL = 2 FOR OLD CEB (NORMAL AND RECOVERY PART)

	AF5(47:48) =  0.0D0 !RLXK(J)		 

	K = 185
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = 0.0D0  !ALFA(I,J)   !CREEP COEF  (186-215) ALLOCATE 30 VALUES  
	ENDDO
	ENDDO

	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	AF5(K) = 0.0D0  !APHI(I,J)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO


C	========================================================	
	CASE(1) !STEEL
C	========================================================
	DKAK= 0.0D0 !NOTHING FOR STEEL

	ENDSELECT

C	STORE WORKING ARRAY
C	AF5-R(NON-MECH STRAIN, RIGIDITY)------------------------
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')    
C	--------------------------------------------------------

C1000	CONTINUE
	
	DEALLOCATE(AF5)
	
	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
      SUBROUTINE GUCOPT23(NPM,NPI,ISEC,SA,AAX,NDX5)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	LTMOPT(2)  WORKING ON NON-MECHANICAL DATA
C	 0 = NO UPDATE DATA
C	 1 = UPDATE DATA AT START OF TIME STEP
C	 2 = UPDATE DATA AT END OF TIME STEP
C	 3 = UPDATE DATA FOREACH ITERATION STEP

C	UPDATE NON-MECHANICAL DATA DURING ITERATION
C	IM = CREEP INTEGRATION OPTION  



	DIMENSION SA(6),NPM(1),NPI(1),AAX(1)
	DIMENSION RLXK(2)
	DIMENSION EMP(6),ENP(6),EPCP(6),SCP(6),SINC(6)
	DIMENSION DSIG(6),DMATX(6,6)
	DIMENSION DECP(6),DICP(6)

C	FOR FLOW VECTOR
	DIMENSION AVECT(6)

	ALLOCATABLE GTJ(:,:),APHI(:,:)
	ALLOCATABLE AF5(:)

	N5 = NPM(1)
	ALLOCATE(AF5(N5))

	CALL RELFILL('$SOC',FMAT,1,ISEC,0)           !MAT NO
	CALL RELFILL('$SOC',FTYP,2,ISEC,0)           !MAT TYPE  0=CONCRETE 1=STEEL
	CALL RELFILL('$SOC',SRTO,3,ISEC,0)           !STEEL RATIO	
	CALL RELFILL('$SOC',FATS,4,ISEC,0)           !MAT No. FOR STEEL
	CALL RELFILL('$SOC',SPAT,5,ISEC,0)           !SHRINKAGE PATTERN	

	MATN = INT(FMAT)
	MTYP = INT(FTYP)	
	MATS = INT(FATS)
	ISKP = INT(SPAT)
	
	!DO 1000 ILAY = 1,NLAY
	
C     AF5-R(SECTION DATA)----------------------------------------
	NPI(2) = 1  ! NEED TO BE MODIFIED
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')    !CALLING DATA OF THIS PART
C     -----------------------------------------------------------
C	CALLING DATA FROM WORKING ARRAY HERE
	EMOD     = AF5(2)
	POISN    = AF5(3)
	FUNLO    = AF5(7)         !UNLOADING FLAG  0.0 FOR UNLOAD  1.0 FOR LOAD	
	
	RLXK(1:2)= AF5(47:48)     !R Factor for creep calculation 3.3.14b 3.3.19b
	
	EMP(1:6) = AF5(11:16)     !TOTAL MECHANICAL STRAIN
	ENP(1:6) = AF5(17:22)     !TOTAL NON-MECHANICAL STRAIN
	SCP(1:6) = AF5(23:28)     !OLD STRESS
	EPCP(1:6)= AF5(35:40)     !TOTAL CREEP STRAIN
	SINC(1:6)= AF5(41:46)     !INCREMENTAL EFF. STRESS FOR CREEP CALCULATION

	DICP(1:6)= AF5(50:55)     !OLD ITERATIVE CREEP STRAIN
C     -----------------------------------------------------------
C	MATERIAL STATE & STRESS & STIFFNESS
	SELECTCASE(MTYP)

C	===========================================================
	CASE(0) !CONC	
C	===========================================================
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
      CALL INTFILL('%CPT',IM ,7,MATN,0)  !CREEP INTG FLAG
C     -----------------------------------------------------------      
	ALLOCATE(GTJ(6,NRT),APHI(NRT,KL))	
	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	GTJ(I,J)  = AF5(K)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO
	
	K = 215
	DO I = 1,NRT
	DO J = 1,KL
	K = K + 1
	APHI(I,J) = AF5(K)   !CREEP COEF  (216-245) ALLOCATE 30 VALUES 
	ENDDO
	ENDDO
	
C     -----------------------------------------------------------
	DO 100 ISIG = 1,6

	DECP(ISIG) = 0.0D0
C	HERE FOR C-C
      LL = 1
      IF(KL.EQ.2) THEN
	IF(SINC(ISIG).LT. 0.D0) LL = 2
      ENDIF
	
C     -----------------------------------------------------------
C                       INCREMENTAL CREEP STRAIN
C     -----------------------------------------------------------
      IF(IM.LT.2) GOTO 90
      LL = 1
      IF(KL.EQ.2) THEN
	CSTS = SINC(ISIG) * RLXK(1)
	IF(CSTS .LT. 0.D0) LL = 2
      ENDIF
      
      DECP(ISIG) = RLXK(LL) * SINC(ISIG) !INCREMENTAL CREEP STRAIN

C	RRR = MOD(DECP,1.0E-6)
C	DECP = DECP - RRR
90	CONTINUE

	DO J = 1,NRT
	GTJ(ISIG,J) = GTJ(ISIG,J) + APHI(J,LL)*SINC(ISIG)    !UPDATE g(tj)
	ENDDO
	

100	CONTINUE

C	----------------------------------------------------------	
C     DISPLACEMENTS 
      !WRITE(IWRIT,*) 'SINC23 = ',SINC(1)         
C	----------------------------------------------------------	

C     -----------------------------------------------------------	
C	DRUCKER PRAGER YIELD SURFACE
C     -----------------------------------------------------------
	ALPI1 = 0.0D0
	CALL SOYSURF(SCP,SN2,AVECT,ALPI1)  !GET SN1 HERE

	SUMW  = SCP(1)*DECP(1)  + SCP(2)*DECP(2)  + SCP(3)*DECP(3) + 
	1		  SCP(4)*DECP(4)  + SCP(5)*DECP(5)  + SCP(6)*DECP(6)
	SUME  = SCP(1)*AVECT(1) + SCP(2)*AVECT(2) + SCP(3)*AVECT(3) + 
	1		  SCP(4)*AVECT(4) + SCP(5)*AVECT(5) + SCP(6)*AVECT(6)
	FLAMD = 0.0D0
	IF(SUME .NE.0.0D0) FLAMD = SUMW/SUME

C     -----------------------------------------------------------
C                   INCREMENTAL NON-MECHANICAL STRAIN
C     -----------------------------------------------------------
C	FOR UNLOADING
C	IF(FUNLO.EQ.0.0D0) THEN
	FLAMD = 1.0D0
	AVECT(1:6) = DECP(1:6)
C	ENDIF

	DECP(1:6) = FLAMD*AVECT(1:6)

	IF(IM.GE.2) THEN
C	COREECTION DUE TO CREEP CONVERGENCE ACCELERATION
	CR = 0.71 !0.71
	DECP(1:6) = (DECP(1:6) + CR*DICP(1:6))/(1.0D0+CR)
	ENDIF

C	UPDATE TOTAL CREEP STRAIN
      EPCP(1:6) = EPCP(1:6) + DECP(1:6)  
	    
C	UPDATE TOTAL NON-MECHANICAL STRAIN
	ENP(1:6)  = ENP(1:6)  + DECP(1:6)

C     -----------------------------------------------------------
C                 INCREMENTAL NON-MECHANICAL STRESS
C     -----------------------------------------------------------
C	INCREMENTAL STRESS DUE TO NON-MECHA. STRAIN
	CALL SCMODUC(DMATX,EMODP,POISN)
	DO I = 1,6
	DSIG(I) = 0.0D0
	DO J = 1,6
	DSIG(I) = DSIG(I) + DMATX(I,J)*DECP(J)
	ENDDO
	ENDDO
	DSIG = MATMUL(DMATX,DECP)


C     -----------------------------------------------------------
C	UPDATE WORKING ARRAY
C     -----------------------------------------------------------
	AF5(17:22) = ENP(1:6)   !TOTAL NON-MECHANICAL STRAIN
	AF5(35:40) = EPCP(1:6)  !TOTAL CREEP STRAIN
	
	AF5(246:251) = AF5(246:251) + DSIG(1:6) !STORE TOTAL NON-MECHANICAL STRESS
	
	K = 85
	DO I = 1,6
	DO J = 1,NRT
	K = K + 1
	AF5(K) = GTJ(I,J)   !CREEP COEF  (86-185) ALLOCATE 100 VALUES  
	ENDDO
	ENDDO

	AF5(50:55) = DECP(1:6)		    !UPDATE ITERATIVE CREEP STRAIN

C	WRITE(*,222) IFIB,SINC
C	WRITE(110,222) IFIB,SINC,APHI(1,1),GTJ(1)
C222	FORMAT(I5,20E15.6)

	DEALLOCATE(GTJ,APHI)


C	===========================================================
	CASE(1) !STEEL
C	===========================================================	
	DSIG = 0.0D0 !NOTHING FOR STEEL
	ENDSELECT
C	--------------------------------------------

C	STRESS CONTRIBUTION
	DO I = 1,6
	!SA(I) = SA(I) + DSIG(I)
	SA(I) = DSIG(I)
	ENDDO

C	STORE WORKING ARRAY
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')  
	
C1000	CONTINUE
	DEALLOCATE(AF5)

	
	RETURN
	END	
	
C	=======================================================================
C	=======================================================================
C	=======================================================================	