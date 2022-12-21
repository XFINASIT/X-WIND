C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE GRLSTS
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	CALLING NEG
	CALL INTFILL('%NUB',NEG,1,5,0)


C	----------------------------------------
C	LOOP OVER ELEMENT GROUP 
C	----------------------------------------	
	DO 5000 IEG = 1,NEG


	CALL ELLSTS(IEG)

C	----------------------------------------
5000	CONTINUE
C	----------------------------------------



	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ELLSTS(IEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	ELEMENT FORCE AND STIFFNESS
	DIMENSION NPM(10),NPI(10),FOCT(12)
	ALLOCATABLE IF1(:),IF2(:),IF3(:),AF6(:),AF8(:),AF9(:)
	ALLOCATABLE COORD(:,:)
	ALLOCATABLE LMEQ(:),LREQ(:),RE(:),RL(:),SE(:),STLD(:),EXDI(:)
	ALLOCATABLE IAX(:),AAX(:)



      LOPT = 10            !LAST UPDATE FOR TOTAL COMPONENT
	IF(LOPT.LE.0) RETURN

	NNN= LOPT - 1


      CALL LOCATN('-GSF',KGSF,NUM,NWK ,2) 
	CALL LOCATN('-LMA',KLMA,NUM,NEQ1,1)  


	CALL INTFILL('%NUB',NSN,1,1,0)
	CALL INTFILL('%NUB',NSF,1,2,0)
	CALL INTFILL('%NUB',NSC,1,4,0)


	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	CALL INTFILL('%GRP',NNM  ,8 ,IEG,0) !
	
	CALL INTFILL('%GRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('%GRF',N2   ,2 ,IEG,0) !
	CALL INTFILL('%GRF',N3   ,3 ,IEG,0) !
	CALL INTFILL('%GRF',N6   ,6 ,IEG,0) !
	CALL INTFILL('%GRF',N8   ,8 ,IEG,0) !
	CALL INTFILL('%GRF',N9   ,9 ,IEG,0) !

	CALL INTFILL('%GRF',NDX1 ,11,IEG,0) !
	CALL INTFILL('%GRF',NDX2 ,12,IEG,0) !
	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX6 ,16,IEG,0) !
	CALL INTFILL('%GRF',NDX8 ,18,IEG,0) !
	CALL INTFILL('%GRF',NDX9 ,19,IEG,0) !
	
	CALL INTFILL('%GRF',NFLI ,21,IEG,0) !
	CALL INTFILL('%GRF',NFLR ,22,IEG,0) !
	
	CALL INTFILL('%GRF',LASI ,23,IEG,0) !
	CALL INTFILL('%GRF',LASR ,24,IEG,0) !


	NEK = (NEF*NEF+NEF)/2

	ALLOCATE(LMEQ(NEF),LREQ(NEF),RE(NEF),RL(NEF),SE(NEK))
	ALLOCATE(COORD(3,NNM))
	ALLOCATE(IF1(N1),IF2(N2),IF3(N3),AF6(N6),AF8(N8),AF9(N9))
	ALLOCATE(STLD(NEF),EXDI(NEF))
	ALLOCATE(IAX(LASI),AAX(LASR))

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
      
	SE(1:NEK) = 0.0D0
	RE(1:NEF) = 0.0D0
	RL(1:NEF) = 0.0D0
	EXDI(1:NEF) = 0.0D0
	STLD(1:NEF) = 0.0D0

	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ') !ELEMENT DATA
	IXIST = IF3(4)     !TOTAL EXISTING FLAG 
	IF(IXIST.EQ.0) GOTO 5000
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

	NPM(1) = N2
	NPI(1) = N2
	CALL CALLINTG(IF2,IAX,NDX2,NPM,NPI,1,'READ')
	DO IEF = 1,NEF
	LREQ(IEF) = IF2(IEF)       !ORIGINAL EQUA
	LMEQ(IEF) = IF2(IEF+NEF)   !IF2(NEF+1:NEF+NEF) = CURRENT EQAUATION
	ENDDO


C	LFOC  0 = NO COMPUTE OF ELM FORCE   1 = COMPUTE ELM FORCE
C	------------------------------------
	LSTF = 0
	LFOC = 2
	CALL MDUOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)
C	--------------------------
C	BACK UP NODAL FORCE FROM MECHA
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ') 
	NPM(1) = N8
	NPI(1) = N8
	CALL CALLREAL(AF8,AAX,NDX8,NPM,NPI,1,'READ') 
	NPM(1) = N9
	NPI(1) = N9
	CALL CALLREAL(AF9,AAX,NDX9,NPM,NPI,1,'READ')         
	NUM1 = 1   + NEF + NEF + NEF
	NUM2 = NEF + NEF + NEF + NEF
	RL(1:NEF) = AF6(NUM1:NUM2)
	NUM1 = 1   + NEF*NNN
	NUM2 = NEF + NEF*NNN
	AF8(NUM1:NUM2) = RL(1:NEF)     !FIXEND 
	NPM(1) = N8
	NPI(1) = N8
	CALL CALLREAL(AF8,AAX,NDX8,NPM,NPI,1,'WRIT') 
	NUM1 = 1   + NEF 
	NUM2 = NEF + NEF 
	RL(1:NEF) = AF6(NUM1:NUM2)
	NUM1 = 1   + NEF*NNN
	NUM2 = NEF + NEF*NNN
	AF9(NUM1:NUM2) = RL(1:NEF)     !INTERNAL FORCE
	NPM(1) = N9
	NPI(1) = N9
	CALL CALLREAL(AF9,AAX,NDX9,NPM,NPI,1,'WRIT')         
C	--------------------------
C	------------------------------------



C	--------------------------
	NPM(1) = N8
	NPI(1) = N8
	CALL CALLREAL(AF8,AAX,NDX8,NPM,NPI,1,'READ') 
	NPM(1) = N9
	NPI(1) = N9
	CALL CALLREAL(AF9,AAX,NDX9,NPM,NPI,1,'READ')  
	NUM1 = 1   + NEF*NNN
	NUM2 = NEF + NEF*NNN  
	NUM3 = 1   + NEF*(NNN+10)
	NUM4 = NEF + NEF*(NNN+10)
	NUM5 = 1   + NEF*(NNN+20)
	NUM6 = NEF + NEF*(NNN+20)
	AF8(NUM5:NUM6) = AF8(NUM1:NUM2) - AF8(NUM3:NUM4)  !GET INC
	AF9(NUM5:NUM6) = AF9(NUM1:NUM2) - AF9(NUM3:NUM4)  !GET INC
	AF8(NUM3:NUM4) = AF8(NUM1:NUM2)  !UPDATE
	AF9(NUM3:NUM4) = AF9(NUM1:NUM2)  !UPDATE
	NPM(1) = N8
	NPI(1) = N8
	CALL CALLREAL(AF8,AAX,NDX8,NPM,NPI,1,'WRIT') 
	NPM(1) = N9
	NPI(1) = N9
	CALL CALLREAL(AF9,AAX,NDX9,NPM,NPI,1,'WRIT')    
C	--------------------------
	
	IF(ITYPE.EQ.5) CALL CSASTSF(IEG,IEL,NNN,IAX,AAX)  !FRAME STRESS
	
      WRITE(NFLI,REC=IEL) IAX
      WRITE(NFLR,REC=IEL) AAX
      
C	------------------------------------
5000	CONTINUE
C	------------------------------------
      
	DEALLOCATE(IAX,AAX)
	DEALLOCATE(LMEQ,LREQ,RE,RL,SE)
	DEALLOCATE(COORD)
	DEALLOCATE(IF1,IF2,IF3,AF6,AF8,AF9)
	DEALLOCATE(STLD,EXDI)


	RETURN
	END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CSASTSF(IEG,IEL,NNN,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	INITIALIZE WORKING ARRAY
	DIMENSION NPM(10),NPI(10),PROPA(20),IAX(1),AAX(1)
C	STRESS CALCULATION
	DIMENSION STC(4,12),STCON(4,4),SSCON(4,4)
	DIMENSION ABC(4),PMT(4),STS(12)

	ALLOCATABLE IF3(:),AF4(:),AF5(:),AF7(:),AF8(:),AF9(:)
	ALLOCATABLE RL(:)
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
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !
	CALL INTFILL('%GRF',N8   ,8 ,IEG,0) !
	CALL INTFILL('%GRF',N9   ,9 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
	CALL INTFILL('%GRF',NDX8 ,18,IEG,0) !
	CALL INTFILL('%GRF',NDX9 ,19,IEG,0) !
	
C	----------------------------------------
	ALLOCATE(IF3(N3),AF4(N4),AF5(N5),AF7(N7),AF8(N8),AF9(N9))
	ALLOCATE(RL(NEF))

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
	DO 1000 NIJ = 1,2  !LOOP OVER 2 ENDS
      IGR = 1								!GET SECTION PROPERTIES AT FIRST GAUSS
	IF(NIJ.EQ.2) IGR = NGR				!GET SECTION PROPERTIES AT LAST  GAUSS

	NPM(1:2) = [N7,NGR]
	NPI(1:2) = [N7,IGR]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')   !CALLING RIGIDITY AND STRESS (MCHA & NON-MCHA)
	N7_R = 1           !POINTER OF RIGIDITY

	DO I = 1,20
	PROPA(I) = AF7(I+N7_R-1) !RIGIDITY
	ENDDO
	A    = PROPA(1)
	QS   = PROPA(2)
	QT   = PROPA(3)
	SIS  = PROPA(4)
	SIT  = PROPA(5)
	SIST = PROPA(6)
	FWWP = PROPA(9)
	FWSP = PROPA(10)
	FWTP = PROPA(11)
	FWP  = PROPA(12)

C     ----------------------------------------------------------
	DO 300 IPRT = 1,NPRT
C     ----------------------------------------------------------

	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'READ')
	NPIB    = INT(AF4(2))   !NUM FIBER FOR THIS PART

	SELECTCASE(LSTY)

	CASE(0) !READY PROP
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   1,IGR ,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')      !CALLING DATA OF THIS PART
	EMOD  = AF5(26)  

	CASE(1) !FIBER PROP
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,   1,IGR ,IPRT]                    !SET NFIB TO BE 1  (USE YOUNG MODULUS At FIRST FIBER OF PART)
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'READ')      !CALLING DATA OF THIS PART
	EMOD  = AF5(26)  

	ENDSELECT

C	-----------------------------------------
	DO IST = 1,12
	STC(1,IST) = 1.0
	INDEX = 4*(NIJ-1)+4*2*(IST-1)+4*2*12*(IPRT-1)
	CALL RELFILL('@FST',FIS       ,1+INDEX,ISEC,0)
	CALL RELFILL('@FST',STC(2,IST),2+INDEX,ISEC,0)
	CALL RELFILL('@FST',STC(3,IST),3+INDEX,ISEC,0)
	CALL RELFILL('@FST',STC(4,IST),4+INDEX,ISEC,0)
	ENDDO

C	A
	STCON(1,1) = A
C	Qs
	STCON(2,1) = QT
C	Qt
	STCON(3,1) = QS
C	Qw
	STCON(4,1) = FWP

C	Qt
	STCON(1,2) = QT
C	It
	STCON(2,2) = SIT
C	Ist
	STCON(3,2) = SIST
C	Isw
	STCON(4,2) = FWSP

C	Qs
	STCON(1,3) = QS
C	Ist
	STCON(2,3) = SIST
C	Is
	STCON(3,3) = SIS
C	Itw
	STCON(4,3) = FWTP

C	Qw
	STCON(1,4) = FWP
C	Isw
	STCON(2,4) = FWSP
C	Itw
	STCON(3,4) = FWTP
C	Iw
	STCON(4,4) = FWWP
C	--------------------------------------

	NPM(1) = N8
	NPI(1) = N8
	CALL CALLREAL(AF8,AAX,NDX8,NPM,NPI,1,'READ')  
	NPM(1) = N9
	NPI(1) = N9
	CALL CALLREAL(AF9,AAX,NDX9,NPM,NPI,1,'READ') 
	NUM5 = 1   + NEF*(NNN+20)  !INC END FORCE
	NUM6 = NEF + NEF*(NNN+20)  !INC END FORCE
	RL(1:NEF) = AF9(NUM5:NUM6) -  AF8(NUM5:NUM6)  !TOTAL END FORCE

	PMT(1) = RL(1+7*(NIJ-1))
	PMT(2) =-RL(6+7*(NIJ-1))
	PMT(3) = RL(5+7*(NIJ-1))
	PMT(4) = RL(7+7*(NIJ-1))


c      IF(IEL.EQ.40.AND.NNN.EQ.9) THEN
c      WRITE(*,*) NIJ
c      WRITE(*,*) A/EMOD,QT/EMOD,SIT/EMOD

c      WRITE(*,*) QT/A,3.5+QT/A
c      WRITE(*,*) PMT(1),PMT(2)
      
c      WRITE(*,*) PMT(2)+PMT(1)*QT/A
c      PAUSE
c      ENDIF
C	--------------------------------------
	STCON(1:4,1:4) = STCON(1:4,1:4)/EMOD
	IF(ABS(STCON(4,4)).LT.1.0E-6) THEN
	STCON(4,1:4) = 0.0D0
	STCON(1:4,4) = 0.0D0
	STCON(4,4) = 1.0D0
	PMT(4) = 0.0D0
	ENDIF
C	--------------------------------------
	N = 4
	CALL INVMATF(STCON,SSCON,N,IB)
C	--------------------------------------

	ABC = MATMUL(SSCON,PMT)
	STS = MATMUL(ABC,STC)
	MXSTRP = 12
	DO IST = 1,12
	INDEX = 4*(NIJ-1)+4*2*(IST-1)+4*2*12*(IPRT-1)
	CALL RELFILL('@FST',FIS,1+INDEX,ISEC,0)
	IF(FIS.EQ.0.0D0) STS(IST) = 0.0D0
	INDEX = IST+MXSTRP*(NIJ-1)+MXSTRP*2*(IPRT-1)+MXSTRP*2*MPRT*NNN
	CALL RELFILL('+FCS',STS(IST),INDEX,IEL,2)   !ACCUMULATE STRESS AT 2 END I-J
	ENDDO

C	--------------------------------------


C	-----------------------------------------

C     ----------------------------------------------------------
300	CONTINUE
C     ----------------------------------------------------------


1000	CONTINUE


	DEALLOCATE(IF3,AF4,AF5,AF7,AF8,AF9)
	DEALLOCATE(RL)


	RETURN

	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================