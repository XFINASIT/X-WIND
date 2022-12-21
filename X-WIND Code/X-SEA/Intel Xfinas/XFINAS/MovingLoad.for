C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE LANREAD(LEST)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	=============================================================
      DIMENSION LEST(1),ECLAN(3000),FILAN(3000),LANEL(3000),LANNS(3000)
	DIMENSION EGDATA(100000),PLNLG(10),LANEN(2,3000),SPANF(101)
	DIMENSION ECCLN(3000)
C	RETURN


	NLANE = 0

	READ(ITI,*)
	READ(ITI,*) NLANE

	IF(NLANE.LE.0) RETURN


	CALL DEFNINT('LNMX',KFILE,1,200)   !MOVING LOAD FILE      
	CALL INTZERO('LNMX')
	
	CALL DEFNINT('LANE',KLANE,1,20 )   !LANE DATA    
	CALL INTZERO('LANE')


	CALL INTFILL('LANE',NLANE,1 ,1,1)  !NUMBER OF LANE

	NFIL = 600
	CALL INTFILL('LNMX',NFIL ,1 ,1 ,1)  !LANE DATA									NLANE
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWS')  !SEQUENTIAL

	NFIL = 650
	CALL INTFILL('LNMX',NFIL ,1 ,2 ,1)  !LANE FIXEND FORCE DATA						NLANE
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWS')  !SEQUENTIAL

	NFIL = 700
	CALL INTFILL('LNMX',NFIL ,1 ,3 ,1)  !RESPONSE DUE TO UNIT LOAD AT EACH STATION	NPL*NLANE
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	NFIL = 701
	CALL INTFILL('LNMX',NFIL ,1 ,4 ,1)  !MAX RESPONSE DUE TO VEHICLE				NVEH*NLANE
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	NFIL = 702
	CALL INTFILL('LNMX',NFIL ,1 ,5 ,1)  !MAX RESPONSE DUE TO UNIT UNITFORM LOAD		NLANE   (INTEGRAT FROM UNIT LAD)
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	NFIL = 703
	CALL INTFILL('LNMX',NFIL ,1 ,6 ,1)  !MAX RESPONSE DUE TO UNIT LOAD		        SPAN*NLANE (AT EACH SPAN)
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	NFIL = 704
	CALL INTFILL('LNMX',NFIL ,1 ,7 ,1)  !MAX RESPONSE DUE TO VEHICLE CLASS			NVEHC*NLANE
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	NFIL = 705
	CALL INTFILL('LNMX',NFIL ,1 ,8 ,1)  !VAHICLE CLASS DATA     
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWS')  !SEQUENTIAL
    
	NFIL = 706
	CALL INTFILL('LNMX',NFIL ,1 ,9 ,1)  !DUMMY FILE 1       
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	NFIL = 707
	CALL INTFILL('LNMX',NFIL ,1 ,10,1)  !DUMMY FILE 2        
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	NFIL = 708
	CALL INTFILL('LNMX',NFIL ,1 ,11,1)  !DUMMY FILE 3    
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS
                 
	NFIL = 709
	CALL INTFILL('LNMX',NFIL ,1 ,12,1)  !DUMMY FILE 4    
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS


	NEDATA    = 65
	CALL INTFILL('LANE',NEDATA,1 ,2,1)  !DATA FOR LANE ELEMENT

	NSPANX = 0  !MAXIMUM NUMBER OF SPAN

	READ(ITI,*)
	DO 1000 ILANE = 1,NLANE

	CALL INTFILL('LNMX',NFIL1,1 ,1 ,0) 
	LP1 = NFIL1 + ILANE
	OPEN (LP1,FORM='UNFORMATTED',STATUS='SCRATCH')

	CALL INTFILL('LNMX',NFIL2,1 ,2 ,0) 
	LP2 = NFIL2 + ILANE
	OPEN (LP2,FORM='UNFORMATTED',STATUS='SCRATCH')


	READ(ITI,*) IL,NLNS,IDREC,NPL,IOPT,NSPAN,WTHLN		!IDREC  1= NORMAL  2= REVERSE  3= BOTH---NPL NUMBER OF INFLUENCE POINT ON THE LANE
	IF(NPL.LE.2) NPL = 3
	READ(ITI,*) (SPANF(J),J=1,NSPAN+1)					!SPAN/LEN RATIO
	NEM = 0
	DO ILNS = 1,NLNS 
	READ(ITI,*)  NE
	IF(IOPT.EQ.0) READ(ITI,*) (LANEL(NEM+J),J=1,NE)    !ELEMENT NUMBER OF THE LANE	
	IF(IOPT.EQ.1) THEN
	READ(ITI,*) (LANEL(J),J=1,NE+1)					   !NODE NUMBER OF THE LANE	
	DO J = 1,NE
	LANEN(1,NEM+J) = LANEL(J)
	LANEN(2,NEM+J) = LANEL(J+1)
	ENDDO
	ENDIF

	READ(ITI,*) (ECLAN(NEM+J),J=1,NE)  !SUB LANE FACTOR
	READ(ITI,*) (ECCLN(NEM+J),J=1,NE)  !SUB LANE FACTOR
	READ(ITI,*) (FILAN(NEM+J),J=1,NE)  !IMPACT FACTOR
	DO J = 1,NE
	LANNS(NEM+J) = ILNS                !SUB LANE NUMBER
	ENDDO
	NEM = NEM + NE	
	ENDDO

	PLNLG = 0.0D0
	CALL LANDATE1(LEST,NEM,IDREC,NPL,LANEL,ECLAN,ECCLN,FILAN,
	1			  LANNS,LP2,EGDATA,NEDATA,PLNLG,LANEN,IOPT)

	WRITE(LP1)  IL,NLNS,NEM,IDREC,NPL,IOPT,NSPAN,WTHLN
	WRITE(LP1)  (SPANF(J),J=1,NSPAN+1)
	WRITE(LP1)  (PLNLG(J),J=1,NLNS)


	CALL LANLFP(NEM,LP1,EGDATA,NEDATA)

	IF(NSPAN.GT.NSPANX) NSPANX = NSPAN

1000	CONTINUE


	CALL INTFILL('LANE',NSPANX,1 ,3,1)  !MAX NUMBER OF SPAN

	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE LANDATE1(LEST,NE,IDREC,NPL,LANEL,ECLAN,ECCLN,FILAN,
	1			        LANNS,LP2,EGDATA,NEDATA,PLNLG,LANEN,IOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN
	COMMON /GiDEle/ LGID
	COMMON /REACT/ LRC,LRCT,MFQ,LRID
	COMMON A(9000000),IA(9000000)
C
      DIMENSION LEST(1),EGDATA(1),LANEL(1),ECLAN(1),FILAN(1),LANNS(1)
	DIMENSION PLNLG(1),LANEN(2,1),ECCLN(1)

C
	IF(IOPT.EQ.1) GOTO 950
C	----------------------------------------------------------
      DO 900  IEG=1,NEG
      NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     

      CALL MOVLEV (2)
	
	KEG = IEG
	CALL LANDATE2(IA(LCN),IA(LLM),A(LCO),IA(LGS),
	1			  A(LGP),IA(LGID),IA(LHG),IA(LHGN),
	2			  NE,LANEL,ECLAN,FILAN,LANNS,
     3			  EGDATA,NEDATA,PLNLG,ECCLN)
	
      REWIND NELEMI
      REWIND NELEMA
C      WRITE (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      WRITE (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      WRITE (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      WRITE (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     
 900  CONTINUE
	
	GOTO 1000

C	----------------------------------------------------------
950	CONTINUE

	CALL LANDATN(IA(LID),IA(LRID),A(LXY),NSF,NSN,
	2			 NE,LANEL,ECLAN,FILAN,LANNS,
     3			 EGDATA,NEDATA,PLNLG,LANEN,ECCLN)	
	
C	----------------------------------------------------------
1000	CONTINUE
	CALL LANLOADF(NE,IDREC,NPL,IA(LDS),A(LDC),NLS,LP2,PLNLG,
	1			  EGDATA,NEDATA,IOPT)
C
      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANDATE2(LMN,LM,XYZ,IGSET,PROPG,IGIDM,LRPIN,IHSET,
	1					NE,LANEL,ECLAN,FILAN,LANNS,
	2					EGDATA,NEDATA,PLNLG,ECCLN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	=============================================================
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	=============================================================
	DIMENSION  LMN(NEF,1),LM(NEF,1),IGIDM(1),XYZ(6,1)
	DIMENSION  IGSET(1),PROPG(NGP,1),LRPIN(14,1),IHSET(1),IPIN(14)
	DIMENSION  VR(3),RG(24)
	DIMENSION  LANEL(1),ECLAN(1),FILAN(1),LANNS(1)
	DIMENSION  EGDATA(NEDATA,1)
	DIMENSION  PLNLG(1),ECCLN(1)


	ITEST = 0 !INITIALIZE LANE DIRECTION FLAG

	DO IE = 1,NE

	IGM = LANEL(IE)
	FAL = ECLAN(IE)
	FIM = FILAN(IE)
	NLN = LANNS(IE)
	ECC = ECCLN(IE)
	MLE = IGM
	CALL GIDORDER(IGIDM,NELE,MLE)
	IF(MLE.EQ.0) GOTO 100

	CALL XFSECTION(KEG,MLE,1)

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !ROTATION ANGLE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)
	NN1 = LMN(1,MLE)
	NN2 = LMN(2,MLE)
	PLNLG2 = PLNLG(NLN) + ELN
C	----------------------
	EGDATA(1,IE) = FLOAT(IGM)
	EGDATA(2,IE) = FLOAT(KEG)
	EGDATA(3,IE) = FLOAT(MLE)
	EGDATA(4,IE) = FLOAT(ITEST)
	EGDATA(5,IE) = XYZ(1,MLE)
	EGDATA(6,IE) = XYZ(2,MLE)
	EGDATA(7,IE) = XYZ(3,MLE)
	EGDATA(8,IE) = XYZ(4,MLE)
	EGDATA(9,IE) = XYZ(5,MLE)
	EGDATA(10,IE)= XYZ(6,MLE)
	EGDATA(11,IE)= VR(1)
	EGDATA(12,IE)= VR(2)
	EGDATA(13,IE)= VR(3)
	EGDATA(14,IE)= ELN
	EGDATA(15,IE)= RANG
	EGDATA(16,IE)= PLNLG(NLN)	
	EGDATA(17,IE)= PLNLG2
	EGDATA(18,IE)= FAL
	EGDATA(19,IE)= FIM
	EGDATA(20,IE)= FLOAT(NN1)
	EGDATA(21,IE)= FLOAT(NN2)
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
	DO IEF = 1,14
	LEQ = LM(IEF,MLE)
	EGDATA(21+IEF,IE) = FLOAT(LEQ)
	LRP = IPIN(IEF)
	EGDATA(35+IEF,IE) = FLOAT(LRP)
	ENDDO
	
	EGDATA(50,IE)= FLOAT(NLN)
	EGDATA(65,IE)= ECC
C	----------------------
	PLNLG(NLN)  = PLNLG2

100	CONTINUE
	ENDDO

	RETURN

	END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANDATN(ID,IDRCT,XYZ,NSF,NSN,
	1				   NE,LANEL,ECLAN,FILAN,LANNS,
	2				   EGDATA,NEDATA,PLNLG,LANEN,ECCLN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	DIMENSION  XYZ(NSN,1),ID(NSF,NSN),IDRCT(NSF,NSN)
	DIMENSION  VR(3),RG(24)
	DIMENSION  LANEL(1),ECLAN(1),FILAN(1),LANNS(1)
	DIMENSION  EGDATA(NEDATA,1),LANEN(2,1)
	DIMENSION  PLNLG(1),ECCLN(1)


	ITEST = 0 !INITIALIZE LANE DIRECTION FLAG

	DO IE = 1,NE

	NN1 = LANEN(1,IE)
	NN2 = LANEN(2,IE)
	FAL = ECLAN(IE)
	FIM = FILAN(IE)
	NLN = LANNS(IE)
	ECC = ECCLN(IE)

	VR(1) = XYZ(NN2,1)-XYZ(NN1,1)
	VR(2) = XYZ(NN2,2)-XYZ(NN1,2)
	VR(3) = XYZ(NN2,3)-XYZ(NN1,3)

	CALL SCALEN(VR,VR,ELN,3)
	PLNLG2 = PLNLG(NLN) + ELN
C	----------------------
	EGDATA(1,IE) = 0.0D0 !FLOAT(IGM)
	EGDATA(2,IE) = 0.0D0 !FLOAT(KEG)
	EGDATA(3,IE) = 0.0D0 !FLOAT(MLE)
	EGDATA(4,IE) = FLOAT(ITEST)
	EGDATA(5,IE) = XYZ(NN1,1) !XYZ(1,MLE)
	EGDATA(6,IE) = XYZ(NN1,2) !XYZ(2,MLE)
	EGDATA(7,IE) = XYZ(NN1,3) !XYZ(3,MLE)
	EGDATA(8,IE) = XYZ(NN2,1) !XYZ(4,MLE)
	EGDATA(9,IE) = XYZ(NN2,2) !XYZ(5,MLE)
	EGDATA(10,IE)= XYZ(NN2,3) !XYZ(6,MLE)
	EGDATA(11,IE)= VR(1)
	EGDATA(12,IE)= VR(2)
	EGDATA(13,IE)= VR(3)
	EGDATA(14,IE)= ELN
	EGDATA(15,IE)= 0.0D0      !RANG
	EGDATA(16,IE)= PLNLG(NLN)	
	EGDATA(17,IE)= PLNLG2
	EGDATA(18,IE)= FAL
	EGDATA(19,IE)= FIM
	EGDATA(20,IE)= FLOAT(NN1)
	EGDATA(21,IE)= FLOAT(NN2)
	DO IEF = 1,14
	EGDATA(21+IEF,IE) = 0.0D0
	EGDATA(50+IEF,IE) = 0.0D0
	EGDATA(35+IEF,IE) = 0.0D0 !FLOAT(LRP)
	ENDDO
	DO IEF = 1,7
	LEQ = ID(IEF,NN1)
	EGDATA(21+IEF,IE) = FLOAT(LEQ)
	LEQ = ID(IEF,NN2)
	EGDATA(28+IEF,IE) = FLOAT(LEQ)
	ENDDO
	EGDATA(50,IE)= FLOAT(NLN)
	DO IEF = 1,7
	LEQ = IDRCT(IEF,NN1)
	EGDATA(50+IEF,IE) = FLOAT(LEQ)
	LEQ = IDRCT(IEF,NN2)
	EGDATA(57+IEF,IE) = FLOAT(LEQ)
	ENDDO
	EGDATA(65,IE)= ECC
C	----------------------
	PLNLG(NLN)  = PLNLG2

100	CONTINUE
	ENDDO

	RETURN

	END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANLFP(NE,LP,EGDATA,NEDATA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	=============================================================
	DIMENSION  EGDATA(NEDATA,1)

	DO IE = 1,NE
	WRITE(LP) (EGDATA(JJ,IE),JJ=1,NEDATA)
	ENDDO

	RETURN

	END



C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE LANLOADF(NE,IDREC,NPL,IDSET,DIRCOS,NLS,LP,PLNLG,
	1				    EGDATA,NEDATA,IOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	COMMON /MGRAV/ NGRAV
	COMMON /BF_SMOTH/ NP_SMH
C	=============================================================
	DIMENSION  LM(14),LRPIN(14)
	DIMENSION  VR(3),XYZ(6)
	DIMENSION  LANEP(5,1)
	DIMENSION  IDSET(1),DIRCOS(1)
	DIMENSION  EGDATA(NEDATA,1)
	DIMENSION  LED(100),RG(1000),RL(1000),RS(1000)
	DIMENSION  PLNLG(1),LEDR(100),VG(3),VE(3),VME(3)


	VG(1:3)   = 0.0D0
	VG(NGRAV) =-1.0D0

	DO 2000 IPL = 1,NPL

	NEFX = 0

	DO IE = 1,NE


	FLN1  = EGDATA(16,IE)
	FLN2  = EGDATA(17,IE)

	NLN   = INT(EGDATA(50,IE))

	DLEN  = PLNLG(NLN)/FLOAT(NPL-1)
	FLEN  = DLEN*(IPL-1)
	IF(IPL.NE.1.AND.IPL.NE.NPL) FLEN = FLEN-DLEN+0.999*DLEN


	NE1 = IE
	NE2 = IE + 1
	IF(IE.NE.NE) THEN
	NLN1 = INT(EGDATA(50,IE+1))
	IF(NLN.NE.NLN1) NE1 = IE - 1
	IF(NLN.NE.NLN1) NE2 = IE
	ENDIF
	IF(IE.EQ.NE) NE1 = IE - 1
	IF(IE.EQ.NE) NE2 = IE


	IF(IPL.NE.NPL) THEN
	IF(FLEN.LT.FLN1.OR.FLEN.GT.FLN2) GOTO 500
	ELSE
	IF(FLEN.LT.FLN1.OR.FLEN.GT.FLN2*1.001) GOTO 500
	ENDIF

	NEFX = NEFX + 1
	NUMA = 1 + 14*(NEFX-1)
	NUMI = 1 + 17*(NEFX-1)

	IGM   = INT(EGDATA(1,IE))
	KEG   = INT(EGDATA(2,IE))
	MLE   = INT(EGDATA(3,IE))
	VR(1) = EGDATA(11,IE)
	VR(2) = EGDATA(12,IE)
	VR(3) = EGDATA(13,IE)
	ELN   = EGDATA(14,IE)
	RANG  = EGDATA(15,IE)
	FAL   = EGDATA(18,IE)
	FIM   = 1.0D0 + EGDATA(19,IE)
	ECT   = EGDATA(65,IE)

	LED(NUMI+0) = IGM
	LED(NUMI+1) = KEG
	LED(NUMI+2) = MLE
	DO IEF = 1,14
	 LED(IEF+NUMI+3) = INT(EGDATA(21+IEF,IE))   !LM
	LEDR(IEF+NUMI+3) = INT(EGDATA(50+IEF,IE))   !LMRCT
	LRPIN(IEF)       = INT(EGDATA(35+IEF,IE))
	ENDDO


	PM   = -1.0D0*FAL*FIM   ! MINUS FOR NEGATIVE GRAVITY
	CALL LANDIRC(NE1,NE2,VR,ITEST,EGDATA,NEDATA)
	EGDATA(4,IE) = FLOAT(ITEST)                     !STORE ITEST

	IF(ITEST.EQ.0) CALL VECPRD(VG,VR,VE)
	IF(ITEST.EQ.1) CALL VECPRD(VR,VG,VE)
	VE = VE*ECT

	IF(ITEST.EQ.0) AL = FLEN-FLN1	
	IF(ITEST.EQ.1) AL = ELN - (FLEN-FLN1)
	BL = ELN - AL
	IDR  = NGRAV
 
	SELECTCASE(IOPT)
	CASE(0) !LANE ON FRAME ELEM
	CALL LANFFIX(PM,AL,ECT,IDR,VR,ELN,RG(NUMA),RANG,LRPIN,VE)
	CALL FRMMOV (RL(NUMA),RG(NUMA),VR,RANG)
	
	W2 = 0.0D0
	BL = AL
	NP = NP_SMH
	NPH = 1 + 6*NP_SMH*(NEFX-1)
	CALL LNUNIF(IDR,PM,W2,AL,BL,ECT,VR,ELN,RANG,NP,RS(NPH),VE) 


	CASE(1) !LANE ON NODE
	RG(NUMA:NUMA+13) = 0.0D0
	RG(NUMA+IDR-1+0) = PM*BL/ELN
	RG(NUMA+IDR-1+7) = PM*AL/ELN
	CALL VECPRD(VE,RG(NUMA+1-1+0),VME)
	DO IJ = 1,3
	RG(NUMA+IJ-1+3)  = RG(NUMA+IJ-1+3)  + VME(IJ) 
	ENDDO
	CALL VECPRD(VE,RG(NUMA+1-1+7),VME)
	DO IJ = 1,3
	RG(NUMA+IJ-1+10) = RG(NUMA+IJ-1+10) + VME(IJ) 
	ENDDO

	ENDSELECT


	NN1 = INT(EGDATA(20,IE)) 
	NN2 = INT(EGDATA(21,IE)) 
      IF (NLS.NE.0) CALL LOPREB (IDSET,DIRCOS,RG(NUMA),NN1,NN2)


500	CONTINUE
 
	ENDDO



	WRITE(LP) IPL,NEFX,IOPT
	DO IEFX = 1,NEFX

	NUMA  = 1 + 14*(IEFX-1)
	NUMI  = 1 + 17*(IEFX-1)
	IGM = LED(NUMI+0)
	KEG = LED(NUMI+1)
	MLE = LED(NUMI+2)
	WRITE(LP) IGM,KEG,MLE
	DO IEF = 1,14 !FORCE IN LOCAL SUPPORT FOR RHS
	IEQ = LED(IEF+NUMI+3)
	WRITE(LP) IEQ,RG(IEF+NUMA-1)
	ENDDO
	DO IEF = 1,14 !FORCE IN ELEMENT LOCAL SYSTEM FOR FIXEND
	IF(IOPT.EQ.0) IEQ = LED(IEF+NUMI+3)
	IF(IOPT.EQ.0) WRITE(LP) IEQ,RL(IEF+NUMA-1) !ON FRAME -- PUT TO ELEMEENT FORCE -- IT WILL TRANSFER TO REACTION LATER IN ELLOOP
	IF(IOPT.EQ.1) IEQ = LEDR(IEF+NUMI+3)
	IF(IOPT.EQ.1) WRITE(LP) IEQ,RG(IEF+NUMA-1) !ON NODE -- PUT DIRECTLY TO GLOBAL REACTION -- IN ELLOOP
	ENDDO

	IF(IOPT.EQ.0) THEN
	DO IEF = 1,NP_SMH !ELEMENT FORCE FOR SMOOTH LINE DIAGRAM
	II = 6*(IEF-1) + 6*NP_SMH*(IEFX-1)
	WRITE(LP) (RS(II+J),J=1,6)
	ENDDO
	ENDIF

	ENDDO

2000	CONTINUE



	
	RETURN

	END



C	=============================================================
C	=============================================================
C	=============================================================

	SUBROUTINE LANFIXF1(LP,NPL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	
	COMMON /LANEFIX/ ALF1(50000),ALF2(50000),ALF3(50000),ILF(10000)
	COMMON /BF_SMOTH/ NP_SMH

	DIMENSION VV(6)

	REWIND(LP)

	NMI = 0
	NMA = 0
	NMS = 0

	DO 2000 ILPL = 1,NPL

	READ(LP) IPL,NEFX,IOPT

	DO 999 IEFX = 1,NEFX
	NUMI = 1 +  4*(IEFX-1) + NMI
	READ(LP) IGM,KEG,MLE
	ILF(NUMI+0) = IGM
	ILF(NUMI+1) = KEG
	ILF(NUMI+2) = MLE
	ILF(NUMI+3) = NEFX
	DO IEF = 1,14 !FORCE IN LOCAL SUPPORT FOR RHS
	N1 = 2*IEF-1 + 14*(IEFX-1)*2 + NMA
	N2 = 2*IEF-0 + 14*(IEFX-1)*2 + NMA
	READ(LP) IEQ,VALV
	ALF1(N1) = FLOAT(IEQ)
	ALF1(N2) = VALV
	ENDDO
	DO IEF = 1,14 !FORCE IN ELEMENT LOCAL SYSTEM FOR FIXEND
	N1 = 2*IEF-1 + 14*(IEFX-1)*2 + NMA
	N2 = 2*IEF-0 + 14*(IEFX-1)*2 + NMA
	READ(LP) IEQ,VALV
	ALF2(N1) = FLOAT(IEQ)
	ALF2(N2) = VALV
	ENDDO

	IF(IOPT.EQ.0) THEN		!FOR LANE ON FRAME ELEMENT IOPT = LAN_OPT = 0
	DO IEF = 1,NP_SMH		!SMOOTH LINE DIAGRAHM
	N1 = 7*(IEF-1) + 7*NP_SMH*(IEFX-1) + NMS
	READ(LP) VV(1:6)
	DO J = 1,6
	ALF3(N1+J) = VV(J)
	ENDDO
	ENDDO
	ENDIF

999	CONTINUE

	NMI = NMI +  4*NEFX
	NMA = NMA + 14*NEFX*2
	NMS = NMS + 7*NP_SMH*NEFX
	 
2000	CONTINUE


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANFIXFV(NEFX,RG,RL,LED,LER,IOPT,RS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	
	COMMON /LANEFIX/ ALF1(50000),ALF2(50000),ALF3(50000),ILF(10000)
	COMMON /BF_SMOTH/ NP_SMH

	DIMENSION LED(1),LER(1),RG(1),RL(1),RS(1)


	DO IEFX = 1,NEFX
	NUMI = 1 +  4*(IEFX-1) 
	NA   = 1 + 14*(IEFX-1)
	NI   = 1 + 17*(IEFX-1)
	IGM = LED(NI+0)
	KEG = LED(NI+1)
	MLE = LED(NI+2)

	ILF(NUMI+0) = IGM
	ILF(NUMI+1) = KEG
	ILF(NUMI+2) = MLE
	ILF(NUMI+3) = NEFX

	DO IEF = 1,14 !FORCE IN LOCAL SUPPORT FOR RHS
	N1 = 2*IEF-1 + 14*(IEFX-1)*2 
	N2 = 2*IEF-0 + 14*(IEFX-1)*2 

	IEQ  = LED(IEF+NI+3)
	VALV =  RG(IEF+NA-1)

	ALF1(N1) = FLOAT(IEQ)
	ALF1(N2) = VALV
	ENDDO
	DO IEF = 1,14 !FORCE IN ELEMENT LOCAL SYSTEM FOR FIXEND
	N1 = 2*IEF-1 + 14*(IEFX-1)*2 
	N2 = 2*IEF-0 + 14*(IEFX-1)*2 

	IF(IOPT.EQ.0) IEQ  = LED(IEF+NI+3)
	IF(IOPT.EQ.0) VALV =  RL(IEF+NA-1)

	IF(IOPT.EQ.1) IEQ  = LER(IEF+NI+3)
	IF(IOPT.EQ.1) VALV =  RG(IEF+NA-1)

	ALF2(N1) = FLOAT(IEQ)
	ALF2(N2) = VALV
	ENDDO

	IF(IOPT.EQ.0) THEN		!FOR LANE ON FRAME ELEMENT IOPT = LAN_OPT = 0
	DO IEF = 1,NP_SMH		!SMOOTH LINE DIAGRAHM
	N1 = 7*(IEF-1) + 7*NP_SMH*(IEFX-1)
	N2 = 6*(IEF-1) + 6*NP_SMH*(IEFX-1)
	DO J = 1,6
	ALF3(N1+J) = RS(N2+J)
	ENDDO
	ENDDO
	ENDIF

	ENDDO


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANFIXF2(RHS,ILPL,NPL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	

	COMMON /LANEFIX/ ALF1(50000),ALF2(50000),ALF3(50000),ILF(10000)
	DIMENSION RHS(1)

	NMI = 0
	NMA = 0
	DO 2000 IPL = 1,NPL

	NEFX = ILF(NMI+4)

	IF(IPL.EQ.ILPL) THEN
	DO IEFX = 1,NEFX
	DO IEF = 1,14
	N1 = 2*IEF-1 + 14*(IEFX-1)*2 + NMA
	N2 = 2*IEF-0 + 14*(IEFX-1)*2 + NMA
	IEQ  = INT(ALF1(N1))
	VALV = ALF1(N2)
	IF(IEQ.NE.0) RHS(IEQ) = RHS(IEQ) + VALV
	ENDDO
	ENDDO
	RETURN
	ENDIF


	NMI = NMI +  4*NEFX
	NMA = NMA + 14*NEFX*2
	
2000	CONTINUE


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANFIXF3(FIXEN,IEL,IEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	
	COMMON /LANPRIN/  LAN_PRN,ILPL,MLANE,NPL,LAN_OPT
	COMMON /LANEFIX/ ALF1(50000),ALF2(50000),ALF3(50000),ILF(10000)
	DIMENSION FIXEN(1)


	IF(LAN_OPT.EQ.1) RETURN

	NMI = 0
	NMA = 0
	DO 2000 IPL = 1,NPL

	NEFX = ILF(NMI+4)

	IF(IPL.EQ.ILPL) THEN
	DO IEFX = 1,NEFX
	NUMI = 1 +  4*(IEFX-1) + NMI
	IGM = ILF(NUMI+0)
	KEG = ILF(NUMI+1)
	MLE = ILF(NUMI+2)

	IF(IEL.EQ.MLE.AND.IEG.EQ.KEG) THEN
	DO IEF = 1,14
	N1 = 2*IEF-1 + 14*(IEFX-1)*2 + NMA
	N2 = 2*IEF-0 + 14*(IEFX-1)*2 + NMA
	VALV = ALF2(N2)
	FIXEN(IEF) = FIXEN(IEF) + VALV
	ENDDO
	ENDIF
	ENDDO
	RETURN
	ENDIF


	NMI = NMI +  4*NEFX
	NMA = NMA + 14*NEFX*2
	
2000	CONTINUE


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANFIXF4(RHS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	

	COMMON /LANPRIN/  LAN_PRN,ILPL,MLANE,NPL,LAN_OPT
	COMMON /LANEFIX/ ALF1(50000),ALF2(50000),ALF3(50000),ILF(10000)
	DIMENSION RHS(1)

	IF(LAN_OPT.EQ.0) RETURN

	NMI = 0
	NMA = 0
	DO 2000 IPL = 1,NPL

	NEFX = ILF(NMI+4)

	IF(IPL.EQ.ILPL) THEN
	DO IEFX = 1,NEFX
	DO IEF = 1,14
	N1 = 2*IEF-1 + 14*(IEFX-1)*2 + NMA
	N2 = 2*IEF-0 + 14*(IEFX-1)*2 + NMA
	IEQ  = INT(ALF2(N1))
	VALV = ALF2(N2)
	IF(IEQ.NE.0) RHS(IEQ) = RHS(IEQ) - VALV
	ENDDO
	ENDDO
	RETURN
	ENDIF


	NMI = NMI +  4*NEFX
	NMA = NMA + 14*NEFX*2
	
2000	CONTINUE


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANDIRC(IE1,IE2,VR,ITEST,EGDATA,NEDATA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	=============================================================
	DIMENSION  LM(14),LRPIN(14),VK(3)
	DIMENSION  VR(3),RG(24),XYZ(6)
	DIMENSION  VTS(3)
	DIMENSION  EGDATA(NEDATA,1)


	XYZ(1) = EGDATA(5 ,IE1) 
	XYZ(2) = EGDATA(6 ,IE1) 
	XYZ(3) = EGDATA(7 ,IE1) 
	XYZ(4) = EGDATA(8 ,IE1) 
	XYZ(5) = EGDATA(9 ,IE1) 
	XYZ(6) = EGDATA(10,IE1)
		
	X1 = 0.5*(XYZ(1) + XYZ(4))
	Y1 = 0.5*(XYZ(2) + XYZ(5))
	Z1 = 0.5*(XYZ(3) + XYZ(6))

	XYZ(1) = EGDATA(5 ,IE2) 
	XYZ(2) = EGDATA(6 ,IE2) 
	XYZ(3) = EGDATA(7 ,IE2) 
	XYZ(4) = EGDATA(8 ,IE2) 
	XYZ(5) = EGDATA(9 ,IE2) 
	XYZ(6) = EGDATA(10,IE2)

	X2 = 0.5*(XYZ(1) + XYZ(4))
	Y2 = 0.5*(XYZ(2) + XYZ(5))
	Z2 = 0.5*(XYZ(3) + XYZ(6))

	VTS(1) = X2 - X1
	VTS(2) = Y2 - Y1
	VTS(3) = Z2 - Z1

	CALL SCALEN(VTS,VTS,DUM,3)

	TEST = VTS(1)*VR(1) + VTS(2)*VR(2) + VTS(3)*VR(3)


	ITEST = 1
	IF(TEST.GT.0.0D0) ITEST = 0 

C	IF(ILANE.EQ.2) THEN
C	WRITE(*,*) VR,VTS,ITEST
C	PAUSE
C	ENDIF
	
	RETURN

	END



C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE LANFFIX(PM,AL,ECT,IDR,VR,ELN,RG,ANG,LREAS,VE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	==============================================================
	DIMENSION VR(3),VS(3),VT(3),RG(7,2)
	DIMENSION COEF(6),TRANS(14,14)
	DIMENSION FIXG(14),FIXD(14)
	DIMENSION TRANH(14,14),LREAS(14)
	DIMENSION VE(3),VME(3)


	RANG = ANG

	ECT  = 0.0D0
	ECS  = 0.0D0

	CALL FMVEVR(VR,VS,VT)
	CALL ROMBAC(VR,VS,VT,RANG)
	CALL TRANLG(VR,VS,VT,TRANS)


	PR = VR(IDR)*PM
	PS = VS(IDR)*PM
	PT = VT(IDR)*PM

	PMR = -ECT*PS
	PMS =  ECT*PR
	PMT =  0.0D0

	FIXD(1:14) = 0.0D0

C	FROM CONCENTRIC LOAD

C	LOCAL AXIAL FORCE
	CALL FXLANE(PR,AL,ELN,COEF,0)
	FIXD(1) = COEF(5)
	FIXD(8) = COEF(6)

C	SHEAR IN S-AXIS AND MOMENT IN T-AXIS
	CALL FXLANE(PS,AL,ELN,COEF,0)
	FIXD(2)  = COEF(1)
	FIXD(6)  = COEF(2)
	FIXD(9)  = COEF(3)
	FIXD(13) = COEF(4)	

C	SHEAR IN T-AXIS AND MOMENT IN S-AXIS
	CALL FXLANE(PT,AL,ELN,COEF,0)
	FIXD(3)  = COEF(1)
	FIXD(5)  =-COEF(2)
	FIXD(10) = COEF(3)
	FIXD(12) =-COEF(4)	

C	FROM ECCENTRICITY MOMENT

C	LOCAL AXIAL MOMENT
	CALL FXLANE(PMR,AL,ELN,COEF,1)
	FIXD(4)  = FIXD(4)  + COEF(5)
	FIXD(11) = FIXD(11) + COEF(6)


C	SHEAR IN T-AXIS AND MOMENT IN S-AXIS
	CALL FXLANE(PMS,AL,ELN,COEF,1)
	FIXD(3)  = FIXD(3)  - COEF(1)
	FIXD(5)  = FIXD(5)  + COEF(2)
	FIXD(10) = FIXD(10) - COEF(3)
	FIXD(12) = FIXD(12) + COEF(4)


C	SHEAR IN S-AXIS AND MOMENT IN T-AXIS
	CALL FXLANE(PMT,AL,ELN,COEF,1)
	FIXD(2)  = FIXD(2)  + COEF(1)
	FIXD(6)  = FIXD(6)  + COEF(2)
	FIXD(9)  = FIXD(9)  + COEF(3)
	FIXD(13) = FIXD(13) + COEF(4)

C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
	CALL TRNHIG(TRANH,ELN,LREAS)
	CALL TRNMUL(TRANH,FIXD,2)

	FIXG = MATMUL(TRANS,FIXD)

	DO I = 1,6
	RG(I,1) = 0.0D0
	RG(I,2) = 0.0D0
	RG(I,1) = FIXG(I)
	RG(I,2) = FIXG(I+7)
	ENDDO

	CALL VECPRD(VE,RG(1,1),VME)
	DO I = 1,3
	RG(3+I,1) = RG(3+I,1) + VME(I) 
	ENDDO

	CALL VECPRD(VE,RG(1,2),VME)
	DO I = 1,3
	RG(3+I,2) = RG(3+I,2) + VME(I) 
	ENDDO



      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================

	SUBROUTINE FXLANE(W,A,ELN,COEF,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION C(6),COEF(6)
C	=====================================================================
	
	ELN2 = ELN*ELN
	ELN3 = ELN*ELN2

C	---------------------------------------------------------
	DA = A
	DA2 = DA*DA
	DA3 = DA*DA2
	DA4 = DA*DA3

	DB = ELN - A
	DB2 = DB*DB
	DB3 = DB*DB2
	DB4 = DB*DB3

	SELECTCASE(IND)

	CASE(0)
	C(1) =  DB3/ELN3 + 3.0*DB2*DA/ELN3
	C(2) =  DB2*DA/ELN2
	C(3) =  DA3/ELN3 + 3.0*DA2*DB/ELN3
	C(4) = -DA2*DB/ELN2 
	C(5) =  DB/ELN
	C(6) =  DA/ELN
	CASE(1)
	C(1) =  -6.0*DA*DB/ELN3
	C(2) =  -2.0*DA*DB/ELN2 + DB2/ELN2
	C(3) =   6.0*DA*DB/ELN3
	C(4) =   2.0*DA*DB/ELN2 - DA2/ELN2
	C(5) =   DB/ELN
	C(6) =   DA/ELN

	ENDSELECT
C	---------------------------------------------------------

	DO I = 1,6
	COEF(I) = C(I)*W
	ENDDO


	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FRMMOV (RL,RG,VR,ANG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION VR(3),VS(3),VT(3),TRANS(7,7)
	DIMENSION RL(7,2),RG(7,2)

	
	RANG = ANG

	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

	TRANS = 0.0D0
	DO I = 1,3
	TRANS(I,1) = VR(I)
	TRANS(I,2) = VS(I)
	TRANS(I,3) = VT(I)
	TRANS(I+3,4) = VR(I)
	TRANS(I+3,5) = VS(I)
	TRANS(I+3,6) = VT(I)
	ENDDO

	RL = 0.0D0
	DO INM = 1,2

	DO IFF = 1,6
	DO JFF = 1,6
	RL(IFF,INM) = RL(IFF,INM) + TRANS(JFF,IFF)*RG(JFF,INM)
	ENDDO
	ENDDO

	ENDDO


	RETURN

	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE INFLCAL(AA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	

      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
	COMMON A(9000000),IA(9000000)
	DIMENSION AA(1)

	CALL INFLCAL1(IA(LID),IA(LMA),A(LDK),NEQ,AA)


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE INFLCAL1(ID,MAXA,D,NEQ,AA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

	COMMON /LANPRIN/  LAN_PRN,ILPL,MLANE,NPL,LAN_OPT


	COMMON A(9000000),IA(9000000)

	DIMENSION RHS(NEQ)
	DIMENSION D(NEQ),MAXA(NEQ+1),ID(NSF,NSN)
	DIMENSION DD(6)
	DIMENSION PLNLG(10)
	DIMENSION SPANF(101)
	DIMENSION AA(1)
	ALLOCATABLE EDATA(:)


	INDPD = 0 !POSITIVE DEFINITE STIFF

	LAN_PRN = 1

	CALL INTFILL('LANE',NLANE ,1 ,1,0)  !NUMBER OF LANE
	CALL INTFILL('LANE',NEDATA,1 ,2,0)  !NUMBER OF DATA FOR ELEMENT
	CALL INTFILL('LANE',NSPANX,1 ,3,0)  !MAX NUMBER OF SPAN

	ALLOCATE(EDATA(NEDATA))

	CALL INTFILL('LNMX',NFIL1 ,1 ,1 ,0)  !LANE DATA	
	CALL INTFILL('LNMX',NFIL2 ,1 ,2 ,0)  !LANE DATA	
	CALL INTFILL('LNMX',NFIL3 ,1 ,3 ,0)  !LANE DATA	
	CALL INTFILL('LNMX',NFIL5 ,1 ,5 ,0)  !LANE DATA	
	CALL INTFILL('LNMX',NFIL6 ,1 ,6 ,0)  !LANE DATA	

      
	DO 1000 ILANE = 1,NLANE

      IPROG = 0
      
	WRITE(*,*) 'PERFORM INFLUENCE LINE CALCULATION OF LANE',ILANE

	MLANE = ILANE
	LP1 = NFIL1 + ILANE
	REWIND(LP1)
	READ(LP1) IL,NLNS,NE,IDREC,NPL,IOPT,NSPAN,WTHLN		!IDREC  1= NORMAL  2= REVERSE  3= BOTH---NPL NUMBER OF INFLUENCE POINT ON THE LANE
	LAN_OPT = IOPT
	READ(LP1) (SPANF(J),J=1,NSPAN+1)
	READ(LP1) (PLNLG(J),J=1,NLNS)
	PLNAV = 0.0D0
	DO J = 1,NLNS
	PLNAV = PLNAV + PLNLG(J)
	ENDDO
	PLNAV = PLNAV/FLOAT(NLNS)

	LP5 = ILANE
	FF = 0.0D0
	CALL INIOPER(NFIL5,LP5,FF,'CLEF','OLDF')  !CLEAR FILE FOR INFLUENCE LINE INTEGRATION (MAX RESPONSE FROM UNIFORM LAOD)

	DO IE = 1,NE
	READ(LP1) (EDATA(JJ),JJ=1,NEDATA)
	ENDDO

	LP2 = NFIL2 + ILANE

	CALL LANFIXF1(LP2,NPL)

	DO 600 ISPAN = 1,NSPAN
	SPF1 = SPANF(ISPAN+0)*PLNAV 
	SPF2 = SPANF(ISPAN+1)*PLNAV*0.999
	IF(ISPAN.EQ.NSPAN) SPF2 = SPANF(ISPAN+1)*PLNAV*1.001

	LP6 = ISPAN + NSPANX*ILANE
	FF = 0.0D0
	CALL INIOPER(NFIL6,LP6,FF,'CLEF','OLDF')  !CLEAR FILE FOR MAX RESPONSE FOR EACH SPAN DUE TO UNIT LOAD 

	DLEN = PLNAV/FLOAT(NPL-1)
	DO 500 ILPL = 1,NPL	
	FLEN = DLEN*(ILPL-1)

	IF(FLEN.LT.SPF1.OR.FLEN.GT.SPF2) GOTO 500 !JUMP HERE IF THE POINT NOT FALLING IN CURRENT SPAN 

C	WRITE(*,*) 'INFLUENCE POINT NO.',ILPL

	CALL CLEARA(RHS,NEQ)
	CALL LANFIXF2(RHS,ILPL,NPL)
	
	CALL COLSOL (MAXA,AA,D,RHS,2,INDPD,'TEMP','TEMP')        !GETTING DISPLACEMENT


	CALL MOVE (RHS,A(LDL),NEQ)
	CALL MOVE (A(LDL),A(LDT),NEQ)


	CALL CLEROUT

C	NEW OUTPUT SONGSAK JUL2007
	CALL DISOUT(ID,RHS)

      ITASK = 3
      IFREF = 1
      ISPRI = 0
      CALL GRLOOP (IA(LEL),KSC)

	CALL INFINTG(NFIL5,LP5,PLNAV,ILPL,NPL)        !INFLUENCE LINE INTEGRATION

	FF = 1.0D0
	CALL INIOPER(NFIL6,LP6,FF,'SORT','OLDF')      !SORT TO FILE (MAX RESPONSE FOR EACH SPAN DUE TO UNIT LOAD)

	FF = 1.0D0
	LP3 = ILPL + NPL*(ILANE-1)
	CALL INIOPER(NFIL3,LP3,FF,'REPC','OLDF')      !RECORD RESPONSE DUE TO UNIT LOAD (AT EACH STATION) TO FILE    

      MPROG = NPL
      IPROG = IPROG + 1
      CALL PROGBAR(IPROG,MPROG)

C	WRITE(101,8000) MLANE,ILPL
C
C	DO I = 1,NSN
C	DD(1:6) = 0.0D0
C	DO J = 1,6
C	IEQ = ID(J,I)
C	IF(IEQ.NE.0) DD(J) = DD(J) + RHS(IEQ)
C	ENDDO
C	WRITE(101,8100) I,DD(1:6)
C	ENDDO
C	WRITE(101,8200)

500	CONTINUE


600	CONTINUE

1000	CONTINUE


 8000 FORMAT ('Result "LaneDisplacement"',2x,'"MovingLoad Lane ',I3,'"',
     1         2X,I5,5X,'Matrix',2X,'OnNodes'/
     2        'ComponentNames  "X-Displacement", "Y-Displacement"',
     3        ', "Z-Displacement", "X-Rotation", "Y-Rotation", 
	1"Z-Rotation"'/'Values') 

 8100	FORMAT (I5,3X,7E12.4)

 8200 FORMAT ('End Values'/)


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE INTERPOL(XX,VV,DAT,NUMT,IT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	PRODUCED BY SONGSAK JUN2007
C	LINEAR INTERPOLATION OF DATA
C	=============================================================
	DIMENSION DAT(2,1)

	NP = NUMT/4       !NUM 1ST PICK
	IF(NP.LE.1) NP = 1

	IF(NP.GT.NUMT) NP = NUMT
	I1 = 1
	I2 = NP

	DO I = 1,NUMT
	
	IF(I2.GT.NUMT) THEN
	I2 = NUMT
	I1 = I2 - NP + 1
	GOTO 5
	ENDIF

	TT1 = DAT(1,I1)
	TT2 = DAT(1,I2)

	IF(XX.GE.TT1.AND.XX.LE.TT2) GOTO 5

	I1 = I2
	I2 = I2 + NP

	ENDDO

	I1 = 1
	I2 = NUMT

5	CONTINUE


	DO I = I1,I2
      IF(I.EQ.I2.AND.I2.EQ.NUMT) EXIT
	TT1 = DAT(1,I  )
	TT2 = DAT(1,I+1)

	IF(XX.GE.TT1.AND.XX.LE.1.0000001*TT2) THEN
	IT = I
	GOTO 10
	ENDIF 
	ENDDO

	VV = 0.0D0
	GOTO 20

10	CONTINUE

	DT = XX -TT1
	TT = TT2-TT1
	H1 = 1.0 - DT/TT
	H2 = DT/TT

	VV = DAT(2,IT)*H1 + DAT(2,IT+1)*H2

20	CONTINUE


	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE BETWEEN(XX,XL,NUMT,IT1,IT2,AL,EL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	PRODUCED BY SONGSAK JUN2007
C	SEARCH THE POINT THAT DATA FALLING BETWEEN
C	=============================================================
	DIMENSION XL(1)

	NP = NUMT/4       !NUM 1ST PICK
	IF(NP.LE.1) NP = 1

	IF(NP.GT.NUMT) NP = NUMT
	I1 = 1
	I2 = NP

	DO I = 1,NUMT
	
	    IF(I2.GT.NUMT) THEN
	        I2 = NUMT
	        I1 = I2 - NP + 1
	        GOTO 5
	    ENDIF

	    TT1 = XL(I1)
	    TT2 = XL(I2)

	    IF(XX.GE.TT1.AND.XX.LE.TT2) GOTO 5

	    I1 = I2
	    I2 = I2 + NP

	ENDDO

	I1 = 1
	I2 = NUMT

5	CONTINUE


	IT1 = 0
	IT2 = 0
	AL = 0.0D0
	EL = 0.0D0

	DO I = I1,I2
	    IF(I.EQ.I2.AND.I2.EQ.NUMT) EXIT
	    TT1 = XL(I  )
	    TT2 = XL(I+1)

	    IF(XX.GE.TT1.AND.XX.LE.1.0000001*TT2) THEN
	          IT1 = I	
	          IT2 = I	+ 1
                AL = XX -TT1
                EL = TT2-TT1
	        RETURN
	    ENDIF 
	ENDDO



	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE INTGFACT(XP,NP,TL,NL,NPF,IPF,FAC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	PRODUCED BY SONGSAK JUN2007
C	SEARCH THE POINT THAT DATA FALLING BETWEEN
C	=============================================================
	DIMENSION XP(1),XL(NL),IPF(1),FAC(1)
	DIMENSION IFD(2,NP),DATL(2,NP)


      XL(1)= 0.0D0
      DO IL = 2,NL
      IF(IL.EQ.2) DL = TL/(NL-1) 
      XL(IL) = DL*(IL-1)
      ENDDO
      
      
      IFD  = 0
      DATL = 0.0D0
      
      IFOUND = 0
      DO 1000 IP = 1,NP
      
      XX = XP(IP)
      CALL BETWEEN(XX,XL,NL,IT1,IT2,AL,EL)
      IF(IT1.NE.0) THEN
        IFOUND = IFOUND + 1
        IFD(1,IP) = IT1
        IFD(2,IP) = IT2
        DATL(1,IP) = AL
        DATL(2,IP) = EL
      ENDIF
      
1000  CONTINUE
      
      
      NPF = 0
      IPF(1:NPL) = 0
      FAC(1:NPL) = 0.0D0
      
      IF(IFOUND.EQ.0) RETURN
      
      
      IF(NP.EQ.1) THEN
          AL = DATL(1,1)
          EL = DATL(2,1)
          IF(EL.NE.0.0D0) THEN 
              NPF = 2
              IPF(1) = IFD(1,1)
              IPF(2) = IFD(2,1)
              FAC(1) = (EL-AL)/EL
              FAC(2) = AL/EL
          ENDIF
          RETURN
      ENDIF
      
      
      IF(NP.EQ.2) THEN
      
        IF(IFD(1,1).EQ.0) THEN
            IFD(1,1) = NL-1
            IFD(2,1) = NL
            EL = XL(NL) - XL(NL-1)
            AL = EL
            DATL(1,1) = AL
            DATL(2,1) = EL     
        ENDIF
      
        IF(IFD(1,2).EQ.0) THEN
            IFD(1,2) = 1
            IFD(2,2) = 2
            EL = XL(2) - XL(1)
            AL = 0.0D0
            DATL(1,2) = AL
            DATL(2,2) = EL     
        ENDIF
        
        IP1I = IFD(1,1)
        IP1J = IFD(2,1)
        IP2I = IFD(1,2)
        IP2J = IFD(2,2)
        NER = IP1I - IP2J
        NEA = 2 + NER
        
        IF(NEA.EQ.1) THEN
            NPF = 0
            
            AL = DATL(1,2)
            BL = DATL(1,1)
            EL = DATL(2,1)
            IF(EL.NE.0.0D0) THEN
                FC = 0.5*(BL-AL)/EL*(2.0*EL-AL-BL)
                IF(FC.NE.0.0D0) THEN
                    NPF = NPF + 1
                    IPF(NPF) = IFD(1,1)
                    FAC(NPF) = FC
                ENDIF
            ENDIF
            
            AL = DATL(1,2)
            BL = DATL(1,1)
            EL = DATL(2,1)
            IF(EL.NE.0.0D0) THEN
                FC = 0.5*(BL-AL)/EL*(AL+BL)
                IF(FC.NE.0.0D0) THEN
                    NPF = NPF + 1
                    IPF(NPF) = IFD(2,1)
                    FAC(NPF) = FC
                ENDIF
            ENDIF
            
            RETURN
        ENDIF
        
        IPI = IP2J
        IPJ = IP2J + 1
        
        NPF = 0
        DO 90 IEA = 1,NEA
            IF(IEA.EQ.1) THEN
                
                AL = DATL(1,2)
                EL = DATL(2,2)
                BL = EL
                IF(EL.NE.0.0D0) THEN
                    FC = 0.5*(BL-AL)/EL*(2.0*EL-AL-BL)
                    IF(FC.NE.0.0D0) THEN
                        NPF = NPF + 1
                        IPF(NPF) = IFD(1,2)
                        FAC(NPF) = FC
                    ENDIF
                ENDIF
                
                AL = DATL(1,2)
                EL = DATL(2,2)
                BL = EL
                IF(EL.NE.0.0D0) THEN
                    FC = 0.5*(BL-AL)/EL*(AL+BL)
                    IF(FC.NE.0.0D0) THEN
                        NPF = NPF + 1
                        IPF(NPF) = IFD(2,2)
                        FAC(NPF) = FC
                    ENDIF
                ENDIF
                
                GOTO 90
            ENDIF
            IF(IEA.EQ.NEA) THEN
                
                AL = 0.0
                BL = DATL(1,1)
                EL = DATL(2,1)
                IF(EL.NE.0.0D0) THEN
                    FC = 0.5*(BL-AL)/EL*(2.0*EL-AL-BL)
                    IF(FC.NE.0.0D0) THEN
                        NPF = NPF + 1
                        IPF(NPF) = IFD(1,1)
                        FAC(NPF) = FC
                    ENDIF
                ENDIF
                
                AL = 0.0
                BL = DATL(1,1)
                EL = DATL(2,1)
                IF(EL.NE.0.0D0) THEN
                    FC = 0.5*(BL-AL)/EL*(AL+BL)
                    IF(FC.NE.0.0D0) THEN
                        NPF = NPF + 1
                        IPF(NPF) = IFD(2,1)
                        FAC(NPF) = FC
                    ENDIF
                ENDIF
                
                GOTO 90
            ENDIF
            
        
            AL = 0.0
            EL = XL(IPJ) - XL(IPI)
            BL = EL
            IF(EL.NE.0.0D0) THEN
                FC = 0.5*(BL-AL)/EL*(2.0*EL-AL-BL)
                IF(FC.NE.0.0D0) THEN
                    NPF = NPF + 1
                    IPF(NPF) = IPI
                    FAC(NPF) = FC
                ENDIF
            ENDIF
            
            AL = 0.0
            EL = XL(IPJ) - XL(IPI)
            BL = EL
            IF(EL.NE.0.0D0) THEN
                FC = 0.5*(BL-AL)/EL*(AL+BL)
                IF(FC.NE.0.0D0) THEN
                    NPF = NPF + 1
                    IPF(NPF) = IPJ
                    FAC(NPF) = FC
                ENDIF
            ENDIF
            
            IPI = IPI + 1
            IPJ = IPJ + 1
            
            
90      CONTINUE
        
        DO IP = 1,NPF
            IF(IPF(IP).NE.0) THEN
                DO JP = IP+1,NPF
                    IF(IPF(IP).EQ.IPF(JP))  THEN
                        IPF(JP) = 0
                        FAC(IP) = FAC(IP) + FAC(JP)
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        
        MPF = 0
        DO IP = 1,NPF
            IF(IPF(IP).NE.0) THEN
                MPF = MPF + 1
                IPF(MPF) = IPF(IP)
                FAC(MPF) = FAC(IP)
            ENDIF
        ENDDO
        NPF = MPF
        
      ENDIF 
      
      
      
      
      
      
	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE INFINTG(LP,ISET,PLNLG,IPL,NPL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
	DIMENSION FILAR(1)


	DLEN = PLNLG/FLOAT(NPL-1)

	CON1 = 0.5*DLEN
	CON2 = 1.0D0
	IF(IPL.NE.1.AND.IPL.NE.NPL) CON2 = 2.0D0  !TRAPIZOIDAL
	CON = CON1*CON2
	CALL INIOPER(LP,ISET,CON,'COML','OLDF')   !ACCUMULATE TO FILE FOR ISET 


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LANELCAL(NFIL4,LP4,NFIL5,LP5,FAC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------

	CALL CLEROUT
	FACTOR = 1.0D0
	CALL INIOPER(NFIL5,LP5,FACTOR,'CALL','OLDF')  !CLEAR DUMMY FILE FOR 2 SET 
	FACTOR = FAC
	CALL INIOPER(NFIL4,LP4,FACTOR,'COMB','OLDF')  !CLEAR DUMMY FILE FOR 2 SET 

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE CLEROUT
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------

	LP    = 0
	NSTEP = 0
	FACT  = 0.0D0 
	CALL INIOPER(LP,NSTEP,FACT,'CLER','OLDF')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LNUNIF(IDIR,W1,W2,AL,BL,ECT,VR,ELN,ANG,NP,RESUL,VE) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
C	=======================================================================
	DIMENSION VR(3),VS(3),VT(3),VE(3),VME(3)
	DIMENSION RESUL(6*NP),ECCMT(3,3)


	RANG = ANG

	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

	CALL CLEARA(RESUL,6*NP)
C	-----------------------------------------------------------
C	ECCENTRICITY EFFECT
C	-----------------------------------------------------------
	ECT   = 0.0D0
	ECS   = 0.0D0
	ECCMT = 0.0D0
	ECCMT(1,1) =  0.0
	ECCMT(1,2) = -ECT
	ECCMT(1,3) =  ECS
	ECCMT(2,1) =  ECT
	ECCMT(3,1) = -ECS
C	-----------------------------------------------------------
C	-----------------------------------------------------------
C	LOAD VALUE IN LOCAL SYSTEM
C	-----------------------------------------------------------
	W1R = VR(IDIR)*W1 !FORCE
	W1S = VS(IDIR)*W1
	W1T = VT(IDIR)*W1

	W1MR = ECCMT(1,1)*W1R + ECCMT(1,2)*W1S + ECCMT(1,3)*W1T !MOMENT
	W1MS = ECCMT(2,1)*W1R + ECCMT(2,2)*W1S + ECCMT(2,3)*W1T 
	W1MT = ECCMT(3,1)*W1R + ECCMT(3,2)*W1S + ECCMT(3,3)*W1T 

	W2R = VR(IDIR)*W2 !FORCE
	W2S = VS(IDIR)*W2
	W2T = VT(IDIR)*W2

	W2MR = ECCMT(1,1)*W2R + ECCMT(1,2)*W2S + ECCMT(1,3)*W2T !MOMENT
	W2MS = ECCMT(2,1)*W2R + ECCMT(2,2)*W2S + ECCMT(2,3)*W2T 
	W2MT = ECCMT(3,1)*W2R + ECCMT(3,2)*W2S + ECCMT(3,3)*W2T 
C	-----------------------------------------------------------

	DL = BL - AL

	X   = 0.0
	DS  = ELN/(NP+1)
	DO IP = 1,NP

	X = X + DS

	XM1 = VMCAULY(X,AL)
	XM2 = VMCAULY(X,BL)

	IF(AL.EQ.BL) GOTO 100 !FOR POINT LOAD


C	FROM WMR
	CALL TRAPIZ(XM1,DL,W1MR,W2MR,W1MR,DMR1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MR,W2MR,W2MR,DMR2,DUM1)
	RESUL(4+6*(IP-1)) = RESUL(4+6*(IP-1)) + DMR1-DMR2  !TORSION

C	FROM WMS
	CALL TRAPIZ(DL ,DL,W1MS,W2MS,W1MS,D1,DUM1)
	RT = -D1/ELN
	CALL TRAPIZ(XM1,DL,W1MS,W2MS,W1MS,D1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MS,W2MS,W2MS,D2,DUM1)
	RMS = RT*X + D1 - D2 
	RESUL(3+6*(IP-1)) = RESUL(3+6*(IP-1))       !SHEAR T
	RESUL(5+6*(IP-1)) = RESUL(5+6*(IP-1)) + RMS !BENDING S

C	FROM WMT
	CALL TRAPIZ(DL ,DL,W1MT,W2MT,W1MT,D1,DUM1)
	RS = -D1/ELN
	CALL TRAPIZ(XM1,DL,W1MT,W2MT,W1MT,D1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MT,W2MT,W2MT,D2,DUM1)
	RMT = RS*X + D1 - D2
	RESUL(2+6*(IP-1)) = RESUL(2+6*(IP-1))       !SHEAR S
	RESUL(6+6*(IP-1)) = RESUL(6+6*(IP-1)) - RMT !BENDING T



C	FROM WR
	CALL TRAPIZ(XM1,DL,W1R,W2R,W1R,RF1,DUM1)
	CALL TRAPIZ(XM2,DL,W1R,W2R,W2R,RF2,DUM1)
	RESUL(1+6*(IP-1)) = RESUL(1+6*(IP-1)) + RF1 - RF2  !AXIAL

C	FROM WS
	CALL TRAPIZ(DL ,DL,W1S,W2S,W1S,D1,DIS)
	D2 = D1*(ELN-AL-DIS)
	RS =-D2/ELN
	CALL TRAPIZ(XM1,DL,W1S,W2S,W1S,D1,DIS1)
	CALL TRAPIZ(XM2,DL,W1S,W2S,W2S,D2,DIS2)
	D3 = D1*(X  -AL-DIS1)
	D4 = D2*(X  -BL-DIS2)
	RMT = RS*X + D3 - D4
	RESUL(2+6*(IP-1)) = RESUL(2+6*(IP-1)) + D1 - D2		 !SHEAR S
	RESUL(6+6*(IP-1)) = RESUL(6+6*(IP-1)) + RMT          !BENDING T


C	FROM WT
	CALL TRAPIZ(DL ,DL,W1T,W2T,W1T,D1,DIS)
	D2 = D1*(ELN-AL-DIS)
	RT =-D2/ELN
	CALL TRAPIZ(XM1,DL,W1T,W2T,W1T,D1,DIS1)
	CALL TRAPIZ(XM2,DL,W1T,W2T,W2T,D2,DIS2)
	D3 = D1*(X  -AL-DIS1)
	D4 = D2*(X  -BL-DIS2)
	RMS = RT*X + D3 - D4
	RESUL(3+6*(IP-1)) = RESUL(3+6*(IP-1)) + D1 - D2		 !SHEAR T
	RESUL(5+6*(IP-1)) = RESUL(5+6*(IP-1)) - RMS          !BENDING S

	GOTO 200

100	CONTINUE

	DUM = 0.0
	IF(XM1.GT.0.0) DUM = 1.0 

C	FROM P MR
	RESUL(4+6*(IP-1)) = RESUL(4+6*(IP-1)) + W1MR*DUM  !TORSION

C	FROM P MS
	RT = -W1MS/ELN
	RMS = RT*X + W1MS*DUM
	RESUL(3+6*(IP-1)) = RESUL(3+6*(IP-1))       !SHEAR T
	RESUL(5+6*(IP-1)) = RESUL(5+6*(IP-1)) + RMS !BENDING S

C	FROM P MT
	RS = -W1MT/ELN
	RMT = RS*X + W1MT*DUM
	RESUL(2+6*(IP-1)) = RESUL(2+6*(IP-1))       !SHEAR S
	RESUL(6+6*(IP-1)) = RESUL(6+6*(IP-1)) - RMT !BENDING T


C	FROM P R
	RESUL(1+6*(IP-1)) = RESUL(1+6*(IP-1)) + W1R*DUM  !AXIAL

C	FROM P S
	D2 = W1S*(ELN-AL)
	RS =-D2/ELN
	RMT = RS*X   + W1S*XM1
	RESUL(2+6*(IP-1)) = RESUL(2+6*(IP-1)) + W1S*DUM       !SHEAR S
	RESUL(6+6*(IP-1)) = RESUL(6+6*(IP-1)) + RMT           !BENDING T

C	FROM P T
	D2 = W1T*(ELN-AL)
	RT =-D2/ELN
	RMS = RT*X   + W1T*XM1
	RESUL(3+6*(IP-1)) = RESUL(3+6*(IP-1)) + W1T*DUM       !SHEAR T
	RESUL(5+6*(IP-1)) = RESUL(5+6*(IP-1)) - RMS           !BENDING S


200	CONTINUE


	CALL VECPRD(VE,RESUL(1+6*(IP-1)),VME)                 !ACCOUNT FOR OFFSET
	DO I = 1,3
	RESUL(3+6*(IP-1)+I) = RESUL(3+6*(IP-1)+I) + VME(I) 
	ENDDO

	ENDDO
C	-----------------------------------------------------------	

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
