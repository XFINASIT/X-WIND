C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE INPLINK(NCON,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000) 
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000) 
C	-------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	-------------------------------------------------------------------
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
C	-------------------------------------------------------------------
	DIMENSION NCON(NSN,1)
	DIMENSION LEFLNK(NSF),IDSP(NSF,NSN)

	IF(IND.EQ.1) GOTO 100

	! LNKLEN = 31  !FOR CWR ANALYSIS
       LNKLEN = 30 !TOEY MODIFY 04/2015

	READ(ITI,*) 
	READ(ITI,*) NGLSET,NGLELM
	IF(NGLSET.EQ.0) RETURN


C	SET NOLIN TO BE 1  SO THAT IT WILL ACTIVATE AN ITERATION OPTION IN LINEAR ANALYSIS
	IF(ISOLOP.EQ.1) NOLIN = 1


	READ(ITI,*)
	DO ISET = 1,NGLSET
	NUM = LNKLEN*(ISET-1)
	READ(ITI,*) (DATLNK(J+NUM),J=1,LNKLEN)
	ENDDO

	READ(ITI,*)
	DO IELE = 1,NGLELM
	NUM = 3*(IELE-1)   
	READ(ITI,*) ID,(LEKNOD(J+NUM),J=1,3)   !ID  NOD-I NOD-J  ISET
	DO I = 1,2
	IN = LEKNOD(I+NUM)
	NCON(IN,1) = ID
	ENDDO
	ENDDO

	LNKWRK = 18+18

	DO IELE = 1,NGLELM
	NUM = LNKWRK*(IELE-1)
	DO J = 1,LNKWRK
	WATLNK(J+NUM) = 0.0D0
	ENDDO
	ENDDO


C	THIS BLOCK IS FOR BACKUP OF WORKING ARRAY
	DO IELE = 1,NGLELM
	NUM = LNKWRK*(IELE-1) + LNKWRK*NGLELM
	DO J = 1,LNKWRK
	WATLNK(J+NUM) = 0.0D0
	ENDDO
	ENDDO


	LNKFRE = 0
	DO I = 1,NSF
	IF(IDOF(I).LE.6) LNKFRE = LNKFRE + 1
	ENDDO

C     ------------------------------------------------------
      REWIND(8)
      READ(8) IDSP
      
      LEFLNK = 0
      II = 0
	DO I = 1,NSF
	    IF(IDOF(I).LE.6) THEN
	        II = II + 1
	        LEFLNK(II) = IDOF(I)
	    ENDIF
	ENDDO
      
      DO IELE = 1,NGLELM
      	  NUM = 3*(IELE-1)   
          DO INM = 1,2
              NOD = LEKNOD(INM+NUM)
              IF(NOD.LE.0) GOTO 300
              DO INF = 1,LNKFRE
                  ISF = 0
                  DO JSF = 1,NSF
                      IF(LEFLNK(INF).EQ.IDOF(JSF)) THEN
                          ISF = JSF
                          EXIT
                      ENDIF
                  ENDDO
                  IF(ISF.NE.0) IDSP(ISF,NOD) = 1
              ENDDO
300           CONTINUE
          ENDDO
      ENDDO
  
      REWIND(8)
      WRITE(8) IDSP
C     ------------------------------------------------------
      
      
	RETURN

100	CONTINUE
C	UPDATE STIFFNESS MATRIX COLUMN HEIGHT
	CALL LNKCOLH

C	CASE(18)
C	KDOF   = 0
C	NDFREE = 0
C	DO I = 1,NSF
C	IF(IDOF(I).LE.6) NDFREE = NDFREE + 1
C	IF(IDOF(I).LE.6) KDOF(NDFREE) = IDOF(I)
C	ENDDO
C	NFLINK = 0
C	DO I = 1,NDFREE
C	NN = NDFREE-I
C	NFLINK = NFLINK + KDOF(I)*10**NN
C	ENDDO


      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LINKELM(RE,REAC,SETIF,IPRN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000)
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000) 
C	-------------------------------------------------------------------
	DIMENSION RE(1),REAC(1),SETIF(1)

	
	DO IELE = 1,NGLELM
	CALL LINKTRN(DATLNK,WATLNK,LEKNOD,RE,REAC,SETIF,IELE)
	ENDDO

C	IF(IND.EQ.1) CALL OUTLINK
	IF(IPRN.EQ.0) CALL LIKOUT


      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LINKWNI
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000) 
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000)  
C	-------------------------------------------------------------------
C	INITIALIZE WORKING ARRAY FOR LINK ELEMENT

	
	DO IELE = 1,NGLELM
	NUM = LNKWRK*(IELE-1)
	DO J = 1,LNKWRK
	WATLNK(J+NUM) = 0.0D0
	ENDDO
	ENDDO


C	THIS BLOCK IS FOR BACKUP OF WORKING ARRAY
	DO IELE = 1,NGLELM
	NUM = LNKWRK*(IELE-1) + LNKWRK*NGLELM
	DO J = 1,LNKWRK
	WATLNK(J+NUM) = 0.0D0
	ENDDO
	ENDDO


      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LINKTRN(PROP,WA,LNKDT,RE,REAC,SETIF,IELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

	COMMON /REACT/ LRC,LRCT,MFQ,LRID

      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON A(9000000),IA(9000000)
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000) 
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000) 
      
      COMMON /CWROPT/ ICWR !FOR CWR ANALYSIS
      COMMON /LNKCWR/ IELECWR !FOR CWR ANALYSIS
C	-------------------------------------------------------------------
	DIMENSION PROP(LNKLEN,1),WA(LNKWRK,1),LNKDT(3,1)
C	-------------------------------------------------------------------
	DIMENSION RE(1),REAC(1)
	DIMENSION COORD(6),EDIS(12),LM(12),LMRCT(12)
	DIMENSION S(78),ELOD(12),FIN(12),SETIF(NEQ,LCS)
	
	IELECWR = IELE !FOR CWR ANALYSIS
      
      NN = 10**5
	IF(NSF.LE.6) NN = 10**(NSF-1)
	LINK = 0
	DO I = 1,NSF
	IF(IDOF(I).LE.6.AND.IDOF(I).NE.0) THEN
	LINK = LINK + NN*IDOF(I)
	NN = NN/10
	ENDIF
	ENDDO

	NNO = 2
	NNF = LNKFRE
	NEF = NNF*NNO


	NOD1 = LNKDT(1,IELE)
	NOD2 = LNKDT(2,IELE)
	ISET = LNKDT(3,IELE)

	EDIS(1:12) = 0.0D0


	CALL LINKERT (IA(LID),IA(LRID),A(LXY),A(LDT),A(LDI),COORD,EDIS,
	1			  LM,LMRCT,NOD1,NOD2)


      IF (NLS.GT.0) THEN
	CALL LIKRES (IA(LID),IA(LDS),A(LDC),LM,S,EDIS,
	1             ELOD,NNF,NNO,NEF,LINK,1)
	ENDIF


	IF(ITASK.NE.1) THEN
	CALL SETELM (EDIS,LMRCT,IA(LRID),NEF,LM,IA(LID))
	ENDIF


C	CALLING WORKING ARRAY FROM BACKUP
	WA(1:LNKWRK,IELE) = WA(1:LNKWRK,IELE+NGLELM)

      IF(ICWR.EQ.0) THEN
	CALL ELINK2(PROP(1,ISET),WA(1,IELE),S,COORD,EDIS,ELOD,FIN,
	1			NNF,NEF)
      ELSEIF(ICWR.EQ.1) THEN
      CALL ELINK_CWR(PROP(1,ISET),WA(1,IELE),S,COORD,EDIS,ELOD,FIN,
	1			NNF,NEF) ! Coupling Interface Element by ANAPHAT for CWR
      ENDIF

      IF (NLS.GT.0) THEN
	CALL LIKRES (IA(LID),IA(LDS),A(LDC),LM,S,EDIS,
	1             ELOD,NNF,NNO,NEF,LINK,2)
	ENDIF


C	============================================================
C	ASSEMBLE THE RESTRAINED STIFFNESS FOR SETTLEMENT LOAD CALC.
C	============================================================
	IF(IFREF.EQ.0) THEN
	IF(ITASK.LE.4) THEN
	CALL SETLODR(SETIF,IA(LID),IA(LRID),LM,LMRCT,S,NEF)
	ENDIF
	ENDIF
C	============================================================


	IF (ITASK.LE.2) GOTO 300
	
	IF(ITASK.EQ.5.AND.IFEIG.EQ.0) THEN
	  S(1:78) = 0.0D0
	  CALL ADSLINK (LM,IA(LMA),S,NEF)
	RETURN
	ENDIF
	
	IF(ITASK.EQ.6.AND.IFEIG.EQ.0) THEN
	  CALL ADSLINK (LM,IA(LMA),S,NEF)
	RETURN
	ENDIF
	
      IF (IFEIG.EQ.0) RETURN

300	CONTINUE

	DO 400  IEF=1,NEF
      IEQ = LM(IEF)
	IRC = LMRCT(IEF)
	IF (IRC.NE.0) REAC(IRC) = REAC(IRC) + ELOD(IEF)
400	IF (IEQ.NE.0) RE(IEQ)   = RE(IEQ)   - ELOD(IEF)



	IF(IFREF.EQ.0) THEN
	IF(ITASK.LE.4.OR.ITASK.EQ.6) THEN
	CALL ADSLINK (LM,IA(LMA),S,NEF)
	ENDIF
	ENDIF


C	BACKUP WORKING ARRAY FOR ITASK = 1 
	IF(ITASK.EQ.1) THEN
	WA(1:LNKWRK,IELE+NGLELM) = WA(1:LNKWRK,IELE)
	ENDIF


      !WRITE (300,'(E12.5)') S(1:78)
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LINKERT (ID,IDRCT,XYZ,DISP,DISPI,COORD,EDIS,LM,LMRCT,
	1					NOD1,NOD2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000) 
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000) 
C	-------------------------------------------------------------------
	DIMENSION ID(NSF,NSN),XYZ(NSN,1),DISP(1),DISPI(1),IDRCT(NSF,NSN)
C	-------------------------------------------------------------------
	DIMENSION COORD(3,2),EDIS(LNKFRE,2),LM(LNKFRE,2),LMRCT(LNKFRE,2)
C	-------------------------------------------------------------------

	EDIS = 0.0D0
	DO I = 1,LNKFRE
	IQ1 = ID(I,NOD1)
	IQ2 = ID(I,NOD2)
	LM(I,1) = IQ1
	LM(I,2) = IQ2
	IR1 = IDRCT(I,NOD1)
	IR2 = IDRCT(I,NOD2)
	LMRCT(I,1) = IR1
	LMRCT(I,2) = IR2
	IF(IQ1.NE.0) EDIS(I,1) = DISP(IQ1)
	IF(IQ2.NE.0) EDIS(I,2) = DISP(IQ2)
      IF (ITASK.NE.2) GOTO 200
      IF(IQ1.NE.0) EDIS(I,1)  = EDIS(I,1) + DISPI(IQ1)
      IF(IQ2.NE.0) EDIS(I,2)  = EDIS(I,2) + DISPI(IQ2)
 200  CONTINUE
	ENDDO


	DO I = 1,3
	COORD(I,1) = XYZ(NOD1,I)
	COORD(I,2) = XYZ(NOD2,I)
	ENDDO



      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE ELINK2(PROPM,WA,S,COORD,EDIS,RE,FIN,NNF,NEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     ------------------------------------------------------------- 
C
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      DIMENSION PROPM(1),EDIS(NEF)
      DIMENSION COORD(3,1),S(1)
	DIMENSION FIN(NEF),RE(NEF),WA(1)

	DIMENSION B(6,12),D(6,6),TRANS(12,12)
	DIMENSION S1(NEF,NEF),BB(6,NEF)
	DIMENSION VR(3),VS(3),VT(3)

	DIMENSION STRAIN(6),STRESS(6),KDOF(9)

	DIMENSION UVW(6),FOC(6),STF(6)
	DIMENSION ST1(6),ST2(6),FSY(6),GAP(6),HOK(6),EPT(6),DAP(6)
	DIMENSION LAT(6)



	KDOF = 0
	NN   = 0
	DO I = 1,NNF
	IF(IDOF(I).LE.6) NN = NN + 1
	IF(IDOF(I).LE.6) KDOF(NN) = IDOF(I)
	ENDDO

	JND		 = PROPM(2)  !FLAG FOR OPTION TYPE
	IND		 = PROPM(3)  !FLAG FOR LINEAR OR BILINEAR
	HOK(1:6) = PROPM(4)
	GAP(1:6) = PROPM(5)
	RANG     = PROPM(6)
	ST1(1:6) = PROPM(7 :12)
	ST2(1:6) = PROPM(13:18)
	FSY(1:6) = PROPM(19:24)
	DAP(1:6) = PROPM(25:30)


	NUMW = 0
	IF(ITASK.EQ.6) THEN
	IND = 0
	ST1(1:6) = DAP(1:6)
	NUMW = 18
	ENDIF

C	-------------------
C	INITIALIZATION 
C	-------------------
	B(1:6,1:12) = 0.0D0
	D(1:6,1:6 ) = 0.0D0

C	B MATRIX
	DO I = 1,6
	J = I + 6
	B(I,I) = -1.0D0
	B(I,J) =  1.0D0
	ENDDO

	VR = 0.0D0
C	LOCAL VECTOR
	DO I = 1,3
	VR(I) = COORD(I,2) - COORD(I,1) 
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)
	IF(ELN.EQ.0) VR(1:3) = [1.0D0,0.0D0,0.0D0]
	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

C	ACCOUNT FOR ELEMENT LENGTH
	B(3,5)  =  ELN*0.5
	B(3,11) =  ELN*0.5
	B(2,6)  = -ELN*0.5
	B(2,12) = -ELN*0.5


C	TRANSFORMATION MATRIX
	TRANS(1:12,1:12) = 0.0D0
	TRANS(1+0,1:3) = VR(1:3)
	TRANS(2+0,1:3) = VS(1:3)
	TRANS(3+0,1:3) = VT(1:3)
	TRANS(1+3,4:6) = VR(1:3)
	TRANS(2+3,4:6) = VS(1:3)
	TRANS(3+3,4:6) = VT(1:3)
	TRANS(1+6,7:9) = VR(1:3)
	TRANS(2+6,7:9) = VS(1:3)
	TRANS(3+6,7:9) = VT(1:3)
	TRANS(1+9,10:12) = VR(1:3)
	TRANS(2+9,10:12) = VS(1:3)
	TRANS(3+9,10:12) = VT(1:3)	

C	GLOBAL B MATRIX
	B = MATMUL(B,TRANS)

C	REDUCED B MATRIX
	DO 100 I = 1,2
	NN = NNF*(I-1)
	MM =   6*(I-1)
	DO 100 J = 1,NNF
	KK = KDOF(J)
	DO 100 K = 1,6
100	BB(K,J+NN) = B(K,KK+MM)

C	-------------------
C	STRAIN
C	-------------------
	STRAIN = MATMUL(BB,EDIS)

C	D MATRIX AND STRESS
	SELECT CASE(JND)
C	---------------------------
	CASE(1)
	LAT(1:6) = 0
	LAT(1:6) = 1  !ACTIVE ALL COMPONENT
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)
	CALL ETMOD1(FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1		    GAP,LAT,IND)
C	---------------------------
	CASE(2)
	LAT(1:6) = 0
	LAT(1)   = 1  !ACTIVE ONLY R COMPONENT
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)
	CALL ETMOD2(FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1			GAP,HOK,LAT,IND)
C	---------------------------
	CASE(3)
	LAT(1:6) = 0
	LAT(1)   = 1  !ACTIVE ONLY R COMPONENT
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)
	CALL ETMOD3(FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1			HOK,LAT,IND)
C	---------------------------
	CASE(4)
	LAT(1:6) = 0
	LAT(1)   = 1  !ACTIVE ONLY R COMPONENT
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)

	CALL ETMOD4(FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1		    GAP,LAT,IND)
C	---------------------------
	END SELECT
	
	J = NUMW
	WA(1+J:6+J) = FOC(1:6) 
	J = J + 6
	WA(1+J:6+J) = UVW(1:6)
	J = J + 6
	WA(1+J:6+J) = EPT(1:6)


C	RIGIDITY MATRIX
	DO I = 1,6
	D(I,I) = STF(I)
	ENDDO

C	ELEMENT STIFFNESS
	S1 = MATMUL(TRANSPOSE(BB),MATMUL(D,BB))

C	------------------------
C	FILL INTO UPPER-TRIANGLE
C	------------------------
	K = 0
	DO I=1,NEF
	DO J=I,NEF
	K = K + 1
	S(K) = S1(I,J)
	ENDDO
	ENDDO


	FIN = MATMUL(TRANSPOSE(BB),STRESS)  !RE(I)


C	RESISTING FORCE
	RE = MATMUL(TRANSPOSE(BB),STRESS)


C	STORE STRESS INTO WORKING ARRAY
	WA(1:6) = STRESS(1:6)



	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ELINK_CWR(PROPM,WA,S,COORD,EDIS,RE,FIN,NNF,NEF) ! Coupling Interface Element by ANAPHAT for CWR
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     ------------------------------------------------------------- 
C
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG,IFLAGCWR
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /LNKCWR/ IELECWR

      DIMENSION PROPM(1),EDIS(NEF)
      DIMENSION COORD(3,2),S(1)
	DIMENSION FIN(NEF),RE(NEF),WA(1)

	DIMENSION B(6,12),D(6,6),TRANS(12,12)
	DIMENSION S1(NEF,NEF),BB(6,NEF)
	DIMENSION VR(3),VS(3),VT(3)

	DIMENSION STRAIN(6),STRESS(6),KDOF(9)

	DIMENSION UVW(6),FOC(6),STF(6)
C	DIMENSION ST1(6),ST2(6),FSY(6),GAP(6),HOK(6),EPT(6),DAP(6)
	DIMENSION LAT(6)

C     CHANGE BY ANAPHAT CWR 07/2015
C     ST(1:6)     =   LINEAR CONDITION
C     ST1         =   UNDER CRITICAL REL DISP
C     ST2         =   OVER CRITICAL REL DISP
C     ST(13:18)   =   BI-LINEAR UNLOADED CONDITION
C     ST(7:12)    =   BI-LINEAR LOADED CONDITION
      
      DIMENSION ST1(18),ST2(18),FSY(18),GAP(6),HOK(6),EPT(6),DAP(6)
      
C     CHANGE BY ANAPAHT CWR      
c      if (IELECWR.EQ.270) then
c      write(*,*) 'T1',EDIS(2)
c      endif
      
c      WRITE(*,*) "TEST", ST1(1),ST1(2), ST2(1),ST2(2),FSY(1),FSY(2)
      
      
c      ST1(13) = 3600000000
c      ST1(14) = 72000000
c      ST1(14) = 36000000
c      ST1(15) = 0.0
c      ST1(16) = 0.0
c      ST1(17) = 0.0
c      ST1(18) = 0.0
      
c      ST2(13) = 3600000000
c      ST2(14) = 0.00001
c      ST2(15) = 0.0
c      ST2(16) = 0.0
c      ST2(17) = 0.0
c      ST2(18) = 0.0
      
c      FSY(13) = 3600000
c      FSY(14) = 250
c      FSY(14) = 18000
c      FSY(15) = 0.0
c      FSY(16) = 0.0
c      FSY(17) = 0.0
c      FSY(18) = 0.0
      
c      ST1(7) = 7200000000
c      ST1(8) = 72000000
c      ST1(9) = 0.0
c      ST1(10) = 0.0
c      ST1(11) = 0.0
c      ST1(12) = 0.0
      
c      ST2(7) = 7200000000
c      ST2(8) = 0.00001
c      ST2(9) = 0.0
c      ST2(10) = 0.0
c      ST2(11) = 0.0
c      ST2(12) = 0.0
      
c      FSY(7) = 7200000
c      FSY(8) = 36000
c      FSY(9) = 0.0
c      FSY(10) = 0.0
c      FSY(11) = 0.0
c      FSY(12) = 0.0
      
      READ(11114,*)
      READ(11114,*) (ST1(I),I=13,18)
      READ(11114,*) (ST2(I),I=13,18)
      READ(11114,*) (FSY(I),I=13,18)
      READ(11114,*)
      READ(11114,*) (ST1(I),I=7,12)
      READ(11114,*) (ST2(I),I=7,12)
      READ(11114,*) (FSY(I),I=7,12)

      REWIND (11114)
      
c      ST1(13) = 3000000000
c      ST1(14) = 30000000
c      ST1(15) = 0.0
c      ST1(16) = 0.0
c      ST1(17) = 0.0
c      ST1(18) = 0.0
      
c      ST2(13) = 3000000000
c      ST2(14) = 0.00001
c      ST2(15) = 0.0
c      ST2(16) = 0.0
c      ST2(17) = 0.0
c      ST2(18) = 0.0
      
c      FSY(13) = 3000000
c      FSY(14) = 50
c      FSY(15) = 0.0
c      FSY(16) = 0.0
c      FSY(17) = 0.0
c      FSY(18) = 0.0

	KDOF = 0
	NN   = 0
	DO I = 1,NNF
	IF(IDOF(I).LE.6) NN = NN + 1
	IF(IDOF(I).LE.6) KDOF(NN) = IDOF(I)
	ENDDO

	JND		 = PROPM(2)  !FLAG FOR OPTION TYPE
	IND		 = PROPM(3)  !FLAG FOR LINEAR OR BILINEAR
	HOK(1:6) = PROPM(4)
	GAP(1:6) = PROPM(5)
	RANG     = PROPM(6)
	ST1(1:6) = PROPM(7 :12)
	ST2(1:6) = PROPM(13:18)
	FSY(1:6) = PROPM(19:24)
	DAP(1:6) = PROPM(25:30)
      
      IFLAGCWR = PROPM(31)
      
c      if (IELECWR.EQ.270) then
c       WRITE(*,*) COORD(1:3,1)
c       WRITE(*,*) COORD(1:3,2)
c       pause
c      endif
      

c      ST1(1) = 7200000000
c      ST1(2) = 72000000
c      ST1(8) = 30000000
c      ST1(3) = 0.0
c      ST1(4) = 0.0
c      ST1(5) = 0.0
c      ST1(6) = 0.0
      
c      ST2(1) = 7200000000
c      ST2(2) = 0.00001
c      ST2(3) = 0.0
c      ST2(4) = 0.0
c      ST2(5) = 0.0
c      ST2(6) = 0.0
      
c      FSY(1) = 7200000
c      FSY(2) = 250
c      FSY(8) = 200
c      FSY(8) = 36000
c      FSY(3) = 0.0
c      FSY(4) = 0.0
c      FSY(5) = 0.0
c      FSY(6) = 0.0
      
c      WRITE(*,*) 'FLAGCWR', FLAGCWR2, FLAGCWR

	NUMW = 0
	IF(ITASK.EQ.6) THEN
	IND = 0
	ST1(1:6) = DAP(1:6)
	NUMW = 18
	ENDIF

C	-------------------
C	INITIALIZATION 
C	-------------------
	B(1:6,1:12) = 0.0D0
	D(1:6,1:6 ) = 0.0D0

C	B MATRIX
	DO I = 1,6
	J = I + 6
	B(I,I) = -1.0D0
	B(I,J) =  1.0D0
	ENDDO

	VR = 0.0D0
C	LOCAL VECTOR
	DO I = 1,3
	VR(I) = COORD(I,2) - COORD(I,1) 
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)
	IF(ELN.EQ.0) VR(1:3) = [1.0D0,0.0D0,0.0D0]
	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

C	ACCOUNT FOR ELEMENT LENGTH
c	B(3,5)  =  ELN*0.5
c	B(3,11) =  ELN*0.5
c	B(2,6)  = -ELN*0.5
c	B(2,12) = -ELN*0.5


C	TRANSFORMATION MATRIX
	TRANS(1:12,1:12) = 0.0D0
	TRANS(1+0,1:3) = VR(1:3)
	TRANS(2+0,1:3) = VS(1:3)
	TRANS(3+0,1:3) = VT(1:3)
	TRANS(1+3,4:6) = VR(1:3)
	TRANS(2+3,4:6) = VS(1:3)
	TRANS(3+3,4:6) = VT(1:3)
	TRANS(1+6,7:9) = VR(1:3)
	TRANS(2+6,7:9) = VS(1:3)
	TRANS(3+6,7:9) = VT(1:3)
	TRANS(1+9,10:12) = VR(1:3)
	TRANS(2+9,10:12) = VS(1:3)
	TRANS(3+9,10:12) = VT(1:3)	

C	GLOBAL B MATRIX
	B = MATMUL(B,TRANS)

C	REDUCED B MATRIX
	DO 100 I = 1,2
	NN = NNF*(I-1)
	MM =   6*(I-1)
	DO 100 J = 1,NNF
	KK = KDOF(J)
	DO 100 K = 1,6
100	BB(K,J+NN) = B(K,KK+MM)

C	-------------------
C	STRAIN
C	-------------------
	STRAIN = MATMUL(BB,EDIS)
      
c      WRITE(*,*) 'ST1-2', ST1(1)

C     CHECK BY ANAPHAT CWR 07/2015
c      IF (COORD(1,1).EQ.O) THEN
c          WRITE (*,*) "COORD X Y Z"
c          WRITE (*,*) COORD(1,1), COORD(2,1), COORD(3,1)
c          WRITE (*,*) STRAIN(1), ST1(1), ST2(1)
c      ENDIF
      
C	D MATRIX AND STRESS
	SELECT CASE(JND)
C	---------------------------
	CASE(1)
	LAT(1:6) = 0
	LAT(1:6) = 1  !ACTIVE ALL COMPONENT
c      write(*,*) 'LAT',LAT(1)
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)

	CALL ETMOD1_CWR(IFLAGCWR,COORD,FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1		    GAP,LAT,IND)
C	---------------------------
	CASE(2)
	LAT(1:6) = 0
	LAT(1)   = 1  !ACTIVE ONLY R COMPONENT
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)
	CALL ETMOD2(FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1			GAP,HOK,LAT,IND)
C	---------------------------
	CASE(3)
	LAT(1:6) = 0
	LAT(1)   = 1  !ACTIVE ONLY R COMPONENT
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)
	CALL ETMOD3(FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1			HOK,LAT,IND)
C	---------------------------
	CASE(4)
	LAT(1:6) = 0
	LAT(1)   = 1  !ACTIVE ONLY R COMPONENT
	J = NUMW
	FOC(1:6) = WA(1+J:6+J)
	J = J + 6
	UVW(1:6) = WA(1+J:6+J)
	J = J + 6
	EPT(1:6) = WA(1+J:6+J)

	CALL ETMOD4(FOC,UVW,STRAIN,STRESS,ST1,ST2,FSY,EPT,STF,
	1		    GAP,LAT,IND)
C	---------------------------
	END SELECT
	
	J = NUMW
	WA(1+J:6+J) = FOC(1:6) 
	J = J + 6
	WA(1+J:6+J) = UVW(1:6)
	J = J + 6
	WA(1+J:6+J) = EPT(1:6)


C	RIGIDITY MATRIX
	DO I = 1,6
	D(I,I) = STF(I)
	ENDDO

C	ELEMENT STIFFNESS
	S1 = MATMUL(TRANSPOSE(BB),MATMUL(D,BB))
      
C     CHECK BY ANAPHAT CWR 07/2015
c      IF (COORD(1,1).EQ.O) THEN
c          WRITE (*,*) "COORD X Y Z"
c          WRITE (*,*) COORD(1,1), COORD(2,1), COORD(3,1)
c          WRITE (*,*) "Local STIFF", STF(1)
c      ENDIF      

C	------------------------
C	FILL INTO UPPER-TRIANGLE
C	------------------------
	K = 0
	DO I=1,NEF
	DO J=I,NEF
	K = K + 1
	S(K) = S1(I,J)
	ENDDO
	ENDDO


	FIN = MATMUL(TRANSPOSE(BB),STRESS)  !RE(I)


C	RESISTING FORCE
	RE = MATMUL(TRANSPOSE(BB),STRESS)


C	STORE STRESS INTO WORKING ARRAY
	WA(1:6) = STRESS(1:6)
C     CHANGE BY ANAPAHT CWR 
      if (ITASK.NE.5.OR.ITASK.NE.6) then
          if (IELECWR.EQ.270.AND.RE(2).NE.0.0) then
c          write(11110,*) 'T2',EDIS(1),RE(1)
          endif
      endif


	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================

	SUBROUTINE ETMOD1(FOC,UVW,STN,STS,ST1,ST2,FSY,EPT,STF,
	1				  GAP,LAT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	GENERAL
	DIMENSION STN(6),STS(6)
	DIMENSION FOC(6),ST1(6),ST2(6),FSY(6),EPT(6),UVW(6)
	DIMENSION STF(6),GAP(6),LAT(6)

	STS(1:6) = 0.0D0
	STF(1:6) = 0.0D0

	DO I = 1,6
	IF(LAT(I).EQ.0  ) GOTO 100

	CALL LIKFOC(UVW(I),FOC(I),STN(I),STS(I),ST1(I),ST2(I),
	1			FSY(I),STF(I),EPT(I),IND) 

100	CONTINUE
	ENDDO

      RETURN
      END
C
C	===============================================================
C	===============================================================
C	===============================================================
	SUBROUTINE ETMOD2(FOC,UVW,STN,STS,ST1,ST2,FSY,EPT,STF,
	1				  GAP,HOK,LAT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	HOOK&GAP
	DIMENSION STN(6),STS(6)
	DIMENSION FOC(6),ST1(6),ST2(6),FSY(6),EPT(6),UVW(6)
	DIMENSION STF(6),GAP(6),HOK(6),LAT(6)

	STS(1:6) = 0.0D0
	STF(1:6) = 0.0D0

	DO I = 1,6
	IF(LAT(I).EQ.0  ) GOTO 200

	IF(STN(I).GT.0.0) GOTO 100
	GP = STN(I) + GAP(I)
	IF(GP.GT.0.0) GOTO 100

	CALL LIKFOC(UVW(I),FOC(I),GP,STS(I),ST1(I),ST2(I),
	1			FSY(I),STF(I),EPT(I),IND) 

100	CONTINUE

	IF(STN(I).LT.0.0) GOTO 200
	HK = STN(I) - HOK(I)
	IF(HK.LT.0.0) GOTO 200

	CALL LIKFOC(UVW(I),FOC(I),HK,STS(I),ST1(I),ST2(I),
	1			FSY(I),STF(I),EPT(I),IND) 

200	CONTINUE

	ENDDO

      RETURN
      END
C
C	===============================================================
C	===============================================================
C	===============================================================
	SUBROUTINE ETMOD3(FOC,UVW,STN,STS,ST1,ST2,FSY,EPT,STF,
	1				  HOK,LAT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	HOOK
	DIMENSION STN(6),STS(6)
	DIMENSION FOC(6),ST1(6),ST2(6),FSY(6),EPT(6),UVW(6)
	DIMENSION STF(6),HOK(6),LAT(6)

	STS(1:6) = 0.0D0
	STF(1:6) = 0.0D0

	DO I = 1,6
	IF(LAT(I).EQ.0  ) GOTO 100
	IF(STN(I).LT.0.0) GOTO 100
	HK = STN(I) - HOK(I)
	IF(HK.LT.0.0) GOTO 100

	CALL LIKFOC(UVW(I),FOC(I),HK,STS(I),ST1(I),ST2(I),
	1			FSY(I),STF(I),EPT(I),IND) 

100	CONTINUE
	ENDDO


      RETURN
      END
C
C	===============================================================
C	===============================================================
C	===============================================================
	SUBROUTINE ETMOD4(FOC,UVW,STN,STS,ST1,ST2,FSY,EPT,STF,
	1				  GAP,LAT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	GAP
	DIMENSION STN(6),STS(6)
	DIMENSION FOC(6),ST1(6),ST2(6),FSY(6),EPT(6),UVW(6)
	DIMENSION STF(6),GAP(6),LAT(6)

	STS(1:6) = 0.0D0
	STF(1:6) = 0.0D0

	DO I = 1,6
	IF(LAT(I).EQ.0  ) GOTO 100
	IF(STN(I).GT.0.0) GOTO 100
	GP = STN(I) + GAP(I)
	IF(GP.GT.0.0) GOTO 100

	CALL LIKFOC(UVW(I),FOC(I),GP,STS(I),ST1(I),ST2(I),
	1			FSY(I),STF(I),EPT(I),IND) 

100	CONTINUE
	ENDDO

      RETURN
      END
C
C	===============================================================
C	===============================================================
C	===============================================================
	SUBROUTINE ETMOD1_CWR(IFLAGCWR,COORD,FOC,UVW,STN,STS,ST1,ST2,FSY,EPT,STF,
	1				  GAP,LAT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	GENERAL
	DIMENSION STN(6),STS(6)
C	DIMENSION FOC(6),ST1(6),ST2(6),FSY(6),EPT(6),UVW(6)
	DIMENSION STF(6),GAP(6),LAT(6)
      
C     CHANGE BY ANAPHAT 07/2015 ..........
      DIMENSION FOC(6),ST1(18),ST2(18),FSY(18),EPT(6),UVW(6),COORD(3,2)
      
      DUMMY = STN(2)
      
C ........................................      

	STS(1:6) = 0.0D0
	STF(1:6) = 0.0D0

	DO I = 1,6

	IF(LAT(I).EQ.0  ) GOTO 100
      
C     CHANGE BY ANAPHAT 07/2015 FOR BI-LINEAR FIX IN X-COORDINATE
      J = I
c      WRITE(*,*) 'STN', ST1(1:2), STN(1:6), FSY(2)

c .....BiLinear Case.....
      IF (IFLAGCWR.EQ.2) THEN
          IF(DUMMY.GT.0.00000000000001) THEN
              J = I+12
          ENDIF
      
          IF(DUMMY.EQ.0.00000000000001) THEN
              J = I+6
          ENDIF
      
          IF(DUMMY.LT.0.00000000000001) THEN
              J = I+6
          ENDIF
      ENDIF
      
	CALL LIKFOCWR(UVW(I),FOC(I),STN(I),STS(I),ST1(J),ST2(J),
	1			FSY(J),STF(I),EPT(I),IND) 


100	CONTINUE
	ENDDO

      RETURN
      END
C
C	===============================================================
C	===============================================================
C	===============================================================      
	SUBROUTINE LIKFOC(EPSP,SIGP,EPS,SIG,YOUNG,YOUN2,YIELD,STIFF,
	1				  EPTTN,IND) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	BILINEAR KINEMATIC HARDENING MATERIAL
C	=============================================================
	HARDS = YOUNG*YOUN2/(YOUNG-YOUN2)
	DEPS  = EPS - EPSP
	DSIG  = DEPS*YOUNG

	SIGO = SIGP

	SIG  = SIGO + DSIG 

	PREY  = YIELD 

	IF(IND.EQ.0) GOTO 40  !FLAG FOR LINEAR

	IF(ABS(SIGO - HARDS*EPTTN).GE.PREY) GO TO 20
	ESCUR = ABS(SIG - HARDS*EPTTN) - PREY
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(DSIG)
	GO TO 30


20	IF(SIG - HARDS*EPTTN.GT.0.0.AND.DSIG.LT.0.0) GO TO 40
	IF(SIG - HARDS*EPTTN.LT.0.0.AND.DSIG.GT.0.0) GO TO 40

	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	SIGO  = SIGO + REDUC*DSIG + 
	1	    RFACT*DSIG*(1.0-YOUNG/(YOUNG+HARDS))
	EPTTN = EPTTN + RFACT*DSIG/(YOUNG+HARDS)
	STIFF = YOUN2

	GO TO 50	

40	STIFF = YOUNG
	SIGO  = SIGO + DSIG
50	SIGP  = SIGO

	EPSP = EPS

	SIG  = SIGP

	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE LIKFOCWR(EPSP,SIGP,EPS,SIG,YOUNG,YOUN2,YIELD,STIFF,
	1				  EPTTN,IND) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	BILINEAR KINEMATIC HARDENING MATERIAL BY ANAPHAT
C	=============================================================
	
      
      IF(YOUNG.EQ.0.0D0) THEN
	HARDS = 0.0D0
	ELSE
	HARDS = YOUNG*YOUN2/(YOUNG-YOUN2)
      ENDIF
	 
	DEPS  = EPS - EPSP
	DSIG  = DEPS*YOUNG

	SIGO = SIGP

	SIG  = SIGO + DSIG 

	PREY  = YIELD 

	IF(IND.EQ.0) GOTO 40  !FLAG FOR LINEAR
	IF(YOUNG.EQ.YOUN2) GOTO 40  !FLAG FOR LINEAR

	IF(ABS(SIGO - HARDS*EPTTN).GE.PREY) GO TO 20
	ESCUR = ABS(SIG - HARDS*EPTTN) - PREY
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(DSIG)
	GO TO 30


20	IF(SIG - HARDS*EPTTN.GT.0.0.AND.DSIG.LT.0.0) GO TO 40
	IF(SIG - HARDS*EPTTN.LT.0.0.AND.DSIG.GT.0.0) GO TO 40

	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	SIGO  = SIGO + REDUC*DSIG + 
	1	    RFACT*DSIG*(1.0-YOUNG/(YOUNG+HARDS))
	EPTTN = EPTTN + RFACT*DSIG/(YOUNG+HARDS)
	STIFF = YOUN2

	GO TO 50	

40	STIFF = YOUNG
	SIGO  = SIGO + DSIG
50	SIGP  = SIGO

	EPSP = EPS

	SIG  = SIGP

	RETURN

      END
      
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE LNKCOLH
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON A(9000000),IA(9000000)
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000) 
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000) 
C	-------------------------------------------------------------------

	DO IELE = 1,NGLELM
	CALL LKCOLHI (IA(LID),IA(LCH),LEKNOD,IELE)
	ENDDO

      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LKCOLHI(ID,MHT,LNKDT,IELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000) 
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000) 
C	-------------------------------------------------------------------
	DIMENSION ID(NSF,NSN),LNKDT(3,1),MHT(1)
C	-------------------------------------------------------------------
	DIMENSION LM(12)

	NNF = LNKFRE
	NEF = NNF*2

	NOD1 = LNKDT(1,IELE)
	NOD2 = LNKDT(2,IELE)
	ISET = LNKDT(3,IELE)

	DO I = 1,LNKFRE
	LM(I)	  = ID(I,NOD1)
	LM(I+NNF) = ID(I,NOD2)
	ENDDO

      MEQ = 0
      DO 390  IEF=1,NEF
      IEQ = LM(IEF)
      IF (IEQ)     310,390,310
 310  IF (MEQ.EQ.0) MEQ = IEQ
      IF (MEQ-IEQ) 320,390,390
 320  MEQ = IEQ
 390  CONTINUE
 
      DO 400  IEF=1,NEF
      IEQ = LM(IEF)
      IF (IEQ.EQ.0) GOTO 400
      KHT = MEQ-IEQ
      IF (KHT.GT.MHT(IEQ)) MHT(IEQ) = KHT
 400	CONTINUE

      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE ADSLINK (LM,MAXA,S,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ASSEMBLES UPPER TRIANGULAR ELEMENT STIFFNESS INTO COMPACTED
C     GLOBAL STIFFNESS BLOCK (CONTRIBUTIONS BETWEEN NEQF AND NEQL)
C	------------------------------------------------------------
C     LM(NEF)       = EQUATION NUMBERS FOR ELEMENT DEGREE OF FREEDOMS
C     MAXA(NEQ1)    = ADDRESSES OF DIAGONAL ELEMENTS IN A
C     A(ISTOR)      = GLOBAL COMPACTED STIFFNESS BLOCK
C     S(NWS)        = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     NEQF,NEQL     = FIRST AND LAST EQUATION CONTAINED IN BLOCK
C     NPRE          = NUMBER OF PREVIOUS ELEMENTS IN A
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ----------------------------------------------------------------
C
      DIMENSION LM(1),MAXA(1),S(1)

	CALL MESTIF(LM,S,0,NEF,'WRT')
	RETURN
C
      NDI = 0
      DO 200 IEF=1,NEF
      II = LM(IEF)
      MI = MAXA(II)
      KS = IEF
      DO 220  JEF=1,NEF
      JJ = LM(JEF)
      IF (JJ) 220,220,110
 110  IJ = II-JJ
      IF (IJ) 220,210,210
 210  KK = MI+IJ
      KSS = KS
      IF (JEF.GE.IEF) KSS = JEF+NDI
C	A(KK) = A(KK)+S(KSS)
 220  KS = KS+NEF-JEF
 200  NDI = NDI+NEF-IEF

C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE OUTLINK
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------------
      COMMON /MGLENK/ NGLELM,NGLSET,LNKLEN,LNKWRK,LNKFRE,LEKNOD(5000) 
      COMMON /DTLENK/ DATLNK(50000),WATLNK(50000) 
C	-------------------------------------------------------------------	
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET


	DIMENSION SIG(6)

	IF(NGLELM.EQ.0) RETURN

	IF(NSTEP.EQ.1) IGIDSTEP = KSTEP/NPRIN
	IF(NSTEP.GT.1) IGIDSTEP = KSTEP-1/NPRIN

	
	WRITE (101,1000) IGIDSTEP
	DO IELE = 1,NGLELM
	NUM = LNKWRK*(IELE-1)
	DO J = 1,6
	SIG(J) = WATLNK(J+NUM)
	ENDDO
	NUM = 3*(IELE-1)   
	IN1 = LEKNOD(1+NUM)
	IN2 = LEKNOD(2+NUM)
	WRITE (101,1500) IN1,SIG(1:3)
	WRITE (101,1500) IN2,SIG(1:3)
	ENDDO
	WRITE (101,2000)
	
	WRITE (101,1100) IGIDSTEP
	DO IELE = 1,NGLELM
	NUM = LNKWRK*(IELE-1)
	DO J = 1,6
	SIG(J) = WATLNK(J+NUM)
	ENDDO
	NUM = 3*(IELE-1)   
	IN1 = LEKNOD(1+NUM)
	IN2 = LEKNOD(2+NUM)
	WRITE (101,1500) IN1,SIG(4:6)
	WRITE (101,1500) IN2,SIG(4:6)
	ENDDO
	WRITE (101,2000)



1000	FORMAT (//'Result "Link Force"',2x,'"XFinas"',2x,i5,5x,
     1 'Vector',2x,'OnNodes'/
	2 'ComponentNames "F-R", "F-S", "F-T"'/'Values')

1100	FORMAT ('Result "Link Moment"',2x,'"XFinas"',2x,i5,5x,
     1 'Vector',2x,'OnNodes'/
	2 'ComponentNames "M-R", "M-S", "M-T"'/'Values')

1500	FORMAT(X,I5,2X,3E17.6)

2000	FORMAT('End Values'//)


	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE LIKRES (ID,IDSET,DIRCOS,LM,ESTI,EDIS,ELOD,
     1                   MNF,NNO,NEF,LKLINK,MODE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIt INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     SCANS NODES TO DETERMINE WHICH HAVE LOCAL (SKEWED) COORDINATES
C     AND CALLS ROUTINES TO PERFORM APPROPRIATE TRANSFORMATIONS
C	---------------------------------------------------------
C     MODE = 1  TRANSFORMS ELEMENT DISPLACEMENT VECTOR INTO GLOBAL
C     MODE = 2  TRANSFORMS ELEMENT LOAD VECTOR INTO LOCAL  AND
C        IF IFREF = 0  TRANSFORMS ELEMENT STIFFNESS MATRIX INTO GLOBAL
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON A(9000000),IA(9000000)

      DIMENSION ID(NSF,1),IDSET(1),DIRCOS(9,1),LM(MNF,1),ESTI(1)
      DIMENSION EDIS(1),ELOD(1)

	DIMENSION COSN(81),TCOS(81),S(81),ST(81)

      DIMENSION IGPOS(9),IDPOS(9)
C
C
      IF (ITASK.EQ.5) RETURN
	IF ((NOLIN.EQ.0).AND.(ITASK.NE.3).AND.(MODE.EQ.1)) RETURN
C
	NNF = MNF
      LNF = MNF
C     ----------------------------------------------------
C     IGPOS = GLOBAL D.O.F. NOS. POSSIBLE FOR THIS ELEMENT
C     IDPOS = POSITIONS OF THESE D.O.F.'S IN ID
C     ----------------------------------------------------
  20  LINK = LKLINK
      IGF = 0
      DO 300  INF = 1,NNF
      IGPOS(INF) = LINK/10**(NNF-INF)
      IF(IGPOS(INF).GT.6) LNF = LNF-1
 320  IGF = IGF+1
	IF(IGF.GT.9) GOTO 301
      IF (IDOF(IGF)-IGPOS(INF)) 320,340,340
 340  LINK = LINK-10**(NNF-INF)*IGPOS(INF)
 300  IDPOS(INF) = IGF
 301	CONTINUE

      LINK = LKLINK
      DO INF=1,NNF
      IGPOS(INF) = LINK/10**(NNF-INF)
	LINK = LINK-10**(NNF-INF)*IGPOS(INF)
	ENDDO
C     ---------------------------------------------------
C     SCAN LM FOR FIRST EQN. NO. FOR EACH NODE OF ELEMENT
C     ---------------------------------------------------
      DO 100  INO = 1,NNO
      DO 120  INF = 1,NNF
      IF (LM(INF,INO).EQ.0) GOTO 120
      IEQ   = LM(INF,INO)
	IGF   = IGPOS(INF)
	IDROW = IDOFCALL(IDOF,IGF)
      GOTO 140
 120  CONTINUE
	GOTO 100
C     --------------------------------
C     SCAN ID FOR AXES SET NO. OF NODE
C     --------------------------------
 140  DO 200  ISN = 1,NSN
      IF (ID(IDROW,ISN).EQ.IEQ) GOTO 220
 200  CONTINUE
 220  ISET = IDSET(ISN)
      IF (ISET.EQ.0) GOTO 100
C     -------------------------------------------
C     EXTRACT APPROPRIATE DIRECTION COSINE MATRIX
C     -------------------------------------------
      CALL COSMAT (DIRCOS(1,ISET),COSN,IGPOS,LNF)
C
      NADD = (INO-1)*MNF+1
      IF (MODE-2) 400,500,100
C     ----------------------------------------------------------
C     TRANSFORM ELEMENT DISPLACEMENT VECTOR INTO GLOBAL (MODE=1)
C     ----------------------------------------------------------
 400  CALL MATRAN (COSN,TCOS,LNF,LNF,1)
      CALL TREVEC (EDIS(NADD),TCOS,LNF)
      GOTO 100
C     ------------------------------------------------------
C     TRANSFORM ELEMENT LOAD VECTOR INTO GLOBAL AND
C     ELEMENT STIFFNESS OR INITIAL STRESS MATRIX IF REQUIRED
C     ------------------------------------------------------
 500  CONTINUE
	CALL TREVEC (ELOD(NADD),COSN,LNF)  !IF (NOLIN.NE.0) CALL TREVEC (ELOD(NADD),COSN,LNF)   CHANGED BY SONGSAK JUL2006 TO ACTIVATE THE REACTION
      IF (ITASK.GE.3) GOTO 520
      IF (IFREF.NE.0) GOTO 100
      GOTO 540
 520  IF (ITASK.EQ.3 .AND. IFEIG.NE.0) GOTO 100
 540  CALL MATRAN (COSN,TCOS,LNF,LNF,1)
      CALL TRESTI (ESTI,COSN,TCOS,MNF,LNF,NNO,INO,S,ST)

	IF(LSYMM.EQ.1) THEN
	NEF2 = (NEF*NEF+NEF)/2 + 1 
	CALL TRESTI (ESTI(NEF2),COSN,TCOS,MNF,LNF,NNO,INO,S,ST)
	ENDIF

      GOTO 100
C     ------------------------------------

C
 100  CONTINUE
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================