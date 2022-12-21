C	=======================================================================
C	=======================FRAME THERMAL LOAD==============================
C	=======================================================================
	SUBROUTINE TEMPFRM(LM,XYZ,IGSET,PROPG,MTSET,PROPM,NTMPLD,IGIDM) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	=======================================================================
C	LAXROT = 1  !ROTATE ABOUT S AXIS
C	LAXROT = 2  !ROTATE ABOUT T AXIS
C	=======================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

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

	COMMON /LCSS/ ILCN,ILCC
C	==================================================================
	COMMON A(9000000),IA(9000000)

      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION VR(3),VS(3),VT(3),TRANS(7,7)
	DIMENSION RG(7,3),XYZ(NCO*NNM,NELE),LM(NEF,NELE)
	DIMENSION BA(14),BB(14),FIXD(14),FIXG(14),IGIDM(NELE)
	DIMENSION BWG(10),BPG(10),BBX(14)
	DIMENSION AA(10)

      ALLOCATE(RAL(NEF),MAL(NEF))

	DO ITEMP = 1,NTMPLD

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	READ(ITI,*) MLE,LAXROT,LOPTMP,NUMPAT,CONS,DIFS,DTEMP,TH,ILCN,ILCC

C     AXIS OF ROTATION OPTION       LAXROT   1 = ROTATE ABOUT S AXIS, 2 = ROTATE ABOUT T AXIS
C     TEMPERATURE PATTERN OPTION    LOPTMP   0 = STANDARD INPUT, 1 = ADVANCE INPUT
C     TEMPERATURE PATTERN NO.       NUMPAT
C     CONSTANT TEMPERATURE CHANGE   CONS
C     DIFF.    TEMPERATURE CHANGE   DIFS
C     NORMINAL DIMENSION            TH
      MLEGID = MLE
      
	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)


	CALL XFSECTION(KEG,MLE,1)

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FTYP ,1,1 ,0) !SECTION TYPE	0=STANDARD 1=FIBER
	CALL MRELFIL(NAMEI,FGAS ,1,3 ,0) !NUMBER OF STATION POINT ALONG ELEMENT
	CALL MRELFIL(NAMEI,FFIB ,1,4 ,0) !NUMBER OF FIBER
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !ROTATION ANGLE
	ITYP = INT(FTYP)
	NGAS = INT(FGAS)
	NFIB = INT(FFIB)

	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

C	----------------------------------------------------------
	FIXD = 0.0D0
C	----------------------------------------------------------
C	CALLING ELEMENT GAUSS DATA
	CALL FRGAUS (NGAS,ELN,BPG,BWG)

C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET STRESS CONT TO ELEMENT FORCE VECTOR
C     ----------------------------------------------------------
      DO 500 IPT = 1,NGAS

C	CALL SECTION PROPERTIES
	CALL XFSECTION(KEG,MLE,IPT)
	
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
C	==========================================================
	CALL MRELFIL(NAMEI,BNORM,1,111,0) !BNORM NORMINAL WIDTH
	CALL MRELFIL(NAMEI,HNORM,1,112,0) !HNORM NORMINAL HIGHT
	CALL MRELFIL(NAMEI,SMAX ,1,113,0) !
	CALL MRELFIL(NAMEI,SMIN ,1,114,0) !
	CALL MRELFIL(NAMEI,TMAX ,1,115,0) !
	CALL MRELFIL(NAMEI,TMIN ,1,116,0) !
C	==========================================================
	CALL MRELFIL(NAMEI,EAREA,1,41,0) !Equivalent AREA-E-ALPHA
	CALL MRELFIL(NAMEI,EQS  ,1,42,0) !Equivalent  QS-E-ALPHA
	CALL MRELFIL(NAMEI,EQT  ,1,43,0) !Equivalent  QT-E-ALPHA
	CALL MRELFIL(NAMEI,ESIT ,1,44,0) !Equivalent ISS-E-ALPHA
	CALL MRELFIL(NAMEI,ETIT ,1,45,0) !Equivalent ITT-E-ALPHA
C	==========================================================

C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
	BXD = BPG(IPT)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
	CALL BBXFRM(BBX,ELN,BXD)

C	----------------------
	IF(TH.EQ.0.0D0) THEN
	    IF(LAXROT.EQ.1) THEN  !ROTATE ABOUT S
	      TH = BNORM
	    ENDIF
	    IF(LAXROT.EQ.2) THEN  !ROTATE ABOUT T
	      TH = HNORM
	    ENDIF
	ENDIF
C	----------------------
	IF(TH.EQ.0.0D0) THEN
	    FIXD = 0.0D0
	    GOTO 1000
	ENDIF
C	----------------------

	ACON = CONS           !0.5*(TPR+TPL)
	BCON = -DIFS/TH       !(TPR-TPL)/TH

C	----------------------
	BA = 0.0 !B-MATRIX FOR AXIAL
	BB = 0.0 !B-MATRIX FOR BENDING
      IF(LAXROT.EQ.1) THEN  !ROTATE ABOUT S
	    BA(1) =  BBX(1)		!EDIS(1)    !SEE ALSO FRAMENL
	    BA(8) =  BBX(8)		!EDIS(8)
	    BB(3) =  BBX(3)		!EDIS(3)
	    BB(5) =  BBX(4)		!EDIS(5)
	    BB(10)=  BBX(10)	    !EDIS(10)
	    BB(12)=  BBX(11)	    !EDIS(12)
      ENDIF
      IF(LAXROT.EQ.2) THEN  !ROTATE ABOUT T
	    BA(1) =  BBX(1)		!EDIS(1)    !SEE ALSO FRAMENL
	    BA(8) =  BBX(8)		!EDIS(8)
	    BB(2) =  BBX(2)		!EDIS(2)
	    BB(6) =  BBX(5)		!EDIS(6)
	    BB(9) =  BBX(9)		!EDIS(9)
	    BB(13)=  BBX(12)	    !EDIS(13)
      ENDIF	
C	----------------------

C	====================================================================
      IF(ITYP.EQ.0) THEN !STANDARD SECTION
      
        SELECTCASE(LOPTMP) !LOPTMP   0 = STANDARD INPUT, 1 = ADVANCE INPUT
        
            CASE(0) !0 = STANDARD INPUT --- INPUT CONSTANT AND DIFFERENTIAL TEMP.   
              IF(LAXROT.EQ.1) THEN  !ROTATE ABOUT S
	              VP = ACON*EAREA + BCON*EQS
	              VM = ACON*EQS   + BCON*ESIT
              ENDIF
              IF(LAXROT.EQ.2) THEN  !ROTATE ABOUT T
	              VP = ACON*EAREA + BCON*EQT
	              VM = ACON*EQT   + BCON*ETIT
              ENDIF	 
              
            CASE(1) !1 = ADVANCE INPUT --- INPUT TEMPERATURE PATTERN OVER SECTION
              WRITE(*,*) '!!WARNING ---- IN FRAME ELEMENT NUMBER:',MLEGID
              WRITE(*,*) '     TEMP. PATTERN CAN NOT BE APPLIED TO FRAME WITH STANDARD OR GENERAL SECTION' 
              GOTO 1000    
        ENDSELECT      
      ENDIF
C	====================================================================



C	====================================================================
      IF(ITYP.EQ.1) THEN !STANDARD SECTION
      
        SELECTCASE(LOPTMP) !LOPTMP   0 = STANDARD INPUT, 1 = ADVANCE INPUT
        
        CASE(0) !0 = STANDARD INPUT --- INPUT CONSTANT AND DIFFERENTIAL TEMP.   
          IF(LAXROT.EQ.1) THEN  !ROTATE ABOUT S
	          VP = ACON*EAREA + BCON*EQS
	          VM = ACON*EQS   + BCON*ESIT
          ENDIF
          IF(LAXROT.EQ.2) THEN  !ROTATE ABOUT T
	          VP = ACON*EAREA + BCON*EQT
	          VM = ACON*EQT   + BCON*ETIT
          ENDIF	 
          
        CASE(1) !1 = ADVANCE INPUT --- INPUT TEMPERATURE PATTERN OVER SECTION
C	---------------------------------------------------------

            ISET = IGSET(MLE) !GEOMETRIC SET
	      III = 1
	      INAME(1:4) = [5,1,ISET,KEG] !XSEH
	      CALL ICONC(INAME,NAMEI)
	      DO II = 1,7
	          CALL MRELFIL(NAMEI,AA(II),III,II,0) !Store Section Props 
	      ENDDO
	      NSTS  = INT(AA(2))
	      NTOUR = INT(AA(6))
	      NSEGT = INT(AA(7))

            VP = 0.0D0 ; VM = 0.0D0
C           INTEGRATE OVER SECTION
	      DO 300 IFIB = 1,NFIB
C	---------------------------------------------------------

	          III = 1 + NSEGT + (1+NTOUR+NFIB+NSTS)*(IPT-1) + 1 + NTOUR + IFIB
	          DO II = 1,7
	              CALL MRELFIL(NAMEI,AA(II),III,II,0) !Store Section Props 
	          ENDDO

	          DA = AA(1) ; SC = AA(2) ; TC = AA(3) ; FMT = AA(4) ; WP = AA(6) ; TR = AA(7)

                SS = SC - SMIN ; TT = TC - TMIN !DISTANCE FROM BOTTOM OF SECTION
                IF(LAXROT.EQ.1) THEN  !ROTATE ABOUT S
	              CALL VALTEMPPAT(NUMPAT,TH,TT,FTEMP,'TMPT')
                ENDIF
                IF(LAXROT.EQ.2) THEN  !ROTATE ABOUT T
	              CALL VALTEMPPAT(NUMPAT,TH,SS,FTEMP,'TMPT')
                ENDIF	
          
	          MSET = INT(FMT)
	          YOUNG = PROPM(1 ,MSET)
	          ALPHA = PROPM(13,MSET)
	
                VP = VP + FTEMP*DTEMP*ALPHA*YOUNG*DA
                IF(LAXROT.EQ.1) THEN  !ROTATE ABOUT S
                    VM = VM - TT*FTEMP*DTEMP*ALPHA*YOUNG*DA
                ENDIF
                IF(LAXROT.EQ.2) THEN  !ROTATE ABOUT T
                    VM = VM + SS*FTEMP*DTEMP*ALPHA*YOUNG*DA*-1.0D0 !MULTIPLY -1 here to make positive deflection
                ENDIF	
      	
C	---------------------------------------------------------
300         CONTINUE
C	--------------------------------------------------------- 
        ENDSELECT          
      ENDIF
C	====================================================================


C     INTEGRATE OVER LENGTH
      DO I = 1,14
        FIXD(I) = FIXD(I) + (BA(I)*VP + BB(I)*VM)*BWG(IPT) 
      ENDDO

      
500	CONTINUE
C	----------------------------------------------------------
C	----------------------------------------------------------
C	TRANSFORMATION
1000	TRANS = 0.0
	DO I = 1,3
	    TRANS(I,1)   = VR(I)
	    TRANS(I,2)   = VS(I)
	    TRANS(I,3)   = VT(I)
	    TRANS(I+3,4) = VR(I)
	    TRANS(I+3,5) = VS(I)
	    TRANS(I+3,6) = VT(I)
	ENDDO
	FIXG = 0.0
	DO I = 1,6
	    DO J = 1,6
	      FIXG(I) = FIXG(I) + TRANS(I,J)*FIXD(J)
	    ENDDO
	ENDDO
	DO I = 1,6
	    DO J = 1,6
	      FIXG(I+7) = FIXG(I+7) + TRANS(I,J)*FIXD(J+7)
	    ENDDO
	ENDDO

C	----------------------------------------------------------
C	----------------------------------------------------------


	RG = 0.0D0
	DO IK = 1,7
	  RG(IK,1) = FIXG(IK)
	  RG(IK,2) = FIXG(IK+7)
	ENDDO

C	TRANSFORM WITH END RELEASE CONDITION
	CALL FIXHNG (RG,VR,RANG,ELN,IA(LHG),IA(LHGN),MLE,NELE)

	CALL FRMBEP (RG,VR,KEG,MLE,RANG)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1            A(LEI),RG,NSF,NNF,5)
711	CONTINUE

	IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RG(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	 
      
      
	ENDDO


      DEALLOCATE(RAL,MAL)	 
      
	RETURN
	END
	
C	=======================================================================
C	=======================FRAME THERMAL LOAD==============================
C	=======================================================================
	SUBROUTINE VALTEMPPAT(N,T,D,FTEMP,NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 NAME
C     INPUT      
C     N = TEMP. PATTERN GROUP NUMBER
C     T = TOTAL DIMENSION OF THE SECTION
C     D = POSITION IN SECTION TO BE CALCULATED FOR TEMP. FACTOR (Measure from bottom of section)
C     OUTPUT
C     FTEMP = OUTPUT OF TEMPERATURE FACTOR AT SPECIFIC POINT D
      
      ALLOCATABLE AA(:),TFACT(:,:),DEPT(:),LFLAG(:),TREV(:,:)
      
C     GET MAXIMUM MEMORY REQUIRE TO STORE DATA ---- NMAX           
	CALL LOCATN (NAME,KTMPT,NMAX,NTMPT,2)  
      ALLOCATE(AA(NMAX))
C     AA(1) = NUMBER OF TEMP. BLOCK OVER THE SECTION
C     AA(2) = NUMBER OF THE BLOCK THAT DEPTH IS VARY
C     AA(2+1:2+NEPT) = DEPTH OF EACH TEMP. BLOCK OVER THE SECTION
C     AA(2+NEPT+1:2+NEPT+NEPT) = FLAG FOR EACH TEMP. BLOCK OVER THE SECTION (0 = REAL VALUE, 1 = NORMALIZED VALUE)
C     AA(2+NEPT+NEPT+1:2+NEPT+NEPT+NPT) = TEMP. FACTOR OF EACH STATION POINT OVER THE SECTION  !NPT IS NUMBER OF TEMP. POINT = NEPT+1 

C     CALLING DATA 
      DO I = 1,NMAX
	CALL RELFILL(NAME,AA(I),I,N,0)
	ENDDO

      NEPT = INT(AA(1))
      NVARY= INT(AA(2))
      NPT  = NEPT+1      
      
      ALLOCATE(TFACT(2,NPT),DEPT(NEPT),LFLAG(NEPT),TREV(2,NPT))

C     DEPTH OF EACH TEMP. BLOCK OVER THE SECTION
      DO IEPT = 1,NEPT
      DEPT(IEPT) = AA(2+IEPT)
      ENDDO

C     FLAG FOR EACH TEMP. BLOCK OVER THE SECTION (0 = REAL VALUE, 1 = NORMALIZED VALUE)
      DO IEPT = 1,NEPT
      LFLAG(IEPT) = AA(2+NEPT+IEPT)
      ENDDO
            
C     TEMP. FACTOR AT EACH LEVEL      
      DO IPT = 1,NPT
      TFACT(2,IPT) = AA(2+NEPT+NEPT+IPT)
      ENDDO

C     FOR MORMALIZED BLOCK, WE NEED TO MULTIPLY IT WITH DEPTH OF SECTION
      DO IEPT = 1,NEPT
      LF = LFLAG(IEPT)
      IF(LF.EQ.1) DEPT(IEPT) = DEPT(IEPT)*T
      ENDDO      
     
C     CALCULATE FOR THE VARY DEPTH   (OVERALL DEPTH MUST EQUAL TO REAL DEPTH)  
      IF(NVARY.GT.0) THEN
          SIZE = 0.0D0 
          DO IEPT = 1,NEPT
          IF(IEPT.NE.NVARY) SIZE = SIZE + DEPT(IEPT)
          ENDDO
          DT = T - SIZE
          DEPT(NVARY) = DT
      ENDIF
      
      
C     FILL IN THE STATION ALONG THE SECTION DEPTH (FROM BOTTOM TO TOP)
      TFACT(1,1) = 0.0
      DO IEPT = 1,NEPT
      TFACT(1,IEPT+1) = TFACT(1,IEPT) + DEPT(IEPT)
      ENDDO
      
      JPT = NPT
      DO IPT = 1,NPT
      TREV(1,IPT) = T - TFACT(1,JPT)
      TREV(2,IPT) = TFACT(2,JPT)
      JPT = JPT - 1
      ENDDO

C     INTERPOLATE FOR FTEMP    
	CALL INTERPOL(D,FTEMP,TREV,NPT,IV)
      
      DEALLOCATE(AA,TFACT,DEPT,LFLAG,TREV)
      
	RETURN
	END
C	=======================================================================
C	=======================FRAME THERMAL LOAD==============================
C	=======================================================================
	SUBROUTINE PATTERN(NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 NAME
C     READ TEMPERATURE PATTERN DATA
      
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     
C     NEPT = NUMBER OF TEMP. BLOCK OVER THE SECTION
C     NVARY= NUMBER OF THE BLOCK THAT DEPTH IS VARY

      ALLOCATABLE NN(:,:),TFACT(:,:,:),AA(:)
      
      READ(ITI,*)
      READ(ITI,*) NTMPT
      
      ALLOCATE(NN(2,NTMPT),TFACT(20,3,NTMPT)) !MAX IS 20 DATA POINTS
      
      NMAX = 0
      DO I = 1,NTMPT
      READ(ITI,*) NEPT,NVARY
      NN(1,I) = NEPT ; NN(2,I) = NVARY
      NPT = NEPT+1
      READ(ITI,*) TFACT(1:NEPT,1,I) !DEPTH OF EACH TEMP. BLOCK OVER THE SECTION
      READ(ITI,*) TFACT(1:NEPT,2,I) !FLAG FOR EACH TEMP. BLOCK OVER THE SECTION (0 = REAL VALUE, 1 = NORMALIZED VALUE)
      READ(ITI,*) TFACT(1: NPT,3,I) !TEMP. FACTOR OF EACH STATION POINT OVER THE SECTION  !NPT IS NUMBER OF TEMP. POINT = NEPT+1 
      MAX = 2+NEPT+NEPT+NPT
      IF(MAX.GT.NMAX) NMAX = MAX
      ENDDO
      
      
C     ALLOCATE STORE DATA ---- NMAX,NTMPT           
	CALL DEFNREL(NAME,KTMPT,NMAX,NTMPT)
      ALLOCATE(AA(NMAX))
C     AA(1) = NUMBER OF TEMP. BLOCK OVER THE SECTION
C     AA(2) = NUMBER OF THE BLOCK THAT DEPTH IS VARY
C     AA(2+NEPT+1:2+NEPT+NEPT) = FLAG FOR EACH TEMP. BLOCK OVER THE SECTION (0 = REAL VALUE, 1 = NORMALIZED VALUE)
C     AA(2+NEPT+NEPT+1:2+NEPT+NEPT+NPT) = TEMP. FACTOR OF EACH STATION POINT OVER THE SECTION  !NPT IS NUMBER OF TEMP. POINT = NEPT+1 

C     --------------------------------
C     STORE DATA 
      DO N = 1,NTMPT
      AA(1:NMAX) = 0.0D0
      NEPT = NN(1,N) ; NVARY= NN(2,N)
      AA(1) = FLOAT(NEPT)
      AA(2) = FLOAT(NVARY)
      AA(2+1:2+NEPT) = TFACT(1:NEPT,1,N)
      AA(2+NEPT+1:2+NEPT+NEPT) = TFACT(1:NPT,2,N)
      AA(2+NEPT+NEPT+1:2+NEPT+NEPT+NPT) = TFACT(1:NPT,3,N)
            
      DO I = 1,NMAX
	CALL RELFILL(NAME,AA(I),I,N,1)
	ENDDO
	ENDDO
C     --------------------------------

      
      DEALLOCATE(NN,TFACT,AA)
      
	RETURN
	END
C	=======================================================================
C	=======================FRAME THERMAL LOAD==============================
C	=======================================================================
	SUBROUTINE TEMPLBM(LM,XYZ,IGSET,PROPG,MTSET,PROPM,NTMPLD,IGIDM) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	=======================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

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

	COMMON /LCSS/ ILCN,ILCC
C	==================================================================
	COMMON A(9000000),IA(9000000)
      
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION VR(3),VS(3),VT(3),TRANS(7,7)
	DIMENSION RG(7,3),XYZ(NCO*NNM,NELE),LM(NEF,NELE)
	DIMENSION BA(14),BB(14),FIXD(14),FIXG(14),IGIDM(NELE)
	DIMENSION BWG(10),BPG(10),BBX(14)

      ALLOCATE(RAL(NEF),MAL(NEF))

	DO ITEMP = 1,NTMPLD

	READ(ITI,*) MLE,CONT,DIFT,CONS,DIFS,TH,ALPHA,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)


	CALL XFSECTION(KEG,MLE,1)

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FTYP ,1,1 ,0) !SECTION TYPE	0=STANDARD 1=FIBER
	CALL MRELFIL(NAMEI,FGAS ,1,3 ,0) !NUMBER OF STATION POINT ALONG ELEMENT
	CALL MRELFIL(NAMEI,FFIB ,1,4 ,0) !NUMBER OF FIBER
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !ROTATION ANGLE
	ITYP = INT(FTYP)
	NGAS = INT(FGAS)
	NFIB = INT(FFIB)

	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

C	----------------------------------------------------------
	FIXD = 0.0D0
C	----------------------------------------------------------
C	CALLING ELEMENT GAUSS DATA
	CALL FRGAUS (NGAS,ELN,BPG,BWG)

C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET STRESS CONT TO ELEMENT FORCE VECTOR
C     ----------------------------------------------------------
      DO 50 II = 1,NGAS


C	CALL SECTION PROPERTIES
	CALL XFSECTION(KEG,MLE,II)
	
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
C	==========================================================
	CALL MRELFIL(NAMEI,BNORM,1,111,0) !BNORM NORMINAL WIDTH
	CALL MRELFIL(NAMEI,HNORM,1,112,0) !HNORM NORMINAL HIGHT
C	==========================================================
	CALL MRELFIL(NAMEI,EAREA,1,41,0) !Equivalent AREA-E-ALPHA
	CALL MRELFIL(NAMEI,EQS  ,1,42,0) !Equivalent  QS-E-ALPHA
	CALL MRELFIL(NAMEI,EQT  ,1,43,0) !Equivalent  QT-E-ALPHA
	CALL MRELFIL(NAMEI,ESIT ,1,44,0) !Equivalent ISS-E-ALPHA
	CALL MRELFIL(NAMEI,ETIT ,1,45,0) !Equivalent ITT-E-ALPHA
C	==========================================================

C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
	BXD = BPG(II)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
	CALL BBXFRM(BBX,ELN,BXD)

C	----------------------
	IF(TH.EQ.0.0D0) THEN

	IF(CONS.EQ.0.0.AND.DIFS.EQ.0.0) THEN  !ROTATE ABOUT T
	TH = HNORM
	ENDIF
	IF(CONT.EQ.0.0.AND.DIFT.EQ.0.0) THEN  !ROTATE ABOUT S
	TH = BNORM
	ENDIF

	ENDIF
C	----------------------
	IF(TH.EQ.0.0D0) THEN
	FIXD = 0.0D0
	GOTO 100
	ENDIF
C	----------------------

	ACON1 = CONT        !0.5*(TPT+TPB)
	BCON1 = -DIFT/TH     !(TPT-TPB)/TH
	
	ACON2 = CONS        !0.5*(TPR+TPL)
	BCON2 = -DIFS/TH     !(TPR-TPL)/TH

	BA = 0.0
	BB = 0.0
	
	IF(CONS.EQ.0.0.AND.DIFS.EQ.0.0) THEN  !ROTATE ABOUT T
	BA(1) =  BBX(1)		!EDIS(1)    !SEE ALSO FRAMENL
	BA(8) =  BBX(8)		!EDIS(8)
	BB(2) =  BBX(2)		!EDIS(2)
	BB(6) =  BBX(5)		!EDIS(6)
	BB(9) =  BBX(9)		!EDIS(9)
	BB(13)=  BBX(12)	!EDIS(13)
	DO I = 1,14
	FIXD(I) = FIXD(I) + (BA(I)*(ACON1*EAREA+BCON1*EQT) + BB(I)*(ACON1*EQT+BCON1*ETIT))*BWG(II)
	ENDDO
	ELSEIF(CONT.EQ.0.0.AND.DIFT.EQ.0.0) THEN
	BA(1) =  BBX(1)		!EDIS(1)    !SEE ALSO FRAMENL
	BA(8) =  BBX(8)		!EDIS(8)
	BB(3) =  BBX(3)		!EDIS(3)
	BB(5) =  BBX(4)		!EDIS(5)
	BB(10)=  BBX(10)	!EDIS(10)
	BB(12)=  BBX(11)	!EDIS(12)
	DO I = 1,14
	FIXD(I) = FIXD(I) + (BA(I)*(ACON2*EAREA+BCON2*EQS) + BB(I)*(ACON2*EQS+BCON2*ESIT))*BWG(II)
	ENDDO
	ENDIF

50	CONTINUE
C	----------------------------------------------------------
C	----------------------------------------------------------
C	TRANSFORMATION
100	TRANS = 0.0
	DO I = 1,3
	TRANS(I,1)   = VR(I)
	TRANS(I,2)   = VS(I)
	TRANS(I,3)   = VT(I)
	TRANS(I+3,4) = VR(I)
	TRANS(I+3,5) = VS(I)
	TRANS(I+3,6) = VT(I)
	ENDDO
	FIXG = 0.0
	DO I = 1,6
	DO J = 1,6
	FIXG(I) = FIXG(I) + TRANS(I,J)*FIXD(J)
	ENDDO
	ENDDO
	DO I = 1,6
	DO J = 1,6
	FIXG(I+7) = FIXG(I+7) + TRANS(I,J)*FIXD(J+7)
	ENDDO
	ENDDO

C	----------------------------------------------------------
C	----------------------------------------------------------


	RG = 0.0D0
	DO IK = 1,7
	RG(IK,1) = FIXG(IK)
	RG(IK,2) = FIXG(IK+7)
	ENDDO


C	TRANSFORM WITH END RELEASE CONDITION
	CALL FIXHNG (RG,VR,RANG,ELN,IA(LHG),IA(LHGN),MLE,NELE)

	CALL FRMBEP (RG,VR,KEG,MLE,RANG)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1            A(LEI),RG,NSF,NNF,5)
711	CONTINUE

      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RG(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	      
      
      
      ENDDO

      
      DEALLOCATE(RAL,MAL)	  

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE TEMPLBM_OLD(LM,XYZ,IGSET,PROPG,MTSET,PROPM,NTMPLD,IGIDM) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	=======================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

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

	COMMON /LCSS/ ILCN,ILCC
C	==================================================================
	COMMON A(9000000),IA(9000000)
      
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION VR(3),VS(3),VT(3),TRANS(7,7)
	DIMENSION RG(7,3),XYZ(NCO*NNM,NELE),LM(NEF,NELE)
	DIMENSION BA(14),BB(14),FIXD(14),FIXG(14),IGIDM(NELE)

      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO ITEMP = 1,NTMPLD

	READ(ITI,*) MLE,CONT,DIFT,CONS,DIFS,TH,ALPHA,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)


	CALL XFSECTION(KEG,MLE,1)

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !ROTATION ANGLE



	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

C	==========================================================
	CALL MRELFIL(NAMEI,EAREA  ,1,31,0) !Equivalent AREA-E-ALPHA
	CALL MRELFIL(NAMEI,ESIT   ,1,34,0) !Equivalent ISS-E-ALPHA
	CALL MRELFIL(NAMEI,ETIT   ,1,35,0) !Equivalent ITT-E-ALPHA
C	==========================================================


C	IF(TH.EQ.0) TH = SQRT(AREA)

	ACON1 = CONT        !0.5*(TPT+TPB)
	BCON1 = DIFT/TH     !(TPT-TPB)/TH
	
	ACON2 = CONS        !0.5*(TPR+TPL)
	BCON2 = DIFS/TH     !(TPR-TPL)/TH

	BA = 0.0
	BB = 0.0
	
	IF(CONS.EQ.0.0.AND.DIFS.EQ.0.0) THEN
	BA(1) = -1.0
	BA(8) =  1.0
	BB(2) =  0.0
	BB(6) = +1.0
	BB(9) =  0.0
	BB(13)= -1.0
	DO I = 1,14
	FIXD(I) =  (BA(I)*ACON1*EAREA + BB(I)*BCON1*ETIT)
	ENDDO
	ELSEIF(CONT.EQ.0.0.AND.DIFT.EQ.0.0) THEN
	BA(1) = -1.0
	BA(8) =  1.0
	BB(3) =  0.0
	BB(5) = -1.0
	BB(10)=  0.0
	BB(12)= +1.0
	DO I = 1,14
	FIXD(I) =  (BA(I)*ACON2*EAREA + BB(I)*BCON2*ESIT)
	ENDDO
	ENDIF


	TRANS = 0.0
	DO I = 1,3
	TRANS(I,1)   = VR(I)
	TRANS(I,2)   = VS(I)
	TRANS(I,3)   = VT(I)
	TRANS(I+3,4) = VR(I)
	TRANS(I+3,5) = VS(I)
	TRANS(I+3,6) = VT(I)
	ENDDO
	FIXG = 0.0
	DO I = 1,6
	DO J = 1,6
	FIXG(I) = FIXG(I) + TRANS(I,J)*FIXD(J)
	ENDDO
	ENDDO
	DO I = 1,6
	DO J = 1,6
	FIXG(I+7) = FIXG(I+7) + TRANS(I,J)*FIXD(J+7)
	ENDDO
	ENDDO

	RG = 0.0

	DO IK = 1,7
	RG(IK,1) = FIXG(IK)
	RG(IK,2) = FIXG(IK+7)
	ENDDO


C	TRANSFORM WITH END RELEASE CONDITION
	CALL FIXHNG (RG,VR,RANG,ELN,IA(LHG),IA(LHGN),MLE,NELE)

	CALL FRMBEP (RG,VR,KEG,MLE,RANG)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1            A(LEI),RG,NSF,NNF,5)
711	CONTINUE

	IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RG(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	 
      

	ENDDO


      DEALLOCATE(RAL,MAL)	
      
	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE TEMPLSH(LM,XYZ,IGSET,PROPG,MTSET,PROPM,NTMPLD,
	1				   IGIDM,NODEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)
	COMMON /LCSS/ ILCN,ILCC
C	================================================================== 
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
	DIMENSION STS(8)
      DIMENSION PROPG(NGP,1),IGSET(1),PROPM(NMP,1),MTSET(1)
	DIMENSION XYZ(NCO*NNM,1),NODEX(NEX,1),LM(NEF,1)
      DIMENSION H(9),P(2,9),XJI(4),FA(4)
      DIMENSION VR(3),VS(3),VT(3),IGIDM(1)
	DIMENSION RR(3),SS(3),STR(8),BB(8,NEF),DR(8,8)
	DIMENSION RHS(NEF)
	ALLOCATABLE STINI(:)

C	CALL FOR NUMBER OF INITIAL STRESS
	CALL MINTFIL('GWOK',NSTH,KEG,4,0)
	ALLOCATE(STINI(NSTH))
	  
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 1000 ITEMP = 1,NTMPLD

	READ(ITI,*) MLE,CONT,DIFT,ALPHA,ILCN,ILCC

	CALL ELEREODER(IGIDM,NELE,MLE)

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
      ISET = IGSET(MLE)
	TH   = PROPG(2,ISET)
	MSET = MTSET(MLE)	

	ALPHA = 0.0D0
	IF(MTMOD.NE.2)  ALPHA = PROPM(13,MSET)

C	SETTING MATERIAL PARAMETER
	CALL HOKLAW (PROPM(1,MSET),PROPG(1,ISET),2)

C	DETERMINE THE STRAIN COMPONENT DUE TO TEMPERATURE CHANGED
	STR(1) = CONT*ALPHA
	STR(2) = CONT*ALPHA
	STR(3) = 0.0D0
	STR(4) = DIFT*ALPHA/TH
	STR(5) = DIFT*ALPHA/TH
	STR(6) = 0.0D0		
	STR(7) = 0.0D0
	STR(8) = 0.0D0

	SELECT CASE(NNO)
	CASE(3)
	MMGR = 3
	MMGS = 1
	CASE(4,8)
	MMGR = 2
	MMGS = 2
	CASE(9)
	MMGR = 3
	MMGS = 3
	END SELECT


	RHS(1:NEF) = 0.0D0
	IPT = 0
C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
	IPT = IPT + 1
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	SELECT CASE(NNO)
	CASE(3)
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,MLE),P,VR,VS,VT,FA,XJI,DET,IELE,NNO)
	CASE(4,8)
	CALL SHAP2D (RI,SI,H,P,NODEX(1,MLE),NNO)
	CALL SHJACO (NNO,XYZ(1,MLE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	CASE(9)
	CALL SHAP2D9(RI,SI,H,P,NODEX(1,MLE),NNO)            
	CALL SHJACO (NNO,XYZ(1,MLE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	END SELECT
C     --------------------------------------------------------
C	CALLING STRAIN-DISPLACEMENT MATRIX
	CALL LINBATS(XJI,H,P,VR,VS,VT,BB,NNO,NNF,NEF) 
C     --------------------------------------------------------
	DVOL = DET*WT
C     -----------------------------------------
C     DETERMINE ELASTIC STRESS
C     -----------------------------------------
	CALL MPSIGA (STR,STS)
C     --------------------------------------------------------
C	INTEGRATE OVER THE SURFACE
	RHS = RHS + MATMUL(TRANSPOSE(BB),STS)*DVOL
C	--------------------------------------------------------
	NN = (NSTH/NPT)*(IPT-1)
	DO ISTH = 1,NSTH/NPT
	STINI(ISTH+NN) = STS(ISTH)
	ENDDO
	CALL ADRINI(KEG,MLE,STINI,ILCN,'ADD','VARY')  !ADD TO INITIAL STRESS VARY LOAD
	CALL ADRINI(KEG,MLE,STINI,ILCC,'ADD','CONT')  !ADD TO INITIAL STRESS CONT LOAD
C	POINTER FOR ELEMENT GROUP FOR TEMPERATURE LOAD (0 = NO TEMP LOAD ON THIS ELEMENT,1 = TEMP LOAD ON THIS ELEMENT)
	CALL MINTFIL('GWOK',1,KEG,5,1)
C	--------------------------------------------------------
800	CONTINUE        !END OF GAUSS LOOP
C	--------------------------------------------------------

	CALL FIXSHE (RHS,KEG,MLE,NNM,6) !INCLUDED MOMENT

C	--------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES
C	--------------------------------------------------------
      IF (NLS.EQ.0) GOTO 850
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1            A(LEI),RHS,NSF,NNF,4)

850	IEFL = 0
      DO I = 1,NEF
      IEQ = LM(I,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RHS(I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)
      
C	--------------------------------------------------------
1000  CONTINUE    !END TEMPERATURE LOAD
C	--------------------------------------------------------

      DEALLOCATE(RAL,MAL)	 

	DEALLOCATE(STINI)


      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE TEMPLSO(LM,XYZ,IGSET,PROPG,MTSET,PROPM,NTMPLD,
	1				   IGIDM,NODEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)
	COMMON /LCSS/ ILCN,ILCC
C	==================================================================
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
	DIMENSION STS(6)
      DIMENSION PROPG(NGP,1),IGSET(1),PROPM(NMP,1),MTSET(1)
	DIMENSION XYZ(NCO*NNM,1),NODEX(NEX,1),LM(NEF,1)
      DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)
      DIMENSION VR(3),VS(3),VT(3),IGIDM(1)
	DIMENSION STR(6),DR(6,6)
	DIMENSION RHS(NEF),B(NEF),BB(6,NEF)

	ALLOCATABLE STINI(:)

C	CALL FOR NUMBER OF INITIAL STRESS
	CALL MINTFIL('GWOK',NSTH,KEG,4,0)
	ALLOCATE(STINI(NSTH))

      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 1000 ITEMP = 1,NTMPLD

	READ(ITI,*) MLE,CONT,DIFT,ALPHA,ILCN,ILCC

	CALL ELEREODER(IGIDM,NELE,MLE)

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
      ISET = IGSET(MLE)
	MSET = MTSET(MLE)
		
	ALPHA = 0.0D0
	ALPHA = PROPM(13,MSET)

C	SETTING MATERIAL PARAMETER
	CALL HOKLAW_S (PROPM(1,MSET),PROPG(1,ISET),1)

C	DETERMINE THE STRAIN COMPONENT DUE TO TEMPERATURE CHANGED
	STR(1) = CONT*ALPHA
	STR(2) = CONT*ALPHA
	STR(3) = CONT*ALPHA
	STR(4) = 0.0D0
	STR(5) = 0.0D0
	STR(6) = 0.0D0		


	SELECT CASE(NNO)
	CASE(4)
	MGR = 4
	MGS = 1
	MGT = 1
	CASE(8)
	MGR = 2
	MGS = 2
	MGT = 2
	CASE(10)
	MGR = 4
	MGS = 1
	MGT = 1
	END SELECT


	RHS(1:NEF) = 0.0D0
	IPT = 0
C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 800  IGR=1,MGR
      RI = GLOC(IGR,MGR)
      DO 800  IGS=1,MGS
      SI = GLOC(IGS,MGS)
	DO 800  IGT=1,MGS
      TI = GLOC(IGT,MGT)
      WT = GWT(IGR,MGR)*GWT(IGS,MGS)*GWT(IGT,MGT)
	IPT = IPT + 1
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	SELECT CASE(NNO)

	CASE(4)
	CALL GAUSST (RI,SI,TI,WT,IPT,MGR,1)
	CALL SHAP3DT(RI,SI,TI,H,P,NNM)
	CALL JACO3D(XYZ(1,MLE),P,XJ,XJI,DET,MLE,NNM)
	DET = DET/6.0

	CASE(8)
	CALL SHAP3D8 (RI,SI,TI,H,P) 
	CALL JACO3D(XYZ(1,MLE),P,XJ,XJI,DET,MLE,NNM)

	CASE(10)
	CALL GAUSST (RI,SI,TI,WT,IPT,MGR,1)
	CALL SHAP3DT(RI,SI,TI,H,P,NNM)
	CALL JACO3D(XYZ(1,MLE),P,XJ,XJI,DET,MLE,NNM)
	DET = DET/6.0

	END SELECT

C     --------------------------------------------------------
C	CALLING STRAIN-DISPLACEMENT MATRIX
	CALL BMATSLT (P,XJI,B,BB,NNM)
C     --------------------------------------------------------
	DVOL = DET*WT
C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES DR
C     -----------------------------------------
	CALL DMATSLD_S(DR)
C     --------------------------------------------------------
C	INTEGRATE OVER THE VOLUME
	STS = MATMUL(DR,STR)
	RHS = RHS + MATMUL(TRANSPOSE(BB),STS)*DVOL
C	--------------------------------------------------------
	NN = (NSTH/NPT)*(IPT-1)
	DO ISTH = 1,NSTH/NPT
	STINI(ISTH+NN) = STS(ISTH)
	ENDDO
	CALL ADRINI(KEG,MLE,STINI,ILCN,'ADD','VARY')  !ADD TO INITIAL STRESS VARY LOAD
	CALL ADRINI(KEG,MLE,STINI,ILCC,'ADD','CONT')  !ADD TO INITIAL STRESS CONT LOAD
C	POINTER FOR ELEMENT GROUP FOR TEMPERATURE LOAD (0 = NO TEMP LOAD ON THIS ELEMENT,1 = TEMP LOAD ON THIS ELEMENT)
	CALL MINTFIL('GWOK',1,KEG,5,1)
C	--------------------------------------------------------
800	CONTINUE        !END OF GAUSS LOOP
C	--------------------------------------------------------

	CALL FIXSOL (RHS,KEG,MLE,NEF,NNM)

C	--------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES
C	--------------------------------------------------------
      IF (NLS.EQ.0) GOTO 850
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1            A(LEI),RHS,NSF,NNF,4)

850	IEFL = 0
      DO I = 1,NEF
      IEQ = LM(I,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RHS(I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)
      
C	--------------------------------------------------------
1000  CONTINUE    !END TEMPERATURE LOAD
C	--------------------------------------------------------

      DEALLOCATE(RAL,MAL)	

	DEALLOCATE(STINI)

      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE TEMPLTS(LM,XYZ,IGSET,PROPG,MTSET,PROPM,NTMPLD,
	1				   IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)
	COMMON /LCSS/ ILCN,ILCC
C	=================================================================
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
      DIMENSION PROPG(NGP,1),IGSET(1),PROPM(NMP,1),MTSET(1)
	DIMENSION XYZ(NCO*NNM,1),LM(NEF,1)
      DIMENSION H(NNM),P(NNM)
      DIMENSION VR(3),IGIDM(1)
	DIMENSION RHS(NEF),BB(NEF)

	ALLOCATABLE STINI(:)

C	CALL FOR NUMBER OF INITIAL STRESS
	CALL MINTFIL('GWOK',NSTH,KEG,4,0)
	ALLOCATE(STINI(NSTH))

      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 1000 ITEMP = 1,NTMPLD

	READ(ITI,*) MLE,CONT,DIFT,ALPHA,ILCN,ILCC

	CALL ELEREODER(IGIDM,NELE,MLE)

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
      ISET = IGSET(MLE)
      AREA = PROPG(2,ISET)
	MSET = MTSET(MLE)
		
	ALPHA = 0.0D0
	ALPHA = PROPM(13,MSET)

C	SETTING MATERIAL PARAMETER
	YG = PROPM(1,MSET)

C	DETERMINE THE STRAIN COMPONENT DUE TO TEMPERATURE CHANGED
	STR = CONT*ALPHA

	IF(ISTYP.EQ.4) STR = -1.0*ABS(STR)         !ONLY CONTRACTION FOR CABLE ELEMENT

	MGR = 1
	RHS(1:NEF) = 0.0D0
C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
	IPT = 0
      DO 800  IGR=1,MGR
      RI = GLOC(IGR,MGR)
      WT = GWT(IGR,MGR)
	IPT = IPT +1
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	CALL SHAP1D(RI,H,P,NNM)
	CALL JACO1D(XYZ(1,MLE),P,VR,DET,MLE,NNM)
C     --------------------------------------------------------
C	CALLING STRAIN-DISPLACEMENT MATRIX
	CALL LINBATT(P,DET,VR,BB,NNM,NNF,NEF) 
C     --------------------------------------------------------
	DVOL = DET*WT
C     --------------------------------------------------------
C	INTEGRATE OVER THE SURFACE
	STS = YG*STR*AREA
	RHS = RHS + BB*STS*DVOL
C	--------------------------------------------------------
	NN = (NSTH/NPT)*(IPT-1)
	DO ISTH = 1,NSTH/NPT
	STINI(ISTH+NN) = STS
	ENDDO
	CALL ADRINI(KEG,MLE,STINI,ILCN,'ADD','VARY')  !ADD TO INITIAL STRESS VARY LOAD
	CALL ADRINI(KEG,MLE,STINI,ILCC,'ADD','CONT')  !ADD TO INITIAL STRESS CONT LOAD
C	POINTER FOR ELEMENT GROUP FOR TEMPERATURE LOAD (0 = NO TEMP LOAD ON THIS ELEMENT,1 = TEMP LOAD ON THIS ELEMENT)
	CALL MINTFIL('GWOK',1,KEG,5,1)
C	--------------------------------------------------------
800	CONTINUE        !END OF GAUSS LOOP
C	--------------------------------------------------------

	CALL FIXTRS (RHS,KEG,MLE,NEF,NNM)
		
C	--------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES
C	--------------------------------------------------------
      IF (NLS.EQ.0) GOTO 850
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1            A(LEI),RHS,NSF,NNF,4)

850	IEFL = 0
      DO I = 1,NEF
      IEQ = LM(I,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RHS(I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)
      
C	--------------------------------------------------------
1000  CONTINUE    !END TEMPERATURE LOAD
C	--------------------------------------------------------

      DEALLOCATE(RAL,MAL)	  

	DEALLOCATE(STINI)


      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE TEMPLME(LM,XYZ,IGSET,PROPG,MTSET,PROPM,NTMPLD,
	1				   IGIDM,NODEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)
	COMMON /LCSS/ ILCN,ILCC
C	==================================================================
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
	DIMENSION STS(4)
      DIMENSION PROPG(NGP,1),IGSET(1),PROPM(NMP,1),MTSET(1)
	DIMENSION XYZ(NCO*NNM,1),NODEX(NEX,1),LM(NEF,1)
      DIMENSION H(9),P(2,9),XJ(4),XJI(4)
      DIMENSION IGIDM(1),STR(4),BB(4,NEF)
	DIMENSION TSM(6,4),DR(6,6),DM(4,4)
	DIMENSION RHS(NEF)
	ALLOCATABLE STINI(:)

C	CALL FOR NUMBER OF INITIAL STRESS
	CALL MINTFIL('GWOK',NSTH,KEG,4,0)
	ALLOCATE(STINI(NSTH))

	TSM = 0.0D0
	TSM(1,1) = 1.0D0
	TSM(2,2) = 1.0D0
	TSM(3,3) = 1.0D0
	TSM(6,4) = 1.0D0

      ALLOCATE(RAL(NEF),MAL(NEF))

	DO 1000 ITEMP = 1,NTMPLD

	READ(ITI,*) MLE,CONT,DIFT,ALPHA,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	

	CALL ELEREODER(IGIDM,NELE,MLE)

      ISET = IGSET(MLE)
	TH   = PROPG(2,ISET)
	MSET = MTSET(MLE)

		
	ALPHA = 0.0D0
	ALPHA = PROPM(13,MSET)

C	DETERMINE THE STRAIN COMPONENT DUE TO TEMPERATURE CHANGED
	STR(1) = CONT*ALPHA
	STR(2) = CONT*ALPHA
	STR(3) = 0.0D0
	STR(4) = 0.0D0


	MMGR = SQRT(1.0*NPT)
	MMGS = SQRT(1.0*NPT)

	ISTYM = ISTYP
	IF(ISTYP.GE.3) ISTYM  = ISTYP - 3
	IF(ISTYM.EQ.0) STR(4) = CONT*ALPHA !AXISYMETRIC
	IF(ITYPE.EQ.8) ISTYM  = 1          !PLANE STRAIN CONSOLIDATION ELEMENT

C	SETTING MATERIAL PARAMETER
	CALL HOKLAW (PROPM(1,MSET),PROPG(1,ISET),ISTYM)

	RHS(1:NEF) = 0.0D0
	IPT = 0
C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
	IPT = IPT + 1
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	CALL SHAP2D (RI,SI,H,P,NODEX(1,MLE),NNO)
	CALL JACOB2D(MLE,NNM,XYZ(1,MLE),P,XJ,XJI,DET,NCO)
C     --------------------------------------------------------
C	CALLING STRAIN-DISPLACEMENT MATRIX
	CALL MEBMAT (XYZ(1,MLE),H,P,XJI,BB,ISTYM,XBAR,NNM)
C     --------------------------------------------------------
	IF (ISTYM.NE.0) XBAR = TH
      DVOL = WT*DET*XBAR
C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES DR
C     -----------------------------------------
	CALL DELA3D (DR,ISTYM)
	DM = MATMUL(TRANSPOSE(TSM),MATMUL(DR,TSM))
C     --------------------------------------------------------
C	INTEGRATE OVER THE SURFACE
	STS = MATMUL(DM,STR)
	RHS = RHS + MATMUL(TRANSPOSE(BB),STS)*DVOL
C	--------------------------------------------------------
	NN = (NSTH/NPT)*(IPT-1)
	DO ISTH = 1,NSTH/NPT
	STINI(ISTH+NN) = STS(ISTH)
	ENDDO
	CALL ADRINI(KEG,MLE,STINI,ILCN,'ADD','VARY')  !ADD TO INITIAL STRESS VARY LOAD
	CALL ADRINI(KEG,MLE,STINI,ILCC,'ADD','CONT')  !ADD TO INITIAL STRESS CONT LOAD
C	POINTER FOR ELEMENT GROUP FOR TEMPERATURE LOAD (0 = NO TEMP LOAD ON THIS ELEMENT,1 = TEMP LOAD ON THIS ELEMENT)
	CALL MINTFIL('GWOK',1,KEG,5,1)
C	--------------------------------------------------------
800	CONTINUE        !END OF GAUSS LOOP
C	--------------------------------------------------------

	CALL FIXMRE (RHS,KEG,MLE,NEF,NNM)

C	--------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES
C	--------------------------------------------------------
      IF (NLS.EQ.0) GOTO 850
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1            A(LEI),RHS,NSF,NNF,4)

850   IEFL = 0
      DO I = 1,NEF
      IEQ = LM(I,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RHS(I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)
      
C	--------------------------------------------------------
1000  CONTINUE    !END TEMPERATURE LOAD
C	--------------------------------------------------------

      DEALLOCATE(RAL,MAL)	  

	DEALLOCATE(STINI)


      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE LINBATS(XJI,H,P,VR,VS,VT,BM,NNM,NNF,NEF) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION XJI(2,2),H(NNM),P(2,NNM),VR(3),VS(3),VT(3)
	DIMENSION F(NNM),G(NNM)
	DIMENSION BM(8,1)


	DO I = 1,NNM
	F(I) = XJI(1,1)*P(1,I) + XJI(1,2)*P(2,I)
	G(I) = XJI(2,1)*P(1,I) + XJI(2,2)*P(2,I)
	ENDDO
	
	BM(1:8,1:NEF) = 0.0D0

	DO I = 1,NNM
	NN = NNF*(I-1)
	DO J = 1,3
	JJ = J  + NN
	KK = JJ + 3
	BM(1,JJ) = F(I)*VR(J)                !MEMBRANE
	BM(2,JJ) = G(I)*VS(J)
	BM(3,JJ) = F(I)*VS(J) + G(I)*VR(J)
	BM(4,KK) = F(I)*VS(J)                !BENDING
	BM(5,KK) =-G(I)*VR(J)
	BM(6,KK) = G(I)*VS(J) - F(I)*VR(J)
	BM(7,JJ) = F(I)*VT(J)                !SHEAR
	BM(7,KK) = H(I)*VR(J)
	BM(8,JJ) = G(I)*VT(J)                
	BM(8,KK) =-H(I)*VS(J)
	ENDDO
	
	ENDDO	


	RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE JACO1D(XY,P,VR,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(3,NNO),P(NNO),XJ(3),VR(3)

C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  J=1,3
      DUM = 0.0
      DO 90   K=1,NNO
 90   DUM = DUM + P(K)*XY(J,K)
 100  XJ(J) = DUM
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = SQRT(XJ(1)*XJ(1) + XJ(2)*XJ(2) + XJ(3)*XJ(3)) 
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')   

	VR(1:3) = XJ(1:3)/DET
C
      RETURN
      END


C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE LINBATT(P,DET,VR,BM,NNM,NNF,NEF) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION P(NNM),VR(3)
	DIMENSION BM(1)


	BM(1:NEF) = 0.0D0
	
	DO I = 1,NNM

	NN = NNF*(I-1)
	DO J = 1,3
	JJ = J  + NN
	BM(JJ) = P(I)*VR(J)/DET
	ENDDO
	
	ENDDO	


	RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE INISUBT(ITYPE,WOREL,IEG,IEL,NPT,NWG,ILC,RHO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	==================================================================
	DIMENSION WOREL(1)
	ALLOCATABLE STINI(:),STINO(:)


	IF(ILC.LE.0) RETURN

C	CALL POINTER FOR ELEMENT GROUP FOR TEMPERATURE LOAD (0 = NO TEMP LOAD ON THIS ELEMENT,1 = TEMP LOAD ON THIS ELEMENT)
	CALL MINTFIL('GWOK',NSTHG,IEG,5,0)
	IF(NSTHG.EQ.0) RETURN


	CALL MINTFIL('GWOK',NSTH,IEG,4,0)

	MSTH = NSTH/NPT

	ALLOCATE(STINI(NSTH),STINO(NSTH))
	CALL ADRINI(IEG,IEL,STINI,ILC,'RED','VARY')	!INITIAL STRESS VARY LOAD
	CALL ADRINI(IEG,IEL,STINO,ILC,'RED','CONT')	!INITIAL STRESS VARY LOAD


	SELECT CASE(ITYPE)	

	CASE(2)
	WOREL(1) = WOREL(1) - STINI(1)*RHO - STINO(1)

	CASE(6,8,9,10,11)
	DO IPT = 1,NPT
	NUM = NWG *(IPT-1)
	MUM = MSTH*(IPT-1)
	DO ISTH = 1,MSTH
	WOREL(ISTH + NUM) = WOREL(ISTH + NUM) - STINI(ISTH + MUM)*RHO
	1									  - STINO(ISTH + MUM)
	ENDDO
	ENDDO


	ENDSELECT



	DEALLOCATE(STINI,STINO)


	RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
