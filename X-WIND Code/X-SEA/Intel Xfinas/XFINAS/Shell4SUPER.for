      SUBROUTINE SUPELLOOP (PROPM,PROPG,PROPO,IPINS,MTSET,IGSET,IOSET,IHSET,
     1           NODEX,LM,RE,DISLI,DISP,S,COORD,EDIS,EDISI,ELOD,
     2           AMV,IAX,MMP,MGP,MEX,MWA,LMRCT,REAC,
     3		   MFQ,SETIF,IDRCT,CABFZ,CABFX,CABFF,TAMBT,TCHET,ISFAC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     SETS UP ELEMENT LOOP
C	----------------------------------------------------------------
C     CALLS ELCORD TO SET ELEMENT DISPLACEMENTS AND CO-ORDINATES
C     CALLS MODULE TO CALCULATE EQ.LOADS AND STIFFNESS MATRIX(IFREF=0)
C     CALLS ADDBAN TO ASSEMBLE NEW TANGENTIAL STIFFNESS (IFREF=0)
C     CALLS ELPRIN TO PRINT STRESSES (IFPRI=0)
C	----------------------------------------------------------------
C     FOR INPUT AND OUTPUT VARIABLES SEE ROUTINE LOCATI
C     ----------------------------------------------------------------
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /DYNA/ CDEN,IMASS
      COMMON /FTIM/ TIM(20),IDATE,ITIME

      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(7000000),IA(6000000)
C
	COMMON /GiDEle/ LGID 

C	COMMON BLOCK FOR EAS SONGSAK FEB2006
	COMMON /MMENH/ MM,MM1,MM2,NDIMC


C	LOADCASE
	COMMON /STCAS/ ILC

C	==================================================================
C	CABLE PRETENSION OPTIMIZATION
	COMMON /CB556/ LCPZ,NCABZ,KBOPZ,NCOBZ
C	==================================================================
C	-------------------------------------------------------------
	COMMON /CB557I/ MCSUM,MCFOC,MCDEF,MCMOM,LS557,LF557,LD557,LM557,
	1			    LR557,NM556
	COMMON /CB557R/ TPOPZ(10000)
C	-------------------------------------------------------------
C	CONSOLIDATION SACHARUCK DEC2006
	COMMON /CONSO/ NCONSO

C	COMMON BLOCK FOR HEAT SONGSAK MAR2007
	COMMON /SHEAT/ NHAEF,LHEAT1,LHEAT2,LHEAT3,LHEAT4
C	==================================================================
C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN
C	==================================================================
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
      COMMON /TALTA/ TAS
C	==================================================================

C	MOVING LOAD PRINTING FLAG SONGSAK JUL2007
	COMMON /LANPRIN/  LAN_PRN,ILPL,MLANE,NPL,LAN_OPT
     
C     SAVE IEL
      COMMON /INELE/ IEL

C     FOR POINTLOAD CASE CONTROL DIAPHRAGM CONSTRATINS BY BJ
      COMMON /JLCN/ ILCNM(1000),IL2 
      COMMON /IJEG/ IEG,JEG,JSTYP

C     FOR SUPER ELEMENT BY BJ
      COMMON /SUPERSTIF/ SM(24,24)       
      COMMON /STOREN/ IELEMENT(100000),N1(100000),N2(100000),N3(100000),N4(100000),N5(100000),N6(100000)
     1               ,N7(100000),N8(100000)  
CC      COMMON /CONNE/ INCINEW(4,100000) !TO SAVE THE CONNECTIVITY INTO INCINEW MATRIX BY BJ
      COMMON /JNELE/ JEL
      COMMON /IBO/ IBOUND(5000,10) ! BOUNDARY CONDITIONS WERE SAVED  
      COMMON /HEI/ HEIGHT ! HEIGHT OF EACH SUPER ELEMENT
      COMMON /LINKSTIFF/ S2(6,6)
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV 
      COMMON /KITEC/ KITE2      
      
      COMMON /DISO/ DISPO(1000)
      
      COMMON /STYPSU/ISTYPSUPER
      
      ALLOCATABLE ISUPPT(:,:)
      
      ALLOCATABLE LNNM(:,:) !ELEMENT CONNECTIVITY
C	==================================================================        

      DIMENSION PROPM(MMP,1),PROPG(MGP,1),MTSET(1),IGSET(1),NODEX(MEX,1)
	DIMENSION PROPO(6,1),IOSET(1),PROPOF(6)
	DIMENSION IPINS(14,1),IHSET(1),IPINSF(14)
      DIMENSION LM(NEF,1),WA(MWA,1),RE(1),DISLI(1),DISP(1),S(1)
      DIMENSION COORD(1),EDIS(1),EDISI(1),ELOD(1)
	DIMENSION AMV(3,NMV),IAX(NELE)
	DIMENSION FIN(NEF)


C	EAS VARIABLES SONGSAK MAR2006
	DIMENSION ASEL(MM,MM1),APHA(MM),
	1		  AINF(MM),ADIS(MM1+MM2),
	2		  ASEQ(MM,MM2),AINP(MM),
	3		  APHI(MM)
	DIMENSION SEDI(MM,MM),SEL(MM,MM1),SEQ(MM,MM2)
	DIMENSION ALPHA(MM),ALPHI(MM)

	DIMENSION HINFC(MM),HINFP(MM)
C	============================

      DIMENSION LMRCT(NEF,NELE),REAC(MFQ),SETIF(NEQ,LCS),
	1		  IDRCT(NSF,NSN)
	DIMENSION FIXEO(NEF),FIXLO(NEF),FIXEN(NEF),FIXLR(NEF)	

C	FOR CABLE PRETENSION OPTIMIZATION SONGSAK NOV2006     
	DIMENSION CABFZ(NEQ,NCABZ),CABFX(NEQ,NCABZ),CABFF(NCABZ),CABVEC(6)

C	WORKING ARRAY IN THE NEW STORAGE BLOCK SONGSAK JAN2007
	DIMENSION WOREL(MWA)
	DIMENSION TCHET(NHAEF*NELE,LCS),TAMBT(2,LCS)  !HEAT SONGSAK MAR2007
	DIMENSION PMATRL(NMP,NMPS)                    !STORE WHOLE MATERIAL PROPERTIES SONGSAK JUN2007

C     SHELL WITH NO ROTATION
C     NEFC = NUMBER OF ELEMENT FACE = NEFC  (SEE ELEMIN)
      DIMENSION ISFAC(2*NEFC,1)   !STORE SUPPORT FLAG FOR EACH ELEMENT FACE (SEE SHELL 3 NODE ONATE)
      
C     FOR SUPER ELEMENT BY BJ
      DIMENSION XYZM(25,3),ICONN(4,16),COORDMESH(12),IDSUPER(25,6),IEDOF(24),SKS(150,150),SKSSUP(150,150),
     1          SKSC(24,24),SKSCSUP(24,24), !CONDENSED STIFFNESS MATRIX 
     2          ELODM(24,24),ELODAM1(126,1),ELODAM(150,150),ELODAMC(24,24), !ASSEMBLED INTERNAL FORCES MATRIX 
     3          SAS(300),TEDIS(150,1),
     4          TMACON(150,24),TMACONI(24,150),TMACONM(126,24),TMACONIM(150,150),
     5          EDIS2(NEF),ELOD2(NEF),ELODCON(NEF),DELOD(NEF),EDISC(NEF),EDISIC(NEF),
     6          WOREL2(MWA),WORELM(MWA)!,
C     7          SUPPT(6,7) 
      DIMENSION NSAVE(5000,2)
      DIMENSION INTERMAT(NELE,8), ! SAVE THE NODES IN INTERFACE SIDE OF ELEMENTS (8 IS TOTAL NUMBER OF SIDES IN ONE ELEMENT)
     1          INTERELE(NELE,8),INTERELE2(NELE,8),INTERELE3(NELE,8) ! SAVE INTERFACE SIDE NUMBER
      DIMENSION COORDMAT(NSN,3) ! SAVE THE COORDINATE INTO COORDMAT MATRIX
      DIMENSION STIFFA(30,30),STIFFBEAM(12,12),NODEMAT(5),IEDOFB(30),IEDOFC(12),STIFFAMAT(150,150),STIFFAMATC(24,24),
     1          SKSCONE(300) !SAVE SKSC MATRIX INTO ONE LINE      
      DIMENSION AA1(3,1),AA2(1,3)
      DIMENSION REACNEW(24),EDISNEW(24),SKSAVE(150,150),SKCSAVE(24,24),SKSCMEM(24,24),DDELOD(150,1),DISPSUPER(24)
       
 
      ISTYPSUPER = ISTYP
C	-----------------------
C	SONGSAK ADD THID CONDITION SINCE THERE IS NO ELEMENT DAMPING INSIDE
C	-----------------------
	IF(ITASK.EQ.6) RETURN

C      IF(ITASK.EQ.3) THEN
C          IFREF = 0 !FLAG TO GET STIFFNESS MATRIX FOR ITASK = 3
C      ENDIF      
C	------------------------------------------------------------------
C     ELEMENT LOOP,INITIALISATION
C     ---------------------------
      DO I = 1,NMPS
	DO J = 1,NMP
	PMATRL(J,I) = PROPM(J,I)
	ENDDO
	ENDDO

	CALL CPU_TIME (TIM3)
      
	NEQF = 1
      NEQL = NEQ
      NPRE = 0

	IF(ITYPE.EQ.2.AND.ISTYP.EQ.4.AND.KBOPZ.EQ.1) THEN  !CABLE OPTIMIZATION
	CABFX = 0.0D0
	CABFZ = 0.0D0
      ENDIF

      
C     TO FIND THE INTERFACE ELEMENT IN ONE SUPER ELEMENT ----------------------
      CALL CLEARMATI(INTERMAT,NELE,8) !CLEARMATI : CLEAR INTEGER MATRIX
      CALL CLEARMATI(INTERELE,NELE,8) !CLEARMATI : CLEAR INTEGER MATRIX
      CALL CLEARMATI(INTERELE2,NELE,8) !CLEARMATI : CLEAR INTEGER MATRIX
      CALL CLEARMATI(INTERELE3,NELE,8) !CLEARMATI : CLEAR INTEGER MATRIX
      
      CALL INTFILL('OGRF',NFL4 ,14,KEG,0) !CALL FILE STORING THE CONNECTIVITY
	CALL INTFILL('OGRF',NNNM ,4 ,KEG,0) !LENGTH OF ELEMENT CONNECTIVITY (EQUAL TO NNM)
      ALLOCATE(LNNM(NNNM,2))
      
      ! TO GET INTERFACE SIDE IN ONE SUPER ELEMENT
      DO IEL1 = 1,NELE  
          
	    IGM1 = IA(LGID+IEL1-1)
          CALL READCON(IGM1,LNNM(1,1),NFL4)
      
          INT = 0
          INTR = 0
          DO INC1 = 1,4
             IF(INC1.EQ.4) THEN
                 NO1 = LNNM(4,1) !INCINEW(4,IEL1)
                 NO2 = LNNM(1,1) !INCINEW(1,IEL1)
             ELSE    
                 NO1 =  LNNM(INC1+0,1) !INCINEW(INC1,IEL1)
                 NO2 =  LNNM(INC1+1,1) !INCINEW(INC1+1,IEL1)
             ENDIF
             
             DO IEL2 = 1,NELE
                 IF(IEL1.EQ.IEL2) GO TO 1

                 IGM2 = IA(LGID+IEL2-1)
                 CALL READCON(IGM2,LNNM(1,2),NFL4)
                 
                 DO INC2 = 1,4
                     IF(INC2.EQ.4) THEN
                         NO3 = LNNM(4,2) !INCINEW(4,IEL2)
                         NO4 = LNNM(1,2) !INCINEW(1,IEL2)
                     ELSE    
                         NO3 =  LNNM(INC2+0,2) !INCINEW(INC2,IEL2)
                         NO4 =  LNNM(INC2+1,2) !INCINEW(INC2+1,IEL2)
                     ENDIF
                     
                     IF(NO1.EQ.NO3.OR.NO1.EQ.NO4) THEN
                         IF(NO2.EQ.NO3.OR.NO2.EQ.NO4) THEN
                             INT = INT + 1
                             INTERMAT(IEL1,2*(INT-1)+1) = NO1
                             INTERMAT(IEL1,2*(INT-1)+2) = NO2
                             INTERELE(IEL1,2*(INT-1)+1) = INC1
                             IF(INC1.EQ.1) THEN
                                 INC11 = 3
                             ELSEIF(INC1.EQ.2) THEN
                                 INC11 = 4 
                             ELSEIF(INC1.EQ.3) THEN
                                 INC11 = 1
                             ELSEIF(INC1.EQ.4) THEN
                                 INC11 = 2    
                             ENDIF    
                             INTERELE(IEL1,2*(INT-1)+2) = INC11
                         ENDIF    
                     ENDIF    
                 ENDDO
1            ENDDO
          ENDDO
      ENDDO        

      
      ! CASE OF TOTAL ONE SHELL ELEMENT
      ! : NEED TO CLASSIFY THE SIDE THAT PUT THE FICTITOUS BEAM
      IF(INT.EQ.0) THEN
          !GO TO 2
          
          IEL1 = 1
	    IGM1 = IA(LGID+IEL1-1)
          CALL READCON(IGM1,LNNM(1,1),NFL4)
          
          DO II = 1,4
              INT2 = 0
              IF(II.EQ.4) THEN
                  INCI1 = LNNM(II,1) !INCINEW(II,1)
                  INCI2 = LNNM(1 ,1) !INCINEW(1,1)
              ELSE    
                  INCI1 = LNNM(II+0,1) !INCINEW(II,1)
                  INCI2 = LNNM(II+1,1) !INCINEW(II+1,1)
              ENDIF
              DO JJ = 1,NBS-1
                   IBNUM1 = IBOUND(JJ,1)
                   IBNUM2 = IBOUND(JJ+1,1)
                   
                   IF(INCI1.EQ.IBNUM1.OR.INCI1.EQ.IBNUM2) THEN
                       IF(INCI2.EQ.IBNUM1.OR.INCI2.EQ.IBNUM2) THEN
                           INT2 = INT2 + 1
                           INTERELE(1,2*(INT2-1)+1) = II
                             IF(II.EQ.1) THEN
                                 II2 = 3
                             ELSEIF(II.EQ.2) THEN
                                 II2 = 4 
                             ELSEIF(II.EQ.3) THEN
                                 II2 = 1
                             ELSEIF(II.EQ.4) THEN
                                 II2 = 2    
                             ENDIF    
                             INTERELE(1,2*(INT2-1)+2) = II2                           
                       ENDIF
                   ENDIF    
              ENDDO
          ENDDO
      ENDIF          

      DEALLOCATE(LNNM)
      
      
!2     !IF(INT.EQ.0) THEN
         ! DO JJ =1,NELE
              !INTERELE3(JJ,1) = 2
              !INTERELE3(JJ,2) = 4
          !ENDDO  
      !ELSE    
          DO JJ =1,NELE
              INTERELE3(JJ,1:8) = INTERELE(JJ,1:8)
          ENDDO
      !ENDIF
      
      DO II = 1,NELE
      INTERELE(II,1) = 1
      INTERELE(II,2) = 2
      INTERELE(II,3) = 3
      INTERELE(II,4) = 4    
      
      !INTERELE(II,5) = 1
      !INTERELE(II,6) = 2
      !INTERELE(II,7) = 3
      !INTERELE(II,8) = 4       
      ENDDO      
      
C     ---------------------------------------------------------------------------                 
      
      TAS = 0.0	
C     FOR READING T MATRIX BY BJ      
      KREC = 5002
      REWIND(KREC)
      KREC1 = 5003
      REWIND(KREC1)
      KREC2 = 5004
      REWIND(KREC2)
       
      
	DO 900  IEL = 1,NELE

      CALL CLEARA (S   ,NWS)
      CALL CLEARA (ELOD,NEF)
      CALL CLEARA (FIN ,NEF)

      MSET   = MTSET(IEL)
      ISET   = IGSET(IEL)
      
      PROPOF(1:6 ) = 0.0D0
      IPINSF(1:14) = 0
      IF(ITYPE.EQ.5) THEN
	    IOET   = IOSET(IEL)
	    IHET   = IHSET(IEL)
	    IF(IOET.NE.0) PROPOF(1:6 ) = PROPO(1:6 ,IOET)
	    IF(IHET.NE.0) IPINSF(1:14) = IPINS(1:14,IHET)
      ENDIF

      MEL    = IEL

C	FIXEND FORCE IN GLOBAL SYSTEM FOR REACTION (FRAME ELEMENT)
C	SONGSAK JUN2006
      CALL CLEARA (FIXLR,NEF)

C     --------------------------------------------------------
C     TOTAL DISPLACEMENTS (EDIS) AND INCREMENTAL DISP.(EDISI)
C     NODAL CO-ORDINATES FOR TOTAL OR UPDATED LAGRANGIAN FORM.
C     --------------------------------------------------------
CB      IF(ITASK.EQ.3) DISP(1:NEQ) = DISPO(1:NEQ)
      CALL ELCORD (A(LCO),LM,A(LDI),DISP,COORD,EDIS,EDISI,NEF)
 
	 
C	ADDED THE PRESCRIBE VALUE OF NODAL DISPLACEMENT TO THE ELEMENT DOF
C	SONGSAK
	IF(ITASK.NE.1) THEN
	CALL SETELM (EDIS,LMRCT(1,IEL),IDRCT,NEF,LM(1,IEL),IA(LID))
	ENDIF

C     -----------------------------------------
C     COMPUTATION OF EQUILIBRIUM LOADS ELOD(NEF)
C     AND STIFFNESS MATRIX S(NWS) (IFREF=0)
C     -----------------------------------------
	NV = IAX(IEL)  !FOR BRIDGE LOAD ANALYSIS

C ================================ START FOR SUPER ELEMENT ==================================================      

C     DIVIDE MESH IN EACH ONE ELEMENT- GET COORDINATE AND CONNECTIVI
      CALL MESHONE(COORD,XYZM,ICONN,IDSUPER)

      !INITIALIAZE ASSEMBLED MATRIX OF ONE SUPER ELEMENT 
      CALL CLEARMAT(SKS,150,150)
      !INITIALIAZE ASSEMBLED MATRIX OF THE INTERNAL FORCES
      CALL CLEARMAT(ELODAM,150,150)      

C	===================================================================
C	CALLING WORKING ARRAY FROM NEW STORAGE 
	CALL ADREWT(KEG,IEL,WOREL,'RED')
      ! SAVE WOREL INTO WOREL2 MATRIX
        WOREL2(1:MWA) = WOREL(1:MWA)
C	===================================================================
	CALL CLEARA(FIXEN,NEF)
	CALL CLEARA(FIXEO,NEF)


	SELECTCASE(LAN_PRN)
	CASE(0)
	CALL FIXSTOR(ILC,KEG,IEL,FIXEN,NEF,1.0D0,'RED','VARY')	!FORCE IN GLOBAL SYSTEM NOT LOCAL SUPPORT (IN ELEMENT LOCAL FOR FRAME ELEMENT) VARY
	CALL FIXSTOR(ILC,KEG,IEL,FIXEO,NEF,1.0D0,'RED','CONT')	!FORCE IN GLOBAL SYSTEM NOT LOCAL SUPPORT (IN ELEMENT LOCAL FOR FRAME ELEMENT) CONSTANT
	CASE(1)
	CALL LANFIXF3(FIXEN,IEL,KEG)							!MOVING LOAD PRINTING FLAG SONGSAK JUL2007
	ENDSELECT

	IF(ITYPE.EQ.17) THEN !TENDON CALL DATA
	CALL TDINIT(IA(LTDN),WOREL,IEL,ILC,KSTEP,MSET,ISET,'CALL') 
      ENDIF

C     ------------------------------------------------------------------------      
C     DIVIDE MESH IN EACH ONE ELEMENT- GET COORDINATE AND CONNECTIVI
      CALL MESHONE(COORD,XYZM,ICONN,IDSUPER)

      !INITIALIAZE ASSEMBLED MATRIX OF ONE SUPER ELEMENT
      CALL CLEARMAT(SKSC,24,24)
      CALL CLEARMAT(SKSCSUP,24,24)
      CALL CLEARMAT(SKSAVE,150,150)
      CALL CLEARMAT(SKCSAVE,24,24)
      !INITIALIAZE ASSEMBLED MATRIX OF THE INTERNAL FORCES
      CALL CLEARMAT(ELODAM,150,150)   
      CALL CLEARMAT(STIFFAMAT,150,150)

      CALL  MATRIX_CONDEN_TEST
      
C     SAVE EDIS IN  EDISNEW MATRIX
C      IF(ITASK.EQ.3) THEN
      IF(ITASK.EQ.2.OR.ITASK.EQ.3) THEN
         EDISNEW(1:24) = EDIS(1:24)
      ENDIF

C     SAVE THE T MATRIX INTO THE FILE (KREC=5002)
C      IF(ITASK.EQ.3) THEN
       IF(ITASK.EQ.2.OR.ITASK.EQ.3.OR.ITASK.EQ.5) THEN
C      IF(ITASK.EQ.2.OR.ITASK.EQ.3) THEN
          READ(KREC,1000) TMACON(1:150,1:24)
          IF(ITASK.EQ.5) THEN
              GO TO 4
          ELSE
              CALL CALTEDIS(TMACON,EDIS,TEDIS)
          ENDIF
4      ENDIF
 
C     -------------------------------------------------------------------------------------

      DO JEL = 1,16 !===================================== START OF LOOP FOR 4X4 MESHES IN ONE SUPER ELEMENT   
      MEL = JEL
      ! RETURN WOREL FROM WOREL2 MATRIX
        WOREL(1:MWA) = WOREL2(1:MWA)
           
      !INITIALIZATION OF MATRICES    
      CALL CLEARA (SAS   ,NWS)
      CALL CLEARA (ELOD,NEF)
      CALL CLEARA (ELOD2,NEF)
      CALL CLEARMAT(ELODM,NEF,NEF)
      CALL CLEARMAT(SM,24,24)
      CALL CLEARA(WORELM,MWA)

          
      ! FOR EDIS FROM TEDIS         
C      IF(ITASK.EQ.3) THEN   
      IF(ITASK.EQ.2.OR.ITASK.EQ.3) THEN
          DO KEL = 1,4
              NT = ICONN(KEL,JEL) !READ NODE BY CONNECTIVITY FOR EACH MESH IN ONE SUPER ELEMENT
               
              NTD = 6*(NT-1) + 1
              ND  = 6*(KEL-1) + 1
              EDISC(ND:ND+5) = TEDIS(NTD:NTD+5,1)            
              
          ENDDO 
              EDISIC(1:NEF) = EDISC(1:NEF)
      ENDIF
          
          JC = 0
          DO KEL = 1,4
              NT = ICONN(KEL,JEL) !READ NODE BY CONNECTIVITY FOR EACH MESH IN ONE SUPER ELEMENT
              JC = 3*(KEL-1) + 1
              COORDMESH(JC:JC+2) = XYZM(NT,1:3)   
              
              DO IS = 1,6
                  IDO = IDSUPER(NT,IS)
                  IEDOF(6*KEL-(6-IS)) = IDO 
              ENDDO    
          ENDDO
 
      SELECT CASE(ITYPE)
	CASE(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20,21,17)
	CALL MODULE (PROPM(1,MSET),PROPG(1,ISET),PROPOF,
	1              NODEX(1,IEL),WOREL,AMV(1,NV),SAS,
	2	            COORDMESH,EDISC,EDISIC,ELOD,ALPHA,SEL,SEDI,FIN,
	3	            IPINSF,ISET,MSET,HINFC,FIXEN,FIXLR,
	4	            FIXEO,FIXLO,CABFF,SEQ,HINFP,TCHET(1,ILC),
	5	            TAMBT(1,ILC),PMATRL,ISFAC(1,IEL))
      END SELECT
      
C      DELOD = DELOD + ELOD
      
5     NEID = 0
      DO I = 1,NEQ
          EDI = EDIS(I)
          IF(EDI.GT.0) NEDI = NEID + 1
      ENDDO          
      
C     ASSEMBLE STRESS
      WORELM(1:MWA) = WORELM(1:MWA) + WOREL
C     --------------------------------------
C     FOR ASSEMBLING THE EACH ELEMENT IN ONE SUPER ELEMENT
      DO IEF = 1,24
          IEQ = IEDOF(IEF)
          DO JEF = 1,24
              JEQ = IEDOF(JEF)
              SKS(IEQ,JEQ) = SKS(IEQ,JEQ) + SM(IEF,JEF)
          ENDDO    
      ENDDO              
          
C     ---------------------------------------------------- 
      ENDDO !=============================================END OF LOOP FOR 4X4 MESHES IN ONE SUPER ELEMENT     

      
C     STIFFNESS REDUCTION FACTOR      
      FACSR = 1.0D0!0.98
      SKS = FACSR*SKS      

      IF(ITASK.EQ.5) GO TO 10 !SKIP THE FICTITIOUS BEAM ROUTINE FOR MASS MATRIX
      
C     FOR CALCULATION OF HEIGHT OF FICTITIOUS STIFF BEAMS 
      IL = 0
      DO I = 1,4
         INEL = INTERELE(IEL,I)
         IF(INEL.NE.0) IL = IL + 1
      ENDDO          
      
       ! CALCULATE HEIGTH OF EACH ELEMENT
       HEIGHT1 = ABS(XYZM(1,NGRAV)-XYZM(2,NGRAV))
       HEIGHT2 = ABS(XYZM(2,NGRAV)-XYZM(3,NGRAV))
       
       IF(HEIGHT1.GT.HEIGHT2) THEN
           HEIGHT = HEIGHT1
       ELSEIF(HEIGHT2.GT.HEIGHT1) THEN           
           HEIGHT = HEIGHT2
       ELSEIF(HEIGHT1.EQ.HEIGHT2) THEN           
           HEIGHT = HEIGHT1           
       ENDIF           
       HEIGHT = HEIGHT
      DO JEL = 1,IL !===================================== START OF LOOP FOR FICTITIOUS STIFF BEAMS
          INEL = INTERELE(IEL,JEL)
          CALL FICBEAM(INEL,XYZM,PROPM(1,MSET),PROPG(1,ISET),STIFFA,NODEMAT)
          
          DO KEL = 1,5
              NT = NODEMAT(KEL) !READ NODE BY CONNECTIVITY FOR EACH MESH IN ONE SUPER ELEMENT
              
              DO IS = 1,6
                  IDO = IDSUPER(NT,IS)
                  IEDOFB(6*KEL-(6-IS)) = IDO 
              ENDDO    
          ENDDO
          
          !COMBINE STIFFNESS MATRIX OF SUPER ELEMENT WITH THAT OF FICTITIOUS SIFF BEAMS
          DO IEF = 1,30
              IEQ = IEDOFB(IEF)
              DO JEF = 1,30
                  JEQ = IEDOFB(JEF)
                    SKS(IEQ,JEQ) = SKS(IEQ,JEQ) + STIFFA(IEF,JEF)
CB                  STIFFAMAT(IEQ,JEQ) = STIFFAMAT(IEQ,JEQ) +  STIFFA(IEF,JEF)
              ENDDO
          ENDDO
  
      ENDDO !===================================== END OF LOOP FOR FICTITIOUS STIFF BEAMS      

      SKSSUP(1:30,1:30) = SKS(1:30,1:30) ! STIFFNESS MATRIX THAT BOUNDARY CONDITIONS ARE APPLIED
      
C     FOR APPLYING BOUNDARY CONDITIONS TO THE NODES OF FINE MESH  ---------------------------------------         
          ALLOCATE(ISUPPT(3,7))
          
          CALL CLEARMATI(ISUPPT,3,7)
          
          LBC = 0
          DO INC = 1,5
              KBC = 0
CC              NODEA = INCINEW(INC,IEL)   !SONGSAK SUPPRESS HERE  ... NEED TO CHECK BECAUSE MEMORY WILL OVER STORAGE FRO INC GT 4
CC              NODEB = INCINEW(INC+1,IEL) !SONGSAK SUPPRESS HERE  ... NEED TO CHECK BECAUSE MEMORY WILL OVER STORAGE FRO INC GT 4
              JNC1 = INC
              JNC2 = INC+1
              IF(INC.EQ.5) THEN
                  JNC = 1
CC                  NODEA = INCINEW(JNC+3,IEL) !SONGSAK SUPPRESS HERE  ... NEED TO CHECK BECAUSE MEMORY WILL OVER STORAGE FRO INC GT 4
CC                  NODEB = INCINEW(JNC,IEL)   !SONGSAK SUPPRESS HERE  ... NEED TO CHECK BECAUSE MEMORY WILL OVER STORAGE FRO INC GT 4
                  JNC1 = JNC+3
                  JNC2 = JNC
              ENDIF    
              DO IBS = 1,NBS
                  NODEC = IBOUND(IBS,1)
                  IF(NODEA.EQ.NODEC) THEN
                      KBC = KBC + 1
                      GO TO 6   
                  ENDIF
              ENDDO            
              
6            DO JBS = 1,NBS
                  NODED = IBOUND(JBS,1)
                  IF(NODEB.EQ.NODED) THEN
                      KBC = KBC + 1
                      IF(KBC.EQ.2) THEN
                          LBC = LBC + 1
                          NSAVE(LBC,1) = JNC1
                          NSAVE(LBC,2) = JNC2
                      ENDIF    
                      GO TO 7
                  ENDIF
              ENDDO    
7          ENDDO

          DO ILB = 1,LBC
              NOD1 = NSAVE(ILB,1)
              NOD2 = NSAVE(ILB,2)
              
              ! SAVE CONSTRAINTED DOF IN SUPPT MATRIX
              ISUPPT(1,2:7) = IBOUND(ILB,2:7)  !THIS ONE SHOULD BE READ FROM IBOUND
                                                   !IT SHOULD BE CONSIDERED FOR DIFFERENT BOUNDARY CONDITIONS OF TWO NODES
                                                      ! EX) NODE 4 : ROLLER, NODE 1:HINGE (#1)
              ISUPPT(2,2:7) = IBOUND(ILB,2:7)
              ISUPPT(3,2:7) = IBOUND(ILB,2:7)
              
              !CHECK THE TWO BOUNDARY CONDITIONS IN ONE SUPER ELEMENT ---- (#1)
              IF(NOD1.EQ.1.AND.NOD2.EQ.2) THEN
                  ISUPPT(1,1) = 9
                  ISUPPT(2,1) = 5
                  ISUPPT(3,1) = 10
              ELSEIF(NOD1.EQ.2.AND.NOD2.EQ.3) THEN
                  ISUPPT(1,1) = 11
                  ISUPPT(2,1) = 6
                  ISUPPT(3,1) = 12
              ELSEIF(NOD1.EQ.3.AND.NOD2.EQ.4) THEN
                  ISUPPT(1,1) = 13
                  ISUPPT(2,1) = 7
                  ISUPPT(3,1) = 14                   
              ELSEIF(NOD1.EQ.4.AND.NOD2.EQ.1) THEN
                  ISUPPT(1,1) = 15
                  ISUPPT(2,1) = 8
                  ISUPPT(3,1) = 16                   
              ENDIF    
          ENDDO
C     ------------------------------------------------------------------------------------------------------------------
               
C     FOR APPLYING BOUNDARY CONDITIONS ASSEMBLED STIFFNESS MATRIX OF FINE MESHES ---------------------------------------          
           DO IASUPPT = 1,3
               IPP = ISUPPT(IASUPPT,1)
               DO IPG = 1,6
                  M = ISUPPT(IASUPPT,IPG+1)                          
                  NPP = IDSUPER(IPP,IPG)
                  IF(M.EQ.1) THEN
                      SKSSUP(NPP,1:150) = 0.0
                      SKSSUP(1:150,NPP) = 0.0
                      SKSSUP(NPP,NPP)   = 1.0
                  ENDIF    
              ENDDO    
           ENDDO 
C      ENDIF       
      DEALLOCATE(ISUPPT)
C     -----------------------------------------------------------------------------------------------------------------        
C     CALL FOF STATIC CONDENSATION OF STIFFNESS MATRIX ----------------------------------------------------------------
10    LLVC = LLV+NEQ*(ILC-1) !INDICATOR FOR ASSEMBLED LOAD A(LLVC)      
      IF(ITASK.EQ.5) THEN
         SKSC = MATMUL(TRANSPOSE(TMACON),MATMUL(SKS,TMACON))
      ELSE    
          CALL MATRIX_CONDEN(SKS,SKSC,A(LLVC),TMACON,150,NEQ,ICONN,0)
          CALL MATRIX_CONDEN(SKSSUP,SKSCSUP,A(LLVC),TMACON,150,NEQ,ICONN,0)

      ENDIF          
      
      IF(ITASK.EQ.5) GO TO 20 !SKIP THE FICTITIOUS BEAM ROUTINE FOR MASS MATRIX
C     ------------------------------------------------------------------------------    
      
C     SUBSTRACT THE STIFFNESS OF FICTITIOUS BEAM -----------------------------------------------------       
      DO JEL = 1,IL !===================================== START OF LOOP FOR FICTITIOUS STIFF BEAMS
          INEL = INTERELE(IEL,JEL)
          CALL FICBEAM_SUB(INEL,XYZM,PROPM(1,MSET),PROPG(1,ISET),STIFFBEAM,NODEMAT)
          
          DO KEL = 1,2
              NT = NODEMAT(KEL) !READ NODE BY CONNECTIVITY FOR EACH MESH IN ONE SUPER ELEMENT
              DO IS = 1,6
                  IDO = IDSUPER(NT,IS)
                  IEDOFC(6*KEL-(6-IS)) = IDO 
              ENDDO    
          ENDDO
          
          !SUBSTRACT STIFFNESS MATRIX OF FICTIOUS BEAM ELEMENT 
          DO IEF = 1,12
              IEQ = IEDOFC(IEF)
              DO JEF = 1,12
                  JEQ = IEDOFC(JEF)
                  SKSC(IEQ,JEQ) = SKSC(IEQ,JEQ) - STIFFBEAM(IEF,JEF)
                  SKSCSUP(IEQ,JEQ) = SKSCSUP(IEQ,JEQ)- STIFFBEAM(IEF,JEF)
              ENDDO
          ENDDO
  
      ENDDO      

C     ------------------------------------------------------------------------------   
C     ADD THE STIFFNESS OF BEAM THAT MOMENT ABOUT Z-AXIS IS RELESED -----------------------------------------------------    
      DO II = 1,NELE   
          INTERELE2(II,1) = INTERELE3(II,1) 
          INTERELE2(II,2) = INTERELE3(II,2) 
          INTERELE2(II,3) = INTERELE3(II,1) 
          INTERELE2(II,4) = INTERELE3(II,2)      
          
         !THE CASE THAT ONLY THETA Z MOMENT APPLIED
          !INTERELE2(II,1) = 1
          !INTERELE2(II,2) = 2
          !INTERELE2(II,3) = 3
          !INTERELE2(II,4) = 4
          !INTERELE2(II,5) = 5
          !INTERELE2(II,6) = 6
          !INTERELE2(II,7) = 7
          !INTERELE2(II,8) = 8         
      ENDDO  
      
      DO JEL = 1,4 !===================================== START OF LOOP FOR STIFF BEAMS THAT MOMENT ABOUT Z-AXIS IS RELESED
         !DO JEL = 1,8 !THE CASE THAT ONLY THETA Z MOMENT APPLIED
          !IRLELEASE = 4
          
          INEL = INTERELE2(IEL,JEL)
          
          IF(JEL.EQ.1.OR.JEL.EQ.5) THEN    
              IRELEASE = 1
          ELSEIF(JEL.EQ.2.OR.JEL.EQ.6) THEN
              IRELEASE = 1
          ELSEIF(JEL.EQ.3.OR.JEL.EQ.7) THEN
              IRELEASE = 2
          ELSEIF(JEL.EQ.4.OR.JEL.EQ.8) THEN
              IRELEASE = 2              
          ENDIF

          CALL FICBEAM_SUB2(INEL,XYZM,PROPM(1,MSET),PROPG(1,ISET),STIFFBEAM,NODEMAT,IRELEASE)
          
          DO KEL = 1,2
              NT = NODEMAT(KEL) !READ NODE BY CONNECTIVITY FOR EACH MESH IN ONE SUPER ELEMENT
              DO IS = 1,6
                  IDO = IDSUPER(NT,IS)
                  IEDOFC(6*KEL-(6-IS)) = IDO 
              ENDDO    
          ENDDO
          
          !ADD STIFFNESS MATRIX OF BEAM
          DO IEF = 1,12
              IEQ = IEDOFC(IEF)
              DO JEF = 1,12
                  JEQ = IEDOFC(JEF)
                  SKSC(IEQ,JEQ) = SKSC(IEQ,JEQ) + STIFFBEAM(IEF,JEF)
                  SKSCSUP(IEQ,JEQ) = SKSCSUP(IEQ,JEQ) + STIFFBEAM(IEF,JEF)
              ENDDO
          ENDDO
  
      ENDDO          
C     --------------------------------------------------------------------------------------------         
      !SAVE CONDENSED STIFFNESS MATRIX TO CALCULATE THE REACTION FORCES
C    IF(ITASK.EQ.1.OR.ITASK.EQ.5) THEN
      IF(ITASK.EQ.1) THEN          
          IONE = 0
          DO I = 1,24
              DO J = I,24
                  IONE = IONE + 1
                  SKSCONE(IONE) = SKSC(I,J)
              ENDDO    
          ENDDO                
          KREC2 = 5004
            WRITE(KREC2,1000) SKSCONE(1:300)
      ELSEIF(ITASK.EQ.2.OR.ITASK.EQ.3) THEN   !READ CONDENSED STIFFNESS MATRIX TO CALCULATE THE REACTION FORCES
            KREC2 = 5004
            READ(KREC2,1000) SKSCONE(1:300)
          IONE = 0
          DO I = 1,24
              DO J = I,24
                  IONE = IONE + 1
                  SKCSAVE(I,J) = SKSCONE(IONE)  
              ENDDO    
          ENDDO

          DO I = 1,24
              DO J = 1,24
                  IF(I.NE.J) SKCSAVE(J,I) = SKCSAVE(I,J)  
              ENDDO    
          ENDDO
      ENDIF  
C     -------------------------------------------------------------------------------------------------       
C     DIRECT CALCULATION OF ELOD (MEMEBER FORCES)
          ALPHA = 1.0D0
          BETA  = 1.0D0
            
          M = 24
          N = 1
          K = 24                        
             
          CALL CLEARMAT(DELOD,M,N)
          CALL DGEMM('N','N',M,N,K,ALPHA,SKCSAVE,M,EDISNEW,K,BETA,DELOD,M)
              ELOD(1:NEF) = DELOD(1:NEF)   
C     --------------------------------------------------------------------------------         
      
        IF(ITASK.EQ.1) THEN
          WRITE(KREC,1000) TMACON(1:150,1:24)     
      ENDIF   

      
C     ASSEMBLE STRESS
      WOREL(1:MWA) = WORELM(1:MWA) 
C     --------------------------------------
  
C     SAVE CONDENSED STIFFNESS MATRIX IN TO S MATRIX      
20    IF(ITASK.EQ.5) THEN
          CALL STIFFONE(SKSC,S)
      ELSE
         CALL STIFFONE(SKSC,S)!CALL STIFFONE(SKSCSUP,S)    
      ENDIF

C ============================== END FOR SUPER ELEMENT ==============================            
      
C     FOR LINK CONSTRAINTS BY BJ
      IF(IFLOOR.EQ.3) THEN           
        DO I = 1,NWK
          TAS = TAS + ABS(S(I))
        ENDDO	   
      ENDIF
C     ----------------------------  
      
 
	IF(ITYPE.EQ.17) THEN !TENDON STORE DATA
	CALL TDINIT(IA(LTDN),WOREL,IEL,ILC,KSTEP,MSET,ISET,'RECD') 
	ENDIF

C	===================================================================
C	HERE SUBTRACT INITIAL STRESS FROM THERMAL LOAD SONGSAK MAR2007
	CALL INISUBT(ITYPE,WOREL,KEG,IEL,NPT,NWG,ILC,RHO)
C	===================================================================
C	STORE WORKING ARRAY TO NEW STORAGE 
	CALL ADREWT(KEG,IEL,WOREL,'WRT')
C	===================================================================

C	====================================================
C	UPDATE EAS VARIABLE SONGSAK MAR2006
C	MODIFIED JAN2007
C	====================================================
	IF(IESPT.GT.0.AND.ITASK.EQ.1) THEN 
	CALL EASSTOR(SEDI,SEL,SEQ,HINFC,HINFP,ASEL,ASEQ,AINF,AINP,
	1			 ADIS,EDIS)
	CALL ADREAS(KEG,IEL,APHI,ASEL,APHA,AINF,ADIS,ASEQ,AINP,
	1			'WRT')
	ENDIF 

C	============================================================
C	TRANSFER THE INTERNAL FORCE AND STIFFNESS TO LOCAL SUPPORT
C	============================================================
	IF (NLS.GT.0) CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IEL),
     1                           S,EDIS,EDISI,ELOD,NSF,NNF,2)


C	============================================================
C	ASSEMBLE THE RESTRAINED STIFFNESS FOR SETTLEMENT LOAD CALC.
C	============================================================
	IF(IFREF.EQ.0) THEN
	IF(ITASK.LE.4) THEN
	CALL SETLODR(SETIF,IA(LID),IDRCT,LM(1,IEL),LMRCT(1,IEL),S,NEF)
	ENDIF
	ENDIF

C	============================================================
C	ASSEMBLE THE REACTION FORCE (IN THE SKEW SUPPORT FOR NLS NE 0)
C	============================================================

	SELECTCASE (ITYPE)
	
	CASE(5)  !FIXEND IN THe ELEMENT LOCAL COORDINATE

	IF (NLS.GT.0) CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IEL),
     1                           S,EDIS,EDISI,FIXLR,NSF,NNF,5)       !VARY LOAD
	IF (NLS.GT.0) CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IEL),   
     1                           S,EDIS,EDISI,FIXLO,NSF,NNF,5)       !CONT LOAD

	DO IEF=1,NEF
	FIN(IEF) = FIN(IEF) - RHO*FIXEN(IEF) - FIXEO(IEF)					!ELEMENT LOCAL
	IRC = LMRCT(IEF,IEL)
	IF (IRC.NE.0) REAC(IRC) = REAC(IRC) + ELOD(IEF) 
	1						- RHO*FIXLR(IEF) - FIXLO(IEF)				!SKEW SUPPORT
	ENDDO
	
	CASE(2,6,9,8,10,11,15,16)  !FIXEND IN GLOBAL SYSTEM

	FIXLR(1:NEF) = FIXEN(1:NEF)
	IF (NLS.GT.0) CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IEL),
     1                           S,EDIS,EDISI,FIXLR,NSF,NNF,5)
	FIXLO(1:NEF) = FIXEO(1:NEF)
	IF (NLS.GT.0) CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IEL),
     1                           S,EDIS,EDISI,FIXLO,NSF,NNF,5)

	DO IEF=1,NEF
	FIN(IEF) = FIN(IEF) - RHO*FIXEN(IEF) - FIXEO(IEF)                   !GLOBAL SYSTEM
	IRC = LMRCT(IEF,IEL)
	IF (IRC.NE.0) REAC(IRC) =  REAC(IRC) + ELOD(IEF) 
	1						  -RHO*FIXLR(IEF) - FIXLO(IEF)              !SKEW SUPPORT
	ENDDO

	CASE DEFAULT

	DO IEF=1,NEF
	IRC = LMRCT(IEF,IEL)
	IF (IRC.NE.0) THEN
	REAC(IRC) = REAC(IRC) + ELOD(IEF)
	ENDIF
	ENDDO
	
	ENDSELECT 

C	============================================================
C	CABLE PRETENSION OPTIMIZATION
C	============================================================
	IF(ITYPE.EQ.2.AND.ISTYP.EQ.4.AND.KBOPZ.EQ.1) THEN
	CABDX = COORD(4) - COORD(1)
	CABDY = COORD(5) - COORD(2)
	CABDZ = COORD(6) - COORD(3)
	CABLEN = SQRT(CABDX*CABDX + CABDY*CABDY + CABDZ*CABDZ)
	CABVEC(1) = CABDX/CABLEN
	CABVEC(2) = CABDY/CABLEN
	CABVEC(3) = CABDZ/CABLEN
	CABVEC(4) =-CABDX/CABLEN
	CABVEC(5) =-CABDY/CABLEN
	CABVEC(6) =-CABDZ/CABLEN

C	ALL FORCE COMPONENT
	DO IEF = 1,NEF
	IEQ = LM(IEF,IEL)
	IF (IEQ.NE.0) CABFZ(IEQ,IEL) = CABVEC(IEF)
	ENDDO

C	ONLY X-COMPONENT
	IEQ = LM(1,IEL)
	IF (IEQ.NE.0) CABFX(IEQ,IEL) = CABVEC(1)
	IEQ = LM(4,IEL)
	IF (IEQ.NE.0) CABFX(IEQ,IEL) = CABVEC(4)
	ENDIF

	IF(KSTEP.EQ.1.AND.KBOPZ.EQ.1) THEN
C	DEFINED MEMBER FORCE CONSTRAINT ELEMENT GROUP
	CALL MOMGRP(IA(LGID),TPOPZ(LM557),MEL,KEG,NELE)
	ENDIF

	IF(KSTEP.GT.1.AND.KBOPZ.EQ.1) THEN
C	MEMBER FORCE CONSTRAINT
	CALL MOMCON1(IA(LGID),FIN,NEF,NELE,TPOPZ(LM557),
	1	        TPOPZ(LM557+4*MCMOM),MEL,KEG)
	ENDIF

C	============================================================
	IF(ITYPE.EQ.5.AND.ISPRI.EQ.0) THEN		
	IF(LAN_PRN.EQ.0) CALL SMHFAC(RHO,FIN,IEL)
	IF(LAN_PRN.EQ.1) CALL LNHFAC(FIN,KEG,RHO,IEL,NELE,LCS)
	ENDIF
C	============================================================
      IF (ITASK-IFPR(5).GT.0 .OR. IFPRI.NE.0) GOTO 300
      IF (ITASK.LE.3 .AND. IFREF.NE.0) GOTO 300
      CALL MATOUT (S,NEF,NEF,5,22,'E',15,10,3,' STIFFNESS')


C
 300  IF (ITASK.LE.2) GOTO 310
      IF (IFEIG.EQ.0) GOTO 510

      GOTO 900

 310  CONTINUE
C     -----------------------------------------------------
C     ASSEMBLE THE RESISTING FORCE
C     -----------------------------------------------------
      DO 400  IEF=1,NEF
      IEQ = LM(IEF,IEL)
 400  IF (IEQ.NE.0) RE(IEQ)   = RE(IEQ) - ELOD(IEF)

C     -----------------------------------------------------
C     ASSEMBLE TANGENTIAL STIFFNESS (IFREF=0)
C     ASSEMBLE GEOMETRIC STIFFNESS MATRIX (ITASK=4,IFEIG=0)
C     -----------------------------------------------------
 500  IF (IFREF) 900,510,900
C     ---------------------------
C     ASSEMBLE LUMPED MASS MATRIX
C     ---------------------------
 510  CALL CPU_TIME (TIM1)

C     ------------------------------------------------------------
C     ASSEMBLE TANGENTIAL STIFFNESS, GEOMETRIC STIFFNESS
C     OR CONSISTENT MASS MATRIX (NBLOCK=1),OR DAMPING MATRIX
C	Damping Matrix is assembled here(ITASK=6),(14Nov03,NguyenDV)
C     ------------------------------------------------------------
 520  CALL ADDBAN (LM(1,IEL),IA(LMA),S,NEQF,NEQL,NPRE,NEF,NWS)
	

c	K = 0
c	DO I =1,NEF
c	DO J =I,NEF
c	K = K + 1
c	IF(I.EQ.11.and.j.eq.11) WRITE(*,*) MEL,K,S(K)
c	IF(I.EQ.11.and.j.eq.14) WRITE(*,*) MEL,K,S(K)
c	IF(I.EQ.14.and.j.eq.14) WRITE(*,*) MEL,K,S(K)
c
c	IF(I.EQ.4.and.j.eq.4) WRITE(*,*) MEL,K,S(K)
c	IF(I.EQ.4.and.j.eq.7) WRITE(*,*) MEL,K,S(K)
c	IF(I.EQ.7.and.j.eq.7) WRITE(*,*) MEL,K,S(K)
c	IF(I.EQ.J) WRITE(*,*) MEL,K,S(K)
c	ENDDO
c	ENDDO
C	WRITE(*,*) S(1:(NEF*NEF+NEF)/2)

 590  CALL CPU_TIME (TIM2)
      TIM(13) = TIM(13) + TIM2-TIM1

 900  CONTINUE

C	STOP
C     ----------------------------------------------
C     STRESS OUTPUT (ISPRI=0)
C     ----------------------------------------------
      CALL CPU_TIME (TIM1)
      TIM(11) = TIM(11) + (TIM1-TIM3)

      IF (ISPRI.NE.0) RETURN

C	PRINT OUT ELEMENT STRESS
	CALL ELPOUT


      CALL CPU_TIME (TIM2)
      TIM(17) = TIM(17) + (TIM2-TIM1)

1000  FORMAT(F20.5)    
      RETURN
      END
C
C=====================================================================      
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SHELL4SUPER(PROPM,PROPG,NODEX,WA,AMV,S,COORD,
	1					EDIS,EDISI,RE,MWG,MSET,SEDI,SEL,ALPHA,HINFC)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /DYNA/  CDEN,IMASS
	COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST

	COMMON /MMENH/ MM,MM1,MM2,NDIMC
      
      
C     FOR SUPER ELEMENT BY BJ
      COMMON /SUPERSTIF/ SM(24,24) 
      COMMON /JNELE/ JEL
      COMMON /REE/RE1(24),RE2(24)
      
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL   
      
      DIMENSION NCOLM(24)     

      DIMENSION PROPM(1),PROPG(1),NODEX(1),WA(MWG,1),S(300),COORD(3,4)
      DIMENSION EDIS(24),EDISI(24),RE(24),REDIS(24),COORDI(3,4)
      DIMENSION DR(8,8),H(8),HD(2,8),XJI(4),HR(8),HS(8),B(240)
      DIMENSION BDRL(48),DISD(12),EPS(8),EPSQ(8),SIGR(8)
      DIMENSION VR(3),VS(3),VT(3),SL(4,2),FF(6),NOD(4)
      DIMENSION COVR(3),COVS(3)
	DIMENSION BA(4,120),FJ(4),RRN(4),SSN(4)
	DIMENSION AMV(3)
      
C      DIMENSION RE2(NEF)
C     --------------
C       EAS METHOD
C     --------------

	DIMENSION GE(8,MM),DRM(8,8),REAS(24),SEAS(24,24)
	DIMENSION XJO(4),FJO(4),TTO(8,8),BB(8,24)
	DIMENSION EAS(8),ALPHA(MM)
      DIMENSION SED(MM,MM),SEL(MM,24),SEDI(MM,MM),RH(MM)
	DIMENSION HINFC(MM)
           
	
      DATA SL /2.5E5,.1334,.1334,.1334,1.E20,19.2,4.00,4.00/
C
C     EAS FLAG
      LEAS = 0
      IF(NLOPT.EQ.3) LEAS = 0
      IF(MTMOD.NE.1) LEAS = 0

      B(1:240) = 0.0D0

C     DRILLING RIGIDITY FACTOR
      FACTD = 0.0D0!0.0D0!10.0
      IF(LEAS.EQ.1) FACTD = 1.0D0!0.01D0!5.0E-05
            
      PI=3.1415926535898
      CALL CLEARA (BDRL,48)
C     ---------------------------------
C     OBTAIN LINEAR STRESS - STRAIN LAW
C     mtmode=1, ipel=1 for linear analysis
	IPEL=1
      If (MTMOD.NE.2) CALL HOKLAW (PROPM,PROPG,2)
C	here 2 indicates plane stress parameters

C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C     ------------------------------------------------------
      CALL SHINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)
C	this subroutine updated coordinate adding displacement for nonlinear 
C	analysis Loop over sampling points, four points for transverse strains
      RRN(1)=1.0D0
	SSN(1)=0.0D0
	RRN(2)=-1.0D0
	SSN(2)=0.0D0
	RRN(3)=0.0D0
	SSN(3)=1.0D0
	RRN(4)=0.0D0
	SSN(4)=-1.0D0
C	make all assumed strain Bs values zero at the start
      DO 5 I=1,4
	FJ(I)=0.0D0
	XJI(I)=0.0D0
	DO 5 J=1,120
 5	BA(I,J)=0.0D0
 	NSP=4
      DO 6 I=1,NSP
	CALL SHAP2D (RRN(I),SSN(I),H,HD,NODEX,NNO) 
C	nodex IS NOT USED here for four node element, 
      CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FJ)
C	unit vectors and jacobians are calculated at sampling points
	CALL SHBMATS (NNO,H,HD,VR,VS,VT,XJI,HR,HS,BA,I,FJ)
 6    CONTINUE
C	BA matrix at 4 nodes are completed

C     GET TRANSFORMATION MATRIX AT CENTER (EAS)
      IF(LEAS.EQ.1) THEN
	RN = 0.0D0
	SN = 0.0D0
	CALL SHAP2D (RN,SN,H,HD,NODEX,NNO) 
      CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJO,DETO,RR,SS,SNA,1,FJO)
	CALL TNEAS (FJO,TTO)
      ENDIF
      
C     INITIALIZE EAS TERMS
	SED = 0.0D0
	SEDI = 0.0D0
	SEL = 0.0D0
	RH  = 0.0D0
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
      IPT=0
      DO 900  IGR=1,NGR
      RN = GLOC(IGR,NGR)
      DO 900  IGS=1,NGS
      SN = GLOC(IGS,NGS)
      WT = GWT(IGR,NGR)*GWT(IGS,NGS)
      IPT= IPT+1

C     -----------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C     -----------------------------------------------------
      CALL SHAP2D (RN,SN,H,HD,NODEX,NNO)
C     ----------------------------------------------------
C     INVERSE JACOBIAN (XJI) , JACOBIAN coefficient F, 
C	DETERMINANT (DET) AND STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------------------
      CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FJ)
      DVOL=WT*DET
	
	CALL SHSTORL1 (COORD,PROPG,WA(1,IPT),NWG,VR,VS,VT)	!STORE SHELL LAX FOR STRESS CALCULATION JAN09

C     GET GE MATRIX IN PHYSICAL COORDINATE (EAS)
	IF(LEAS.EQ.1) CALL MTEAS(RN,SN,TTO,DETO,DET,GE,MM)
	
	
      TH=PROPG(2)
      
      CALL RELFILL('TWGH',DVOL*TH*PROPM(5),1,KEG,2)
      
      
	IF (ITASK.NE.5) GOTO 140
C     ---------------------
C     MASS MATRIX (ITASK=5)
C     ---------------------
C	CONMSS ADDED BY SONGSAK FOR RC SHELL
	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6) THEN
	    CALL CONMSS(TH,MSET,RHORC)
	    CALL SHMASS (S,H,VR,VS,DVOL,RHORC,TH,IMASS,NNO,NEF)
	ELSEIF(MTMOD.EQ.2) THEN
	    PROPM(5) = CDEN
          CALL SHMASS (S,H,VR,VS,DVOL,PROPM(5),TH,IMASS,NNO,NEF)
	ELSE
          CALL SHMASS (S,H,VR,VS,DVOL,PROPM(5),TH,IMASS,NNO,NEF)
	ENDIF
      GOTO 900

      
  140 IF(IFLOOR.EQ.3) THEN !FOR SHEAR WALL ELEMENT FROM SUPER ELEMENT
          CALL SHBMAT11 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,RN,SN,BA,NEF)
          CALL SHBMAT22 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B)
      ELSEIF(IFLOOR.NE.3) THEN !FOR JUST SUPER ELEMENT
          CALL SHBMAT1 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,RN,SN,BA,NEF)
          CALL SHBMAT2 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B)          
      ENDIF
      
      BB(1:8,1:24) = 0.0D0
      K=1
      M=1
      DO I = 1,4
      DO J = 1,3
      L=K+3
      BB(1,K) = B(M+0)
      BB(2,K) = B(M+1)
      BB(3,K) = B(M+2)
      BB(7,K) = B(M+3)
      BB(8,K) = B(M+4)
      BB(4,L) = B(M+15)
      BB(5,L) = B(M+16)
      BB(6,L) = B(M+17)
      BB(7,L) = B(M+18)
      BB(8,L) = B(M+19)
      K=K+1
      M=M+5
      ENDDO
      K=K+3
      M=M+15
      ENDDO
      
C     -------------------------------------------
C     SHEAR LOCKING CONSTRAINT FACTORS (SLR,SLS), 
C	MODIFIED WITH ASSUME STRAIN
C     -------------------------------------------
      SLR=1.
	SLS=1.
C	--------------------------------
C	DETERMINE INITIAL MATERIAL ANGLE
C	added by gislon - sept2002
C	--------------------------------
      SELECTCASE(MTMOD)
      CASE(2,5,6)
	CALL CANGLE(VR,VS,AMV,ANG)
	ENDSELECT

C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES (IFREF.EQ.0)
C     -----------------------------------------
      IF (MTMOD.NE.2) CALL SHDELA (DR,DVOL,SLR,SLS)

C     EAS MATRIX      
      IF(LEAS.EQ.1) THEN
	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(DR,GE))
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(DR,BB))
      ENDIF
      
      K=1
      M=1
      DO I = 1,4
      DO J = 1,3
      L=K+3
C      SEL(1,K) = SEL(1,K)
C      SEL(2,K) = SEL(2,K)
C      SEL(3,K) = SEL(3,K)
C      SEL(7,K) = 0.0D0
C      SEL(8,K) = 0.0D0

C      SEL(4,L) = 0.0D0
C      SEL(5,L) = 0.0D0
C      SEL(6,L) = 0.0D0
C      SEL(7,L) = 0.0D0
C      SEL(8,L) = 0.0D0
      
      K=K+1
      M=M+5
      ENDDO
      K=K+3
      M=M+15
      ENDDO      
      
C      IF (NLOPT+ITASK.EQ.1) GOTO 700
C     -------------------------------------------------
C     REMOVE RIGID BODY TRANSLATIONS AND ROTATIONS FROM
C     TOTAL DISPLACEMENT VECTOR (NLOPT>1)
C     -------------------------------------------------
      IF (NLOPT.LE.1) GOTO 180
      CALL SHMDSP (COORD,COORDI,EDIS,REDIS,H,HD,VR,VS,VT,NNO)
C     ------------------------
C     DISPLACEMENT DERIVATIVES
C     ------------------------
  180 CALL CLEARA (DISD,12)
      K=1
      M=1
      DO 200 I=1,NNO
      DO 190 J=1,3
      L=K+3
      DISD(1)=DISD(1)+B(M)*REDIS(K)	
      DISD(2)=DISD(2)+B(M+1)*REDIS(K)
      DISD(3)=DISD(3)+(B(M+2)-HR(I)*VS(J))*REDIS(K)
      DISD(4)=DISD(4)+(B(M+2)-HS(I)*VR(J))*REDIS(K)
      DISD(5)=DISD(5)+B(M+15)*REDIS(L)
      DISD(6)=DISD(6)+B(M+16)*REDIS(L)
      DISD(7)=DISD(7)+B(M+17)*REDIS(L)
      DISD(9)=DISD(9)+B(M+3)*REDIS(K)
      DISD(10)=DISD(10)+B(M+18)*REDIS(L)
      DISD(11)=DISD(11)+B(M+4)*REDIS(K)
      DISD(12)=DISD(12)+B(M+19)*REDIS(L)
      K=K+1
 190  M=M+5
      K=K+3
 200  M=M+15
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      EPS(1) = DISD(1)
      EPS(2) = DISD(2)
      EPS(3) = DISD(3)+DISD(4)
      EPS(4) = DISD(5)
      EPS(5) = DISD(6)
      EPS(6) = DISD(7)
      EPS(7) = DISD(9) +DISD(10)
      EPS(8) = DISD(11)+DISD(12)
      
      IF(LEAS.EQ.1) THEN
	EAS = MATMUL(GE,ALPHA)
	DO I=1,8
        EPS(I)= EPS(I) - EAS(I)
	ENDDO
	ENDIF
	
C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400
      EPSQ(1)=0.5D0*(DISD(1)*DISD(1)+DISD(4)*DISD(4)+DISD(9)*DISD(9))
      EPSQ(2)=0.5D0*(DISD(2)*DISD(2)+DISD(3)*DISD(3)+DISD(11)*DISD(11))
      EPSQ(3)=     DISD(3)*DISD(1)+DISD(2)*DISD(4)+DISD(11)*DISD(9)
      EPSQ(7)=     DISD(10)*DISD(1)+DISD(12)*DISD(4)
      EPSQ(8)=     DISD(10)*DISD(3)+DISD(12)*DISD(2)
      DO 300  I=1,3
 300  EPS(I)= EPS(I)-EPSQ(I)
      EPS(7)= EPS(7)-EPSQ(7)
      EPS(8)= EPS(8)-EPSQ(8)
 400  EPS(7)= EPS(7)
      EPS(8)= EPS(8)
      
C     ---------------------------------------
C     GAUSS POINT STRESS - RESULTANTS SIGR(8)
C     ---------------------------------------
      GO TO (410,465,450,460,467,468), MTMOD
 410  CALL MPSIGA (EPS,SIGR)

	CALL SHSTRS(VR,VS,VT,EPS,WA(9,IPT),PROPM,TH)

       DO 420 I=1,8
 420  WA(I,IPT)=SIGR(I)
      GO TO 500

c-----COMPOSITE 
 465  CALL COMRGD(EPS,SIGR,DR,PROPM,TH,ANG,DVOL,'SHELL')
      DO 426 I=1,8
	WA(I  ,IPT) =  EPS(I)
 426	WA(I+8,IPT) = SIGR(I)
      GO TO 500

 450	IF (HP.EQ.0.) THEN
	CALL IVANOV (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,SIGR,DR,DVOL)
      ELSE
	EPCA = 0.
	CALL IVANOVH (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,
	1					SIGR,DR,DVOL,WA(18,IPT),EPCA)
	END IF
	  APEL = WA(17,IPT)
	  IPEL = 2
      GOTO 500

 460	IF (HP.EQ.0.) THEN
	CALL MLAYER (WA(1,IPT),EPS,SIGR,DR,DVOL)
      IPEL=2
      DR(7,7)=DR(7,7)*SLR
      DR(8,8)=DR(8,8)*SLS
	GOTO 500
	ELSE
	CALL MLAYERH (WA(1,IPT),EPS,SIGR,DR,DVOL)
      IPEL=2
      DR(7,7)=DR(7,7)*SLR
      DR(8,8)=DR(8,8)*SLS
	GO TO 500
	END IF

 467	CALL RCLAYR (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
      DR(7,7)=DR(7,7)*SLR
      DR(8,8)=DR(8,8)*SLS
	GO TO 500

 468	CALL RCLAYRE (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
      DR(7,7)=DR(7,7)*SLR
      DR(8,8)=DR(8,8)*SLS
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510 I=1,8
 510  SIGR(I)=SIGR(I)*DVOL


C     FOR BUCKLINGANALYSIS -- GET GEOMETRIX STIFFNESS	
      IF (ISOLOP.EQ.4.AND.IFEIG.EQ.0) GOTO 800
      
      
 520  K=1
      M=1
      DO 560 I=1,NNO
      DO 540 J=1,3
      L=K+3
      RE(K)=RE(K)+B(M)*SIGR(1)+B(M+1)*SIGR(2)+B(M+2)*SIGR(3)
     1                      +B(M+3)*SIGR(7)+B(M+4)*SIGR(8)
      RE(L)=RE(L)+B(M+15)*SIGR(4)+B(M+16)*SIGR(5)+B(M+17)*SIGR(6)
     1                         +B(M+18)*SIGR(7)+B(M+19)*SIGR(8)
      K=K+1
 540  M=M+5
      K=K+3
 560  M=M+15

C     EAS FORCE
	IF(LEAS.EQ.1) RH  = RH  + MATMUL(TRANSPOSE(GE),SIGR)
	
C	---------------------------------------
C	FORCES FROM DRILLING DOF SONGSAK JAN09
      IF(NLOPT.EQ.0) THEN
      CALL SHBDRL (NNO,H,HR,HS,VR,VS,VT,BDRL)
	FDRL=FACTD*DR(3,3)
	DO IEF = 1,NEF
	DO JEF = 1,NEF
	RE(IEF) = RE(IEF) + BDRL(IEF)*FDRL*BDRL(JEF)*REDIS(JEF)
	ENDDO
	ENDDO
	ENDIF
C	---------------------------------------

C     ------------------------------------------------------
C     FIND LINEAR CONTRIBUTION TO STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------------------
      IF (IFREF ) 900,700,900
c 700  CALL SECOND(T1,TIM1)
C     -----------------------------------------------------------
C     REMOVE SINGULARITY FROM STIFFNESS OF EXACTLY PLANE ELEMENTS
C     -----------------------------------------------------------
 700  IJ=1
      N=3*NEF-2
      LROW=NEF-3
      FAC=DMIN1(DR(3,3),DR(4,4)/DET)*1.E-6
      DO 720 I=1,3
      FI=FAC*VT(I)
      DO 720 J=I,3
      FF(IJ)=FI*VT(J)
 720  IJ=IJ+1
      DO 740 I=1,NNO
      HH=H(I)*H(I)
      S(N)=S(N)+HH*FF(1)
      S(N+1)=S(N+1)+HH*FF(2)
      S(N+2)=S(N+2)+HH*FF(3)
      N=N+LROW
      S(N)=S(N)+HH*FF(4)
      S(N+1)=S(N+1)+HH*FF(5)
      N=N+LROW-1
      S(N)=S(N)+HH*FF(6)
      N=N+4*LROW-14
 740  LROW=LROW-6
      CALL SHBDRL (NNO,H,HR,HS,VR,VS,VT,BDRL)
C      CALL SHKLINEAS_SUPER (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD,FACTD)
      CALL SHKLINEAS (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD,FACTD)

C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 810
 800  CALL SHGEO1 (S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)
 810  CONTINUE 
 
 900  CONTINUE
               

      IF (ISOLOP.EQ.4.AND.IFEIG.EQ.0) GOTO 905 !NO NEED TO SUBTRACT EAS FOR GEOMETRIC STIFFNESS FOR BUCKLING ANALYSIS
C     ------------------------------------------------------
C     CORRECTION OF STIFFNESS AND INTERNAL FORCES DUE TO EAS
      IF(LEAS.EQ.1) THEN
	    HINFC(1:MM) = RH(1:MM)
    	
	    CALL INVMATF(SED,SEDI,MM,IB)
	    IF(IB.EQ.1) GOTO 905

	    REAS = MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,RH))
	    RE(1:24) = RE(1:24) - REAS(1:24)
    	
	    SEAS = MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,SEL))
	    K = 0
	    DO I = 1,24
	        DO J = I,24
	            K = K + 1
	            S(K) = S(K) - SEAS(I,J)
	        ENDDO
	    ENDDO
	ENDIF
C     ------------------------------------------------------
905   CONTINUE	

C	KEF = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	ENDDO
C	ENDDO
C	MEF = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	MEF = MEF + 1
C	S(KEF) = S(MEF) 
C	ENDDO
C	ENDDO

C     SAVE THE STIFFNESS MATRIX INTO FULL MATRIX FORM
      II = 25
      DO I = 1,24  
          II = II-1
          NCOLM(I) = II
      ENDDO

      KK = 0
      IJ = 0
      DO II = 1,24
          IJ = IJ + 1
         DO JJ = II,24
             KK = KK + 1
             SM(II,JJ) = S(KK)
         ENDDO    
      ENDDO      
      
      DO II = 1,24
         DO JJ = 1,24
          IF(II.NE.JJ) SM(JJ,II) = SM(II,JJ)
         ENDDO
      ENDDO          
      
      
      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SHKLINEAS_SUPER (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD,FACT)
C
C     ------------------------------------------------------------
C     EVALUATES LINEAR CONTRIBUTION TO TANGENTIAL STIFFNESS MATRIX
C
C     S(1176)  = STIFFNESS MATRIX STORED UPPER-TRIANGULAR ROW-WISE
C     DR(64)   = ELASTO-PLASTIC RIGIDITY MATRIX STORED COLUMN-WISE
C     B(240)   = STRAIN-DISPLACEMENT MATRIX
C     BDRL(48) = DERIVATIVE OPERATORS FOR DRILLING STRAIN
C     NNO      = NUMBER OF NODES FOR ELEMENT
C     NEF      = NUMBER OF D.O.F FOR ELEMENT
C     IPEL     = SECTION PLASTICITY INDICATOR (1=EL,2=EL-PL)
C     ------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(*),DR(64),B(240),BDRL(48)
C
      N=1
      I=4
      FDRL=FACT*DR(19)
C     -----------------------------
C     TRANSVERSE SHEAR CONTRIBUTION
C     -----------------------------
      DO 30 NROW=1,NEF
      A1=DR(55)*B(I) + DR(56)*B(I+1)
      A2=DR(63)*B(I) + DR(64)*B(I+1)
      A3=FDRL*BDRL(NROW)
      J=I
      DO 20 NCOL=NROW,NEF
      S(N)=S(N)+A1*B(J)+A2*B(J+1)+A3*BDRL(NCOL)
      N=N+1
   20 J=J+5
   30 I=I+5
      LROW=NEF
      I=1
      N1=1
      N2=3*NEF-2
      IF (IPEL.EQ.2.OR.MTMOD.EQ.2) GO TO 120
C     -----------------------
C     ELASTIC RIGIDITY MATRIX
C     -----------------------
      DO 100 IR=1,NNO
      DO 90 IRR=1,3
      A1=DR(1)*B(I)+DR(2)*B(I+1)
      A2=DR(2)*B(I)+DR(10)*B(I+1)
      A3=DR(19)*B(I+2)
      A4=DR(28)*B(I+15)+DR(29)*B(I+16)
      A5=DR(29)*B(I+15)+DR(37)*B(I+16)
      A6=DR(46)*B(I+17)
      J=I
      IF=4-IRR
C     ------------------------------------------------
C     UPPER TRIANGULAR PART OF DIAGONAL 3*3 PARTITIONS
C     ------------------------------------------------
      DO 40 JR=1,IF
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
   40 J=J+5
C     ----------------------------------------
C     INTERMEDIATE OFF-DIAGONAL 3*3 PARTITIONS
C     ----------------------------------------
      N1=N1+3
      IF (IR.EQ.NNO) GO TO 70
      NB=NNO-IR
      N2=N2+3
      J=J+15
      DO 60 JB=1,NB
      DO 50 JR=1,3
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
   50 J=J+5
      N1=N1+3
      N2=N2+3
   60 J=J+15
      N2=N2-3
   70 I=I+5
   90 LROW=LROW-1
      I=I+15
      N1=N2
      LROW=LROW-3
  100 N2=N2+3*LROW-3
      RETURN
C     --------------------------------------
C     MULTILAYERED ANISOTROPIC COMPOSITE AND
C     ELASTO-PLASTIC RIGIDITY MATRIX
C     --------------------------------------
  120 DO 200 IR=1,NNO
      DO 190 IRR=1,3
      A1=DR(1)*B(I)+DR(2)*B(I+1)+DR(3)*B(I+2)
      A2=DR(2)*B(I)+DR(10)*B(I+1)+DR(11)*B(I+2)
      A3=DR(3)*B(I)+DR(11)*B(I+1)+DR(19)*B(I+2)
      A4=DR(28)*B(I+15)+DR(29)*B(I+16)+DR(30)*B(I+17)
      A5=DR(29)*B(I+15)+DR(37)*B(I+16)+DR(38)*B(I+17)
      A6=DR(30)*B(I+15)+DR(38)*B(I+16)+DR(46)*B(I+17)
      A7=DR(4)*B(I)+DR(12)*B(I+1)+DR(20)*B(I+2)
      A8=DR(5)*B(I)+DR(13)*B(I+1)+DR(21)*B(I+2)
      A9=DR(6)*B(I)+DR(14)*B(I+1)+DR(22)*B(I+2)
      J=I
      IF=4-IRR
C     ------------------------------------------------
C     UPPER TRIANGULAR PART OF DIAGONAL 6*6 PARTITIONS
C     ------------------------------------------------
      DO 140 JR=1,IF
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
  140 J=J+5
      S(N1)=S(N1)+A7*B(J)+A8*B(J+1)+A9*B(J+2)
      S(N1+1)=S(N1+1)+A7*B(J+5)+A8*B(J+6)+A9*B(J+7)
      S(N1+2)=S(N1+2)+A7*B(J+10)+A8*B(J+11)+A9*B(J+12)
C     ---------------------------
C     OFF-DIAGONAL 6*6 PARTITIONS
C     ---------------------------
      N1=N1+3
      IF (IR.EQ.NNO) GO TO 170
      NB=NNO-IR
      N2=N2+3
      J=J+15
      A10=DR( 4)*B(I+15)+DR( 5)*B(I+16)+DR( 6)*B(I+17)
      A11=DR(12)*B(I+15)+DR(13)*B(I+16)+DR(14)*B(I+17)
      A12=DR(20)*B(I+15)+DR(21)*B(I+16)+DR(22)*B(I+17)
      DO 160 JB=1,NB
      DO 150 JR=1,3
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      S(N1+3)=S(N1+3)+A7*B(J+15)+A8*B(J+16)+A9*B(J+17)
      S(N2-3)=S(N2-3)+A10*B(J)+A11*B(J+1)+A12*B(J+2)
      N1=N1+1
      N2=N2+1
  150 J=J+5
      N1=N1+3
      N2=N2+3
  160 J=J+15
      N2=N2-3
  170 I=I+5
  190 LROW=LROW-1
      I=I+15
      N1=N2
      LROW=LROW-3
  200 N2=N2+3*LROW-3
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE MESHONE (COORD,XYZM,ICONN,IDSUPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     SUPER ELEMENT MESH
C     25 NODES, 16 ELEMENTS      
C     ------------------------------------------------------------
      DIMENSION COORD(12),XYZM(25,3),DXYZ(3) 
      DIMENSION IELEM(5),IELEMO(35),ICONS(4,16),ICONN(4,16),IDSUPER(25,6)
C     MESH MAPPING
      DATA IELEMO /1,2,10,5,9, 
     2             2,3,12,6,11, 
     3             3,4,14,7,13, 
     4             4,1,16,8,15,
     5             16,11,23,18,22, 
     6             8,6,19,17,21,
     7             15,12,24,20,25/ 
C     CONNECTIVITY IN ONE SUPER ELEMENT
      DATA ICONS /  1,9,22,16,
     2              9,5,18,22,
     3              5,10,23,18,
     4              10,2,11,23,
     5              16,22,21,8,
     6              22,18,17,21,
     7              18,23,19,17,
     8              23,11,6,19,
     9              8,21,25,15,
     1              21, 17, 20, 25, 
     2              17, 19, 24, 20,
     3              19, 6, 12, 24,
     4              15, 25, 14, 4,
     5              25, 20, 7, 14,
     6              20, 24, 13, 7,
     7              24, 12, 3, 13/
                   
      ICONN(1:4,1:16) = ICONS(1:4,1:16)
      

      !GET THE COORDINATE IN EACH NODE FOLLOWING CONNECTIVITY ORDER      
      
      J = 0
      DO I = 1,4
       J = 3*(I-1)+1
      XYZM(I,1:3) = COORD(J:J+2)
      ENDDO
      
      KK = 0
      DO K = 1,7
       KK = 5*(K-1)+1    
      IELEM(1:5) = 0 
      IELEM(1:5) = IELEMO(KK:KK+4)    
      DO II = 1,3
          DXYZ(II) = XYZM(IELEM(1),II) - XYZM(IELEM(2),II)
          XYZM(IELEM(3),II) = DXYZ(II)*0.25 + XYZM(IELEM(2),II)
          XYZM(IELEM(4),II)  = DXYZ(II)*0.50 + XYZM(IELEM(2),II)
          XYZM(IELEM(5),II)  = DXYZ(II)*0.75 + XYZM(IELEM(2),II)            
      ENDDO
      ENDDO
      
      DO INODE = 1,25
          IDSUPER(INODE,1) = 6*INODE-5
          IDSUPER(INODE,2) = 6*INODE-4
          IDSUPER(INODE,3) = 6*INODE-3
          IDSUPER(INODE,4) = 6*INODE-2
          IDSUPER(INODE,5) = 6*INODE-1
          IDSUPER(INODE,6) = 6*INODE-0
      ENDDO
      
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE MATRIX_CONDEN(SKS,SKSC,R,T,NNEQ,NEQ,ICONN,ILO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     STATIC CONDENSTAION OF MATRIX
C     SKS : STIFFNESS MATRIX TO BE REDUCED
C     R : LOAD VECTOR
C     NEQ : TOTAL NUMBER OF EQUATION
C     LEQ : NUMBER OF EQUATION TO BE ELIMINATED      
C     ------------------------------------------------------------
      CHARACTER*3 OPER
         
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C     SAVE IEL
      COMMON /INELE/ IEL   
      COMMON /JNELE/ JEL
      
C     STIFFNESS MATRIX PART
      DIMENSION SKS(NNEQ,NNEQ)
      DIMENSION SKUU(24,24),SKUC(24,126),SKCC(126,126),SKCU(126,24),
     1          SKCCI(126,126), ! INVERSION OF SKCC MATRIX
     2          SKSC(24,24),  ! CONDENSED STIFFNESS MATRIX
     3          SKCIU(24,126),SKSM(24,24), 
     4          TB(126,24),T(150,24),TI(24,24),TP(126,1)
      
      DIMENSION ICONN(4,16)

C     LOAD MATRIX PART       
      DIMENSION RELEM(NNEQ),R(NEQ),
     1          FCS(24), !CONDENSED LOAD MATRIX      
     2          FU(24),FC(126),FUC(24)
      
      DIMENSION AKK(4,4),AKK2(4,4)
      
      DIMENSION FMASS(150,150),FMASSC(24,24),FMASSCU(150,24)
      
      DATA AKK /10,5,3,2,
     1          2,6,7,5,
     2          1,5,32,7,
     3          2,4,7,5/

      
C     1. REARRANGE THE STIFFNESS MATRIX (SUPER NODES AND ELIMINATED NODES)
      
      NEQU = 24
      NEQC = 126
      
      CALL CLEARMAT(SKUU,24,24)
      CALL CLEARMAT(SKUC,24,126)
      CALL CLEARMAT(SKCC,126,126)
      CALL CLEARMAT(SKCU,126,24)
      CALL CLEARMAT(SKCCI,126,126)
      CALL CLEARMAT(SKSC,24,24)
      CALL CLEARMAT(SKCIU,24,126)
      CALL CLEARMAT(SKSM,24,24)
      CALL CLEARMAT(SKSC,24,24)
      CALL CLEARMAT(TB,126,24)

      
      SKUU(1:24,1:24)   = SKS(1:24,1:24)
      SKUC(1:24,1:126)  = SKS(1:24,25:150)
      SKCC(1:126,1:126) = SKS(25:150,25:150)
      SKCU(1:126,1:24)  = SKS(25:150,1:24)
      
      
C     2. INVERSE SKCC MATRIX
         
      CALL INVMATRIX(SKCC,SKCCI,126)
      CALL INVMATRIX(AKK,AKK2,4)

C     3. CONDENSATION OF STIFFNESS MATRIX  
C     =========================================================      
C     MULTIPLYING MATRICES USING 'DGEMM'
C     ---------------------------------------------------------
C     C = ALPHA*A*B + BETA*C   
C     A(M,K), B(K,N), C(M,N)   
C     "INITAILIZE C MATRIX BEFORE DOING CALL DGEMM"
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)     
C     =========================================================
      ALPHA = 1.0
      BETA  = 1.0
            
      M = NEQU
      N = NEQC
      K = NEQC
CC      CALL DGEMM('N','N',M,N,K,ALPHA,SKUC,M,SKCCI,K,BETA,SKCIU,M)
      
      SKCIU = MATMUL(SKUC,SKCCI)
      
      M = NEQU
      N = NEQU
      K = NEQC     
CC      CALL DGEMM('N','N',M,N,K,ALPHA,SKCIU,M,SKCU,K,BETA,SKSM,M)    
      
      SKSM = MATMUL(SKCIU,SKCU)
      
      SKSC(1:24,1:24) = SKUU(1:24,1:24) - SKSM(1:24,1:24)
      
      IF(ITASK.EQ.2.OR.ITASK.EQ.3.OR.ITASK.EQ.5) RETURN
C      !SAVE THE MATRIX
C	KREC = 5002
C      CALL MESTIF_CON(KREC,SKUU,SKUC,SKCU,SKCC,OPER)
      
      SKCCI(1:126,1:126) = -1.0D0*SKCCI(1:126,1:126)
C     TB-MATRIX          
      M = NEQC
      N = NEQU
      K = NEQC     
CC      CALL DGEMM('N','N',M,N,K,ALPHA,SKCCI,M,SKCU,K,BETA,TB,M) 
      TB = MATMUL(SKCCI,SKCU)
      

C     T-MATRIX    

      DO IT = 1,24
          DO JT = 1,24
              IF(IT.EQ.JT) THEN
                  TI(IT,JT) = 1.0
              ELSE
                  TI(IT,JT) = 0.0
              ENDIF    
          ENDDO
      ENDDO          
      T(1:24,1:24) = TI(1:24,1:24)      
      T(25:150,1:24) = TB(1:126,1:24)
      
      DO I = 1,150
          DO J = 1,150
              IF(I.EQ.J) FMASS(I,J) = 3.0D0
          ENDDO    
      ENDDO    
      
      FMASSCU = MATMUL(FMASS,T)
      FMASSC = MATMUL(TRANSPOSE(T),FMASSCU)
      
      
      IF(ILO.EQ.0) RETURN    
      
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================      
      SUBROUTINE STIFFONE(SKSC,S)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     SAVE CONDENSED STIFFNESS MATRIX IN TO S MATRIX
      DIMENSION SKSC(24,24), S(1)
      
      K = 0
      DO I = 1,24
          DO J = I,24
             K = K + 1 
             S(K) = SKSC(I,J) 
          ENDDO
      ENDDO          
      

      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================         
      SUBROUTINE SHELL4SUPER0(PROPM,PROPG,NODEX,WA,AMV,S,COORD,
	1					EDIS,EDISI,RE,MWG,FIN,MSET)
C	FIN - ADDED TO PREVIOUS LINE BY GILSON - JUL2003 (INT FORCE)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	modified by Hari December 2000
C     ----------------------------------------------------------------
C     EVALUATES THE TANGENTIAL STIFFNESS, GAUSS POINT
C     STRAINS, GAUSS POINT STRESS-RESULTANTS AND NODAL STRESS-
C     RESULTANTS FOR 4 TO 8 NODED MINDLIN TYPE ISOPAROMETRIC
C     SHELL ELEMENTS WITH 6 DOF AT EACH NODE
C	--------------------------------------
C     INPUT VARIABLES
C	---------------
C     PROPM(NMP)    = MATERIAL PROPERTIES (YM,PR,YLD,HP,DEN)
C     PROPG(NGP)    = GEOMETRIC PROPERTIES (NNO)
C     NODEX(NEX)    = LOCATIONS OF EXCESS NODES (MIDSIDE NODES)
C     WA(MWG,NPT)   = WORKING ARRAY (8 STRESSES + (8 STRAINS,YLD,IPEL))
C     COORD(3,NNO)  = CURRENT NODAL COORDINATES X,Y,Z
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     EDISI(NEF)    = CURRENT NODAL DISPLACEMENT INCREMENTS
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     S(NWS) 1176   = ELEMENT STIFFNESS MATRIX (UPPER TRIANG.ROW-WISE)
C     RE(NEF)       = EQUILIBRIUM LOADS AT ELEMENT NODES
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/
C	--------------------------------
C     NAME(2)       = NAME OF ELEMENT MODULE
C     ITYPE         = CODE NUMBER FOR ELEMENT MODULE
C     NLOPT         = CODE FOR NONLINEAR OPTION
C     NLOPT=0         LINEAR ANALYSIS
C     NLOPT=1         MATERIALLY NONLINEAR ONLY
C     NLOPT=2,3       TOTAL LAGRANGIAN,UPDATED LAGRANGIAN
C     MTMOD         = CODE FOR MATERIAL MODULE
C     MTMOD=1         LINEAR ELASTIC,ISOTROPIC
C     MTMOD=2         LINEAR ELASTIC,ORTHOTROPIC
C     MTMOD=3         ELASTO-PLASTIC (IVANOV)
C     MTMOD=4         ELASTO-PLASTIC (MULTI-LAYER)
C     MTMOD=5         CONCRETE WITH CRACKING
C     NSINC         = FACTOR CONTROLLING NUMBER OF SUBINCREMENTS
C     ITOLEY         = TOLERANCE ON YIELD FUNCTION
C     NELE          = NUMBER OF ELEMENTS IN THIS GROUP
C     NMPS          = NUMBER OF MATERIAL PROPERTY SETS
C     NGPS          = NUMBER OF GEOMETRIC PROPERTY SETS
C     NMP           = NUMBER OF MATERIAL PROPERTIES PER SET
C     NGP           = NUMBER OF GEOMETRIC PROPERTIES PER SET
C     NNM           = MAXIMUM NUMBER OF NODES FOR ANY ONE ELEMENT
C     NEX           = MAXIMUM NUMBER OF EXCESS NODES
C     NCO           = NUMBER OF NODAL COORDINATES
C     NNF           = NUMBER OF NODAL DEGREES OF FREEDOM
C     NEF           = MAXIMUM NUMBER OF ELEMENT DEGREES OF FREEDOM
C     NWG           = NUMBER OF STORAGE LOCATIONS AT EACH GAUSS POINT
C     NPT           = NUMBER OF GAUSS POINTS
C     NWA           = SIZE OF WORKING ARRAY
C     NWS           = SIZE OF ELEMENT STIFFNESS MATRIX
C     MEL           = CURRENT ELEMENT NUMBER
C     NNO           = NUMBER OF NODES FOR THIS ELEMENT
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR THIS ELEMENT
C     NELTOT        = TOTAL NUMBER OF ELEMENTS (ALL GROUPS)
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /GAUS/
C	--------------------------------
C     GLOC(10,10)   = NATURAL GAUSS POINT COORDINATES (1*1 TO 10*10)
C     GWT (10,10)   = GAUSS POINT WEIGHTS
C     NGR,NGS,NGT   = NUMBER OF GAUSS POINTS IN RN,SN,TN DIRECTION
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /FLAG/
C	--------------------------------
C     IFPRI,ISPRI   = FLAG FOR PRINTING DISPL.OR STRESSES (ISPRI=0)
C     IFPLO         = FLAG FOR PLOT OUTPUT
C     IFREF         = FLAG FOR REFORMATION OF STIFFNESS (IFREF=0)
C     IFEIG         = FLAG FOR EIGENVALUE SOLUTION (IFEIG=0)
C     ITASK = 1       FIRST ENTRY INTO ELEMENT MODULE
C     ITASK = 2       ENTRY DURING EQUILIBRIUM ITERATIONS
C     ITASK = 3       ENTRY TO WORK OUT STRESSES (LAST STEP ONLY)
C     ITASK = 4       ENTRY TO DETERMINE GEOMETRIC STIFF.MATRIX ONLY
C     KSTEP           CURRENT STEP NUMBER
C     KITE            CURRENT ITERATION NUMBER
C	---------------
C     LOCAL VARIABLES
C	---------------
C     COORDI(3,8)   = INITIAL NODAL COORDINATES
C     REDIS(48)     = COROTATIONAL FORM OF EDIS
C     DISD(12)      = DISPLACEMENT DERIVATIVES
C     EPS(8)        = GAUSS POINT STRAINS
C     EPSQ(8)       = QUADRATIC PART OF GAUSS POINT STRAINS
C     SIGR(8)       = GAUSS POINT STRESS-RESULTANTS
C     DR(64)        = ELASTO-PLASTIC RIGIDITY MATRIX STORED COLUMN-WISE
C     PR            = POISSON'S RATIO
C     TH            = GAUSS POINT THICKNESS
C     RN,SN         = NATURAL (NON-DIMENSIONAL) COORDINATES
C     SLR,SLS       = SHEAR LOCKING CONSTRAINT FACTORS ALONG RN,SN
C     NGR,NGS       = NUMBER OF GAUSS PONTS ALONG RN,SN
C     RWT(4),SWT(4) = WEIGHTING FACTORS FOR GAUSS POINTS ALONG RN,SN
C     H(8)          = SHAPE FUNCTIONS
C     HD(2,8)       = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)   = SHAPE FUNCTION DERIVATIVES W.R.T R,S
C     XJI(4)        = INVERSE JACOBIAN STORED COLUMN-WISE
C     DET           = JACOBIAN DETERMINANT
C     DVOL          = INTEGRATION FACTOR (DR*DS=DVOL*DRN*DSN)
C     VT(3)         = DIRECTION COSINE VECTOR ALONG OUTWARD NORMAL TO
C                     LOCAL RN/SN SURFACE (TRUE GEOMETRIC NORMAL)
C     VR(3),VS(3)   = DIRECTION COSINE VECTORS WITH VR TANGENTIAL TO RN
C                     AND VS NORMAL TO VR/VT PLANE (VR,VS,VT FORM A
C                     RIGHT-HANDED ORTHOGONAL SET)
C     RR,SS         = SQUARED BASE VECTOR LENGTHS
C     SNA           = SIN OF ANGLE SUBTENDED BY BASE VECTORS
C     IPEL          = SECTION PLASTICITY INDICATOR (1=EL,2=EL-PL)
C     ------------------------------------------------------------------
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /DYNA/  CDEN,IMASS
C	NEXT COMMON ADDED BY GILSON - SEPT2002
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
C	ADDED BY DE SILVA
	COMMON /HARD/  HP,DEN

C     FOR SUPER ELEMENT BY BJ
      COMMON /SUPERSTIF/ SM(24,24) 
      DIMENSION NCOLM(24)
C
      DIMENSION PROPM(5),PROPG(*),NODEX(*),WA(MWG,*),S(1176),COORD(3,8)
      DIMENSION EDIS(48),EDISI(48),RE(48),REDIS(48),COORDI(3,8)
      DIMENSION DR(64),H(8),HD(2,8),XJI(4),HR(8),HS(8),B(240)
      DIMENSION BDRL(48),DISD(12),EPS(8),EPSQ(8),SIGR(8)
      DIMENSION VR(3),VS(3),VT(3),SL(4,2),FF(6),NOD(4)
      DIMENSION COVR(3),COVS(3)
	DIMENSION BA(4,120),FJ(4),RRN(4),SSN(4)
C	NEXT LINE ADDED BY GILSON - SEPT2002
	DIMENSION AMV(3)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
C
      EQUIVALENCE (APEL,IPEL)
C
      DATA SL /2.5E5,.1334,.1334,.1334,1.E20,19.2,4.00,4.00/
C
      PI=3.1415926535898
      CALL CLEARA (BDRL,48)
C     ---------------------------------
C     OBTAIN LINEAR STRESS - STRAIN LAW
C     mtmode=1, ipel=1 for linear analysis
	IPEL=1
      If (MTMOD.NE.2) CALL HOKLAW (PROPM,PROPG,2)
C	here 2 indicates plane stress parameters

C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C     ------------------------------------------------------
      CALL SHINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)
C	this subroutine updated coordinate adding displacement for nonlinear 
C	analysis Loop over sampling points, four points for transverse strains
      RRN(1)=1.0
	SSN(1)=0.0
	RRN(2)=-1.0
	SSN(2)=0.0
	RRN(3)=0.0
	SSN(3)=1.0
	RRN(4)=0.0
	SSN(4)=-1.0
C	make all assumed strain Bs values zero at the start
      DO 5 I=1,4
	FJ(I)=0.00
	XJI(I)=0.0
	DO 5 J=1,120
 5	BA(I,J)=0.00
 	NSP=4
      DO 6 I=1,NSP
	CALL SHAP2D (RRN(I),SSN(I),H,HD,NODEX,NNO) 
C	nodex IS NOT USED here for four node element, 
      CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FJ)
C	unit vectors and jacobians are calculated at sampling points
	CALL SHBMATS (NNO,H,HD,VR,VS,VT,XJI,HR,HS,BA,I,FJ)
 6    CONTINUE
C	BA matrix at 4 nodes are completed
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
      IPT=0
      DO 900  IGR=1,NGR
      RN = GLOC(IGR,NGR)
      DO 900  IGS=1,NGS
      SN = GLOC(IGS,NGS)
      WT = GWT(IGR,NGR)*GWT(IGS,NGS)
      IPT= IPT+1

C     -----------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C     -----------------------------------------------------
      CALL SHAP2D (RN,SN,H,HD,NODEX,NNO)
C     ----------------------------------------------------
C     INVERSE JACOBIAN (XJI) , JACOBIAN coefficient F, 
C	DETERMINANT (DET) AND STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------------------
      CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FJ)
      DVOL=WT*DET
	
	CALL SHSTORL1 (COORD,PROPG,WA(1,IPT),NWG,VR,VS,VT)	!STORE SHELL LAX FOR STRESS CALCULATION JAN09

      CALL RELFILL('TWGH',DVOL*TH*PROPM(5),1,KEG,2)
      
	GOTO 120      
C      IF (NGP.LE.2) GOTO 120
C     ----------------------------------------
C     INTERPOLATE NODAL THICKNESSES (NGP.GT.2)
C     ----------------------------------------
      TH=0.
      DO 100 I=1,NNO
  100 TH=TH+H(I)*PROPG(I+1)   
  
	GOTO 121
	  
  120 TH=PROPG(2)

121	IF (ITASK.NE.5) GOTO 140

C     ---------------------
C     MASS MATRIX (ITASK=5)
C     ---------------------
	IF (MTMOD.EQ.2) PROPM(5) = CDEN
C	CONMSS ADDED BY SONGSAK FOR RC SHELL
	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6) THEN
	CALL CONMSS(TH,MSET,RHORC)
	CALL SHMASS (S,H,VR,VS,DVOL,RHORC,TH,IMASS,NNO,NEF)
	ELSE
      CALL SHMASS (S,H,VR,VS,DVOL,PROPM(5),TH,IMASS,NNO,NEF)
	ENDIF
      GOTO 900


  140 CALL SHBMAT1 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,RN,SN,BA,NEF)
      CALL SHBMAT2 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B)
C     -------------------------------------------
C     SHEAR LOCKING CONSTRAINT FACTORS (SLR,SLS), 
C	MODIFIED WITH ASSUME STRAIN
C     -------------------------------------------
      SLR=1.
	SLS=1.
C	--------------------------------
C	DETERMINE INITIAL MATERIAL ANGLE
C	added by gislon - sept2002
C	--------------------------------
      SELECTCASE(MTMOD)
      CASE(2,5,6)
	CALL CANGLE(VR,VS,AMV,ANG)
	ENDSELECT

C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES (IFREF.EQ.0)
C     -----------------------------------------
      IF (MTMOD.NE.2) CALL SHDELA (DR,DVOL,SLR,SLS)
      IF (NLOPT+ITASK.EQ.1) GOTO 700
C     -------------------------------------------------
C     REMOVE RIGID BODY TRANSLATIONS AND ROTATIONS FROM
C     TOTAL DISPLACEMENT VECTOR (NLOPT>1)
C     -------------------------------------------------
      IF (NLOPT.LE.1) GOTO 180
      CALL SHMDSP (COORD,COORDI,EDIS,REDIS,H,HD,VR,VS,VT,NNO)
C     ------------------------
C     DISPLACEMENT DERIVATIVES
C     ------------------------
  180 CALL CLEARA (DISD,12)
      K=1
      M=1
      DO 200 I=1,NNO
      DO 190 J=1,3
      L=K+3
      DISD(1)=DISD(1)+B(M)*REDIS(K)	
      DISD(2)=DISD(2)+B(M+1)*REDIS(K)
      DISD(3)=DISD(3)+(B(M+2)-HR(I)*VS(J))*REDIS(K)
      DISD(4)=DISD(4)+(B(M+2)-HS(I)*VR(J))*REDIS(K)
      DISD(5)=DISD(5)+B(M+15)*REDIS(L)
      DISD(6)=DISD(6)+B(M+16)*REDIS(L)
      DISD(7)=DISD(7)+B(M+17)*REDIS(L)
      DISD(9)=DISD(9)+B(M+3)*REDIS(K)
      DISD(10)=DISD(10)+B(M+18)*REDIS(L)
      DISD(11)=DISD(11)+B(M+4)*REDIS(K)
      DISD(12)=DISD(12)+B(M+19)*REDIS(L)
      K=K+1
 190  M=M+5
      K=K+3
 200  M=M+15
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      EPS(1) = DISD(1)
      EPS(2) = DISD(2)
      EPS(3) = DISD(3)+DISD(4)
      EPS(4) = DISD(5)
      EPS(5) = DISD(6)
      EPS(6) = DISD(7)
      EPS(7) = DISD(9) +DISD(10)
      EPS(8) = DISD(11)+DISD(12)
C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400
      EPSQ(1)=0.5*(DISD(1)*DISD(1)+DISD(4)*DISD(4)+DISD(9)*DISD(9))
      EPSQ(2)=0.5*(DISD(2)*DISD(2)+DISD(3)*DISD(3)+DISD(11)*DISD(11))
      EPSQ(3)=     DISD(3)*DISD(1)+DISD(2)*DISD(4)+DISD(11)*DISD(9)
      EPSQ(7)=     DISD(10)*DISD(1)+DISD(12)*DISD(4)
      EPSQ(8)=     DISD(10)*DISD(3)+DISD(12)*DISD(2)
      DO 300  I=1,3
 300  EPS(I)= EPS(I)-EPSQ(I)
      EPS(7)= EPS(7)-EPSQ(7)
      EPS(8)= EPS(8)-EPSQ(8)
 400  EPS(7)= EPS(7)
      EPS(8)= EPS(8)
C     ---------------------------------------
C     GAUSS POINT STRESS - RESULTANTS SIGR(8)
C     ---------------------------------------
      GO TO (410,465,450,460,467,468), MTMOD
 410  CALL MPSIGA (EPS,SIGR)

	CALL SHSTRS(VR,VS,VT,EPS,WA(9,IPT),PROPM,TH)

       DO 420 I=1,8
 420  WA(I,IPT)=SIGR(I)
      GO TO 500

c-----COMPOSITE 
 465  CALL COMRGD(EPS,SIGR,DR,PROPM,TH,ANG,DVOL,'SHELL')
      DO 426 I=1,8
	WA(I  ,IPT) =  EPS(I)
 426	WA(I+8,IPT) = SIGR(I)
      GO TO 500

 450	IF (HP.EQ.0.) THEN
	CALL IVANOV (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,SIGR,DR,DVOL)
      ELSE
	EPCA = 0.
	CALL IVANOVH (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,
	1					SIGR,DR,DVOL,WA(18,IPT),EPCA)
	END IF
	  APEL = WA(17,IPT)
	  IPEL = 2
      GOTO 500

 460	IF (HP.EQ.0.) THEN
	CALL MLAYER (WA(1,IPT),EPS,SIGR,DR,DVOL)
      IPEL=2
      DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GOTO 500
	ELSE
	CALL MLAYERH (WA(1,IPT),EPS,SIGR,DR,DVOL)
      IPEL=2
      DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GO TO 500
	END IF

 467	CALL RCLAYR (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GO TO 500

 468	CALL RCLAYRE (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510 I=1,8
 510  SIGR(I)=SIGR(I)*DVOL

C	NEXT LINE CHANGED BY GILSON - JUL2003 (INT FORCE)
C      IF (ITASK.LE.2) GOTO 520
      IF (ITASK.LE.3.AND.ISOLOP.NE.4) GOTO 520 

c      CALL SECOND(T1,TIM1)
      IF (IFEIG.EQ.0) GOTO 800

      GOTO 900
 520  K=1
      M=1
      DO 560 I=1,NNO
      DO 540 J=1,3
      L=K+3
      RE(K)=RE(K)+B(M)*SIGR(1)+B(M+1)*SIGR(2)+B(M+2)*SIGR(3)
     1                      +B(M+3)*SIGR(7)+B(M+4)*SIGR(8)
      RE(L)=RE(L)+B(M+15)*SIGR(4)+B(M+16)*SIGR(5)+B(M+17)*SIGR(6)
     1                         +B(M+18)*SIGR(7)+B(M+19)*SIGR(8)
      K=K+1
 540  M=M+5
      K=K+3
 560  M=M+15

C	---------------------------------------
C	FORCES FROM DRILLING DOF SONGSAK JAN09
      IF(NLOPT.EQ.0) THEN
      CALL SHBDRL (NNO,H,HR,HS,VR,VS,VT,BDRL)
	FDRL=10.*DR(19)
	DO IEF = 1,NEF
	DO JEF = 1,NEF
	RE(IEF) = RE(IEF) + BDRL(IEF)*FDRL*BDRL(JEF)*REDIS(JEF)
	ENDDO
	ENDDO
	ENDIF
C	---------------------------------------

C     ------------------------------------------------------
C     FIND LINEAR CONTRIBUTION TO STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------------------
      IF (IFREF ) 900,700,900
c 700  CALL SECOND(T1,TIM1)
C     -----------------------------------------------------------
C     REMOVE SINGULARITY FROM STIFFNESS OF EXACTLY PLANE ELEMENTS
C     -----------------------------------------------------------
 700     IJ=1
      N=3*NEF-2
      LROW=NEF-3
      FAC=DMIN1(DR(19),DR(28)/DET)*1.E-6
      DO 720 I=1,3
      FI=FAC*VT(I)
      DO 720 J=I,3
      FF(IJ)=FI*VT(J)
 720  IJ=IJ+1
      DO 740 I=1,NNO
      HH=H(I)*H(I)
      S(N)=S(N)+HH*FF(1)
      S(N+1)=S(N+1)+HH*FF(2)
      S(N+2)=S(N+2)+HH*FF(3)
      N=N+LROW
      S(N)=S(N)+HH*FF(4)
      S(N+1)=S(N+1)+HH*FF(5)
      N=N+LROW-1
      S(N)=S(N)+HH*FF(6)
      N=N+4*LROW-14
 740  LROW=LROW-6
      CALL SHBDRL (NNO,H,HR,HS,VR,VS,VT,BDRL)
      CALL SHKLIN (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD)
C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 810
C 800  IF (ISTYP.EQ.1) GOTO 805
 800  CALL SHGEO1 (S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)
C      GOTO 810
C 805  CALL SHGEO2 (S,SIGR,HR,HS,NNO)
c 810  CALL SECOND(T1,TIM2)
 810  CONTINUE !TIM(12)=TIM(12)+TIM2
 900  CONTINUE


C	KEF = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	ENDDO
C	ENDDO
C	MEF = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	MEF = MEF + 1
C	S(KEF) = S(MEF) 
C	ENDDO
C	ENDDO



C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
      ENDIF

C     SAVE THE STIFFNESS MATRIX INTO FULL MATRIX FORM
      II = 25
      DO I = 1,24  
          II = II-1
          NCOLM(I) = II
      ENDDO

      KK = 0
      IJ = 0
      DO II = 1,24
          IJ = IJ + 1
         DO JJ = II,24
             KK = KK + 1
             SM(II,JJ) = S(KK)
         ENDDO    
      ENDDO      
      
      DO II = 1,24
         DO JJ = 1,24
          IF(II.NE.JJ) SM(JJ,II) = SM(II,JJ)
         ENDDO
      ENDDO        
      
      RETURN
      END
C
C======================================================================   
C======================================================================
C======================================================================      
      
C======================================================================
C	START 4 NODE SHELL ELEMENT ROUTINES - MAIN SUBROUTINE
C======================================================================
      SUBROUTINE XSHELL24SUPER(PROPM,PROPG,NODEX,WA,AMV,S,COORD,
	1					EDIS,EDISI,RE,MWG,FIN,MSET)
C	FIN - ADDED TO PREVIOUS LINE BY GILSON - JUL2003 (INT FORCE)
C	ANG - ADDED TO PREVIOUS LINE - SEPT2002
C	=================================================================
C	FOUR NODE SHELL ELEMENT WITH 24 DEGREES OF FREEDOM
C	THE DRILLING DEGREE IS BASED ON A CONTINUUM MECHANICS DEFINITION
C	THE GEOMETRIC STIFFNESS IS BASED ON A FULL GREEN STRAIN EXPANSION
C	=================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	integer apel
C
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/ A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
      COMMON /DYNA/ CDEN,IMASS
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
C	ADDED BY DE SILVA
	COMMON /HARD/  HP,DEN
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL      

C     FOR SUPER ELEMENT BY BJ
      COMMON /SUPERSTIF/ SM(24,24) 
      DIMENSION NCOLM(24)      
C
      DIMENSION PROPM(5),PROPG(*),NODEX(*),WA(MWG,*),S(*),COORD(3,8)
      DIMENSION EDIS(48),EDISI(48),RE(48),REDIS(48),COORDI(3,8)
      DIMENSION DR(64),H(8),HD(2,8),XJI(4),HR(8),HS(8)!,B(240)
      DIMENSION DISD(12),EPS(8),EPSQ(8),SIGR(8),FF1(4)
      DIMENSION VR(3),VS(3),VT(3),SL(4,2)!,NOD(4),EC(3),FF(6)
	DIMENSION TEDIS(48),RETMP(24),ST(300)

C	COMMON VARIABLES
	DIMENSION RS(2,4),TRS(4),SLN(4),RN(4),SN(4),BL(4)
	DIMENSION FLR(4),FLS(4),FLSS(4),FLRR(4),ANT(4)
	DIMENSION DWR(12),DWRS(12),DWS(12),DWSR(12),AKBS(300)
	DIMENSION TT(24,24)

C	BENDING VARIABLES
	DIMENSION AB1(4,4),AB2(3,3),CB(11,24),ACB(11,24)!,PDPB(11,11)
	DIMENSION APB(11)!,TEMB1(11,24),TEMB2(24,11),AKB(24,24)
C	DIMENSION ABB2(11,11),ABBI2(11,11),ACB2(11,24)
	DIMENSION PBPM(11,5),PBPB(11,11)

C	MEMBRANE VARIABLES
	DIMENSION CM(5,24),ACM(5,24)!,PSP(5,5),ALM(2,24)
	DIMENSION APM(6)
C	DIMENSION AMM2(5,5),TEMM1(5,24),ACM2(5,24),TEMM2(24,5)
C	DIMENSION AMMI2(5,5),AKM(24,24)
	DIMENSION PMPM(5,5),PMPB(5,11)

C	TRANSVERSE SHEAR VARIABLES
	DIMENSION CS(2,24),CSC(4,24),ACS(2,24),ACSC(4,24)
	DIMENSION APS(4),PTP(2,2)
C	DIMENSION ASS2(2,2),ASSI2(2,2),ACS2(2,24)
C	DIMENSION TEMS1(2,24),TEMS2(24,2),AKS(24,24)

C	FOR GEOMETRIC STIFFNESS
	DIMENSION CG(15,24),ACG(15,24),PFP(15,15)

C	INTERNAL FORCE VECTOR
	DIMENSION FM(24),FB(24),FS(24),FD(24),APSH(2),APMT(5)
	DIMENSION FMB(24),FBM(24)
	DIMENSION PAPM(5),PAPB(11),PAPSH(2)

C	FOR DRILLING STIFFNESS
	DIMENSION CD(24),ACD(24)!,ACT(24)

C	FOR PLASTICITY - ADDED BY GILSON - NOV2002
	DIMENSION ANTR(4),ANTS(4),ANTRS(4)
	DIMENSION ANTRR(4),ANTSS(4),ANTRRS(4),ANTSSR(4),ANTRRSS(4)

C	FOR CHECKING - TEMPORARY VARIABLES
C	DIMENSION STIFF(24,24)
C
C	NEXT LINE ADDED BY GILSON - SEPT2002
	DIMENSION AMV(3)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
C
      EQUIVALENCE (APEL,IPEL)
C
      DATA SL /2.5E5,.1334,.1334,.1334,1.E20,19.2,4.00,4.00/
C
      PI=3.1415926535898
	DVOL=1.0
	SLR=1.0
	SLS=1.0
	KST = 0
C     -----------------------------------
C     OBTAIN LINEAR STRESS - STRAIN LAW
C     MTMOD=1, IPEL=1 FOR LINEAR ANALYSIS
C	HOKLAW - FLAG=2 FOR PLANE STRESS
C	-----------------------------------
	IPEL=1
      IF (MTMOD.NE.2) CALL HOKLAW (PROPM,PROPG,2)

C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C     ------------------------------------------------------
      CALL SHINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)

C	-----------------------------------------------
C	DETERMINE (1)LOCAL CORDINATE VECTORS - VR,VS,VT
C	          (2)LOCAL NODAL COORDINATES - RS
C	-----------------------------------------------
	CALL VRST(COORD,VR,VS,VT,RS,TRS,NNO)
	IF (ITASK.NE.5) GOTO 1333

C	----------------
C	FORM MASS MATRIX
C	----------------
	IF (MTMOD.EQ.2) PROPM(5) = CDEN


	IF (IMASS.EQ.1) THEN
C	  LUMPED MASS MATRIX
C	  ------------------
	  AREA = (-RS(2,4)/2+RS(2,2)/2)*RS(1,1)+(-RS(2,1)/2+RS(2,3)/2)*R
     #         S(1,2)+(RS(2,4)/2-RS(2,2)/2)*RS(1,3)+(RS(2,1)/2-
     #         RS(2,3)/2)*RS(1,4)
        CALL SHMASS (S,H,VR,VS,AREA,PROPM(5),TH,IMASS,NNO,NEF)
	ELSE
C	  CONSISTENT MASS MATRIX
C	  ----------------------
C       LOOP OVER GAUSS POINTS
        DO 911  IGR=1,NGR
          RN = GLOC(IGR,NGR)
          DO 911  IGS=1,NGS
            SN = GLOC(IGS,NGS)
            WT = GWT(IGR,NGR)*GWT(IGS,NGS)
C           SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
            CALL SHAP2D (RN,SN,H,HD,NODEX,NNO)
C           INVERSE JACOBIAN (XJI) , JACOBIAN DETERMINANT (DET)
C           AND STRAIN-DISPLACEMENT MATRIX (B)
            CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FF1)
            DVOL=WT*DET
            TH = PROPG(2)

	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6) THEN
	CALL CONMSS(TH,MSET,RHORC)
	CALL SHMASS (S,H,VR,VS,DVOL,RHORC,TH,IMASS,NNO,NEF)
	ELSE
      CALL SHMASS (S,H,VR,VS,DVOL,PROPM(5),TH,IMASS,NNO,NEF)
	ENDIF

911	  CONTINUE
	ENDIF
	RETURN

C	---------------------------------------
C	FORM TRANSFORMATION MATRIX TRANSPOSE(T)
C	---------------------------------------
1333	DO 30 I = 1,3
	  TT(1,I)=VR(I)
	  TT(2,I)=VS(I)
	  TT(3,I)=VT(I)
30	CONTINUE

C	----------------------------------------
C	ADD RIGID-LINK TERMS FOR WARPED GEOMETRY
C	----------------------------------------
C	NEXT IF BLOCK ADDED BY GILSON - APRIL2003
	IF (NLOPT.GE.2) THEN
	  WFAC = 0.0         !SUPPRESS FOR NONLINEAR
	ELSE
	  WFAC = 1.0
	ENDIF

	IROW = 1
	ICOL = 4
	DO 200 I=1,4
	  DO 220 K=1,3
C		'WFAC' ADDED TO NEXT TWO LINES BY GILSON - APRIL2003
C	    TT(IROW,  ICOL+K-1)=-TRS(I)*VS(K)
C	    TT(IROW+1,ICOL+K-1)=+TRS(I)*VR(K)
	    TT(IROW,  ICOL+K-1)=-TRS(I)*VS(K)*WFAC
	    TT(IROW+1,ICOL+K-1)=+TRS(I)*VR(K)*WFAC
220	  CONTINUE
	  IROW = IROW+6
	  ICOL = ICOL+6
200	CONTINUE

C	---------------------------------------------
C	COMPUTE FOR ELEMENT BOUNDARY 
C		(1)LENGTHS  - SLN
C	    (2)NORMALS  - RN
C	    (3)TANGENTS - SN
C	    (4)SHEAR DEF PAR, LAMBDA - BL, BLR, BLS
C	---------------------------------------------
	CALL LNT(RS,PROPG(2),PROPM(2),SLN,RN,SN,BL,BLR,BLS)

C	--------------------------
C	COMPUTE FOR COMMON FACTORS
C	--------------------------
	CALL COMFACT(RS,BL,BLR,BLS,SLN,RN,SN,FLR,FLS,FLSS,FLRR,
	1                   DWR,DWRS,DWS,DWSR,ANT)
	IF (MTMOD.GT.2.AND.(NLOPT+ITASK).NE.1) THEN
	  CALL COMFACT2(RS,ANTR,ANTS,ANTRS)
	  IF (IFREF.EQ.0) THEN
	    CALL COMFACT3(RS,ANTRR,ANTSS,ANTRRS,ANTSSR,ANTRRSS)
	  ENDIF
	ENDIF

C	-------------------------------------------------------------
C	COMPUTE FOR THE FLEXURAL CONTRIBUTION TO THE LINEAR STIFFNESS
C	-------------------------------------------------------------
C	Cb FOR BENDING STRAINS - CB
	CALL CBEND(RS,SLN,RN,SN,BL,BLR,BLS,FLR,FLS,FLRR,FLSS,
	1                 DWR,DWRS,DWS,DWSR,CB)
C	Ab FOR BENDING STRAINS - AB1,AB2
	CALL ABEND(RS,SLN,RN,SN,AB1,AB2)
C	MULTIPLY (Ab-1)(Cb)  - ACB
	CALL ABCB(AB1,AB2,CB,ACB)

      CB(1:11,1:24) = 0.0D0
      ACB(1:11,1:24) = 0.0D0
      
      
C	----------------------------------------------------------
C	COMP FOR THE MEMBRANE CONTRIBUTION TO THE LINEAR STIFFNESS
C	----------------------------------------------------------
C	Cm FOR MEMBRANE STRAINS - CM
	CALL CMEM20(RS,SLN,RN,SN,CM)
C	MULTIPLY (Am-1)(Cm)  - ACM
	CALL ABCM20(AB1,CM,ACM)

C	------------------------------------------------------------------
C	COMP FOR THE TRANSVERSE SHEAR CONTRIBUTION TO THE LINEAR STIFFNESS
C	------------------------------------------------------------------
C	Cs FOR SHEAR STRAINS - CS
	CALL CSHR(SLN,RN,SN,FLR,FLS,DWR,DWS,CSC,BLR,BLS)
C	MULTIPLY (As-1)(Cs)  - ACS
	CALL ABCS(AB1,CSC,VR,VS,VT,ACS,ACSC)
      
      CS(1:2,1:24) = 0.0D0
      ACS(1:2,1:24) = 0.0D0

C	-----------------------------------------------------------
C	COMP FOR THE TORSIONAL CONTRIBUTION TO THE LINEAR STIFFNESS
C	-----------------------------------------------------------
C	Cd FOR DRILLING STRAINS - CD
	CALL CDRILL(ANT,CM,CD)
C	MULTIPLY (Ad-1)(Cd)  - ACD
	DO 400 IK=1,19,6
	  ACD(IK)   = CD(IK  )/AB1(1,1)
	  ACD(IK+1) = CD(IK+1)/AB1(1,1)
400   ACD(IK+5) = CD(IK+5)/AB1(1,1)
      
      CD(1:24) = 0.0D0
      ACD(1:24) = 0.0D0
C     -------------------------------------------------
C     REMOVE RIGID BODY TRANSLATIONS AND ROTATIONS FROM
C     TOTAL DISPLACEMENT VECTOR (NLOPT>1)
C     -------------------------------------------------
      IF (NLOPT.LE.1) GOTO 1800
      CALL SHMDSP2(COORD,COORDI,EDIS,VR,VS,VT,REDIS,NNO)

C	--------------------------------
C	DETERMINE STRAIN PARAMETER ALPHA
C	--------------------------------
1800	IF (NLOPT+ITASK.NE.1) THEN
C	  TRANSFORM REDIS TO LOCAL - TEDIS
	  CALL TREDIS(TT,REDIS,TEDIS)
C	  MEMBRANE
	  CALL APHAM20(ACM,TEDIS,APM)
C	  BENDING
	  CALL APHAB(ACB,TEDIS,APB)
C	  TRANSVERSE SHEAR
	  CALL APHAS(ACSC,TEDIS,APS)
	ENDIF

C	----------------------------------------------------------
C	LOOP OVER NODES TO DET STRESS CONT TO ELEMENT FORCE VECTOR
C	----------------------------------------------------------
	IPT=0
	IPOINT = 4
	DO 1000 II=1,IPOINT
	IPT=IPT+1
C	--------------------------------
C	DETERMINE INITIAL MATERIAL ANGLE
C	added by gislon - sept2002
C	--------------------------------
      SELECTCASE(MTMOD)
      CASE(2,5,6)
	CALL CANGLE(VR,VS,AMV,ANG)
	ENDSELECT

	
	CALL SHSTORL1 (COORD,PROPG,WA(1,IPT),NWG,VR,VS,VT)	!STORE SHELL LAX FOR STRESS CALCULATION JAN09

      
C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES (IFREF.EQ.0)
C	DVOL=1, SLA=1, SLB=1
C     -----------------------------------------
	TH=PROPG(2)

C	IF(NGP.GT.2) THEN
C	TH = PROPG(IPT+1)
C	ENDIF

C	WRITE(*,*) MEL,PROPG(1),PROPG(2),PROPG(3),PROPG(4),PROPG(5),TH
C	PAUSE

C      IF(IFLOOR.EQ.3.AND.ITASK.EQ.3.AND.NEG.GT.1) THEN
         IF (MTMOD.NE.2) CALL SHDELA_WALL (DR,DVOL,SLR,SLS)
C      ELSEIF(IFLOOR.NE.3.OR.NEG.EQ.1) THEN
C          IF (MTMOD.NE.2) CALL SHDELA (DR,DVOL,SLR,SLS)    
C      ENDIF    
      
      
      IF (NLOPT+ITASK.EQ.1) GOTO 6767

C	------------------------
C     DISPLACEMENT DERIVATIVES
C	------------------------
	CALL CLEARA (DISD,12)
C	COMPONENTS OF MEMBRANE STRAINS
	DISD(1) =APM(1)+RS(2,II)*APM(2)
	DISD(2) =APM(3)+RS(1,II)*APM(4)
	DISD(3) =APM(5)
	DISD(4) =APM(6)
C	COMPONENTS OF BENDING STRAINS
	DISD(5) =APB(1)+APB(2)*RS(1,II)+APB(3)*RS(2,II)+
	1         APB(4)*RS(1,II)*RS(2,II)
	DISD(6) =APB(5)+APB(6)*RS(1,II)+APB(7)*RS(2,II)+
	1         APB(8)*RS(1,II)*RS(2,II)
	DISD(7) =APB(9)+APB(10)*RS(1,II)+APB(11)*RS(2,II)
C	COMPONENTS OF TRANSVERSE SHEAR STRAINS
	DISD(9) =APS(1)
	DISD(10)=APS(2)
	DISD(11)=APS(3)
	DISD(12)=APS(4)

C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      EPS(1) = DISD(1)
      EPS(2) = DISD(2)
      EPS(3) = DISD(3)+DISD(4)
      EPS(4) = DISD(5)
      EPS(5) = DISD(6)
      EPS(6) = DISD(7)
      EPS(7) = DISD(9) +DISD(10)
      EPS(8) = DISD(11)+DISD(12)

C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 4000
      EPSQ(1)=0.5*(DISD(1)*DISD(1)+DISD(4)*DISD(4)+DISD(9)*DISD(9))
      EPSQ(2)=0.5*(DISD(2)*DISD(2)+DISD(3)*DISD(3)+DISD(11)*DISD(11))
      EPSQ(3)=     DISD(3)*DISD(1)+DISD(2)*DISD(4)+DISD(11)*DISD(9)
      EPSQ(7)=     DISD(10)*DISD(1)+DISD(12)*DISD(4)
      EPSQ(8)=     DISD(10)*DISD(3)+DISD(12)*DISD(2)
      DO 3000  I=1,3
3000	EPS(I)= EPS(I)-EPSQ(I)
      EPS(7)= EPS(7)-EPSQ(7)
      EPS(8)= EPS(8)-EPSQ(8)

C     ---------------------------------
C     NODAL STRESS - RESULTANTS SIGR(8)
C     ---------------------------------
4000	GOTO (4100,4650,4500,4600,4605,4606), MTMOD

C	--------------
C	LINEAR ELASTIC
C	--------------
4100	CALL MPSIGA (EPS,SIGR)
	
	CALL SHSTRS(VR,VS,VT,EPS,WA(9,IPT),PROPM,TH)

4110	  DO 4200 I=1,8
4200	  WA(I,IPT)=SIGR(I)
        GOTO 1000

C	--------------
C	COMPOSITE
C	---------
4650	CALL COMRGD(EPS,SIGR,DR,PROPM,TH,ANG,DVOL,'SHELL')
	DO 4654 I=1,8
		WA(I  ,IPT)= EPS(I)
4654		WA(I+8,IPT)=SIGR(I)
      GOTO 1000

C	--------------
C	ELASTO-PLASTIC
C	--------------

4500	IF (HP.EQ.0.) THEN
	CALL IVANOV (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,SIGR,DR,DVOL)
      ELSE
	EPCA = 0.
	CALL IVANOVH (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,
	1			  SIGR,DR,DVOL,WA(18,IPT),EPCA)
	END IF
	APEL = WA(17,IPT)
	IPEL = 2
	GOTO 4700

4600	IF (HP.EQ.0.) THEN
	CALL MLAYER  (WA(1,IPT),EPS,SIGR,DR,DVOL)
	ELSE 
	CALL MLAYERH (WA(1,IPT),EPS,SIGR,DR,DVOL)
	END IF
      IPEL=2
	GOTO 4700

4605	CALL RCLAYR (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	GO TO 4700

4606	CALL RCLAYRE (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	GO TO 4700

4700	CALL PSIGR(IPT,SIGR,ANT(II),ANTR(II),ANTS(II),ANTRS(II),
	1           PAPM,PAPB,PAPSH)
	CALL PDRP(IPT,DR,ANT(II),ANTR(II),ANTS(II),ANTRS(II),ANTRR(II),
	1          ANTSS(II),ANTRRS(II),ANTSSR(II),ANTRRSS(II),
	2          PMPM,PMPB,PBPB,PTP)

1000	CONTINUE

C	NEXT LINE CHANGED BY GILSON - JUL2003 (INT FORCE)
C	IF (ITASK.LE.2) GOTO 6767
	IF (ITASK.LE.3.AND.ISOLOP.NE.4) GOTO 6767
      IF (IFEIG.EQ.0) GOTO 8000

C	---------------------------------
C	DETEREMINE TRANSPOSE(P)*MODULUS*P
C	---------------------------------
6767	IF (IPEL.EQ.1.AND.MTMOD.NE.2) THEN

C	ELASTIC CASE
C	  MEMBRANE - Nm*S*Pm
	  CALL NMSPM20(AB1,DR,PMPM)
C	  BENDING - Nb*D*Pb
	  CALL NBDPB(AB1,DR,PBPB)
C	  TRANSVERSE SHEAR - Ns*T*Ps
	  PTP(1,1) = AB1(1,1)*DR(55)
	  PTP(2,2) = PTP(1,1)
	ELSE IF (MTMOD.EQ.2) THEN
C	COMPOSITE CASE
C	  MEMBRANE-MEMBRANE - Nm*A*Pm
	  CALL SPMPM2(AB1,DR,PMPM)
C	  MEMBRANE-BENDING  - Nm*B*Pb
	  CALL SPMPB2(AB1,DR,PMPB)
C	  BENDING-MEMBRANE  - Nb*B*Pm
	  PBPM = TRANSPOSE(PMPB)
C	  BENDING-BENDING   - Nb*D*Pb
	  CALL SPBPB(AB1,DR,RS,PBPB)
C	  TRANSVERSE SHEAR  - Ns*G*Ps
	  PTP(1,1) = AB1(1,1)*DR(55)
	  PTP(1,2) = AB1(1,1)*DR(56)
	  PTP(2,1) = AB1(1,1)*DR(63)
	  PTP(2,2) = AB1(1,1)*DR(64)
	ELSE
C	INELASTIC CASE - FILL IN LOWER-TRIANGULAR PART
C	  MEMBRANE-MEMBRANE - Nm*A*Pm
	  DO 700 I = 1,5
	   DO 705 J = I,5
	    PMPM(J,I) = PMPM(I,J)
705	   CONTINUE
700	  CONTINUE
C	  MEMBRANE-BENDING  - Nm*B*Pb
C		Determined in the Stress Resultants Routine
C	  BENDING-MEMBRANE  - Nb*B*Pm
	  PBPM = TRANSPOSE(PMPB)
C	  BENDING-BENDING   - Nb*D*Pb
	  DO 710 I = 1,11
	   DO 715 J = I,11
	    PBPB(J,I) = PBPB(I,J)
715	   CONTINUE
710	  CONTINUE
C	  TRANSVERSE SHEAR  - Ns*G*Ps
	  PTP(2,1) = PTP(1,2)
	ENDIF

C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
      IF (NLOPT+ITASK.EQ.1) GOTO 7000

C	MEMBRANE CONTRIBUTION
C	---------------------
	IF (IPEL.EQ.1.OR.MTMOD.EQ.2) THEN
	  APMT(1) = APM(1)
	  APMT(2) = APM(2)
	  APMT(3) = APM(3)
	  APMT(4) = APM(4)
	  APMT(5) = APM(5)+APM(6)
	  PAPM  = MATMUL(PMPM,APMT)
	ENDIF
	DO 7100 IR = 1,19,6
	  IS = IR+1
	  FM(IR) = ACM(1,IR)*PAPM(1)+ACM(2,IR)*PAPM(2)+ACM(5,IR)*PAPM(5)
	  FM(IS) = ACM(3,IS)*PAPM(3)+ACM(4,IS)*PAPM(4)+ACM(5,IS)*PAPM(5)
7100	CONTINUE
C	BENDING CONTRIBUTION
C	--------------------
	IF (IPEL.EQ.1.OR.MTMOD.EQ.2) PAPB = MATMUL(PBPB,APB)
	DO 7200 IR = 3,21,6
	  IS = IR+1
	  IT = IR+2
	  FB(IR) = ACB(1,IR)*PAPB(1)+ACB(2,IR)*PAPB(2)+ACB(3,IR)*PAPB(3)
     #          +ACB(4,IR)*PAPB(4)+ACB(5,IR)*PAPB(5)+ACB(6,IR)*PAPB(6)
     #          +ACB(7,IR)*PAPB(7)+ACB(8,IR)*PAPB(8)+ACB(9,IR)*PAPB(9)
     #          +ACB(10,IR)*PAPB(10)+ACB(11,IR)*PAPB(11)
	  FB(IS) = ACB(1,IS)*PAPB(1)+ACB(2,IS)*PAPB(2)+ACB(3,IS)*PAPB(3)
     #          +ACB(4,IS)*PAPB(4)+ACB(5,IS)*PAPB(5)+ACB(6,IS)*PAPB(6)
     #          +ACB(7,IS)*PAPB(7)+ACB(8,IS)*PAPB(8)+ACB(9,IS)*PAPB(9)
     #          +ACB(10,IS)*PAPB(10)+ACB(11,IS)*PAPB(11)
	  FB(IT) = ACB(1,IT)*PAPB(1)+ACB(2,IT)*PAPB(2)+ACB(3,IT)*PAPB(3)
     #          +ACB(4,IT)*PAPB(4)+ACB(5,IT)*PAPB(5)+ACB(6,IT)*PAPB(6)
     #          +ACB(7,IT)*PAPB(7)+ACB(8,IT)*PAPB(8)+ACB(9,IT)*PAPB(9)
     #          +ACB(10,IT)*PAPB(10)+ACB(11,IT)*PAPB(11)
7200	CONTINUE
C	TRANSVERSE SHEAR CONTRIBUTION
C	-----------------------------
	IF (IPEL.EQ.1.OR.MTMOD.EQ.2) THEN
	  APSH(1) = APS(1)+APS(2)
	  APSH(2) = APS(3)+APS(4)
	  PAPSH= MATMUL(PTP,APSH)
	ENDIF
	DO 7300 IR = 3,21,6
	  IS = IR+1
	  IT = IR+2
	  FS(IR) = ACS(1,IR)*PAPSH(1)+ACS(2,IR)*PAPSH(2)
	  FS(IS) = ACS(1,IS)*PAPSH(1)+ACS(2,IS)*PAPSH(2)
	  FS(IT) = ACS(1,IT)*PAPSH(1)+ACS(2,IT)*PAPSH(2)
7300	CONTINUE

	IF (MTMOD.EQ.2) THEN
C	 MEMBRANE-BENDING CONTRIBUTION
C	 -----------------------------
	 PAPM = MATMUL(PMPB,APB)
	 DO 7400 IR = 1,19,6
	   IS = IR+1
	   FMB(IR) = ACM(1,IR)*PAPM(1)+ACM(2,IR)*PAPM(2)+ACM(5,IR)*PAPM(5)
	   FMB(IS) = ACM(3,IS)*PAPM(3)+ACM(4,IS)*PAPM(4)+ACM(5,IS)*PAPM(5)
7400	 CONTINUE
C	 BENDING-MEMBRANE CONTRIBUTION
C	 -----------------------------
	 PAPB = MATMUL(PBPM,APMT)
	 DO 7500 IR = 3,21,6
	   IS = IR+1
	   IT = IR+2
	   FBM(IR) = ACB(1,IR)*PAPB(1)+ACB(2,IR)*PAPB(2)+ACB(3,IR)*PAPB(3)
     #            +ACB(4,IR)*PAPB(4)+ACB(5,IR)*PAPB(5)+ACB(6,IR)*PAPB(6)
     #            +ACB(7,IR)*PAPB(7)+ACB(8,IR)*PAPB(8)+ACB(9,IR)*PAPB(9)
     #            +ACB(10,IR)*PAPB(10)+ACB(11,IR)*PAPB(11)
	   FBM(IS) = ACB(1,IS)*PAPB(1)+ACB(2,IS)*PAPB(2)+ACB(3,IS)*PAPB(3)
     #            +ACB(4,IS)*PAPB(4)+ACB(5,IS)*PAPB(5)+ACB(6,IS)*PAPB(6)
     #            +ACB(7,IS)*PAPB(7)+ACB(8,IS)*PAPB(8)+ACB(9,IS)*PAPB(9)
     #            +ACB(10,IS)*PAPB(10)+ACB(11,IS)*PAPB(11)
	   FBM(IT) = ACB(1,IT)*PAPB(1)+ACB(2,IT)*PAPB(2)+ACB(3,IT)*PAPB(3)
     #            +ACB(4,IT)*PAPB(4)+ACB(5,IT)*PAPB(5)+ACB(6,IT)*PAPB(6)
     #            +ACB(7,IT)*PAPB(7)+ACB(8,IT)*PAPB(8)+ACB(9,IT)*PAPB(9)
     #            +ACB(10,IT)*PAPB(10)+ACB(11,IT)*PAPB(11)
7500	 CONTINUE
	ENDIF


C	---------------------------------------
C	FORCES FROM DRILLING DOF SONGSAK JAN09
      IF(NLOPT.EQ.0) THEN
	FD = 0.0D0
	DFAC = 10.0*DR(19)*AB1(1,1)
	DO IEF = 1,NEF
	DO JEF = 1,NEF
	FD(IEF) = FD(IEF) + ACD(IEF)*DFAC*ACD(JEF)*TEDIS(JEF)
	ENDDO
	ENDDO
	ENDIF
C	---------------------------------------


C	SUM CONTRIBUTIONS
C	-----------------
	IF (MTMOD.NE.2) THEN
	  DO 333 IDF = 1,NEF
333	    RETMP(IDF) = FM(IDF)+FB(IDF)+FS(IDF)+FD(IDF)   !FD ADDED BY SONGSAK FOR DRILLING DOF JAN09
	ELSE
	  DO 335 IDF = 1,NEF
335	    RETMP(IDF) = FM(IDF)+FMB(IDF)+FBM(IDF)+FB(IDF)+FS(IDF)
	ENDIF

C	TRANSFORM TO GLOBAL
C	-------------------
	DO 9017 I = 1,22,3
	  IX = I+1
	  IY = I+2
	  RE(I)  = TT(1,1)*RETMP(I)+TT(2,1)*RETMP(IX)+TT(3,1)*RETMP(IY)
	  RE(IX) = TT(1,2)*RETMP(I)+TT(2,2)*RETMP(IX)+TT(3,2)*RETMP(IY)
	  RE(IY) = TT(1,3)*RETMP(I)+TT(2,3)*RETMP(IX)+TT(3,3)*RETMP(IY)
9017	CONTINUE
	DO 9018 I = 4,22,6
	  IX = I+1
	  IY = I+2
	  IZ = I-3
	  IW = I-2
	  RE(I ) = RE(I )+TT(IZ,I )*RETMP(IZ)+TT(IW,I )*RETMP(IW)
	  RE(IX) = RE(IX)+TT(IZ,IX)*RETMP(IZ)+TT(IW,IX)*RETMP(IW)
	  RE(IY) = RE(IY)+TT(IZ,IY)*RETMP(IZ)+TT(IW,IY)*RETMP(IW)
9018	CONTINUE


C     ------------------------------------------------------
C     FIND LINEAR CONTRIBUTION TO STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------------------
      IF (IFREF ) 9000,7000,9000

7000	ST = 0.0
	KST = 1
C	ADD FLEXURAL PART OF STIFFNESS MATRIX - S
	 CALL KBEND(ACB,PBPB,ST)
C	ADD MEMBRANE PART OF STIFFNESS MATRIX - S
	 CALL KMEM20(ACM,PMPM,ST)
C	ADD TRANS SHEAR PART OF STIFFNESS MATRIX - S
	 CALL KSHR(ACS,PTP,ST)
C	MEMBRANE-BENDING + BENDING-MEMBRANE PART OF STIFFNESS MATRIX - S
	 IF (IPEL.NE.1.OR.MTMOD.EQ.2) CALL KMEBE(ACM,ACB,PMPB,ST)


C	ADD TORSION PART OF THS STIFFNESS MATRIX K2 - DRILLING STRAIN
	DFAC = 0.0D0!10.0*DR(19)*AB1(1,1)
	CALL KDRLL(ACD,DFAC,ST)

	
C	ADD TORSION PART OF THS STIFFNESS MATRIX K3 - SOFT SPRING
      SFAC=DMIN1(DR(19),DR(28)/(AB1(1,1)/4.0)) !NEW - FEB2004
C	DFAC=DR(19)*1.E-6
	DFAC=0.0D0!AB1(1,1)*SFAC*1.E-6 !CHANGED - FEB2004
	ST(111) = ST(111)+DFAC
	ST(210) = ST(210)+DFAC
	ST(273) = ST(273)+DFAC
	ST(300) = ST(300)+DFAC

C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 8100
C	-- -- -- -- -- -- -- -- -- -- -- -- -
C	inverse(Ag)*Cg FOR NON-LINEAR STRAINS
C	-- -- -- -- -- -- -- -- -- -- -- -- -
8000	IF (KST.EQ.0) ST = 0.0
	CALL CGEO20(AB1,ANT,RN,SN,CB,CM,CS,CSC,RS,SLN,BL,BLR,BLS,
	1                DWR,DWS,DWRS,DWSR,FLR,FLS,FLRR,FLSS,ACG)

C	-- -- -- -- -- -- -- -- -- -- --
C	DETERMINE INT(Ng*F*Pg)drds - PFP
C	-- -- -- -- -- -- -- -- -- -- --
5555	IF (IPEL.EQ.1.AND.MTMOD.NE.2) THEN
	  CALL NGFPG20(AB1,APM,APB,APS,DR,PFP)
	ELSE IF (MTMOD.EQ.2) THEN
	  CALL NFPPLAS2(AB1,APM,APB,APS,DR,PFP)
	ELSE
	  CALL NFPPLAS3(PAPM,PAPB,PAPSH,PFP)
	ENDIF

C	-- -- -- -- -- -- -- -- -- --
C	DETERMINE GEOMETRIC STIFFNESS
C	-- -- -- -- -- -- -- -- -- --
	CALL KGEO20(ACG,PFP,ST)

8100	CALL KGLOBAL(ST,TT,S)

C	TIM(12)=TIM(12)+TIM2
9000	CONTINUE


C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
      ENDIF

C     SAVE THE STIFFNESS MATRIX INTO FULL MATRIX FORM
      II = 25
      DO I = 1,24  
          II = II-1
          NCOLM(I) = II
      ENDDO

      KK = 0
      IJ = 0
      DO II = 1,24
          IJ = IJ + 1
         DO JJ = II,24
             KK = KK + 1
             SM(II,JJ) = S(KK)
         ENDDO    
      ENDDO      
      
      DO II = 1,24
         DO JJ = 1,24
          IF(II.NE.JJ) SM(JJ,II) = SM(II,JJ)
         ENDDO
      ENDDO        
      
      RETURN
      END
C
C======================================================================   
C======================================================================
C======================================================================     
C	An Assumed Natural Transverse Shear Strain Natural-Based Lagrangian 4-Node Shell Element
c	Enhanced Assumed Strain (EAS )Method for Membrane Strain
c	Assumed Strain Method (ANS) Used in Transverse Shear Strains
C	-------------------------------------------------------------
      SUBROUTINE SHELL43SUPER(PROPM,PROPG,WA,S,COORD,
	1					EDIS,EDISI,RE,MWG,ALPHAM,SEDIM,SELM,RHM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(7),LCS,ISOLOP
C
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NDF,NWG,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT

	COMMON /MMENH/ MM,MM1,MM2,NDIMC
C
C     FOR SUPER ELEMENT BY BJ
      COMMON /SUPERSTIF/ SM2(24,24) 
      
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL      
      
      DIMENSION NCOLM(24)        
      
	DIMENSION EAS(3),ALPHAM(MM),RHM(MM),REHM(24)
C
	DIMENSION PROPM(NMP),PROGP(NGP)
	DIMENSION EDIS(24),EDISI(24),RE(24)
C
C	--------		------------
C	VARIABLE		DESCRIPTIONS
C	--------		------------
C
C	COORD		ELEMENT CO-ORDINATES
C	H			SHAPE FUNCTIONS
C	P			SHAPE FUNCTION DERIVATIVE W.R.T. xi,eta
C	DI		    INITIAL UNDEFORMED SURFACE DIRECTORS OF SHELL SURFACE
C	G1			COVARIANT BASIS-1
C	G2			COVARIANT BASIS-2
C	G3			COVARIANT BASIS-3
C	G01			COVARIANT BASIS-1 AT THE ELEMENT CENTER (0,0)
C	G02			COVARIANT BASIS-2 AT THE ELEMENT CENTER (0,0)
C	G03			COVARIANT BASIS-3 AT THE ELEMENT CENTER (0,0)
C	D1			DIRECTOR GRADIENT-1
C	D2			DIRECTOR GRADIENT-2
C	D01			DIRECTOR GRADIENT-1 THE ELEMENT CENTER (0,0)
C	D02			DIRECTOR GRADIENT-2 THE ELEMENT CENTER (0,0)
C	D			DIRECTORS AT DEFINED SURFACES
C	D0			DIRECTORS AT DEFINED ELEMENT CENTER (0,0) SURFACES
C
	DIMENSION COORD(3,4),H(4),P(2,4),XJ(2,2),XJI(2,2),WA(MWG,1)
	DIMENSION VR(3),VS(3),VT(3)
C	
	DIMENSION DI(3,4),EPS(8)
	DIMENSION G1(3),G2(3),G3(3),D1(3),D2(3),D(3)
C	---------------------------------------------
C	NATURAL-BASED COVARIANT OF THE ELEMENT CENTER
C	---------------------------------------------
	DIMENSION G01(3),G02(3),G03(3),D01(3),D02(3),D0(3)
C	-----------------------
C	NATURAL MEMBRANE STRAIN
C	-----------------------
	DIMENSION B11(1,24)
	DIMENSION B22(1,24)
	DIMENSION B12(1,24)
C	---------------------------------------
C	ASSUMED NATURAL TRANSVERSE SHEAR STRAIN
C	---------------------------------------
C	BS13A	TR. SHEAR STRAIN AT SAMPLING POINT A
C	BS13B	TR. SHEAR STRAIN AT SAMPLING POINT B
C	BS23A	TR. SHEAR STRAIN AT SAMPLING POINT C
C	BS23B	TR. SHEAR STRAIN AT SAMPLING POINT D
C
	DIMENSION BS13A(1,24),BS13B(1,24)
	DIMENSION BS23A(1,24),BS23B(1,24)
	DIMENSION BSA(4,24)
	DIMENSION ANS(2,4)
C
	DIMENSION BM(3,24)
	DIMENSION BS(2,24)
	DIMENSION BD(3,24)

	DIMENSION BS13(1,24),BS23(1,24)
C	-----------------------------------------
C	COVATRARIANT BASIS AT INTERGRATION POINTS 
C	-----------------------------------------
	DIMENSION GT1(3),GT2(3),GT3(3)
C	-------------------------------------
C	CONTRAVARIANT BASIS OF ELEMENT CENTER 
C	-------------------------------------
	DIMENSION GT01(3),GT02(3),GT03(3)
C	-------------------------------------------------
C	ASSUMED TRANSVERSE SHEAR STRAIN AT SAMPLING POINT
C	-------------------------------------------------
	DIMENSION RSP(4),SSP(4)
C	--------------------------------
C	RESULTANT CONSTITUTIVE RELATIONS
C	--------------------------------
	DIMENSION DM(3,3),DB(3,3),DS(2,2)
C	--------------------------------------------------------
C	RESULTANT ENHANCED ASSUMED STRAIN CONSTITUTIVE RELATIONS
C	--------------------------------------------------------
	DIMENSION DM0(3,3),DM1(3,3)
	DIMENSION DB0(3,3),DB1(3,3)
C	---------------------------------------
C	ENHANCED ASSUMED STRAIN MEMBRANE MODELS
C	---------------------------------------
	DIMENSION GC(3,MM)
	DIMENSION SEDM(7,7),SELM(7,24),SEDIM(7,7),SEM(24,24)
C	------------------------
C	LINEAR ELEMENT STIFFNESS
C	------------------------
	DIMENSION SM(24,24),SB(24,24),SH(24,24),ST(24,24)
	DIMENSION SD(24,24),BDRL(1,24),DR(1,1)
	DIMENSION S(300)
C	-------------------------------------------------
C	STRESSES, STRAINS AND RESULTANT STRESS COMPONENTS
C	-------------------------------------------------
	DIMENSION EM(3,1),EB(3,1),GAMMAS(2,1)
	DIMENSION EDISD(24,1)
C
C	RESULTANTS STRESSES AND RESULTANT EQUILIBRIUM FORCES
C
C	RN(3,1)			RESULTANT NORMAL STRESS ACCORDING TO MEMBRANE STRAIN
C	RM(3,1)			RESULTANT MOMENT STRESS ACCORDING TO BENDING STRAIN
C	RS(2,1)			RESULTANT TRANSVERSE SHEAR STRESS ACCORDING TO TRANSVERSE SHEAR STRAIN
C	FN(24,1)		EQUILIBRIUM RESULTANT NORMAL FORCES
C	FM(24,1)		EQUILIBRIUM RESULTANT MOMENT FORCES
C	FS(24,1)		EQUILIBRIUM RESULTANT TRANSVERSE SHEAR FORCES
C
C	EAS PARAMETERS
C
C	FHM(MM,1)		ENHANCED EQUILIBRIUM FORCE
C	REHM(24,1)		ENHANCED RESULTANT EQUILIBRIUM FORCE 
C
	DIMENSION RN(3,1),RM(3,1),RS(2,1)
	DIMENSION FN(24,1),FM(24,1),FS(24,1)
C
	DIMENSION FHM(7,1)
C
C	----------------------------------------------
C	STRATING SUBROUTINE XSHELL43 (LINEAR ANALYSIS)
C	----------------------------------------------
	RHM(1:MM) = 0.0D0
	EAS(1:3)  = 0.0D0
C	--------------------------------------
C	INITILIZED EQUILIBRIUM RESULTANT FORCE	
C	--------------------------------------
	FN(1:24,1) = 0.0D0
	FM(1:24,1) = 0.0D0
	FS(1:24,1) = 0.0D0
C
	FHM(1:7,1) = 0.0D0
C	-----------------------------------------------------
C	ASSIGN DISPLACEMENTS AND ROTATIONS TO THE MATRIX FORM
C	-----------------------------------------------------
	EDISD(1:24,1) = 0.0D0
	EDISD(1:24,1) = EDIS(1:24)
C
C	---------------------------
C	ENHANCED MEMBRANE STIFFNESS
C	---------------------------
	SEDM  = 0.0D0
	SELM  = 0.0D0
	SEDIM = 0.0D0
	SEM   = 0.0D0
C	---------------------------------
C	LINEAR ELEMENT STIFFNESS MATRICES
C	---------------------------------
	SM(1:24,1:24) = 0.0D0
	SB(1:24,1:24) = 0.0D0
	SH(1:24,1:24) = 0.0D0
	ST(1:24,1:24) = 0.0D0
	SD(1:24,1:24) = 0.0D0
C
	S(1:300) = 0.0D0
C	-----------------------------------
C	DEFINE UNDEFORMED SURFACE DIRECTORS
C	-----------------------------------
	CALL SURFACEDIRECTOR(PROPG,COORD,DI,TH)
C	---------------------------------------
C	ASSUMED NATURAL TRANSVERSE SHEAR STRAIN
C	---------------------------------------
C	TRANSVERSE SHEAR 13	SAMPLING POINTS A
C	-------------------------------------
	RSP = 0.0D0
	SSP = 1.0D0
	H(1:4)     = 0.0D0
	P(1:2,1:4) = 0.0D0
	BS13A = 0.0D0
	CALL SHAPE4EAS4N(RSP,SSP,H,P,4)
	CALL COVBASE(COORD,DI,P,H,G1,G2,G3,D1,D2,D)
	CALL TRSHEAR13(G1,G2,DI,D,H,P,BS13A)
C	-------------------------------------
C	TRANSVERSE SHEAR 13	SAMPLING POINTS B
C	-------------------------------------
	RSP =  0.0D0
	SSP = -1.0D0
	H(1:4)     = 0.0D0
	P(1:2,1:4) = 0.0D0
	BS13B = 0.0D0
	CALL SHAPE4EAS4N(RSP,SSP,H,P,4)
	CALL COVBASE(COORD,DI,P,H,G1,G2,G3,D1,D2,D)
	CALL TRSHEAR13(G1,G2,DI,D,H,P,BS13B)
C	-------------------------------------
C	TRANSVERSE SHEAR 23 SAMPLING POINTS C
C	-------------------------------------
	RSP = 1.0D0
	SSP = 0.0D0
	H(1:4)     = 0.0D0
	P(1:2,1:4) = 0.0D0
	BS23A = 0.0D0
	CALL SHAPE4EAS4N(RSP,SSP,H,P,4)
	CALL COVBASE(COORD,DI,P,H,G1,G2,G3,D1,D2,D)
	CALL TRSHEAR23(G1,G2,DI,D,H,P,BS23A)
C	-------------------------------------
C	TRANSVERSE SHEAR 23 SAMPLING POINTS D
C	-------------------------------------
	RSP = -1.0D0
	SSP =  0.0D0
	H(1:4)     = 0.0D0
	P(1:2,1:4) = 0.0D0
	BS23B = 0.0D0
	CALL SHAPE4EAS4N(RSP,SSP,H,P,4)
	CALL COVBASE(COORD,DI,P,H,G1,G2,G3,D1,D2,D)
	CALL TRSHEAR23(G1,G2,DI,D,H,P,BS23B)
C	--------------------------------------
C	AN ASSUMED NATURAL SHEAR STRAIN MATRIX
C	--------------------------------------
	BSA(1:4,1:24) = 0.0D0
	DO I=1,24
	BSA(1,I) = BS13A(1,I)
	BSA(2,I) = BS13B(1,I)
	BSA(3,I) = BS23A(1,I)
	BSA(4,I) = BS23B(1,I)
	END DO
C	---------------------------------
C	COVARIANT BASES AT ELEMENT CENTER 
C	---------------------------------
	RI0 = 0.0D0
	SI0 = 0.0D0
	CALL SHAPE4EAS4N(RI0,SI0,H,P,4)
	CALL COVBASE(COORD,DI,P,H,G01,G02,G03,D01,D02,D0)
C	------------------
C	INTEGRATION DOMAIN
C	------------------
C
      IPT = 0
	NGR = 2
	NGS = 2
C
      DO 900  IGR=1,NGR
      RI = GLOC(IGR,NGR)
      DO 900  IGS=1,NGS
      SI = GLOC(IGS,NGS)
C
      WT = GWT(IGR,NGR)*GWT(IGS,NGS)
      IPT= IPT+1
C	------------------------------------
C	INITIALIZED NATURAL BASIS PARAMETERS
C	------------------------------------
	G1 = 0.0D0
	G2 = 0.0D0
	G3 = 0.0D0
	D1 = 0.0D0
	D2 = 0.0D0
	D  = 0.0D0
	P(1:2,1:4) = 0.0D0
	H(1:4)     = 0.0D0
	DET= 0.0D0
	DM = 0.0D0
	DB = 0.0D0
	DS = 0.0D0
	CALL SHAPE4EAS4N(RI,SI,H,P,4)
	CALL COVBASE(COORD,DI,P,H,G1,G2,G3,D1,D2,D)
	CALL CTVBASE(G1,G2,G3,D,GT1,GT2,GT3,DET)      ! CONTRAVIRIANT BASIS
	CALL NATURALCONS(PROPM,GT1,GT2,GT3,DM,DB,DS)  ! THE NATURAL-BASED CONSTITUTIVE RELATION
C
	DAREA = WT*DET
C
	CALL SHAPE4EAS4N(RI,SI,H,P,4)
	CALL COVBASE(COORD,DI,P,H,G1,G2,G3,D1,D2,D)
	
	CALL SHSTORL1 (COORD,PROPG,WA(1,IPT),NWG,VR,VS,VT)	!STORE SHELL LAX FOR STRESS CALCULATION JAN09
	
C	---------------------------------------------
C	NATURAL MEMBRANE STRAIN GRADIENT DISPLACEMENT
C	---------------------------------------------
	B11(1,1:24) = 0.0D0
	B22(1,1:24) = 0.0D0
	B12(1,1:24) = 0.0D0
	CALL STRAIN11(G1,G2,P,B11)
	CALL STRAIN22(G1,G2,P,B22)
	CALL GAMMA12 (G1,G2,P,B12)
	BM(1,1:24) = B11(1,1:24)
	BM(2,1:24) = B22(1,1:24)
	BM(3,1:24) = B12(1,1:24)
C	-----------------------------------------------------
C	ENHANCED ASSUMED STRAIN MEMBRANE 7 EAS TERMS ARE USED 
C	-----------------------------------------------------
	GC(1:3,1:7) = 0.0D0
	GC(1,1) = RI
	GC(2,2) = SI
	GC(3,3) = SI
	GC(3,4) = RI
	GC(1,5) = RI*SI
	GC(2,6) = RI*SI
	GC(3,7) = RI*SI
C	----------------------------
C	ASSUMED NATURAL SHEAR STRAIN
C	MODIFIED TRANSVERSE SHEAR AT SAMPLING POINTS
C	--------------------------------------------
	ANS(1:2,1:4) = 0.0D0
	ANS(1,1) = (1.0D0/2.0D0)*(1.0D0+SI)
	ANS(1,2) = (1.0D0/2.0D0)*(1.0D0-SI)
	ANS(2,3) = (1.0D0/2.0D0)*(1.0D0+RI)
	ANS(2,4) = (1.0D0/2.0D0)*(1.0D0-RI)
C
	BS(1:2,1:24) = 0.0D0
	BS = MATMUL(ANS,BSA)
C
C	-----------------------------------
C	COMPATIBLE TRANSVERSE SHEAR STRAINS
C	-----------------------------------	
C	BS13 = 0.0D0
C	BS23 = 0.0D0
C	CALL TRSHEAR13(G1,G2,DI,D,H,P,BS13)
C	CALL TRSHEAR23(G1,G2,DI,D,H,P,BS23)
C	BS(1,1:24) = BS13(1,1:24)
C	BS(2,1:24) = BS23(1,1:24)
C	--------------
C	BENDING STRAIN
C	--------------
	BD = 0.0D0
	CALL CURVATURE(DI,D1,D2,D,G1,G2,P,BD)
C	------------------
C	INTEGRAL OVER AREA
C	------------------
	SM = SM + MATMUL(TRANSPOSE(BM),MATMUL(DM,BM))*DAREA
	SB = SB + MATMUL(TRANSPOSE(BD),MATMUL(DB,BD))*DAREA
	SH = SH + MATMUL(TRANSPOSE(BS),MATMUL(DS,BS))*DAREA
C	----------------------------------------
C	ENHANCED ASSUMED STRAIN BASIC PARAMETERS
C	----------------------------------------
	CALL CTVBASE(G01,G02,G03,D0,GT01,GT02,GT03,DET0)         ! CONTRAVIRIANT BASED
	CALL ENHANCEDCON(PROPM,GT01,GT02,GT03,DM0,DB0)           ! CENTER OF JACOBIAN
	CALL COUPLEMEM(PROPM,GT1,GT2,GT3,GT01,GT02,GT03,DM1,DB1) ! CENTER OF JACOBIAN
C	------------------------
C	ENHANCED MEMBRANE STRAIN
C	------------------------
	SEDM =SEDM+MATMUL(TRANSPOSE(GC),MATMUL(DM0,GC))*DAREA
	1              *(DET0/DET)*(DET0/DET)
	SELM =SELM+MATMUL(TRANSPOSE(GC),MATMUL(DM1,BM))*DAREA*(DET0/DET)
C
C	----------------------------------------
C	UPDATING STRAIN LINEAR STRAIN COMPONENTS
C	----------------------------------------
	EM     = 0.0D0
	EB     = 0.0D0
	GAMMAS = 0.0D0
C
	EM     = MATMUL(BM,EDISD)
	EB     = MATMUL(BD,EDISD)
	GAMMAS = MATMUL(BS,EDISD)
C	-----------------------------------
C	ADDITIONAL ENHANCED STRAIN MEMBRANE
C	-----------------------------------
	EAS = 0.0D0
	EAS = MATMUL(GC,ALPHAM)
	EM(1,1) = EM(1,1)+EAS(1)	
	EM(2,1) = EM(2,1)+EAS(2)	
	EM(3,1) = EM(3,1)+EAS(3)	
	
      EPS(1:3) = EM(1:3,1)
      EPS(4:6) = EB(1:3,1)
      EPS(7:8) = GAMMAS(1:2,1)
	CALL SHSTRS(VR,VS,VT,EPS,WA(9,IPT),PROPM,TH)	
	

C
C	-----------------------------------------------
C	COMPUTE NATURAL-BASED CAUCHY RESULTANT STRESSES
C	------------------------------------
C	RESULTANT MEMBRANE STRESS (NORMAL STRESS RESULTANTS)
C	RESULTANT BENDING STRESS  (MOMENT STRESS RESULTANTS)
C	RESULTANT TRANSVERSE SHEAR STRESS (TRANSVERSE SHEAR STRESS RESULTANTS)
C	----------------------------------------------------------------------
	RN(1:3,1) = 0.0D0
	RM(1:3,1) = 0.0D0
	RS(1:2,1) = 0.0D0
C
	RN = MATMUL(DM,EM)
	RM = MATMUL(DB,EB)
	RS = MATMUL(DS,GAMMAS)
	
	
      WA(1:3,IPT)=RN(1:3,1)
      WA(4:6,IPT)=RM(1:3,1)
      WA(7:8,IPT)=RS(1:2,1)
 
C	------------------------------------
C	COMPATIBLE ELEMENT EQULIBRIUM FORCES
C	------------------------------------
	FN = FN + MATMUL(TRANSPOSE(BM),RN)*DAREA
	FM = FM + MATMUL(TRANSPOSE(BD),RM)*DAREA
	FS = FS + MATMUL(TRANSPOSE(BS),RS)*DAREA
C	---------------------------------------------------
C 	ENHANCED ASSUMED MEMBRANE STRAIN EQUILIBRIUM FORCES
C	---------------------------------------------------
	FHM = FHM + MATMUL(TRANSPOSE(GC),RN)*DAREA
 900  CONTINUE
C	----------------------------------------
C	ASSEMBLY THE RESULTANT EQUILIBRIUM FORCE
C	----------------------------------------
	RE(1:24) = FN(1:24,1)+FM(1:24,1)+FS(1:24,1)
C	------------------------
C	ENHANCED MEMBRANE STRAIN
C	------------------------
	CALL INVMAT1(SEDM,SEDIM,7)
	SEM = 0.0D0
	SEM = MATMUL(TRANSPOSE(SELM),MATMUL(SEDIM,SELM))
	SM = SM-SEM
C	---------------------------------
C	EQUILIBRIUM DUE TO ENHANCED FORCE
C	---------------------------------
	RHM(1:7) = FHM(1:7,1)
	REHM = MATMUL(MATMUL(TRANSPOSE(SELM),SEDIM),RHM)
C	-------------------------
C	FIND THE EQUILIBRIUM LOAD
C	-------------------------
	DO I=1,24
	RE(I) = RE(I)-REHM(I)
	END DO
C	-------------------------------------------------------------
C	ASSEMBLY ELEMENT STIFFNESS WITHOUT DRILLING DEGREE OF FREEDOM
C	-------------------------------------------------------------
	ST = SM+SB+SH
C	------------------------------------------------------------
C	SPLITING ENERGY OF DRILLING DEGREES OF FREEDOM CONTRIBUTIONS
C	------------------------------------------------------------
C
	NGR = 2
	NGS = 2
C
C	CHOICE OF PENALTY TESTING
C	PENALTY = 0.5E-5
C	PENALTY = 0.01D0
C	PENALTY = 1.0D-05
	PENALTY = 0.5E-10 !LINEAR ANALYSIS FOR THIS VALUE
C
      DO 901  IGR=1,NGR
      RI = GLOC(IGR,NGR)
      DO 901  IGS=1,NGS
      SI = GLOC(IGS,NGS)
C
      WT = GWT(IGR,NGR)*GWT(IGS,NGS)
      IPT= IPT+1
C
C	------------------------
C	INITILAZED WORKING ARRAY
C	SHAPE FUNCTION
C	SHAPE FUNCTION DERIVATIVES
C	VECTOR COMPONENT VR,VS,VT
C	-------------------------

	H(1:4)     = 0.0D0
	P(1:2,1:4) = 0.0D0 
	VR(1:3)	   = 0.0D0
	VS(1:3)	   = 0.0D0
	VT(1:3)	   = 0.0D0

C	--------------------------------------------
C	COMPUTE JACOBIAN MATRIX ,INVERSE OF JACOBIAN
C	VECTOR COMPONENTS AND DETERMINANT
C	---------------------------------
	CALL SHAPE4EAS4N (RI,SI,H,P,4)
	CALL JACOBSH4(COORD,P,XJ,XJI,DET,VR,VS,VT)
      DVOL=WT*DET
C	----------------------------------------
C	DRILLING DEGREES OF FREEDOM CONTRIBUTION
C	----------------------------------------
	BDRL(1:1,1:24) = 0.0D0
	CALL DRILING(XJI,P,H,VR,VS,VT,BDRL)
	DR(1,1) = PENALTY*(1./2.)*(PROPM(1)/(1.+ PROPM(2)))
C
	SD = SD + MATMUL(TRANSPOSE(BDRL),MATMUL(DR,BDRL))*DVOL
C
 901  CONTINUE
C
C	ADDED DRILLING CONTRIBUTIONS FOR SINGULARITY FREE
C
	ST = ST + SD
C
C	--------------------------------
C	MAKING UPPER TRIANGULAR ROW-WISE
C	-------------------------------- 

	KDX = 0

	DO I =1,24
	 DO J=I,24
	   KDX = KDX + 1
	  S(KDX) = ST(I,J)
	END DO
      END DO
C
C     SAVE THE STIFFNESS MATRIX INTO FULL MATRIX FORM
      II = 25
      DO I = 1,24  
          II = II-1
          NCOLM(I) = II
      ENDDO

      KK = 0
      IJ = 0
      DO II = 1,24
          IJ = IJ + 1
         DO JJ = II,24
             KK = KK + 1
             SM2(II,JJ) = S(KK)
         ENDDO    
      ENDDO      
      
      DO II = 1,24
         DO JJ = 1,24
          IF(II.NE.JJ) SM2(JJ,II) = SM2(II,JJ)
         ENDDO
      ENDDO   
      
      
1500	FORMAT(50F12.5)
	RETURN
      END
C
C-----------------------------------------------------------------------------      
      SUBROUTINE MESTIF_CON(KREC,SUU,SUC,SCU,SCC,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*3 OPER
C	=============================================================
C	DIMENSION S(1),SUU(1,1),SUC(1,1),SCU(1,1),SCC(1,1)
      DIMENSION S(1),SUU(24,24),SUC(24,126),SCU(126,24),SCC(126,126)


	IND = -1
	IF(OPER.EQ.'RED') IND = 0 
	IF(OPER.EQ.'WRT') IND = 1
	IF(OPER.EQ.'NOO') IND = 2 !PASS
	IF(IND.EQ.-1) RETURN
	
      NEQU = 24
      NEQC = 126
      
      N1 = NEQU*NEQU
      N2 = NEQU*NEQC
      N3 = NEQC*NEQU
      N4 = NEQC*NEQC
      
C      OPEN(UNIT=1563,FILE="KAK.TXT")
      
	SELECT CASE(IND)

C	---------------------------------------
      CASE(0)     
           READ(1563,100) SUU(1:NEQU,1:NEQU)!,SUC(1:NEQU,1:NEQC),SCU(1:NEQC,1:NEQU),SCC(1:NEQC,1:NEQC)
 
C	---------------------------------------
      CASE(1)
          WRITE(1563,200) SUU(1:NEQU,1:NEQU)!,SUC(1:NEQU,1:NEQC),SCU(1:NEQC,1:NEQU),SCC(1:NEQC,1:NEQC)
C	---------------------------------------
      CASE(2)
          CONTINUE
C	---------------------------------------

	ENDSELECT	 


	RETURN
100   FORMAT (F20.4) 
200   FORMAT (F20.4)
      END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE INVMATRIX_ONE(A,B,N,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------
C     MATRIX INVERSION BY ELIMMINATION WITH PARTIAL PIVOTING
C     AND DETERMINANT OF MATRIX A
C     ROGINAL MATRIX =A, INVERSE MATRIX=B
C
C     A: INPUT SQUARE MATRIX
C     B: INVERSE MATRIX OF A
C     N: SIZE OF A AND B
C     EPS: CONTROL VARIABLE
C     DET: DETERMINANT OF A
C     --------

      DIMENSION A(N,N2),B(N,N2),C(N,N2)

C     CONSTRUCT IDENTITY MATRIX B(I,J)=I
      EPS = 1.0E-12
	DO 88 I=1,N
	DO 88 J=1,N2
	C(I,J)=A(I,J)
  88  CONTINUE	
 	DO 6 I=1,N
	DO 5 J=1,N2
	IF(I-J) 4,3,4
   3  B(I,J)=1.0
      GOTO 5 
   4  B(I,J)=0.0
   5  CONTINUE
   6  CONTINUE

C     LOCATE MAXIMUM MAGNITUDE A(I,K) ON OR BELOW MAIN DIAGONAL
      DET=1.0
	DO 45 K=1,N
	IF (K-N) 12,30,30
   12 IMAX=K
      AMAX=ABS(C(K,1))
	KP1=K+1
	DO 20 I=KP1,N
	IF (AMAX-ABS(C(I,1))) 15,20,20
   15 IMAX=I
      AMAX=ABS(C(I,1))
   20 CONTINUE

C     INTERCHANGE ROWS IMAX AND K IF IMAX NOT EQUAL TO K
      IF (IMAX-K) 25,30,25
   25 DO 29 J=1,N2
      ATMP=C(IMAX,J)
	C(IMAX,J)=C(K,J)
	C(K,J)=ATMP
	BTMP=B(IMAX,J)
	B(IMAX,J)=B(K,J)
   29 B(K,J)=BTMP
      DET=-DET
   30 CONTINUE

C     TEST FOR SINGULAR MATRIX
C      IF (ABS(C(K,K))-EPS) 33,33,35
C   33 WRITE(*,*) 'SINGULAR MATRIX EPS- INSIDE SUB INVMATRIX'
C	STOP
   35 DET=C(K,1)*DET

C     DIVIDE PIVOT ROW BY ITS MAIN DIAGONAL ELEMENT
      DIV=C(K,1)
	
	IF(DIV.EQ.0.0) DIV = EPS
	IF(ABS(DIV).LT.EPS) DIV = EPS*DIV/ABS(DIV)

	DO 38 J=1,N2
	C(K,J)=C(K,J)/DIV
   38 B(K,J)=B(K,J)/DIV

C     REPLACE EACH ROW BY LINEAR COMBINATION WITH PIVOT ROW
      DO 43 I=1,N
	AMULT=C(I,1)
	IF (I-K) 39,43,39
   39 DO 42 J=1,N2
      C(I,J)=C(I,J)-AMULT*C(K,J)
   42 B(I,J)=B(I,J)-AMULT*B(K,J)
   43 CONTINUE
   45 CONTINUE
C
      RETURN
      END
C	=========================================================      
      SUBROUTINE MATRIX_CONDEN_TEST
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     STATIC CONDENSTAION OF MATRIX
C     SKS : STIFFNESS MATRIX TO BE REDUCED
C     R : LOAD VECTOR
C     NEQ : TOTAL NUMBER OF EQUATION
C     LEQ : NUMBER OF EQUATION TO BE ELIMINATED      
C     ------------------------------------------------------------
      CHARACTER*3 OPER
       
C     SAVE IEL
      COMMON /INELE/ IEL   
      COMMON /JNELE/ JEL
      
C     STIFFNESS MATRIX PART
      DIMENSION SKUU(2,2),SKUC(2,1),SKCC(1,1),SKCU(1,2),
     1          SKCCI(1,1), ! INVERSION OF SKCC MATRIX
     2          SKSC(2,2),  ! CONDENSED STIFFNESS MATRIX
     3          SKSM(2,2),SKCIU(2,1),TB(1,2),T(3,2)
      DIMENSION FMASS(3,3),FMASSC(2,2)
      
C     1. REARRANGE THE STIFFNESS MATRIX (SUPER NODES AND ELIMINATED NODES)
      
      NEQU = 2
      NEQC = 1
      
      CALL CLEARMAT(SKUU,2,2)
      CALL CLEARMAT(SKUC,2,1)
      CALL CLEARMAT(SKCC,1,1)
      CALL CLEARMAT(SKCU,1,2)
      CALL CLEARMAT(SKCCI,1,1)
      CALL CLEARMAT(SKCIU,2,1)
      CALL CLEARMAT(SKSC,1,1)
      CALL CLEARMAT(SKSM,2,2)

      CALL CLEARMAT(FMASS,3,3)
      CALL CLEARMAT(FMASSC,2,2)
      CALL CLEARMAT(T,3,2)
      
      SKUU(1,1) = 17500
      SKUU(1,2) = -10000
      SKUU(2,1) = -10000
      SKUU(2,2) = 10000
      
      SKUC(1,1) = 0.0
      SKUC(2,1) = 0.0
      
      SKCU(1,1) = -0.25
      SKCU(1,2) = 0.0
      
      SKCC(1,1) = 1
        
      FMASS(1,1) =  100
      FMASS(2,2) =  50
      FMASS(3,3) =  25
      
      
C     2. INVERSE SKCC MATRIX
         
      CALL INVMATRIX(SKCC,SKCCI,1)
C      CALL INVMATRIX(AKK,AKK2,4)

C     3. CONDENSATION OF STIFFNESS MATRIX  
C     =========================================================      
C     MULTIPLYING MATRICES USING 'DGEMM'
C     ---------------------------------------------------------
C     C = ALPHA*A*B + BETA*C   
C     A(M,K), B(K,N), C(M,N)   
C     "INITAILIZE C MATRIX BEFORE DOING CALL DGEMM"
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)     
C     =========================================================
      ALPHA = 1.0
      BETA  = 1.0
            
      M = NEQU
      N = NEQC
      K = NEQC
      CALL DGEMM('N','N',M,N,K,ALPHA,SKUC,M,SKCCI,K,BETA,SKCIU,M)
      
      M = NEQU
      N = NEQU
      K = NEQC     
      CALL DGEMM('N','N',M,N,K,ALPHA,SKCIU,M,SKCU,K,BETA,SKSM,M)    
      
      SKSC = SKUU - SKSM
      
      CALL INVMATRIX(SKCC,SKCCI,1)      
      
C     TB-MATRIX          
      M = NEQC
      N = NEQU
      K = NEQC     
      CALL DGEMM('N','N',M,N,K,ALPHA,-SKCCI,M,SKCU,K,BETA,TB,M)  
      
      T(1,1) = 1.0D0
      T(2,2) = 1.0D0
      T(3,1) = TB(1,1)
      T(3,2) = TB(1,2)
      
      FMASSC = MATMUL(TRANSPOSE(T),MATMUL(FMASS,T))
    
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================           
      SUBROUTINE SAVET(T,IEL,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     STATIC CONDENSTAION OF MATRIX
C     SKS : STIFFNESS MATRIX TO BE REDUCED
C     R : LOAD VECTOR
C     NEQ : TOTAL NUMBER OF EQUATION
C     LEQ : NUMBER OF EQUATION TO BE ELIMINATED      
C     ------------------------------------------------------------
      CHARACTER*3 OPER
      
c      COMMON /TMA/ TMATRIX()
      
      DIMENSION T(1)
      
c      TMATRIX(IEL,1) = T

    
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================       
      SUBROUTINE CALTEDIS(T,EDIS,TEDIS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     STATIC CONDENSTAION OF MATRIX
C     SKS : STIFFNESS MATRIX TO BE REDUCED
C     R : LOAD VECTOR
C     NEQ : TOTAL NUMBER OF EQUATION
C     LEQ : NUMBER OF EQUATION TO BE ELIMINATED      
C     ------------------------------------------------------------
          
C     SAVE IEL
      COMMON /INELE/ IEL   
      COMMON /JNELE/ JEL
      
C     STIFFNESS MATRIX PART
      DIMENSION TB(126,24),T(150,24),TP(126,1),
     1          TEDIS(150,1),EDIS(24)


      NEQU = 24
      NEQC = 126
C     =========================================================      
C     MULTIPLYING MATRICES USING 'DGEMM'
C     ---------------------------------------------------------
C     C = ALPHA*A*B + BETA*C   
C     A(M,K), B(K,N), C(M,N)   
C     "INITAILIZE C MATRIX BEFORE DOING CALL DGEMM"
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)     
C     =========================================================
      ALPHA = 1.0
      BETA  = 1.0 
      
      TB(1:126,1:24) = T(25:150,1:24)
      
C     DISPLACEMENT OF 25 NODES IN ONE SUPER ELEMENT
      M = NEQC
      N = 1
      K = NEQU     
      CALL DGEMM('N','N',M,N,K,ALPHA,TB,M,EDIS,K,BETA,TP,M) 
      
      TEDIS(1:24,1) = EDIS(1:24)
      TEDIS(25:150,1) = TP(1:126,1)
   
C      M = NEQC+NEQU
C      N = 1
C      K = NEQU     
C      CALL DGEMM('N','N',M,N,K,ALPHA,T,M,DISP,K,BETA,TEDIS,M)     
      
      
      IF(ILO.EQ.0) RETURN    
      
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================       
      SUBROUTINE COORDBEAM(COORDMAT)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     READ COORDINATE FROM COORDNEW FOR DIAPHRAGM CONSTRAINTS BY BJ
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NODEDATA      = SAVE NODE DATA IN EACH FLOOR
C     -----------------------------------------------------      
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      
      COMMON A(7000000),IA(6000000)
      
      
      DIMENSION COORDMAT(NSN,3)
     

        COORDMAT(1:NSN,1) = A(LXY:LXY+(NSN-1))
        COORDMAT(1:NSN,2) = A(LXY+NSN:LXY+(2*NSN-1))
        COORDMAT(1:NSN,3) = A(LXY+(2*NSN):LXY+(3*NSN-1))
        
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================      
      SUBROUTINE FICBEAM(INEL,XYZM,PROPM,PROPG,STIFFA,NODEMAT)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     CALCULATE THE STIFFNESS MATRIX OF FINE MESHED FICTITIOUS BEAM
C	------------------------------------------------------------------------    
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      
      COMMON A(7000000),IA(6000000)
      
      
      DIMENSION XYZM(25,3),NODEBEAM(4,2),XYZMESH(2,3)
      DIMENSION PROPM(1),PROPG(1),STIFFBEAM(12,12),STIFFA(30,30)
      DIMENSION IDBEAM(5,6),IEDOFB(12),NODEMAT(5)
     
      IF(INEL.EQ.1) THEN  !SIDE NUMBER IN ONE SUPER ELEMENT
          NOD1  =  1
          NOD2  =  9
          NOD3  =  5
          NOD4  = 10
          NOD5  =  2
      ELSEIF(INEL.EQ.2) THEN    
          NOD1  =  2
          NOD2  = 11
          NOD3  =  6
          NOD4  = 12
          NOD5  =  3          
      ELSEIF(INEL.EQ.3) THEN    
          NOD1  =  3
          NOD2  = 13 
          NOD3  =  7
          NOD4  = 14
          NOD5  =  4  
      ELSEIF(INEL.EQ.4) THEN    
          NOD1  =  4
          NOD2  = 15 
          NOD3  =  8
          NOD4  = 16
          NOD5  =  1
      ELSEIF(INEL.EQ.5) THEN           
          NOD1  =  2
          NOD2  =  10
          NOD3  =  5
          NOD4  =  9
          NOD5  =  1
      ELSEIF(INEL.EQ.6) THEN    
          NOD1  =  3
          NOD2  = 12
          NOD3  =  6
          NOD4  = 11
          NOD5  =  2          
      ELSEIF(INEL.EQ.7) THEN    
          NOD1  = 4
          NOD2  = 14
          NOD3  = 7
          NOD4  = 13
          NOD5  = 3 
      ELSEIF(INEL.EQ.8) THEN    
          NOD1  =  1
          NOD2  =  16
          NOD3  =  8
          NOD4  = 15
          NOD5  =  4         
      ENDIF    
      
      
          NODEMAT(1) = NOD1 !SAVE THE NODE NUMBER INTO NODEMAT MATRIX    
          NODEMAT(2) = NOD2
          NODEMAT(3) = NOD3
          NODEMAT(4) = NOD4
          NODEMAT(5) = NOD5
      
          NODEBEAM(1,1) = NOD1 !1ST NODE NUMBER OF ELEMENT NO.1 IN INEL SIDE
          NODEBEAM(1,2) = NOD2 !2ND NODE NUMBER OF ELEMENT NO.1 IN INEL SIDE
          
          NODEBEAM(2,1) = NOD2 !1ST NODE NUMBER OF ELEMENT NO.2 IN INEL SIDE
          NODEBEAM(2,2) = NOD3 !2ND NODE NUMBER OF ELEMENT NO.2 IN INEL SIDE
          
          NODEBEAM(3,1) = NOD3 !1ST NODE NUMBER OF ELEMENT NO.3 IN INEL SIDE
          NODEBEAM(3,2) = NOD4 !2ND NODE NUMBER OF ELEMENT NO.3 IN INEL SIDE
          
          NODEBEAM(4,1) = NOD4 !1ST NODE NUMBER OF ELEMENT NO.4 IN INEL SIDE
          NODEBEAM(4,2) = NOD5 !2ND NODE NUMBER OF ELEMENT NO.4 IN INEL SIDE
      
      CALL CLEARMAT(STIFFA,30,30)

      DO INODE = 1,5
          IDBEAM(INODE,1) = 6*INODE-5
          IDBEAM(INODE,2) = 6*INODE-4
          IDBEAM(INODE,3) = 6*INODE-3
          IDBEAM(INODE,4) = 6*INODE-2
          IDBEAM(INODE,5) = 6*INODE-1
          IDBEAM(INODE,6) = 6*INODE-0
      ENDDO
      
      NT = 0
      DO IFL = 1,4 !================ START OF LOOP FOR FICTITIOUS BEAM BY BJ
          NODE1 = NODEBEAM(IFL,1)
          NODE2 = NODEBEAM(IFL,2)
          
          XYZMESH(1,1:3) = XYZM(NODE1,1:3)
          XYZMESH(2,1:3) = XYZM(NODE2,1:3)
          
          CALL MODULEBEAM(XYZMESH,PROPM,PROPG,STIFFBEAM) ! MODULE FOR CALCULATING STIFFNESS MATRIX OF FINE MESH ELEMENT OF FICITIOUS BEAM 

          DO KEL = 1,2
                  NT = NT + 1
                  DO IS = 1,6
                      IDO = IDBEAM(NT,IS)
                      IEDOFB(6*KEL-(6-IS)) = IDO 
                  ENDDO    
          ENDDO 
          NT = NT - 1 
          
          !ASSEMBLING THE STIFFNESS MATRIX IN ONE SIDE OF FICTITIOUS BEAM
          DO IEF = 1,12
                IEQ = IEDOFB(IEF)
              DO JEF = 1,12
                  JEQ = IEDOFB(JEF)
                  STIFFA(IEQ,JEQ) = STIFFA(IEQ,JEQ) + STIFFBEAM(IEF,JEF)
              ENDDO
          ENDDO  
          
      ENDDO !================ END OF LOOP FOR FICTITIOUS BEAM BY BJ   
          
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================       
      SUBROUTINE MODULEBEAM(XYZMESH,PROPM,PROPG,STIFFBEAM)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     CALCULATE THE STIFFNESS MATRIX OF FINE MESHED FICTITIOUS BEAM
C	------------------------------------------------------------------------
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      
      COMMON A(7000000),IA(6000000)
      
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV      

      COMMON /HEI/ HEIGHT ! HEIGHT OF EACH SUPER ELEMENT
      
      DIMENSION XYZMESH(2,3),STIFFBEAM(12,12),TRANS(12,12),STIFFBEAM1(12,12)
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION PROPM(1),PROPG(1)
     
      X1 = XYZMESH(1,1)
      Y1 = XYZMESH(1,2)
      Z1 = XYZMESH(1,3)

      X2 = XYZMESH(2,1)
      Y2 = XYZMESH(2,2)
      Z2 = XYZMESH(2,3)     
      
      ELENGTH = SQRT((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
      
      CALL CLEARMAT(STIFFBEAM,12,12)
      
      E  = PROPM(1)
      PO = PROPM(2)
      TH = PROPG(2)
      
      G = E/(2D0*(1D0+PO))
C      TH = HEIGHT

      AREA = TH*HEIGHT
      EIZ = E*(TH*(HEIGHT**3D0))/12D0 ! E*IZ
      EIY = E*((TH**3D0)*HEIGHT)/12D0 ! E*IY
      
      PM = (EIZ + EIY)/E  !POLAR MOMENT OF INERTIA
      
      STIFFBEAM(1,1)   =  E*AREA/ELENGTH
      STIFFBEAM(2,2)   =  12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(3,3)   =  12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(4,4)   =  G*PM/ELENGTH
      STIFFBEAM(5,5)   =  4D0*EIY/ELENGTH
      STIFFBEAM(6,6)   =  4D0*EIZ/ELENGTH
      
      STIFFBEAM(7,7)   =  E*AREA/ELENGTH
      STIFFBEAM(8,8)   =  12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(9,9)   =  12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(10,10) =  G*PM/ELENGTH
      STIFFBEAM(11,11) =  4D0*EIY/ELENGTH
      STIFFBEAM(12,12) =  4D0*EIZ/ELENGTH      
      
      STIFFBEAM(7,1)   =  -E*AREA/ELENGTH
      
      STIFFBEAM(6,2)   =   6D0*EIZ/(ELENGTH**2D0)
      STIFFBEAM(8,2)   =  -12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(12,2)  =   6D0*EIZ/(ELENGTH**2D0)
      
      STIFFBEAM(5,3)   =  -6D0*EIY/(ELENGTH**2D0)
      STIFFBEAM(9,3)   =  -12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(11,3)  =  -6D0*EIY/(ELENGTH**2D0)
      
      STIFFBEAM(10,4)  =  -G*PM/ELENGTH
      
      STIFFBEAM(9,5)   =   6D0*EIY/(ELENGTH**2D0)
      STIFFBEAM(11,5)  =   2D0*EIY/ELENGTH
      
      STIFFBEAM(8,6)   =   -6D0*EIZ/(ELENGTH**2D0)
CC      STIFFBEAM(8,6)   =   -6*EIY/(LENGTH**2)
      STIFFBEAM(12,6)  =   2D0*EIZ/ELENGTH
      
      STIFFBEAM(12,8)  =  -6D0*EIZ/(ELENGTH**2D0)
CC      STIFFBEAM(12,8)  =  -6*EIY/(LENGTH**2)
      
      STIFFBEAM(11,9)  =   6D0*EIY/(ELENGTH**2D0)
       
      DO I = 1,12
          DO J = 1,12
              IF(I.NE.J) THEN
               STIFFBEAM(I,J) = STIFFBEAM(J,I)
              ENDIF
          ENDDO
      ENDDO    
          
      STIFFBEAM(3:6,1:12) = 0.0D0
      STIFFBEAM(9:12,1:12) = 0.0D0
      
      STIFFBEAM(1:12,3:6) = 0.0D0
      STIFFBEAM(1:12,9:12) = 0.0D0
      
             
C      STIFFBEAM(1:12,1:12) = STIFFBEAM(1:12,1:12)*10**(20)
      
C     TRANSFORMATION MATRIX 
	
	!LOCAL VECTOR
      VR = 0.0D0
	DO I = 1,3
	VR(I) = XYZMESH(2,I) - XYZMESH(1,I) 
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)
	IF(ELN.EQ.0) VR(1:3) = [1.0D0,0.0D0,0.0D0]
	CALL FMVEVR (VR,VS,VT)
      CALL ROMBAC (VR,VS,VT,RANG)
      

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
      
C     TRANSFORM OF STIFFNESS MATRIX 
      !Ke = T'*K*T
      
      STIFFBEAM1 = MATMUL(STIFFBEAM,TRANS)
      
      STIFFBEAM = MATMUL(TRANSPOSE(TRANS),STIFFBEAM1)
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================         
      SUBROUTINE FICBEAM_SUB(INEL,XYZM,PROPM,PROPG,STIFFBEAM,NODEMAT)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     --------------------------------------------------------------------------
C     CALCULATE THE STIFFNESS MATRIX OF FICTITIOUS BEAM THAT WILL BE SUBSTRACTED
C	--------------------------------------------------------------------------    
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      
      COMMON A(7000000),IA(6000000)

C     SAVE IEL
      COMMON /INELE/ IEL
      
      
      DIMENSION XYZM(25,3),NODEBEAM(4,2),XYZMESH(2,3),NODEMAT(2)
      DIMENSION PROPM(1),PROPG(1),STIFFBEAM(12,12)
      DIMENSION IDBEAM(5,6),IEDOFB(12)
     
      IF(INEL.EQ.1) THEN  !SIDE NUMBER IN ONE SUPER ELEMENT
          NOD1  =  1
          NOD2  =  2
      ELSEIF(INEL.EQ.2) THEN    
          NOD1  =  2
          NOD2  =  3
      ELSEIF(INEL.EQ.3) THEN    
          NOD1  =  3
          NOD2  =  4  
      ELSEIF(INEL.EQ.4) THEN    
          NOD1  =  4
          NOD2  =  1
      ELSEIF(INEL.EQ.5) THEN              
          NOD1  =  2
          NOD2  =  1
      ELSEIF(INEL.EQ.6) THEN    
          NOD1  =  3
          NOD2  =  2
      ELSEIF(INEL.EQ.7) THEN    
          NOD1  =  4
          NOD2  =  3  
      ELSEIF(INEL.EQ.8) THEN    
          NOD1  =  1
          NOD2  =  4          
      ENDIF       
      
          NODEMAT(1) = NOD1 !SAVE THE NODE NUMBER INTO NODEMAT MATRIX    
          NODEMAT(2) = NOD2
      
          NODE1 = NOD1
          NODE2 = NOD2
          
          XYZMESH(1,1:3) = XYZM(NODE1,1:3)
          XYZMESH(2,1:3) = XYZM(NODE2,1:3)
          
          CALL MODULEBEAM_SUB(XYZMESH,PROPM,PROPG,STIFFBEAM) ! MODULE FOR CALCULATING STIFFNESS MATRIX OF FICITIOUS BEAM 
                   
          
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================        
      SUBROUTINE MODULEBEAM_SUB(XYZMESH,PROPM,PROPG,STIFFBEAM)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     CALCULATE THE STIFFNESS MATRIX OF FINE MESHED FICTITIOUS BEAM
C	------------------------------------------------------------------------
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      
      COMMON A(7000000),IA(6000000)
      
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV      

      COMMON /HEI/ HEIGHT ! HEIGHT OF EACH SUPER ELEMENT
      
      DIMENSION XYZMESH(2,3),STIFFBEAM(12,12),TRANS(12,12),STIFFBEAM1(12,12)
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION PROPM(1),PROPG(1)
     
      X1 = XYZMESH(1,1)
      Y1 = XYZMESH(1,2)
      Z1 = XYZMESH(1,3)

      X2 = XYZMESH(2,1)
      Y2 = XYZMESH(2,2)
      Z2 = XYZMESH(2,3)     
      
      ELENGTH = SQRT((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
      
      CALL CLEARMAT(STIFFBEAM,12,12)
      
      E  = PROPM(1)
      PO = PROPM(2)
      TH = PROPG(2)
      
      G = E/(2D0*(1D0+PO))
C      TH = HEIGHT

      AREA = TH*HEIGHT
      EIZ = E*(TH*(HEIGHT**3D0))/12D0 ! E*IZ
      EIY = E*((TH**3D0)*HEIGHT)/12D0 ! E*IY
      
      PM = (EIZ + EIY)/E  !POLAR MOMENT OF INERTIA
      
      STIFFBEAM(1,1)   =  E*AREA/ELENGTH
      STIFFBEAM(2,2)   =  12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(3,3)   =  12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(4,4)   =  G*PM/ELENGTH
      STIFFBEAM(5,5)   =  4D0*EIY/ELENGTH
      STIFFBEAM(6,6)   =  4D0*EIZ/ELENGTH
      
      STIFFBEAM(7,7)   =  E*AREA/ELENGTH
      STIFFBEAM(8,8)   =  12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(9,9)   =  12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(10,10) =  G*PM/ELENGTH
      STIFFBEAM(11,11) =  4D0*EIY/ELENGTH
      STIFFBEAM(12,12) =  4D0*EIZ/ELENGTH      
      
      STIFFBEAM(7,1)   =  -E*AREA/ELENGTH
      
      STIFFBEAM(6,2)   =   6D0*EIZ/(ELENGTH**2D0)
      STIFFBEAM(8,2)   =  -12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(12,2)  =   6D0*EIZ/(ELENGTH**2D0)
      
      STIFFBEAM(5,3)   =  -6D0*EIY/(ELENGTH**2D0)
      STIFFBEAM(9,3)   =  -12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(11,3)  =  -6D0*EIY/(ELENGTH**2D0)
      
      STIFFBEAM(10,4)  =  -G*PM/ELENGTH
      
      STIFFBEAM(9,5)   =   6D0*EIY/(ELENGTH**2D0)
      STIFFBEAM(11,5)  =   2D0*EIY/ELENGTH
      
      STIFFBEAM(8,6)   =   -6D0*EIZ/(ELENGTH**2D0)
CC      STIFFBEAM(8,6)   =   -6*EIY/(LENGTH**2)
      STIFFBEAM(12,6)  =   2D0*EIZ/ELENGTH
      
      STIFFBEAM(12,8)  =  -6D0*EIZ/(ELENGTH**2D0)
CC      STIFFBEAM(12,8)  =  -6*EIY/(LENGTH**2)
      
      STIFFBEAM(11,9)  =   6D0*EIY/(ELENGTH**2D0)
       
      DO I = 1,12
          DO J = 1,12
              IF(I.NE.J) THEN
               STIFFBEAM(I,J) = STIFFBEAM(J,I)
              ENDIF
          ENDDO
      ENDDO    
          
      STIFFBEAM(3:6,1:12) = 0.0D0
      STIFFBEAM(9:12,1:12) = 0.0D0
      
      STIFFBEAM(1:12,3:6) = 0.0D0
      STIFFBEAM(1:12,9:12) = 0.0D0
      
          
C      STIFFBEAM(1:12,1:12) = STIFFBEAM(1:12,1:12)*10**(20)
      
C     TRANSFORMATION MATRIX 
	
	!LOCAL VECTOR
      VR = 0.0D0
	DO I = 1,3
	VR(I) = XYZMESH(2,I) - XYZMESH(1,I) 
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)
	IF(ELN.EQ.0) VR(1:3) = [1.0D0,0.0D0,0.0D0]
	CALL FMVEVR (VR,VS,VT)
      CALL ROMBAC (VR,VS,VT,RANG)



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
      
C     TRANSFORM OF STIFFNESS MATRIX 
      !Ke = T'*K*T
      
      STIFFBEAM1 = MATMUL(STIFFBEAM,TRANS)
      
      STIFFBEAM = MATMUL(TRANSPOSE(TRANS),STIFFBEAM1)
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================          
      SUBROUTINE MEM_CONDEN(SKS,ELODAM1,ELODCON,R,NNEQ,NEQ,ICONN,ILO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     STATIC CONDENSTAION OF MATRIX
C     SKS : STIFFNESS MATRIX TO BE REDUCED
C     R : LOAD VECTOR
C     NEQ : TOTAL NUMBER OF EQUATION
C     LEQ : NUMBER OF EQUATION TO BE ELIMINATED      
C     ------------------------------------------------------------
      CHARACTER*3 OPER
         
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C     SAVE IEL
      COMMON /INELE/ IEL   
      COMMON /JNELE/ JEL
      
C     STIFFNESS MATRIX PART
      DIMENSION SKS(NNEQ,NNEQ)
      DIMENSION SKUU(24,24),SKUC(24,126),SKCC(126,126),SKCU(126,24),
     1          SKCCI(126,126), ! INVERSION OF SKCC MATRIX
     2          SKSC(24,24),  ! CONDENSED STIFFNESS MATRIX
     3          SKCIU(24,126),SKSM(24,24), 
     4          TB(126,24),T(150,24),TI(24,24),TP(126,1)
      
      DIMENSION ELODAM1(150),ELODC(126),ELODCON(24),
     1          SKUUI(24,24),SKEC(24),
     2          RP(24),SKUR(24)
      
      DIMENSION ICONN(4,16)

C     LOAD MATRIX PART       
      DIMENSION RELEM(NNEQ),R(NEQ),
     1          FCS(24), !CONDENSED LOAD MATRIX      
     2          FU(24),FC(126),FUC(24)
      
      
C     1. REARRANGE THE STIFFNESS MATRIX (SUPER NODES AND ELIMINATED NODES)
      
      NEQU = 24
      NEQC = 126
      
      CALL CLEARMAT(SKUU,24,24)
      CALL CLEARMAT(SKUC,24,126)
      CALL CLEARMAT(SKCC,126,126)
      CALL CLEARMAT(SKCU,126,24)
      CALL CLEARMAT(SKCCI,126,126)
      CALL CLEARMAT(SKSC,24,24)
      CALL CLEARMAT(SKCIU,24,126)
      CALL CLEARMAT(SKSM,24,24)
      CALL CLEARMAT(SKSC,24,24)
      CALL CLEARMAT(SKUUI,24,24)
      CALL CLEARMAT(ELODCON,24,1)

      
      SKUU(1:24,1:24)   = SKS(1:24,1:24)
      SKUC(1:24,1:126)  = SKS(1:24,25:150)
      SKCC(1:126,1:126) = SKS(25:150,25:150)
      SKCU(1:126,1:24)  = SKS(25:150,1:24)
      
      
C     2. INVERSE SKCC MATRIX
         
      CALL INVMATRIX(SKCC,SKCCI,126)
C     3. CONDENSATION OF STIFFNESS MATRIX  
C     =========================================================      
C     MULTIPLYING MATRICES USING 'DGEMM'
C     ---------------------------------------------------------
C     C = ALPHA*A*B + BETA*C   
C     A(M,K), B(K,N), C(M,N)   
C     "INITAILIZE C MATRIX BEFORE DOING CALL DGEMM"
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)     
C     =========================================================
      ALPHA = 1.0
      BETA  = 1.0
            
      M = NEQU
      N = NEQC
      K = NEQC
      CALL DGEMM('N','N',M,N,K,ALPHA,SKUC,M,SKCCI,K,BETA,SKCIU,M)
      
      M = NEQU
      N = NEQU
      K = NEQC     
      CALL DGEMM('N','N',M,N,K,ALPHA,SKCIU,M,SKCU,K,BETA,SKSM,M)    
      
      SKSC = SKUU - SKSM
      
      
C     NEW T- MATRIX ---------------------------------------------   
      
      DO I = 1,126
          II = I+24
          ELODC(I) = ELODAM1(II)
      ENDDO          
      
      CALL INVMATRIX(SKUU,SKUUI,24)
      
      !KUC*YC
      M = NEQU
      N = 1
      K = NEQC
      CALL DGEMM('N','N',M,N,K,ALPHA,SKUC,M,ELODC,K,BETA,SKEC,M) 
      
      M = NEQU
      N = 1
      K = NEQU
      CALL DGEMM('N','N',M,N,K,ALPHA,-SKUUI,M,SKEC,K,BETA,ELODCON,M)
      
      
      IF(ILO.EQ.0) RETURN    
      
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================        
      SUBROUTINE SHDELA_WALL (DR,DVOL,SLA,SLB)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATES ELASTIC RIGIDITIES FOR MEMBRANE-PLATE OR SHELL ELEM.
C	---------------------------------------------------------------
C     DR(8,8) = ELASTIC MEMBRANE,BENDING AND SHEAR RIGIDITIES
C     DVOL    = WT*DET
C     SLA,SLB = SHEAR LOCKING CONSTRAINT FACTORS
C	----------------------------------------------------------------
C	PLEASE NOT SLA,SLB ARE 1.0 MODIFIED ASSUMED STRAIN METHOD, 
C	NO FACTOR REQUIRED
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION DR(64)
C
      CALL CLEARA (DR,64)

C	WRITE(*,*) PR
C	PAUSE

      FAM = 1.0!1.8!2.82!1.8
      FAB = 1.0!1.05!1.14!1.6
      FAS = 1.0!1.5!1.8

      TH3  = TH*TH*TH/12.
      FACM = DVOL*TH*FAM
      FACB = DVOL*TH3*FAB
      FACS = DVOL*TH*FAS
C
      DR(1)  = FACM*A1
      DR(10) = FACM*A1
      DR(2)  = FACM*B1
      DR(9)  = FACM*B1
      DR(19) = FACM*C1
      DR(28) = FACB*A1
      DR(37) = FACB*A1
      DR(29) = FACB*B1
      DR(36) = FACB*B1
      DR(46) = FACB*C1
      DR(55) = FACS*SLA*D1
      DR(64) = FACS*SLB*D1

C
      RETURN
      END
C
C=====================================================================      
      SUBROUTINE FICBEAM_SUB2(INEL,XYZM,PROPM,PROPG,STIFFBEAM,NODEMAT,IRELEASE)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     --------------------------------------------------------------------------
C     CALCULATE THE STIFFNESS MATRIX OF FICTITIOUS BEAM THAT WILL BE SUBSTRACTED
C	--------------------------------------------------------------------------    
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      
      COMMON A(7000000),IA(6000000)

C     SAVE IEL
      COMMON /INELE/ IEL
      
      
      DIMENSION XYZM(25,3),NODEBEAM(4,2),XYZMESH(2,3),NODEMAT(2)
      DIMENSION PROPM(1),PROPG(1),STIFFBEAM(12,12)
      DIMENSION IDBEAM(5,6),IEDOFB(12)
     
      IF(INEL.EQ.1) THEN  !SIDE NUMBER IN ONE SUPER ELEMENT
          NOD1  =  1
          NOD2  =  2
      ELSEIF(INEL.EQ.2) THEN    
          NOD1  =  2
          NOD2  =  3
      ELSEIF(INEL.EQ.3) THEN    
          NOD1  =  3
          NOD2  =  4  
      ELSEIF(INEL.EQ.4) THEN    
          NOD1  =  4
          NOD2  =  1
      ELSEIF(INEL.EQ.5) THEN              
          NOD1  =  2
          NOD2  =  1
      ELSEIF(INEL.EQ.6) THEN    
          NOD1  =  3
          NOD2  =  2
      ELSEIF(INEL.EQ.7) THEN    
          NOD1  =  4
          NOD2  =  3  
      ELSEIF(INEL.EQ.8) THEN    
          NOD1  =  1
          NOD2  =  4          
      ENDIF    
      
          NODEMAT(1) = NOD1 !SAVE THE NODE NUMBER INTO NODEMAT MATRIX    
          NODEMAT(2) = NOD2
      
          NODE1 = NOD1
          NODE2 = NOD2
          
          XYZMESH(1,1:3) = XYZM(NODE1,1:3)
          XYZMESH(2,1:3) = XYZM(NODE2,1:3)
          
          CALL MODULEBEAM_SUB2(XYZMESH,PROPM,PROPG,STIFFBEAM,IRELEASE) ! MODULE FOR CALCULATING STIFFNESS MATRIX OF FICITIOUS BEAM 
                   
          
      RETURN

      END
C	=======================================================================
C	======================================================================= 
C	=======================================================================        
      SUBROUTINE MODULEBEAM_SUB2(XYZMESH,PROPM,PROPG,STIFFBEAM,IRELEASE)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     CALCULATE THE STIFFNESS MATRIX OF FINE MESHED FICTITIOUS BEAM
C	------------------------------------------------------------------------
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      
      COMMON A(7000000),IA(6000000)
      
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV      

      COMMON /HEI/ HEIGHT ! HEIGHT OF EACH SUPER ELEMENT
      
      DIMENSION XYZMESH(2,3),STIFFBEAM(12,12),TRANS(12,12),STIFFBEAM1(12,12)
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION PROPM(1),PROPG(1)
      
      DIMENSION TRANH(12,12),LREAS(12) !FOR BEAM RELEASE
     
      X1 = XYZMESH(1,1)
      Y1 = XYZMESH(1,2)
      Z1 = XYZMESH(1,3)

      X2 = XYZMESH(2,1)
      Y2 = XYZMESH(2,2)
      Z2 = XYZMESH(2,3)     
      
      ELENGTH = SQRT((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
      
      CALL CLEARMAT(STIFFBEAM,12,12)
      
      E  = PROPM(1)
      PO = PROPM(2)
      TH = PROPG(2)
      
      G = E/(2D0*(1D0+PO))
      !!TH = HEIGHT
 
      AREA = TH*HEIGHT
      EIZ = E*(TH*(HEIGHT**3D0))/12D0 ! E*IZ
      EIY = E*((TH**3D0)*HEIGHT)/12D0 ! E*IY
      
      PM = (EIZ + EIY)/E  !POLAR MOMENT OF INERTIA
      
      STIFFBEAM(1,1)   =  E*AREA/ELENGTH
      STIFFBEAM(2,2)   =  12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(3,3)   =  12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(4,4)   =  G*PM/ELENGTH
      STIFFBEAM(5,5)   =  4D0*EIY/ELENGTH
      STIFFBEAM(6,6)   =  4D0*EIZ/ELENGTH
      
      STIFFBEAM(7,7)   =  E*AREA/ELENGTH
      STIFFBEAM(8,8)   =  12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(9,9)   =  12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(10,10) =  G*PM/ELENGTH
      STIFFBEAM(11,11) =  4D0*EIY/ELENGTH
      STIFFBEAM(12,12) =  4D0*EIZ/ELENGTH      
      
      STIFFBEAM(7,1)   =  -E*AREA/ELENGTH
      
      STIFFBEAM(6,2)   =   6D0*EIZ/(ELENGTH**2D0)
      STIFFBEAM(8,2)   =  -12D0*EIZ/(ELENGTH**3D0)
      STIFFBEAM(12,2)  =   6D0*EIZ/(ELENGTH**2D0)
      
      STIFFBEAM(5,3)   =  -6D0*EIY/(ELENGTH**2D0)
      STIFFBEAM(9,3)   =  -12D0*EIY/(ELENGTH**3D0)
      STIFFBEAM(11,3)  =  -6D0*EIY/(ELENGTH**2D0)
      
      STIFFBEAM(10,4)  =  -G*PM/ELENGTH
      
      STIFFBEAM(9,5)   =   6D0*EIY/(ELENGTH**2D0)
      STIFFBEAM(11,5)  =   2D0*EIY/ELENGTH
      
      STIFFBEAM(8,6)   =   -6D0*EIZ/(ELENGTH**2D0)
CC      STIFFBEAM(8,6)   =   -6*EIY/(LENGTH**2)
      STIFFBEAM(12,6)  =   2D0*EIZ/ELENGTH
      
      STIFFBEAM(12,8)  =  -6D0*EIZ/(ELENGTH**2D0)
CC      STIFFBEAM(12,8)  =  -6*EIY/(LENGTH**2)
      
      STIFFBEAM(11,9)  =   6D0*EIY/(ELENGTH**2D0)

      DO I = 1,12
          DO J = 1,12
              IF(I.NE.J) THEN
               STIFFBEAM(I,J) = STIFFBEAM(J,I)
              ENDIF
          ENDDO
      ENDDO    
      
      
      IF(IRELEASE.EQ.1) THEN
          
          STIFFBEAM(3:5,1:12) = 0.0D0
          STIFFBEAM(1:12,3:5) = 0.0D0

          STIFFBEAM(9:12,1:12) = 0.0D0
          STIFFBEAM(1:12,9:12) = 0.0D0  
          
          !STIFFBEAM(6,1:12) = 0.0D0
          

      ELSEIF(IRELEASE.EQ.2) THEN      
          STIFFBEAM(3:6,1:12) = 0.0D0
          STIFFBEAM(1:12,3:6) = 0.0D0
          
          STIFFBEAM(9:11,1:12) = 0.0D0
          STIFFBEAM(1:12,9:11) = 0.0D0
        
          !STIFFBEAM(12,1:12) = 0.0D0      
      ELSEIF(IRELEASE.EQ.3) THEN
          STIFFBEAM(3:5,1:12) = 0.0D0
          STIFFBEAM(1:12,3:5) = 0.0D0

          STIFFBEAM(9:11,1:12) = 0.0D0
          STIFFBEAM(1:12,9:11) = 0.0D0
      ELSEIF(IRELEASE.EQ.4) THEN ! CASE THAT ONLY THETA T STIFFNESS AT ONE NODE
          STIFFBEAM(1:11,1:12) = 0.0D0
          STIFFBEAM(1:12,1:11) = 0.0D0
      ELSE          
          STIFFBEAM(1:6,1:12) = 0.0D0
          STIFFBEAM(1:12,1:6) = 0.0D0
          
          STIFFBEAM(7:11,1:12) = 0.0D0
          STIFFBEAM(1:12,7:11) = 0.0D0
      ENDIF          
      
       
              
C     TRANSFORMATION MATRIX 
	
	!LOCAL VECTOR
      VR = 0.0D0
	DO I = 1,3
	VR(I) = XYZMESH(2,I) - XYZMESH(1,I) 
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)
	IF(ELN.EQ.0) VR(1:3) = [1.0D0,0.0D0,0.0D0]
	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

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
           

C     TRANSFORM OF STIFFNESS MATRIX 
      !Ke = T'*K*T
      
      STIFFBEAM1 = MATMUL(STIFFBEAM,TRANS)
      
      STIFFBEAM = MATMUL(TRANSPOSE(TRANS),STIFFBEAM1)
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================            