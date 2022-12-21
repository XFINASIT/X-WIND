C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE TENSERH (LEST,NFLINK,NDFREE,XYZ,JEG,IEGT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C     -----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /GiDEle/ LGID 
	
      COMMON A(9000000),IA(9000000)
      
C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN

C	COMMON BLOCK FOR TEND PRAMIN  NOV2010 	
	COMMON /BOX/ NBOX,NLEVL      	
	

      DIMENSION XYZ(NSN,1)
      
      DIMENSION LEST(1),NFLINK(30,30),NDFREE(30,30)

      DIMENSION COORT(3),COORN(3),TOFFS(3),HF(27),NODE(27),VR(3),VXF(3),VXB(3)

      DIMENSION COORN1(3),COORN2(3),TOFFS1(3),TOFFS2(3),NODE1(27),NODE2(27),HF1(27),HF2(27)
      
      ALLOCATABLE TPT(:,:),NST(:),TSTP(:,:),SLT(:)
      ALLOCATABLE MCONT(:),COORD(:,:),NODEX(:)

 
C     BOX SEARCH ---------------------------------------
C     FOR BOX SEARCH ADDED BY PRAMIN NOV2010     
      ALLOCATABLE BOUD(:,:),COORTT(:,:),MCONTT(:,:)  
      ALLOCATABLE EDATA(:,:),ICOUT(:,:),MAPGP(:,:) 
      ALLOCATABLE MAPTEN(:),TENELM(:,:)
C     --------------------------------------------------   

      SELECTCASE(ISTYP)
      CASE(1) !NODAL TENDON
        ITDYP = 0
      CASE(2) !FRAME TENDON
        ITDYP = 5
      CASE(3) !SHELL TENDON
        ITDYP = 9
      CASE(4) !SOLID TENDON
        ITDYP = 10
      ENDSELECT  


      NFIL = 3500	
	CALL SEQOPEN(NFIL) !TEMPORARY FILE FOR ELEMENT GROUP DATA USING FOR SEARCHING
	
      IFIL = NFIL + 2*IEGT - 1
      IFILJ= NFIL + 2*IEGT - 0
	CALL DIROPEN(IFIL ,1218)   !518 BYTE PER ACCESS
	CALL DIROPEN(IFILJ,1218)   !518 BYTE PER ACCESS
	
      CALL DEFNINT('#TED',KTED,1,20)

      ITDTYP = ITYPE
      ISTDYP = ISTYP
      CALL INTFILL('#TED',ITYPE,1,1 ,1)
      CALL INTFILL('#TED',ISTYP,1,2 ,1)
      CALL INTFILL('#TED',NLOPT,1,3 ,1)
      CALL INTFILL('#TED',MTMOD,1,4 ,1)
      CALL INTFILL('#TED',NELE ,1,5 ,1)
      CALL INTFILL('#TED',NMPS ,1,6 ,1)
      CALL INTFILL('#TED',NGPS ,1,7 ,1)
      CALL INTFILL('#TED',NOPS ,1,8 ,1)
      CALL INTFILL('#TED',NHIGE,1,9 ,1)
      CALL INTFILL('#TED',NMV  ,1,10,1)
      CALL INTFILL('#TED',NMP  ,1,11,1)
      CALL INTFILL('#TED',NGP  ,1,12,1)
      CALL INTFILL('#TED',NNM  ,1,13,1)
      CALL INTFILL('#TED',NGR  ,1,14,1)
      CALL INTFILL('#TED',NGS  ,1,15,1)
      CALL INTFILL('#TED',NGT  ,1,16,1)
      


      NELEGT = 0
      REWIND(NFIL)
C     ----------------------------------------------
	DO IEG = 1,JEG-1
	KEG = IEG
	NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     
      CALL MOVLEV (2)
      
      IF(ITYPE.EQ.ITDYP) THEN
      NELEGT = NELEGT + 1
      
      WRITE(NFIL) IEG,ITYPE,ISTYP,NELE,NCO,NEF,NEX,NNM
      
      ALLOCATE(NODEX(NEX),MCONT(NEF),COORD(NCO,NNM))
C     ------------------------------	
	DO IEL = 1,NELE
      CALL TFILLDATA(IA(LCN),A(LCO),IA(LEX),MCONT,COORD,NODEX,IEL,NEF,NCO,NNM,NEX)
	
      SELECTCASE(ITYPE)
      CASE(5)  !FRAME
        NPR = 1
        CALL SLEN_FRAME(IEG,IEL,SMAX) !MAXIMUM ALLOABLE OFFSET DISTANCE
      CASE(9)  !SHELL
        NPR = 2
        CALL SLEN_SHELL(IA(LGS),A(LGP),NGP,IEG,IEL,SMAX)  !MAXIMUM ALLOABLE OFFSET DISTANCE
      CASE(10) !SOLID
        NPR = 3
        SMAX = 0.0D0
      ENDSELECT
      
      WRITE(NFIL) IEL,NPR,SMAX,NODEX(1:NEX),MCONT(1:NEF),COORD(1:NCO,1:NNM)
      ENDDO
C     ------------------------------	     
      DEALLOCATE(MCONT,COORD,NODEX)
      ENDIF
          
      ENDDO
C     ----------------------------------------------

	      
 
C	NTEND = NUMBER OF TENDON
	READ (ITI,*) 
	READ (ITI,*) NTEND
	
C     BACK UP HERE BECAUSE IT WILL RESET TO ZERO AFTER CALL MOVLEV
	NTENDO = NTEND
	
      CALL DEFNINT('LTDN',KLTDN,10,NTEND)
      
      MTNOD = 0
      MTSEG = 0
      MTREC = 0
      
      
C     =====================================================================
C                       BOX DIVISION FOR SOLID ELEMENT LOOP           
C     =====================================================================  
       
C     SELECT TENDON-SOLID ELEMENT SEACHING PROGRAMS (1 = OLD, 2 = NEW) 
      ISEARCH = 2       
                  
C     FOR BOX SEARCH ADDED BY PRAMIN NOV2010 
C     CREATE BOXS 
      NLEVL = 3  ! The Optimized Value  
      ALLOCATE(BOUD(NLEVL,6))
      BOUD = 0.0D0 
      CALL BOXDATA(XYZ,BOUD)
      
C     CALL SOLID ELEMENT COORDINATES AND CONNECTIVITIES
      NGG = NEG-1
	DO 100 IEG = 1,NEG-1
	KEG = IEG
	NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)  
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)         
      CALL MOVLEV (2) ! CHANGE ALL DATA TO BE OTHER ELEMENTS

      IF(ITYPE.NE.10) GOTO 100      

C     FIND MAX AND MIN LIMITS OF EACH ELEMENT
      ALLOCATE(EDATA(NELE,6))
      EDATA = 0.0D0 
C     NUMBER OF ELEMENTS IN EACH BOX      
      ALLOCATE(ICOUT(NBOX,1))    
      ICOUT = 0
C     ELEMENT NUMBERS BEING IN EACH BOX      
      ALLOCATE(MAPGP(NBOX,NELE)) 
      MAPGP = 0.0D0 
C     BACKUP COORDINATES AND CONNECTIVITIES      
      ALLOCATE(COORTT(NCO*NNM*NELE,NEG-1)) 
      COORTT = 0.0D0     
      NSELE = NELE
      ALLOCATE(MCONTT(NEF*NELE,NEG-1))
      MCONTT = 0
C     ELEMENT BOX MAPPING     
      CALL ELEMENT_BOX_MAP(IA(LCN),A(LCO),BOUD,EDATA,ICOUT,MAPGP,COORTT,MCONTT,NGG,IEG) 
      MAXIC = MAXVAL(ICOUT(1:NBOX,1))       
       
100   CONTINUE
      
C     ---------------------------------------------------------------------
     
      
150   CONTINUE
      IIO = 0     
C     =====================================================================
C                                 TENDON LOOP 
C     =====================================================================
      DO 9000 ITEND = 1,NTENDO

C	READ & STORE TENDON NAME
	CALL TENDONMAP(ITEND,ITI,IEGT)
	READ(ITI,*) MSPAN,IMAT,IGEO !NUMBER OF SPAN, MATERIAL No., GEOMETRY No.

C     ----------------------------------------------------------------
      LFIRSTC = 0
      LFIRSTE = 0
      IFOUNDL = 0
      MNOD = 0
      
C     ----------------------------------------------------------------       
	DO 8500 ISPAN = 1,MSPAN
		      
	CALL FREBUF
      CALL FREINT('P',NTPOIN,1) !INPUT NUMBER OF TENDON POINT
      IF(ITDYP.EQ.0) THEN !FOR TENDON ON NODE
        CALL FREINT('S',NSNOD,1) !INPUT NUMBER OF SPAN NODE
        ALLOCATE(TPT(3,NTPOIN),NST(NSNOD),TSTP(2,NTPOIN),SLT(NSNOD))
      ENDIF    
          
C     ------------------------------------------------------------
C     FOR BOX SEARCH ADDED BY PRAMIN NOV2010  
      ALLOCATE(MAPTEN(NTPOIN))
      MAPTEN = 0  
C     TENDON BOX SEARCH
      ALLOCATE(TENELM(NTPOIN,7))
      TENELM = 0.0D0        
C     ------------------------------------------------------------ POINT 

      DO 8000 ITPOIN = 1,NTPOIN

	CALL FREBUF
      CALL FREINT('N',ITP,1)	
      CALL FREREL('P',COORT(1),3) !READ TENDON POINT COORDINATE
      
C     CONTROL READ TENDON POINT VALUES 
      DECI = 0.0000000010D0
      CX = ABS(COORT(1))
      IF (CX <= DECI) THEN 
      COORT(1) = 0.0D0 
      ENDIF
      CY = ABS(COORT(2))
      IF (CY <= DECI) THEN
      COORT(2) = 0.0D0 
      ENDIF
      CZ = ABS(COORT(3))
      IF (CZ <= DECI) THEN
      COORT(3) = 0.0D0    
      ENDIF  
C     --------------------------------       

      SELECTCASE(ITDYP)
      CASE(0) !NODAL TENDON
        GOTO 201
      CASE DEFAULT !FRAME TENDON, SHELL TENDON, SOLID TENDON
        GOTO 202
      ENDSELECT 
  
C     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      
201   CONTINUE !FOR TENDON ON NODE
        TPT(1:3,ITPOIN) = COORT(1:3) !STORE TENDON POINT COORDINATE
 	  NDFREE(ITDTYP,ISTDYP) = 6
	  NFLINK(ITDTYP,ISTDYP) = 123456
      GOTO 8000 
C     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII         
  
  
C     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII          
202   CONTINUE !FOR TENDON ON ELEMENT        

      XT = COORT(1) ; YT = COORT(2) ; ZT = COORT(3) 
      IF(LFIRSTC.EQ.0) THEN
        LFIRSTC = 1
      ELSE
        DIFF = SQRT( (XTO-XT)**2.0 + (YTO-YT)**2.0 + (ZTO-ZT)**2.0 )
C        IF(DIFF.LT.0.5.AND.IFOUNDL.EQ.1) GOTO 8000
        IF(XTO.EQ.XT.AND.YTO.EQ.YT.AND.ZTO.EQ.ZT.AND.IFOUNDL.EQ.1) 
	1GOTO 8000 !PREVENT TENDON NODE FALLING INTO SAME POSITION
      ENDIF
      XTO = XT ; YTO = YT ; ZTO = ZT
      
C     ---------------------------------------------- 
C     PUT TENDONS IN TO BOXES           
      CALL TENDON_BOX(BOUD,MAPTEN,NTPOIN,XT,YT,ZT,ITPOIN)                  
      Xnod = XT 
      Ynod = YT
      Znod = ZT                    
C     ----------------------------------------------         
      
     
      REWIND(NFIL)            
C     ----------------------------------------------
	DO 1000 IELEGT = 1,NELEGT
	
      READ(NFIL) IEG,ITYPE,ISTYP,NELE,NCO,NEF,NEX,NNM
      
      ALLOCATE(NODEX(NEX),MCONT(NEF),COORD(NCO,NNM))
      
C     ----------------------------------------------
C     FOR BOX SEARCH ADDED BY PRAMIN NOV2010 
C     SELECT TENDON-SOLID ELEMENT SEACHING PROGRAM
      IF (ITYPE.NE.10) GOTO 400 
      IF (ITYPE.EQ.10) THEN      
      IF(ISEARCH.EQ.1) GOTO 400   ! FOR OLD APPROACH
      IF(ISEARCH.EQ.2) GOTO 500   ! FOR NEW APPORACH 
      ENDIF     
C     ----------------------------------------------        
      
C     -----------------------------------------------------------------------	
C     -----------------------------------------------------------------------	
400	CONTINUE

	DO 450 IELE = 1,NELE
	
      READ(NFIL) IEL,NPR,SMAX,NODEX(1:NEX),MCONT(1:NEF),COORD(1:NCO,1:NNM)
      CALL TPSEARCH(MCONT,COORD,NODEX,XT,YT,ZT,COORN,TOFFS,NODE,HF,IFOUND,NCO,NEF,NEX,NNM,NPR,SMAX)
      
C     ============================    
      IF(IFOUND.EQ.1) THEN !FOUND THE CORRESPONDING ELEMENT

      IF(LFIRSTE.EQ.0) THEN
        LFIRSTE = 1
      ELSE
C        IF(IELO.EQ.IEL.AND.IFOUNDL.EQ.1) GOTO 8000 !PREVENT TENDON NODE FALLING INTO SAME ELEMENT
      ENDIF
      IELO = IEL
      
      MNOD = MNOD + 1      
      
      WRITE(IFIL,REC=MTREC+MNOD) IEG,IEL,ITYPE,ISTYP,NNM,COORN(1:3),XT,YT,ZT,NODE(1:NNM),TOFFS(1:3),HF(1:NNM)

 	  NDFREE(ITDTYP,ISTDYP) = NDFREE(ITYPE,ISTYP)
	  NFLINK(ITDTYP,ISTDYP) = NFLINK(ITYPE,ISTYP)
	
	IFOUNDL = 1
	
	DEALLOCATE(MCONT,COORD,NODEX)  
      GOTO 8000
      ENDIF
C     ============================ 
450	CONTINUE
C     ------------------------------
      DEALLOCATE(MCONT,COORD,NODEX)
      GOTO 1000 
C     -----------------------------------------------------------------------	
C     -----------------------------------------------------------------------	


C     -----------------------------------------------------------------------	
C     -----------------------------------------------------------------------	
500	CONTINUE
C     TENDON POINT SEARCH   
C     ------------------------------  
C     CALL THE BOX NUMBER IN WHICH THE TENDON NODE IS     
      IBOX  = MAPTEN(ITPOIN)      
      IF (IBOX == 0) IBOX = 1                
      NSEL = ICOUT(IBOX,1)                      
C     Element Loop 
      IFOUND = 0
      ICHECK = 0
      DO 550 IE = 1,NSEL       
      IEL = MAPGP(IBOX,IE)
      XEmin = EDATA(IEL,1)
      XEmax = EDATA(IEL,2)   
      YEmin = EDATA(IEL,3)
      YEmax = EDATA(IEL,4) 
      ZEmin = EDATA(IEL,5)
      ZEmax = EDATA(IEL,6) 
      IF ((Xnod >= XEmin).AND.(Xnod <= XEmax)) THEN                   
      IF ((Ynod >= YEmin).AND.(Ynod <= YEmax)) THEN               
      IF ((Znod >= ZEmin).AND.(Znod <= ZEmax)) THEN          
C     CHECKING THE POINTS CORRESPONDING WITH ALL CONDITIONS                      
      CALL TENDON_SEARCH(IA(LCN),A(LCO),IA(LEX),XT,YT,ZT,COORN,NODE,NCO,NEF,NNM,
     1                         TENELM,ITPOIN,NTPOIN,NN,NELE,IEL,COORTT,NSELE,NEX,
     1                         MCONTT,NGG,IG,NPR,TOFFS,HF)
      ICHECK = INT(TENELM(ITPOIN,1))
      IF (ICHECK.NE.0) THEN
      IFOUND = 1 
      ENDIF               
C     ------------------------------         
      IF(IFOUND.EQ.1) THEN !FOUND THE CORRESPONDING ELEMENT

      IF(LFIRSTE.EQ.0) THEN
      LFIRSTE = 1
      ELSE
C     IF(IELO.EQ.IEL.AND.IFOUNDL.EQ.1) GOTO 8000 !PREVENT TENDON NODE FALLING INTO SAME ELEMENT
      ENDIF      
      IELO = IEL
             
      MNOD = MNOD + 1      
        
      WRITE(IFIL,REC=MTREC+MNOD) IEG,IEL,ITYPE,ISTYP,NNM,COORN(1:3),XT,YT,ZT,NODE(1:NNM),TOFFS(1:3),HF(1:NNM)
	
      IIO = IIO + 1
C      WRITE(*,15) IIO,INT(TENELM(ITPOIN,1)),TENELM(ITPOIN,2:7)
C15    FORMAT(I5,1X,I5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5)
     

 	NDFREE(ITDTYP,ISTDYP) = NDFREE(ITYPE,ISTYP)
	NFLINK(ITDTYP,ISTDYP) = NFLINK(ITYPE,ISTYP)
    	
	IFOUNDL = 1    	
    	
	DEALLOCATE(MCONT,COORD,NODEX) 
      GOTO 8000
      ENDIF
C     ------------------------------ 
      ENDIF 
      ENDIF
      ENDIF  
      
                   
C     ------------------------------                    
550   CONTINUE
C     ------------------------------
      DEALLOCATE(MCONT,COORD,NODEX) 
      GOTO 1000
C     -----------------------------------------------------------------------	
C     -----------------------------------------------------------------------	  


C     ------------------------------
1000	CONTINUE
C     ----------------------------------------------
      GOTO 8000
C     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    


8000	CONTINUE
      DEALLOCATE(MAPTEN)
      DEALLOCATE(TENELM)
C     ------------------------------------------------------------

      SELECTCASE(ITDYP)
      CASE(0)
	  READ(ITI,*) NST(1:NSNOD)    !READ SPAN NODES
        GOTO 8400
      CASE(5,9,10)
	  READ(ITI,*)                 !READ NOTHING
        GOTO 8500
      ENDSELECT
      
C     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      
8400  CONTINUE !FOR TENDON ON NODE

	TELN = 0.0D0 !LENGTH ALONG THE TENDON AT EACH TENDON POINT
	TSTP(1,1) = 0.0D0 !LENTGH AT THE BIGINING OF TENDON EQUAL TO ZERO
      DO ITPOIN = 1,NTPOIN-1
	    DO I = 1,3
            VR(I) = TPT(I,ITPOIN+1) - TPT(I,ITPOIN)
	    ENDDO
	    CALL SCALEN(VR,VR,SELN,3) !GET SEGMENT LENGTH HERE
	    TELN = TELN + SELN
	    TSTP(1,ITPOIN+1) = TELN
	ENDDO

	PELN = 0.0D0 !LENGTH ALONG THE SPAN NODE AT EACH TENDON POINT
	SLT(1) = 0.0D0 !LENTGH AT THE BIGINING OF TENDON EQUAL TO ZERO
	DO ISNOD = 1,NSNOD-1
	N1 = NST(ISNOD)
	N2 = NST(ISNOD+1)
	VR(1:3) = 0.0D0
	DO ISC = 1,NSC
	C1 = XYZ(N1,ISC)  !GETTING HERE NODAL COORDINATE
	C2 = XYZ(N2,ISC)  !GETTING HERE NODAL COORDINATE
	VR(ISC) = C2 - C1
	ENDDO
	CALL SCALEN(VR,VR,SELN,3)            !GET LENGTH OF NODE SPACING HERE
	PELN = PELN + SELN
	SLT(ISNOD+1) = PELN
	ENDDO	

	FACTOR = TELN/PELN
	DO 8450 ISNOD = 1,NSNOD
	NNS = NST(ISNOD)
	SELN = SLT(ISNOD)*FACTOR
	TSTP(2,1:NTPOIN) = TPT(1,1:NTPOIN)   !X
	CALL INTERPOL(SELN,VALX,TSTP,NTPOIN,IV)
	TSTP(2,1:NTPOIN) = TPT(2,1:NTPOIN)   !Y
	CALL INTERPOL(SELN,VALY,TSTP,NTPOIN,IV)
	TSTP(2,1:NTPOIN) = TPT(3,1:NTPOIN)   !Z
	CALL INTERPOL(SELN,VALZ,TSTP,NTPOIN,IV)
	
      XT = VALX ; YT = VALY ; ZT = VALZ
      IF(LFIRSTC.EQ.0) THEN
        LFIRSTC = 1
      ELSE
        DIFF = SQRT( (XTO-XT)**2.0 + (YTO-YT)**2.0 + (ZTO-ZT)**2.0 )
C        IF(DIFF.LT.0.5) GOTO 8450
        IF(XTO.EQ.XT.AND.YTO.EQ.YT.AND.ZTO.EQ.ZT) GOTO 8450 !PREVENT TENDON NODE FALLING INTO SAME POSITION
      ENDIF
      XTO = XT ; YTO = YT ; ZTO = ZT

      IF(LFIRSTE.EQ.0) THEN
        LFIRSTE = 1
      ELSE
        IF(NNSO.EQ.NNS) GOTO 8450 !PREVENT TENDON NODE FALLING INTO SAME NODE
      ENDIF
      NNSO = NNS
 
	VR(1:3) = 0.0D0
	DO ISC = 1,NSC
	VR(ISC) = XYZ(NNS,ISC)  !GETTING HERE NODAL COORDINATE
	ENDDO         	
	OFFX = VALX - VR(1) ; OFFY = VALY - VR(2) ; OFFZ = VALZ - VR(3)    	
     	
      MNOD = MNOD + 1      
      WRITE(IFIL,REC=MTREC+MNOD) NNS,XT,YT,ZT,OFFX,OFFY,OFFZ
      
8450  CONTINUE
      
			
	DEALLOCATE(TPT,NST,TSTP,SLT)
      GOTO 8500
C     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII   

8500  CONTINUE
C     ----------------------------------------------------------------
      MSEG = MNOD - 1 !NUMBER OF TENDON SEGMENT
      IF(MSEG.LT.0) MSEG = 0
      
      SELECTCASE(ITDYP)
      CASE(0)	
      WRITE(ISO,800) ITEND,MSEG
	WRITE(ISO,810) 
      CASE(5)	
      WRITE(ISO,801) ITEND,MSEG
	WRITE(ISO,830)
      CASE(9)	
      WRITE(ISO,802) ITEND,MSEG
	WRITE(ISO,830)
      CASE(10)	
      WRITE(ISO,803) ITEND,MSEG
	WRITE(ISO,830)
	ENDSELECT      
      
      DO 8600 ISEG = 1,MSEG
      ISEGI = ISEG + 0
      ISEGJ = ISEG + 1
      SELECTCASE(ITDYP)
      CASE(0)
          READ(IFIL,REC=MTREC+ISEGI) NNS1,XT1,YT1,ZT1,OFFX1,OFFY1,OFFZ1
          READ(IFIL,REC=MTREC+ISEGJ) NNS2,XT2,YT2,ZT2,OFFX2,OFFY2,OFFZ2
      CASE(5,9,10) !SOLID
          READ(IFIL,REC=MTREC+ISEGI) IEG1,IEL1,ITYPE1,ISTYP1,NNM1,COORN1(1:3),XT1,YT1,ZT1,NODE1(1:NNM1),TOFFS1(1:3),HF1(1:NNM1)
          READ(IFIL,REC=MTREC+ISEGJ) IEG2,IEL2,ITYPE2,ISTYP2,NNM2,COORN2(1:3),XT2,YT2,ZT2,NODE2(1:NNM2),TOFFS2(1:3),HF2(1:NNM2)
      ENDSELECT
      
C	SEGMENT LENGTH AND VECTOR
	SELN = 0.0D0
	VR(1) = XT2 - XT1
	VR(2) = YT2 - YT1
	VR(3) = ZT2 - ZT1
	CALL SCALEN(VR,VR,SELN,3)  !GET LENGTH OF SEGMENT HERE

C	REFERENCE VECTOR FORWARD
	ANGF = 0.0D0 !Segment Angle Forward
	IF(ISEG.GT.1) THEN
	    COST = VR(1)*VXF(1)+VR(2)*VXF(2)+VR(3)*VXF(3)
	    IF(ABS(COST).GT.1.0D0) COST = 1.0D0*COST/ABS(COST)
	    ANGF = ACOS(COST)
	ENDIF
	VXF(1:3) = VR(1:3)

		
C	REFERENCE VECTOR BACKWARD
      ISEGBI = (MSEG + 1) - (ISEG - 1)
      ISEGBJ = (MSEG + 1) - (ISEG - 0)
      SELECTCASE(ITDYP)
      CASE(0)
          READ(IFIL,REC=MTREC+ISEGBI) NNS1,XT1,YT1,ZT1,OFFX1,OFFY1,OFFZ1
          READ(IFIL,REC=MTREC+ISEGBJ) NNS2,XT2,YT2,ZT2,OFFX2,OFFY2,OFFZ2
      CASE(5,9,10) !SOLID
          READ(IFIL,REC=MTREC+ISEGBI) IEG1,IEL1,ITYPE1,ISTYP1,NNM1,COORN1(1:3),XT1,YT1,ZT1,NODE1(1:NNM1),TOFFS1(1:3),HF1(1:NNM1)
          READ(IFIL,REC=MTREC+ISEGBJ) IEG2,IEL2,ITYPE2,ISTYP2,NNM2,COORN2(1:3),XT2,YT2,ZT2,NODE2(1:NNM2),TOFFS2(1:3),HF2(1:NNM2)
      ENDSELECT
      
C	SEGMENT LENGTH AND VECTOR
	SELNB = 0.0D0
	VR(1) = XT2 - XT1
	VR(2) = YT2 - YT1
	VR(3) = ZT2 - ZT1
	CALL SCALEN(VR,VR,SELNB,3)  !GET LENGTH OF SEGMENT HERE
	
	ANGB = 0.0D0 !Segment Angle Backward
	IF(ISEG.GT.1) THEN
	    COST = VR(1)*VXB(1)+VR(2)*VXB(2)+VR(3)*VXB(3)
	    IF(ABS(COST).GT.1.0D0) COST = 1.0D0*COST/ABS(COST)
	    ANGB = ACOS(COST)
	ENDIF
	VXB(1:3) = VR(1:3)


      SELECTCASE(ITDYP)
      CASE(0)
          READ(IFIL,REC=MTREC+ISEGI) NNS1,XT1,YT1,ZT1,OFFX1,OFFY1,OFFZ1
          READ(IFIL,REC=MTREC+ISEGJ) NNS2,XT2,YT2,ZT2,OFFX2,OFFY2,OFFZ2
      CASE(5,9,10) !SOLID
          READ(IFIL,REC=MTREC+ISEGI) IEG1,IEL1,ITYPE1,ISTYP1,NNM1,COORN1(1:3),XT1,YT1,ZT1,NODE1(1:NNM1),TOFFS1(1:3),HF1(1:NNM1)
          READ(IFIL,REC=MTREC+ISEGJ) IEG2,IEL2,ITYPE2,ISTYP2,NNM2,COORN2(1:3),XT2,YT2,ZT2,NODE2(1:NNM2),TOFFS2(1:3),HF2(1:NNM2)
      ENDSELECT
         
C     RECORD SEGMENT DATA TO FILE
      ISEGR = MNOD + ISEG
      SELECTCASE(ITDYP)
      CASE(0)
      WRITE(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB,
     1                      NNS1,XT1,YT1,ZT1,OFFX1,OFFY1,OFFZ1,
     2                      NNS2,XT2,YT2,ZT2,OFFX2,OFFY2,OFFZ2
     	WRITE(ISO,820) SELN,ANGF,ANGB,
     1               NNS1,XT1,YT1,ZT1,OFFX1,OFFY1,OFFZ1,
     2               NNS2,XT2,YT2,ZT2,OFFX2,OFFY2,OFFZ2
      NNM = 2
      CASE(5,9,10) !SOLID
      WRITE(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB,
     1                      IEG1,IEL1,ITYPE1,ISTYP1,NNM1,COORN1(1:3),XT1,YT1,ZT1,NODE1(1:NNM1),TOFFS1(1:3),HF1(1:NNM1),
     2                      IEG2,IEL2,ITYPE2,ISTYP2,NNM2,COORN2(1:3),XT2,YT2,ZT2,NODE2(1:NNM2),TOFFS2(1:3),HF2(1:NNM2)
	WRITE(ISO,840) SELN,ANGF,ANGB,
	1               XT1,YT1,ZT1,TOFFS1(1:3),
     2               XT2,YT2,ZT2,TOFFS2(1:3)
      NNM = NNM1+NNM2
      ENDSELECT
      
8600  CONTINUE

      CALL INTFILL('LTDN',MNOD ,1,ITEND,1)
      CALL INTFILL('LTDN',MSEG ,2,ITEND,1)
      CALL INTFILL('LTDN',IFIL ,3,ITEND,1)
      CALL INTFILL('LTDN',MTREC,4,ITEND,1)
      CALL INTFILL('LTDN',ITDYP,5,ITEND,1)
      CALL INTFILL('LTDN',IMAT ,6,ITEND,1)
      CALL INTFILL('LTDN',IGEO ,7,ITEND,1)
      CALL INTFILL('LTDN',IFILJ,8,ITEND,1)
      CALL INTFILL('LTDN',IEGT ,9,ITEND,1)
           
      CALL INTFILL('LTDN',MNOD ,1,ITEND,0)
      
      MTNOD = MTNOD + MNOD
      MTSEG = MTSEG + MSEG
      MTREC = MTREC + MNOD + MSEG

      
9000	CONTINUE
C     ---------------------------------------------------------------------

      CALL INTFILL('#TED',ITYPE,1,1 ,0)
      CALL INTFILL('#TED',ISTYP,1,2 ,0)
      CALL INTFILL('#TED',NLOPT,1,3 ,0)
      CALL INTFILL('#TED',MTMOD,1,4 ,0)
      CALL INTFILL('#TED',NELE ,1,5 ,0)
      CALL INTFILL('#TED',NMPS ,1,6 ,0)
      CALL INTFILL('#TED',NGPS ,1,7 ,0)
      CALL INTFILL('#TED',NOPS ,1,8 ,0)
      CALL INTFILL('#TED',NHIGE,1,9 ,0)
      CALL INTFILL('#TED',NMV  ,1,10,0)
      CALL INTFILL('#TED',NMP  ,1,11,0)
      CALL INTFILL('#TED',NGP  ,1,12,0)
C     CALL INTFILL('#TED',NNM  ,1,13,0)
      CALL INTFILL('#TED',NGR  ,1,14,0)
      CALL INTFILL('#TED',NGS  ,1,15,0)
      CALL INTFILL('#TED',NGT  ,1,16,0)
      
C     RECALL HERE 
      NTEND = NTENDO
      
      NELE = MTSEG
      NMPS = NTEND
      NGPS = NTEND
 
      CALL DELTINT('#TED')

      CLOSE(NFIL)

800	FORMAT(/5X,'Nodal TENDON No.',I5,2X,'NUMBER OF TENDON SEGMENT =',I5/)
801	FORMAT(/5X,'Frame TENDON No.',I5,2X,'NUMBER OF TENDON SEGMENT =',I5/)
802	FORMAT(/5X,'Shell TENDON No.',I5,2X,'NUMBER OF TENDON SEGMENT =',I5/)
803	FORMAT(/5X,'Solid TENDON No.',I5,2X,'NUMBER OF TENDON SEGMENT =',I5/)

810	FORMAT( 25X,'LENGTH',13X,'ANGLE F.',10X,'ANGLE B.',5X,'NODE-i',19X,'COORDINATE-i',35X,'OFFSET-i',
     1                                                  14X,'NODE-j',19X,'COORDINATE-j',35X,'OFFSET-j'/)
820	FORMAT( 20X,E14.5,5X,E14.5,5X,E14.5,5X,I4,3X,3E14.5,3X,3E14.5,5X,I4,3X,3E14.5,3X,3E14.5)
830	FORMAT( 25X,'LENGTH',13X,'ANGLE F.',10X,'ANGLE B.',26X,'COORDINATE-i',35X,'OFFSET-i',
     1                                                   40X,'COORDINATE-j',35X,'OFFSET-j'/)
840	FORMAT( 20X,E14.5,5X,E14.5,5X,E14.5,5X,3X,3E14.5,3X,3E14.5,5X,3X,3E14.5,3X,3E14.5)

      
      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE TENDATA(LTPROP) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN
           
      DIMENSION LTPROP(10,NTEND)
      
      DO ITEND = 1,NTEND
          CALL INTFILL('LTDN',LTPROP(1 ,ITEND),1 ,ITEND,0)
          CALL INTFILL('LTDN',LTPROP(2 ,ITEND),2 ,ITEND,0)
          CALL INTFILL('LTDN',LTPROP(3 ,ITEND),3 ,ITEND,0)
          CALL INTFILL('LTDN',LTPROP(4 ,ITEND),4 ,ITEND,0)
          CALL INTFILL('LTDN',LTPROP(5 ,ITEND),5 ,ITEND,0)
          CALL INTFILL('LTDN',LTPROP(6 ,ITEND),6 ,ITEND,0)
          CALL INTFILL('LTDN',LTPROP(7 ,ITEND),7 ,ITEND,0)
          CALL INTFILL('LTDN',LTPROP(8 ,ITEND),8 ,ITEND,0)  
          CALL INTFILL('LTDN',LTPROP(9 ,ITEND),9 ,ITEND,0)  
      ENDDO 
      CALL DELTINT('LTDN')
      
      
      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================

	SUBROUTINE TENDNODE(LM,LTPROP,LGID,MSET,ISET,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN

      DIMENSION LTPROP(10,1),LM(NEF,1),MSET(1),ISET(1),LGID(1)
      DIMENSION COORN1(3),COORN2(3),TOFFS1(3),TOFFS2(3),NODE1(27),NODE2(27),HF1(27),HF2(27)

      IEL = 0
      DO ITEND = 1,NTEND
      
	NNOD = LTPROP(1,ITEND)	!Number of tendon node
	NSEG = LTPROP(2,ITEND)	!Number of tendon segment
	IFIL = LTPROP(3,ITEND)	!file number
	MTREC= LTPROP(4,ITEND)	!last record number (just before this tendon)
	ITDYP= LTPROP(5,ITEND)	!ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)
	IMAT = LTPROP(6,ITEND)	!MAT No.
	IGEO = LTPROP(7,ITEND)	!GEO No.
	IFILJ= LTPROP(8,ITEND)	!file number (store force after jacking)

      DO ISEG = 1,NSEG
      
      IEL = IEL + 1  
      
      LGID(IEL) = 0
      MSET(IEL) = ITEND
      ISET(IEL) = ITEND
      
	ISEGR = NNOD + ISEG
      SELECTCASE(ITDYP)
      CASE(0)
      READ(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB,
     1                      NNS1,XT1,YT1,ZT1,OFFX1,OFFY1,OFFZ1,
     2                      NNS2,XT2,YT2,ZT2,OFFX2,OFFY2,OFFZ2
      LM(1:2,IEL) = [NNS1,NNS2]
      CASE(5,9,10) !SOLID
      READ(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB,
     1                      IEG1,IEL1,ITYPE1,ISTYP1,NNM1,COORN1(1:3),XT1,YT1,ZT1,NODE1(1:NNM1),TOFFS1(1:3),HF1(1:NNM1),
     2                      IEG2,IEL2,ITYPE2,ISTYP2,NNM2,COORN2(1:3),XT2,YT2,ZT2,NODE2(1:NNM2),TOFFS2(1:3),HF2(1:NNM2)
      
      LM(1:NNM1,IEL) = NODE1(1:NNM1)
      LM(NNM1+1:NNM1+NNM2,IEL) = NODE2(1:NNM2)
      ENDSELECT
      	
      ENDDO
      
      ENDDO
      


	RETURN
	END


C	=====================================================================
C	=====================================================================
C	=====================================================================
      

	SUBROUTINE TPSEARCH(MCONT,COORD,NODEX,XT,YT,ZT,POSN,OFFS,NODE,H,IFOUND,NCO,NEF,NEX,NNM,NPR,SMAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      DIMENSION MCONT(NEF),COORD(NCO,NNM),NODEX(NEX)
	DIMENSION POSN(1),OFFS(1),NODE(1),COORP(3),H(NNM)

C      CALL TPITER(MCONT(1,IEL),COORD(1,IEL),XT,YT,ZT,COORN,NODE,IFOUND,NNM,IEL)
      
	COORP(1:3) = [XT,YT,ZT]
	CALL TSEARCH(POSN,COORD,NODEX,COORP,NPR,NCO,NNM,NEX,SLEN)
      NODE(1:NNM) = MCONT(1:NNM)

C     CALCULATE OFFSET
	CALL TDOFFSET(POSN,COORD,COORP,OFFS,H,NPR,NCO,NNM,NODEX)
	
      TOL   = 1.001
	ITEST = 0
	DO IPR = 1,NPR
	  IF(ABS(POSN(IPR)).GT.TOL) ITEST = 1 !NOT YET FOUND
	ENDDO

      SELECTCASE(NPR)
      CASE(1) !FRAME
	  IF(SLEN.GT.SMAX) ITEST = 1 !NOT YET FOUND
      CASE(2) !SHELL
	  IF(SLEN.GT.SMAX) ITEST = 1 !NOT YET FOUND
      CASE(3) !SOLID
C	  NOTHING FOR SOLID
      ENDSELECT

	IFOUND = 1
	IF(ITEST.EQ.1) IFOUND = 0
	
	
	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE TSEARCH(POS,COORD,NODEX,COORP,NPR,NCO,NNM,NEX,SLEN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION COORD(NCO,NNM),COORP(NCO),NODEX(NEX)
	DIMENSION F(NPR),HH(NPR,NPR),HHI(NPR,NPR),POS(NPR)

	NTER = 10

C	ASSUMED THE INITIAL VALUE OF RI,SI,TI
	POS(1:NPR) = 0.0D0
	 
	DO 100 ITER = 1,NTER

	CALL PDERIV(POS,COORD,COORP,F,HH,SLEN,NPR,NCO,NNM,NODEX)

	CALL INVMATRIX(HH,HHI,NPR)

	DO I = 1,NPR
	DO J = 1,NPR
	POS(I) = POS(I) - HHI(I,J)*F(J)
	ENDDO
	ENDDO

100	CONTINUE


	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================


	SUBROUTINE TDOFFSET(POS,COORD,COORP,OFFS,H,NPR,NCO,NNM,NODEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION H(NNM),P(NPR,NNM),COORD(NCO,NNM),COORP(NCO),POS(NPR)
	DIMENSION OFFS(NCO),NODEX(1)


	SELECTCASE(NPR)

	CASE(1)
	RI = POS(1)
	CALL POIN1D (RI,H,P)
      
	CASE(2)
	RI = POS(1)
	SI = POS(2)
	IF(NNM.NE.3) THEN
		IF(NNM.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX,NNM)
		IF(NNM.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX,NNM)				!9 NODE ELEMENT
	ELSEIF(NNM.EQ.3) THEN											!3 NODE ELEMENT
		CALL SHAP2D3(RI,SI,H,P,NNM)
	ENDIF

	CASE(3)
	RI = POS(1)
	SI = POS(2)
	TI = POS(3)
	IF(NNM.EQ.4.OR.NNM.EQ.10) GOTO 250
		CALL SHAP3D8 (RI,SI,TI,H,P) 
		GOTO 251
250	CONTINUE
		CALL SHAP3DT (RI,SI,TI,H,P,NNM)
251	CONTINUE

	ENDSELECT
	
	
	DO ICO = 1,3
	  OFFS(ICO) = 0.0D0
	  DO INM = 1,NNM
	      OFFS(ICO) = OFFS(ICO) + H(INM)*COORD(ICO,INM)
	  ENDDO
	  OFFS(ICO) = COORP(ICO) - OFFS(ICO)
	ENDDO
	
	
	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================

	SUBROUTINE TPITER(MCONT,COORD,XT,YT,ZT,COORN,NODE,IFOUND,NNM,IEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION MCONT(8),COORD(3,8),XYZ(8,3),ANS(3,100),R(3,1),ER(3,100),
     1          TDIFFK(3,3),DIFFK(3,3),RDIFFK(3,3),SH(1,8),P(3,8),
     2          TFREAL(1,3),FREAL(3,1),F(3,1),DELTA(3,1)	
      DIMENSION COORN(3),NODE(1)
C	-----------------------------------------------------------------
      IFOUND = 0
      
	R(1,1) = XT
      R(2,1) = YT
      R(3,1) = ZT
C     Set the initial values in the Natural Coordinates
      xi = 0.0
      eta = 0.0
      zeta = 0.0

      ANS(1,1) = 0.0
      ANS(2,1) = 0.0
      ANS(3,1) = 0.0            

      ER(1,1) = 100.0
      ER(2,1) = 100.0
      ER(3,1) = 100.0
            
      XYZ = TRANSPOSE(COORD) 
            
      DO II = 2,100  !.....................................200 
C     Find the Shape function and Shape function derivative     
      CALL SHAPE8(xi,eta,zeta,SH,P)  
      TDIFFK = MATMUL(P,XYZ)     
      DIFFK = TRANSPOSE(TDIFFK) 
      TFREAL = MATMUL(SH,XYZ)     
      FREAL = TRANSPOSE(TFREAL) 
      F = R - FREAL
      CALL INVMATRIX(DIFFK,RDIFFK,3)  
      DELTA = MATMUL(RDIFFK,F)
C     Update the initial value ---------------------
      ANS(1,II) = ANS(1,II-1) + DELTA(1,1)    
      ANS(2,II) = ANS(2,II-1) + DELTA(2,1)   
      ANS(3,II) = ANS(3,II-1) + DELTA(3,1)  
               
      xi = ANS(1,II)
      eta = ANS(2,II)
      zeta = ANS(3,II)
C       Control Zero Variables  -----------------------               
        IF (ABS(xi) <= 0.00001) THEN
          xi = 0.0
        ENDIF
        IF (ABS(eta) <= 0.00001) THEN
          eta = 0.0
        ENDIF    
        IF (ABS(zeta) <= 0.00001) THEN
          zeta = 0.0
        ENDIF
C       Control Limitation of Natural Coordinates ------               
        IF ((xi < -1.0).OR.(xi > 1.0)) THEN
          EXIT
        ELSEIF ((eta < -1.0).OR.(eta > 1.0)) THEN
          EXIT 
        ELSEIF ((zeta < -1.0).OR.(zeta > 1.0)) THEN
          EXIT                 
        ENDIF
C       Control ERROR 1 -------------------------------            
        IF (xi == 0.0) THEN
          ER(1,II) = 0.0
        ELSE
          ER(1,II) = ABS((ANS(1,II)-ANS(1,II-1))
     1               /ANS(1,II)*100.0)   
        ENDIF
               
        IF (eta == 0.0) THEN
          ER(2,II) = 0.0
        ELSE
          ER(2,II) = ABS((ANS(2,II)-ANS(2,II-1))
     1               /ANS(2,II)*100.0)   
        ENDIF
               
        IF (zeta == 0.0) THEN
          ER(3,II) = 0.0
        ELSE
          ER(3,II) = ABS((ANS(3,II)-ANS(3,II-1))
     1               /ANS(3,II)*100.0)   
        ENDIF
C       Control ERROR 2 ------------------------------- 
        IF (ER(1,II) <= 0.001) THEN
          IF (ER(2,II) <= 0.001) THEN
            IF (ER(3,II) <= 0.001) THEN                         
              ! Keep the CORRECT Results                   
              COORN(1:3) = [xi,eta,zeta]
               NODE(1:8) = MCONT(1:8)   
              ! Step for jumping to the next element
              IFOUND = 1                     
              EXIT 
            ENDIF
          ENDIF
        ENDIF  
                               
      ENDDO !........................................200         

C	----------------------------------------------------------	
	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE SLEN_FRAME(IEG,IEL,SLEN) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
      
    
	CALL XFSECTION(IEG,IEL,1)

	INAME(1:4) = [5,0,1,IEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FTYP ,1,1 ,0) !SECTION TYPE	0=STANDARD 1=FIBER
	CALL MRELFIL(NAMEI,FGAS ,1,3 ,0) !NUMBER OF STATION POINT ALONG ELEMENT
	CALL MRELFIL(NAMEI,FFIB ,1,4 ,0) !NUMBER OF FIBER
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !ROTATION ANGLE
	ITYP = INT(FTYP)
	NGAS = INT(FGAS)
	NFIB = INT(FFIB)

      SLEN = 0.0D0
      ARMX = 0.0D0
C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO GET MAXIMUM CHARACTERISTIC LENGTH SLEN
C     ----------------------------------------------------------
      DO 50 II = 1,NGAS

C	CALL SECTION PROPERTIES
	CALL XFSECTION(IEG,IEL,II)
	
	INAME(1:4) = [5,0,1,IEG] !XSEC
	CALL ICONC(INAME,NAMEI)
C	==========================================================
	CALL MRELFIL(NAMEI,BNORM,1,111,0) !BNORM NORMINAL WIDTH
	CALL MRELFIL(NAMEI,HNORM,1,112,0) !HNORM NORMINAL HIGHT
C	==========================================================
	CALL MRELFIL(NAMEI, AREA,1,21 ,0) !SECTIONAL AREA
      IF(AREA.GT.ARMX) ARMX = AREA
C	==========================================================
	
      CLEN = SQRT(BNORM*BNORM + HNORM*HNORM)
      IF(CLEN.GT.SLEN) SLEN = CLEN
C     ----------------------------------------------------------
50    CONTINUE
C     ----------------------------------------------------------

C     IF SLEN = 0 , LET MAKE IT EQUAL TO SQRT OF AREA
      IF(SLEN.EQ.0.0D0) SLEN = SQRT(ARMX)
      
      
      RETURN
      END
C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE SLEN_SHELL(IGSET,PROPG,NGP,IEG,IEL,SLEN) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION PROPG(NGP,1),IGSET(1)
	
	ISET = IGSET(IEL)
	SLEN = PROPG(2,ISET)
      
      RETURN
      END
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE TFILLDATA(MCONT,COORD,NODEX,MCONT1,COORD1,NODEX1,IEL,NEF,NCO,NNM,NEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      DIMENSION MCONT(NEF,1),COORD(NCO*NNM,1),NODEX(NEX,1)
      DIMENSION MCONT1(NEF),COORD1(NCO*NNM),NODEX1(NEX)

      NODEX1(1:NEX)     = NODEX(1:NEX,IEL)
      MCONT1(1:NEF)     = MCONT(1:NEF,IEL)
      COORD1(1:NCO*NNM) = COORD(1:NCO*NNM,IEL)

      RETURN
      END


C	=====================================================================
C	=====================================================================
C	===================================================================== 
      SUBROUTINE BOXDATA(XYZ,BOUD)   
C    ------------------------------------------------------------------------ 
C     Clear all default name      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /BOX/ NBOX,NLEVL
                
          
      DIMENSION XYZ(NSN,NSC),BOUD(NLEVL,6)     
C    ------------------------------------------------------------------------     
C     Find the maximum and minimum points on the structure
      XXmin = MINVAL(XYZ(1:NSN,1))
      XXmax = MAXVAL(XYZ(1:NSN,1))
      YYmin = MINVAL(XYZ(1:NSN,2))
      YYmax = MAXVAL(XYZ(1:NSN,2))
      ZZmin = MINVAL(XYZ(1:NSN,3)) 
      ZZmax = MAXVAL(XYZ(1:NSN,3))   
      
C    ------------------------------------------------------------------ 
      VAL = REAL(NLEVL)
      DX = (XXmax - XXmin)/VAL
      DY = (YYmax - YYmin)/VAL
      DZ = (ZZmax - ZZmin)/VAL   
      
      XLmin = XXmin
      YLmin = YYmin
      ZLmin = ZZmin

      NBOX = NLEVL*NLEVL*NLEVL           
      
      JJ = 0
C    ------------------------------------------------------------------  
C     Lelvel storage ---------------      
      DO ILX = 1,NLEVL 
      XLmax = XLmin + DX  
      YLmax = YLmin + DY
      ZLmax = ZLmin + DZ
    
      BOUD(ILX,1) = XLmin
      BOUD(ILX,2) = XLmax   
      BOUD(ILX,3) = YLmin
      BOUD(ILX,4) = YLmax   
      BOUD(ILX,5) = ZLmin
      BOUD(ILX,6) = ZLmax   
    
      XLmin = XLmax    
      YLmin = YLmax
      ZLmin = ZLmax      
      ENDDO
      
      RETURN
	END


C	=====================================================================
C	=====================================================================
C	=====================================================================      
	
	SUBROUTINE ELEMENT_BOX_MAP(MCONT,COORD,BOUD,EDATA,ICOUT,MAPGP,COORTT,MCONTT,NGG,IG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT   
      COMMON /BOX/ NBOX,NLEVL 
      COMMON /SELEM/ NSELE
      
      ALLOCATABLE MAPEM(:,:)
      
	DIMENSION MCONT(NEF,1),COORD(NCO*NNM,1),COORTT(NCO*NNM*NELE,NGG),CONTT(NEF*NELE,NGG)
	DIMENSION BOUD(NLEVL,6)
      DIMENSION EDATA(NELE,6),ICOUT(NBOX,1),MAPGP(NBOX,NELE)
      
      ALLOCATE(MAPEM(NBOX,NELE))
      MAPEM = 0.0D0  
      
      DO IEL = 1,NELE
      CALL ELEM_BOX(MCONT(1,IEL),COORD(1,IEL),BOUD,EDATA,ICOUT,MAPEM,MAPGP,IEL,COORTT,MCONTT,NGG,IG)      
      ENDDO
      
      DEALLOCATE(MAPEM)
      
	RETURN
	END


C	=====================================================================
C	=====================================================================
C	=====================================================================  
      SUBROUTINE ELEM_BOX(MCONT,COORD,BOUD,EDATA,ICOUT,MAPEM,MAPGP,IEL,COORTT,MCONTT,NGG,IG)  
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT   
      COMMON /BOX/ NBOX,NLEVL             
           
      DIMENSION MCONT(8),COORD(3,8),COORTT(NCO*NNM*NELE,NGG),MCONTT(NEF*NELE,NGG)
      DIMENSION BOUD(NLEVL,6),MAPEM(NBOX,NELE)
      DIMENSION EDATA(NELE,6),ICOUT(NBOX,1),MAPGP(NBOX,NELE)     

C    ------------------------------------------            
C     BACK UP CONNECTIVITIES & COORDINATES FOR BOX SEARCH IN TENDON LOOP           
      DO INOD = 1,NNM
C     CONNECTIVITIES
      ICON = NEF*(IEL-1) + INOD
      MCONTT(ICON,IG) = MCONT(INOD)
      
C     COORDINATES
      ICOX = NEF*(IEL-1) + NCO*INOD - 2
      ICOY = NEF*(IEL-1) + NCO*INOD - 1
      ICOZ = NEF*(IEL-1) + NCO*INOD - 0
    
      COORTT(ICOX,IG) = COORD(1,INOD)
      COORTT(ICOY,IG) = COORD(2,INOD)
      COORTT(ICOZ,IG) = COORD(3,INOD)
      ENDDO 
C    ------------------------------------------      
   
      
C    ------------------------------------------
C    Mapping Elements to Blocks
C    ------------------------------------------
C     Set the maximum and minimum for each element
C     EDATA = [ELEM Xmin Ymin Zmin Xmax Ymax Zmax]      
C     Element data --------------------------      
      EDATA(IEL,1) = MINVAL(COORD(1,1:8))
      EDATA(IEL,2) = MAXVAL(COORD(1,1:8))
      EDATA(IEL,3) = MINVAL(COORD(2,1:8))
      EDATA(IEL,4) = MAXVAL(COORD(2,1:8))
      EDATA(IEL,5) = MINVAL(COORD(3,1:8))  
      EDATA(IEL,6) = MAXVAL(COORD(3,1:8))   
      ! ----------------------------------------
      XEmin = EDATA(IEL,1)
      XEmax = EDATA(IEL,2)    
      YEmin = EDATA(IEL,3)
      YEmax = EDATA(IEL,4)  
      ZEmin = EDATA(IEL,5)
      ZEmax = EDATA(IEL,6)  
      JJ = 0
      DO ILX = 1,NLEVL        
      XLmin = BOUD(ILX,1)
      XLmax = BOUD(ILX,2)                
      DO ILY = 1,NLEVL           
      YLmin = BOUD(ILY,3)
      YLmax = BOUD(ILY,4)
      DO ILZ = 1,NLEVL                
      ZLmin = BOUD(ILZ,5)
      ZLmax = BOUD(ILZ,6)               
C     ------------------------------------- 
      IBOX = NLEVL*NLEVL*(ILX-1) + NLEVL*(ILY-1) + ILZ
      JJ = JJ + 1
      
C     X - DIRECTION -----------------------           
      IF ((XEmin <= XLmin).AND.(XEmax >= XLmax)) THEN   
      XLmin = XEmin
      XLmax = XEmax 
      CALL CHECKCONDITION(IBOX,IEL,XEmin,XEmax,YEmin,YEmax,
     1           ZEmin,ZEmax,XLmin,XLmax,YLmin,YLmax,ZLmin,ZLmax,MAPEM,NELE)
      ELSE
      CALL CHECKCONDITION(IBOX,IEL,XEmin,XEmax,YEmin,YEmax,
     1           ZEmin,ZEmax,XLmin,XLmax,YLmin,YLmax,ZLmin,ZLmax,MAPEM,NELE)
      ENDIF 
      
C     Y - DIRECTION -----------------------  
      IF ((YEmin <= YLmin).AND.(YEmax >= YLmax)) THEN              
      YLmin = YEmin
      YLmax = YEmax 
      CALL CHECKCONDITION(IBOX,IEL,XEmin,XEmax,YEmin,YEmax,
     1           ZEmin,ZEmax,XLmin,XLmax,YLmin,YLmax,ZLmin,ZLmax,MAPEM,NELE)
      ELSE
      CALL CHECKCONDITION(IBOX,IEL,XEmin,XEmax,YEmin,YEmax,
     1           ZEmin,ZEmax,XLmin,XLmax,YLmin,YLmax,ZLmin,ZLmax,MAPEM,NELE)
      ENDIF 
      
C     Z - DIRECTION -----------------------
      IF ((ZEmin <= ZLmin).AND.(ZEmax >= ZLmax)) THEN            
      ZLmin = ZEmin
      ZLmax = ZEmax 
      CALL CHECKCONDITION(IBOX,IEL,XEmin,XEmax,YEmin,YEmax,
     1           ZEmin,ZEmax,XLmin,XLmax,YLmin,YLmax,ZLmin,ZLmax,MAPEM,NELE)
      ELSE
      CALL CHECKCONDITION(IBOX,IEL,XEmin,XEmax,YEmin,YEmax,
     1           ZEmin,ZEmax,XLmin,XLmax,YLmin,YLmax,ZLmin,ZLmax,MAPEM,NELE)
      ENDIF             
                                                                
      A = MAPEM(IBOX,IEL)
      IF (A /= 0) THEN
      ICOUT(IBOX,1) = ICOUT(IBOX,1) + 1
      JJ = ICOUT(IBOX,1)
      MAPGP(IBOX,JJ) = MAPEM(IBOX,IEL) 
      ENDIF                
                 
C     -------------------------------------           
      ENDDO
      ENDDO        
      ENDDO  
            
      
      RETURN
	END
	
	
C	=====================================================================
C	=====================================================================
C	===================================================================== 
      SUBROUTINE CHECKCONDITION(IB,IE,XEmin,XEmax,YEmin,YEmax,
     1           ZEmin,ZEmax,XLmin,XLmax,YLmin,YLmax,ZLmin,ZLmax,MAPEM,NELE)
     
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
            
     
      COMMON /BOX/ NBOX,NLEVL

      DIMENSION MAPEM(NBOX,NELE) 
     
C     ------------------------------------------------------------      
      IF ((XEmin >= XLmin).AND.(XEmin <= XLmax)) THEN
        IF ((YEmin >= YLmin).AND.(YEmin <= YLmax)) THEN 
            IF ((ZEmin >= ZLmin).AND.(ZEmin <= ZLmax)) THEN   
                MAPEM(IB,IE) = IE       
            ENDIF    
            IF ((ZEmax >= ZLmin).AND.(ZEmax <= ZLmax)) THEN    
                MAPEM(IB,IE) = IE         
            ENDIF
        ENDIF    
        IF ((YEmax >= YLmin).AND.(YEmax <= YLmax)) THEN 
            IF ((ZEmin >= ZLmin).AND.(ZEmin <= ZLmax)) THEN   
                MAPEM(IB,IE) = IE           
            ENDIF    
            IF ((ZEmax >= ZLmin).AND.(ZEmax <= ZLmax)) THEN    
                MAPEM(IB,IE) = IE         
            ENDIF
        ENDIF 
      ENDIF   
C     ------------------------------------------------------------        
      IF ((XEmax >= XLmin).AND.(XEmax <= XLmax)) THEN 
        IF ((YEmin >= YLmin).AND.(YEmin <= YLmax)) THEN 
            IF ((ZEmin >= ZLmin).AND.(ZEmin <= ZLmax)) THEN   
                MAPEM(IB,IE) = IE         
            ENDIF    
            IF ((ZEmax >= ZLmin).AND.(ZEmax <= ZLmax)) THEN    
                MAPEM(IB,IE) = IE       
            ENDIF
        ENDIF    
        IF ((YEmax >= YLmin).AND.(YEmax <= YLmax)) THEN 
            IF ((ZEmin >= ZLmin).AND.(ZEmin <= ZLmax)) THEN   
                MAPEM(IB,IE) = IE           
            ENDIF    
            IF ((ZEmax >= ZLmin).AND.(ZEmax <= ZLmax)) THEN    
                MAPEM(IB,IE) = IE         
            ENDIF
        ENDIF     
      ENDIF 
C     ------------------------------------------------------------ 
      RETURN
	END


C	=====================================================================
C	=====================================================================
C	=====================================================================  
      SUBROUTINE TENDON_BOX(BOUD,MAPTEN,NTPOIN,XT,YT,ZT,IP)         
C    ---------------------------------------------------------------------- 
C     Clear all default name      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	COMMON BLOCK FOR TEND PRAMIN  NOV2010 	
	COMMON /BOX/ NBOX,NLEVL
           
      DIMENSION BOUD(NLEVL,6),MAPTEN(NTPOIN)      
C    --------------------------------------
C    Mapping Tendon nodes to Blocks
C    --------------------------------------          
      Xnod = XT    
      Ynod = YT 
      Znod = ZT
         
      DO ILX = 1,NLEVL        
      XLmin = BOUD(ILX,1)
      XLmax = BOUD(ILX,2)
      IF ((Xnod >= XLmin).AND.(Xnod <= XLmax)) THEN            
      DO ILY = 1,NLEVL           
      YLmin = BOUD(ILY,3)
      YLmax = BOUD(ILY,4)            
      IF ((Ynod >= YLmin).AND.(Ynod <= YLmax)) THEN     
      DO ILZ = 1,NLEVL                
      ZLmin = BOUD(ILZ,5)
      ZLmax = BOUD(ILZ,6)               
      !-------------------------------------                         
      IF ((Znod >= ZLmin).AND.(Znod <= ZLmax)) THEN  
      IBOX = NLEVL*NLEVL*(ILX-1) + NLEVL*(ILY-1) + ILZ 
      MAPTEN(IP) = IBOX 
      EXIT           
      ENDIF          
      !-------------------------------------           
      ENDDO
      EXIT
      ENDIF
      ENDDO
      EXIT 
      ENDIF
      ENDDO           
C    ----------------------------------------------------------------------
      RETURN
	END
      
      
C	=====================================================================
C	=====================================================================
C	=====================================================================      
	
	SUBROUTINE TENDON_SEARCH(MCONT,COORD,NODEX,XT,YT,ZT,COORN,NODE,NCO,NEF,NNM,
     1                         TENELM,ITPOIN,NTPOIN,NN,NELE,IEL,COORTT,NSELE,NEX,
     1                         MCONTT,NGG,IG,NPR,OFFS,H)
     
      
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)  
      
C	COMMON BLOCK FOR TEND PRAMIN  NOV2010 	
	COMMON /BOX/ NBOX,NLEVL         
      
	DIMENSION MCONT(NEF,1),COORD(NCO*NNM,1),COORN(3),NODEX(NEX,1)
	DIMENSION COORTT(NCO*NNM*NELE,NGG),MCONTT(NEF*NELE,NGG)
	DIMENSION TENELM(NTPOIN,7) 
	DIMENSION POSN(1),OFFS(1),NODE(1),COORP(3),H(NNM)
	 
	COORP(1:3) = [XT,YT,ZT] 
C     CHECKING THE POINTS CORRESPONDING WITH ALL CONDITIONS                      
      CALL LOCATEPOINT(MCONT(1,IEL),COORD(1,IEL),NODEX(1,IEL),XT,YT,ZT,COORN,NODE,TENELM,
     1                 ITPOIN,IEL,NTPOIN,COORTT,MCONTT,NSELE,NGG,IG,H) 
     
C     CALCULATE OFFSET     
      OFFS(1) = 0.0D0
      OFFS(2) = 0.0D0
      OFFS(3) = 0.0D0   
      
      SELECTCASE(NPR)
      CASE(1) !FRAME
	  IF(SLEN.GT.SMAX) ITEST = 1 !NOT YET FOUND
      CASE(2) !SHELL
	  IF(SLEN.GT.SMAX) ITEST = 1 !NOT YET FOUND
      CASE(3) !SOLID
C	  NOTHING FOR SOLID
      ENDSELECT	
      
           
	RETURN
	END	
	
C	=======================================================================	
C	=======================================================================       
C	======================================================================= 
      SUBROUTINE LOCATEPOINT(MCONT,COORD,NODEX,XT,YT,ZT,COORN,NODE,TENELM,
     1                       IP,IEL,NTPOIN,COORTT,MCONTT,NSELE,NGG,IG,HR)                 
C     -----------------------------------------------------------------------     
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT      
      COMMON /BOX/ NBOX,NLEVL
     
      DIMENSION TENELM(NTPOIN,7),COORTT(NCO*NNM*NELE,NGG),MCONTT(NEF*NELE,NGG)  
                 
      DIMENSION MCONT(8),COORD(3,8),XYZ(8,3),ANS(3,100),R(3,1),ER(3,100),
     1          TDIFFK(3,3),DIFFK(3,3),RDIFFK(3,3),SH(1,8),P(3,8),
     2          TFREAL(1,3),FREAL(3,1),F(3,1),DELTA(3,1)	
      DIMENSION COORN(3),NODE(1),HR(NNM),NODEX(1)
             
C     ----------------------------------------     
      DO II = 1,NEX
      NODEX(II) = 0
      ENDDO  
C     2.2) Create the nodal displacements in each element
            
      R(1,1) = XT
      R(2,1) = YT
      R(3,1) = ZT
C     Guess the initial values in the Natural Coordinate 
      xi   = 0.0D0
      eta  = 0.0D0
      zeta = 0.0D0

      ANS(1,1) = xi
      ANS(2,1) = eta
      ANS(3,1) = zeta

C     ---------------------------------------- 
C     ITERATIONS                
C     ---------------------------------------- 
      ER(1,1) = 100.0D0
      ER(2,1) = 100.0D0
      ER(3,1) = 100.0D0      
          
C     ------------------------------------------            
C     BACK UP CONNECTIVITIES & COORDINATES FOR BOX SEARCH IN TENDON LOOP           
C      DO INOD = 1,NNM
C     CONNECTIVITIES
C      ICON = NEF*(IEL-1) + INOD
C      MCONT(INOD) = MCONTT(ICON,IG) 
      
C     COORDINATES
C      ICOX = NEF*(IEL-1) + NCO*INOD - 2
C      ICOY = NEF*(IEL-1) + NCO*INOD - 1
C      ICOZ = NEF*(IEL-1) + NCO*INOD - 0
    
C      COORD(1,INOD) = COORTT(ICOX,IG) 
C      COORD(2,INOD) = COORTT(ICOY,IG)
C      COORD(3,INOD) = COORTT(ICOZ,IG) 
C      ENDDO 
C     ------------------------------------------       
      
      XYZ = TRANSPOSE(COORD)       
C     ---------------------------------------------
      DEC = 1.0010D0     
      
      DO 200 II = 2,100  
      
C     Shape Function
      CALL SHAPE8(xi,eta,zeta,SH,P) 
C     ---------------------------------------------
      DIFFK = MATMUL(TRANSPOSE(XYZ),TRANSPOSE(P))     
     
      DO I = 1,3
      FREAL(I,1) = 0.0D0
      DO J = 1,8
      FREAL(I,1) = FREAL(I,1) + XYZ(J,I)*SH(1,J)
      ENDDO
      ENDDO
      
      F = R - FREAL       
      CALL INVMATRIX(DIFFK,RDIFFK,3)          
      DELTA = MATMUL(RDIFFK,F) 
        
C     Control Values ---------------------------
      CC = ABS(DELTA(1,1))
      IF (CC <= 0.000010D0) THEN
      DELTA(1,1) = 0.0D0
      ENDIF 
      CC = ABS(DELTA(2,1))   
      IF (CC <= 0.00001) THEN
      DELTA(2,1) = 0.0D0
      ENDIF 
      CC = ABS(DELTA(3,1))
      IF (CC <= 0.000010D0) THEN
      DELTA(3,1) = 0.0D0
      ENDIF 
C     ---------------------------------------------

C     Update the initial value      
      ANS(1,II) = ANS(1,II-1) + DELTA(1,1)
      ANS(2,II) = ANS(2,II-1) + DELTA(2,1)
      ANS(3,II) = ANS(3,II-1) + DELTA(3,1)
      
      xi   = ANS(1,II)
      eta  = ANS(2,II)
      zeta = ANS(3,II)            

C     Control Values ---------------------------
      CC = ABS(xi)
      IF (CC <= 0.00010D0) THEN
      xi = 0.0D0
      ENDIF  
      CC = ABS(eta)  
      IF (CC <= 0.00010D0) THEN
      eta = 0.0D0
      ENDIF 
      CC = ABS(zeta)
      IF (CC <= 0.00010D0) THEN
      zeta = 0.0D0
      ENDIF 

C     Control ERROR ----------------------------    
      IF (xi == 0.0) THEN
      ER(1,II) = 0.0D0
      ELSE
      ER(1,II) = ABS((ANS(1,II)-ANS(1,II-1))/ANS(1,II)*100.0D0)
      ENDIF
      
      IF (eta == 0.0) THEN
      ER(2,II) = 0.0D0
      ELSE
      ER(2,II) = ABS((ANS(2,II)-ANS(2,II-1))/ANS(2,II)*100.0D0)
      ENDIF
      
      IF (zeta == 0.0) THEN
      ER(3,II) = 0.0D0
      ELSE
      ER(3,II) = ABS((ANS(3,II)-ANS(3,II-1))/ANS(3,II)*100.0D0)
      ENDIF         

      EX = ER(1,II)
      EY = ER(2,II)
      EZ = ER(3,II)
      
C     Keep all values matching the conditions       
      IF (EX <=0.00010D0) THEN
      IF (EY <=0.00010D0) THEN
      IF (EZ <=0.00010D0) THEN                                             
C     ---------------------------------------------
      IF ((xi >= -DEC).AND.(xi <= DEC)) THEN 
      IF ((eta >= -DEC).AND.(eta <= DEC)) THEN 
      IF ((zeta >= -DEC).AND.(zeta <= DEC)) THEN
      TENELM(IP,1) = REAL(IEL) 
               
      CALL SHAPE8(xi,eta,zeta,SH,P)
      TENELM(IP,2) = SH(1,1)*XYZ(1,1) + SH(1,2)*XYZ(2,1) 
     1             + SH(1,3)*XYZ(3,1) + SH(1,4)*XYZ(4,1)
     1             + SH(1,5)*XYZ(5,1) + SH(1,6)*XYZ(6,1)
     1             + SH(1,7)*XYZ(7,1) + SH(1,8)*XYZ(8,1)
      TENELM(IP,3) = SH(1,1)*XYZ(1,2) + SH(1,2)*XYZ(2,2)
     1             + SH(1,3)*XYZ(3,2) + SH(1,4)*XYZ(4,2)
     1             + SH(1,5)*XYZ(5,2) + SH(1,6)*XYZ(6,2)
     1             + SH(1,7)*XYZ(7,2) + SH(1,8)*XYZ(8,2)
      TENELM(IP,4) = SH(1,1)*XYZ(1,3) + SH(1,2)*XYZ(2,3)
     1             + SH(1,3)*XYZ(3,3) + SH(1,4)*XYZ(4,3)
     1             + SH(1,5)*XYZ(5,3) + SH(1,6)*XYZ(6,3)
     1             + SH(1,7)*XYZ(7,3) + SH(1,8)*XYZ(8,3)
                        
      TENELM(IP,5) = xi 
      TENELM(IP,6) = eta 
      TENELM(IP,7) = zeta 
      
      COORN(1:3) = [xi,eta,zeta]
      NODE(1:8) = MCONT(1:8)   
      !Step for jumping to the next element  
      
C     RETURN SHAPE FUNCTION
      DO JJ = 1,NNM
      HR(JJ) = SH(1,JJ)
      ENDDO      
                
      ENDIF                            
      ENDIF                             
      ENDIF       
C     ---------------------------------------------
      EXIT
      
      ENDIF
      ENDIF    
      ENDIF    

200   CONTINUE  

 
      
      RETURN
	END
C	=======================================================================	
C	=======================================================================       
C	======================================================================= 


   
