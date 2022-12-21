      SUBROUTINE ILINKCON
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     READ LINK CONSTRAINTS DATA BY BJ      
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NEDATA        = SAVE NODE DATA IN EACH FLOORF
C     ICLO          = THE NUMBER OF CENTER LOAD DATA
C     ILDEC         = DIRECTION NUMBER OF CENTER LOAD
C     CLOAD         = VALUE OF CENTER LOAD
C     -----------------------------------------------------
C   
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
     
!      COMMON /DIAPH/ NUMF,IFNUM(500),NUMN(1000),NUMVT,NEDATA(1000,1000)
      
      COMMON /DAPH/ DIAPH(700000),IDIAPH(600000) ! STORAGE FOR CONSTRAINTS BY BJ
      
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
          
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      
      COMMON /NUMFLO/ NUMF,NUMVT
      COMMON /IEQNUMFLO/ IEQNUM  
        
C      DIMENSION IFNUM(1),NUMN(1),FNEDATA(1)     
     
	READ (ITI,*)
      READ (ITI,*) NUMF,IEQNUM
     
     
	IF(NUMF.EQ.0) RETURN
	
C	    DO INUMF = 1,NUMF
	    JIF = 1
          JIFD = NUMF
          CALL UPDIA(JNUM,JIF,JIFD)
          
          JNUMD = NUMF
          CALL UPDIA(JNED,JNUM,JNUMD)

          CALL READNOD(IDIAPH(JIF),IDIAPH(JNUM))	
C          CALL READNOD(IFNUM,NUMN,NEDATA)	
      
C      ENDDO
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================       
      SUBROUTINE READNOD(IFNUM,NUMN)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     READ DIAPHRAGM CONSTRAINTS DATA BY BJ
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NEDATA        = SAVE NODE DATA IN EACH FLOORF
C     ICLO          = THE NUMBER OF CENTER LOAD DATA
C     ILDEC         = DIRECTION NUMBER OF CENTER LOAD
C     CLOAD         = VALUE OF CENTER LOAD
C     -----------------------------------------------------
C   
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
     
C      COMMON /DIAPH/ NUMF,IFNUM(500),NUMN(1000),NUMVT,NEDATA(1000,1000)

      COMMON /DAPH/ DIAPH(700000),IDIAPH(600000) ! STORAGE FOR CONSTRAINTS BY BJ
      
      COMMON /DIAPH_DATA/ NSLAVE_NODE(500,100),NSLAVE(500)
   
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ       

      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      
      COMMON /NUMFLO/ NUMF,NUMVT
      
      DIMENSION IFNUM(NUMF),NUMN(NUMF)
C      DIMENSION IFNUM(NUMF),NUMN(NUMF),NEDATA(1)
      
      
      TOTOAL_LOOP = NUMF
      DO 100 INUMF = 1,TOTOAL_LOOP
     
        READ (ITI,*) IFNUM(INUMF),OPT_CON,NUMN(INUMF),NUMVT
        
        IF (OPT_CON.EQ.0) THEN
        NUMF =  NUMF - 1   
        JNODE = NUMN(INUMF)
        READ (ITI,*) NSLAVE_NODE(INUMF,1:JNODE)
        NSLAVE(INUMF)             = NUMN(INUMF)
        GOTO 100
        ENDIF
        
        JNODE = NUMN(INUMF)
        NSLAVE(INUMF)             = NUMN(INUMF)
        
        CALL READNE(INUMF,JNODE,IDIAPH(JNED))
C        CALL READNE(INUMF,JNODE,NEDATA)        
        
100   CONTINUE
      
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================       
      SUBROUTINE READNE(INUMF,JNODE,NEDATA)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     READ DIAPHRAGM CONSTRAINTS DATA BY BJ
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NEDATA        = SAVE NODE DATA IN EACH FLOORF
C     ICLO          = THE NUMBER OF CENTER LOAD DATA
C     ILDEC         = DIRECTION NUMBER OF CENTER LOAD
C     CLOAD         = VALUE OF CENTER LOAD
C     -----------------------------------------------------
C   
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
     
C      COMMON /DIAPH/ NUMF,IFNUM(500),NUMN(1000),NUMVT,NEDATA(1000,1000)

      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /DIAPH_DATA/ NSLAVE_NODE(500,100),NSLAVE(500)
      
      COMMON /NUMFLO/ NUMF,NUMVT
  
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      DIMENSION NEDATA(NUMF,JNODE)
      
              
        READ (ITI,*) NEDATA(INUMF,1:JNODE)
        !NSLAVE(INUMF)             = INUMF
        NSLAVE_NODE(INUMF,1:JNODE) = NEDATA(INUMF,1:JNODE)
        
        IF(INUMF.EQ.1) JNEDD = 0
        
        IF(INUMF.EQ.1) JMAX = 0
          
          JNEDD = JNEDD + JNODE
          
          IF(JNODE.GT.JMAX) JMAX = JNODE
        
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================    
      SUBROUTINE UPDIA(JA,JB,JC)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     UPDATE INDICATOR OF DIA MATRIX FOR DIAPHRAGM CONSTRAINTS BY BJ
C	-----------------------     

      JA = JB + JC

      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE CALCOORD(COORDNEW)
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
      
      COMMON /LMN/ LMNEW(4,1000) !TO SAVE THE LM MATRIX INTO LMS FOR DIAPHRAGM CONSTRAITNS BY BJ  
      
      COMMON A(9000000),IA(9000000)
      
      COMMON /DAPH/ DIAPH(700000),IDIAPH(600000) ! STORAGE FOR CONSTRAINTS BY BJ
       
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ   
      
C      COMMON /DIAPH/ NUMF,IFNUM(500),NUMN(1000),NUMVT,NEDATA(1000,1000)
      
      
      DIMENSION COORDNEW(NSN,3)
     
        
        JCOD = NSN*3

        COORDNEW(1:NSN,1) = A(LXY:LXY+(NSN-1))
        COORDNEW(1:NSN,2) = A(LXY+NSN:LXY+(2*NSN-1))
        COORDNEW(1:NSN,3) = A(LXY+(2*NSN):LXY+(3*NSN-1))
        
        CALL UPDIA(JCXY,JCO,JCOD)
 
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE CALCENTRO(COORDNEW,CXY,IFNUM,NUMN,NEDATA)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     CALCULATE CENTROID POINT OF MASTER NODE FOR DIAPHRAGM CONSTRAINTS BY BJ
C     (GEOMETRIC CENTER)      
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NODEDATA      = SAVE NODE DATA IN EACH FLOORF
C     -----------------------------------------------------
C    
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
       
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
	
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
	
 
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      COMMON A(9000000),IA(9000000)
             
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ             
C      COMMON /DIAPH/ NUMF,IFNUM(500),NUMN(1000),NUMVT,NEDATA(1000,1000)
      COMMON /NUMFLO/ NUMF,NUMVT
      
      COMMON /LMN/ LMNEW(4,1000) !TO SAVE THE LM MATRIX INTO LMS FOR DIAPHRAGM CONSTRAITNS BY BJ  
              
      COMMON /NEL2/ NELE2 ! FOR DIAPHRAGM CONSTRAINTS BY BJ
      
      COMMON /AREASA/ ARAA(1000,2)
      
      
      COMMON /DGEOMET/ IGSETD(5000,2) !STORE ASSIGNED GEO.PROP.
      
	COMMON /MEMBA/ JMEMBT(5000,3),JMEMBS(5000,3),JMEMBF(5000,3), !SAVE ELEMENT NUMBER FOR TRUSS & SHELL FOR DIAPHRAGM BY BJ
     1	              AMEMBT(5000,2),AMEMBS(5000,2),AMEMBF(5000,2),
     2                NELET,NELES,NELEF 
      COMMON /KTYP/ KTYPEM(3)  ! FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      
      DIMENSION IFNUM(NUMF),NUMN(NUMF),NEDATA(NUMF,JMAX)
      
      DIMENSION COORDNEW(NSN,3),CXY(NUMF,3)!,IFELEM(NUMF,10000),IFELEM2(NUMF,10000),MNEO(NUMF),NELEDATA(NUMF,10000)

      DIMENSION JIFEM(NUMF)
      
C      NELE = NELE2
      NELE = NELET + NELES + NELEF
 
      IT = NUMVT ! VERTICAL DIRECTION NUMBER OF FLOOR
        
        IF(IT.EQ.1) THEN
            IR = 2
            IS = 3
        ELSEIF(IT.EQ.2) THEN   
            IR = 1
            IS = 3             
        ELSEIF(IT.EQ.3) THEN   
            IR = 1
            IS = 2             
        ENDIF 
      
      
      DO INUMF = 1,NUMF 
          
          IFL = 0
          IEL3 = 0
          
          MXRX = 0.0
          MNRX = 0.0

          MXSY = 0.0
          MNSY = 0.0
          
          JNUMN = NUMN(INUMF)
        
        DO IJ = 1,JNUMN
           JNODE = NEDATA(INUMF,IJ)
           
           RX = COORDNEW(JNODE,IR)
           SY = COORDNEW(JNODE,IS)
           
           Z1 = COORDNEW(JNODE,IT) ! COORDINATE IN VERTICAL DIRECTION
           
C           RX = ABS(RX)
C           SY = ABS(SY)
           
           IF(RX.GT.MXRX) MXRX = RX !MAXIMUM COORDINATE OF X
           IF(RX.LT.MNRX) MNRX = RX !MINIMUM COORDINATE OF X

           IF(SY.GT.MXSY) MXSY = SY !MAXIMUM COORDINATE OF X
           IF(SY.LT.MNSY) MNSY = SY !MINIMUM COORDINATE OF X
         
CS           CX = (MXRX-MNRX)/2.0
CS           CY = (MXSY-MNSY)/2.0
           CX = (MXRX+MNRX)/2.0
           CY = (MXSY+MNSY)/2.0
           CZ = Z1
           
        ENDDO
        
           CXY(INUMF,IR) = CX !SAVE THE CENTROID COORDINATE IN TO CXY MATRIX                                         
           CXY(INUMF,IS) = CY
           CXY(INUMF,IT) = CZ

      ENDDO
       
      JCXYD = NUMF*3
      
C      CALL UPDIA(JDTR,JCXY,JCXYD)
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================   
      SUBROUTINE CALCENTROMASS(COORDNEW,CXY,IFNUM,NUMN,NEDATA)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------------------------
C     CALCULATE CENTROID POINT OF MASTER NODE FOR DIAPHRAGM CONSTRAINTS BY BJ
C     (CENTER OF MASS)      
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NODEDATA      = SAVE NODE DATA IN EACH FLOORF
C     -----------------------------------------------------
C    
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
       
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC
	
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
	
 
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      
      COMMON A(9000000),IA(9000000)
             
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ             
C      COMMON /DIAPH/ NUMF,IFNUM(500),NUMN(1000),NUMVT,NEDATA(1000,1000)
      COMMON /NUMFLO/ NUMF,NUMVT
      
      COMMON /LMN/ LMNEW(4,1000) !TO SAVE THE LM MATRIX INTO LMS FOR DIAPHRAGM CONSTRAITNS BY BJ  
              
      COMMON /NEL2/ NELE2 ! FOR DIAPHRAGM CONSTRAINTS BY BJ
      
      COMMON /AREASA/ ARAA(1000,2)
      
      COMMON /DGEOMET/ IGSETD(5000,2) !STORE ASSIGNED GEO.PROP.
      
	COMMON /MEMBA/ JMEMBT(5000,3),JMEMBS(5000,3),JMEMBF(5000,3), !SAVE ELEMENT NUMBER FOR TRUSS & SHELL FOR DIAPHRAGM BY BJ
     1	              AMEMBT(5000,2),AMEMBS(5000,2),AMEMBF(5000,2),
     2                NELET,NELES,NELEF 
      COMMON /KTYP/ KTYPEM(3)  ! FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      
      DIMENSION IFNUM(NUMF),NUMN(NUMF),NEDATA(NUMF,JMAX)
      
      DIMENSION COORDNEW(NSN,3),CXY(NUMF,3)!,IFELEM(NUMF,10000),IFELEM2(NUMF,10000),MNEO(NUMF),NELEDATA(NUMF,10000)

      DIMENSION JIFEM(NUMF)
      
C      NELE = NELE2
      NELE = NELET + NELES + NELEF
 
      IT = NUMVT ! VERTICAL DIRECTION NUMBER OF FLOOR
        
        IF(IT.EQ.1) THEN
            IR = 2
            IS = 3
        ELSEIF(IT.EQ.2) THEN   
            IR = 1
            IS = 3             
        ELSEIF(IT.EQ.3) THEN   
            IR = 1
            IS = 2             
        ENDIF 
      
      
      DO INUMF = 1,NUMF 
          
          IFL = 0
          IEL3 = 0
          
          MXRX = 0.0
          MNRX = 0.0

          MXSY = 0.0
          MNSY = 0.0
          
          JNUMN = NUMN(INUMF)
        
        DO IJ = 1,JNUMN
           JNODE = NEDATA(INUMF,IJ)
           
           RX = COORDNEW(JNODE,IR)
           SY = COORDNEW(JNODE,IS)
           
           Z1 = COORDNEW(JNODE,IT) ! COORDINATE IN VERTICAL DIRECTION
           
C           RX = ABS(RX)
C           SY = ABS(SY)
           
           IF(RX.GT.MXRX) MXRX = RX !MAXIMUM COORDINATE OF X
           IF(RX.LT.MNRX) MNRX = RX !MINIMUM COORDINATE OF X

           IF(SY.GT.MXSY) MXSY = SY !MAXIMUM COORDINATE OF X
           IF(SY.LT.MNSY) MNSY = SY !MINIMUM COORDINATE OF X
         
CS           CX = (MXRX-MNRX)/2.0
CS           CY = (MXSY-MNSY)/2.0
           CX = (MXRX+MNRX)/2.0
           CY = (MXSY+MNSY)/2.0
           CZ = Z1
           
        ENDDO
        
           CXY(INUMF,IR) = CX !SAVE THE CENTROID COORDINATE IN TO CXY MATRIX                                         
           CXY(INUMF,IS) = CY
           CXY(INUMF,IT) = CZ

      ENDDO
       
      JCXYD = NUMF*3
      
C      CALL UPDIA(JDTR,JCXY,JCXYD)
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================         
      SUBROUTINE CONECT2 (NCON,LM,MSN,MEF,MAXN,NELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     --------------------------------------------------------------
C     SETS UP NODAL CONNECTIVITY TABLE, WHICH IS LATER NEEDED FOR
C     BANDWIDTH OPTIMISATION. NOTE ZERO ENTRIES IN INCIDENCES ARRAY
C     (LM) INDICATE THAT THIS NODE WAS NOT USED TO DEFINE AN ELEMENT
C	--------------------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     NCON(NSN,50)  = NODAL CONNECTIVITY TABLE (MAX.CONNECTIONS = 50)
C                     NCON(I,J) = NODE NUMBER OF JTH CONNECTION TO
C                     NODE NUMBER I. A CONNECTION OF A NODE TO ITSELF
C                     IS NOT LISTED.
C     LM(NEF,NELE)  = ELEMENT CONNECTIVITY ARRAY
C     MSN = NSN     = TOTAL NUMBER OF NODES IN STRUCTURE
C     MEF = NEF     = NUMBER OF ELEMENT DEGREES OF FREEDOM
C     MAXN          = MAX. PERMITTED NUM OF NODES FOR ANY ONE ELEMENT
C     NELE          = NUMBER OF ELEMENTS IN THIS GROUP
C     MAXCON        = MAXIMUM PERMITTED NUMBER OF COUPLING NODES
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /BOPT/ IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
      DIMENSION NCON(NSNL,1),LM(MEF,1)

C
      IF (IREDU.GT.0) GOTO 200
C     -------------------------------------------
C     NO BANDWIDTH OPTIMISATION REQUESTED
C     FLAG ACTIVE NODES BY SETTING NCON(IACT) = 1
C     -------------------------------------------
      DO 100  IEL=1,NELE
      DO 100  INO=1,MAXN
      NODE = LM(INO,IEL)
      IF (NODE.EQ.0) GOTO 100
      NCON(NODE,1) = 1
 100  CONTINUE
      NCON(MSN+1:NSNL,1) = 0
      RETURN
C     ---------------------------------
C     LOOP OVER ALL ELEMENTS IN GROUP
C     LOOP OVER NUMBER OF ELEMENT NODES
C     ---------------------------------
 200  DO 900  IEL=1,NELE
      DO 900  INO=1,MAXN
      NOMAS = LM(INO,IEL)
      IF (NOMAS.EQ.0) GOTO 900
C     ----------------------------------------------------------
C     LOOP OVER NODES WHICH ARE COUPLED WITH MASTER NODE (NOMAS)
C     (CONNECTIONS OF A NODE TO ITSELF ARE NOT LISTED)
C     ----------------------------------------------------------
      DO 800  JNO=1,MAXN
      NOCOP = LM(JNO,IEL)
      IF (NOCOP.EQ.NOMAS) GOTO 800
      IF (NOCOP.EQ.0) GOTO 800
C     -----------------------------------------
C     FIND POSITION OF COUPLING NODE (NOCOP) IN
C     CONNECTIVITY TABLE (ROW NOMAS)
C     -----------------------------------------
      DO 500  JPOS=1,MAXCON
      IF (NCON(NOMAS,JPOS).EQ.0) NCON(NOMAS,JPOS) = NOCOP
      IF (NOCOP-NCON(NOMAS,JPOS)) 410,800,500
C     ------------------------------------------------------------
C     INSERT NEW COUPLING NODE AND SHIFT HIGHER COUPLING POSITIONS
C     ------------------------------------------------------------
 410  IF (JPOS.EQ.MAXCON) GOTO 510
      DO 490  KPOS=JPOS,MAXCON
      NOTEMP = NCON(NOMAS,KPOS)
      NCON(NOMAS,KPOS) = NOCOP
      NOCOP = NOTEMP
      IF (KPOS.GT.MCON) MCON = KPOS
 490  IF (NOTEMP.EQ.0) GOTO 800
C
 500  CONTINUE
 510  CALL ERRORS (20,IEL,MAXCON,'CONECT    ')
C     -------------------------------------------
C     END OF LOOP OVER COUPLING NODES
C     END OF LOOP OVER ELEMENT NODES AND ELEMENTS
C     -------------------------------------------
 800  CONTINUE
 900  CONTINUE
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================      
	SUBROUTINE INPLINK2(NCON,IND)
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
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
      
      COMMON /NUMFLO/ NUMF,NUMVT
      
      COMMON /DAPH/ DIAPH(700000),IDIAPH(600000) ! STORAGE FOR CONSTRAINTS BY BJ
      
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /CNODE/ NCNODE(5000) !SAVE CENTER NODE NUMBER
      
      COMMON /TOMAT/ TEMAT,EMAXM,TEMAT3
      COMMON /ITOMAT/ IEMAT
C	-------------------------------------------------------------------            
      
	DIMENSION NCON(NSNL,1)
	DIMENSION LEFLNK(NSF),IDSP(NSF,NSN)
	
	DIMENSION NUMN(NUMF), INUMAT(NUMF), LNODE(50000), ICNODE(50000)

	IF(IND.EQ.1) GOTO 100

	LNKLEN = 30
	
	READ(ITI,*) 
	READ(ITI,*) NGLSET,NGLELM
	
		! SET NGLSET & NGLELM -------------------
	NGLSET = NFLOOR
	
	NGLELM = 0
	
	NUMN(1:NUMF) = IDIAPH(JNUM:JNUM+NUMF-1)
      
   	DO INUMF = 1,NUMF
        INUM = NUMN(INUMF)
        INUM2 = INUM*2
        NGLELM0 = NGLELM0 + INUM
	  NGLELM = NGLELM + INUM2
        
      ENDDO          
	
	
      KNUM = 0
	DO INUMF = 1,NUMF
	  INUM = NUMN(INUMF)

	  DO IN = 1,INUM
	      INK = KNUM + 3*(IN)
	      LEKNOD(INK) = INUMF	! SAVE ISET INTO LEKNOD MATRIX (THIRD COLUMN) 
C           LEKNOD(INK+(NGLELM*3)) = INUMF !LEKNOD(SECOND COLUMN)  
            LEKNOD(INK+3*NGLELM0) = INUMF
	  ENDDO
	  KNUM = INK 
	ENDDO
      !----------------------------------------------------------------------

	IF(NGLSET.EQ.0) RETURN


C	SET NOLIN TO BE 1  SO THAT IT WILL ACTIVATE AN ITERATION OPTION IN LINEAR ANALYSIS
	IF(ISOLOP.EQ.1) NOLIN = 1


      ! MAKE DATLNK MATRIX -------------
      DATLNK(1:5000) = 0.0
      DO INUMF = 1,NUMF
        JNUM = LNKLEN*(INUMF-1)
        DATLNK(1+JNUM) = INUMF
        DATLNK(2+JNUM) = 1
C       DATLNK(3+JNUM:6+JNUM) = 0
     
      ! LINK STIFFNESS VALUES
        
c        DATLNK(7+JNUM)  =  TEMAT3  !r
c        DATLNK(8+JNUM)  =  2.1e-20 !s = Zero value
c        DATLNK(9+JNUM)  =  TEMAT3  !t
c        DATLNK(10+JNUM) =  2.1e-20!theta r = Zero value
c        DATLNK(11+JNUM) =  TEMAT3  !theta s 
c        DATLNK(12+JNUM) =  TEMAT3   !theta t
        
C        DATLNK(7+JNUM)  =  TEMAT3 !r
C        DATLNK(8+JNUM)  =  TEMAT3 !s = Zero value
C        DATLNK(9+JNUM)  =  2.1E-20 !t
C        DATLNK(10+JNUM) =  2.1E-20 !theta r = Zero value
C        DATLNK(11+JNUM) =  2.1E-20 !theta s 
C        DATLNK(12+JNUM) =  TEMAT3  !theta t
        
C        DATLNK(7+JNUM)  =  TEMAT3 !r
C        DATLNK(8+JNUM)  =  2.1E-20 !s = Zero value
C        DATLNK(9+JNUM)  =  TEMAT3!t
C        DATLNK(10+JNUM) =  2.1E-20 !theta r = Zero value
C        DATLNK(11+JNUM) =  TEMAT3 !theta s 
C        DATLNK(12+JNUM) =  2.1E-20  !theta t
                
        DATLNK(7+JNUM)  = TEMAT3 !r
        DATLNK(8+JNUM)  = 2.1E-20 !s = Zero value
        DATLNK(9+JNUM)  = TEMAT3 !t
        DATLNK(10+JNUM) = 2.1E-20 !theta r = Zero value
        DATLNK(11+JNUM) = TEMAT3 !theta s 
        DATLNK(12+JNUM) = 2.1E-20  !theta t
        
        
c       DATLNK(7+JNUM)  =  0 !r
c        DATLNK(8+JNUM)  =  0  !s = Zero value
c        DATLNK(9+JNUM)  =  0     !t
c        DATLNK(10+JNUM) =  0  !theta r = Zero value
c        DATLNK(11+JNUM) =  0 !theta s 
c        DATLNK(12+JNUM) =  0     !theta t
       
      ENDDO
      !---------------------------------
      
      ! MAKE NCON,LEKNOD MATRIX --------------------------------------------
C        = LEKNOD MATRIX (NODE-i, NODE-j, ISET)
      
      IL = 1
      DO INUMF = 1,NUMF
          
        INUM = NUMN(INUMF)
C        IL = INUM*(INUMF-1) + IL
        DO INU = 1,INUM
C            JNU = JNED+(INUM*(INUMF-1))+(INU-1)
            JNU = JNED+(INUMF-1)+(NUMF*(INU-1))
            LNODE(IL) = IDIAPH(JNU)  ! SAVE THE NODE IN EACH FLOOR
            ICNODE(IL) = NSN0 + INUMF ! SAVE THE CENTER NODE IN EACH FLOOR
            
C            LNODE(IL+INUM) = IDIAPH(JNU)  ! SAVE THE NODE IN EACH FLOOR
C            ICNODE(IL+INUM) = NSN0 + INUMF ! SAVE THE CENTER NODE IN EACH FLOOR             

            LNODE(IL+NGLELM0) = IDIAPH(JNU)  ! SAVE THE NODE IN EACH FLOOR
            ICNODE(IL+NGLELM0) = NSN0 + INUMF ! SAVE THE CENTER NODE IN EACH FLOOR            
            
            IL = IL + 1
        ENDDO   
            NCNODE(INUMF) = ICNODE(IL-1) 
        ENDDO  
      
        
      DO INGLE = 1,NGLELM
         JNUM = 3*(INGLE-1)
         
         NGLELM2 = NGLELM/2
         
         IF(INGLE.LT.NGLELM2.OR.INGLE.EQ.NGLELM2) THEN
           LEKNOD(1+JNUM) = LNODE(INGLE) !LEKNOD(FIRST COLUMN) 
           LEKNOD(2+JNUM) = ICNODE(INGLE) !LEKNOD(SECOND COLUMN)
         ELSEIF(INGLE.GT.NGLELM2) THEN
           LEKNOD(1+JNUM) = ICNODE(INGLE) !LEKNOD(FIRST COLUMN) 
           LEKNOD(2+JNUM) = LNODE(INGLE)  !LEKNOD(SECOND COLUMN)             
      
CC         LEKNOD(1+(NGLELM*3)+JNUM) = ICNODE(INGLE) !LEKNOD(SECOND COLUMN)
CC         LEKNOD(2+(NGLELM*3)+JNUM) = LNODE(INGLE) !LEKNOD(FIRST COLUMN) 
         ENDIF
         
        DO I = 1,2
	      IN = LEKNOD(I+JNUM)
	      NCON(IN,1) = INGLE
	  ENDDO
      ENDDO
      
	!------------------------------------------------ 2014.08.05 FINISHED.              
	
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
      
      DO INUMF = 1,NUMF
         IDSP(1:NSF,NSN0+INUMF) = IDSP(1:NSF,NSN0)
      ENDDO
      
      
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
      SUBROUTINE NODEIN2 (ID,SD,DFACS,IDSET,DIRCOS,XY2,MSF,MSN,SD2,
	1				   LCON1,LCON2)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     GENERATION OF STRUCTURE NODES AND BOUNDARY CONDITIONS
C	-----------------------
C     ID(NSF,NSN)   = BOUNDARY CONDITIONS
C     IDSET(NSN)    = LOCAL AXES SET NUMBERS
C     DIRCOS(9,NLA) = DIRECTION COSINES OF LOCAL AXES
C     XY(NSN,NSC)   = STRUCTURE COORDINATES
C	DFACS(NSN)    = ARRAY CONTAINING NODAL ISOSLATOR DAMPING/STIFF FACTOR
C     -----------------------------------------------------
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /NUMB2/ IDOFB(9)  
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	ADDED BY DE SILVA NOV2003
	COMMON /SPBC/ NSS,NLSS

C	==================================================================
C	JOINT GLOBAL CONSTRAINT SONGSAK JAN2007 
	COMMON /LC007/ LC107,LC207,LC307,LC407
C	==================================================================
     
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
      DIMENSION XY2(NSNL,1)
C      ALLOCATABLE XY2(:,:)
      
C	--------------------------------------------------------------
      DIMENSION ID(MSF,1),IDSET(1),DIRCOS(9,1),XY(MSN,NSC),COORD(3,21)
      
      
C	CHANGED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
C	DIMENSION XYZ(3),ISD(3),IBC(7),IPERM(22),DCOS(9)
	DIMENSION XYZ(3),ISD(3),IBC(9),IPERM(22),DCOS(9)

C	NEXT ADDED BY VSAYAK  11AUG03
	DIMENSION NOD(NSN) !Node number

C	ADDED BY DE SILVA NOV2003
C	CHANGED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
C	DIMENSION SD(6,MSN),SBC(6),STL(6,MSN),STC(6),VR(3),VS(3),VT(3)
	DIMENSION SD(6,MSN),SBC(6),STC(NSF),VR(3),VS(3),VT(3)

	DIMENSION TT1(6,6),DFACS(MSN)

C	BILINEAR SPRING SUPPORT SONGSAK NOV2006
	DIMENSION SD2(7*6,MSN)  !ACTUAL SIZE IS SD2(7*6,NSS)

	DIMENSION LCON1(11,LC107),LCON2(LC207,LC107)
	
C	--------------------------------------------------------------
      DATA IPERM  /1,4*0,1,0,1,0,0,1,1,0,5*1,4*0/

CC      NSN = NSN + NFLOOR
C     --------------
C     INITIALISATION
C     --------------
      CALL CLEARA (XY,NSN*NSC)
      CALL CLEARI (IDSET,NSN)
C     ----------------------
C     DEFINE LOCAL AXES SETS
C     ----------------------
      IF (NLA.EQ.0) GOTO 50
      ISET = 0
      WRITE (ITO,6000) NLA
      WRITE (10,6000) NLA

C	REIMPLIMENTATION BY DE SILVA 29/3/2004
c	read DIRECTION COSINES ***********************************
	READ (ITI,*)
 630  READ (ITI,*,END=631) NSET,(DCOS(N),N=1,6)
 631	DCOS(7) = DCOS(2)*DCOS(6)-DCOS(3)*DCOS(5)
      DCOS(8) = DCOS(3)*DCOS(4)-DCOS(1)*DCOS(6)
      DCOS(9) = DCOS(1)*DCOS(5)-DCOS(2)*DCOS(4)
      IF (NSET.LE.NLA) GOTO 610
      CALL ERRORS (21,ISET,NSET,'LOCAL AXES')
      GOTO 600
 610  UNITL = DCOS(7)*DCOS(7)+DCOS(8)*DCOS(8)+DCOS(9)*DCOS(9)
      IF (DABS(UNITL-1.0).LT.0.0001)  GOTO 620
      CALL ERRORS (22,ISET,NSET,'LOCAL AXES')
      GOTO 600
 620  DO 640  N=1,9
 640  DIRCOS(N,NSET) = DCOS(N)
      ISET = ISET + 1
 600  IF (ISET.LT.NLA) GOTO 630
C     -----------------------
C     PRINT DIRECTION COSINES
C     -----------------------
      WRITE (ISO,6100)
      DO 660 ISET=1,NLA
 660  WRITE (ISO,6200) ISET,(DIRCOS(N,ISET),N=1,9)


C     ----------------------------------------------
C     LOOP OVER NUMBER OF REGIONS, READ REGION CODE,
C     NODEL COORDINATES AND NUMBER OF INTERVALS
C     ----------------------------------------------
  50  NFIRST = 1
      IRE    = 0
 100  WRITE (ITO,1000) IRE+1
      WRITE (10,1000) IRE+1

	READ (ITI,*)
      READ (ITI,*) ICODE
C     ---------------------
C     READ INDIVIDUAL NODES
C     ---------------------
 200  WRITE (ITO,2000)
      WRITE (10,2000)

	READ (ITI,*)
      READ (ITI,*) NNO
	NLAST = NFIRST + NNO - 1
      IF (NLAST.LE.NSN) GOTO 210
      CALL ERRORS (2,IRE,NLAST,'INDIV.NODE')

210	READ (ITI,*)
	DO INO = NFIRST,NLAST
	    READ(ITI,*) NOD(INO),XY(INO,1:NSC)
      ENDDO
      
C     -------------------------------------------------------
C     END OF REGION LOOP, CHECK WHETHER ALL NODES ARE DEFINED
C     -------------------------------------------------------
 280  NFIRST = NLAST+1
      IRE    = IRE+1
 290  IF (IRE.LT.NRE) GOTO 100
      NLAST = NSN
      IF (NLAST.LT.NSN) CALL ERRORS (19,0,NLAST,'COORDINATE')

      !SAVE COORDINATE FROM XY MATRIX INTO XY2 BY BJ      
       XY2(1:MSN-1,1:NSC) = XY(1:MSN-1,1:NSC)
       XY2(MSN:NSNL,1:NSC) = 0.0

C     -----------------------------------------------------
C     -------------------------
C     GLOBAL SUPPORT CONDITIONS
C     -------------------------
	READ (ITI,*)
	READ (ITI,*)

      CALL CLEARI (ID,NSF*NSN)
      
      IF (NBS.EQ.0) GOTO 395
      IF (MODEX.EQ.3) GOTO 500
      WRITE (ITO,3000) NSF,NBS
      WRITE (10,3000) NSF,NBS

      NN  =  0
 300	READ (ITI,*,END=307) NODEA,(IBC(ISF),ISF=1,NSF)
 307	NODEA = IABS(NODEA)
      NODEB = NODEA
	INC = 1

      IF (NODEA.LE.NSN) GOTO 310
      CALL ERRORS (2,NN,NODEA,'SUPPORT C.')
      GOTO 390
 310  IF (NODEB.LE.NSN) GOTO 320
      CALL ERRORS (2,NN,NODEB,'SUPPORT C.')
      GOTO 390
 320  CONTINUE
C
      IF (INC.LE.0)  INC = 1
      DO 350  ISN=NODEA,NODEB,INC
      NN = NN+1
      DO 350  ISF=1,NSF
 350  ID(ISF,ISN) = ID(ISF,ISN) + IBC(ISF)
 390  IF (NN.LT.NBS) GOTO 300

C     -------------------------
C     SKEWED SUPPORT CONDITIONS
C     -------------------------
C 395  IF (NLS.EQ.0) GOTO 500
C	THIS LINE CHANGED BY DE SILVA NOV2003
  395  IF (NLS.EQ.0) GOTO 1036
      WRITE (ITO,4000) NSF,NLS
      WRITE (10,4000) NSF,NLS
C
      NN = 0
 400	READ (ITI,*,END=408) NODEA,NSET,(IBC(ISF),ISF=1,NSF)
 408	NODEA = IABS(NODEA)
      NODEB = NODEA
	INC = 1

      IF (NODEA.LE.NSN) GOTO 410
      CALL ERRORS (2,NN,NODEA,'LOCAL SUP.')
      GOTO 490
 410  IF (NODEB.LE.NSN) GOTO 420
      CALL ERRORS (2,NN,NODEB,'LOCAL SUP.')
      GOTO 490
 420  IF (NSET.LE.NLA) GOTO 425
      CALL ERRORS (21,NN,NSET,'LOCAL SUP.')
      GOTO 490
 425  CONTINUE

      IF (INC.LE.0)  INC = 1
      DO 450  ISN=NODEA,NODEB,INC
      NN = NN+1
      DO 460  ISF=1,NSF
	IF ((ID(ISF,ISN).EQ.0).OR.(IDSET(ISN).EQ.NSET)) GOTO 460
      IF (ISF.EQ.NSF .AND. IDOF(NSF).GE.7) GOTO 460
      CALL ERRORS (23,NN,ISN,'LOCAL SUP.')
      GOTO 490
 460  ID(ISF,ISN) = ID(ISF,ISN) + IBC(ISF)
 450  IDSET(ISN) = NSET
 490  IF (NN.LT.NLS) GOTO 400

C
C	ADDED BY DE SILVA NOV2003
C     ---------------------------------------------------
C	SPRING SUPPORT CONDITIONS ADDED BY DE SILVA NOV2003
C     ---------------------------------------------------
1036  CALL CLEARA (SD,6*NSN)
C	Next line added 12Dec03 by NguyenDV DAMPING FACTOR
	CALL CLEARA (DFACS,NSN)
      IF (NSS.EQ.0) GOTO 1395

      NN  =  0
 1300	READ (ITI,*,END=1307) NODEA,DFAC,(SBC(ISF),ISF=1,6),IOPT2  !GiD input
 1307	NODEA = IABS(NODEA)
      NODEB = NODEA
	INC = 1

      IF (NODEA.LE.NSN) GOTO 1310
      CALL ERRORS (2,NN,NODEA,'SUPPORT C.')
      GOTO 1390
 1310 IF (NODEB.LE.NSN) GOTO 1320
      CALL ERRORS (2,NN,NODEB,'SUPPORT C.')
      GOTO 1390
 1320 CONTINUE
C
      IF (INC.LE.0)  INC = 1

      DO 1350  ISN=NODEA,NODEB,INC
		NN = NN+1
		DFACS(ISN) = DFAC
		DO 1350  ISF=1,6
			SD(ISF,ISN) = SD(ISF,ISN) + SBC(ISF)
 1350 CONTINUE

C	ADDED BY SONGSAK FOR BILINEAR BASE ISOLATION NOV2006
	CALL INPSD2(SD2(1,NN),NODEA,IOPT2,1)
	IF(IOPT2.GT.0) THEN
	CALL INPSD2(SD2(1,NN),NODEA,IOPT2,2)
	ENDIF

 1390 IF (NN.LT.NSS) GOTO 1300
C	-----------------------------------------
C	SKIWED SPRING SUPPORT BOUNDARY CONDITIONS (Local!)
C	-----------------------------------------
 1395 IF (NLSS.EQ.0) GOTO 2036

      NN = 0
 1400	READ (ITI,*,END=1408) NODEA,NSET,(SBC(ISF),ISF=1,6),IOPT2  !GiD input
 1408	NODEA = IABS(NODEA)
      NODEB = NODEA
	INC = 1

      IF (NODEA.LE.NSN) GOTO 1410
      CALL ERRORS (2,NN,NODEA,'LOCAL SUP.')
      GOTO 1490
 1410 IF (NODEB.LE.NSN) GOTO 1420
      CALL ERRORS (2,NN,NODEB,'LOCAL SUP.')
      GOTO 1490
 1420 IF (NSET.LE.NLA) GOTO 1425
      CALL ERRORS (21,NN,NSET,'LOCAL SUP.')
      GOTO 1490
 1425 CONTINUE

C	CHANGE THE SBC MATRIX FROM LOCAL TO GLOBAL
C	----------------------
C	SET THE LOCAL VECTORS
	DO I = 1,3
	VR(I) = DIRCOS(I,NSET)
	VS(I) = DIRCOS(I+3,NSET)
	VT(I) = DIRCOS(I+6,NSET)
	END DO
	CALL TT1B (VR,VS,VT,TT1)
	SBC = MATMUL(TT1,SBC)

      IF (INC.LE.0)  INC = 1
      DO 1460  ISN=NODEA,NODEB,INC
      NN = NN+1
      DO 1460  ISF=1,6
 1460 SD(ISF,ISN) = SD(ISF,ISN) + SBC(ISF)

C	ADDED BY SONGSAK FOR BILINEAR BASE ISOLATION NOV2006
	CALL INPSD2(SD2(1,NN+NSS),NODEA,IOPT2,1)
	IF(IOPT2.GT.0) THEN
	CALL INPSD2(SD2(1,NN+NSS),NODEA,IOPT2,2)
	ENDIF

 1490 IF (NN.LT.NLSS) GOTO 1400

C	----------------------------------------
C	GLOBAL JOINT CONSTRAINT SONGSAK JAN2007
C	----------------------------------------
2036  IF(LC107.EQ.0) GOTO 500

	NCJX = 0
	READ(ITI,*)
3002	NCJX = NCJX + 1
	READ(ITI,*) LCON1(1:11,NCJX)
	LCJX = LCON1(2,NCJX)
	READ(ITI,*) LCON2(1:LCJX,NCJX) 

	IF (NCJX.LT.LC107) GOTO 3002
C     ----------------------------------------
C     PRINT COORDINATES AND SUPPORT CONDITIONS
C     ----------------------------------------
 500  WRITE (ISO,5000)
	
	WRITE(112,*) 'NUMBER OF NODE'    !ADDED BY SONGSAK MAR2006 3D
	WRITE(112,6300) NSN              !ADDED BY SONGSAK MAR2006 3D
	WRITE(112,*) 'NODAL COORDINATES' !ADDED BY SONGSAK MAR2006 3D

      CALL CLEARA (XYZ,3)
      DO 510  I=1,9                    !CHANGED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR) FROM 7 TO 9 DOFS
 510  IBC(I) = 1
      DO 590  ISN=1,NSN
      DO 550  ISC=1,NSC
 550  XYZ(ISC) = XY(ISN,ISC)
      DO 560  ISF=1,NSF
      IPO = IDOF(ISF)
 560  IBC(IPO) = ID(ISF,ISN)
      IF (NLS.EQ.0) GOTO 580
      NSET = IDSET(ISN)
      IF (NSET.EQ.0) GOTO 580
      WRITE (ISO,5200) ISN,XYZ,NSET,IBC
	WRITE (112,5100) ISN,XYZ         !ADDED BY SONGSAK MAR2006 3D
      GOTO 590
 580  WRITE (ISO,5100) ISN,XYZ,IBC
	WRITE (112,5100) ISN,XYZ         !ADDED BY SONGSAK MAR2006 3D
 590  CONTINUE

      DO ISN = 1,NSN 
      IDLC = IDSET(ISN)
      DO ISF = 1,NSF
      IF(IDLC.EQ.0) ID(ISF,ISN) = ID(ISF,ISN) + IDOFB(ISF) 
      IF(IDLC.NE.0.AND.ISF.GT.6) ID(ISF,ISN) = ID(ISF,ISN) + IDOFB(ISF) 
      ENDDO  
      ENDDO
      
C
      RETURN
C
 1000 FORMAT (1X,'READ REGION CODE FOR GENERATION OF STRUCTURE',
	1		1X,'COORDINATES; FOR REGION',I5)
 1100 FORMAT (1X,'READ ',I1,' COORDINATES FOR EACH',I3,' REGION NODES',
	1		13X,'; FOR REGION',I5)
 1200 FORMAT (1X,'READ NUMBER OF DIVISIONS IN ',I2,' DIRECTIONS',
	1		15X,'; FOR REGION',I5)
 1500 FORMAT (26X,28H* INADMISSIBLE REGION CODE *//
     1        1X,20HPOSSIBLE OPTIONS ARE/
     2        1X,50HINDIVIDUAL NODES  . . . . . . . . ICODE = -1/1    /
     3        1X,50HSTRAIGHT LINE . . . . . . . . . . ICODE = -2/2    /
     4        1X,50HPARABOLIC/CIRCULAR LINE . . . . . ICODE = -3/3    /
     5        1X,50H2-D STRAIGHT SIDED QUADRILATERAL  ICODE = -4/4    /
     6        1X,50H2-D PARABOLIC QUADRILATERAL . . . ICODE = -5 TO -8/
     7        1X,50H2-D CIRCULAR IN 1 DIRECTION . . . ICODE =  6      /
     8        1X,50H2-D CIRCULAR IN 2 DIRECTIONS  . . ICODE =  9      /
     9        1X,50H3-D STRAIGHT SIDED HEXAHEDRON . . ICODE =  8      /
     1        1X,50H3-D PARABOLIC HEXAHEDRON  . . . . ICODE = -9 TO-21/
     2        1X,50H3-D CIRCULAR IN 1 DIRECTION . . . ICODE = 12      /
     3        1X,50H3-D CIRCULAR IN 2 DIRECTIONS  . . ICODE = 18      )
 2000 FORMAT (1X'READ NUMBER OF INDIVIDUAL NODES')
 2100 FORMAT (I6/I6/(3E20.10))
 3000 FORMAT (1X,'READ GLOBAL SUPPORT CONDITIONS; ',I1,1X,6HDOF'S,,I5,
	1		1X,'SUPPORTS')
 4000 FORMAT (1X,'READ LOCAL  SUPPORT CONDITIONS; ',I1,1X,6HDOF'S,,I5,
	1		1X,'SUPPORTS')
C	===================================================
C	CHANGED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
 5000 FORMAT (//,28X,23(1H*)/28X,1H*,21X,1H*/
     1        28X,23H*  NODAL COORDINATES  */
     2        28X,23H* SUPPORT  CONDITIONS */
     3        28X,1H*,21X,1H*/28X,23(1H*)//
     4        5H NODE,2X,12HX-COORDINATE,2X,12HY-COORDINATE,2X
     5        ,12HZ-COORDINATE,3X,6HSET NO,5X,
     6         34HDX  DY  DZ  RX  RY  RZ  Wp  T.  P./1X,94(1H-))
 5100 FORMAT (I5,3E14.5,12X,9I4)
 5200 FORMAT (I5,3E14.5,8H  LOCAL ,I2,2X,9I4)
C	===================================================
 6000 FORMAT (1X,'READ DIRECTION COSINES OF',I5,' LOCAL AXES SETS')
 6100 FORMAT (//30X,19(1H*)/,30X,1H*,17X,1H*,/
     1        30X,19H* LOCAL AXES SETS */,30X,1H*,17X,1H*/30X,19(1H*)//
     2        34X,12HLOCAL  1-2-3,/34X,12HGLOBAL X-Y-Z,//
     3        1X,6HSET NO,4X,3H1-X,5X,3H1-Y,5X,3H1-Z,5X,3H2-X,5X,3H2-Y
     4        ,5X,3H2-Z,5X,3H3-X,5X,3H3-Y,5X,3H3-Z,/,1X,79(1H-)/)
 6200 FORMAT (I5,3X,9F8.4)
 6300 FORMAT (X,I8)
C
      END
C
C=====================================================================
      SUBROUTINE NODEADD (XY2,CXY)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     --------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM      
C     ----------------------------------
C     ADD THE CENTER NODE OF EACH FLOOR
C     ----------------------------------
      COMMON /NUMFLO/ NUMF,NUMVT

C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
      DIMENSION XY2(NSNL,1),CXY(NUMF,3)

      DO INUMF = 1,NUMF
        CX = CXY(INUMF,1)
        CY = CXY(INUMF,2)
        CZ = CXY(INUMF,3)
       
        XY2(NSN0+INUMF,1) = CX
        XY2(NSN0+INUMF,2) = CY
        XY2(NSN0+INUMF,3) = CZ
      ENDDO
      

      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE INCENLO
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     READ DIAPHRAGM CONSTRAINTS DATA BY BJ
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NEDATA        = SAVE NODE DATA IN EACH FLOORF
C     ILFLO         = THE NUMBER OF LOAD     
C     ICLO          = THE NUMBER OF CENTER LOAD DATA
C     ILDEC         = DIRECTION NUMBER OF CENTER LOAD
C     CLOAD         = VALUE OF CENTER LOAD
C     -----------------------------------------------------
C   
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
     
     
      COMMON /LOADCEN/ ICLOAD,ICLO ! POINT LOAD IN CENTROID OF EACH FLOOR
      
      COMMON /DAPH/ DIAPH(700000),IDIAPH(600000) ! STORAGE FOR CONSTRAINTS BY BJ
        
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD,JLCS ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ  
       
        READ (ITI,*)     
C        READ (ITI,*) ICLO
        
C        IF(ICLO.EQ.0) THEN
C            JCENL = 0
C            RETURN
C        ENDIF 
        
        JCENL = 1
        
        CALL UPDIA(JILF,JNED,JNEDD)
        
        JILFD = ICLO
        
        CALL UPDIA(JILD,JILF,JILFD)
        
        JILDD = ICLO
        
       CALL UPDIA(JLCS,JILD,JILDD)        
        
        DO ICL = 1,ICLO
        
        CALL READCENLO(IDIAPH(JILF),IDIAPH(JILD),IDIAPH(JLCS),DIAPH(JCENL),ICL)

        ENDDO
         
        JCENLD = ICLO
        
C        JLCSD = ICLO
C      CALL UPDIA(   ,JLCS,JLCSD)         

C      CALL CALCENTRO
C       CALL CALCOORD
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	======================================================================= 
      SUBROUTINE READCENLO(ILFLO,ILDEC,LCSCEN,CENLOAD,ICL)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     READ DIAPHRAGM CONSTRAINTS DATA BY BJ
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NEDATA        = SAVE NODE DATA IN EACH FLOORF
C     ILFLO         = THE NUMBER OF LOAD     
C     ICLO          = THE NUMBER OF CENTER LOAD DATA
C     ILDEC         = DIRECTION NUMBER OF CENTER LOAD
C     CLOAD         = VALUE OF CENTER LOAD
C     LCSCEN        = LOADCASE FOR CENTER LOAD
C     -----------------------------------------------------
C   
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
     
        
      COMMON /LOADCEN/ ICLOAD,ICLO ! POINT LOAD IN CENTROID OF EACH FLOOR     
C      COMMON /LOADFLOOR/ ICLO,ILFLO(1000),ILDEC(1000),CENLOAD(1000) 
      DIMENSION ILFLO(ICLO),ILDEC(ICLO),CENLOAD(ICLO) ,LCSCEN(ICLO)
        
        READ (ITI,*) JILFLO,JILDEC,CENLOAD2,JLCSCEN 
        
         ILFLO(ICL) = JILFLO
         ILDEC(ICL) = JILDEC
         CENLOAD(ICL) = CENLOAD2
         LCSCEN(ICL) = JLCSCEN
        
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================   
      SUBROUTINE NODELO2 (ID,MSF,NNL,ICLO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     READS AND GENERATES NODAL LOADS, ADDS CONTRIBUTIONS TO LOAD V
C     -------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

	COMMON /LCSS/ ILCN,ILCC
	
	COMMON /CNODE/ NCNODE(5000) !SAVE CENTER NODE NUMBER
      
C     FOR DIAPHRAGM CONSTRATINS BY BJ
      COMMON /JLCN/ ILCNM(1000),IL2      
C
      DIMENSION ID(MSF,1),R(NEQ),TNO(NSN),FLOAD(MSF,NSN),FORCE(9)


      IF(NNL.EQ.0) THEN
C       NNL = ICLO
      GO TO 180 
      ENDIF
      
C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
 150  WRITE (6,1100) NNL

	TNO = 0.0

C     ---------------
C     INITIALISATION
C     --------------
 180  CALL CLEARA (R,NEQ)
      CALL CLEARA (FLOAD,NSF*NSN)
      
      IF(NNL.GT.0.AND.ICLO.GT.0) THEN
        NN0 = NNL
        NNL0 = NNL
        NNL  = NNL + ICLO
      ELSEIF(NNL.EQ.0) THEN
        NN0 = NNL
        NNL = ICLO
        NNL0 = ICLO
      ENDIF  
     

      INL = 0
      
C     FOR DIAPHRAGM CONSTRATINS BY BJ
      IL = 1
      JCN = 0
      ILCNM(1:1000) = 0.0
C     -------------------------------- 
      

      READ (ITI,*)
      
 190 	IF(INL.EQ.NNL0) READ (ITI,*)

C 200	READ (ITI,*,END=202) NODEA,ISF,RLOAD,ILCN,ILCC
200   IF(NN0.NE.0.AND.INL.LT.NNL0) THEN
        READ (ITI,*,END=202) NODEA,ISF,RLOAD,ILCN,ILCC
      ELSE
        READ (ITI,*,END=202) NFLOOR,ISF,RLOAD,ILCN,ILCC
        !READ CENTER NODE NUMBER
        NODEA = NCNODE(NFLOOR) 
      ENDIF    
CC 200	READ (ITI,*,END=202) NFLOOR,ISF,RLOAD,ILCN,ILCC
      
      

 202  NODEA = IABS(NODEA)
      NODEB = NODEA
	INC = 1

C     FOR DIAPHRAGM CONSTRATINS BY BJ
C      JLCN = ILCNM(IL)
      IF(JLCN.NE.ILCN) THEN
       ILCNM(IL) = ILCN
       IL = IL + 1
       JLCN = ILCN
      ENDIF
C     --------------------------------        
      
	ISF = IDOFCALL(IDOF,ISF)

      IF (NODEA.LE.NSN) GOTO 210
      CALL ERRORS (2,INL,NODEA,'JOINT LOAD')
      GOTO 290
 210  IF (NODEB.LE.NSN) GOTO 220
      CALL ERRORS (2,INL,NODEB,'JOINT LOAD')
      GOTO 290
C

220	 DO 250  INO=NODEA,NODEB,INC
		TNO(INO) = 1.0
		FLOAD(ISF,INO) = FLOAD(ISF,INO) + RLOAD
		IEQ = ID(ISF,INO)
		INL = INL+1
		IF (IEQ) 240,240,230
 230		R(IEQ) = R(IEQ) + RLOAD
		GOTO 250
 240		CALL ERRORS (7,INL,INO,'JOINT LOAD')
 250    CONTINUE

	CALL LDASEM (R)


	IF (INL.GT.NNL) CALL ERRORS (24,INL,NNL,'JOINT LOAD')
CC 290  IF (INL.LT.NNL) GOTO 200
 290  IF (INL.LT.NNL) GOTO 190
C     -----------------
C     PRINT NODAL LOADS
C     -----------------
C	FORCES
C      WRITE (ISO,2100)
C      CALL CLEARA (FORCE,6)
C      DO 490  ISN=1,NSN
C      IF (TNO(ISN).EQ.0.) GOTO 490
C      DO 450  ISF=1,NSF
C      IPO = IDOF(ISF)
C 450  FORCE(IPO) = FLOAD(ISF,ISN)
C      WRITE (ISO,4000) ISN,FORCE
C 490  CONTINUE
      
      WRITE (ISO,2100)
      CALL CLEARA (FORCE,9)
      
       DO ISN=1,NSN
         DO ISF=1,NSF
          IPO = IDOF(ISF)
          FORCE(IPO) = FLOAD(ISF,ISN)
         ENDDO
          SUMF = SUM(ABS(FORCE(1:9)))
          IF(SUMF.GT.0.0) WRITE (ISO,4000) ISN,FORCE
490   ENDDO
C
      
C     FOR DIAPHRAGM CONSTRATINS BY BJ
      IL2 = IL-1
C     --------------------------------         

C
      RETURN
C
 1100 FORMAT (1X,'READ & GENERATE NODAL LOADS;',I5,' NODAL LOADS')
 2100 FORMAT (///28X,24HCONCENTRATED NODAL LOADS//
     +        1X,'JOINT',8X,'FORCE-X',6X,'FORCE-Y',6X,'FORCE-Z',5X,
     +        'MOMENT-X',5X,'MOMENT-Y',5X,'MOMENT-Z',5X,'WARPING ',5X,'  TEMP  ',4X,'PORE PRS'/1X,80(1H-))     
 2200 FORMAT (///30X,19HNODAL DISPLACEMENTS//
     +        1X,'JOINT',8X,'DISPL-X',5X,'DISPL-Y',5X,'DISPL-Z',4X,
     +        'ROTATI-X',4X,'ROTATI-Y',4X,'ROTATI-Z'/1X,80(1H-))
 4000 FORMAT (1X,I5,3X,9E13.5)
C
      END
C
C=====================================================================
      SUBROUTINE NODELO22 (ID,MSF,ICLO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     READS AND GENERATES NODAL LOADS, ADDS CONTRIBUTIONS TO LOAD V
C     -------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

	COMMON /LCSS/ ILCN,ILCC
	
	COMMON /CNODE/ NCNODE(5000) !SAVE CENTER NODE NUMBER
C
      DIMENSION ID(MSF,1),R(NEQ),TNO(NSN),FLOAD(MSF,NSN),FORCE(9)

      
	TNO = 0.0

C     ---------------
C     INITIALISATION
C     --------------
C 180  CALL CLEARA (R,NEQ)
C      CALL CLEARA (FLOAD,NSF*NSN)
      
      NNL = ICLO

      INL = 0

 	READ (ITI,*)

C 200	READ (ITI,*,END=202) NODEA,ISF,RLOAD,ILCN,ILCC
 200	READ (ITI,*,END=202) NFLOOR,ISF,RLOAD,ILCN,ILCC
      
      !READ CENTER NODE NUMBER
      NODEA = NCNODE(NFLOOR) 
      
 202  NODEA = IABS(NODEA)
      NODEB = NODEA
	INC = 1

	ISF = IDOFCALL(IDOF,ISF)

      IF (NODEA.LE.NSN) GOTO 210
      CALL ERRORS (2,INL,NODEA,'JOINT LOAD')
      GOTO 290
 210  IF (NODEB.LE.NSN) GOTO 220
      CALL ERRORS (2,INL,NODEB,'JOINT LOAD')
      GOTO 290
C

220	 DO 250  INO=NODEA,NODEB,INC
		TNO(INO) = 1.0
		FLOAD(ISF,INO) = FLOAD(ISF,INO) + RLOAD
		IEQ = ID(ISF,INO)
		INL = INL+1
		IF (IEQ) 240,240,230
 230		R(IEQ) = R(IEQ) + RLOAD
		GOTO 250
 240		CALL ERRORS (7,INL,INO,'JOINT LOAD')
 250    CONTINUE

	CALL LDASEM (R)


	IF (INL.GT.NNL) CALL ERRORS (24,INL,NNL,'JOINT LOAD')
 290  IF (INL.LT.NNL) GOTO 200
C     -----------------
C     PRINT NODAL LOADS
C     -----------------
C	FORCES
      WRITE (ISO,2100)
      CALL CLEARA (FORCE,9)
      DO 490  ISN=1,NSN
      IF (TNO(ISN).EQ.0.) GOTO 490
      DO 450  ISF=1,NSF
      IPO = IDOF(ISF)
 450  FORCE(IPO) = FLOAD(ISF,ISN)
      WRITE (ISO,4000) ISN,FORCE
 490  CONTINUE
C
C
      RETURN
C
           
 1100 FORMAT (1X,'READ & GENERATE NODAL LOADS;',I5,' NODAL LOADS')
 2100 FORMAT (///28X,24HCONCENTRATED NODAL LOADS//
     +        1X,'JOINT',8X,'FORCE-X',6X,'FORCE-Y',6X,'FORCE-Z',5X,
     +        'MOMENT-X',5X,'MOMENT-Y',5X,'MOMENT-Z',5X,'WARPING ',5X,'  TEMP  ',4X,'PORE PRS'/1X,80(1H-))     
 2200 FORMAT (///30X,19HNODAL DISPLACEMENTS//
     +        1X,'JOINT',8X,'DISPL-X',5X,'DISPL-Y',5X,'DISPL-Z',4X,
     +        'ROTATI-X',4X,'ROTATI-Y',4X,'ROTATI-Z'/1X,80(1H-))
 4000 FORMAT (1X,I5,3X,9E13.5)

C
      END
C
C=====================================================================
      SUBROUTINE CENTERMASS(CXY,CXYMASS)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     READ LINK CONSTRAINTS DATA BY BJ      
C	-----------------------
C     NUMF          = THE NUMBER OF FLOOR
C     IFNUM         = FLOOR NUMBER
C     NUMN          = THE NUMBER OF NODE IN EACH FLOOR
C     NUMVT         = VERTICAL DIRECTION NUMBER OF FLOOR
C     NEDATA        = SAVE NODE DATA IN EACH FLOORF
C     ICLO          = THE NUMBER OF CENTER LOAD DATA
C     ILDEC         = DIRECTION NUMBER OF CENTER LOAD
C     CLOAD         = VALUE OF CENTER LOAD
C     -----------------------------------------------------
C   
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
     
C      COMMON /DIAPH/ NUMF,IFNUM(500),NUMN(1000),NUMVT,NEDATA(1000,1000)
      
      COMMON /DAPH/ DIAPH(700000),IDIAPH(600000) ! STORAGE FOR CONSTRAINTS BY BJ
      
      COMMON /IDCAT/ JCENL,JCENLD,JCO,JCOD,JCXY,JMCO ! REAL INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
      
      COMMON /IDCATI/ JIF,JNUM,JNED,JNEDD,JMAX,JILF,JILFD,JILD,JILDD,JLC ! INTEGER INDICATOR FOR DIAPHRAGM CONSTRAINTS BY BJ 
          
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      
      COMMON /NUMFLO/ NUMF,NUMVT
      COMMON /NLCASE/ NLCS,LCEQ(1000)
C      COMMON /IEQNUMFLO/ IEQNUM,ILOC
C      COMMON /EQPO/ POMAS  
        
C      DIMENSION IFNUM(1),NUMN(1),FNEDATA(1)            
C      DIMENSION ILOC(NUMF), POMAS(NUMF,2)
      DIMENSION CXY(NUMF,3),CXYMASS(NUMF,3)
      
      READ (ITI,*)
      READ (ITI,*) NLCS
      READ (ITI,*) LCEQ(1:NLCS)
    
      
      IT = NUMVT ! VERTICAL DIRECTION NUMBER OF FLOOR
        
        IF(IT.EQ.1) THEN
            IR = 2
            IS = 3
        ELSEIF(IT.EQ.2) THEN   
            IR = 1
            IS = 3             
        ELSEIF(IT.EQ.3) THEN   
            IR = 1
            IS = 2             
        ENDIF 
      
      DO INUMF = 1,NUMF
          CXYMASS(INUMF,1:3) = CXY(INUMF,1:3)
      ENDDO
      !FLOOR NUMBER, COORDINATES OF EACH DIRECTIONS(X,Z) IN EACH FLOOR
      DO INUMF = 1,NUMF
C       READ (ITI,*) ILOC(INUMF),POMAS(INUMF,1:2)
       READ (ITI,*) ILOC,POMASX,POMASZ
          CXYMASS(ILOC,IR) = POMASX
          CXYMASS(ILOC,IS) = POMASZ
      ENDDO
      
      
      RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================       
    