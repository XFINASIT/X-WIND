C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE DM3DM(IDATM,AA,BB,CC)  
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	------------------------------------------------------------
	CHARACTER*200 NAME
	CHARACTER*200 STROHED,TRNOHED
C	------------------------------------------------------------
      CHARACTER*80 TITLE
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL
     
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG, !Added 7Mar06
     +              ISOLV,IVPRT

	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
	
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS

      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG !added 5Nov2005

      COMMON A(9000000),IA(9000000)

	COMMON /MEMW/ W(7000000),IW(7000000)
C
	DIMENSION AA(1),BB(1),CC(1)
	
	DIMENSION NTRAIN(3)
	
	ALLOCATABLE  NEID(:), PELW(:)
      ALLOCATABLE DMASS(:),DDAMP(:),DSTIF(:)
      ALLOCATABLE VMASS(:),VDAMP(:),VSTIF(:)
      ALLOCATABLE MAXAM(:),WHPOS(:),LWHDOF(:),WHSLOD(:),WCSTIF(:)	
      
C     FOR TRAIN
      ALLOCATABLE PCAB(:,:),CADI(:,:),ADDI(:,:)  
	ALLOCATABLE PRSS(:,:),PBOG(:,:),ALSU(:,:)
	ALLOCATABLE PRPS(:,:),PRWH(:,:),WRCO(:,:) 
	ALLOCATABLE TRNDIM(:,:) 

C     OUTPUT CONTROL	      
      ALLOCATABLE KOST(:,:),STROHED(:),KOVH(:),TRNOHED(:)
      
C	------------------------------------------------------------------------
C	LOADWAY DEFINE
C	--------------
C	Read number of elements to define loaways of bridge
	CALL RENELW

	
C	Read element number array to define loaways of bridge
	ALLOCATE(NEID(2*NELW),PELW(4*NELW))
	CALL REELELW(NEID)

C	-------------------------------
C	RECALL MATRICES OF  STRUCTURES:
C	-------------------------------
	KSTEP = 1

      NUMITE = 0
      NUMREF = 0
      INDPD  = KPOSD

C     Form mass matrix:
C     -----------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 5
      CALL GRLOOP (IA(LEL),KSC)

C
C     Form damping matrix: 
C     --------------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 6
      CALL GRLOOP (IA(LEL),KSC)
C
C     Form tangential stiffness matrix: (IFREF=0)
C     ------------------------------------------
      IFEIG = 1
	IFREF = 0
      ISPRI = 1
      ITASK = 1
      CALL GRLOOP (IA(LEL),KSC)

      NUMREF = NUMREF + 1

C	Read additional parameters of bridge for 3D model (eccentricity)
      CALL  REBRI3D
C
C	FORM DAMPING MATRIX OF BRIDGE
C	-----------------------------
C	Solve for bridge frequencies & Calculate Rayleigh damping coefficients
C	By Subspace
	IF(ISOLV.EQ.1) THEN
		CALL STABIL(AA,BB,'STIF','MASS')
	  CALL MPFCAL(IA(LMA),IA(LID),A(LER),A(LEV),A(LDIM),
	1			    A(LFRE),NROOT,NITEM,BB,'MASS',KSC)
	    CALL RDAMCO(A(LEV),NROOT,NITEM)  
C     By Lanczos
	ELSE IF (ISOLV.EQ.2) THEN
	    CALL LANC(W,IA(LID),IA(LMA),N11,N10,AA,BB,'STIF','MASS')
	  CALL MPFCAL(IA(LMA),IA(LID),W(N10),W(N11),A(LDIM),
	1			   A(LFRE),NROOT,NITEM,BB,'MASS',KSC)
		CALL RDAMCO(W(N11),NROOT,NITEM) 
	ENDIF

C     Form total damping matrix (RDM OR/AND RDK .NE.0.)
	CALL RDAMBR (IA(LMA),RDM,RDK,NWK,NWM,NEQ,AA,BB,CC,
	1			 'STIF','MASS','DAMP')

C	READ PARAMETERS FOR RAIL/TRACK IRREGULARTIES AND GENERATE
	CALL REIRRE3D
	
C	----------------------------
C	FORMULATION OF TRAIN SYSTEM:
C	----------------------------
C	Read control variables for 3D train/truck model:
	CALL REDMASS(NDOFM,NWHEEL,LSTATC,IVTYPE,NTRAIN,NBOGI,IWRIN,NOUTST,NOUTVH,SPACING,NCOPY)

      ALLOCATE(MAXAM(NDOFM+1),WHPOS(NWHEEL),LWHDOF(3*NWHEEL),WHSLOD(NWHEEL),WCSTIF(3*2*NWHEEL))
      
C     ------------------      
      SELECTCASE(IVTYPE)
C     ------------------ 
      
      CASE(0) !DIRECT INPUT OF VEHICLE MATRICES
C	Generate vehicle matrices
      ALLOCATE(DMASS(NDOFM*NDOFM),DDAMP(NDOFM*NDOFM),DSTIF(NDOFM*NDOFM))
      
      CALL REDVMAT(NDOFM,NWHEEL,MAXAM,DMASS,DDAMP,DSTIF,WHPOS,LWHDOF,WHSLOD,WCSTIF,SPACING,NCOPY)
      
C     MAXIMUM STORAGE NEEDED FOR VEHICLE MATRIX
      MAXM = MAXAM(NDOFM+1) - 1
      ALLOCATE(VMASS(MAXM),VDAMP(MAXM),VSTIF(MAXM))
      VMASS(1:MAXM) = DMASS(1:MAXM)
      VDAMP(1:MAXM) = DDAMP(1:MAXM)
      VSTIF(1:MAXM) = DSTIF(1:MAXM)
      
      DEALLOCATE(DMASS,DDAMP,DSTIF)
      
      
      CASE(1) !TRAIN
      
C	FORM VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS OF TRAIN STIFFNESS/DAMPING MATRIX
	CALL MAXATRAIN(NTRAIN,NDOFM,MAXAM,LWHDOF)
	
C     MAXIMUM STORAGE NEEDED FOR TRAIN MATRIX
      MAXM = MAXAM(NDOFM+1) - 1
      ALLOCATE(VMASS(MAXM),VDAMP(MAXM),VSTIF(MAXM))      
      VMASS(1:MAXM) = 0.0
      VDAMP(1:MAXM) = 0.0
      VSTIF(1:MAXM) = 0.0
      
      
C     TOTAL NUMBER OF CAR BODY
      NCARB = NTRAIN(1) + NTRAIN(2) + NTRAIN(3)
C     NUMBER OF ARTIGULATED CAR BODY
      NACAB = NTRAIN(2)
      
      ALLOCATE( PCAB(NCARB,4) ,CADI(NCARB,7) ,ADDI(NACAB,5)  )
	ALLOCATE( PRSS(NBOGI,9) ,PBOG(NBOGI,4) ,ALSU(NACAB-1,6))
	ALLOCATE( PRPS(NWHEEL,9),PRWH(NWHEEL,6),WRCO(NWHEEL,6) )
	ALLOCATE( TRNDIM(3*8,NCARB+NBOGI+NWHEEL) )

C     READ TRAIN PROPERTIES	
      CALL TRAINPROP(VEL,NCARB,NACAB,NBOGI,NWHEEL,IWRIN,WHPOS,
	1               PCAB,CADI,ADDI,ALSU,PRSS,PBOG,PRPS,PRWH,WRCO)  
	
C	cal. static axle forces 
	CALL TRAINAXLEL(NTRAIN,NCARB,NBOGI,NWHEEL,PCAB,CADI,PBOG,PRWH,SCALE,WHSLOD)

C	cal. wheel-rail contac mechanism if required --> IWRIN=1
	IF(IWRIN.EQ.1) CALL RAILCONTCT(NWHEEL,PRWH,WHSLOD,WRCO,VEL) 

C     TRAIN MATRICES
      CALL TRAINMAT(NTRAIN,NCARB,NACAB,NBOGI,NWHEEL,VMASS,VDAMP,VSTIF,
     1              PCAB,CADI,ADDI,ALSU,PRSS,PBOG,PRPS,PRWH,WRCO)  

C     TRANSFER DIAGONAL LUMPED MASS INTO SKYLINE FORMAT (MAKe CONSISTANCE WITH ANOTHER MATRICES)
      ALLOCATE(DMASS(MAXM))      
          DMASS(1:MAXM) = VMASS(1:MAXM)
          VMASS(1:MAXM) = 0.0
          DO IDOFM = 1,NDOFM
            INDX = MAXAM(IDOFM)
            VMASS(INDX) = DMASS(IDOFM)
          ENDDO
      DEALLOCATE(DMASS)

C      CALL PRINSTIFC(VSTIF,MAXAM,NDOFM)
C      STOP
      
C     TRANSFER WHEEL RAIL CONTACT PROPETY TO WORKING MATRIX
      CALL CONTCTRAN(NWHEEL,PRWH,WRCO,WCSTIF)

C      CALL TRAINRENDER(NTRAIN,NCARB,NBOGI,NWHEEL,WHPOS,CADI,PRSS,PRPS,PRWH,TRNDIM,'NODE')         
C      CALL TRAINRENDER(NTRAIN,NCARB,NBOGI,NWHEEL,WHPOS,CADI,PRSS,PRPS,PRWH,TRNDIM,'ELEM')
     
C     ------------------ 
      ENDSELECT
C     ------------------ 


C     ------------------------------------------------------
C     READ OUTPUT CONTROL OF VEHICLE AND STRUCTURE    
C     ------------------------------------------------------     
      ALLOCATE(KOST(2,NOUTST),STROHED(NOUTST),KOVH(NOUTVH),TRNOHED(NOUTVH))
C     ------------------------------------------------------
C     VEHICLE
      READ(ITI,*)
      DO IOUTVH = 1,NOUTVH
	    READ(ITI,'(A200)') NAME
	    READ(ITI,*) II,KDF
             KOVH(IOUTVH) = KDF
          TRNOHED(IOUTVH) = TRIM(NAME)
      ENDDO
C     ------------------------------------------------------
C     STRUCTURE
      READ(ITI,*)
      DO IOUTST = 1,NOUTST
	    READ(ITI,'(A200)') NAME
	    READ(ITI,*) KSN,KDF
           KOST(1,IOUTST) = KSN
           KOST(2,IOUTST) = KDF
          STROHED(IOUTST) = TRIM(NAME)
      ENDDO
C     ------------------------------------------------------
C     ------------------------------------------------------
      
      
      
C     DIRECTION VECTOR OF EACH LANE SEGMENT, SEGMENT LENGTH AND TOTAL LENGTH
      CALL PATHMOV(NELW,NEID,PELW,A(LXY),TLENGTH)


C	-------------------------------------------------
C	SOLVE TRAIN-STRUCTURE INTERACTION IN 3D-MODELLING
C	-------------------------------------------------
      CALL VEHSOL3D(MAXAM,VMASS,VDAMP,VSTIF,WHPOS,LWHDOF,WHSLOD,WCSTIF,NOUTST,NOUTVH,
     1              LSTATC,NDOFM,NWHEEL,NEID,PELW,IA(LMA),KOST,STROHED,KOVH,TRNOHED,
     2              IVTYPE,TRNDIM,NTRAIN,NCARB,NBOGI,AA,BB,CC,'STIF','MASS','DAMP')		  
C	------------------------------------------------------------------------

      DEALLOCATE(KOST,STROHED,KOVH,TRNOHED)
      
	DEALLOCATE(NEID,PELW)
      DEALLOCATE(VMASS,VDAMP,VSTIF)
      DEALLOCATE(MAXAM,WHPOS,LWHDOF,WHSLOD,WCSTIF)

      IF(IVTYPE.EQ.1) THEN !TRAIN
        DEALLOCATE( PCAB,CADI,ADDI)
	  DEALLOCATE( PRSS,PBOG,ALSU)
	  DEALLOCATE( PRPS,PRWH,WRCO)
	  DEALLOCATE( TRNDIM)
	ENDIF
	
      RETURN
	END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE VEHSOL3D(MAXAM,VMASS,VDAMP,VSTIF,WHPOS,LWHDOF,WHSLOD,WCSTIF,NOUTST,NOUTVH,
     1                    LSTATC,NDOFM,NWHEEL,NEID,PELW,MAXA,KOST,STROHED,KOVH,TRNOHED,
     2                    IVTYPE,TRNDIM,NTRAIN,NCARB,NBOGI,AK,AM,AC,STIF,MASS,DAMP) 
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 STIF,MASS,DAMP
C	------------------------------------------------------------------------
	CHARACTER*200 NAME
	CHARACTER*200 STROHED,TRNOHED
C	------------------------------------------------------------------------
C     LSTATC = FLAG FOR STATIC ANALYSIS   0=STATIC 1=DYNAMIC
C	------------------------------------------------------------------------
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
	
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
     
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
	
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM
	
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS

      COMMON A(9000000),IA(9000000)
      
C	-------------------------------------------------------------------

      DIMENSION MAXA(1),AK(1),AM(1),AC(1)
      DIMENSION MAXAM(1),VMASS(1),VDAMP(1),VSTIF(1)
      DIMENSION WHPOS(1),LWHDOF(3,1),WHSLOD(1),WCSTIF(1)
      
      DIMENSION WHPOSN(NWHEEL),IDWHPOS(NWHEEL)  
      DIMENSION XI(NWHEEL)
      DIMENSION RSW(3*2*NWHEEL),FWR(3*NWHEEL)
      DIMENSION FC(NWHEEL),FCO(NWHEEL),FCMAX(2,NWHEEL)          !CONTACT FORCES
      DIMENSION FS(NEQ),RS(NEQ,3),RSN(NEQ,3)
      DIMENSION FT(NDOFM),RT(NDOFM,3),RTN(NDOFM,3)
      
      DIMENSION DDK(NEQ)
      
      DIMENSION RRX(NFFT),RRE(NFFT),RRA(NFFT),RRR(NFFT),TIRAIL(3,NWHEEL)
      
      DIMENSION NEID(1),PELW(1)

      DIMENSION KOST(2,NOUTST),STROHED(NOUTST),KOVH(NOUTVH),TRNOHED(NOUTVH)
      
      DIMENSION VOUTMAXS(6,NOUTST),VOUTMAXV(6,NOUTVH)
      
      DIMENSION NTRAIN(3),TRNDIM(3*8,1)
      
      ALLOCATABLE DMASS(:),DDAMP(:),DSTIF(:)
      
            
C	CALCULATE NUMBER OF SOLUTION STEPS NSOLT
c	Calculate the train length
	TRLENG = DABS(WHPOS(NWHEEL)-WHPOS(1))		
c	BIDIS can be (+): a distance before entering
c				   (-): the train starts within bridge spans
      SOLT=(BIDIS + TLENGTH + EXTDIS + TRLENG)/VEL
      NSOLT=IDNINT(SOLT/DELT)


C	Total time considered for rail irregularities
	TIRR=(TRLENG + BIDIS + TLENGTH + EXTDIS + TRLENG)/VEL
c	Generate 3-D rail irregularities by random functions (using FRA PSD),added 31Jan2007
	IF(IRIN.EQ.3.AND.IRAN.EQ.1)THEN
	  POSLW = WHPOS(NWHEEL) ! POSITION OF LAST WHEEL RELATIVE WITH FRONT WHEEL (MINUS SIGN)
		CALL RAILFRA(POSLW,TIRR,VEL,RRX,RRE,RRA,RRR)
	ENDIF


C     MAXIMUM STORAGE NEEDED FOR VEHICLE MATRIX
      MAXM = MAXAM(NDOFM+1) - 1
      ALLOCATE(DMASS(MAXM),DDAMP(MAXM),DSTIF(MAXM))
       
C	BACK UP VEHICLE MATRICES
	CALL MOVE(VMASS,DMASS,MAXM)
	CALL MOVE(VDAMP,DDAMP,MAXM)
	CALL MOVE(VSTIF,DSTIF,MAXM)

C	***********************************************
C 	*		      ZERO STEP                       *                            
C     ***********************************************
C	Asumming Responses of structure at nodal DOFs: RS(NEQ,3)[DIS,VEL,ACCE] 
	CALL CLEARMAT(RS,NEQ,3)

C	Initialize Responses of structure at wheel position
	CALL CLEARA(RSW,3*2*NWHEEL)  

C     Initialize contact forces
      CALL CLEARA(FC ,NWHEEL)	
      CALL CLEARA(FCO,NWHEEL)
      
C     Initialize Responses of vehicle 
      CALL CLEARMAT(RT ,NDOFM,3)	
      CALL CLEARMAT(RTN,NDOFM,3)	

C     Initialize Responses of bridge 
      CALL CLEARMAT(RS ,NEQ,3)	
      CALL CLEARMAT(RSN,NEQ,3)	

C	Initialize Responses of structure at wheel position
	CALL CLEARA(FWR,  3*NWHEEL)  

C	Identify element number correspoding to wheel position at zero step
      CALL CLEARI(IDWHPOS,NWHEEL)
      CALL IDWHEEL(NELW,NEID,PELW,WHPOS,IDWHPOS,XI,NWHEEL)


C     -----------------------------
      SELECTCASE(LSTATC)
      
      CASE(0) !STATIC
      
C	DECOMPOSE EFFECTIVE STIFFNESS MATRIX INTO [L]*[D]*[L]**T   !VEHICLE
      IF(NDOFM.EQ.1) THEN
        VSTIF(1) = 1.0D0/VSTIF(1)
      ELSE
        ISH=0
        CALL DECOMP(VSTIF,MAXAM,NDOFM,ISH,IOUT)
      ENDIF   
C	DIRECTLY COPY STIFFNESS TO EFFECTIVE STIFFNESS            !STRUCTURE     
      CALL MDMOVE('STIF','EFTF')
      
      
      CASE(1) !DYNAMIC
      
C	CALCULATE EFFECTIVE STIFFNESS MATRIX OF TRAIN
C	 - Before EFFSTIF, AKTR = Original stiffness matrix
C	 - After  EFFSTIF  AKTR = Decomposed Effective stiffness matrix
      CALL EFFSTIF3(VMASS,VDAMP,VSTIF,MAXAM,DELT,ALPHA,BITA,NDOFM,MAXM,MAXM)	

C	CALCULATE EFFECTIVE STIFFNESS MATRIX OF STRUCTURE
C	 - Before EFFSTIF, AK = Original stiffness matrix
C	 - After  EFFSTIF  AK = Decomposed Effective stiffness matrix
      CALL EFFSTIFB(AM,AC,AK,MAXA,DELT,ALPHA,BITA,NEQ,NWK,DDK,
	1			 'STIF','MASS','DAMP','EFTF','TEMP')
	
	
      CASE(2) !DYNAMIC WITH ONLY MOVING FORCES
      
C	DECOMPOSE EFFECTIVE STIFFNESS MATRIX INTO [L]*[D]*[L]**T   !VEHICLE
      IF(NDOFM.EQ.1) THEN
        VSTIF(1) = 1.0D0/VSTIF(1)
      ELSE
        ISH=0
        CALL DECOMP(VSTIF,MAXAM,NDOFM,ISH,IOUT)
      ENDIF      
      	

C	CALCULATE EFFECTIVE STIFFNESS MATRIX OF STRUCTURE
C	 - Before EFFSTIF, AK = Original stiffness matrix
C	 - After  EFFSTIF  AK = Decomposed Effective stiffness matrix
      CALL EFFSTIFB(AM,AC,AK,MAXA,DELT,ALPHA,BITA,NEQ,NWK,DDK,
	1			 'STIF','MASS','DAMP','EFTF','TEMP')
	
      ENDSELECT
C     -----------------------------
      

C     GATHER IRREGULARITY AT THe WHEEL POSITION
      IF(IRIN.EQ.2) THEN !IRREGULARITY GENERATION
        CALL GENIRRE3D(WHPOS,TLENGTH,1,NWHEEL,TIRAIL)			
      ELSEIF(IRIN.EQ.0.OR.IRIN.EQ.1.OR.IRIN.EQ.3) THEN
        CALL CALIRR(RRX,RRE,RRA,RRR,WHPOS,1,NWHEEL,TIRAIL) 
      ENDIF
      
C     PRINT OUTJOB & OUTTRAIN HEADER      
      CALL TITLEBTI(SOLT,TRLENG,NOUTST,NOUTVH,KOST,STROHED,KOVH,TRNOHED) 
      	          
C     ---------------------------------------------------------
C 	ITERATION PROCESS OF ZERO STEP: 
C	In order to cal. Responses of structure at nodal DOFs: RSN (later, check whether RS)                                   
C     ---------------------------------------------------------
      DO 20 JTE=1,NMAX

C     -------------
C     VEHICLE PHASE
C     -------------

C     CALCULATE FORCE APPLIED TO VEHICLE'DOFs AT ZERO STEP
      CALL CALFORV3D(NDOFM,NWHEEL,IDWHPOS,LWHDOF,RSW,FWR,FT,WCSTIF,TIRAIL) 
      
C     -----------------------------
      SELECTCASE(LSTATC)
      CASE(0,2) !STATIC & DYNAMIC WITH ONLY MOVING FORCES
        CALL STASOLV(VSTIF,RTN,FT,MAXAM,NDOFM,MAXM)
      CASE(1) !DYNAMIC
C     USING NEWMARK BETA METHOD TO CAL. DISP,VELO & ACCE AT ZERO STEP [RTN(NDOFM,3)]
C			VSTIF  = Decomposed Effective stiffness matrix
        CALL NMSOL(VMASS,VDAMP,VSTIF,RT,RTN,FT,MAXAM,DELT,ALPHA,BITA,NDOFM,MAXM,MAXM)     
      ENDSELECT
C     -----------------------------
      
C     ----------------
C     STRUCTURAL PHASE
C     ----------------
C     CALCULATE INTERACTION FORCE APPLIED TO STRUCTURE'DOFs AT ZERO STEP
      CALL CALFORS3D(FS,FC,WHSLOD,FWR,RTN,IDWHPOS,LWHDOF,NEID,PELW,XI,IA(LID),
     1               NDOFM,NWHEEL,NEQ,WCSTIF)


C     -----------------------------
      SELECTCASE(LSTATC)
      CASE(0) !STATIC
        CALL STASOLB(AK,RSN,FS,MAXA,NEQ,NWK,DDK,'EFTF')
      CASE(1,2) !DYNAMIC & DYNAMIC WITH ONLY MOVING FORCES
C     USING NEWMARK BETA METHOD TO CAL. DISP,VELO & ACCE OF ZERO STEP [RSN(NEQ,3)]
C           Decompositon & Solving is doing in this subroutine
	  CALL NMSOLB(AM,AC,AK,RS,RSN,FS,MAXA,DELT,ALPHA,BITA,NEQ,NWK,NWK,
	1			  DDK,'EFTF','MASS','DAMP')
      ENDSELECT
C     -----------------------------	
      
	CALL CALDISW(IDWHPOS,NEID,PELW,XI,IA(LID),RSN,RSW,NWHEEL,NEQ)
	
C	GET CONTACT FORCE NORM & CHECK CONVERGENCE
	DNORMN = DINORM(FC ,NWHEEL)
	DNORMO = DINORM(FCO,NWHEEL)
	IF(DNORMN.NE.0.0D0) THEN
	  DCONV = ABS((DNORMN-DNORMO)/DNORMN*100.0D0)
	ELSE
	  WRITE(*,*) 'ALL WHEEL LOSS THEIR CONTACT!'
	  STOP
	ENDIF

      FCO(1:NWHEEL) = FC(1:NWHEEL) !UPDATE OLD CONTACT FORCE
      
	IF(DCONV.LT.TOLER) EXIT  !CONVERGENCE ACHIEVED

C     ---------------------------------------------------------
   20 CONTINUE
C     ---------------------------------------------------------

C     UPDATE RESPONSE OF VEHICLE & STRUCTURE
      CALL MOVMAT(RTN,RT,NDOFM,3)
      CALL MOVMAT(RSN,RS,NEQ  ,3)       

C
C     ******************************************************************
C     *							FOLLOWING STEP                       *                                       
C     ******************************************************************
      DO 140 ITE=1,NSOLT

C	  CALCULATE (ITE+1)TH WHEEL POSITION
        TTIMC = ITE*DELT
	  CALL CLEARA(WHPOSN,NWHEEL)	 
        DO 175 J=1,NWHEEL
          WHPOSN(J)=WHPOS(J) + TTIMC*VEL
  175   CONTINUE

C	  IDENTIFY ELEMENT NUMBER CORRESPONDING TO WHEEL POSITION AT (ITE+1)TH STEP
	  CALL CLEARI(IDWHPOS,NWHEEL)
        CALL IDWHEEL(NELW,NEID,PELW,WHPOSN,IDWHPOS,XI,NWHEEL)


C     GATHER IRREGULARITY AT THe WHEEL POSITION
      IF(IRIN.EQ.2) THEN !IRREGULARITY GENERATION
        CALL GENIRRE3D(WHPOSN,TLENGTH,1,NWHEEL,TIRAIL)			
      ELSEIF(IRIN.EQ.0.OR.IRIN.EQ.1.OR.IRIN.EQ.3) THEN
        CALL CALIRR(RRX,RRE,RRA,RRR,WHPOSN,1,NWHEEL,TIRAIL) 
      ENDIF
      
C     -----------------------------------
C 	ITERATION PROCESS OF FOLLOWING STEP                                
C     -----------------------------------
        DO 150 JTE=1,NMAX

C     -------------
C     VEHICLE PHASE
C     -------------

C     CALCULATE FORCE APPLIED TO VEHICLE'DOFs AT ZERO STEP
      CALL CALFORV3D(NDOFM,NWHEEL,IDWHPOS,LWHDOF,RSW,FWR,FT,WCSTIF,TIRAIL) 

C     -----------------------------
      SELECTCASE(LSTATC)
      CASE(0,2) !STATIC & DYNAMIC WITH ONLY MOVING FORCES
        CALL STASOLV(VSTIF,RTN,FT,MAXAM,NDOFM,MAXM)
      CASE(1) !DYNAMIC
C     USING NEWMARK BETA METHOD TO CAL. DISP,VELO & ACCE AT ZERO STEP [RTN(NDOFM,3)]
C			VSTIF  = Decomposed Effective stiffness matrix
        CALL NMSOL(VMASS,VDAMP,VSTIF,RT,RTN,FT,MAXAM,DELT,ALPHA,BITA,NDOFM,MAXM,MAXM)
      ENDSELECT
C     -----------------------------

      
C     ----------------
C     STRUCTURAL PHASE
C     ----------------
C     CALCULATE INTERACTION FORCE APPLIED TO STRUCTURE'DOFs AT ZERO STEP
      CALL CALFORS3D(FS,FC,WHSLOD,FWR,RTN,IDWHPOS,LWHDOF,NEID,PELW,XI,IA(LID),
     1               NDOFM,NWHEEL,NEQ,WCSTIF)


C     -----------------------------
      SELECTCASE(LSTATC)
      CASE(0) !STATIC
        CALL STASOLB(AK,RSN,FS,MAXA,NEQ,NWK,DDK,'EFTF')
      CASE(1,2) !DYNAMIC & DYNAMIC WITH ONLY MOVING FORCES
C     USING NEWMARK BETA METHOD TO CAL. DISP,VELO & ACCE OF ZERO STEP [RSN(NEQ,3)]
C           Decompositon & Solving is doing in this subroutine
	  CALL NMSOLB(AM,AC,AK,RS,RSN,FS,MAXA,DELT,ALPHA,BITA,NEQ,NWK,NWK,
	1			  DDK,'EFTF','MASS','DAMP')
      ENDSELECT
C     -----------------------------	
	
	CALL CALDISW(IDWHPOS,NEID,PELW,XI,IA(LID),RSN,RSW,NWHEEL,NEQ) 

	
C	GET CONTACT FORCE NORM & CHECK CONVERGENCE
	DNORMN = DINORM(FC ,NWHEEL)
	DNORMO = DINORM(FCO,NWHEEL)
	IF(DNORMN.NE.0.0D0) THEN
	  DCONV = ABS((DNORMN-DNORMO)/DNORMN*100.0D0)
	ELSE
	  WRITE(*,*) 'ALL WHEEL LOSS THEIR CONTACT!'
	  STOP
	ENDIF

      FCO(1:NWHEEL) = FC(1:NWHEEL) !UPDATE OLD CONTACT FORCE
      
	IF(DCONV.LT.TOLER) EXIT  !CONVERGENCE ACHIEVED
	
                             
C     -----------------------------------
  150   CONTINUE                             
C     -----------------------------------> end of iteration for time step 

C     STORE MAXIMUM WHEEL CONTACT FORCE
      DO IWHEEL = 1,NWHEEL
          IF(ITE.EQ.1) THEN
              FCMAX(1,IWHEEL) = FCO(IWHEEL)
              FCMAX(2,IWHEEL) = FCO(IWHEEL)
          ELSE
              IF(FCMAX(1,IWHEEL).GT.FCO(IWHEEL)) FCMAX(1,IWHEEL) = FCO(IWHEEL)
              IF(FCMAX(2,IWHEEL).LT.FCO(IWHEEL)) FCMAX(2,IWHEEL) = FCO(IWHEEL)
          ENDIF
      ENDDO
      
C     UPDATE RESPONSE OF VEHICLE & STRUCTURE
      CALL MOVMAT(RTN,RT,NDOFM,3)
      CALL MOVMAT(RSN,RS,NEQ  ,3)   
     
      CALL PRNOUTBTI(IA(LID),NEQ,NDOFM,RS,RT,NOUTST,NOUTVH,KOST,KOVH,ITE,JTE,
     1               NSOLT,DCONV,TTIMC,WHPOSN(1),VOUTMAXS,VOUTMAXV) 
      

C      IF(IVTYPE.EQ.1) THEN
C        CALL TRAINRENRES(NTRAIN,NCARB,NBOGI,NWHEEL,RT,ITE,VEL,TTIMC,TRNDIM)
C      ENDIF      

      GOTO 140
C	--------------------------------------------------------------
C	ADD THIS BLOCK FOR PRINTING STRESS OF ALL STEP SONGSAK JUL2007
C	NEW OUTPUT SONGSAK JUL2007
	CALL DISOUT(IA(LID),RS(1,1))
	A(LDT:LDT+NEQ-1) = RS(1:NEQ,1)

      IFPRI = 0
      IFPLO = 0
	IFEIG = 1
      ITASK = 3
      IFREF = 1
      ISPRI = 0
      CALL GRLOOP (IA(LEL),KSC)

C	NEW OUTPUT SONGSAK JUL2007
	CALL PRNFLAG('NONE','NONE','GSUP','LSUP','DISP','NONE','NONE')
	CALL  PRNOUT('STND','PONE','NONE',ITE)
C	--------------------------------------------------------------

C	--------------------------------------------------------------
  140 CONTINUE !Continue next time step <<------------------>>
C	--------------------------------------------------------------

      CALL PRNMAXBTI(NOUTST,NOUTVH,VOUTMAXS,VOUTMAXV,WHPOS,WHSLOD,FCMAX,NWHEEL) 
        
  
      DEALLOCATE(DMASS,DDAMP,DSTIF)

      
      RETURN  
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE REDMASS(NDOFM,NWHEEL,LSTATC,IVTYPE,NTRAIN,NBOGI,IWRIN,NOUTST,NOUTVH,SPACING,NCOPY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     LSTATC = FLAG FOR STATIC ANALYSIS   0=STATIC 1=DYNAMIC 2=DYNAMIC WITH ONLY MOVING FORCE

C     *** FOR TRAIN ***
C     NPCARF = NUMBER OF INDEPENDENT CAR OR FRONT POWER CAR !FOR INDEPENDENT TRAIN  THERE WITH BE ONLY -- NPCARF
C     NACARB  = NUMBER OF ARTIGULATED CAR
C     NPCARR = NUMBER OF INDEPENDENT CAR OR REAR  POWER CAR

C	NBOGI = Total number of bogies     = NPCARF*2+(NACARB+1)+NPCARR*2
C	NWHEEL= Total number of wheels     = NBOGI*2

C	IWRIN : Idenfifier for Input of Wheel-Rail contact Mechanism 
C	   = 0: Input the contact stiffness & damping directly (OLD)
C	   = 1: Input the wheel-rail properties and calculate contact stiffness & damping
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
C	----------------------------------------------------------------
      DIMENSION NTRAIN(3)
      
C     ------------------
C     READ CONTROL DATA:
C     ------------------
      READ(ITI,*)
      READ(ITI,*) ALPHA,BITA,DELT,NMAX,TOLER,SCALE,LSTATC,IVTYPE,NOUTST
      
C     ------------------------------      
      SELECTCASE(IVTYPE)
C     ------------------------------   
      
      CASE(0) !DIRECT INPUT OF VEHICLE MATRICES
          READ(ITI,*)
          READ(ITI,*) NDOFM,NWHEEL,VEL,BIDIS,EXTDIS,NOUTVH,SPACING,NCOPY	
          NCOPY = NCOPY+1
      	  NDOFM =  NDOFM*NCOPY
      	  NWHEEL= NWHEEL*NCOPY
      	
      CASE(1) !TRAIN 
          READ(ITI,*)
          READ(ITI,*) NPCARF,NACARB,NPCARR,IWRIN,VEL,BIDIS,EXTDIS,NOUTVH	 	 
    	
C         IN GENERAL ARTIGULATED TRAIN HAVE AT LEAST 3 PASSENGER CAR (2-TRANSITION + PASSENGER)
          IF(NACARB.EQ.0) THEN
            INDAR = 0     !0=INDEPENDENT TRAIN  1=ARTIGULATED TRAIN
          ELSEIF(NACARB.GE.3) THEN
            INDAR = 1     !0=INDEPENDENT TRAIN  1=ARTIGULATED TRAIN
          ELSE
            WRITE(*,*) 'ERROR IN ARTIGULATED TRAIN DATA'  
            WRITE(*,*) 'CHECK THE NUMBER OF TRANSITIONAL AND PASSENGER CAR' 
            STOP
          ENDIF

	    NTRAIN(1:3) = [NPCARF,NACARB,NPCARR]
          SELECTCASE(INDAR)
          CASE(0) !INDEPENDENT TRAIN
	      NBOGI = NPCARF*2 + NPCARR*2
	      NDOFM = (5+2*(5+3+3))*NPCARF + (5+2*(5+3+3))*NPCARR	
C	      NWKTR = NPCARF*167 + NPCARR*167
          CASE(1) !ARTIGULATED TRAIN
	      NBOGI = NPCARF*2 + (NACARB+1 ) + NPCARR*2
	      NDOFM = (5+2*(5+3+3))*NPCARF + NACARB*(2*3+5+5)+11 + NPCARR*(5+2*(5+3+3))
C	      NWKTR = NPCARF*167 + (83+(NACARB-1)*(78+29)+72) + NPCARR*167	  
          ENDSELECT
	    NWHEEL= NBOGI*2
	
C     ------------------------------   
	ENDSELECT
C     ------------------------------   

	
C     INITIALIZE DATA
		IF(ALPHA.EQ.0.0) ALPHA=0.25
		IF(BITA.EQ.0.0) BITA=0.5
		IF(ALPHA.GT.0.25.AND.BITA.NE.0.5) THEN
			WRITE(*,100) ALPHA,BITA
			STOP
		ENDIF
		IF(NMAX.EQ.0) NMAX=100
		IF(TOLER.EQ.0.0) TOLER=1.0D-06
		
		VEL = VEL*5.0/18.0 !FROM kMpH to m/s
		

  100 FORMAT(//,1X,'INPUT ERROR IN NEWMARK COEFFICIENT',/
     &1X,'RECOMMENDED VALUES ALPHA = 0.15, BITA=0.5 BUT INPUT VALUES',/
     &1X,'ARE AS FOLLOWING',/
     &1X,'  ALPHA = ',F10.4/
     &1X,'  BITA  = ',F10.4)

      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE REDVMAT(NDOFM,NWHEEL,MAXAM,VMASS,VDAMP,VSTIF,WHPOS,LWHDOF,WHSLOD,WCSTIF,SPACING,NCOPY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	----------------------------------------------------------------
      DIMENSION WHPOS(NWHEEL),LWHDOF(3,NWHEEL),WHSLOD(NWHEEL),WCSTIF(3,2,NWHEEL) !WCSTIF  3=component  2=Stif&Damping
      DIMENSION VMASS(1),VDAMP(1),VSTIF(1)
      DIMENSION VMAS(NDOFM,NDOFM),VDAM(NDOFM,NDOFM),VSTF(NDOFM,NDOFM)
      DIMENSION LCOLH(NDOFM),MAXAM(NDOFM+1)
      
      MDOFM =  NDOFM/NCOPY
      MWHEEL= NWHEEL/NCOPY
      
      VMAS(1:NDOFM,1:NDOFM) = 0.0D0 
      VDAM(1:NDOFM,1:NDOFM) = 0.0D0 
      VSTF(1:NDOFM,1:NDOFM) = 0.0D0 
      
C     READ MASS MATRIX      
      READ(ITI,*)
      DO IDOFM = 1,MDOFM
        READ(ITI,*) VMAS(IDOFM,1:MDOFM)
      ENDDO
      
C     READ DAMPING MATRIX      
      READ(ITI,*)
      DO IDOFM = 1,MDOFM
        READ(ITI,*) VDAM(IDOFM,1:MDOFM)
      ENDDO

C     READ STIFFNESS MATRIX      
      READ(ITI,*)
      DO IDOFM = 1,MDOFM
        READ(ITI,*) VSTF(IDOFM,1:MDOFM)
      ENDDO
      
      
      IF(NCOPY.GT.1) THEN
          DO ICOPY = 2,NCOPY
              IDF = MDOFM*(ICOPY-1)
              JDF = MDOFM*(ICOPY-1)
              DO I = 1,MDOFM
                  DO J = 1,MDOFM
                        VMAS(I+IDF,J+JDF) = VMAS(I,J)
                        VDAM(I+IDF,J+JDF) = VDAM(I,J)
                        VSTF(I+IDF,J+JDF) = VSTF(I,J)
                  ENDDO
              ENDDO
          ENDDO
      ENDIF

C     COLUMN HEIGHT
      MAXAM(1) = 1
      DO IDOFM = 1,NDOFM
          ICH = IDOFM
          DO JDOFM = 1,IDOFM
              VALM = VMAS(IDOFM,JDOFM)
              VALC = VDAM(IDOFM,JDOFM)
              VALK = VSTF(IDOFM,JDOFM)
              IF(VALM.NE.0.0D0.OR.VALC.NE.0.0D0.OR.VALK.NE.0.0D0) THEN
                  ICH = IDOFM - JDOFM + 1
                  EXIT
              ENDIF
          ENDDO
          LCOLH(IDOFM) = ICH
          MAXAM(IDOFM+1) = MAXAM(IDOFM) + ICH
      ENDDO

C     RESTORE MATRIX INTO SKYLINE FORMAT
      KDOFM = 0
      DO IDOFM = 1,NDOFM
          ICH = LCOLH(IDOFM)
          JJ = IDOFM
          DO JDOFM = IDOFM-ICH+1,IDOFM
            KDOFM = KDOFM + 1
            VMASS(KDOFM) = VMAS(IDOFM,JJ)
            VDAMP(KDOFM) = VDAM(IDOFM,JJ)
            VSTIF(KDOFM) = VSTF(IDOFM,JJ)
            JJ = JJ - 1
          ENDDO
      ENDDO      



C     READ WHEEL DOF AND STATIC WHEEL LOAD
      READ(ITI,*)
      DO IWHEEL = 1,MWHEEL
        READ(ITI,*) IW,LWHDOF(1:3,IWHEEL),WHSLOD(IWHEEL)  !LWHDOF(1,IWHEEL) = DOF NUMBER (POSITION 1=Lateral 2=Rolling 3=Vertical)
      ENDDO    

C     READ WHEEL RELATIVE POSITION AND THEN GENERATE ABSOLUTE WHELL POSITION
      READ(ITI,*)
      TLEN = 0.0
      DO IWHEEL = 1,MWHEEL
          READ(ITI,*) IW,SPEC
          WHPOS(IWHEEL) = TLEN
          TLEN = TLEN - SPEC
      ENDDO  

C     READ WHEEL CONTACT STIFFNESS AND DAMPING
      READ(ITI,*)
      DO IWHEEL = 1,MWHEEL
        READ(ITI,*) IW,CSTIFL,CDAMPL,CSTIFR,CDAMPR,CSTIFV,CDAMPV
        WCSTIF(1,1,IWHEEL) = CSTIFL !LATERAL
        WCSTIF(1,2,IWHEEL) = CDAMPL !LATERAL
        WCSTIF(2,1,IWHEEL) = CSTIFR !ROLLING
        WCSTIF(2,2,IWHEEL) = CDAMPR !ROLLING
        WCSTIF(3,1,IWHEEL) = CSTIFV !VERTICAL
        WCSTIF(3,2,IWHEEL) = CDAMPV !VERTICAL
C        IDFW = LWHDOF(1,IWHEEL) !LATERAL
C        IF(IDFW.NE.0) THEN
C            INDX = MAXAM(IDFW)
C            VSTIF(INDX) = VSTIF(INDX) + CSTIFL !ADD TO SYSTEM MATRIX
C            VDAMP(INDX) = VDAMP(INDX) + CDAMPL !ADD TO SYSTEM MATRIX
C        ENDIF
C        IDFW = LWHDOF(2,IWHEEL) !ROLLING
C        IF(IDFW.NE.0) THEN
C            INDX = MAXAM(IDFW)
C            VSTIF(INDX) = VSTIF(INDX) + CSTIFR !ADD TO SYSTEM MATRIX
C            VDAMP(INDX) = VDAMP(INDX) + CDAMPR !ADD TO SYSTEM MATRIX
C        ENDIF
C        IDFW = LWHDOF(3,IWHEEL) !VERTICAL
C        IF(IDFW.NE.0) THEN
C            INDX = MAXAM(IDFW)
C            VSTIF(INDX) = VSTIF(INDX) + CSTIFV !ADD TO SYSTEM MATRIX
C            VDAMP(INDX) = VDAMP(INDX) + CDAMPV !ADD TO SYSTEM MATRIX
C        ENDIF
      ENDDO 



      IF(NCOPY.GT.1) THEN
          TL = TLEN
          DO ICOPY = 2,NCOPY
              IWL = MWHEEL*(ICOPY-1)
              
              DO I = 1,MWHEEL
                DO J = 1,3
                    LWHDOF(J,I+IWL) = LWHDOF(J,I)
                    IF(LWHDOF(J,I).NE.0) LWHDOF(J,I+IWL) = LWHDOF(J,I) + MDOFM*(ICOPY-1)
                ENDDO
                WHSLOD(I+IWL) = WHSLOD(I)
              ENDDO
              
              TL = TL - SPACING
              DO I = 1,MWHEEL
                WHPOS(I+IWL) = TL +  WHPOS(I)
              ENDDO
              TL = TL + TLEN
 
              DO I = 1,MWHEEL
                WCSTIF(1:3,1,I+IWL) = WCSTIF(1:3,1,I) !STIF
                WCSTIF(1:3,2,I+IWL) = WCSTIF(1:3,2,I) !DAMP
              ENDDO

          ENDDO
      ENDIF
      
      
                  

      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE CALFORV3D(NDOFM,NWHEEL,IDWHPOS,LWHDOF,RSW,FWR,FT,WCSTIF,TIRAIL) 
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C     TIRAIL(3,1) = TRACK IRREGULARITY  1=LATERAL 2=ROLLING 3=VERTICAL
C	-----------------------------------------------------------------------

      DIMENSION IDWHPOS(1),LWHDOF(3,1),RSW(3,2,1),WCSTIF(3,2,1),TIRAIL(3,1) !WCSTIF  3=component  2=Stif&Damping
      DIMENSION FWR(3,1),FT(1)
      
      
C	FORM FORCE ACTING ON VEHICLE WHEELS:
C	------------------------------------
	CALL CLEARMAT(FWR,3,NWHEEL)

	DO 100 IWH = 1,NWHEEL
	  DO ICOM = 1,3
	    IF(IDWHPOS(IWH).GT.0) THEN
	      FWR(ICOM,IWH) = WCSTIF(ICOM,1,IWH)*RSW(ICOM,1,IWH) + WCSTIF(ICOM,2,IWH)*RSW(ICOM,2,IWH) !VERTICAL FORCE   WHEEL ON BRIDGE
	      FWR(ICOM,IWH) = FWR(ICOM,IWH) + WCSTIF(ICOM,1,IWH)*TIRAIL(ICOM,IWH)                     !EFFECT OF TRACK IRREGULARITY
	    ELSE
	      FWR(ICOM,IWH) = WCSTIF(ICOM,1,IWH)*0.0             + WCSTIF(ICOM,2,IWH)*0.0             !VERTICAL FORCE   WHEEL NOT ON BRIDGE
	      FWR(ICOM,IWH) = FWR(ICOM,IWH) + WCSTIF(ICOM,1,IWH)*TIRAIL(ICOM,IWH)                     !EFFECT OF TRACK IRREGULARITY
	    ENDIF
	  ENDDO  
100   CONTINUE

C     STORE IN OVERALL VEHICLE LOAD VECTOR
	CALL CLEARA(FT ,NDOFM)
	
	DO 200 IWH = 1,NWHEEL
	  DO ICOM = 1,3
	      IDFW = LWHDOF(ICOM,IWH)
            IF(IDFW.NE.0) THEN
	        FT(IDFW) = FWR(ICOM,IWH)
            ENDIF
	  ENDDO  
200   CONTINUE

      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE CALFORS3D(FS,FC,WHSLOD,FWR,RT,IDWHPOS,LWHDOF,NEID,PELW,XI,ID,
     1                     NDOFM,NWHEEL,NEQ,WCSTIF)	
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
     
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
	
C	GRAVTITY DIRECTION 
	COMMON /MGRAV/ NGRAV	
C	-----------------------------------------------------------------------
	
      DIMENSION FS(1),FC(1),FWR(3,1),RT(NDOFM,3)
      DIMENSION IDWHPOS(1),LWHDOF(3,1),WHSLOD(1),WCSTIF(3,2,1) !WCSTIF  3=component  2=Stif&Damping
      DIMENSION NEID(2,1),PELW(4,1),XI(1),ID(NSF,1)
      
      DIMENSION FSS(3,NWHEEL)
	DIMENSION VR(3),VS(3),VT(3),TRANS(14,14),FIXD(14),FIXG(14),FSE(7,2)
      
C	FORM CONTACT FORCE BETWEEN WHEELS & BRIDGE:
C	-----------------------------------------
	CALL CLEARMAT(FSS,3,NWHEEL)

	DO 10 IWH = 1,NWHEEL
	  DO ICOM = 1,3
	    IDFW = LWHDOF(ICOM,IWH)
          IF(IDFW.NE.0) THEN
	        FTOP = WCSTIF(ICOM,1,IWH)*RT(IDFW,1) + WCSTIF(ICOM,2,IWH)*RT(IDFW,2)

              FSS(ICOM,IWH) = FTOP - FWR(ICOM,IWH)
                
	        IF(ICOM.EQ.3) THEN
	          FSS(ICOM,IWH) = FSS(ICOM,IWH) - WHSLOD(IWH) 
	          IF(FSS(ICOM,IWH).GT.0.0D0) FSS(ICOM,IWH) = 0.0D0 !NO TENSION VERTICAL FORCE ON RAIL
	        ENDIF
	    ENDIF
	  ENDDO  
	  FC(IWH) = FSS(3,IWH) !STORE VERTICAL CONTACT FORCE FOR CONVERGENCE CHECKING
10    CONTINUE	

C	-----------------------------------------------------------------------------------
C	EQUIVALENT NODAL FORCE CORRESPONDING TO INTERACTION LUMPED FORCE AT WHEEL POSITION
C	-----------------------------------------------------------------------------------
	CALL CLEARA(FS,NEQ)

      DO 40 I=1,NWHEEL
        IF(IDWHPOS(I).GT.0) THEN

C	CALCULATE EQUIVALENT NODAL FORCES VECTOR OF THE CURRENT BEAM ELEMENT
            IENUM=IDWHPOS(I)		!ID of current element
            
C     LOCAL VECTOR & SEGMENT LENGTH		
            VR(1:3) = PELW(1:3,IENUM)
               SLEN = PELW(  4,IENUM)
               
            CALL FMVEVR(VR,VS,VT)
            RANG = 0.0D0
	      CALL ROMBAC(VR,VS,VT,RANG)  
	      CALL TRANLG(VR,VS,VT,TRANS)
    	
	      WR = VR(NGRAV)*FSS(3,I)
	      WS = VS(NGRAV)*FSS(3,I)
	      WT = VT(NGRAV)*FSS(3,I) + FSS(1,I)
	      MR = -WS*ECC + FSS(2,I)
	          	        
	      CI = 1.0-XI(I) !1.0-AL/ELN
            CJ = XI(I)     !AL/ELN
    
C	LOCAL AXIAL FORCE
	      FIXD(1)  = CI*WR
	      FIXD(8)  = CJ*WR

C	SHEAR IN S-AXIS AND MOMENT IN T-AXIS
	      FIXD(2)  = CI*WS
	      FIXD(9)  = CJ*WS
	
C	SHEAR IN T-AXIS AND MOMENT IN S-AXIS
	      FIXD(3)  = CI*WT
	      FIXD(10) = CJ*WT
	  
C	LOCAL TORQUE
	      FIXD(4)  = CI*MR
	      FIXD(11) = CJ*MR

C     TRANSFORM LOCAL FORCES TO GLOBAL
	      FIXG = 0.0
	      FIXG = MATMUL(TRANS,FIXD)
       
            FSE(1:7,1) = FIXG(1+0:7+0)
            FSE(1:7,2) = FIXG(1+7:7+7)
        
C	ASSEMBLE EQUIVALENT NODAL FORCE OF ELEMENT INTO GLOBAL NODAL FORCE
            DO 50 J=1,2			        !Loops over number of nodes per element
                NDNUM=NEID(J,IENUM)	        !Current node number
                DO 60 ISF=1,NSF		        !Loops over number of DOFs per node
			        IDF = IDOF(ISF)         !DOF NUMBER
			        IF(IDF.GT.0.AND.IDF.LE.7) THEN
				        IEQ=ID(ISF,NDNUM)	!Equation number of current DOF
				        IF(IEQ.GT.0) THEN
				            FS(IEQ)=FS(IEQ)+FSE(IDF,J)
				        ENDIF
			        ENDIF
   60           CONTINUE
   50       CONTINUE
   
        ENDIF			!----->(IDWHPOS(I).GT.0)
   40 CONTINUE		! Loops over NWHEEL

      RETURN
      END

C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE GETSTDIS(ID,NCNOD,NCDOF,NEQ,DISP,VAL,OPTN)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPTN
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
     
      DIMENSION ID(NSF,1),DISP(NEQ,1)
      
      IEQ = ID(NCDOF,NCNOD)
      
      SELECTCASE(OPTN)
          CASE('DIS')
          NC = 1
          CASE('VEL')
          NC = 2
          CASE('ACC')
          NC = 3
      ENDSELECT
      
      VAL = 0.0D0
      IF(IEQ.NE.0) VAL = DISP(IEQ,NC)
      
      IF(ABS(VAL).LT.1.0E-30) VAL = 0.0D0
      
      RETURN
      END

      	
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE CALDISW(IDWHPOS,NEID,PELW,XI,ID,RS,RSW,NWHEEL,NEQ)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	-------------------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
     
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
C
      DIMENSION IDWHPOS(1),NEID(2,1),PELW(4,1),XI(1),ID(NSF,1)
      DIMENSION VR(3),VS(3),VT(3),TRANS(14,14),RL(14,3)
	DIMENSION RSW(3,2,1),RS(NEQ,3)
      DIMENSION RSE(14,3)
C	------------------------------------------------------------------------

      DO 10 I=1,NWHEEL
      
        IF(IDWHPOS(I).GT.0) THEN

C		CALCULATE RESPONSES OF ELEMENT WHERE WHEEL RUNS
        IENUM=IDWHPOS(I)

		CALL CLEARMAT (RSE,14,3) ! 3 for [DIS,VEL,ACCE]

          DO 30 J=1,2
            NDNUM=NEID(J,IENUM)
            DO 40 ISF=1,NSF
			    IDF = IDOF(ISF)
			    IF(IDF.GT.0.AND.IDF.LE.7) THEN !SONGSAK CHANGED NSF TO 7 OCT2007
				    IEQ=ID(ISF,NDNUM)
				    IF(IEQ.GT.0) THEN
				        DO 50 L=1,3
					        RSE((J-1)*7+IDF,L)=RS(IEQ,L)
   50		                CONTINUE
				    ENDIF
			    ENDIF
   40       CONTINUE
   30     CONTINUE


C     TRANSFORM THE NODAL DISP. INTO RAIL COORDINATE
          VR(1:3) = PELW(1:3,IENUM)
	    CALL FMVEVR(VR,VS,VT)
          RANG = 0.0D0
	    CALL ROMBAC(VR,VS,VT,RANG)
	    CALL TRANLG(VR,VS,VT,TRANS)
	    RL = MATMUL(TRANSPOSE(TRANS),RSE)
	
	
C	    INTERPOLATE ELEMENT RESPONSES INTO RESPONSES AT WHEEL POSITIONS (AT SECTION CENTER)
          DO 60 J=1,2 ![DIS,VEL]
         
            CI = 1.0-XI(I) ; CJ = XI(I)
            
			RSW(1,J,I) = CI*RL(3  ,J) + CJ*RL(7+3,J)      !Zwi lateral
			
			RSW(2,J,I) = CI*RL(4  ,J) + CJ*RL(7+4,J)      !Rxwi(Twist) rolling
			
			RSW(3,J,I) = CI*RL(2  ,J) + CJ*RL(7+2,J)      !Ywi vertical
            RSW(3,J,I) = RSW(3,J,I) - RSW(2,J,I)*ECC
 
   60     CONTINUE

        ELSE
        
          DO 70 J=1,2 ![DIS,VEL]
            RSW(1:3,J,I) = 0.0D0
   70     CONTINUE
   
        ENDIF
        
   10 CONTINUE

  100 RETURN
      END
      
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE STASOLV(AK,RN,FN,MAXA,NEQ,NWK)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C     STATIC SOLUTION OF VEHICLE RESPONSE
C	-----------------------------------------------------------------------
      DIMENSION AK(NWK),RN(NEQ,3),FN(NEQ)
      DIMENSION MAXA(NEQ+1)
C	-----------------------------------------------------------------------

      RN(1:NEQ,1) = FN(1:NEQ)
      
      IF(NEQ.EQ.1) THEN
        RN(1,1) = AK(1)*RN(1,1)
      ELSE
        CALL REDBAK(AK,RN(1,1),MAXA,NEQ) !SOLVE FOR STATIC DISP
      ENDIF      
            
      
      RN(1:NEQ,2) = 0.0D0 !NO VELOCITY FOR STATIC
      RN(1:NEQ,3) = 0.0D0 !NO ACCELERATION FOR STATIC

      RETURN
      END

C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE STASOLB(AK,RN,FN,MAXA,NEQ,NWK,DDK,EFTF)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 EFTF
C	----------------------------------------------------------------------
C     STATIC SOLUTION OF STRUCTURAL RESPONSE
C	-----------------------------------------------------------------------
      DIMENSION AK(NWK),RN(NEQ,3),FN(NEQ)
      DIMENSION MAXA(NEQ+1),DDK(1)
C	-----------------------------------------------------------------------

      
      RN(1:NEQ,1) = FN(1:NEQ)
	CALL COLSOL(MAXA,AK,DDK,RN(1,1),1,INDPD,EFTF,'TEMP') !STIFFNESS DECOMPOSITION
	CALL COLSOL(MAXA,AK,DDK,RN(1,1),2,INDPD,EFTF,'TEMP') !SOLVE FOR STATIC DISP

      RN(1:NEQ,2) = 0.0D0 !NO VELOCITY FOR STATIC
      RN(1:NEQ,3) = 0.0D0 !NO ACCELERATION FOR STATIC
      

      RETURN
      END

C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE MAXATRAIN(NTRAIN,NDOFM,MAXAM,LWHDOF)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------------------------------------------
C	MAKING VECTOR CONTANING NUMBER OF ELEMENTS IN COLUMNS OF TRAIN'S STIFF/DAMPING MATRIX
C	--------------------------------------------------------------------------------------
	DIMENSION NTRAIN(3),MAXAM(1),LWHDOF(3,1)
	DIMENSION NCHTR(NDOFM)
	

      INDAR = 0 !0=INDEPENDENT TRAIN  1=ARTIGULATED TRAIN
      IF(NTRAIN(2).NE.0) INDAR = 1
            
C     ------------------      
      SELECTCASE(INDAR)
C     ------------------ 
      CASE(0) !INDEPENDENT TRAIN 
        CALL TRAINCOLHI(NTRAIN,NCHTR,NDOFM,LWHDOF)
      CASE(1) !ARTIGULATED TRAIN
        CALL TRAINCOLHA(NTRAIN,NCHTR,NDOFM,LWHDOF)
C     ------------------ 
      ENDSELECT
C     ------------------ 
      
C	-------------------------------------------------
C	MAKING VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS
C	-------------------------------------------------
      MAXAM(1) = 1
      DO IDOFM = 1,NDOFM
          ICH = NCHTR(IDOFM)
          MAXAM(IDOFM+1) = MAXAM(IDOFM) + ICH
      ENDDO
      

      RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE TRAINCOLHI(NTRAIN,NCHTR,NDOFM,LWHDOF)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	COLUMN HEIGHT FOR INDEPENDENT TRAIN
C	-----------------------------------------------------------------------
      DIMENSION NTRAIN(3),NCHTR(1),LWHDOF(3,1)
      
	CALL CLEARI(NCHTR,NDOFM)
	
	ICOL = 0
      IWHEEL = 0
      
C	Front power & independent cars:
C	-----------------
	DO 60 IP=1,NTRAIN(1)  
		NCHTR(ICOL+1 ) = 1 !!BODY
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 1
		
		NCHTR(ICOL+6 ) = 6 !FRONT BOGIE
		NCHTR(ICOL+7 ) = 7
		NCHTR(ICOL+8 ) = 1
		NCHTR(ICOL+9 ) = 6
		NCHTR(ICOL+10) = 1
		
		NCHTR(ICOL+11) = 6 !FRONT BOGIE WHEEL		
		NCHTR(ICOL+12) = 6
		NCHTR(ICOL+13) = 5
		NCHTR(ICOL+14) = 9	
		NCHTR(ICOL+15) = 9
		NCHTR(ICOL+16) = 8
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+11 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+12 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+13 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+14 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+15 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+16 !VERTICAL
		
		NCHTR(ICOL+17) = 17 !REAR BOGIE
		NCHTR(ICOL+18) = 18
		NCHTR(ICOL+19) = 1
		NCHTR(ICOL+20) = 17
		NCHTR(ICOL+21) = 1
		
		NCHTR(ICOL+22) = 6 !REAR BOGIE WHEEL	
		NCHTR(ICOL+23) = 6
		NCHTR(ICOL+24) = 5
		NCHTR(ICOL+25) = 9
		NCHTR(ICOL+26) = 9
		NCHTR(ICOL+27) = 8

		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+22 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+23 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+24 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+25 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+26 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+27 !VERTICAL
		
		ICOL = ICOL+27
   60 CONTINUE


C	Rear power & independent cars:
C	----------------
	DO 80 IP=1,NTRAIN(3)   
		NCHTR(ICOL+1 ) = 1 !!BODY
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 1
		
		NCHTR(ICOL+6 ) = 6 !FRONT BOGIE
		NCHTR(ICOL+7 ) = 7
		NCHTR(ICOL+8 ) = 1
		NCHTR(ICOL+9 ) = 6
		NCHTR(ICOL+10) = 1
		
		NCHTR(ICOL+11) = 6 !FRONT BOGIE WHEEL		
		NCHTR(ICOL+12) = 6
		NCHTR(ICOL+13) = 5
		NCHTR(ICOL+14) = 9	
		NCHTR(ICOL+15) = 9
		NCHTR(ICOL+16) = 8
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+11 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+12 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+13 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+14 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+15 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+16 !VERTICAL
		
		NCHTR(ICOL+17) = 17 !REAR BOGIE
		NCHTR(ICOL+18) = 18
		NCHTR(ICOL+19) = 1
		NCHTR(ICOL+20) = 17
		NCHTR(ICOL+21) = 1
		
		NCHTR(ICOL+22) = 6 !REAR BOGIE WHEEL	
		NCHTR(ICOL+23) = 6
		NCHTR(ICOL+24) = 5
		NCHTR(ICOL+25) = 9
		NCHTR(ICOL+26) = 9
		NCHTR(ICOL+27) = 8

		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+22 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+23 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+24 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+25 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+26 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+27 !VERTICAL
		
		ICOL = ICOL+27
   80 CONTINUE


      RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE TRAINCOLHA(NTRAIN,NCHTR,NDOFM,LWHDOF)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	COLUMN HEIGHT FOR ARTIGULATED TRAIN
C	-----------------------------------------------------------------------
      DIMENSION NTRAIN(3),NCHTR(1),LWHDOF(3,1)
      
	CALL CLEARI(NCHTR,NDOFM)
	
	ICOL = 0
      IWHEEL = 0
	
C	Front power & independent cars:
C	-----------------
	DO 60 IP=1,NTRAIN(1)   
		NCHTR(ICOL+1 ) = 1 !!BODY
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 1
		
		NCHTR(ICOL+6 ) = 6 !FRONT BOGIE
		NCHTR(ICOL+7 ) = 7
		NCHTR(ICOL+8 ) = 1
		NCHTR(ICOL+9 ) = 6
		NCHTR(ICOL+10) = 1
		
		NCHTR(ICOL+11) = 6 !FRONT BOGIE WHEEL		
		NCHTR(ICOL+12) = 6
		NCHTR(ICOL+13) = 5
		NCHTR(ICOL+14) = 9	
		NCHTR(ICOL+15) = 9
		NCHTR(ICOL+16) = 8

		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+11 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+12 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+13 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+14 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+15 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+16 !VERTICAL
				
		NCHTR(ICOL+17) = 17 !REAR BOGIE
		NCHTR(ICOL+18) = 18
		NCHTR(ICOL+19) = 1
		NCHTR(ICOL+20) = 17
		NCHTR(ICOL+21) = 1
		
		NCHTR(ICOL+22) = 6 !REAR BOGIE WHEEL	
		NCHTR(ICOL+23) = 6
		NCHTR(ICOL+24) = 5
		NCHTR(ICOL+25) = 9
		NCHTR(ICOL+26) = 9
		NCHTR(ICOL+27) = 8

		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+22 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+23 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+24 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+25 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+26 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+27 !VERTICAL
				
		ICOL = ICOL+27
   60 CONTINUE



C	Articulated group (Transition and Passenger Cars)
C	-----------------
C   	 1.The first transition car:
		NCHTR(ICOL+1 ) = 1 !WHEEL
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 2	
		NCHTR(ICOL+6 ) = 1

		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+1 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+2 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+3 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+4 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+5 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+6 !VERTICAL
				
		NCHTR(ICOL+7 ) = 7 !BOGIE
		NCHTR(ICOL+8 ) = 8
		NCHTR(ICOL+9 ) = 9
		NCHTR(ICOL+10) = 8
		NCHTR(ICOL+11) = 9
		
		NCHTR(ICOL+12) = 6 !BODY				
		NCHTR(ICOL+13) = 7
		NCHTR(ICOL+14) = 8 
		NCHTR(ICOL+15) = 6
		NCHTR(ICOL+16) = 7
						 
		ICOL = ICOL+16	

C	 2.The passenger cars (1 <n <N),& last transition car (n = N)
C	   looping over each pair (bogie & car body):
	        DO 70 IA = 2,NTRAIN(2)
		        NCHTR(ICOL+1 ) = 1 !WHEEL
		        NCHTR(ICOL+2 ) = 2
		        NCHTR(ICOL+3 ) = 1
		        NCHTR(ICOL+4 ) = 1
		        NCHTR(ICOL+5 ) = 2
		        NCHTR(ICOL+6 ) = 1

		        IWHEEL = IWHEEL + 1
		        LWHDOF(1,IWHEEL) = ICOL+1 !LATERAL
		        LWHDOF(2,IWHEEL) = ICOL+2 !ROLLING
		        LWHDOF(3,IWHEEL) = ICOL+3 !VERTICAL
        		
		        IWHEEL = IWHEEL + 1
		        LWHDOF(1,IWHEEL) = ICOL+4 !LATERAL
		        LWHDOF(2,IWHEEL) = ICOL+5 !ROLLING
		        LWHDOF(3,IWHEEL) = ICOL+6 !VERTICAL
		        		
		        NCHTR(ICOL+7 ) = 7 + 0 + 5 !BOGIE
		        NCHTR(ICOL+8 ) = 8 + 0 + 5
		        NCHTR(ICOL+9 ) = 9 + 0 + 5
		        NCHTR(ICOL+10) = 8 + 2 + 5
		        NCHTR(ICOL+11) = 9 + 2 + 5	

		        NCHTR(ICOL+12) = 6 !BODY	 			
		        NCHTR(ICOL+13) = 7
		        NCHTR(ICOL+14) = 3
		        NCHTR(ICOL+15) = 6
		        NCHTR(ICOL+16) = 7

		        ICOL = ICOL+16
70            CONTINUE

C      3.The last bogie: (n = N+1) 
		NCHTR(ICOL+1 ) = 6 !BOGIE
		NCHTR(ICOL+2 ) = 7
		NCHTR(ICOL+3 ) = 3 
		NCHTR(ICOL+4 ) = 6
		NCHTR(ICOL+5 ) = 7 

		NCHTR(ICOL+6 ) = 6 !WHEEL
		NCHTR(ICOL+7 ) = 6
		NCHTR(ICOL+8 ) = 5
		NCHTR(ICOL+9 ) = 9
		NCHTR(ICOL+10) = 9
		NCHTR(ICOL+11) = 8

        IWHEEL = IWHEEL + 1
        LWHDOF(1,IWHEEL) = ICOL+6 !LATERAL
        LWHDOF(2,IWHEEL) = ICOL+7 !ROLLING
        LWHDOF(3,IWHEEL) = ICOL+8 !VERTICAL
		
        IWHEEL = IWHEEL + 1
        LWHDOF(1,IWHEEL) = ICOL+9  !LATERAL
        LWHDOF(2,IWHEEL) = ICOL+10 !ROLLING
        LWHDOF(3,IWHEEL) = ICOL+11 !VERTICAL
		        
		ICOL = ICOL+11 

	
C	Rear power & independent cars:
C	----------------
	DO 80 IP=1,NTRAIN(3)   
		NCHTR(ICOL+1 ) = 1 !!BODY
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 1
		
		NCHTR(ICOL+6 ) = 6 !FRONT BOGIE
		NCHTR(ICOL+7 ) = 7
		NCHTR(ICOL+8 ) = 1
		NCHTR(ICOL+9 ) = 6
		NCHTR(ICOL+10) = 1
		
		NCHTR(ICOL+11) = 6 !FRONT BOGIE WHEEL		
		NCHTR(ICOL+12) = 6
		NCHTR(ICOL+13) = 5
		NCHTR(ICOL+14) = 9	
		NCHTR(ICOL+15) = 9
		NCHTR(ICOL+16) = 8

		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+11 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+12 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+13 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+14 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+15 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+16 !VERTICAL
				
		NCHTR(ICOL+17) = 17 !REAR BOGIE
		NCHTR(ICOL+18) = 18
		NCHTR(ICOL+19) = 1
		NCHTR(ICOL+20) = 17
		NCHTR(ICOL+21) = 1
		
		NCHTR(ICOL+22) = 6 !REAR BOGIE WHEEL	
		NCHTR(ICOL+23) = 6
		NCHTR(ICOL+24) = 5
		NCHTR(ICOL+25) = 9
		NCHTR(ICOL+26) = 9
		NCHTR(ICOL+27) = 8
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+22 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+23 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+24 !VERTICAL
		
		IWHEEL = IWHEEL + 1
		LWHDOF(1,IWHEEL) = ICOL+25 !LATERAL
		LWHDOF(2,IWHEEL) = ICOL+26 !ROLLING
		LWHDOF(3,IWHEEL) = ICOL+27 !VERTICAL
		
		ICOL = ICOL+27
   80 CONTINUE


      RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE TRAINPROP(VEL,NCARB,NACAB,NBOGI,NWHEEL,IWRIN,WPOS,
	1				       PCAB,CADI,ADDI,ALSU,PRSS,PBOG,PRPS,PRWH,WRCO)    
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO READ & WRITE PROPERTIES OF TRAIN (CAR BODY,SUSPENSION,BOGIE)
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	NBOGI = Total number of bogies     
C	NWHEEL= Total number of wheels     

C     TOTAL NUMBER OF CAR BODY = NCARB = NTRAIN(1) + NTRAIN(2) + NTRAIN(3)
C     NUMBER OF ARTIGULATED CAR BODY = NACAB = NTRAIN(2)
      
C	-------
C	OUTPUT:
C	-------
C	NTRAIN(3)	   = Train composition stored in an array
C	PCAB(NCARB,4)  = Car body mechanical properties (MASS,Jx,Jy,Jz)		
C	CADI(NCARB,7)  = Dimension parameters of cars (sA,sB,qA,qB,hC,h1A,h1B) 
C	ADDI(NACAB,5)  = Additional Dimension Parameters for Articulated Cars (b3,h5A,h5B,h6A,h6B)
C	ALSU(NACAB-1,6)= Longitudinal Suspensions of Articulated Cars  =(kTH,kTV,cX)
C	PRSS(NBOGI,9)  = Secondary Suspensions (k2H,c2H,k2V,c2V,k2X,c2X,c2RX,b2,h2) 
C	PBOG(NBOGI,4)  = Bogie of whole train (MASS,Jx,Jy,Jz)				
C	PRPS(NWHEEL,9) = Primary Suspensions (k1H,c1H,k1V,c1V,k1X,c1X,b1,h3,t)
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jy,Jz,rW,b0)	
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwY,cwY)
C	WPOS(NWHEEL)   = Wheel position			
C	-------------------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

C     FOR WHEEL-RAIL CONTACT MECHANISM
	COMMON /WRCT/ Ewr,WMU,RDR,RDW,CANG,SH,D11,D22,D23,D33,DHV,DV		
C
      DIMENSION NTRAIN(3),PCAB(NCARB,4),CADI(NCARB,7),ADDI(NACAB,5)
	DIMENSION ALSU(NACAB-1,6),PRSS(NBOGI,9),PBOG(NBOGI,4)
	DIMENSION PRPS(NWHEEL,9),PRWH(NWHEEL,6),WRCO(NWHEEL,6)
	DIMENSION WPOS(NWHEEL)
C	----------------------------------------------------------------------  


C	-------------------------------------------------
      WRITE(ITO,1005)
      WRITE(10,1005)
      WRITE(100,900)
 
 1005 FORMAT(/,1X,'READ TRAIN PROPERTIES (CAR,SUSPENSION,BOGIE,WHEEL)')
  900	FORMAT (//1H#,18X,49(1H*)/1H#,18X,1H*,47X,1H*/
     +1H#,18X,1H*,1X,45HTRAIN PROPERTIES (CAR,SUSPENSION,BOGIE,WHEEL),
     +1X,1H*/1H#,18X,1H*,47X,1H*/1H#,18X,49(1H*))
C	-------------------------------------------------


C	Read Title for Train Properties
	READ(ITI,*)


C	---------------------
C	1b.CAR BODY PROPERTIES 
C	---------------------
	READ(ITI,*)
	WRITE(100,1100)
      WRITE(ITO,*) ('Read & write car body properties') 
      WRITE(10,*) ('Read & write car body properties') 

C	Read Number of Typical Car Bodies
	READ(ITI,*) NTCARB

	CALL GENTRAIN(NCARB,4,NTCARB,PCAB) 

      DO 110 I=1,NCARB
        WRITE(100,1120) I,(PCAB(I,J),J=1,4)
  110 CONTINUE

 1100 FORMAT(//,1X,'PROPERTIES OF CAR BODY',/
     +1X,'Car No.      Mass     Rolling moment   Yawing moment  Nodding 
     +moment')
 1110 FORMAT(4F10.0)
 1120 FORMAT(I3,4X,4E15.6)
 
 
C	-----------------------------------
C	2b.DIMMENSION PARAMETERS OF CAR BODY:
C	-----------------------------------
	READ(ITI,*)
      WRITE(100,1210)
      WRITE(ITO,*) ('Read & write dimension parameters of car bodies')
      WRITE(10,*) ('Read & write dimension parameters of car bodies')

C	Read Number of Typical Sets of Dimmension Parameters 
	READ(ITI,*) NTDIM

	CALL GENTRAIN(NCARB,7,NTDIM,CADI)

      DO 120 I=1,NCARB
        WRITE(100,1220) I,(CADI(I,J),J=1,7)
  120 CONTINUE
  
 1210 FORMAT(//,1X,'DIMENSION PARAMETERS OF CAR BODIES (in m)',/
     +1X,'Distance btw. current car body center & front bogie center. sA 
     +',/
     +1X,'Distance btw. current car body center & rear bogie center. .sB 
     +',/
     +1X,'Distance from the car body center to front side of car. . . qA
     +',/
     +1X,'Distance from the car body center to rear side of car . . . qB
     +',/
     +1X,'Height of car body.	. .	. .	. . . . .	. .	. .	. . . . . . hC
     +',/
     +1X,'Vertical dist. btw.car body center & front 2nd suspension  h1A 
     +',/
     +1X,'Vertical dist. btw.car body center & rear 2nd suspension. .h1B 
     +',//
     +'CarNo.',2X,'sA',7X,'sB',8X,'qA',8X,'qB',8X,'hC',7X,'h1A',7X,'h1B'
     +)

 1220 FORMAT(I3,2X,7(1X,F8.3))  
 
 
  
C	-----------------------------------------------------------
C	3b.ADDITIONAL DIMMENSION PARAMETERS OF ARTICULATED CAR BODY:
C	-----------------------------------------------------------
	IF(NACAB.GE.3) THEN	

	READ(ITI,*)
      WRITE(100,1310)
      WRITE(ITO,*) ('Read & write additional dimension parameters for Articulated Cars')
      WRITE(10,*) ('Read & write additional dimension parameters for Articulated Cars')

C	Read Number of Typical Sets of Additional Dimmension Parameters 
	READ(ITI,*) NTADIM

	CALL GENTRAIN(NACAB,5,NTADIM,ADDI)

      DO 130 I=1,NACAB
        WRITE(100,1320) I,(ADDI(I,J),J=1,5)
  130 CONTINUE

 1310 FORMAT(//,1X,'ADDITIONAL DIMENSIONS OF ARTICULATED CARS (in m)',/
     +1X,'Half lateral span of longitudinal dampers. . . . . . . . . b3
     +'/
     +1X,'Vert.dist.btw.car body center & front upper longi.dampers .h5A 
     +'/
     +1X,'Vert.dist.btw.car body center & back upper longi.dampers. .h5B 
     +'/
     +1X,'Vert.dist.btw.car body center & front lower longi.dampers .h6A 
     +'/
     +1X,'Vert.dist.btw.car body center & back lower longi.dampers. .h6B 
     +'//
     +'ACarNo.',2X,'b3',6X,'h5A',6X,'h5B',6X,'h6A',6X,'h6B')

 1320 FORMAT(I3,2X,5(1X,F8.3))
   
  
C	----------------------------------------------
C	4b.LONGITUDINAL SUSPENSIONS OF ARTICULATED CARS:
C	----------------------------------------------
	READ(ITI,*)
      WRITE(100,1410)
      WRITE(ITO,*) ('Read & write longitudinal suspensions of Articulated Cars')
      WRITE(10,*) ('Read & write longitudinal suspensions of Articulated Cars')

C	Read Number of Typical Sets of Longitudinal Suspensions
	READ(ITI,*) NTLSUS

	CALL GENTRAIN(NACAB-1,3,NTLSUS,ALSU)

      DO 140 I=1,NACAB-1
        WRITE(100,1420) I,(ALSU(I,J),J=1,3)
	  DO 140 J=4,6
		ALSU(I,J) = 0.0  ! Assign the nominal zero values
  140 CONTINUE

 1410 FORMAT(//,1X,'LONGITUDINAL SUSPENSIONS OF ARTICULATED CAR BODY ',/
     +1X,'Transverse stiffness of central hinge . . . . . kTH '/
     +1X,'Vertical stiffness of central hinge. . . .  . . kTV '/
     +1X,'Damping coefficient of longitudinal damper. . . .cX '/
     +1X,'ACarNo.        kTH            kTV           cX')

 1420 FORMAT(I3,4X,3E15.4)
 
	ENDIF !	IF(NACAB.GE.3)
	
	
	
C	-----------------------
C	5b.SECONDARY SUSPENSIONS:
C	-----------------------
	READ(ITI,*)
      WRITE(100,1500)
      WRITE(ITO,*) ('Read & write secondary suspensions')
      WRITE(10,*) ('Read & write secondary suspensions')

C	Read Number of Typical Secondary Suspensions
	READ(ITI,*) NTSSUS

	CALL GENTRAIN(NBOGI,9,NTSSUS,PRSS)
      DO 150 I=1,NBOGI
        WRITE(100,1520) I,(PRSS(I,J),J=1,9)
  150 CONTINUE

 1500 FORMAT(//,1X,'SECONDARY SUSPENSIONS PARAMETERS',/
     +1X,'Lateral Stiffness  of one of 2 springs in each bogie . .k2H '/
     +1X,'Lateral Damping   of one of 2 dampers in each bogie . . c2H '/
     +1X,'Vertical Stiffness  of one of 2 springs in each bogie . k2V '/
     +1X,'Vertical Damping   of one of 2 dampers in each bogie . .c2V '/
     +1X,'Longitudinal spring stiffness . . . . . . . . . . . . . k2X '/
     +1X,'Longitudinal damping  . . . . . .  . . . . . .  . . . . c2X '/
     +1X,'Rotational damping (for rolling velocity). . . . . .  .c2RX '/
     +1X,'Half lateral span of the secondary suspension. . .  . . .b2 '/
     +1X,'Vertical dist. btw. bogie center and second suspension. .h2 ',
     +//1X,'2ndSuspNo.',2X,'k2H',9X,'c2H',12X,'k2V',12X,'c2V',12X,'k2X',
     +12X,'c2X',12X,'c2RX',12x,'b2',12x,'h2')

 1510 FORMAT(9F10.0)
 1520 FORMAT(I3,4X,9E15.6)
   
  
C	-----------------------
C	6b.PROPERTIES OF BOGIES:
C	-----------------------  
      READ(ITI,*)
	WRITE(100,1600)
      WRITE(ITO,*) ('Read & write properties of bogie')
      WRITE(10,*) ('Read & write properties of bogie')

C	Read Number of Typical Properties of Bogie
	READ(ITI,*) NTPBOG

	CALL GENTRAIN(NBOGI,4,NTPBOG,PBOG)

      DO 160 I=1,NBOGI
        WRITE(100,1620) I,(PBOG(I,J),J=1,4)
  160 CONTINUE

 1600 FORMAT(//,1X,'PROPERTIES OF BOGIES',/
     +1X,'BogieNo.',2x,'Mass',10x,'RollingMoment',2x,'YawingMoment',3x,'
     +NoddingMoment')
 1610 FORMAT(4F10.0)
 1620 FORMAT(I3,4X,4E15.6)
   
  
C	----------------------
C	7b.PRIMARY SUSPENSIONS:
C	----------------------
	READ(ITI,*)
      WRITE(100,1700)
      WRITE(ITO,*) ('Read & write primary suspensions')
      WRITE(10,*) ('Read & write primary suspensions')

C	Read Number of Typical Primary Suspensions
	READ(ITI,*) NTPSUS

	CALL GENTRAIN(NWHEEL,9,NTPSUS,PRPS)

      DO 170 I=1,NWHEEL
        WRITE(100,1720)I,(PRPS(I,J),J=1,9)
  170 CONTINUE
  
 1700 FORMAT(//,1X,'PRIMARY SUSPENSION PAREMETERS',/
     +1X,'Lateral stiffness. . . . . . . . . . . . . .k1H',/
     +1X,'Lateral damping coeficient . . . . . . . . .c1H',/
     +1X,'Vertical stiffness . . . . . . . . . . . . .k1V',/
     +1X,'Vertical damping coeficient. . . . . . . . .c1V',/
     +1X,'Longitudinal spring stiffness . . . . . . . k1X '/
     +1X,'Longitudinal damping  . . . . . .  . . . . .c1X '/
     +1X,'Half lateral span of primary springs. . . . .b1',/
     +1X,'Vert.dist.btw. bogie center & axle center . .h3',/
     +1X,'Half lateral span of front bogie . . . . . . .t',//
     +1X,'1stSuspNo.',2X,'k1H',9X,'c1H',9X,'k1V',9X,'c1V',9X,'k1X',
     +9X,'c1X',9X,'b1',9X,'h3',10X,'t')
 1710 FORMAT(9F10.0)
 1720 FORMAT(I3,4X,9E12.4)

  
C	---------------------------
C	8b.WHEEL-AXLE SET PROPERTIES: 
C	---------------------------
      READ(ITI,*)
	WRITE(100,1800)
      WRITE(ITO,*) ('Read & write wheel-axle set properties')
      WRITE(10,*) ('Read & write wheel-axle set properties')

C	Read Number of Typical Wheel-Axle Set Properties
	READ(ITI,*) NTWHPR

	CALL GENTRAIN(NWHEEL,6,NTWHPR,PRWH)

      DO 180 I=1,NWHEEL
        WRITE(100,1820) I,(PRWH(I,J),J=1,6)
  180 CONTINUE

 1800 FORMAT(//,1X,'PROPERTIES OF WHEEL-AXLE SET',/
     +1X,'WheelNo.',2x,'Mass',10x,'RollingMoment',2x,'YawingMoment', 
     +3x,'NoddingMoment',2x,'WheelRadious',3x,'WheelDistance')
 1810 FORMAT(6F10.0)
 1820 FORMAT(I3,4X,6E15.6)
  
  
C	---------------------------------------------------
C	9b.WHEEL-RAIL CONTACT MECHANISM (Three-dimmensional): WRCO(kwX,cwX,kwH,cwH,kwV,cwV)
C	---------------------------------------------------
	READ(ITI,*)
      WRITE(ITO,*) ('Read & write 3D Wheel-Rail contact Mechanism')
      WRITE(10,*) ('Read & write 3D Wheel-Rail contact Mechanism')
C
C	NTWRCO: Number of Typical Wheel-Rail contact Mechanism
C	IWRIN : Idenfifier for Input of Wheel-Rail contact Mechanism
C	   = 0: Input the contact stiffness & damping directly (OLD)
C	   = 1: Input the wheel-rail properties and calculate contact stiffness & damping
c	READ(ITI,*)NTWRCO,IWRIN !changed next 13Mar08 (read IWRIN in REMOVN)
C	READ(ITI,*)NTWRCO	!Moved down 13Apr08

	IF(IWRIN.EQ.0) THEN
	  WRITE(100,1900)
	  READ(ITI,*) NTWRCO		
	CALL GENTRAIN(NWHEEL,6,NTWRCO,WRCO)

      DO 190 I=1,NWHEEL
	    WRCO(I,4) = WRCO(I,4)/VEL   !cwH = {1.022284*E**(1/3)*P**(2/3)}/V(m/sec)
        WRITE(100,1920) I,(WRCO(I,J),J=1,6)
  190 CONTINUE

 1900 FORMAT(//,1X,'WHEEL-RAIL CONTACT MECHANISM PARAMETERS',/
     +1X,'Longitudinal spring stiffness . . . . . . . . . . . . . kwX '/
     +1X,'Longitudinal damping  . . . . . .  . . . . . .  . . . . cwX '/
     +1X,'Lateral Stiffness  of wheel-rail contact. . . . . .  . .kwH '/
     +1X,'Lateral Damping  of wheel-rail contact . .  . . . . . . cwH '/
     +1X,'Vertical Stiffness of wheel-rail contact . . . . . . . .kwV '/
     +1X,'Vertical Damping  of wheel-rail contact . . . . . . . . cwV',/
     +/1X,'WheelNo.',6X,'kwX',9X,'cwX',12X,'kwH',12X,'cwH',12X,'kwV',
     +12X,'cwV')
 1910 FORMAT(6F10.0)
 1920 FORMAT(I3,4X,6E15.6)


	ELSEIF(IWRIN.EQ.1) THEN !automatically calculate contact propeties
c	    -Input the wheel-rail properties (E,v,Rr, static/dynamic wheel force P)
c	    -Calculate contact stiffness & damping kwH,cwH,kwV
		READ(ITI,*) Ewr,WMU,RDR,CANG
C			Ewr = Elastic modulus of rail and wheel material	
C			WMU = Poisson ratio 
C			RDR = Radius of rail head
C			CANG= Contact angle between the wheel and rail
	    WRITE(100,2010) Ewr,WMU,RDR,CANG  
	ENDIF

 2010 FORMAT(//,1X,'MECHANICAL PROPERTIES OF WHEEL-RAIL CONTACTS',/
     +1X,'Elastic modulus of rail and wheel material (Ewr) = ',E11.4/
     +1X,'Poisson ratio . . . . . . . . . . . . . . .(WMU) = ',F11.4/
     +1X,'Radius of rail head . . . . . . . . . . . .(RDR) = ',F11.4/
     +1X,'Contact angle between the wheel and rail. (CANG) = ',F11.4)
     
     
C	----------------------------------------------
C	10b.WHEEL POSITION:
C	---------------------------------------------- 
      READ(ITI,*)
	WRITE(100,1950)
      WRITE(ITO,*) ('Read & write wheel positions')
      WRITE(10,*) ('Read & write wheel positions')
	 
	DO 195 I=1,NWHEEL
c		Read the relative wheel positions
		READ(ITI,*) IWH,WPOS(I)	
C		Write the relative wheel positions		 
		WRITE(100,1970) (I,WPOS(I))

c		Calculate absolute wheel position	
		WPOS(I) = WPOS(I)-BIDIS
  195	CONTINUE
C	----------------------------------------------
C	----------------------------------------------

 1950 FORMAT(//,1X,'RELATIVE WHEEL POSITION',/
     +1X,'WheelNo.	   Position relative to 1st wheel')
 1960 FORMAT(2F10.0)
 1970 FORMAT(I4,6X,E15.4)
 
      
      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE TRAINAXLEL(NTRAIN,NCARB,NBOGI,NWHEEL,PCAB,CADI,PBOG,PRWH,SCALE,FAXL)	
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	------------------------------------------------------------------------------------------------
C	PROGRAM TO CALCULATE STATIC AXLE FORCE (VERTICAL FORCE FROM WHEEL POSITION DUE TO TRAIN WEIGHT)
C	------------------------------------------------------------------------------------------------
C	INPUT:
C	------
C	PCAB(NCARB,4) = Car body mechanical properties (MASS,Jx,Jz,Jy)
C	CADI(NCARB,7) = Dimension parameters of cars (sA,sB,qA,qB,hC,h1A,h1B) 
C	PBOG(NBOGI,4) = Bogie of whole train (MASS,Jx,Jz,Jy)
C	PRWH(NWHEEL,6)= Wheel-axle set properties (MASS,Jx,Jz,Jy,rW,b0)
C	-------
C	OUTPUT:
C	-------
C	FAXL(NWHEEL)  = Static axle force (vertical force due to train weights acting onto bridge)

C	INDEX:
C	------
C	NACAB		    = Number of articulated car bodies = IT2
C	NCARB		    = Total number of car bodies = IT1+IT2+IT3
C	NBOGI		    = Total number of bogies     = IT1*2+(IT2+1)+IT3*2
C	NWHEEL	    = Total number of wheels     = NBOGI*2
C	SCALE		    = Gravity acceleration (g)
C	-------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION NTRAIN(3)
	DIMENSION PCAB(NCARB,4),CADI(NCARB,7),PBOG(NBOGI,4),PRWH(NWHEEL,6)
      DIMENSION FAXL(NWHEEL)
C	-----------------------------------------------------------------------------------------
C	CALCULATE TRAIN/RAIL INTERACTION LUMPED FORCE APPLIED TO STRUCTURE AT WHEEL POSITION{FSS}
C	-----------------------------------------------------------------------------------------
	CALL CLEARA(FAXL,NWHEEL)

      IT1 = NTRAIN(1) !FRONT POWER CAR & INDEPENDENT CARS
      IT2 = NTRAIN(2) !ARTIGULATED CARS
      IT3 = NTRAIN(3) !REAR  POWER CAR & INDEPENDENT CARS
      
      
C	Front Power Cars:
C	-----------------
	IF(IT1.GT.0)THEN
       DO 100 I=1,IT1		 
	  SCN  =  CADI(I,1) + CADI(I,2)		
	    	   
C	  Two Bogies:
	  DO 100 KB=1,2
	    JB = 2*(I-1)+KB			
		IF(KB.EQ.1) THEN
			WCN =  (CADI(I,2)/SCN)*PCAB(I,1)	
		ELSEIF(KB.EQ.2)THEN
			WCN =  (CADI(I,1)/SCN)*PCAB(I,1)	
		ENDIF 
		 
C	    Two wheels of each bogie:
	    DO 100 KW=1,2
		  JW = 2*(JB-1)+KW		
		  FAXL(JW)=(PRWH(JW,1)+0.5*PBOG(JB,1)+0.5*WCN)*SCALE 
  100  CONTINUE
	ENDIF



C	Articulated Cars (for case of Articulated train):
C	-------------------------------------------------
	IF(IT2.GE.3)THEN	
	DO 200 I = IT1+1,IT1+IT2			
	  SCN  =  CADI(I,1) + CADI(I,2)		
	  IF(I.EQ.(IT1+1)) THEN
		WCN =  (CADI(I,2)/SCN)*PCAB(I,1)
	  ELSE
		SCN1=  CADI(I-1,1) + CADI(I-1,2) 
		WCN =  (CADI(I-1,1)/SCN1)*PCAB(I-1,1)
	1		  +(CADI(I,2)/SCN)*PCAB(I,1)	
	  ENDIF

	  JB = IT1 + I			 
C		 Two wheels of front bogie:
	     DO 200 KW=1,2
		   JW = 2*(JB-1) + KW	
		   FAXL(JW)=(PRWH(JW,1)+0.5*PBOG(JB,1)+0.5*WCN)*SCALE 
  200 CONTINUE

C	The last bogie of rear transition car:
C	--------------------------------------
	ICAR = IT1+IT2
	SCN  = CADI(ICAR,1) + CADI(ICAR,2)	
	WCN =  (CADI(ICAR,1)/SCN)*PCAB(ICAR,1) 

	JB = 2*IT1+IT2+1		 
C	   Two wheels of last bogie:
	   DO 250 KW=1,2
		   JW = 2*(JB-1) + KW	
		   FAXL(JW)=(PRWH(JW,1)+0.5*PBOG(JB,1)+0.5*WCN)*SCALE 
  250 CONTINUE
	ENDIF		!<-----IF(IT2.GE.3), for case of Articulated train



C	Rear Power Cars:
C	-----------------
	IF(IT3.GT.0)THEN
       DO 300 I=IT1+IT2+1,IT1+IT2+IT3		
	  SCN  =  CADI(I,1) + CADI(I,2)  

C	  Two Bogies:
	  DO 300 KB=1,2
		JB = 2*I-IT2-1+KB 	
		IF(KB.EQ.1) THEN
			WCN =  (CADI(I,2)/SCN)*PCAB(I,1)
		ELSEIF(KB.EQ.2)THEN
			WCN =  (CADI(I,1)/SCN)*PCAB(I,1)	
		ENDIF 
									  
C	    Two wheels of each bogie:
	    DO 300 KW=1,2
		  JW = 2*(JB-1)+KW		
		  FAXL(JW)=(PRWH(JW,1)+0.5*PBOG(JB,1)+0.5*WCN)*SCALE 

  300  CONTINUE	
	ENDIF


C	WRITE OUT THE CALCULATED STATIC AXLE FORCES
C	-------------------------------------------
	WRITE(100,1100)
	 
	DO 500 I=1,NWHEEL	 
		WRITE(100,1200) I,FAXL(I)
  500	CONTINUE

 1100 FORMAT(//,1X,'CALCULATED STATIC AXLE FORCES',/
     +1X,'Wheel No.	   Vertical Axle Force')
 1200 FORMAT(I4,6X,E15.4)

	RETURN
      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE RAILCONTCT(NWHEEL,PRWH,FAXL,WRCO,VEL)             
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO GENERATE WHEEL-RAIL CONTACT STIFFNESS & DAMPING FROM MECHANICAL PROP
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	FAXL(NWHEEL)   = Static axle force (vertical force due to train weights acting onto bridge)
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jy,Jz,rW,b0)	!updated 19Jan08
C	NWHEEL= Total number of wheels     = NBOGI*2
C	VEL	  =  Running speed (m/s)
C	Ewr	  = Elastic modulus of rail and wheel material	
C	WMU   = Poisson ratio 
C	RDR   = Radius of rail head
C	CANG  = Contact angle between the wheel and rail
C	-------
C	OUTPUT:
C	-------
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwY,cwY)	
C	-------------------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

	COMMON /WRCT/ Ewr,WMU,RDR,RDW,CANG,SH,D11,D22,D23,D33,DHV,DV
C	-------------------------------------------------------------------------------------------
	DIMENSION PRWH(NWHEEL,6),FAXL(NWHEEL),WRCO(NWHEEL,6)
	
C	---------------------------------------------------
C	WHEEL-RAIL CONTACT MECHANISM (Three-dimmensional): 
C	---------------------------------------------------
      WRITE(ISO,1900)
      
C	Define constants for all wheels
	CKV=(1.5**(1./3.))*((Ewr/(1.-WMU**2))**(2./3.))*(RDR**(1./6.))
	CC1=((1.-WMU)**(2./3.))/((1.+WMU)**(1./3.))
	C22=4.18028
	CCH=((3./16.)**(2./3.))*CC1*(Ewr**(1./3.))*C22/VEL

      DO 90 IW=1,NWHEEL
	  WRCO(IW,1) = 0.0
	  WRCO(IW,2) = 0.0
	  WRCO(IW,5) = CKV*(FAXL(IW)**(1./3.))*(PRWH(IW,5)**(1./6.))
	  WRCO(IW,6) = 0.0
	  WRCO(IW,3) = 0.8134*WRCO(IW,5)
	  WRCO(IW,4) = CCH*(FAXL(IW)**(2./3.))*((PRWH(IW,5)*RDR)**(1./3.))
        WRITE(ISO,1920) IW,(WRCO(IW,J),J=1,6)
   90 CONTINUE

 1900 FORMAT(//,1X,'WHEEL-RAIL CONTACT MECHANISM PARAMETERS',/
     +1X,'Longitudinal spring stiffness . . . . . . . . . . . . . kwX '/
     +1X,'Longitudinal damping  . . . . . .  . . . . . .  . . . . cwX '/
     +1X,'Lateral Stiffness  of wheel-rail contact. . . . . .  . .kwH '/
     +1X,'Lateral Damping  of wheel-rail contact . .  . . . . . . cwH '/
     +1X,'Vertical Stiffness of wheel-rail contact . . . . . . . .kwV '/
     +1X,'Vertical Damping  of wheel-rail contact . . . . . . . . cwV',/
     +/1X,'WheelNo.',6X,'kwX',9X,'cwX',12X,'kwH',12X,'cwH',12X,'kwV',
     +12X,'cwV')

 1910 FORMAT(6F10.0)
 1920 FORMAT(I3,4X,6E15.6)

	RETURN
      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE CONTCTRAN(NWHEEL,PRWH,WRCO,WCSTIF)             
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	-----------------------------------------------------------------------
	DIMENSION PRWH(NWHEEL,6),WRCO(NWHEEL,6),WCSTIF(3,2,NWHEEL) !WCSTIF  3=component  2=Stif&Damping

      DO 90 IW=1,NWHEEL
      
          WCSTIF(1,1,IW) = 2.0*WRCO(IW,3)
          WCSTIF(1,2,IW) = 2.0*WRCO(IW,4)
          
          WCSTIF(2,1,IW) = 2.0*(WRCO(IW,5)*(PRWH(IW,6)**2)+WRCO(IW,3)*(PRWH(IW,5)**2))
          WCSTIF(2,2,IW) = 2.0*(WRCO(IW,6)*(PRWH(IW,6)**2)+WRCO(IW,4)*(PRWH(IW,5)**2)) 
          
          WCSTIF(3,1,IW) = 2.0*WRCO(IW,5)
          WCSTIF(3,2,IW) = 2.0*WRCO(IW,6)
      
   90 CONTINUE


	RETURN
      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE TRAINMAT(NTRAIN,NCARB,NACAB,NBOGI,NWHEEL,AMTR,ACTR,AKTR,
     1                    PCAB,CADI,ADDI,ALSU,PRSS,PBOG,PRPS,PRWH,WRCO)     
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO FORM MASS, DAMPING & STIFFNESS MATRICES OF THE TRAIN SYSTEM
C	----------------------------------------------------------------------
C	INPUT:
C	------
C	NTRAIN(3)	   = Train composition stored in an array
C	PCAB(NCARB,4)  = Car body mechanical properties (MASS,Jx,Jy,Jz) !change posi Jy,Jz 19Jan08
C	CADI(NCARB,7)  = Dimension parameters of cars (sA,sB,qA,qB,hC,h1A,h1B) 
C	ADDI(NACAB,5)  = Additional Dimension Parameters for Articulated Cars (b3,h5A,h5B,h6A,h6B)
C	ALSU(NACAB-1,6)= Longitudinal Suspensions of Articulated Cars  =(kTH,kTV,cX)
C	PRSS(NBOGI,9)  = Secondary Suspensions (k2H,c2H,k2V,c2V,k2X,c2X,c2RX,b2,h2) 
C	PBOG(NBOGI,4)  = Bogie of whole train (MASS,Jx,Jy,Jz)			!change posi Jy,Jz 19Jan08
C	PRPS(NWHEEL,9) = Primary Suspensions (k1H,c1H,k1V,c1V,k1X,c1X,b1,h3,t)
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jy,Jz,rW,b0)
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwZ,cwZ)
C
C	NACAB = Number of articulated car bodies     = NTRAIN(2)
C	Parameters for whole train system:
C	NCARB  = Total number of car bodies = NPCAB + NACAB
C	NBOGI  = Total number of bogies     = NPBOG + NABOG
C	NWHEEL = Total number of wheels     = NPWHL + NAWHL
C	NEQTR  = Number of equations of stiffness matrix
C	NWKTR  = Number of elements below skyline of upper tria. stiffness /damping matrix

C	OUTPUT:
C	-------
C	AMTR(NEQTR)  = MASS MATRIX OF TRAIN (Lumped mass matrix)
C	ACTR(NWKTR)  = DAMPING MATRIX OF TRAIN
C	AKTR(NWKTR)  = STIFFNESS MATRIX OF TRAIN
C	------------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C
      DIMENSION NTRAIN(3),PCAB(NCARB,1),CADI(NCARB,1),ADDI(NACAB,1)
	DIMENSION ALSU(NACAB-1,1),PRSS(NBOGI,1),PBOG(NBOGI,1)
	DIMENSION PRPS(NWHEEL,1),PRWH(NWHEEL,1),WRCO(NWHEEL,1)	
			
	DIMENSION AMTR(1),ACTR(1),AKTR(1)
C	---------------------------------------------------------------------------

C	------------------------------------------------------------------------------
C 	A.EVALUATION OF SYSTEM MATRICES(MASS, DAMPING, STIFFNESS) OF FRONT POWER CARS:                                                    
C	------------------------------------------------------------------------------
      IPM=0  ! Position pointer for element of mass matrix stored in array
	IPK=0  ! Position pointer for element of stiffness matrix stored in array

      IF(NTRAIN(1).GT.0) THEN
	  DO 10 I=1,NTRAIN(1)
	    CALL POWCMAT(NTRAIN,NCARB,NACAB,NBOGI,NWHEEL,
	1                PCAB,CADI,PRSS,PBOG,PRPS,PRWH,WRCO,
	1		          AMTR,ACTR,AKTR,IPM,IPK,I)   
	    IPM = IPM + 27
	    IPK = IPK + 167
   10	  CONTINUE
      ENDIF
      
C	------------------------------------------------------------
C     B.EVALUATION OF SYSTEM MATRICES(MASS, DAMPING, STIFFNESS) OF               
C	ARTICULATED GROUP (INTERMEDIATE AND PASSENGER CARS):                                    
C	------------------------------------------------------------       
      IF(NTRAIN(2).GE.3) THEN

C	FORM THE MASS MATRIX:
C	---------------------
C	Form the lumped mass matrix of each pair of articulated car (bogie & car body):
C	------------------------------------------------------------------------------	
	DO 20 IA = NTRAIN(1)+1,NTRAIN(1)+NTRAIN(2)	!loop over number of articulated cars
	      JB = NTRAIN(1)+IA		! Corresponding bogie number

	  JW1=2*JB-1                
        AMTR(IPM+1 )=PRWH(JW1,1)
        AMTR(IPM+2 )=PRWH(JW1,2)
        AMTR(IPM+3 )=PRWH(JW1,1)

	  JW2=2*JB                  
        AMTR(IPM+4 )=PRWH(JW2,1)
        AMTR(IPM+5 )=PRWH(JW2,2)
        AMTR(IPM+6 )=PRWH(JW2,1)

C	  Account for bogie:
        AMTR(IPM+7 )= PBOG(JB,1)
        AMTR(IPM+8 )= PBOG(JB,2)
        AMTR(IPM+9 )= PBOG(JB,3)
        AMTR(IPM+10)= PBOG(JB,1)
        AMTR(IPM+11)= PBOG(JB,4)

C	  Account for car bodies
        AMTR(IPM+12)= PCAB(IA,1)
        AMTR(IPM+13)= PCAB(IA,2)
        AMTR(IPM+14)= PCAB(IA,3)
        AMTR(IPM+15)= PCAB(IA,1)
        AMTR(IPM+16)= PCAB(IA,4)

	  IPM=IPM+16          
   20 CONTINUE
C
C	The last independent bogie and its wheel-axles:
C	-----------------------------------------------
	  JB = NTRAIN(1)*2+NTRAIN(2)+1

C	  Account for bogie:
        AMTR(IPM+1 )= PBOG(JB,1)
        AMTR(IPM+2 )= PBOG(JB,2)
        AMTR(IPM+3 )= PBOG(JB,3)
        AMTR(IPM+4 )= PBOG(JB,1)
        AMTR(IPM+5 )= PBOG(JB,4)

C	  Account for 1st wheel of last bogie:
	  JW1=2*JB-1               
        AMTR(IPM+6 )=PRWH(JW1,1)
        AMTR(IPM+7 )=PRWH(JW1,2)
        AMTR(IPM+8 )=PRWH(JW1,1)

C	  Account for 2nd wheel of last bogie:
	  JW2=2*JB                  
        AMTR(IPM+9 )=PRWH(JW2,1)
        AMTR(IPM+10)=PRWH(JW2,2)
        AMTR(IPM+11)=PRWH(JW2,1)
  
	  IPM=IPM+11						

C	FORM THE STIFF & DAMPING MATRICES:
C	---------------------------------------------------------------------------------------
C	1.The first transition car (sub-matrices of the first pair (car body & bogie) (K11,C11))
C	2.The intermidate passenger cars (1 <n <N) and the last transition car (n=N):	
C	 (Form the tri-diagonal matrices by looping over each pair (bogie & car body))
C	---------------------------------------------------------------------------------------
	CALL ARTGMAT(NTRAIN,NCARB,NACAB,NBOGI,NWHEEL,
	1            CADI,ADDI,ALSU,PRSS,PRPS,PRWH,WRCO,
	1			  ACTR,AKTR,IPK)
C	------------------------------------------------------------------------
C    3.For the last bogie: (n = N+1) and its coupling with the last (Nth) car
C	The matrices must be formulated seperately since its chacteristics
C	------------------------------------------------------------------------
	IA = NTRAIN(1)+NTRAIN(2)
	JB = NTRAIN(1)+IA+1	    
	JW1= 2*JB-1			
	JW2= 2*JB				
	IAC= IA-NTRAIN(1)		

C-----  Stiffness MATRIX:
	AKTR(IPK+1 )= 2.*PRSS(JB,1)+2.*PRPS(JW1,1)+2.*PRPS(JW2,1) !-->bogie, !K(N+1,N+1)
	AKTR(IPK+2 )= 0.									!K(N,N+1)
	AKTR(IPK+3 )= 0.

	AKTR(IPK+4 )= 2.*PRSS(JB,1)*CADI(IA,2)				!K(N,N+1)
	AKTR(IPK+5 )= 2.*PRSS(JB,1)*CADI(IA,7)
	AKTR(IPK+6 )=-2.*PRSS(JB,1)

      AKTR(IPK+7 )= 2.*PRSS(JB,1)*PRSS(JB,9)**2 +			!K(N+1,N+1)
	1			  2.*PRSS(JB,3)*PRSS(JB,8)**2 +
	2		      2.*PRPS(JW1,1)*PRPS(JW1,8)**2 +
     3			  2.*PRPS(JW2,1)*PRPS(JW2,8)**2 +
     4			  2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
     5			  2.*PRPS(JW2,3)*PRPS(JW2,7)**2 

      AKTR(IPK+8 )= 2.*PRSS(JB,1)*PRSS(JB,9) 
	1			 -2.*PRPS(JW1,1)*PRPS(JW1,8)
	2			 -2.*PRPS(JW2,1)*PRPS(JW2,8)

	AKTR(IPK+9 )= 0.										!K(N,N+1)
	AKTR(IPK+10)= 0.

	AKTR(IPK+11)= 2.*PRSS(JB,1)*CADI(IA,2)*PRSS(JB,9)

	AKTR(IPK+12)= 2.*PRSS(JB,1)*CADI(IA,7)*PRSS(JB,9)
	1			 -2.*PRSS(JB,3)*PRSS(JB,9)**2

	AKTR(IPK+13)=-2.*PRSS(JB,1)*PRSS(JB,9)

 
      AKTR(IPK+14)= 2.*PRPS(JW1,1)*PRPS(JW1,9)**2 +			!K(N+1,N+1)
     1			  2.*PRPS(JW2,1)*PRPS(JW2,9)**2
CC      Taken from 1st transition bogie  dof 64-65 
	AKTR(IPK+15)=0.0 !-2.*PRSS(JB,1)*CADI(IA,1)*CADI(IA,6)
c	1			      -ALSU(IAC,1)*CADI(IA,2)*CADI(IA,7) 

	AKTR(IPK+16)=0.0 ! 2.*PRSS(JB,1)*CADI(IA,1) + ALSU(IAC,1)*CADI(IA,2)
CC

	AKTR(IPK+17)= 2.*PRSS(JB,3) + 2.*PRPS(JW1,3) + 2.*PRPS(JW2,3) 

	AKTR(IPK+18)= 0.
	AKTR(IPK+19)= 0.
	AKTR(IPK+20)= 0.

	AKTR(IPK+21)=-2.*PRSS(JB,1)*CADI(IA,2)					!K(N,N+1)
	AKTR(IPK+22)=-2.*PRSS(JB,1)
	

      AKTR(IPK+23)= 2.*PRPS(JW1,3)*PRPS(JW1,9)**2 +			!K(N+1,N+1)
     1			  2.*PRPS(JW2,3)*PRPS(JW2,9)**2	
CC	
      AKTR(IPK+24)= 0.
	AKTR(IPK+25)= 0.
	AKTR(IPK+26)= 0.
	AKTR(IPK+27)= 0.	
CC      Taken from 1st transition bogie  dof 77-78 		
      AKTR(IPK+28)= 0.
	AKTR(IPK+29)= 0.0 !2.*PRSS(JB,3)*CADI(IA,1)	
CC				



	AKTR(IPK+30)= 2.*PRPS(JW1,1) + 2.*WRCO(JW1,3)	!--> 1st wheel

      AKTR(IPK+31)= 0.
      AKTR(IPK+32)= 0.

	AKTR(IPK+33)= 2.*PRPS(JW1,1)*PRPS(JW1,9)

	AKTR(IPK+34)= 2.*PRPS(JW1,1)*PRPS(JW1,8)

	AKTR(IPK+35)=-2.*PRPS(JW1,1)					

      AKTR(IPK+36)= 2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
	1			  2.*WRCO(JW1,5)*PRWH(JW1,6)**2 +
	2			  2.*WRCO(JW1,3)*PRWH(JW1,5)**2

	AKTR(IPK+37)=-2.*WRCO(JW1,3)*PRWH(JW1,5)
	AKTR(IPK+38)= 0.
	AKTR(IPK+39)= 0.
      AKTR(IPK+40)= 0.

	AKTR(IPK+41)=-2.*PRPS(JW1,3)*PRPS(JW1,7)**2		

	AKTR(IPK+42)= 2.*PRPS(JW1,3) + 2.*WRCO(JW1,5)

	AKTR(IPK+43)= 0.
	AKTR(IPK+44)= 0.

	AKTR(IPK+45)=-2.*PRPS(JW1,3)*PRPS(JW1,9)

	AKTR(IPK+46)=-2.*PRPS(JW1,3)

	AKTR(IPK+47)= 2.*PRPS(JW2,1) + 2.*WRCO(JW2,3)	!--> 2nd wheel

	AKTR(IPK+48)= 0.
	AKTR(IPK+49)= 0.

	AKTR(IPK+50)= 0.
	AKTR(IPK+51)= 0.
	AKTR(IPK+52)= 0.					

	AKTR(IPK+53)=-2.*PRPS(JW2,1)*PRPS(JW2,9)
	AKTR(IPK+54)= 2.*PRPS(JW2,1)*PRPS(JW2,8)

	AKTR(IPK+55)=-2.*PRPS(JW2,1)			

	AKTR(IPK+56)= 2.*PRPS(JW2,3)*PRPS(JW2,7)**2 +
	1			  2.*WRCO(JW2,5)*PRWH(JW2,6)**2 +
	2			  2.*WRCO(JW2,3)*PRWH(JW2,5)**2		

	AKTR(IPK+57)=-2.*WRCO(JW2,3)*PRWH(JW2,5)

	AKTR(IPK+58)= 0.
	AKTR(IPK+59)= 0.
	AKTR(IPK+60)= 0.
	AKTR(IPK+61)= 0.
	AKTR(IPK+62)= 0.
	AKTR(IPK+63)= 0.

	AKTR(IPK+64)=-2.*PRPS(JW2,3)*PRPS(JW2,7)**2

	AKTR(IPK+65)= 2.*PRPS(JW2,3) + 2.*WRCO(JW2,5)
	AKTR(IPK+66)= 0.
	AKTR(IPK+67)= 0.
	AKTR(IPK+68)= 0.
	AKTR(IPK+69)= 0.
	AKTR(IPK+70)= 0.

	AKTR(IPK+71)= 2.*PRPS(JW2,3)*PRPS(JW2,9)
	AKTR(IPK+72)=-2.*PRPS(JW2,3)					

C-----  Damping MATRIX:
	ACTR(IPK+1 )= 2.*PRSS(JB,2)+2.*PRPS(JW1,2)+2.*PRPS(JW2,2)	!-->bogie, !C(N+1,N+1)
	ACTR(IPK+2 )= 0.											!C(N,N+1)
	ACTR(IPK+3 )= 0.

	ACTR(IPK+4 )= 2.*PRSS(JB,2)*CADI(IA,2)						!C(N,N+1)
	ACTR(IPK+5 )= 2.*PRSS(JB,2)*CADI(IA,7)
	ACTR(IPK+6 )=-2.*PRSS(JB,2)

      ACTR(IPK+7 )= 2.*PRSS(JB,2)*PRSS(JB,9)**2 +					!C(N+1,N+1)
	1			  2.*PRSS(JB,4)*PRSS(JB,8)**2 +
	2		      2.*PRPS(JW1,2)*PRPS(JW1,8)**2 +
     3			  2.*PRPS(JW2,2)*PRPS(JW2,8)**2 +
     4			  2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
     5			  2.*PRPS(JW2,4)*PRPS(JW2,7)**2 

      ACTR(IPK+8 )= 2.*PRSS(JB,2)*PRSS(JB,9) 
	1			 -2.*PRPS(JW1,2)*PRPS(JW1,8)
	2			 -2.*PRPS(JW2,2)*PRPS(JW2,8)

	ACTR(IPK+9 )= 0.											!C(N,N+1)
	ACTR(IPK+10)= 0.

	ACTR(IPK+11)= 2.*PRSS(JB,2)*CADI(IA,2)*PRSS(JB,9)

	ACTR(IPK+12)= 2.*PRSS(JB,2)*CADI(IA,7)*PRSS(JB,9)
	1			 -2.*PRSS(JB,4)*PRSS(JB,9)**2

	ACTR(IPK+13)=-2.*PRSS(JB,2)*PRSS(JB,9)

      ACTR(IPK+14)= 2.*PRPS(JW1,2)*PRPS(JW1,9)**2 +			!C(N+1,N+1)
     1			  2.*PRPS(JW2,2)*PRPS(JW2,9)**2
     
CC      Taken from 1st transition bogie  dof 64-65 
	ACTR(IPK+15)= 0.0 !-2.*PRSS(JB,2)*CADI(IA,1)*CADI(IA,6)
c	1			      -ALSU(IAC,4)*CADI(IA,2)*CADI(IA,7)
	ACTR(IPK+16)= 0.0 !2.*PRSS(JB,2)*CADI(IA,1) + ALSU(IAC,4)*CADI(IA,2)
CC	
     
     
     

	ACTR(IPK+17)= 2.*PRSS(JB,4) + 2.*PRPS(JW1,4) + 2.*PRPS(JW2,4) 

	ACTR(IPK+18)= 0.
	ACTR(IPK+19)= 0.
	ACTR(IPK+20)= 0.

	ACTR(IPK+21)=-2.*PRSS(JB,2)*CADI(IA,2)					!C(N,N+1)
	ACTR(IPK+22)=-2.*PRSS(JB,2)
      ACTR(IPK+23)= 2.*PRPS(JW1,4)*PRPS(JW1,9)**2 +			!C(N+1,N+1)
     1			  2.*PRPS(JW2,4)*PRPS(JW2,9)**2		
CC	
      ACTR(IPK+24)= 0.
	ACTR(IPK+25)= 0.
	ACTR(IPK+26)= 0.
	ACTR(IPK+27)= 0.	
CC      Taken from 1st transition bogie  dof 77-78 		
      ACTR(IPK+28)= 0.
	ACTR(IPK+29)= 0.0 !2.*PRSS(JB,4)*CADI(IA,1)
CC	
     		

	ACTR(IPK+30)= 2.*PRPS(JW1,2) + 2.*WRCO(JW1,4)	!--> 1st wheel

      ACTR(IPK+31)= 0.
      ACTR(IPK+32)= 0.

	ACTR(IPK+33)= 2.*PRPS(JW1,2)*PRPS(JW1,9)

	ACTR(IPK+34)= 2.*PRPS(JW1,2)*PRPS(JW1,8)

	ACTR(IPK+35)=-2.*PRPS(JW1,2)					

      ACTR(IPK+36)= 2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
	1			  2.*WRCO(JW1,6)*PRWH(JW1,6)**2 +
	2			  2.*WRCO(JW1,4)*PRWH(JW1,5)**2

	ACTR(IPK+37)=-2.*WRCO(JW1,4)*PRWH(JW1,5)
	ACTR(IPK+38)= 0.
	ACTR(IPK+39)= 0.
      ACTR(IPK+40)= 0.

	ACTR(IPK+41)=-2.*PRPS(JW1,4)*PRPS(JW1,7)**2		

	ACTR(IPK+42)= 2.*PRPS(JW1,4) + 2.*WRCO(JW1,5)

	ACTR(IPK+43)= 0.
	ACTR(IPK+44)= 0.

	ACTR(IPK+45)=-2.*PRPS(JW1,4)*PRPS(JW1,9)

	ACTR(IPK+46)=-2.*PRPS(JW1,4)

	ACTR(IPK+47)= 2.*PRPS(JW2,2) + 2.*WRCO(JW2,4)	!--> 2nd wheel

	ACTR(IPK+48)= 0.
	ACTR(IPK+49)= 0.

	ACTR(IPK+50)= 0.
	ACTR(IPK+51)= 0.
	ACTR(IPK+52)= 0.					

	ACTR(IPK+53)=-2.*PRPS(JW2,2)*PRPS(JW2,9)
	ACTR(IPK+54)= 2.*PRPS(JW2,2)*PRPS(JW2,8)

	ACTR(IPK+55)=-2.*PRPS(JW2,2)			

	ACTR(IPK+56)= 2.*PRPS(JW2,4)*PRPS(JW2,7)**2 +
	1			  2.*WRCO(JW2,6)*PRWH(JW2,6)**2 +
	2			  2.*WRCO(JW2,4)*PRWH(JW2,5)**2		

	ACTR(IPK+57)=-2.*WRCO(JW2,4)*PRWH(JW2,5)

	ACTR(IPK+58)= 0.
	ACTR(IPK+59)= 0.
	ACTR(IPK+60)= 0.
	ACTR(IPK+61)= 0.
	ACTR(IPK+62)= 0.
	ACTR(IPK+63)= 0.

	ACTR(IPK+64)=-2.*PRPS(JW2,4)*PRPS(JW2,7)**2

	ACTR(IPK+65)= 2.*PRPS(JW2,4) + 2.*WRCO(JW2,6)
	ACTR(IPK+66)= 0.
	ACTR(IPK+67)= 0.
	ACTR(IPK+68)= 0.
	ACTR(IPK+69)= 0.
	ACTR(IPK+70)= 0.

	ACTR(IPK+71)= 2.*PRPS(JW2,4)*PRPS(JW2,9)
	ACTR(IPK+72)=-2.*PRPS(JW2,4)

	IPK=IPK+72									!STOP HERE

      ENDIF   !(NTRAIN(2).GE.3)
C	-----------------------------------------------------------------------------
C	C.EVALUATION OF SYSTEM MATRICES(MASS, DAMPING, STIFFNESS) OF REAR POWER CARS                                                    .
C	-----------------------------------------------------------------------------
	IF(NTRAIN(3).GT.0) THEN 
       DO 50 I = NTRAIN(1)+NTRAIN(2)+1,NTRAIN(1)+NTRAIN(2)+NTRAIN(3)    !Loop over number of rear power cars

	    CALL POWCMAT(NTRAIN,NCARB,NACAB,NBOGI,NWHEEL,
	1                PCAB,CADI,PRSS,PBOG,PRPS,PRWH,WRCO,
	1		          AMTR,ACTR,AKTR,IPM,IPK,I)          
       
	   IPM = IPM + 27
	   IPK = IPK + 167
   50	 CONTINUE
      ENDIF
C     
      RETURN
 1000 FORMAT(1X,'NUMBER OF INTERMEDIATE AND PASSENGER CARS SHOULD BE',/,
     &1X,'EQUAL TO 0 OR GREATER THAN AND EQUAL TO 3',/,
     &1X,'IT2 = ',I5)
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE POWCMAT(NTRAIN,NCARB,NACAB,NBOGI,NWHEEL,
     1                  PCAB,CADI,PRSS,PBOG,PRPS,PRWH,WRCO,
	1				    AMTR,ACTR,AKTR,IPM,IPK,I)     
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO FORM MASS, DAMPING & STIFFNESS MATRICES OF A POWER CAR
C	----------------------------------------------------------------------
C	INPUT:
c	------
C	PCAB(NCARB,4)  = Car body mechanical properties (MASS,Jx,Jy,Jz)	  !updated 19Jan08
C	PBOG(NBOGI,4)  = Bogie of whole train (MASS,Jx,Jy,Jz)			  !updated 19Jan08
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jy,Jz,rW,b0)  !updated 19Jan08
C      IPM  :Position pointer for element of mass matrix stored in array
C	 IPK  :Position pointer for element of stiffness matrix stored in array
C	 I	  :Current car number
C	OUTPUT:
c	-------
C	AMTR(NEQTR)  = MASS MATRIX OF TRAIN (Lumped mass matrix)
C	ACTR(NWKTR)  = DAMPING MATRIX OF TRAIN
C	AKTR(NWKTR)  = STIFFNESS MATRIX OF TRAIN
C	----------------------------------------------------------------------
      DIMENSION NTRAIN(3)
      DIMENSION PCAB(NCARB,1),CADI(NCARB,1),PRSS(NBOGI,1),PBOG(NBOGI,1)
	DIMENSION PRPS(NWHEEL,1),PRWH(NWHEEL,1),WRCO(NWHEEL,1)	
	DIMENSION AMTR(1),ACTR(1),AKTR(1)
C	--------------------------------------------------------------------------

      IT1 = NTRAIN(1)
      IT2 = NTRAIN(2)
      IT3 = NTRAIN(3)
      

	IF(I.LE.IT1) THEN
C	 POINTERS FOR FRONT POWER CARS
	 JB1=2*I-1                  
	   JW1=2*JB1-1              
	   JW2=2*JB1                
	 JB2=2*I                    
	   JW3=2*JB2-1             
	   JW4=2*JB2                

	ELSE IF(I.GT.(IT1+IT2).AND.I.LE.NCARB) THEN
C	 POINTERS FOR REAR POWER CARS
	 JB1=2*IT1+(IT2+1)+2*(I-IT1-IT2)-1      
	   JW1=2*JB1-1							
	   JW2=2*JB1							
	 JB2=2*IT1+(IT2+1)+2*(I-IT1-IT2)		
	   JW3=2*JB2-1							
	   JW4=2*JB2							
	ENDIF

C	Form the lumped mass matrix of each power car:
C	---------------------------------------------
C	  Account for car bodies
        AMTR(IPM+1 )=PCAB(I,1)
        AMTR(IPM+2 )=PCAB(I,2)
        AMTR(IPM+3 )=PCAB(I,3)
        AMTR(IPM+4 )=PCAB(I,1)
        AMTR(IPM+5 )=PCAB(I,4)

C	  Account for 1st bogie:
        AMTR(IPM+6 )=PBOG(JB1,1)
        AMTR(IPM+7 )=PBOG(JB1,2)
        AMTR(IPM+8 )=PBOG(JB1,3)
        AMTR(IPM+9 )=PBOG(JB1,1)
        AMTR(IPM+10)=PBOG(JB1,4)

C	  Account for 1st wheel of 1st bogie:
        AMTR(IPM+11)=PRWH(JW1,1)
        AMTR(IPM+12)=PRWH(JW1,2)
        AMTR(IPM+13)=PRWH(JW1,1)

C	  Account for 2nd wheel of 1st bogie:
        AMTR(IPM+14)=PRWH(JW2,1)
        AMTR(IPM+15)=PRWH(JW2,2)
        AMTR(IPM+16)=PRWH(JW2,1)

C	  Account for 2nd bogie:
        AMTR(IPM+17)=PBOG(JB2,1)
        AMTR(IPM+18)=PBOG(JB2,2)
        AMTR(IPM+19)=PBOG(JB2,3)
        AMTR(IPM+20)=PBOG(JB2,1)
        AMTR(IPM+21)=PBOG(JB2,4)

C	  Account for 1st wheel of 2nd bogie:
        AMTR(IPM+22)=PRWH(JW3,1)
        AMTR(IPM+23)=PRWH(JW3,2)
        AMTR(IPM+24)=PRWH(JW3,1)

C	  Account for 2nd wheel of 2nd bogie:
        AMTR(IPM+25)=PRWH(JW4,1)
        AMTR(IPM+26)=PRWH(JW4,2)
        AMTR(IPM+27)=PRWH(JW4,1)

C	Form the stiffness matrix of each power car:
C	-------------------------------------------       
        AKTR(IPK+1 )= 2.*(PRSS(JB1,1) + PRSS(JB2,1))
        AKTR(IPK+2 )= 2.*PRSS(JB1,1)*CADI(I,6)**2 + 
	1				2.*PRSS(JB2,1)*CADI(I,7)**2 +
     2	 		    2.*PRSS(JB1,3)*PRSS(JB1,8)**2 + 
	3				2.*PRSS(JB2,3)*PRSS(JB2,8)**2
        AKTR(IPK+3 )=-2.*PRSS(JB1,1)*CADI(I,6)-2.*PRSS(JB2,1)*CADI(I,7)

        AKTR(IPK+4 )= 2.*PRSS(JB1,1)*CADI(I,1)**2 + 
	1				2.*PRSS(JB2,1)*CADI(I,2)**2

        AKTR(IPK+5 )= 2.*(PRSS(JB1,3) + PRSS(JB2,3))

        AKTR(IPK+6 )= 2.*PRSS(JB1,3)*CADI(I,1)**2 + 
	1				2.*PRSS(JB2,3)*CADI(I,2)**2

        AKTR(IPK+7 )= 2.*PRSS(JB1,1) + 2.*PRPS(JW1,1) + 2.*PRPS(JW2,1)

        AKTR(IPK+8 )= 0.
        AKTR(IPK+9 )= 0.

C        AKTR(IPK+10)=-2.*PRSS(JB1,1)*CADI(I,1)  !Changed to next 9Feb07 for consistent coordinate
        AKTR(IPK+10)= 2.*PRSS(JB1,1)*CADI(I,1)

        AKTR(IPK+11)= 2.*PRSS(JB1,1)*CADI(I,6)

        AKTR(IPK+12)=-2.*PRSS(JB1,1)

        AKTR(IPK+13)= 2.*PRSS(JB1,1)*PRSS(JB1,9)**2 +
	1			    2.*PRSS(JB1,3)*PRSS(JB1,8)**2 +
	2				2.*PRPS(JW1,1)*PRPS(JW1,8)**2 +
     3				2.*PRPS(JW2,1)*PRPS(JW2,8)**2 +
     4				2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
     5				2.*PRPS(JW2,3)*PRPS(JW2,7)**2   				 

        AKTR(IPK+14)= 2.*PRSS(JB1,1)*PRSS(JB1,9) 
	1			   -2.*PRPS(JW1,1)*PRPS(JW1,8)
	2			   -2.*PRPS(JW2,1)*PRPS(JW2,8)
        AKTR(IPK+15)= 0.
        AKTR(IPK+16)= 0.

C        AKTR(IPK+17)=-2.*PRSS(JB1,1)*CADI(I,1)*PRSS(JB1,9)  !Changed to next 9Feb07 for consistent coordinate
        AKTR(IPK+17)= 2.*PRSS(JB1,1)*CADI(I,1)*PRSS(JB1,9)

        AKTR(IPK+18)= 2.*PRSS(JB1,1)*CADI(I,6)*PRSS(JB1,9)-
     1			    2.*PRSS(JB1,3)*PRSS(JB1,8)**2

	  AKTR(IPK+19)=-2.*PRSS(JB1,1)*PRSS(JB1,9)

        AKTR(IPK+20)= 2.*PRPS(JW1,1)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,1)*PRPS(JW2,9)**2

        AKTR(IPK+21)= 2.*PRSS(JB1,3) + 2.*PRPS(JW1,3) + 2.*PRPS(JW2,3)
        AKTR(IPK+22)= 0.
        AKTR(IPK+23)= 0.
        AKTR(IPK+24)= 0.
	  
c	  AKTR(IPK+25)= 2.*PRSS(JB1,3)*CADI(I,1)  !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+25)=-2.*PRSS(JB1,3)*CADI(I,1)

        AKTR(IPK+26)=-2.*PRSS(JB1,3)

        AKTR(IPK+27)= 2.*PRPS(JW1,3)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,3)*PRPS(JW2,9)**2		

	  AKTR(IPK+28)= 2.*PRPS(JW1,1) + 2.*WRCO(JW1,3)	!--> 1st wheel

        AKTR(IPK+29)= 0.
        AKTR(IPK+30)= 0.

c	  AKTR(IPK+31)=-2.*PRPS(JW1,1)*PRPS(JW1,9) !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+31)= 2.*PRPS(JW1,1)*PRPS(JW1,9)

	  AKTR(IPK+32)= 2.*PRPS(JW1,1)*PRPS(JW1,8)

	  AKTR(IPK+33)=-2.*PRPS(JW1,1)	

        AKTR(IPK+34)= 2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
	1				2.*WRCO(JW1,5)*PRWH(JW1,6)**2 +
	2				2.*WRCO(JW1,3)*PRWH(JW1,5)**2

	  AKTR(IPK+35)=-2.*WRCO(JW1,3)*PRWH(JW1,5)
	  AKTR(IPK+36)= 0.
	  AKTR(IPK+37)= 0.
        AKTR(IPK+38)= 0.

	  AKTR(IPK+39)=-2.*PRPS(JW1,3)*PRPS(JW1,7)**2		

	  AKTR(IPK+40)= 2.*PRPS(JW1,3) + 2.*WRCO(JW1,5)

	  AKTR(IPK+41)= 0.
	  AKTR(IPK+42)= 0.

c	  AKTR(IPK+43)= 2.*PRPS(JW1,3)*PRPS(JW1,9) !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+43)=-2.*PRPS(JW1,3)*PRPS(JW1,9)

	  AKTR(IPK+44)=-2.*PRPS(JW1,3)

	  AKTR(IPK+45)= 2.*PRPS(JW2,1) + 2.*WRCO(JW2,3)	!--> 2nd wheel

	  AKTR(IPK+46)= 0.
	  AKTR(IPK+47)= 0.
	  AKTR(IPK+48)= 0.
	  AKTR(IPK+49)= 0.
	  AKTR(IPK+50)= 0.					

c	  AKTR(IPK+51)= 2.*PRPS(JW2,1)*PRPS(JW2,9)  !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+51)=-2.*PRPS(JW2,1)*PRPS(JW2,9)

	  AKTR(IPK+52)= 2.*PRPS(JW2,1)*PRPS(JW2,8)

	  AKTR(IPK+53)=-2.*PRPS(JW2,1)			

	  AKTR(IPK+54)= 2.*PRPS(JW2,3)*PRPS(JW2,7)**2 +
	1				2.*WRCO(JW2,5)*PRWH(JW2,6)**2 +
	2				2.*WRCO(JW2,3)*PRWH(JW2,5)**2		

	  AKTR(IPK+55)=-2.*WRCO(JW2,3)*PRWH(JW2,5)

	  AKTR(IPK+56)= 0.
	  AKTR(IPK+57)= 0.
	  AKTR(IPK+58)= 0.
	  AKTR(IPK+59)= 0.
	  AKTR(IPK+60)= 0.
	  AKTR(IPK+61)= 0.

	  AKTR(IPK+62)=-2.*PRPS(JW2,3)*PRPS(JW2,7)**2

	  AKTR(IPK+63)= 2.*PRPS(JW2,3) + 2.*WRCO(JW2,5)
	  AKTR(IPK+64)= 0.
	  AKTR(IPK+65)= 0.
	  AKTR(IPK+66)= 0.
	  AKTR(IPK+67)= 0.
	  AKTR(IPK+68)= 0.

c	  AKTR(IPK+69)=-2.*PRPS(JW2,3)*PRPS(JW2,9) !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+69)= 2.*PRPS(JW2,3)*PRPS(JW2,9)

	  AKTR(IPK+70)=-2.*PRPS(JW2,3)

	  AKTR(IPK+71)= 2.*PRSS(JB2,1) + 2.*PRPS(JW3,1) + 2.*PRPS(JW4,1) !-->2nd bogie
	  AKTR(IPK+72)= 0.
	  AKTR(IPK+73)= 0.
	  AKTR(IPK+74)= 0.
	  AKTR(IPK+75)= 0.
	  AKTR(IPK+76)= 0.
	  AKTR(IPK+77)= 0.
	  AKTR(IPK+78)= 0.
	  AKTR(IPK+79)= 0.
	  AKTR(IPK+80)= 0.
	  AKTR(IPK+81)= 0.
	  AKTR(IPK+82)= 0.
	  AKTR(IPK+83)= 0.
	  AKTR(IPK+84)= 0.

c        AKTR(IPK+85)= 2.*PRSS(JB2,1)*CADI(I,2) !Changed to next 9Feb07 for consistent coordinate
        AKTR(IPK+85)=-2.*PRSS(JB2,1)*CADI(I,2)

        AKTR(IPK+86)= 2.*PRSS(JB2,1)*CADI(I,7)

        AKTR(IPK+87)=-2.*PRSS(JB2,1)

        AKTR(IPK+88)= 2.*PRSS(JB2,1)*PRSS(JB2,9)**2 +
	1			    2.*PRSS(JB2,3)*PRSS(JB2,8)**2 +
	2				2.*PRPS(JW3,1)*PRPS(JW3,8)**2 +
     3				2.*PRPS(JW4,1)*PRPS(JW4,8)**2 +
     4				2.*PRPS(JW3,3)*PRPS(JW3,7)**2 +
     5				2.*PRPS(JW4,3)*PRPS(JW4,7)**2   				 

        AKTR(IPK+89)= 2.*PRSS(JB2,1)*PRSS(JB2,9) 
	1			   -2.*PRPS(JW3,1)*PRPS(JW3,8)
	2			   -2.*PRPS(JW4,1)*PRPS(JW4,8)

	  AKTR(IPK+90)= 0.
	  AKTR(IPK+91)= 0.
	  AKTR(IPK+92)= 0.
	  AKTR(IPK+93)= 0.
	  AKTR(IPK+94)= 0.
	  AKTR(IPK+95)= 0.
	  AKTR(IPK+96)= 0.
	  AKTR(IPK+97)= 0.
	  AKTR(IPK+98)= 0.
	  AKTR(IPK+99)= 0.
	  AKTR(IPK+100)=0.
	  AKTR(IPK+101)=0.
	  AKTR(IPK+102)=0.

c        AKTR(IPK+103)= 2.*PRSS(JB2,1)*CADI(I,2)*PRSS(JB2,9) !Changed to next 9Feb07 for consistent coordinate
        AKTR(IPK+103)=-2.*PRSS(JB2,1)*CADI(I,2)*PRSS(JB2,9)

        AKTR(IPK+104)= 2.*PRSS(JB2,1)*CADI(I,7)*PRSS(JB2,9)-
     1			     2.*PRSS(JB2,3)*PRSS(JB2,8)**2

	  AKTR(IPK+105)=-2.*PRSS(JB2,1)*PRSS(JB2,9)		

        AKTR(IPK+106)= 2.*PRPS(JW3,1)*PRPS(JW3,9)**2 +
     1				 2.*PRPS(JW4,1)*PRPS(JW4,9)**2

        AKTR(IPK+107)= 2.*PRSS(JB2,3) + 2.*PRPS(JW3,3) + 2.*PRPS(JW4,3)

	  AKTR(IPK+108)=0.
	  AKTR(IPK+109)=0.
	  AKTR(IPK+110)=0.
	  AKTR(IPK+111)=0.
	  AKTR(IPK+112)=0.
	  AKTR(IPK+113)=0.
	  AKTR(IPK+114)=0.
	  AKTR(IPK+115)=0.
	  AKTR(IPK+116)=0.
	  AKTR(IPK+117)=0.
	  AKTR(IPK+118)=0.
	  AKTR(IPK+119)=0.
	  AKTR(IPK+120)=0.
	  AKTR(IPK+121)=0.

c	   AKTR(IPK+122)=-2.*PRSS(JB2,3)*CADI(I,2)  !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+122)= 2.*PRSS(JB2,3)*CADI(I,2)

        AKTR(IPK+123)=-2.*PRSS(JB2,3)		

        AKTR(IPK+124)= 2.*PRPS(JW3,3)*PRPS(JW3,9)**2 +		
     1				 2.*PRPS(JW4,3)*PRPS(JW4,9)**2	

	  AKTR(IPK+125)= 2.*PRPS(JW3,1) + 2.*WRCO(JW3,3)	!--> 3rd wheel

        AKTR(IPK+126)= 0.
        AKTR(IPK+127)= 0.

c	   AKTR(IPK+128)=-2.*PRPS(JW3,1)*PRPS(JW3,9) !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+128)= 2.*PRPS(JW3,1)*PRPS(JW3,9)

	  AKTR(IPK+129)= 2.*PRPS(JW3,1)*PRPS(JW3,8)

	  AKTR(IPK+130)=-2.*PRPS(JW3,1)					

        AKTR(IPK+131)= 2.*PRPS(JW3,3)*PRPS(JW3,7)**2 +
	1				 2.*WRCO(JW3,5)*PRWH(JW3,6)**2 +
	2				 2.*WRCO(JW3,3)*PRWH(JW3,5)**2

	  AKTR(IPK+132)=-2.*WRCO(JW3,3)*PRWH(JW3,5)
	  AKTR(IPK+133)= 0.
	  AKTR(IPK+134)= 0.
        AKTR(IPK+135)= 0.

	  AKTR(IPK+136)=-2.*PRPS(JW3,3)*PRPS(JW3,7)**2		

	  AKTR(IPK+137)= 2.*PRPS(JW3,3) + 2.*WRCO(JW3,5)

	  AKTR(IPK+138)= 0.
	  AKTR(IPK+139)= 0.

c 	   AKTR(IPK+140)= 2.*PRPS(JW3,3)*PRPS(JW3,9) !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+140)=-2.*PRPS(JW3,3)*PRPS(JW3,9)

	  AKTR(IPK+141)=-2.*PRPS(JW3,3)

	  AKTR(IPK+142)= 2.*PRPS(JW4,1) + 2.*WRCO(JW4,3)	!--> 4th wheel

	  AKTR(IPK+143)= 0.
	  AKTR(IPK+144)= 0.

	  AKTR(IPK+145)= 0.
	  AKTR(IPK+146)= 0.
	  AKTR(IPK+147)= 0.					

c	   AKTR(IPK+148)= 2.*PRPS(JW4,1)*PRPS(JW4,9) !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+148)=-2.*PRPS(JW4,1)*PRPS(JW4,9) 

	  AKTR(IPK+149)= 2.*PRPS(JW4,1)*PRPS(JW4,8)

	  AKTR(IPK+150)=-2.*PRPS(JW4,1)			

	  AKTR(IPK+151)= 2.*PRPS(JW4,3)*PRPS(JW4,7)**2 +
	1				 2.*WRCO(JW4,5)*PRWH(JW4,6)**2 +
	2				 2.*WRCO(JW4,3)*PRWH(JW4,5)**2		

	  AKTR(IPK+152)=-2.*WRCO(JW4,3)*PRWH(JW4,5)

	  AKTR(IPK+153)= 0.
	  AKTR(IPK+154)= 0.
	  AKTR(IPK+155)= 0.
	  AKTR(IPK+156)= 0.
	  AKTR(IPK+157)= 0.
	  AKTR(IPK+158)= 0.

	  AKTR(IPK+159)=-2.*PRPS(JW4,3)*PRPS(JW4,7)**2

	  AKTR(IPK+160)= 2.*PRPS(JW4,3) + 2.*WRCO(JW4,5)
	  AKTR(IPK+161)= 0.
	  AKTR(IPK+162)= 0.
	  AKTR(IPK+163)= 0.
	  AKTR(IPK+164)= 0.
	  AKTR(IPK+165)= 0.

c 	   AKTR(IPK+166)=-2.*PRPS(JW4,3)*PRPS(JW4,9)  !Changed to next 9Feb07 for consistent coordinate
	  AKTR(IPK+166)= 2.*PRPS(JW4,3)*PRPS(JW4,9)

	  AKTR(IPK+167)=-2.*PRPS(JW4,3)

C	Form the damping matrix of each power car:
C	------------------------------------------- 
        ACTR(IPK+1 )= 2.*(PRSS(JB1,2) + PRSS(JB2,2))	
        ACTR(IPK+2 )= 2.*PRSS(JB1,2)*CADI(I,6)**2 + 
	1				2.*PRSS(JB2,2)*CADI(I,7)**2 +
     2	 		    2.*PRSS(JB1,4)*PRSS(JB1,8)**2 + 
	3				2.*PRSS(JB2,4)*PRSS(JB2,8)**2
        ACTR(IPK+3 )=-2.*PRSS(JB1,2)*CADI(I,6)-2.*PRSS(JB2,2)*CADI(I,7)

        ACTR(IPK+4 )= 2.*PRSS(JB1,2)*CADI(I,1)**2 + 
	1				2.*PRSS(JB2,2)*CADI(I,2)**2

        ACTR(IPK+5 )= 2.*(PRSS(JB1,4) + PRSS(JB2,4))

        ACTR(IPK+6 )= 2.*PRSS(JB1,4)*CADI(I,1)**2 + 
	1				2.*PRSS(JB2,4)*CADI(I,2)**2

        ACTR(IPK+7 )= 2.*PRSS(JB1,2) + 2.*PRPS(JW1,2) + 2.*PRPS(JW2,2)

        ACTR(IPK+8 )= 0.
        ACTR(IPK+9 )= 0.

c        ACTR(IPK+10)=-2.*PRSS(JB1,2)*CADI(I,1) !Changed to next 9Feb07 for consistent coordinate
        ACTR(IPK+10)= 2.*PRSS(JB1,2)*CADI(I,1) 

        ACTR(IPK+11)= 2.*PRSS(JB1,2)*CADI(I,6)

        ACTR(IPK+12)=-2.*PRSS(JB1,2)

        ACTR(IPK+13)= 2.*PRSS(JB1,2)*PRSS(JB1,9)**2 +
	1			    2.*PRSS(JB1,4)*PRSS(JB1,8)**2 +
	2				2.*PRPS(JW1,2)*PRPS(JW1,8)**2 +
     3				2.*PRPS(JW2,2)*PRPS(JW2,8)**2 +
     4				2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
     5				2.*PRPS(JW2,4)*PRPS(JW2,7)**2   				 

        ACTR(IPK+14)= 2.*PRSS(JB1,2)*PRSS(JB1,9) 
	1			   -2.*PRPS(JW1,2)*PRPS(JW1,8)
	2			   -2.*PRPS(JW2,2)*PRPS(JW2,8)
        ACTR(IPK+15)= 0.
        ACTR(IPK+16)= 0.

c        ACTR(IPK+17)=-2.*PRSS(JB1,2)*CADI(I,1)*PRSS(JB1,9)  !Changed to next 9Feb07 for consistent coordinate
        ACTR(IPK+17)= 2.*PRSS(JB1,2)*CADI(I,1)*PRSS(JB1,9) 

        ACTR(IPK+18)= 2.*PRSS(JB1,2)*CADI(I,6)*PRSS(JB1,9)-
     1			    2.*PRSS(JB1,4)*PRSS(JB1,8)**2

	  ACTR(IPK+19)=-2.*PRSS(JB1,2)*PRSS(JB1,9)

        ACTR(IPK+20)= 2.*PRPS(JW1,2)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,2)*PRPS(JW2,9)**2

        ACTR(IPK+21)= 2.*PRSS(JB1,4) + 2.*PRPS(JW1,4) + 2.*PRPS(JW2,4)
        ACTR(IPK+22)= 0.
        ACTR(IPK+23)= 0.
        ACTR(IPK+24)= 0.
	  
c	  ACTR(IPK+25)= 2.*PRSS(JB1,4)*CADI(I,1)  !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+25)=-2.*PRSS(JB1,4)*CADI(I,1)

        ACTR(IPK+26)=-2.*PRSS(JB1,4)

        ACTR(IPK+27)= 2.*PRPS(JW1,4)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,4)*PRPS(JW2,9)**2		

	  ACTR(IPK+28)= 2.*PRPS(JW1,2) + 2.*WRCO(JW1,4)	!--> 1st wheel

        ACTR(IPK+29)= 0.
        ACTR(IPK+30)= 0.

c	  ACTR(IPK+31)=-2.*PRPS(JW1,2)*PRPS(JW1,9)  !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+31)= 2.*PRPS(JW1,2)*PRPS(JW1,9)

	  ACTR(IPK+32)= 2.*PRPS(JW1,2)*PRPS(JW1,8)

	  ACTR(IPK+33)=-2.*PRPS(JW1,2)	

        ACTR(IPK+34)= 2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
	1				2.*WRCO(JW1,6)*PRWH(JW1,6)**2 +
	2				2.*WRCO(JW1,4)*PRWH(JW1,5)**2

	  ACTR(IPK+35)=-2.*WRCO(JW1,4)*PRWH(JW1,5)
	  ACTR(IPK+36)= 0.
	  ACTR(IPK+37)= 0.
        ACTR(IPK+38)= 0.

	  ACTR(IPK+39)=-2.*PRPS(JW1,4)*PRPS(JW1,7)**2		

	  ACTR(IPK+40)= 2.*PRPS(JW1,4) + 2.*WRCO(JW1,6)

	  ACTR(IPK+41)= 0.
	  ACTR(IPK+42)= 0.

c	  ACTR(IPK+43)= 2.*PRPS(JW1,4)*PRPS(JW1,9) !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+43)=-2.*PRPS(JW1,4)*PRPS(JW1,9)

	  ACTR(IPK+44)=-2.*PRPS(JW1,4)

	  ACTR(IPK+45)= 2.*PRPS(JW2,2) + 2.*WRCO(JW2,4)	!--> 2nd wheel

	  ACTR(IPK+46)= 0.
	  ACTR(IPK+47)= 0.
	  ACTR(IPK+48)= 0.
	  ACTR(IPK+49)= 0.
	  ACTR(IPK+50)= 0.					

c	  ACTR(IPK+51)= 2.*PRPS(JW2,2)*PRPS(JW2,9)  !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+51)=-2.*PRPS(JW2,2)*PRPS(JW2,9)

	  ACTR(IPK+52)= 2.*PRPS(JW2,2)*PRPS(JW2,8)

	  ACTR(IPK+53)=-2.*PRPS(JW2,2)			

	  ACTR(IPK+54)= 2.*PRPS(JW2,4)*PRPS(JW2,7)**2 +
	1				2.*WRCO(JW2,6)*PRWH(JW2,6)**2 +
	2				2.*WRCO(JW2,4)*PRWH(JW2,5)**2		

	  ACTR(IPK+55)=-2.*WRCO(JW2,4)*PRWH(JW2,5)

	  ACTR(IPK+56)= 0.
	  ACTR(IPK+57)= 0.
	  ACTR(IPK+58)= 0.
	  ACTR(IPK+59)= 0.
	  ACTR(IPK+60)= 0.
	  ACTR(IPK+61)= 0.

	  ACTR(IPK+62)=-2.*PRPS(JW2,4)*PRPS(JW2,7)**2

	  ACTR(IPK+63)= 2.*PRPS(JW2,4) + 2.*WRCO(JW2,6)
	  ACTR(IPK+64)= 0.
	  ACTR(IPK+65)= 0.
	  ACTR(IPK+66)= 0.
	  ACTR(IPK+67)= 0.
	  ACTR(IPK+68)= 0.

c	  ACTR(IPK+69)=-2.*PRPS(JW2,4)*PRPS(JW2,9) !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+69)= 2.*PRPS(JW2,4)*PRPS(JW2,9)

	  ACTR(IPK+70)=-2.*PRPS(JW2,4)

	  ACTR(IPK+71)= 2.*PRSS(JB2,2) + 2.*PRPS(JW3,2) + 2.*PRPS(JW4,2) !-->2nd bogie
	  ACTR(IPK+72)= 0.
	  ACTR(IPK+73)= 0.
	  ACTR(IPK+74)= 0.
	  ACTR(IPK+75)= 0.
	  ACTR(IPK+76)= 0.
	  ACTR(IPK+77)= 0.
	  ACTR(IPK+78)= 0.
	  ACTR(IPK+79)= 0.
	  ACTR(IPK+80)= 0.
	  ACTR(IPK+81)= 0.
	  ACTR(IPK+82)= 0.
	  ACTR(IPK+83)= 0.
	  ACTR(IPK+84)= 0.

c        ACTR(IPK+85)= 2.*PRSS(JB2,2)*CADI(I,2) !Changed to next 9Feb07 for consistent coordinate
        ACTR(IPK+85)=-2.*PRSS(JB2,2)*CADI(I,2)

        ACTR(IPK+86)= 2.*PRSS(JB2,2)*CADI(I,7)

        ACTR(IPK+87)=-2.*PRSS(JB2,2)

        ACTR(IPK+88)= 2.*PRSS(JB2,2)*PRSS(JB2,9)**2 +
	1			    2.*PRSS(JB2,4)*PRSS(JB2,8)**2 +
	2				2.*PRPS(JW3,2)*PRPS(JW3,8)**2 +
     3				2.*PRPS(JW4,2)*PRPS(JW4,8)**2 +
     4				2.*PRPS(JW3,4)*PRPS(JW3,7)**2 +
     5				2.*PRPS(JW4,4)*PRPS(JW4,7)**2   				 

        ACTR(IPK+89)= 2.*PRSS(JB2,2)*PRSS(JB2,9) 
	1			   -2.*PRPS(JW3,2)*PRPS(JW3,8)
	2			   -2.*PRPS(JW4,2)*PRPS(JW4,8)

	  ACTR(IPK+90)= 0.
	  ACTR(IPK+91)= 0.
	  ACTR(IPK+92)= 0.
	  ACTR(IPK+93)= 0.
	  ACTR(IPK+94)= 0.
	  ACTR(IPK+95)= 0.
	  ACTR(IPK+96)= 0.
	  ACTR(IPK+97)= 0.
	  ACTR(IPK+98)= 0.
	  ACTR(IPK+99)= 0.
	  ACTR(IPK+100)=0.
	  ACTR(IPK+101)=0.
	  ACTR(IPK+102)=0.

c        ACTR(IPK+103)= 2.*PRSS(JB2,2)*CADI(I,2)*PRSS(JB2,9) !Changed to next 9Feb07 for consistent coordinate
        ACTR(IPK+103)=-2.*PRSS(JB2,2)*CADI(I,2)*PRSS(JB2,9)

        ACTR(IPK+104)= 2.*PRSS(JB2,2)*CADI(I,7)*PRSS(JB2,9)-
     1			     2.*PRSS(JB2,4)*PRSS(JB2,8)**2

	  ACTR(IPK+105)=-2.*PRSS(JB2,2)*PRSS(JB2,9)		

        ACTR(IPK+106)= 2.*PRPS(JW3,2)*PRPS(JW3,9)**2 +
     1				 2.*PRPS(JW4,2)*PRPS(JW4,9)**2

        ACTR(IPK+107)= 2.*PRSS(JB2,4) + 2.*PRPS(JW3,4) + 2.*PRPS(JW4,4)

	  ACTR(IPK+108)=0.
	  ACTR(IPK+109)=0.
	  ACTR(IPK+110)=0.
	  ACTR(IPK+111)=0.
	  ACTR(IPK+112)=0.
	  ACTR(IPK+113)=0.
	  ACTR(IPK+114)=0.
	  ACTR(IPK+115)=0.
	  ACTR(IPK+116)=0.
	  ACTR(IPK+117)=0.
	  ACTR(IPK+118)=0.
	  ACTR(IPK+119)=0.
	  ACTR(IPK+120)=0.
	  ACTR(IPK+121)=0.

c	  ACTR(IPK+122)=-2.*PRSS(JB2,4)*CADI(I,2) !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+122)= 2.*PRSS(JB2,4)*CADI(I,2)

        ACTR(IPK+123)=-2.*PRSS(JB2,4)		

        ACTR(IPK+124)= 2.*PRPS(JW3,4)*PRPS(JW3,9)**2 +		
     1				 2.*PRPS(JW4,4)*PRPS(JW4,9)**2	

	  ACTR(IPK+125)= 2.*PRPS(JW3,2) + 2.*WRCO(JW3,4)	!--> 3rd wheel

        ACTR(IPK+126)= 0.
        ACTR(IPK+127)= 0.

c	  ACTR(IPK+128)=-2.*PRPS(JW3,2)*PRPS(JW3,9)  !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+128)= 2.*PRPS(JW3,2)*PRPS(JW3,9)

	  ACTR(IPK+129)= 2.*PRPS(JW3,2)*PRPS(JW3,8)

	  ACTR(IPK+130)=-2.*PRPS(JW3,2)					

        ACTR(IPK+131)= 2.*PRPS(JW3,4)*PRPS(JW3,7)**2 +
	1				 2.*WRCO(JW3,6)*PRWH(JW3,6)**2 +
	2				 2.*WRCO(JW3,4)*PRWH(JW3,5)**2

	  ACTR(IPK+132)=-2.*WRCO(JW3,4)*PRWH(JW3,5)
	  ACTR(IPK+133)= 0.
	  ACTR(IPK+134)= 0.
        ACTR(IPK+135)= 0.

	  ACTR(IPK+136)=-2.*PRPS(JW3,4)*PRPS(JW3,7)**2		

	  ACTR(IPK+137)= 2.*PRPS(JW3,4) + 2.*WRCO(JW3,6)

	  ACTR(IPK+138)= 0.
	  ACTR(IPK+139)= 0.

c	  ACTR(IPK+140)= 2.*PRPS(JW3,4)*PRPS(JW3,9) !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+140)=-2.*PRPS(JW3,4)*PRPS(JW3,9)

	  ACTR(IPK+141)=-2.*PRPS(JW3,4)

	  ACTR(IPK+142)= 2.*PRPS(JW4,2) + 2.*WRCO(JW4,4)	!--> 4th wheel

	  ACTR(IPK+143)= 0.
	  ACTR(IPK+144)= 0.

	  ACTR(IPK+145)= 0.
	  ACTR(IPK+146)= 0.
	  ACTR(IPK+147)= 0.					

c	  ACTR(IPK+148)= 2.*PRPS(JW4,2)*PRPS(JW4,9)  !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+148)=-2.*PRPS(JW4,2)*PRPS(JW4,9)

	  ACTR(IPK+149)= 2.*PRPS(JW4,2)*PRPS(JW4,8)

	  ACTR(IPK+150)=-2.*PRPS(JW4,2)			

	  ACTR(IPK+151)= 2.*PRPS(JW4,4)*PRPS(JW4,7)**2 +
	1				 2.*WRCO(JW4,6)*PRWH(JW4,6)**2 +
	2				 2.*WRCO(JW4,4)*PRWH(JW4,5)**2		

	  ACTR(IPK+152)=-2.*WRCO(JW4,4)*PRWH(JW4,5)

	  ACTR(IPK+153)= 0.
	  ACTR(IPK+154)= 0.
	  ACTR(IPK+155)= 0.
	  ACTR(IPK+156)= 0.
	  ACTR(IPK+157)= 0.
	  ACTR(IPK+158)= 0.

	  ACTR(IPK+159)=-2.*PRPS(JW4,4)*PRPS(JW4,7)**2

	  ACTR(IPK+160)= 2.*PRPS(JW4,4) + 2.*WRCO(JW4,6)
	  ACTR(IPK+161)= 0.
	  ACTR(IPK+162)= 0.
	  ACTR(IPK+163)= 0.
	  ACTR(IPK+164)= 0.
	  ACTR(IPK+165)= 0.

c	  ACTR(IPK+166)=-2.*PRPS(JW4,4)*PRPS(JW4,9)  !Changed to next 9Feb07 for consistent coordinate
	  ACTR(IPK+166)= 2.*PRPS(JW4,4)*PRPS(JW4,9)

	  ACTR(IPK+167)=-2.*PRPS(JW4,4)
C	<<-----------------END OF POWER CAR MATRIX
	RETURN
	END	
C	
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE ARTGMAT(NTRAIN,NCARB,NACAB,NBOGI,NWHEEL,
     1                  CADI,ADDI,ALSU,PRSS,PRPS,PRWH,WRCO,
	1				    ACTR,AKTR,IPK)            
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO FORM DAMPING & STIFFNESS MATRICES OF THE ARTICULATED GROUP
C	----------------------------------------------------------------------
C	INPUT:
C	------
C	NTRAIN(3)	   = Train composition stored in an array
C	PCAB(NCARB,4)  = Car body mechanical properties (MASS,Jx,Jz,Jy)
C	CADI(NCARB,7)  = Dimension parameters of cars (sA,sB,qA,qB,hC,h1A,h1B) 
C	ADDI(NACAB,5)  = Additional Dimension Parameters for Articulated Cars (b3,h5A,h5B,h6A,h6B)
C	ALSU(NACAB-1,6)= Longitudinal Suspensions of Articulated Cars  =(kTH,kTV,cX)
C	PRSS(NBOGI,9)  = Secondary Suspensions (k2H,c2H,k2V,c2V,k2X,c2X,c2RX,b2,h2) 
C	PBOG(NBOGI,4)  = Bogie of whole train (MASS,Jx,Jz,Jy)
C	PRPS(NWHEEL,9) = Primary Suspensions (k1H,c1H,k1V,c1V,k1X,c1X,b1,h3,t)
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jz,Jy,rW,b0)
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwZ,cwZ)
C	IPK	= Position pointer for element of stiffness matrix stored in array, for last term in Power Car

C	OUTPUT:
C	-------
C	ACTR(NWKTR)  = DAMPING MATRIX OF TRAIN
C	AKTR(NWKTR)  = STIFFNESS MATRIX OF TRAIN
C	IPK	= Position pointer for element of stiffness matrix stored in array (last
C	------------------------------------------------------------------------------------

      DIMENSION NTRAIN(3),CADI(NCARB,1),ADDI(NACAB,1),ALSU(NACAB-1,1)
	DIMENSION PRSS(NBOGI,1),PRPS(NWHEEL,1),PRWH(NWHEEL,1)		
	DIMENSION WRCO(NWHEEL,1)			
	DIMENSION ACTR(1),AKTR(1)  	
	
C	---------------------------------------------------------------------------
C	FORM THE STIFF & DAMPING MATRICES:
C	-------------------------------------------------------
C	1.The first transition car:
C	 Matrices of the first pair (car body & bogie) (K11,C11):
C	-------------------------------------------------------
	  IA = NTRAIN(1)+1		
	  JB = NTRAIN(1)+IA	
	  JW1=2*JB-1            
	  JW2=2*JB              
	  IAC= 1				

C	Stiffness K11

	  AKTR(IPK+1 )= 2.*PRPS(JW1,1) + 2.*WRCO(JW1,3)	!--> 1st wheel
        AKTR(IPK+2 )= 2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
	1				2.*WRCO(JW1,5)*PRWH(JW1,6)**2 +
	2				2.*WRCO(JW1,3)*PRWH(JW1,5)**2
	  AKTR(IPK+3 )=-2.*WRCO(JW1,3)*PRWH(JW1,5)
	  AKTR(IPK+4 )= 2.*PRPS(JW1,3) + 2.*WRCO(JW1,5)

	  AKTR(IPK+5 )= 2.*PRPS(JW2,1) + 2.*WRCO(JW2,3)	!--> 2nd wheel
        AKTR(IPK+6 )= 2.*PRPS(JW2,3)*PRPS(JW2,7)**2 +
	1				2.*WRCO(JW2,5)*PRWH(JW2,6)**2 +
	2				2.*WRCO(JW2,3)*PRWH(JW2,5)**2
	  AKTR(IPK+7 )=-2.*WRCO(JW2,3)*PRWH(JW2,5)
	  AKTR(IPK+8 )= 2.*PRPS(JW2,3) + 2.*WRCO(JW2,5)

	  AKTR(IPK+9 )= 2.*PRSS(JB,1) + 2.*PRPS(JW1,1) + 2.*PRPS(JW2,1) !-->bogie
	  AKTR(IPK+10)= 0.
	  AKTR(IPK+11)= 0.

 	  AKTR(IPK+12)=-2.*PRPS(JW2,1)

	  AKTR(IPK+13)= 0.
	  AKTR(IPK+14)= 0.
 	  AKTR(IPK+15)=-2.*PRPS(JW1,1)		

        AKTR(IPK+16)= 2.*PRSS(JB,1)*PRSS(JB,9)**2 +
	1			    2.*PRSS(JB,3)*PRSS(JB,8)**2 +
	2				2.*PRPS(JW1,1)*PRPS(JW1,8)**2 +
     3				2.*PRPS(JW2,1)*PRPS(JW2,8)**2 +
     4				2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
     5				2.*PRPS(JW2,3)*PRPS(JW2,7)**2   				 

        AKTR(IPK+17)= 2.*PRSS(JB,1)*PRSS(JB,9) 
	1			   -2.*PRPS(JW1,1)*PRPS(JW1,8)
	2			   -2.*PRPS(JW2,1)*PRPS(JW2,8)		

        AKTR(IPK+18)= 0.
	  AKTR(IPK+19)=-2.*PRPS(JW2,3)*PRPS(JW2,7)**2
	  AKTR(IPK+20)= 2.*PRPS(JW2,1)*PRPS(JW2,8)

        AKTR(IPK+21)= 0.
	  AKTR(IPK+22)=-2.*PRPS(JW1,3)*PRPS(JW1,7)**2
	  AKTR(IPK+23)= 2.*PRPS(JW1,1)*PRPS(JW1,8)
	  		
        AKTR(IPK+24)= 2.*PRPS(JW1,1)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,1)*PRPS(JW2,9)**2	
	  AKTR(IPK+25)= 0.
	  AKTR(IPK+26)= 0.
	  AKTR(IPK+27)= 0.
	  AKTR(IPK+28)= 0.
 	  AKTR(IPK+29)=-2.*PRPS(JW2,1)*PRPS(JW2,9)
	  AKTR(IPK+30)= 0.
	  AKTR(IPK+31)= 0.
 	  AKTR(IPK+32)= 2.*PRPS(JW1,1)*PRPS(JW1,9)	
 	  
	  AKTR(IPK+33)= 2.*PRSS(JB,3) + 2.*PRPS(JW1,3) + 2.*PRPS(JW2,3) 
	  AKTR(IPK+34)= 0.
	  AKTR(IPK+35)= 0.
	  AKTR(IPK+36)= 0.
	  AKTR(IPK+37)=-2.*PRPS(JW2,3)
	  AKTR(IPK+38)= 0.
	  AKTR(IPK+39)= 0.
	  AKTR(IPK+40)=-2.*PRPS(JW1,3)
	  
        AKTR(IPK+41)= 2.*PRPS(JW1,3)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,3)*PRPS(JW2,9)**2		
        AKTR(IPK+42)= 0.
        AKTR(IPK+43)= 0.
        AKTR(IPK+44)= 0.
	  AKTR(IPK+45)= 0.

	  AKTR(IPK+46)= 2.*PRPS(JW2,3)*PRPS(JW2,9)
        AKTR(IPK+47)= 0.
	  AKTR(IPK+48)= 0.
	  AKTR(IPK+49)=-2.*PRPS(JW1,3)*PRPS(JW2,9)	

	  AKTR(IPK+50)=	2.*PRSS(JB,1) + ALSU(IAC,1)	
	  AKTR(IPK+51)= 0.
	  AKTR(IPK+52)= 0.
	  AKTR(IPK+53)= 0.

	  AKTR(IPK+54)=-2.*PRSS(JB,1)*PRSS(JB,9)
	  AKTR(IPK+55)=-2.*PRSS(JB,1)	
	  		
	  AKTR(IPK+56)= 2.*PRSS(JB,1)*CADI(IA,6)**2 + 
	1				2.*PRSS(JB,3)*PRSS(JB,8)**2 + 
	2				   ALSU(IAC,1)*CADI(IA,7)**2

	  AKTR(IPK+57)=-2.*PRSS(JB,1)*CADI(IA,6) - ALSU(IAC,1)*CADI(IA,7) 
	  AKTR(IPK+58)= 0.
	  AKTR(IPK+59)= 0.
	  AKTR(IPK+60)= 0.

	  AKTR(IPK+61)= 2.*PRSS(JB,1)*CADI(IA,6)* PRSS(JB,9)
	1			   -2.*PRSS(JB,3)*PRSS(JB,8)**2

	  AKTR(IPK+62)= 2.*PRSS(JB,1)*CADI(IA,6)
	  		
	  AKTR(IPK+63)= 2.*PRSS(JB,1)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,1)*CADI(IA,2)**2 +
	2				4.*ALSU(IAC,6)*ADDI(IAC,1)**2

	  AKTR(IPK+64)= 0.0 !-2.*PRSS(JB,1)*CADI(IA,1)*CADI(IA,6)
C	1			      -ALSU(IAC,1)*CADI(IA,2)*CADI(IA,7)

	  AKTR(IPK+65)= 0.0 !2.*PRSS(JB,1)*CADI(IA,1) + ALSU(IAC,1)*CADI(IA,2)	

cc
	  AKTR(IPK+66)= 0.0
	  AKTR(IPK+67)= 0.0
	  AKTR(IPK+68)= 0.0
	  AKTR(IPK+69)= 2.*PRSS(JB,1)*CADI(IA,1)*PRSS(JB,9)
	  AKTR(IPK+70)= 2.*PRSS(JB,1)*CADI(IA,1)
cc

	  AKTR(IPK+71)=	2.*PRSS(JB,3) + ALSU(IAC,2)	

	  AKTR(IPK+72)= 0.
	  AKTR(IPK+73)= 0.
	  AKTR(IPK+74)= 0.
	  AKTR(IPK+75)= 0.
	  AKTR(IPK+76)=-2.*PRSS(JB,3)		

	  AKTR(IPK+77)= 2.*PRSS(JB,3)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,2)*CADI(IA,2)**2 +
	2				2.*ALSU(IAC,6)*ADDI(IAC,3)**2 +
	3				2.*ALSU(IAC,6)*ADDI(IAC,5)**2		

	  AKTR(IPK+78)= 0.0 !-2.*PRSS(JB,3)*CADI(IA,1) - ALSU(IAC,2)*CADI(IA,2)
	   
	  AKTR(IPK+79)= 0.
	  AKTR(IPK+80)= 0.
	  AKTR(IPK+81)= 0.
	  AKTR(IPK+82)= 0.
	  AKTR(IPK+83)= 2.*PRSS(JB,3)*CADI(IA,1)

C	Damping C11
	  ACTR(IPK+1 )= 2.*PRPS(JW1,2) + 2.*WRCO(JW1,4)	!--> 1st wheel
        ACTR(IPK+2 )= 2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
	1				2.*WRCO(JW1,6)*PRWH(JW1,6)**2 +
	2				2.*WRCO(JW1,4)*PRWH(JW1,5)**2

	  ACTR(IPK+3 )=-2.*WRCO(JW1,4)*PRWH(JW1,5)
	  ACTR(IPK+4 )= 2.*PRPS(JW1,4) + 2.*WRCO(JW1,6)

	  ACTR(IPK+5 )= 2.*PRPS(JW2,2) + 2.*WRCO(JW2,4)	!--> 2nd wheel
        ACTR(IPK+6 )= 2.*PRPS(JW2,4)*PRPS(JW2,7)**2 +
	1				2.*WRCO(JW2,6)*PRWH(JW2,6)**2 +
	2				2.*WRCO(JW2,4)*PRWH(JW2,5)**2

	  ACTR(IPK+7 )=-2.*WRCO(JW2,4)*PRWH(JW2,5)
	  ACTR(IPK+8 )= 2.*PRPS(JW2,4) + 2.*WRCO(JW2,6)

	  ACTR(IPK+9 )= 2.*PRSS(JB,2) + 2.*PRPS(JW1,2) + 2.*PRPS(JW2,2) !-->bogie
	  ACTR(IPK+10)= 0.
	  ACTR(IPK+11)= 0.

 	  ACTR(IPK+12)=-2.*PRPS(JW2,2)

	  ACTR(IPK+13)= 0.
	  ACTR(IPK+14)= 0.
 	  ACTR(IPK+15)=-2.*PRPS(JW1,2)		

        ACTR(IPK+16)= 2.*PRSS(JB,2)*PRSS(JB,9)**2 +
	1			    2.*PRSS(JB,4)*PRSS(JB,8)**2 +
	2				2.*PRPS(JW1,2)*PRPS(JW1,8)**2 +
     3				2.*PRPS(JW2,2)*PRPS(JW2,8)**2 +
     4				2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
     5				2.*PRPS(JW2,4)*PRPS(JW2,7)**2   				 

        ACTR(IPK+17)= 2.*PRSS(JB,2)*PRSS(JB,9) 
	1			   -2.*PRPS(JW1,2)*PRPS(JW1,8)
	2			   -2.*PRPS(JW2,2)*PRPS(JW2,8)		

        ACTR(IPK+18)= 0.
	  ACTR(IPK+19)=-2.*PRPS(JW2,4)*PRPS(JW2,7)**2
	  ACTR(IPK+20)= 2.*PRPS(JW2,2)*PRPS(JW2,8)

        ACTR(IPK+21)= 0.
	  ACTR(IPK+22)=-2.*PRPS(JW1,4)*PRPS(JW1,7)**2
	  ACTR(IPK+23)= 2.*PRPS(JW1,2)*PRPS(JW1,8)
	  		
        ACTR(IPK+24)= 2.*PRPS(JW1,2)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,2)*PRPS(JW2,9)**2	

	  ACTR(IPK+25)= 0.
	  ACTR(IPK+26)= 0.
	  ACTR(IPK+27)= 0.
	  ACTR(IPK+28)= 0.

 	  ACTR(IPK+29)=-2.*PRPS(JW2,2)*PRPS(JW2,9)
	  ACTR(IPK+30)= 0.
	  ACTR(IPK+31)= 0.

 	  ACTR(IPK+32)= 2.*PRPS(JW1,2)*PRPS(JW1,9)	
	  ACTR(IPK+33)= 2.*PRSS(JB,4) + 2.*PRPS(JW1,4) + 2.*PRPS(JW2,4) 
	  ACTR(IPK+34)= 0.
	  ACTR(IPK+35)= 0.
	  ACTR(IPK+36)= 0.
	  ACTR(IPK+37)=-2.*PRPS(JW2,4)
	  ACTR(IPK+38)= 0.
	  ACTR(IPK+39)= 0.
	  ACTR(IPK+40)=-2.*PRPS(JW1,4)
        ACTR(IPK+41)= 2.*PRPS(JW1,4)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,4)*PRPS(JW2,9)**2		
        ACTR(IPK+42)= 0.
        ACTR(IPK+43)= 0.
        ACTR(IPK+44)= 0.
	  ACTR(IPK+45)= 0.

	  ACTR(IPK+46)= 2.*PRPS(JW2,4)*PRPS(JW2,9)
        ACTR(IPK+47)= 0.
	  ACTR(IPK+48)= 0.
	  ACTR(IPK+49)=-2.*PRPS(JW1,4)*PRPS(JW2,9)	

	  ACTR(IPK+50)=	2.*PRSS(JB,2) + ALSU(IAC,4)	
	  ACTR(IPK+51)= 0.
	  ACTR(IPK+52)= 0.
	  ACTR(IPK+53)= 0.

	  ACTR(IPK+54)=-2.*PRSS(JB,2)*PRSS(JB,9)
	  ACTR(IPK+55)=-2.*PRSS(JB,2)	
	  		
	  ACTR(IPK+56)= 2.*PRSS(JB,2)*CADI(IA,6)**2 + 
	1				2.*PRSS(JB,4)*PRSS(JB,8)**2 + 
	2				   ALSU(IAC,4)*CADI(IA,7)**2

	  ACTR(IPK+57)=-2.*PRSS(JB,2)*CADI(IA,6) - ALSU(IAC,4)*CADI(IA,7) 
	  ACTR(IPK+58)= 0.
	  ACTR(IPK+59)= 0.
	  ACTR(IPK+60)= 0.

	  ACTR(IPK+61)= 2.*PRSS(JB,2)*CADI(IA,6)* PRSS(JB,9)
	1			   -2.*PRSS(JB,4)*PRSS(JB,8)**2

	  ACTR(IPK+62)= 2.*PRSS(JB,2)*CADI(IA,6)		
	  ACTR(IPK+63)= 2.*PRSS(JB,2)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,4)*CADI(IA,2)**2 +
	2				4.*ALSU(IAC,3)*ADDI(IAC,1)**2

	  ACTR(IPK+64)= 0.0 !-2.*PRSS(JB,2)*CADI(IA,1)*CADI(IA,6)
	1			      -ALSU(IAC,4)*CADI(IA,2)*CADI(IA,7)

	  ACTR(IPK+65)= 0.0 !2.*PRSS(JB,2)*CADI(IA,1) + ALSU(IAC,4)*CADI(IA,2)	

cc
	  ACTR(IPK+66)= 0.
	  ACTR(IPK+67)= 0.
	  ACTR(IPK+68)= 0.
	  ACTR(IPK+69)= 2.*PRSS(JB,2)*CADI(IA,1)*PRSS(JB,9)
	  ACTR(IPK+70)= 2.*PRSS(JB,2)*CADI(IA,1)
cc

	  
	  ACTR(IPK+71)=	2.*PRSS(JB,4) + ALSU(IAC,5)	

	  ACTR(IPK+72)= 0.
	  ACTR(IPK+73)= 0.
	  ACTR(IPK+74)= 0.
	  ACTR(IPK+75)= 0.
	  ACTR(IPK+76)=-2.*PRSS(JB,4)		

	  ACTR(IPK+77)= 2.*PRSS(JB,4)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,5)*CADI(IA,2)**2 +
	2				2.*ALSU(IAC,3)*ADDI(IAC,3)**2 +
	3				2.*ALSU(IAC,3)*ADDI(IAC,5)**2		

	  ACTR(IPK+78)= 0.0 !-2.*PRSS(JB,4)*CADI(IA,1) - ALSU(IAC,5)*CADI(IA,2)
	   
	  ACTR(IPK+79)= 0.
	  ACTR(IPK+80)= 0.
	  ACTR(IPK+81)= 0.
	  ACTR(IPK+82)= 0.
	  ACTR(IPK+83)= 2.*PRSS(JB,4)*CADI(IA,1)
	  
	  IPK=IPK+83

C	---------------------------------------------------------------------------
C	2.The intermidate passenger cars (1 <n <N) and the last transition car (n=N):	
C	 Form the tri-diagonal matrices by looping over each pair (bogie & car body)
C		- The last transition car (n = N), several terms of matrices must be  
C		  formulated seperately since the chacteristics of the transition car
C	---------------------------------------------------------------------------
	DO 30 IA = NTRAIN(1)+2,NTRAIN(1)+NTRAIN(2)
	  JB = NTRAIN(1)+IA		
	  JW1=2*JB-1            
	  JW2=2*JB              
	  IAC= IA-NTRAIN(1)		
C
C	  Stiffness Kn,n & Kn-1,n

	  AKTR(IPK+1 )= 2.*PRPS(JW1,1) + 2.*WRCO(JW1,3)	!--> 1st wheel
        AKTR(IPK+2 )= 2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
	1				2.*WRCO(JW1,5)*PRWH(JW1,6)**2 +
	2				2.*WRCO(JW1,3)*PRWH(JW1,5)**2
	  AKTR(IPK+3 )=-2.*WRCO(JW1,3)*PRWH(JW1,5)
	  AKTR(IPK+4 )= 2.*PRPS(JW1,3) + 2.*WRCO(JW1,5)

	  AKTR(IPK+5 )= 2.*PRPS(JW2,1) + 2.*WRCO(JW2,3)	!--> 2nd wheel
        AKTR(IPK+6 )= 2.*PRPS(JW2,3)*PRPS(JW2,7)**2 +
	1				2.*WRCO(JW2,5)*PRWH(JW2,6)**2 +
	2				2.*WRCO(JW2,3)*PRWH(JW2,5)**2
	  AKTR(IPK+7 )=-2.*WRCO(JW2,3)*PRWH(JW2,5)
	  AKTR(IPK+8 )= 2.*PRPS(JW2,3) + 2.*WRCO(JW2,5)

	  AKTR(IPK+9 )= 2.*PRSS(JB,1) + 2.*PRPS(JW1,1) + 2.*PRPS(JW2,1) !-->bogie
	  AKTR(IPK+10)= 0.
	  AKTR(IPK+11)= 0.

 	  AKTR(IPK+12)=-2.*PRPS(JW2,1)

	  AKTR(IPK+13)= 0.
	  AKTR(IPK+14)= 0.
 	  AKTR(IPK+15)=-2.*PRPS(JW1,1)
 	
CC      Taken from 1st transition bogie  dof 51-55  
 	  AKTR(IPK+16)= 0.
 	  AKTR(IPK+17)= 0.
 	  AKTR(IPK+18)= 2.*PRSS(JB,1)*CADI(IA,1)
 	  AKTR(IPK+19)= 2.*PRSS(JB,1)*CADI(IA,6)
 	  AKTR(IPK+20)=-2.*PRSS(JB,1)
CC 	  
 	  		

        AKTR(IPK+21)= 2.*PRSS(JB,1)*PRSS(JB,9)**2 +
	1			    2.*PRSS(JB,3)*PRSS(JB,8)**2 +
	2				2.*PRPS(JW1,1)*PRPS(JW1,8)**2 +
     3				2.*PRPS(JW2,1)*PRPS(JW2,8)**2 +
     4				2.*PRPS(JW1,3)*PRPS(JW1,7)**2 +
     5				2.*PRPS(JW2,3)*PRPS(JW2,7)**2   				 

        AKTR(IPK+22)= 2.*PRSS(JB,1)*PRSS(JB,9) 
	1			   -2.*PRPS(JW1,1)*PRPS(JW1,8)
	2			   -2.*PRPS(JW2,1)*PRPS(JW2,8)		

        AKTR(IPK+23)= 0.
	  AKTR(IPK+24)=-2.*PRPS(JW2,3)*PRPS(JW2,7)**2
	  AKTR(IPK+25)= 2.*PRPS(JW2,1)*PRPS(JW2,8)

        AKTR(IPK+26)= 0.
	  AKTR(IPK+27)=-2.*PRPS(JW1,3)*PRPS(JW1,7)**2
	  AKTR(IPK+28)= 2.*PRPS(JW1,1)*PRPS(JW1,8)

CC      Taken from 1st transition bogie  dof 58-62
 	  AKTR(IPK+29)= 0.
 	  AKTR(IPK+30)= 0.
 	  AKTR(IPK+31)= 2.*PRSS(JB,1)*CADI(IA,1)*PRSS(JB,9)
 	  AKTR(IPK+32)= 2.*PRSS(JB,1)*CADI(IA,6)* PRSS(JB,9)
	1			   -2.*PRSS(JB,3)*PRSS(JB,8)**2
 	  AKTR(IPK+33)= -2.*PRSS(JB,1)*PRSS(JB,9)	  	
CC 	  	
	  		
	  		
	  		
        AKTR(IPK+34)= 2.*PRPS(JW1,1)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,1)*PRPS(JW2,9)**2	

	  AKTR(IPK+35)= 0.
	  AKTR(IPK+36)= 0.
	  AKTR(IPK+37)= 0.
	  AKTR(IPK+38)= 0.

 	  AKTR(IPK+39)=-2.*PRPS(JW2,1)*PRPS(JW2,9)
	  AKTR(IPK+40)= 0.
	  AKTR(IPK+41)= 0.

 	  AKTR(IPK+42)= 2.*PRPS(JW1,1)*PRPS(JW1,9)	


CC 
 	  AKTR(IPK+43)= 0.
	  AKTR(IPK+44)= 0.
	  AKTR(IPK+45)= 0.
	  AKTR(IPK+46)= 0.
	  AKTR(IPK+47)= 0.
CC

 	  
	  AKTR(IPK+48)= 2.*PRSS(JB,3) + 2.*PRPS(JW1,3) + 2.*PRPS(JW2,3) 
	  AKTR(IPK+49)= 0.
	  AKTR(IPK+50)= 0.
	  AKTR(IPK+51)= 0.
	  AKTR(IPK+52)=-2.*PRPS(JW2,3)
	  AKTR(IPK+53)= 0.
	  AKTR(IPK+54)= 0.
	  AKTR(IPK+55)=-2.*PRPS(JW1,3)

CC      Taken from 1st transition bogie  dof 70-71  
	  AKTR(IPK+56)= 0.
	  AKTR(IPK+57)= 0.	
	  AKTR(IPK+58)=-2.*PRSS(JB,3)*CADI(IA,1)
	  AKTR(IPK+59)=-2.*PRSS(JB,3)	
	  AKTR(IPK+60)= 0.
	  AKTR(IPK+61)= 0.	
	  AKTR(IPK+62)= 0.
CC	  	  
	  
	  
        AKTR(IPK+63)= 2.*PRPS(JW1,3)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,3)*PRPS(JW2,9)**2		
        AKTR(IPK+64)= 0.
        AKTR(IPK+65)= 0.
        AKTR(IPK+66)= 0.
	  AKTR(IPK+67)= 0.

	  AKTR(IPK+68)= 2.*PRPS(JW2,3)*PRPS(JW2,9)
        AKTR(IPK+69)= 0.
	  AKTR(IPK+70)= 0.
	  AKTR(IPK+71)=-2.*PRPS(JW1,3)*PRPS(JW2,9)	

CC      Taken from 1st transition bogie  dof 77-78	  
	  AKTR(IPK+72)= 0.
	  AKTR(IPK+73)= 0.
	  AKTR(IPK+74)= 0.
	  AKTR(IPK+75)= 0.
	  AKTR(IPK+76)= 0.
	  AKTR(IPK+77)= 0.
	  AKTR(IPK+78)= 0.
CC	  
	  

	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars
		AKTR(IPK+79)= ALSU(IAC-1,1) + 2.*PRSS(JB,1) + ALSU(IAC,1)
	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car
		AKTR(IPK+79)= ALSU(IAC-1,1) + 2.*PRSS(JB,1) + 2.*PRSS(JB+1,1)
	  AKTR(IPK+79)=	2.*PRSS(JB,1) + ALSU(IAC,1)
	  ENDIF

	  AKTR(IPK+80)= 0.
	  AKTR(IPK+81)= 0.
	  AKTR(IPK+82)= 0.

	  AKTR(IPK+83)=-2.*PRSS(JB,1)*PRSS(JB,9)
	  AKTR(IPK+84)=-2.*PRSS(JB,1)
	  
  
	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 				
		AKTR(IPK+85)=   ALSU(IAC-1,1)*CADI(IA,6)**2 +
	1				 2.*PRSS(JB,1)*CADI(IA,6)**2 + 
	2				 2.*PRSS(JB,3)*PRSS(JB,8)**2 + 
	3					ALSU(IAC,1)*CADI(IA,7)**2

	    AKTR(IPK+86)=-ALSU(IAC-1,1)*CADI(IA,6)
     1				 - 2.*PRSS(JB,1)*CADI(IA,6)
	2			     -ALSU(IAC,1)*CADI(IA,7) 

	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car
		AKTR(IPK+85)=   ALSU(IAC-1,1)*CADI(IA,6)**2 +
	1				 2.*PRSS(JB,1)*CADI(IA,6)**2 + 
	2				 2.*PRSS(JB,3)*PRSS(JB,8)**2 + 
	3				 2.*PRSS(JB+1,1)*CADI(IA,7)**2 + 
	4				 2.*PRSS(JB+1,3)*PRSS(JB+1,8)**2

	    AKTR(IPK+86)=-ALSU(IAC-1,1)*CADI(IA,6) 
	1				 - 2.*PRSS(JB,1)*CADI(IA,6)
	2			     - 2.*PRSS(JB+1,1)*CADI(IA,7)
	  AKTR(IPK+85)= 2.*PRSS(JB,1)*CADI(IA,6)**2 + 
	1				2.*PRSS(JB,3)*PRSS(JB,8)**2 + 
	2				   ALSU(IAC,1)*CADI(IA,7)**2

	  AKTR(IPK+86)=-2.*PRSS(JB,1)*CADI(IA,6) - ALSU(IAC,1)*CADI(IA,7) 
	  ENDIF

	  AKTR(IPK+87)= 0.
	  AKTR(IPK+88)= 0.
	  AKTR(IPK+89)= 0.

	  AKTR(IPK+90)= 2.*PRSS(JB,1)*CADI(IA,6)* PRSS(JB,9)
	1			   -2.*PRSS(JB,3)*PRSS(JB,8)**2
	  AKTR(IPK+91)= 2.*PRSS(JB,1)*CADI(IA,6)			

	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 	
	    AKTR(IPK+92)=  ALSU(IAC-1,1)*CADI(IA,1)**2 +
     1				 4.*ALSU(IAC-1,6)*ADDI(IAC,1)**2 + 
     2			     2.*PRSS(JB,1)*CADI(IA,1)**2 + 
	3				    ALSU(IAC,1)*CADI(IA,2)**2 +
	4				 4.*ALSU(IAC,6)*ADDI(IAC,1)**2			

	    AKTR(IPK+93)=  0.0 !-ALSU(IAC-1,1)*CADI(IA,1)*CADI(IA,6) +
c	1			    -2.*PRSS(JB,1)*CADI(IA,1)*CADI(IA,6) +
c	2			        ALSU(IAC,1)*CADI(IA,2)*CADI(IA,7)

	    AKTR(IPK+94)=  0.0 !ALSU(IAC-1,1)*CADI(IA,1) + 
c	1				2.*PRSS(JB,1)*CADI(IA,1)
c	2			     - ALSU(IAC,1)*CADI(IA,2)

	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car   
	    AKTR(IPK+92)=   ALSU(IAC-1,1)*CADI(IA,1)**2 +
     1				 4.*ALSU(IAC-1,6)*ADDI(IAC,1)**2 + 
     2			     2.*PRSS(JB,1)*CADI(IA,1)**2 + 
     2			     2.*PRSS(JB+1,1)*CADI(IA,2)**2		

	    AKTR(IPK+93)=   0.0 !-ALSU(IAC-1,1)*CADI(IA,1)*CADI(IA,6) +
c	1			     -2.*PRSS(JB,1)*CADI(IA,1)*CADI(IA,6) +
c	2			      2.*PRSS(JB+1,1)*CADI(IA,2)*CADI(IA,7)

	    AKTR(IPK+94)=    0.0 !ALSU(IAC-1,1)*CADI(IA,1) + 
c	1				  2.*PRSS(JB,1)*CADI(IA,1)
c	2			     -2.*PRSS(JB+1,1)*CADI(IA,2) 

	  		
	  AKTR(IPK+92)= 2.*PRSS(JB,1)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,1)*CADI(IA,2)**2 +
	2				4.*ALSU(IAC,6)*ADDI(IAC,1)**2

	  AKTR(IPK+93)= 0.0 !-2.*PRSS(JB,1)*CADI(IA,1)*CADI(IA,6)
C	1			      -ALSU(IAC,1)*CADI(IA,2)*CADI(IA,7)

	  AKTR(IPK+94)= 0.0 !2.*PRSS(JB,1)*CADI(IA,1) + ALSU(IAC,1)*CADI(IA,2)	
	  ENDIF  			 		

	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 
	    AKTR(IPK+95)= ALSU(IAC-1,2) + 2.*PRSS(JB,3) + ALSU(IAC,2)
	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car 
	    AKTR(IPK+95)= ALSU(IAC-1,2) + 2.*PRSS(JB,3) + 2.*PRSS(JB+1,3)	 
	  AKTR(IPK+95)=	2.*PRSS(JB,3) + ALSU(IAC,2) 
	  ENDIF	  		
			
	  AKTR(IPK+96)= 0.
	  AKTR(IPK+97)= 0.
	  AKTR(IPK+98)= 0.
	  AKTR(IPK+99)= 0.

	  AKTR(IPK+100)=-2.*PRSS(JB,3)

	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 
	    AKTR(IPK+101)=   ALSU(IAC-1,2)*CADI(IA,1)**2 +
	1				  2.*ALSU(IAC-1,6)*ADDI(IAC,2)**2 +
	2				  2.*ALSU(IAC-1,6)*ADDI(IAC,4)**2 +
	3				  2.*PRSS(JB,3)*CADI(IA,1)**2 + 
	4				     ALSU(IAC,2)*CADI(IA,2)**2 +
	5				  2.*ALSU(IAC,6)*ADDI(IAC,3)**2 +
	6				  2.*ALSU(IAC,6)*ADDI(IAC,5)**2

	    AKTR(IPK+102)= 0.0 !-ALSU(IAC-1,2)*CADI(IA,1)
c	1				  -2.*PRSS(JB,3)*CADI(IA,1)
c	2				  -ALSU(IAC,2)*CADI(IA,2)

        ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car 
	    AKTR(IPK+101)=   ALSU(IAC-1,2)*CADI(IA,1)**2 +
	1				  2.*ALSU(IAC-1,6)*ADDI(IAC,2)**2 +
	2				  2.*ALSU(IAC-1,6)*ADDI(IAC,4)**2 +
	3				  2.*PRSS(JB,3)*CADI(IA,1)**2 + 
	4				  2.*PRSS(JB+1,3)*CADI(IA,2)**2

	    AKTR(IPK+102)= 0.0 !-ALSU(IAC-1,2)*CADI(IA,1)
c	1				  -2.*PRSS(JB,3)*CADI(IA,1)
c	2				  +2.*PRSS(JB+1,3)*CADI(IA,2)


	  AKTR(IPK+101)= 2.*PRSS(JB,3)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,2)*CADI(IA,2)**2 +
	2				2.*ALSU(IAC,6)*ADDI(IAC,3)**2 +
	3				2.*ALSU(IAC,6)*ADDI(IAC,5)**2
	
	  AKTR(IPK+102)= 0.0 !-2.*PRSS(JB,3)*CADI(IA,1) - ALSU(IAC,2)*CADI(IA,2)	
	  ENDIF

	  AKTR(IPK+103)= 0.
	  AKTR(IPK+104)= 0.
	  AKTR(IPK+105)= 0.
	  AKTR(IPK+106)= 0.
	  AKTR(IPK+107)= 2.*PRSS(JB,3)*CADI(IA,1)  		

C	  Damping Cn,n & Cn-1,n
	  ACTR(IPK+1 )= 2.*PRPS(JW1,2) + 2.*WRCO(JW1,4)	!--> 1st wheel
        ACTR(IPK+2 )= 2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
	1				2.*WRCO(JW1,6)*PRWH(JW1,6)**2 +
	2				2.*WRCO(JW1,4)*PRWH(JW1,5)**2

	  ACTR(IPK+3 )=-2.*WRCO(JW1,4)*PRWH(JW1,5)
	  ACTR(IPK+4 )= 2.*PRPS(JW1,4) + 2.*WRCO(JW1,6)		!stop here

	  ACTR(IPK+5 )= 2.*PRPS(JW2,2) + 2.*WRCO(JW2,4)	!--> 2nd wheel

        ACTR(IPK+6 )= 2.*PRPS(JW2,4)*PRPS(JW2,7)**2 +
	1				2.*WRCO(JW2,6)*PRWH(JW2,6)**2 +
	2				2.*WRCO(JW2,4)*PRWH(JW2,5)**2

	  ACTR(IPK+7 )=-2.*WRCO(JW2,4)*PRWH(JW2,5)
	  ACTR(IPK+8 )= 2.*PRPS(JW2,4) + 2.*WRCO(JW2,6)

	  ACTR(IPK+9 )= 2.*PRSS(JB,2) + 2.*PRPS(JW1,2) + 2.*PRPS(JW2,2) !-->bogie
	  ACTR(IPK+10)= 0.
	  ACTR(IPK+11)= 0.

 	  ACTR(IPK+12)=-2.*PRPS(JW2,2)

	  ACTR(IPK+13)= 0.
	  ACTR(IPK+14)= 0.
 	  ACTR(IPK+15)=-2.*PRPS(JW1,2)		

CC      Taken from 1st transition bogie  dof 51-55  
 	  ACTR(IPK+16)= 0.
 	  ACTR(IPK+17)= 0.
 	  ACTR(IPK+18)= 0.0 !2.*PRSS(JB,2)*CADI(IA,1)
 	  ACTR(IPK+19)= 0.0 !2.*PRSS(JB,2)*CADI(IA,6)
 	  ACTR(IPK+20)= 0.0 !-2.*PRSS(JB,2)
CC  	  

        ACTR(IPK+21)= 2.*PRSS(JB,2)*PRSS(JB,9)**2 +
	1			    2.*PRSS(JB,4)*PRSS(JB,8)**2 +
	2				2.*PRPS(JW1,2)*PRPS(JW1,8)**2 +
     3				2.*PRPS(JW2,2)*PRPS(JW2,8)**2 +
     4				2.*PRPS(JW1,4)*PRPS(JW1,7)**2 +
     5				2.*PRPS(JW2,4)*PRPS(JW2,7)**2   				 

        ACTR(IPK+22)= 2.*PRSS(JB,2)*PRSS(JB,9) 
	1			   -2.*PRPS(JW1,2)*PRPS(JW1,8)
	2			   -2.*PRPS(JW2,2)*PRPS(JW2,8)		

        ACTR(IPK+23)= 0.
	  ACTR(IPK+24)=-2.*PRPS(JW2,4)*PRPS(JW2,7)**2
	  ACTR(IPK+25)= 2.*PRPS(JW2,2)*PRPS(JW2,8)

        ACTR(IPK+26)= 0.
	  ACTR(IPK+27)=-2.*PRPS(JW1,4)*PRPS(JW1,7)**2
	  ACTR(IPK+28)= 2.*PRPS(JW1,2)*PRPS(JW1,8)

CC      Taken from 1st transition bogie  dof 58-62
 	  ACTR(IPK+29)= 0.
 	  ACTR(IPK+30)= 0.
 	  ACTR(IPK+31)= 2.*PRSS(JB,2)*CADI(IA,1)* PRSS(JB,9)
 	  ACTR(IPK+32)= 2.*PRSS(JB,2)*CADI(IA,6)* PRSS(JB,9)
	1			   -2.*PRSS(JB,4)*PRSS(JB,8)**2
 	  ACTR(IPK+33)=-2.*PRSS(JB,2)*PRSS(JB,9)	  	
CC 	 


	  		
        ACTR(IPK+34)= 2.*PRPS(JW1,2)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,2)*PRPS(JW2,9)**2	

	  ACTR(IPK+35)= 0.
	  ACTR(IPK+36)= 0.
	  ACTR(IPK+37)= 0.
	  ACTR(IPK+38)= 0.

 	  ACTR(IPK+39)=-2.*PRPS(JW2,2)*PRPS(JW2,9)
	  ACTR(IPK+40)= 0.
	  ACTR(IPK+41)= 0.

 	  ACTR(IPK+42)= 2.*PRPS(JW1,2)*PRPS(JW1,9)
 	  

CC      
 	  ACTR(IPK+43)= 0.
 	  ACTR(IPK+44)= 0.
 	  ACTR(IPK+45)= 0.
 	  ACTR(IPK+46)= 0.
 	  ACTR(IPK+47)= 0.	  	
CC 	 

	
	  ACTR(IPK+48)= 2.*PRSS(JB,4) + 2.*PRPS(JW1,4) + 2.*PRPS(JW2,4) 
	  ACTR(IPK+49)= 0.
	  ACTR(IPK+50)= 0.
	  ACTR(IPK+51)= 0.
	  ACTR(IPK+52)=-2.*PRPS(JW2,4)
	  ACTR(IPK+53)= 0.
	  ACTR(IPK+54)= 0.
	  ACTR(IPK+55)=-2.*PRPS(JW1,4)
CC      Taken from 1st transition bogie  dof 70-71  
	  ACTR(IPK+56)= 0.
	  ACTR(IPK+57)= 0.	
	  ACTR(IPK+58)= -2.*PRSS(JB,4)*CADI(IA,1)
	  ACTR(IPK+59)= -2.*PRSS(JB,4)
 	  ACTR(IPK+60)= 0.
 	  ACTR(IPK+61)= 0.
 	  ACTR(IPK+62)= 0.		
CC
	  
	  
        ACTR(IPK+63)= 2.*PRPS(JW1,4)*PRPS(JW1,9)**2 +
     1				2.*PRPS(JW2,4)*PRPS(JW2,9)**2		
        ACTR(IPK+64)= 0.
        ACTR(IPK+65)= 0.
        ACTR(IPK+66)= 0.
	  ACTR(IPK+67)= 0.

	  ACTR(IPK+68)= 2.*PRPS(JW2,4)*PRPS(JW2,9)
        ACTR(IPK+69)= 0.
	  ACTR(IPK+70)= 0.
	  ACTR(IPK+71)=-2.*PRPS(JW1,4)*PRPS(JW2,9)	
CC      Taken from 1st transition bogie  dof 77-78	  
	  ACTR(IPK+72)= 0.
	  ACTR(IPK+73)= 0.
	  ACTR(IPK+74)= 0.
	  ACTR(IPK+75)= 0.
 	  ACTR(IPK+76)= 0.
 	  ACTR(IPK+77)= 0.
 	  ACTR(IPK+78)= 0.
CC	 


	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars
		ACTR(IPK+79)= ALSU(IAC-1,4)	+ 2.*PRSS(JB,2) + ALSU(IAC,4)	
	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car
		ACTR(IPK+79)= ALSU(IAC-1,4) + 2.*PRSS(JB,2) + 2.*PRSS(JB+1,2)
	  ACTR(IPK+79)=	2.*PRSS(JB,2) + ALSU(IAC,4)	
	  ENDIF

	  ACTR(IPK+80)= 0.
	  ACTR(IPK+81)= 0.
	  ACTR(IPK+82)= 0.

	  ACTR(IPK+83)=-2.*PRSS(JB,2)*PRSS(JB,9)
	  ACTR(IPK+84)=-2.*PRSS(JB,2)

	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 				
	    ACTR(IPK+85)=    ALSU(IAC-1,4)*CADI(IA,6)**2 +
	1			      2.*PRSS(JB,2)*CADI(IA,6)**2 + 
	2				  2.*PRSS(JB,4)*PRSS(JB,8)**2 + 
	3				     ALSU(IAC,4)*CADI(IA,7)**2

	    ACTR(IPK+86)=    -ALSU(IAC-1,4)*CADI(IA,6)
	1				 - 2.*PRSS(JB,2)*CADI(IA,6)
	2			         -ALSU(IAC,4)*CADI(IA,7)  

	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car
		ACTR(IPK+85)=   ALSU(IAC-1,4)*CADI(IA,6)**2 +
	1				 2.*PRSS(JB,2)*CADI(IA,6)**2 + 
	2				 2.*PRSS(JB,4)*PRSS(JB,8)**2 + 
	3				 2.*PRSS(JB+1,2)*CADI(IA,7)**2 + 
	4				 2.*PRSS(JB+1,4)*PRSS(JB+1,8)**2

	    AKTR(IPK+86)=    -ALSU(IAC-1,4)*CADI(IA,6)
	1				 - 2.*PRSS(JB,2)*CADI(IA,6)
	2			     - 2.*PRSS(JB+1,2)*CADI(IA,7)
	  ACTR(IPK+85)= 2.*PRSS(JB,2)*CADI(IA,6)**2 + 
	1				2.*PRSS(JB,4)*PRSS(JB,8)**2 + 
	2				   ALSU(IAC,4)*CADI(IA,7)**2

	  ACTR(IPK+86)=-2.*PRSS(JB,2)*CADI(IA,6) - ALSU(IAC,4)*CADI(IA,7) 
	  ENDIF 	  				

	  ACTR(IPK+87)= 0.
	  ACTR(IPK+88)= 0.
	  ACTR(IPK+89)= 0.

	  ACTR(IPK+90)= 2.*PRSS(JB,2)*CADI(IA,6)* PRSS(JB,9)
	1			   -2.*PRSS(JB,4)*PRSS(JB,8)**2

	  ACTR(IPK+91)= 2.*PRSS(JB,2)*CADI(IA,6)	
	  
	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 	
	    ACTR(IPK+92)=    ALSU(IAC-1,2)*CADI(IA,1)**2 +
     1				  4.*ALSU(IAC-1,3)*ADDI(IAC,1)**2 
     2			      2.*PRSS(JB,2)*CADI(IA,1)**2 + 
	3				     ALSU(IAC,4)*CADI(IA,2)**2 +
	4				  4.*ALSU(IAC,3)*ADDI(IAC,1)**2			

	    ACTR(IPK+93)=   0.0 !-ALSU(IAC-1,4)*CADI(IA,1)*CADI(IA,6) +
c	1			     -2.*PRSS(JB,2)*CADI(IA,1)*CADI(IA,6) +
c	2			         ALSU(IAC,4)*CADI(IA,2)*CADI(IA,7)

	    ACTR(IPK+94)=   0.0 ! ALSU(IAC-1,4)*CADI(IA,1) + 
c	1				  2.*PRSS(JB,2)*CADI(IA,1)
c	2			       - ALSU(IAC,4)*CADI(IA,2)

	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car   
	    ACTR(IPK+92)=   ALSU(IAC-1,4)*CADI(IA,1)**2 +
     1				 4.*ALSU(IAC-1,3)*ADDI(IAC,1)**2 + 
     2			     2.*PRSS(JB,2)*CADI(IA,1)**2 + 
     2			     2.*PRSS(JB+1,2)*CADI(IA,2)**2		

	    ACTR(IPK+93)= 0.0 !  -ALSU(IAC-1,4)*CADI(IA,1)*CADI(IA,6) +
c	1			     -2.*PRSS(JB,2)*CADI(IA,1)*CADI(IA,6) +
c	2			      2.*PRSS(JB+1,2)*CADI(IA,2)*CADI(IA,7)

	    ACTR(IPK+94)=  0.0 !  ALSU(IAC-1,4)*CADI(IA,1) + 
c	1				  2.*PRSS(JB,2)*CADI(IA,1)
c	2			     -2.*PRSS(JB+1,2)*CADI(IA,2) 

	  ACTR(IPK+92)= 2.*PRSS(JB,2)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,4)*CADI(IA,2)**2 +
	2				4.*ALSU(IAC,3)*ADDI(IAC,1)**2

	  ACTR(IPK+93)= 0.0 !-2.*PRSS(JB,2)*CADI(IA,1)*CADI(IA,6)
c	1			      -ALSU(IAC,4)*CADI(IA,2)*CADI(IA,7)

	  ACTR(IPK+94)= 0.0 !2.*PRSS(JB,2)*CADI(IA,1) + ALSU(IAC,4)*CADI(IA,2)	
	  ENDIF  				  						

	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 
	    ACTR(IPK+95)= ALSU(IAC-1,5) + 2.*PRSS(JB,4) + ALSU(IAC,5)
	  ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car 
	    ACTR(IPK+95)= ALSU(IAC-1,5) + 2.*PRSS(JB,4) + 2.*PRSS(JB+1,4)	 
	  ACTR(IPK+95)=	2.*PRSS(JB,4) + ALSU(IAC,5)	 
	  ENDIF
	  	
	  ACTR(IPK+96)= 0.
	  ACTR(IPK+97)= 0.
	  ACTR(IPK+98)= 0.
	  ACTR(IPK+99)= 0.

	  ACTR(IPK+100)=-2.*PRSS(JB,4)			

	  IF(IAC.LT.NTRAIN(2)) THEN				!Intermidiate cars 
	    ACTR(IPK+101)=    ALSU(IAC-1,5)*CADI(IA,1)**2 +
	1				   2.*ALSU(IAC-1,3)*ADDI(IAC,2)**2 +
	2				   2.*ALSU(IAC-1,3)*ADDI(IAC,4)**2 +
	3				   2.*PRSS(JB,4)*CADI(IA,1)**2 + 
	4				      ALSU(IAC,5)*CADI(IA,2)**2 +
	5				   2.*ALSU(IAC,3)*ADDI(IAC,3)**2 +
	6				   2.*ALSU(IAC,3)*ADDI(IAC,5)**2

	    ACTR(IPK+102)=0.0 !-ALSU(IAC-1,5)*CADI(IA,1)
c	1				  -2.*PRSS(JB,4)*CADI(IA,1) 
c	2				  -ALSU(IAC,5)*CADI(IA,2)

        ELSE IF (IAC.EQ.NTRAIN(2)) THEN		!Last transition car 
	    ACTR(IPK+101)=   ALSU(IAC-1,5)*CADI(IA,1)**2 +
	1				  2.*ALSU(IAC-1,3)*ADDI(IAC,2)**2 +
	2				  2.*ALSU(IAC-1,3)*ADDI(IAC,4)**2 +
	3				  2.*PRSS(JB,4)*CADI(IA,1)**2 + 
	4				  2.*PRSS(JB+1,4)*CADI(IA,2)**2

	    ACTR(IPK+102)=0.0 !-ALSU(IAC-1,5)*CADI(IA,1)
c	1				  -2.*PRSS(JB,4)*CADI(IA,1)
c	2				  +2.*PRSS(JB+1,4)*CADI(IA,2)

	  ACTR(IPK+101)= 2.*PRSS(JB,4)*CADI(IA,1)**2 + 
	1				   ALSU(IAC,5)*CADI(IA,2)**2 +
	2				2.*ALSU(IAC,3)*ADDI(IAC,3)**2 +
	3				2.*ALSU(IAC,3)*ADDI(IAC,5)**2		

	  ACTR(IPK+102)=0.0 !-2.*PRSS(JB,4)*CADI(IA,1) - ALSU(IAC,5)*CADI(IA,2)
	  ENDIF

	  ACTR(IPK+103)= 0.
	  ACTR(IPK+104)= 0.
	  ACTR(IPK+105)= 0.
	  ACTR(IPK+106)= 0.

	  ACTR(IPK+107)= 2.*PRSS(JB,4)*CADI(IA,1)

	  IPK=IPK+107

   30 CONTINUE

	RETURN
	END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE TITLEBTI(SOLT,TRLENG,NOUTST,NOUTVH,KOST,STROHED,KOVH,TRNOHED) 
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
	CHARACTER*200 NAME
	CHARACTER*200 STROHED,TRNOHED
	CHARACTER*40 HED
	CHARACTER*4000 PHED
C	-----------------------------------------------------------------------
C	PROGRAM TO:
C		- PRINT OUT TITLE FOR BRIDGE/TRAIN IN DYNAMIC INTERACTION ANALYSIS
C	-----------------------------------------------------------------------
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
C	-----------------------------------------------------------------------
      DIMENSION KOST(2,NOUTST),STROHED(NOUTST),KOVH(NOUTVH),TRNOHED(NOUTVH)
C	-----------------------------------------------------------------------

      
c	Velocity in Km/hr
	VELKMH = VEL*3.6 

c	Total analysis length:
	TALENG = SOLT*VEL

C	Time running before 1st wheel enter the bridge
	BITIME = BIDIS/VEL
			
C	Time for train running on bridge (1st wheel enter to last wheel left)
	BRTIME = (TLENGTH + TRLENG)/VEL

c	Extra time after last wheel leaves the path (for free vibration)
	TEXTRA = EXTDIS/VEL	

	WRITE(*,200)	

C	PRINT TITLE FOR BRIDGE OUTPUT -- OUTJOB.OUT
	WRITE(100,500)
	WRITE(100,600)
	WRITE(100,1000) VELKMH,VEL,TALENG,SOLT,BIDIS,BITIME,TLENGTH,
     1				  TRLENG,BRTIME,EXTDIS,TEXTRA


      IF(NOUTST.GT.0) THEN
      
      WRITE(100,1200)

      DO IOUTST = 1,NOUTST
          NAME = STROHED(IOUTST)
	    NAML = LEN_TRIM(NAME)
    	
          NCHA = 40
          HED(1:NAML) = NAME(1:NAML)
          HED(NAML+1:NCHA) = ' '
          
          I1 = 1    + NCHA*(IOUTST-1)
          I2 = NCHA + NCHA*(IOUTST-1)
          PHED(I1:I2) = HED(1:NCHA)
      ENDDO
      WRITE(100,1300) PHED(1:NCHA*NOUTST)  
      
      
      NCHA = 40
      HED(1:4) = 'DIS.' 
      HED(5:11) = ' '
      HED(12:15) = 'VEL.' 
      HED(16:22) = ' '
      HED(23:26) = 'ACC.' 
      HED(27:40) = ' '
      DO IOUTST = 1,NOUTST
          I1 = 1    + NCHA*(IOUTST-1)
          I2 = NCHA + NCHA*(IOUTST-1)
          PHED(I1:I2) = HED(1:NCHA)
      ENDDO
      WRITE(100,1400) PHED(1:NCHA*NOUTST)    
      
      ENDIF


C	PRINT TITLE FOR TRAIN OUTPUT -- OUTRAIN.OUT
		WRITE(500,500)
		WRITE(500,700)
		WRITE(500,1000) VELKMH,VEL,TALENG,SOLT,BIDIS,BITIME,TLENGTH,
     1					TRLENG,BRTIME,EXTDIS,TEXTRA
     
      IF(NOUTVH.GT.0) THEN

      WRITE(500,1500)

      DO IOUTVH = 1,NOUTVH
          NAME = TRNOHED(IOUTVH)
	    NAML = LEN_TRIM(NAME)
    	
          NCHA = 40
          HED(1:NAML) = NAME(1:NAML)
          HED(NAML+1:NCHA) = ' '
          
          I1 = 1    + NCHA*(IOUTVH-1)
          I2 = NCHA + NCHA*(IOUTVH-1)
          PHED(I1:I2) = HED(1:NCHA)
      ENDDO
      WRITE(500,1300) PHED(1:NCHA*NOUTVH)  
      
      
      NCHA = 40
      HED(1:4) = 'DIS.' 
      HED(5:11) = ' '
      HED(12:15) = 'VEL.' 
      HED(16:22) = ' '
      HED(23:26) = 'ACC.' 
      HED(27:40) = ' '
      DO IOUTVH = 1,NOUTVH
          I1 = 1    + NCHA*(IOUTVH-1)
          I2 = NCHA + NCHA*(IOUTVH-1)
          PHED(I1:I2) = HED(1:NCHA)
      ENDDO
      WRITE(500,1400) PHED(1:NCHA*NOUTVH)   
      
      ENDIF
      
      

  200	FORMAT (/,1X,'>>> SOLVING VEHICLE-BRIDGE DYNAMIC INTERACTION')

  500	FORMAT (//1H#,5X,50(1H*)/
     + 1H#,5X,1H*,48X,1H*/
     + 1H#,5X,49H* DYNAMIC ANALYSIS OF VEHICLE-BRIDGE INTERACTION 1H*/ 
     + 1H#,5X,1H*,48X,1H*/
     + 1H#,5X,50(1H*)/)

  600	FORMAT (1H#,15X,27(1H+)/
     + 1H#,15X,25H+ JOB PROGRESS OF BRIDGE 1X,1H+/
     + 1H#,15X,27(1H+)/)

  700	FORMAT (1H#,15X,26(1H@)/
     + 1H#,15X,24H@ JOB PROGRESS OF TRAIN 1X,1H@/
     + 1H#,15X,26(1H@)/)
     
c	changed next 31Oct07
 1000 FORMAT (2X,'TRAIN SPEED  =',F10.4,'(km/hr) =',F10.4,'(m/s)',/
     +2X,'TOTAL ANALYSIS LENGTH (m) . . . . . . . = ',F10.4,'
     +; TOTAL ANALYSIS TIME (s) . . . . . .  = ',F10.4/
     +2X,'Initial distance before bridge (m). . . = ',F10.4,'
     +; Running time before on bridge (s). . = ',F10.4/
     +2X,'Traffic length of bridge (m). . . . . . = ',F10.4/
     +2X,'Train length (distance of wheels)(m)  . = ',F10.4,'
     +; Running time of wheels on bridge (s) = ',F10.4/
     +2X,'Extra running distance after bridge (m) = ',F10.4,'
     +; Extra running time after bridge (s). = ',F10.4/)

 1200 FORMAT (//2X,29(1H-),/2X,'OUTPUT OF STRUCTURAL RESPONSE',/2X,29(1H-))
 
 1300 FORMAT (/2X,12X,'List of Output',11X,6X,A)
 1400 FORMAT (/5X,'TIME(s)',3X,'POS.-1st-Wheel(m)',11X,A)
 
 1500 FORMAT (//2X,26(1H-),/2X,'OUTPUT OF VEHICLE RESPONSE',/2X,26(1H-))
 
	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE PRNOUTBTI(ID,NEQ,NDOFM,RS,RT,NOUTST,NOUTVH,KOST,KOVH,ITE,JTE,
     1                     NSOLT,DCONV,TIME,WPOS,VOUTMAXS,VOUTMAXV) 
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO:
C		- PRINT OUT FOR BRIDGE/TRAIN IN DYNAMIC INTERACTION ANALYSIS
C	-----------------------------------------------------------------------
      DIMENSION ID(1),RS(1),KOST(2,NOUTST),KOVH(NOUTVH),RT(NDOFM,3)
C	-----------------------------------------------------------------------
      DIMENSION BOUT(3,NOUTST),VOUT(3,NOUTVH),VOUTMAXS(6,NOUTST),VOUTMAXV(6,NOUTVH)


      
      VITE  = DFLOAT(ITE)
      VSOLT = DFLOAT(NSOLT)
      PROGS = 100.0*VITE/VSOLT
      WRITE(*,1000) PROGS,JTE,DCONV

      IF(NOUTST.GT.0) THEN
      DO IOUTST = 1,NOUTST
          NCNOD = KOST(1,IOUTST) !NODE
          NCDOF = KOST(2,IOUTST) !DOF
          CALL GETSTDIS(ID,NCNOD,NCDOF,NEQ,RS,BOUT(1,IOUTST),'DIS')
          CALL GETSTDIS(ID,NCNOD,NCDOF,NEQ,RS,BOUT(2,IOUTST),'VEL')
          CALL GETSTDIS(ID,NCNOD,NCDOF,NEQ,RS,BOUT(3,IOUTST),'ACC')
          
          IF(ITE.EQ.1) THEN
              VOUTMAXS(1:3,IOUTST) = BOUT(1:3,IOUTST)
              VOUTMAXS(4:6,IOUTST) = BOUT(1:3,IOUTST)
          ELSE
              DO I = 1,3
                IF(BOUT(I,IOUTST).GT.VOUTMAXS(I+0,IOUTST)) VOUTMAXS(I+0,IOUTST) = BOUT(I,IOUTST)
                IF(BOUT(I,IOUTST).LT.VOUTMAXS(I+3,IOUTST)) VOUTMAXS(I+3,IOUTST) = BOUT(I,IOUTST)
              ENDDO
          ENDIF
      ENDDO
      WRITE(100,1100) TIME,WPOS,(BOUT(1:3,IOUTST),IOUTST=1,NOUTST)
      ENDIF
      

      IF(NOUTVH.GT.0) THEN
      DO IOUTVH = 1,NOUTVH
          NCDOF = KOVH(IOUTVH) !DOF
          VOUT(1,IOUTVH) = RT(NCDOF,1) !DIS
          VOUT(2,IOUTVH) = RT(NCDOF,2) !VEL
          VOUT(3,IOUTVH) = RT(NCDOF,3) !ACC
  
          IF(ABS(VOUT(1,IOUTVH)).LT.1.0E-30) VOUT(1,IOUTVH) = 0.0D0
          IF(ABS(VOUT(2,IOUTVH)).LT.1.0E-30) VOUT(2,IOUTVH) = 0.0D0
          IF(ABS(VOUT(3,IOUTVH)).LT.1.0E-30) VOUT(3,IOUTVH) = 0.0D0
      
          IF(ITE.EQ.1) THEN
              VOUTMAXV(1:3,IOUTVH) = VOUT(1:3,IOUTVH)
              VOUTMAXV(4:6,IOUTVH) = VOUT(1:3,IOUTVH)
          ELSE
              DO I = 1,3
                IF(VOUT(I,IOUTVH).GT.VOUTMAXV(I+0,IOUTVH)) VOUTMAXV(I+0,IOUTVH) = VOUT(I,IOUTVH)
                IF(VOUT(I,IOUTVH).LT.VOUTMAXV(I+3,IOUTVH)) VOUTMAXV(I+3,IOUTVH) = VOUT(I,IOUTVH)
              ENDDO
          ENDIF
      ENDDO
      WRITE(500,1100) TIME,WPOS,(VOUT(1:3,IOUTVH),IOUTVH=1,NOUTVH)
      ENDIF
       
1000  FORMAT(3X,'PROGRESS %',F6.2,2X,'Number of Iteration',I4,3X,'Relative Different %',E12.4)     

1100  FORMAT(3X, E11.4 ,3X, 3X,E12.4,X, 6X, 100(3E12.4,4X))

 
	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE PRNMAXBTI(NOUTST,NOUTVH,VOUTMAXS,VOUTMAXV,WHPOS,WHSLOD,FCMAX,NWHEEL) 
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO:
C		- PRINT OUT FOR BRIDGE/TRAIN IN DYNAMIC INTERACTION ANALYSIS LAST STEP
C         FOR MAXIMUM VALUE
C	-----------------------------------------------------------------------
C	-----------------------------------------------------------------------
      DIMENSION VOUTMAXS(6,NOUTST),VOUTMAXV(6,NOUTVH)
      DIMENSION WHPOS(NWHEEL),WHSLOD(NWHEEL),FCMAX(2,NWHEEL)
      DIMENSION ABSMAXS(3,NOUTST),ABSMAXV(3,NOUTVH)
      
      IF(NOUTST.GT.0) THEN
      WRITE(100,*)
      WRITE(100,1100) (VOUTMAXS(1:3,IOUTST),IOUTST=1,NOUTST)
      WRITE(100,1200) (VOUTMAXS(4:6,IOUTST),IOUTST=1,NOUTST)
          DO IOUTST=1,NOUTST
              DO I = 1,3
                  ABSMAXS(I,IOUTST) = ABS(VOUTMAXS(I,IOUTST))
                  IF(ABS(VOUTMAXS(I+3,IOUTST)).GT.ABSMAXS(I,IOUTST)) 
     1                  ABSMAXS(I,IOUTST) = ABS(VOUTMAXS(I+3,IOUTST))
              ENDDO
          ENDDO
      WRITE(100,1300) (ABSMAXS(1:3,IOUTST),IOUTST=1,NOUTST)
      ENDIF

      IF(NOUTVH.GT.0) THEN
      WRITE(500,*)
      WRITE(500,1100) (VOUTMAXV(1:3,IOUTVH),IOUTVH=1,NOUTVH)
      WRITE(500,1200) (VOUTMAXV(4:6,IOUTVH),IOUTVH=1,NOUTVH)
          DO IOUTVH=1,NOUTVH
              DO I = 1,3
                  ABSMAXV(I,IOUTVH) = ABS(VOUTMAXV(I,IOUTVH))
                  IF(ABS(VOUTMAXV(I+3,IOUTVH)).GT.ABSMAXV(I,IOUTVH)) 
     1                  ABSMAXV(I,IOUTVH) = ABS(VOUTMAXV(I+3,IOUTVH))
              ENDDO
          ENDDO
      WRITE(500,1300) (ABSMAXV(1:3,IOUTVH),IOUTVH=1,NOUTVH)
      ENDIF
      
      WRITE(500,2000)
      WRITE(500,2500)
      DO IWHEEL = 1,NWHEEL
        WRITE(500,3000) IWHEEL,WHPOS(IWHEEL)-WHPOS(1),-WHSLOD(IWHEEL),FCMAX(1:2,IWHEEL)
      ENDDO
      
1100  FORMAT('*MAXIMUM VALUE*', 3X, 12X ,X, 6X, 100(3E12.4,4X))
1200  FORMAT('*MINIMUM VALUE*', 3X, 12X ,X, 6X, 100(3E12.4,4X))
1300  FORMAT('*ABSOLUT VALUE*', 3X, 12X ,X, 6X, 100(3E12.4,4X))

2000  FORMAT(//2X,38(1H-),/2X,'MAXIMUM & MINIMUM WHEEL CONTACT FORCES',/2X,38(1H-))
2500  FORMAT(/2X,'WHEEL No.',2X,'Relative-Position',2X,'Static-Axle-Load',2X,
     1       'Maximum-Contact-Force',2X,'Minimum-Contact-Force') 
3000  FORMAT(X,I6,5X,E15.6,4X,E15.6,6X,E15.6,7X,E15.6)

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================

      SUBROUTINE TRAINRENDER(NTRAIN,NCARB,NBOGI,NWHEL,WHPOS,CADI,PRSS,PRPS,PRWH,TRNDIM,OPTN)
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TRNTYP,OPTN

	DIMENSION NTRAIN(3),WHPOS(1),II(8),TRNDIM(3*8,1)
      DIMENSION CADI(NCARB,1),PRSS(NBOGI,1),PRPS(NWHEL,1),PRWH(NWHEL,1)
      DIMENSION XX(8),YY(8),ZZ(8)
      
      ALLOCATABLE XCC(:),XCB(:),XCW(:),LCARB(:,:),LBOGW(:,:)
      
      ITO = 110
      
      TRNTYP = 'INDP'
      IF(NTRAIN(2).GT.3) TRNTYP = 'ARTG'
            
c      SELECTCASE(TRNTYP)
c          CASE('INDP') !INDEPENDENT TRAIN
c              NCARB = NTRAIN(1)+NTRAIN(3)
c              NBOGI = 2*NCARB
c              NWHEL = 2*NBOGI
c          CASE('ARTG') !ARTIGULATED TRAIN
c              NCARB = NTRAIN(1)+NTRAIN(1)+NTRAIN(3)
c              NBOGI = 2*NTRAIN(1) + (1+1*NTRAIN(2)) + 2*NTRAIN(3)
c              NWHEL = 2*NBOGI
c      ENDSELECT

      ALLOCATE (XCC(NCARB),XCB(NBOGI),XCW(NWHEL),LCARB(2,NCARB),LBOGW(2,NBOGI))
C     OBTAIN THE LONGITUDINAL POSITION OF EACH CAR&BOGIE&WHEEL
C     XCC = POSITION OF THE FRONT SIDE OF EACH CAR BODY	
C     XCB = POSITION OF CENTER OF EACH BOGIE
C     XCW = POSITION OF CENTER OF EACH WHEEL
      CALL TRAINPOSITON(NTRAIN,NBOGI,NWHEL,WHPOS,XCC,XCB,XCW,LCARB,LBOGW,TRNTYP)

      SELECTCASE(OPTN)
      CASE('NODE')
          WRITE(ITO,10000) 
          WRITE(ITO,10100)  
          WRITE(10,10000) 
          WRITE(10,10100)  
      CASE('ELEM')
          WRITE(ITO,10300) 
          WRITE(10,10300) 
      ENDSELECT
      
              
      ICUB = 0
      INOD = 0
C	------------------------------------------------------------             
C	ALL CAR BODY;                                
C	------------------------------------------------------------  
      DO 1000 ICARB = 1,NCARB
          
          ICUB = ICUB + 1
          
          SA = CADI(ICARB,1) !CAR CG TO FRONT BOGIE CG
          SB = CADI(ICARB,2) !CAR CG TO REAR  BOGIE CG
          
          HH = CADI(ICARB,5) !HIEGHT OF CAR BODY
          H2 = 0.5*HH        !HALF F CAR BODY HEIGHT
          
          
          IBOGI = LCARB(1,ICARB)
          IWHEL = LBOGW(1,IBOGI)
          H1A = CADI(ICARB,6)
          H2A = PRSS(IBOGI,9)
          H3A = PRPS(IWHEL,8)
           RR = PRWH(IWHEL,5) !WHEEL RADIUS
          YC = RR+H1A+H2A+H3A
          
          B2A = PRSS(IBOGI,8)
           B2 = B2A


          XX(1:8) = [ SA,-SB,-SB, SA, SA,-SB,-SB, SA]
          YY(1:8) = [ H2, H2, H2, H2,-H2,-H2,-H2,-H2]
          ZZ(1:8) = [-B2,-B2, B2, B2,-B2,-B2, B2, B2]
          
          DO IN = 1,8
              TRNDIM(3*IN-2,ICUB) = XX(IN)
              TRNDIM(3*IN-1,ICUB) = YY(IN)
              TRNDIM(3*IN-0,ICUB) = ZZ(IN)
          ENDDO
          
          XC = XCC(ICARB) - SA
          
          SELECTCASE(OPTN)
          CASE('NODE')
              DO IN = 1,8
                  INOD = INOD + 1
                  XP = XX(IN)+XC
                  YP = YY(IN)+YC
                  ZP = ZZ(IN)
                  WRITE(ITO,6000) INOD,XP,YP,ZP
                  WRITE(10,6000) INOD,XP,YP,ZP
              ENDDO

          CASE('ELEM')
              DO IN = 1,8
                  INOD = INOD + 1
                  II(IN) = INOD
              ENDDO
              WRITE(ITO,7000) ICUB,II(1:8)
              WRITE(10,7000) ICUB,II(1:8)
          ENDSELECT  
C	------------------------------------------------------------  
1000  CONTINUE
C	------------------------------------------------------------  
      
      
C	------------------------------------------------------------             
C	ALL BOGIE;                                
C	------------------------------------------------------------  
      DO 2000 IBOGI = 1,NBOGI
          
          ICUB = ICUB + 1
      
          IWHEL = LBOGW(1,IBOGI)
          
          SA = PRPS(IWHEL,9)
          SB = PRPS(IWHEL,9)
          
          H2A = PRSS(IBOGI,9)
          H3A = PRPS(IWHEL,8)
           RR = PRWH(IWHEL,5) !WHEEL RADIUS
          H2 = H3A - RR + 0.5*H2A
          YC = RR+H3A
          
          B2 = PRPS(IWHEL,7)
                    
          XX(1:8) = [ SA,-SB,-SB, SA, SA,-SB,-SB, SA]
          YY(1:8) = [ H2, H2, H2, H2,-H2,-H2,-H2,-H2]
          ZZ(1:8) = [-B2,-B2, B2, B2,-B2,-B2, B2, B2]

          DO IN = 1,8
              TRNDIM(3*IN-2,ICUB) = XX(IN)
              TRNDIM(3*IN-1,ICUB) = YY(IN)
              TRNDIM(3*IN-0,ICUB) = ZZ(IN)
          ENDDO
                    
          XC = XCB(IBOGI)
          
          SELECTCASE(OPTN)
          CASE('NODE')
              DO IN = 1,8
                  INOD = INOD + 1
                  XP = XX(IN)+XC
                  YP = YY(IN)+YC
                  ZP = ZZ(IN)
                  WRITE(ITO,6000) INOD,XP,YP,ZP
                  WRITE(10,6000) INOD,XP,YP,ZP
              ENDDO
          
          CASE('ELEM')
              DO IN = 1,8
                  INOD = INOD + 1
                  II(IN) = INOD
              ENDDO
              WRITE(ITO,7000) ICUB,II(1:8)
              WRITE(10,7000) ICUB,II(1:8)
          ENDSELECT  
C	------------------------------------------------------------  
2000  CONTINUE
C	------------------------------------------------------------  
    
C	------------------------------------------------------------             
C	ALL WHEEL;                                
C	------------------------------------------------------------  
      DO 3000 IWHEL = 1,NWHEL
          
          ICUB = ICUB + 1
          
          RR = PRWH(IWHEL,5) !WHEEL RADIUS
          SA = RR 
          SB = RR 
          
          H2 = RR            !HALF OF HEIGHT OF WHEEL
          
          B2 = PRWH(IWHEL,6) !HALF OF WIDTH Of WHEEL PAIR
          
          XX(1:8) = [ SA,-SB,-SB, SA, SA,-SB,-SB, SA]
          YY(1:8) = [ H2, H2, H2, H2,-H2,-H2,-H2,-H2]
          ZZ(1:8) = [-B2,-B2, B2, B2,-B2,-B2, B2, B2]
              
          DO IN = 1,8
              TRNDIM(3*IN-2,ICUB) = XX(IN)
              TRNDIM(3*IN-1,ICUB) = YY(IN)
              TRNDIM(3*IN-0,ICUB) = ZZ(IN)
          ENDDO
          
          XC = XCW(IWHEL)
            
          SELECTCASE(OPTN)
          CASE('NODE')
              DO IN = 1,8
                  INOD = INOD + 1
                  XP = XX(IN)+XC
                  YP = YY(IN)
                  ZP = ZZ(IN)
                  WRITE(ITO,6000) INOD,XP,YP,ZP
                  WRITE(10,6000) INOD,XP,YP,ZP
              ENDDO
                
          CASE('ELEM')
              DO IN = 1,8
                  INOD = INOD + 1
                  II(IN) = INOD
              ENDDO
              WRITE(ITO,7000) ICUB,II(1:8)
              WRITE(10,7000) ICUB,II(1:8)
          ENDSELECT    
C	------------------------------------------------------------  
3000  CONTINUE
C	------------------------------------------------------------  


      DEALLOCATE (XCC,XCB,XCW,LCARB,LBOGW)
      

      SELECTCASE(OPTN)
      CASE('NODE')
          WRITE(ITO,10200)  
          WRITE(10,10200) 
      CASE('ELEM')
          WRITE(ITO,10400) 
          WRITE(10,10400) 
      ENDSELECT
      
      
6000  FORMAT(2X,I10,4X,3E20.6)
7000  FORMAT(2X,I6,4X,8(X,I6))


10000 FORMAT('MESH dimension 3 ElemType Hexahedra Nnode 8')
10100 FORMAT('Coordinates')
10200 FORMAT('end coordinates'/)
10300 FORMAT('Elements')
10400 FORMAT('end elements'/)

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================

      SUBROUTINE TRAINRENRES(NTRAIN,NCARB,NBOGI,NWHEL,RT,ITIME,VEL,TIME,TRNDIM)
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TRNTYP,OPTN

	DIMENSION NTRAIN(3),RT(1),ITDOF(5),TRNDIM(3*8,1)
      DIMENSION CUBDIS(6),DISR(3),DIST(3),VARM(3),XX(8),YY(8),ZZ(8)
      
      
      ITO = 110
      
      TRNTYP = 'INDP'
      IF(NTRAIN(2).GT.3) TRNTYP = 'ARTG'



C     DISTANCE TRAVELED BY TRAIN      
      DISX = 0.0 !VEL*TIME


      IF(ITIME.EQ.1) WRITE(ITO,10400) 
      IF(ITIME.EQ.1) WRITE(10,10400) 
      WRITE(ITO,10500) ITIME
      WRITE(ITO,10600) 
      WRITE(ITO,10700) 
      WRITE(10,10500) ITIME
      WRITE(10,10600) 
      WRITE(10,10700) 

      
      ICUB = 0
      INOD = 0
C	------------------------------------------------------------             
C	ALL CAR BODY;                                
C	------------------------------------------------------------  
      DO 1000 ICARB = 1,NCARB
          
          ICUB = ICUB + 1
          
          DO IN = 1,8
              XX(IN) = TRNDIM(3*IN-2,ICUB) 
              YY(IN) = TRNDIM(3*IN-1,ICUB) 
              ZZ(IN) = TRNDIM(3*IN-0,ICUB) 
          ENDDO
          
          CALL TRAINDOFCAR(ICARB,NTRAIN,ITDOF,'CAR',TRNTYP)
          
          I1 = ITDOF(1) !LATERAL    W
          I2 = ITDOF(2) !ROLLING    RX
          I3 = ITDOF(3) !YAWING     RY
          I4 = ITDOF(4) !VERTICAL   V
          I5 = ITDOF(5) !PITCHING   RZ
          
          CUBDIS(1) = DISX   !LONGITUD   U
          CUBDIS(2) = RT(I4) !VERTICAL   V
          CUBDIS(3) = RT(I1) !LATERAL    W
          CUBDIS(4) = RT(I2) !ROLLING    RX
          CUBDIS(5) = RT(I3) !YAWING     RY
          CUBDIS(6) = RT(I5) !PITCHING   RZ
          IF(TRNTYP.EQ.'ARTG') CUBDIS(6) =-RT(I5)
          
          DO IN = 1,8
            INOD = INOD + 1
              DIST(1:3) = CUBDIS(1:3)
              VARM(1:3) = [XX(IN),YY(IN),ZZ(IN)]
              CALL VECPRD(CUBDIS(4),VARM,DISR)
              DIST(1:3) = DIST(1:3) + DISR(1:3)
              
              WRITE(ITO,6000) INOD,DIST(1:3)
              WRITE(10,6000) INOD,DIST(1:3)
          ENDDO
          
C	------------------------------------------------------------  
1000  CONTINUE
C	------------------------------------------------------------  
      
      
C	------------------------------------------------------------             
C	ALL BOGIE;                                
C	------------------------------------------------------------  
      DO 2000 IBOGI = 1,NBOGI
          
          ICUB = ICUB + 1
      
          DO IN = 1,8
              XX(IN) = TRNDIM(3*IN-2,ICUB) 
              YY(IN) = TRNDIM(3*IN-1,ICUB) 
              ZZ(IN) = TRNDIM(3*IN-0,ICUB) 
          ENDDO
      
          CALL TRAINDOFCAR(IBOGI,NTRAIN,ITDOF,'BOG',TRNTYP)
          

          I1 = ITDOF(1) !LATERAL    W
          I2 = ITDOF(2) !ROLLING    RX
          I3 = ITDOF(3) !YAWING     RY
          I4 = ITDOF(4) !VERTICAL   V
          I5 = ITDOF(5) !PITCHING   RZ
          
          CUBDIS(1) = DISX   !LONGITUD   U
          CUBDIS(2) = RT(I4) !VERTICAL   V
          CUBDIS(3) = RT(I1) !LATERAL    W
          CUBDIS(4) = RT(I2) !ROLLING    RX
          CUBDIS(5) = RT(I3) !YAWING     RY
          CUBDIS(6) = RT(I5) !PITCHING   RZ
          
          DO IN = 1,8
            INOD = INOD + 1
              DIST(1:3) = CUBDIS(1:3)
              VARM(1:3) = [XX(IN),YY(IN),ZZ(IN)]
              CALL VECPRD(CUBDIS(4),VARM,DISR)
              DIST(1:3) = DIST(1:3) + DISR(1:3)
              
              WRITE(ITO,6000) INOD,DIST(1:3)
              WRITE(10,6000) INOD,DIST(1:3)
          ENDDO
          
C	------------------------------------------------------------  
2000  CONTINUE
C	------------------------------------------------------------  
    
C	------------------------------------------------------------             
C	ALL WHEEL;                                
C	------------------------------------------------------------  
      DO 3000 IWHEL = 1,NWHEL
          
          ICUB = ICUB + 1
          
          DO IN = 1,8
              XX(IN) = TRNDIM(3*IN-2,ICUB) 
              YY(IN) = TRNDIM(3*IN-1,ICUB) 
              ZZ(IN) = TRNDIM(3*IN-0,ICUB) 
          ENDDO
      
          CALL TRAINDOFCAR(IWHEL,NTRAIN,ITDOF,'WHL',TRNTYP)
          
          I1 = ITDOF(1) !LATERAL    W
          I2 = ITDOF(2) !ROLLING    RX
          I3 = ITDOF(3) !VERTICAL   V
          
          CUBDIS(1) = DISX   !LONGITUD   U
          CUBDIS(2) = RT(I3) !VERTICAL   V
          CUBDIS(3) = RT(I1) !LATERAL    W
          CUBDIS(4) = RT(I2) !ROLLING    RX
          CUBDIS(5) = 0.0D0  !YAWING     RY
          CUBDIS(6) = 0.0D0  !PITCHING   RZ
          
          DO IN = 1,8
            INOD = INOD + 1
              DIST(1:3) = CUBDIS(1:3)
              VARM(1:3) = [XX(IN),YY(IN),ZZ(IN)]
              CALL VECPRD(CUBDIS(4),VARM,DISR)
              DIST(1:3) = DIST(1:3) + DISR(1:3)
              
              WRITE(ITO,6000) INOD,DIST(1:3)
              WRITE(10,6000) INOD,DIST(1:3)
          ENDDO

          
C	------------------------------------------------------------  
3000  CONTINUE
C	------------------------------------------------------------  




      WRITE(ITO,10800) 
      WRITE(10,10800)

      
6000  FORMAT(2X,I10,4X,3E20.6)

10400 FORMAT('GiD Post Results File 1.0'//)
10500 FORMAT('Result "Displacement" "Loadcase 1"',X,I5,X,'Vector OnNodes')
10600 FORMAT('ComponentNames "U" "V" "W"')
10700 FORMAT('Values')
10800 FORMAT('End Values'/)

      
	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================

      SUBROUTINE TRAINDOFCAR(INDX,NTRAIN,ITDOF,TYP,TRNTYP)
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 TYP
	CHARACTER*4 TRNTYP
C     RETURN NUMBER OF TRAIN CAR BODY	
	DIMENSION NTRAIN(3),ITDOF(1)

      
      SELECTCASE(TRNTYP)
      CASE('INDP') !INDEPENDENT TRAIN
      
      IDOF = 0
      JCAR = 0
      JBOG = 0
      JWHL = 0
      DO ICAR = 1,NTRAIN(1)+NTRAIN(3)
          JCAR = JCAR + 1
          IF(TYP.EQ.'CAR'.AND.JCAR.EQ.INDX) GOTO 100
          IDOF = IDOF + 5           !5 DOF PER CAR BODY
          DO IBOG = 1,2
              JBOG = JBOG + 1
              IF(TYP.EQ.'BOG'.AND.JBOG.EQ.INDX) GOTO 100
              IDOF = IDOF + 5           !5 DOF PER BOGIE
              DO IWHL = 1,2
                JWHL = JWHL + 1
                IF(TYP.EQ.'WHL'.AND.JWHL.EQ.INDX) GOTO 100
                IDOF = IDOF + 3           !3 DOF PER WHEEL
              ENDDO
          ENDDO
      ENDDO
      
100   CONTINUE  


      CASE('ARTG') !ARTIGULATED TRAIN
      
      IDOF = 0
      JCAR = 0
      JBOG = 0
      JWHL = 0
      
      DO ICAR = 1,NTRAIN(1)
          JCAR = JCAR + 1
          IF(TYP.EQ.'CAR'.AND.JCAR.EQ.INDX) GOTO 200
          IDOF = IDOF + 5           !5 DOF PER CAR BODY
          DO IBOG = 1,2
              JBOG = JBOG + 1
              IF(TYP.EQ.'BOG'.AND.JBOG.EQ.INDX) GOTO 200
              IDOF = IDOF + 5           !5 DOF PER BOGIE
              DO IWHL = 1,2
                JWHL = JWHL + 1
                IF(TYP.EQ.'WHL'.AND.JWHL.EQ.INDX) GOTO 200
                IDOF = IDOF + 3           !3 DOF PER WHEEL
              ENDDO
          ENDDO
      ENDDO


      DO ICAR = 1,NTRAIN(2)
          DO IBOG = 1,1  !LAST BOGIE OF ARTIGULATED CARS
              DO IWHL = 1,2
                JWHL = JWHL + 1
                IF(TYP.EQ.'WHL'.AND.JWHL.EQ.INDX) GOTO 200
                IDOF = IDOF + 3           !3 DOF PER WHEEL
              ENDDO
              JBOG = JBOG + 1
              IF(TYP.EQ.'BOG'.AND.JBOG.EQ.INDX) GOTO 200
              IDOF = IDOF + 5           !5 DOF PER BOGIE
          ENDDO
          JCAR = JCAR + 1
          IF(TYP.EQ.'CAR'.AND.JCAR.EQ.INDX) GOTO 200  
          IDOF = IDOF + 5           !5 DOF PER CAR BODY        
      ENDDO
      
C     LAST   BOGIE OF TRANSITIONAL CAR
      JBOG = JBOG + 1
      IF(TYP.EQ.'BOG'.AND.JBOG.EQ.INDX) GOTO 200
      IDOF = IDOF + 5           !5 DOF PER BOGIE
C     FIRST  WHEEL OF TRANSITIONAL CAR LAST BOGIE
      JWHL = JWHL + 1
      IF(TYP.EQ.'WHL'.AND.JWHL.EQ.INDX) GOTO 200
      IDOF = IDOF + 3           !3 DOF PER WHEEL
C     SECOND WHEEL OF TRANSITIONAL CAR LAST BOGIE
      JWHL = JWHL + 1
      IF(TYP.EQ.'WHL'.AND.JWHL.EQ.INDX) GOTO 200
      IDOF = IDOF + 3           !3 DOF PER WHEEL


      DO ICAR = 1,NTRAIN(3)
          JCAR = JCAR + 1
          IF(TYP.EQ.'CAR'.AND.JCAR.EQ.INDX) GOTO 200
          IDOF = IDOF + 5           !5 DOF PER CAR BODY
          DO IBOG = 1,2
              JBOG = JBOG + 1
              IF(TYP.EQ.'BOG'.AND.JBOG.EQ.INDX) GOTO 200
              IDOF = IDOF + 5           !5 DOF PER BOGIE
              DO IWHL = 1,2
                JWHL = JWHL + 1
                IF(TYP.EQ.'WHL'.AND.JWHL.EQ.INDX) GOTO 200
                IDOF = IDOF + 3           !3 DOF PER WHEEL
              ENDDO
          ENDDO
      ENDDO

    
200   CONTINUE  
      
      ENDSELECT
      
      
      SELECTCASE(TYP)
          CASE('CAR')
              ITDOF(1) = IDOF+1 
              ITDOF(2) = IDOF+2 
              ITDOF(3) = IDOF+3 
              ITDOF(4) = IDOF+4 
              ITDOF(5) = IDOF+5
          CASE('BOG')
              ITDOF(1) = IDOF+1 
              ITDOF(2) = IDOF+2 
              ITDOF(3) = IDOF+3 
              ITDOF(4) = IDOF+4 
              ITDOF(5) = IDOF+5
          CASE('WHL')
              ITDOF(1) = IDOF+1 
              ITDOF(2) = IDOF+2 
              ITDOF(3) = IDOF+3 
      ENDSELECT  
        
      

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================



      SUBROUTINE TRAINPOSITON(NTRAIN,NBOGI,NWHEL,WHPOS,XCC,XCB,XCW,LCARB,LBOGW,TRNTYP)
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TRNTYP
C     XCC = POSITION OF THE FRONT SIDE OF EACH CAR BODY	
C     XCB = POSITION OF CENTER OF EACH BOGIE
C     XCW = POSITION OF CENTER OF EACH WHEEL
      DIMENSION NTRAIN(3),WHPOS(1)
	DIMENSION XCC(1),XCB(1),XCW(1),LCARB(2,1),LBOGW(2,1)


      DO IWHEL = 1,NWHEL
        XCW(IWHEL) = WHPOS(IWHEL)
      ENDDO
      
      DO IBOGI = 1,NBOGI
        XCB(IBOGI) = 0.5*( XCW(2*IBOGI-1) + XCW(2*IBOGI-0) )
        LBOGW(1,IBOGI) = 2*IBOGI-1
        LBOGW(2,IBOGI) = 2*IBOGI-0
      ENDDO
      

      SELECTCASE(TRNTYP)
      
      CASE('INDP') !INDEPENDENT TRAIN
      
      JCAR = 0
      JBOG = 0
      DO ICAR = 1,NTRAIN(1)+NTRAIN(3)
          JCAR = JCAR + 1
          DO IBOG = 1,2
              JBOG = JBOG + 1
              IF(IBOG.EQ.1) XCC(JCAR) = XCB(JBOG)
              LCARB(IBOG,JCAR) = JBOG
          ENDDO
      ENDDO
      


      CASE('ARTG') !ARTIGULATED TRAIN
      
      JCAR = 0
      JBOG = 0
      
      DO ICAR = 1,NTRAIN(1)
          JCAR = JCAR + 1
          DO IBOG = 1,2
              JBOG = JBOG + 1
              IF(IBOG.EQ.1) XCC(JCAR) = XCB(JBOG)
              LCARB(IBOG,JCAR) = JBOG
          ENDDO
      ENDDO
      

C     FIRST  BOGIE OF TRANSITIONAL CAR
      JBOG = JBOG + 1

      DO ICAR = 1,NTRAIN(2)
          JCAR = JCAR + 1
          XCC(JCAR) = XCB(JBOG) !ARTIGULATED CAR USE THE PREVIOUS BOGIE AS REFERENCE
          DO IBOG = 1,1  !LAST BOGIE OF ARTIGULATED CARS
              JBOG = JBOG + 1
              LCARB(1,JCAR) = JBOG-1
              LCARB(2,JCAR) = JBOG-0
          ENDDO
      ENDDO


      DO ICAR = 1,NTRAIN(1)
          JCAR = JCAR + 1
          DO IBOG = 1,2
              JBOG = JBOG + 1
              IF(IBOG.EQ.1) XCC(JCAR) = XCB(JBOG)
              LCARB(IBOG,JCAR) = JBOG
          ENDDO
      ENDDO

        
      ENDSELECT
      

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================















