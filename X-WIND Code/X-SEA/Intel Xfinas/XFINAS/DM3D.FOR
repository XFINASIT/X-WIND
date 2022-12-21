      SUBROUTINE DM3D(IDATM,AA,BB,CC)  
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
C	MAIN PROGRAM FOR 3D TRAIN-BRIDGE DYNAMIC INTERACTION ANALYSIS
C	Developed by Mr.Nguyen DV, Asian Institute of Technology
C	(New Version 2007, Train & Bridge are in consistent Coordinates)
C	-------------------------------------------------------------
C	IDMA	= DYNAMIC MOVING MODELLING TYPE
C				1: Vehicle modeled as moving loads (inputed axle loads) with 3D bridge model
C				2: 2D model analysis of train interaction with bridge
C				3: 3D model analysis of train interaction with bridge
C
C	IANA	= ANALYSIS METHOD OF MOVING LOAD (=1 if IDMA =1,3) (= following if IDMA =2)
C				1 = Dynamic Analysis of Train with moving mass
C				2 = Dynamic Analysis of Train with moving forces
C				3 = Static Analysis of Train with moving forces
C				4 = Dynamic Analysis of Moving Mass
C				5 = Static Analysis of Moving mass
C
C	IDSO	= SOLUTION METHOD
C			    1 = Direct Integration
C				2 = Modal superpostion (available for only IDMA=1) 

C	STRUCTURE PARAMETERS:
C	--------------------
C	NSN     = NUMBER OF NODES FOR ENTIRE STRUCTURE
C	NEG     = NUMBER OF ELEMENT GROUPS
C	NBS     = NUMBER OF GLOBAL SUPPORTS
C	NLS     = NUMBER OF LOCAL SUPPORTS (SKEWED SUPPORTS)
C	NSS		=  NUMBER OF SPRING SUPPORTED NODES GLOBALY
C	NLSS	= NUMBER OF SPRING SUPPORTED NODES LOCALY
C	NLA     = NUMBER OF LOCAL AXES
C	NSC     = MAXIMUM NUMBER OF COORDINATES AT ANY ONE NODE
C	NSF     = MAXIMUM NUMBER OF DEGREES OF FREEDOM AT ANY ONE NODE

C     A(LXY)  = XY(NSN,NSC): NODAL CO-ORDINATES X,Y,Z
C     IA(LID) = NRL(NSF,NSN): BOUNDARY CODES AND EQUATION NUMBERS
C     A(LDK)  = D(NEQ)     :  DIAGONAL TERMS OF GLOBAL STIFFNESS MATRIX


C	Element Parameters
C	NCO     = NUMBER OF NODAL COORDINATES
C	NNF     = NUMBER OF NODAL DEGREES OF FREEDOM
C	NEF     = NUMBER OF ELEMENT DEGREES OF FREEDOM
C	NNO     = CURRENT NUMBER OF ELEMENT NODES
C	NELE    = NUMBER OF ELEMENTS IN GROUP
C     IA(LLM) = LM(NEF,NELE): ELEMENT CONNECTIVITY FOR CURRENT ELEMENT GROUP
C						    & GLOBAL EQUATION NUMBER FOR ELEMENT DOF

C	NEID(2,NELW)  = NODE NUMBER OF LANE WHERE WHEEL RUNS
C
C	TRAIN PARAMETERS:
C	-----------------
C	NEQTR  = Number of equations of stiffness matrix = 3*(ITRAIN(1)*5) + (2*ITRAIN(2)+1)*5 + 3*(ITRAIN(3)*5)
C	NWKTR  = NUMBER OF ELEMENTS BELOW SKYLINE OF UPPER TRIA. STIFFNESS MATRIX(=ITRAIN(1)*29+(12+ITRAIN(2)*17)+ITRAIN(3)*29)
C	IW(I_NCHTR)= NCHTR(NEQTR) : VECTOR CONTANING COLUMN HEIGHT OF TRAIN'S STIFF/DAMPING MATRIX
C	IW(I_NDTR) = NDTR(NEQTR+1):VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS OF TRAIN STIFFNESS/DAMPING MATRIX
C	IWRIN : Idenfifier for Input of Wheel-Rail contact Mechanism
C	   = 0: Input the contact stiffness & damping directly (OLD)
C	   = 1: Input the wheel-rail properties and calculate contact stiffness & damping

C	OUTPUT CONTROLS	
C	IW(I_PRTR) = IPRTTR(NPRTTR)  : DEGREES OF FREEDOM OF TRAIN TO PRINT
C	IW(I_PRST) = IPRTST(NPRTST,2): NODES NUMBER OF STRUCTURE[NODE #, DOF #]
C	NPRTTR	   = NUMBER OF DEGREES OF FREEDOM TO PRINT OF TRAIN
C	NPRTST	   = NUMBER OF NODES TO PRINT OF STRUCTURE
C	NPRWH      = NUMBER OF DEGREES OF FREEDOM OF WHEEL TO PRINT
C	KPRTR	   = TRAIN PRINT FLAG (in order to control size of train output file)
C				   = 0: displacements of train body is not print, only accelerations are printed
C				   = 1: print
C	
C	--------NOTES (16Apr06)-------
C	+ Current max. number of train DOFs to print = 10, may change if/MAX/AMXTR(10,2)
C	+ Current max. number of bridge nodes to print = 10,may change if/MAX/AMXST(10,2)
C	+ Current max. number of train wheel NWHEEL = 100,may change size if/MAX1/AMXAXL(100,3)
C
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
	COMMON /SPBC/ NSS,NLSS
      COMMON /DMCO/ IDMA,IANA,IDSO,NNBR,NNRA,NEBEBR,NPBEBR,NLRB,
     +			  NEWKTK,NPWKTK,NEBERA,NPBERA,NDOF,NNBE,NNWK
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
	COMMON /BRIW/ I_PELW,I_NEID,I_PRST

	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH

	COMMON /TRW3/ I_PCAB,I_CADI,I_ADDI,I_ALSU,I_PRSS,I_PBOG,I_PRPS,
	1			  I_PRWH,I_WRCO,I_WPOS,I_AMTR,I_ACTR,I_AKTR,I_AXLD

	COMMON /NIAX/ NWHEEL
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG !added 5Nov2005

	COMMON /ITRM/ ITRK	!Created here 25Mar06 Temperarily for ARMAT pointer

      COMMON A(9000000),IA(9000000)

	COMMON /MEMW/ W(7000000),IW(7000000)
C
	DIMENSION AC(1)
	DIMENSION AA(1),BB(1),CC(1)
C	------------------------------------------------------------------------
C	LOADWAY DEFINE
C	--------------
      ITRK = 1
C	Read number of elements to define loaways of bridge
	CALL RENELW

	CALL LOCAW(4,2)
C	Read element number array to define loaways of bridge
	CALL REELELW(IW(I_NEID))

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
  505	IFEIG = 1
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

C	Read additional parameters of bridge for 3D model (eccentricity)
C	CALL  REBRI3D !Moved up 8Mar

C	READ PARAMETERS FOR RAIL/TRACK IRREGULARTIES AND GENERATE
	CALL REIRRE3D					!added Sep06
C	----------------------------
C	FORMULATION OF TRAIN SYSTEM:
C	----------------------------
C	Read control variables for 3D train/truck model:
	CALL REMOVM(LSTATC)

C	Calculate train general parameters (NBOGI,NCARB,NPCAB,NACAB,NABOG,NWHEEL,NEQTR,NWKTR)
C	CALL TRGPARA3D	!Changed next 12mar08
C	IF(IWRIN.EQ.0) THEN !Changed next Dec2008 by Nguyen (new)
	IF(IWRIN.EQ.0.OR.IWRIN.EQ.1) CALL TRGPARA3D
	
C	Allocate storage of train
	CALL LOCAW(4,3)

C	READ TRAIN PROPERTIES:
	CALL RETRN3D (IW(I_ITRAIN),W(I_PCAB),W(I_CADI),W(I_ADDI),W(I_ALSU)
	1			,W(I_PRSS),W(I_PBOG),W(I_PRPS),W(I_PRWH),W(I_WRCO),
     2			 W(I_WPOS),IDATM)
  
C	Next added Dec2008 to cal. static axle forces for wheel-rail contacts
	CALL CALSAXL(W(I_PCAB),W(I_CADI),W(I_PBOG),W(I_PRWH),SCALE,
     1			 W(I_AXLD))

C	Next added Dec2008 to cal. wheel-rail contac mechanism if required
	IF(IWRIN.EQ.1)CALL GENWHPRO(W(I_PRWH),W(I_AXLD),W(I_WRCO),VEL) 
  
C	FORM TRAIN MATRICES
      CALL TRNMAT3D(IW(I_ITRAIN),W(I_PCAB),W(I_CADI),W(I_ADDI),
	1			   W(I_ALSU),W(I_PRSS),W(I_PBOG),W(I_PRPS),W(I_PRWH),
     2			   W(I_WRCO),W(I_AMTR),W(I_ACTR),W(I_AKTR)) 

C	FORM VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS OF TRAIN STIFFNESS/DAMPING MATRIX
	CALL MAXATRN(IW(I_ITRAIN),NEQTR,NWKTR,IW(I_NCHTR),IW(I_NDTR))
      
      CALL PRINSTIFC(W(I_AKTR),IW(I_NDTR),NEQTR)
      
C	DEFINE THE LOADWAY OF BRIDGE, WHERE THE TRAIN IS MOVING ON IW(I_NEID)
C		IPATH = 1 : LOAD RUNS ON BRIDGE, OTHERWISE ON RAIL
      IF(IPATH.EQ.1) THEN
	  CALL PATHMOV(NELW,IW(I_NEID),W(I_PELW),A(LXY),TLENGTH)
      ELSE
	  GOTO 1000
      ENDIF
C	-------------------------------------------------
C	SOLVE TRAIN-STRUCTURE INTERACTION IN 3D-MODELLING
C	-------------------------------------------------
C	READ DOFs OF TRAIN AND NODES NUMBER & DOFs OF STRUCTURE TO PRINT
C	CALL REPRCONT(IW(I_PRTR),IW(I_PRST),NPRTTR,NPRTST) !Changed next 03Mar07
	CALL REPRCO(IW(I_PRTR),IW(I_PRST),IW(I_PTWH),NPRTTR,NPRTST,NPRWH)

 	
      CALL TRSOL3D(IA(LMA),W(I_WPOS),
	1			 IW(I_PRTR),IW(I_PRST),IW(I_PTWH),W(I_AXLD),
     2			 AA,BB,CC,'STIF','MASS','DAMP')		!added IWHPR 16Apr06    
C	------------------------------------------------------------------------
 1000	RETURN
	END
C
C	=======================================================================
      SUBROUTINE REIRRE3D
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------
C     TO READ PARAMETERS FOR 3D RAIL/TRACK IRREGULARITIES
C	Programmed 6Sep2006 by NguyenDV, Upgraded 30Jan07
C	-----------------------------------------------------
C	PARAMETERS:
C	IRIN  = GENERATING METHOD OF IRREGULARITY 
C		    =0: No Irregularity considered
C			=1: Text file (Measurement data)
C			=2: Generate by semi-periodic functions
C			=3: Generate by random functions, where:
C				IRAN = 1: by FRA Power Spectral Density (Federal Railway Administration, USA)
C				IRAN = 2: by SBCF Power Spectral Density (French National Railway)
C			=4: Rail joints 
C     NFFT  = Number of spectrums for total rail length
C     DSEED = Seed value for random number ()     
C     WLMIN = Shortest wave length considered (m), correspoding to upper limit of wave number (1/m)= 1/WLMIN 
C     WLMAX = Longest wave length considered (m), correspoding to lower limit of wave number (1/m) =1/WLMAX
C     NOSIM = Number of simulation			
C	IRCLA = Rail irregularities class, according to FRA
c			= 1,2,3: (Not reccommended)
c			= 4: FRA Class = very poor (IRv=4.05; IRh=5.10; IRc=1.50mm)
c			= 5: FRA Class = poor      (IRv=3.38; IRh=4.25; IRc=1.25mm)
c			= 6: FRA Class = moderate  (IRv=2.70; IRh=3.40; IRc=1.00mm)
C	RLENG = Norminal Length of a rail bar (between two rail joint)
C	SLSTEP= Sleeper step (average distance between two adjacent sleepers)
C	DVMIN = Min. distance btw two adjacent deterioted spots in vertical surface
C	DHMIN = Min. distance btw two adjacent deterioted spots in horizontal surface
C	GAUGE =	Nominal (average or design) gauge of rails = inner distance btw. 2 rails
C	-----------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	COMMON /IR3I/ IRIN,IRCLA,RLENG,SLSTEP,DVMIN,DHMIN,GAUGE
C	Previous seperated into Real and Integer common blocks 30Jan07
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM
	COMMON /RIR3/ RLENG,SLSTEP,DVMIN,DHMIN,GAUGE,DSEED,GAMU,GAML 
c	COMMON /RIR3/ RLENG,SLSTEP,DVMIN,DHMIN,GAUGE,DSEED,WLMIN,WLMAX !transfer wave numbers
C	-----------------------------------------------------
C	Read option for generate method
	READ(ITI,*) 
	READ(ITI,*) IRIN 
	WRITE(ITO,*)'READ PARAMETERS FOR 3D RAIL/TRACK IRREGULARITIES'
      WRITE(10,*)'READ PARAMETERS FOR 3D RAIL/TRACK IRREGULARITIES'

	IF(IRIN.LT.0.OR.IRIN.GT.4) THEN  !Check validity of input data
	  WRITE(*,*) 'PLEASE CHECK INPUT TYPE OF RAIL IRREGULARIRIES IRIN'
	  STOP
	ENDIF

	WRITE(ISO,500)
     	
	SELECT CASE (IRIN)
	CASE (0)   !No irregularity
	  WRITE(ITO,*)'RAIL IRREGULARITIES IS NOT CONSIDERED'
        WRITE(10,*)'RAIL IRREGULARITIES IS NOT CONSIDERED'
	  GOTO 5000

      CASE (1)  !Text file (Measurement data)
	  WRITE(ITO,*)'OPTION CURRENTLY NOT AVAILABLE'
        WRITE(10,*)'OPTION CURRENTLY NOT AVAILABLE'
	  STOP

      CASE (2)  !Generate by semi-periodic functions
	  READ(ITI,*) IRCLA,RLENG,SLSTEP,DVMIN,DHMIN,GAUGE
	  IF(IRCLA.EQ.1.OR.IRCLA.EQ.2.OR.IRCLA.EQ.3) THEN
	    WRITE(*,*) '===> TOO LOW QUALITY OF RAIL, TRY FRA CLASS 4,5,6'
	    STOP
	  ENDIF
	  WRITE(ISO,2000)IRIN,IRCLA,RLENG,SLSTEP,DVMIN,DHMIN,GAUGE

      CASE (3)  !Generate by random functions
c	  READ(ITI,*) IRAN,IRCLA,NFFT,DSEED,GAMU,GAML,NOSIM,GAUGE !Changed next 10Nov07
	  READ(ITI,*) IRAN,IRCLA,NFFT,DSEED,WLMIN,WLMAX,NOSIM,GAUGE

	  IF(IRAN.EQ.1) THEN !by FRA Power Spectral Density
		IF(IRCLA.EQ.0.OR.IRCLA.GT.6) IRCLA = 6
		IF(NFFT.EQ.0)    NFFT = 1024
		IF(NFFT.GT.8192) NFFT = 8192
	    IF(DSEED.LE.0.0D0) DSEED = 5.0
c	    IF(GAMU.LE.0.0D0)  GAMU = 0.65617   !changed next 25June07 (account for shorter wave lengths)
c	    IF(GAMU.LE.0.0D0)  GAMU = 2.0		!correspond to lamda,min = 0.5m
c	    IF(GAML.LE.0.0D0)  GAML =0.0032808	!correspond to lamda,max = 304.8m
c		Previous two lines changed next 10Nov07
	    IF(WLMIN.LE.0.0D0)  WLMIN = 0.5		!correspond to GAMU = 2.0 (1/m)
	    IF(WLMAX.LE.0.0D0)  WLMAX = 50.		!correspond to GAML = 0.0032808 (1/m)
C		next added 10Nov07
		IF(WLMIN.GT.WLMAX) THEN
			WRITE(ITO,*)'ERROR! CHECK INPUT WAVE LENGTH LIMITS'
              WRITE(10,*)'ERROR! CHECK INPUT WAVE LENGTH LIMITS'
			STOP
		ENDIF
	    IF(NOSIM.EQ.0) NOSIM = 1
		IF(GAUGE.LE.0.0.OR.GAUGE.GT.1.45) GAUGE = 1.0
	    WRITE(ISO,3100)IRAN,IRCLA

	  ELSEIF(IRAN.EQ.2) THEN !by SBCF Power Spectral Density
		IF(NFFT.EQ.0)    NFFT = 1024
		IF(NFFT.GT.8192) NFFT = 8192
	    IF(DSEED.LE.0.0D0) DSEED = 5.0
c	    IF(GAMU.LE.0.0D0)  GAMU = 10.0
c	    IF(GAML.LE.0.0D0)  GAML = 0.01
c		Previous two lines changed next 10Nov07
	    IF(WLMIN.LE.0.0D0)  WLMIN = 0.1		!correspond to GAMU = 10.0 (1/m)
	    IF(WLMAX.LE.0.0D0)  WLMAX = 100		!correspond to GAML = 0.01 (1/m)
C		next added 10Nov07
		IF(WLMIN.GT.WLMAX) THEN
			WRITE(ITO,*)'ERROR! CHECK INPUT WAVE LENGTH LIMITS'
              WRITE(10,*)'ERROR! CHECK INPUT WAVE LENGTH LIMITS'
			STOP
		ENDIF
	    IF(NOSIM.EQ.0) NOSIM = 1
	    WRITE(ISO,3100)IRAN
	  ENDIF

C		next added 10Nov07
		GAMU = 1./WLMIN
		GAML = 1./WLMAX

c	  WRITE(ISO,3500) NFFT,DSEED,GAMU,GAML,NOSIM,GAUGE !changed next 10Nov07
	  WRITE(ISO,3500) NFFT,DSEED,WLMIN,WLMAX,NOSIM,GAUGE

      CASE (4)  !Rail joints
	  WRITE(ITO,*)'OPTION CURRENTLY NOT AVAILABLE'
        WRITE(10,*)'OPTION CURRENTLY NOT AVAILABLE'
	  STOP

	END SELECT

  500 FORMAT(//,
     +	21X,'PARAMETERS FOR 3D RAIL/TRACK IRREGULARITIES',/
     +	(1H#,19X,45(1H-)))

 2000 FORMAT(/
     +1X,'GENERATING METHOD: by semi-periodic functions;
     +IRIN     = ',I3/
     +1X,'Rail irregularities class, according to FRA',/
     +5X,'= 0: No Irregulartities (all IRR = 0)',/
     +5X,'= 1,2,3: (Not reccommended)',/
     +5X,'= 4: FRA Class = very poor (IRv=4.05; IRh=5.10; IRc=1.50mm)',/
     +5X,'= 5: FRA Class = poor      (IRv=3.38; IRh=4.25; IRc=1.25mm)',/
     +5X,'= 6: FRA Class = moderate  (IRv=2.70; IRh=3.40; IRc=1.00mm)',/
     +3X,'FRA CLASS    = ',I3/
     +1X,'Nominal Length of a rail bar (between two rail joint)
     +RLENG    = ',F10.4/
     +1X,'Sleeper step (average distance between two adjacent sleepers)
     +SLSTEP   = ',F10.4/
     +1X,'Min. distance bwt. two adjacent deterioted spots in vertical 
     +DVMIN    = ',F10.4/
     +1X,'Min. distance btw two adjacent deterioted spots in horizontal 
     +DHMIN    = ',F10.4/
     +1X,'Nominal rail gauge (inner distance between two rails)
     +GAUGE    = ',F10.4)

 3100 FORMAT(/10X,'GENERATING METHOD: BY RANDOM FUNCTIONS, IRIN = 3',/
     +5X,'PSD functions: by Federal Railway Administration (FRA), USA
     +. . . .  IRAN   = ',I3/
     +5X,'Rail irregularities class (by FRA); IRCLA  = ',I3)

 3200 FORMAT(/10X,'GENERATING METHOD: BY RANDOM FUNCTIONS, IRIN = 3',/
     +5X,'PSD functions: by SBCF (French) . .  IRAN   = ',I3)

 3500 FORMAT(/
     +5X,'Number of spectrums for total rail length (NFFT) = ',I10/
     +5X,'Seed value to generate random number (DSEED)     = ',F10.4/
c     +5X,'Upper limit of wave number (GAMU, 1/m)           = ',F10.4/ 
c     +5X,'Lower limit of wave number (GAML, 1/m)           = ',F10.4/ !changed next 10Nov07
     +5X,'Shortest wave length considered (WLMIN, m)        = ',F10.4/
     +5X,'Longest wave length considered (WLMAX, m)         = ',F10.4/
     +5X,'Number of simulations (NOSIM)					 = ',I10/
     +5X,'Nominal (design) gauge of rails (GAUGE,M)        = ',F10.4)
C
 5000 RETURN	
	END
C
C	=======================================================================
      SUBROUTINE RAILFRA(POSI,TMAX,SPEED,GEOx,GEOe,GEOa,GEOc)

      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------
C     PROGRAM TO GENERATE 3D RAIL GEOMETRY USING FRA PSD:
C	Programmed 30Jan2007 by NguyenDV
C	-----------------------------------------------------
C	INPUT:
C	------
C	POSI :  Last wheel position (in train running direction)
C	TMAX :  Total time considered for rail irregularities (Sec.)
C	SPEED:  Train speed (M/SEC)

C	Track irregularity is calculated from X to X+SPEED*TMAX
C	IRCLA:  Rail grade number based on FRA group
C	NFFT :  Number of spectrums used for total rail length (Default value: 1024)
C	DSEED : Seed value to generate random number (= arbitrary value)

C	GAMMA RANGE USED IN FRA PSD
C	GAMU :  Upper limit of wave number GAMMA,default = 0.65617 Cycle/M(=1.524M=5ft) 
C	GAML :  Lower limit of wave number GAMMA,default = 0.0032808 Cycle/M(=304.8M=1000ft)
C	NOSIM : Number of simulations (= arbitrary value, in general NOSIM = 1)

C	OUTPUT:
C	-------
C	GEOx(NFFT) = X-ordinate array for rail profiles
C	GEOe(NFFT) = Elevation (vertical) irregularity profile
C	GEOa(NFFT) = Alignment (lateral) irregularity profile
C	GEOc(NFFT) = Rotational irregularity profile

C	WORKING
C	RR(NR) : TO CALCULATE PHASE(BETA) AS RANDOM NUMBER GENERATED IN SUBROUTINE RND
C	GAMMA=1/LAMBDA=1/(VT)=OMEGA/(2*PI)
C		GAMMA : WAVE NUMBER(1/M) 
C		LAMBDA : WAVE LENGTH(M)
C		V : VELOCITY(M/SEC)
C		T : PERIOD(SEC)
C		OMEGA : SPATIAL FREQUENCY(RAD/SEC) 
C	GEOc(NFFT) = Cross (superelevation) irregularity
C	GEOg(NFFT) = Gauge irregularity
C	-------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM
	COMMON /RIR3/ RLENG,SLSTEP,DVMIN,DHMIN,GAUGE,DSEED,GAMU,GAML

      DIMENSION GEOx(NFFT),GEOe(NFFT),GEOa(NFFT),GEOc(NFFT),GEOg(NFFT)  !later, change size to (1)?
	DIMENSION Re(NFFT),Ra(NFFT),Rc(NFFT),Rg(NFFT)  !later, change size to (1)?
C	-------------------------------------------------------------------------------
      PI2=2.*3.141592654
      DT=TMAX/DFLOAT(NFFT)
      DGAM=(GAMU-GAML)/DFLOAT(NFFT)

C	Re-compute wave length limits (added 10Nov07)
	WLMIN = 1./GAMU
	WLMAX = 1./GAML
C	WRITE HEADERS & INPUT PARAMETERS
c      WRITE(100,1000)
c	WRITE(100,1100)IRCLA,NFFT,DSEED,GAMU,GAML,NOSIM,GAUGE
c	Previous 2 lines changed next 07Nov07
c	WRITE(100,1000)IRCLA,NFFT,DSEED,GAMU,GAML,NOSIM,GAUGE !Changed next 10Nov07
	WRITE(100,1000)IRCLA,NFFT,DSEED,WLMIN,WLMAX,NOSIM,GAUGE
c	Next added 07Nov07
	WRITE(100,1100)
C
C	GENERATION OF AVERAGE TRACK IRREGULARITIES OF 'NOSIM' SIMULATIONS
C	Initialization:
	CALL CLEARA(Re,NFFT)
	CALL CLEARA(Ra,NFFT)
	CALL CLEARA(Rc,NFFT)
	CALL CLEARA(Rg,NFFT)
	CALL CLEARA(GEOx,NFFT)
	CALL CLEARA(GEOe,NFFT)
	CALL CLEARA(GEOa,NFFT)
	CALL CLEARA(GEOc,NFFT)
	CALL CLEARA(GEOg,NFFT)

C	DEFINITION OF FRA PSD constants, Unit: A[1/m3], O1,O2[1/m]
C	(Reference: Fryba, Dynamics of Railway Bridges, pp125)	 
      IF(IRCLA.EQ.1) THEN
C	  Elevation	  
        Ae  = 15.53E-8
        O1e = 23.3E-3
        O2e = 13.1E-2

c	  Alignment
        Aa  = 9.83E-8
        O1a = 32.8E-3
        O2a = 18.4E-2

c	  Cross
        Ac  = 4.52E-8
        O1c = 23.3E-3
        O2c = 13.1E-2

c	  Gauge
        Ag  = 9.83E-8
        O1g = 29.2E-3
        O2g = 23.3E-2
	  	  
      ELSE IF(IRCLA.EQ.2) THEN
C	  Elevation
        Ae  = 8.85E-8
        O1e = 23.3E-3
        O2e = 13.1E-2

c	  Alignment
        Aa  = 5.51E-8
        O1a = 32.8E-3
        O2a = 18.4E-2

c	  Cross
        Ac  = 3.15E-8
        O1c = 23.3E-3
        O2c = 13.1E-2

c	  Gauge
        Ag  = 5.51E-8
        O1g = 29.2E-3
        O2g = 23.3E-2

      ELSE IF(IRCLA.EQ.3) THEN
C	  Elevation
        Ae  = 4.92E-8
        O1e = 23.3E-3
        O2e = 13.1E-2

c	  Alignment
        Aa  = 3.15E-8
        O1a = 32.8E-3
        O2a = 18.4E-2

c	  Cross
        Ac  = 2.16E-8
        O1c = 23.3E-3
        O2c = 13.1E-2

c	  Gauge
        Ag  = 3.15E-8
        O1g = 29.2E-3
        O2g = 23.3E-2

      ELSE IF(IRCLA.EQ.4) THEN
C	  Elevation
        Ae  = 2.75E-8
        O1e = 23.3E-3
        O2e = 13.1E-2

c	  Alignment
        Aa  = 1.77E-8
        O1a = 32.8E-3
        O2a = 18.4E-2

c	  Cross
        Ac  = 1.38E-8
        O1c = 23.3E-3
        O2c = 13.1E-2

c	  Gauge
        Ag  = 1.77E-8
        O1g = 29.2E-3
        O2g = 23.3E-2

      ELSE IF(IRCLA.EQ.5) THEN
C	  Elevation
        Ae  = 1.57E-8
        O1e = 23.3E-3
        O2e = 13.1E-2

c	  Alignment
        Aa  = 0.98E-8
        O1a = 32.8E-3
        O2a = 18.4E-2

c	  Cross
        Ac  = 0.98E-8
        O1c = 23.3E-3
        O2c = 13.1E-2

c	  Gauge
        Ag  = 0.98E-8
        O1g = 29.2E-3
        O2g = 23.3E-2

      ELSE IF(IRCLA.EQ.6) THEN
C	  Elevation
        Ae  = 0.98E-8
        O1e = 23.3E-3
        O2e = 13.1E-2

c	  Alignment
        Aa  = 0.59E-8
        O1a = 32.8E-3
        O2a = 18.4E-2

c	  Cross
        Ac  = 0.59E-8
        O1c = 23.3E-3
        O2c = 13.1E-2

c	  Gauge
        Ag  = 0.59E-8
        O1g = 29.2E-3
        O2g = 23.3E-2

      ENDIF

      DO 100 IOSIM=1,NOSIM
	
C	GENETATE RANDOM NUMBERS
      CALL RND(DSEED,NFFT,Re)
      CALL RND(DSEED,NFFT,Ra)
      CALL RND(DSEED,NFFT,Rc)
      CALL RND(DSEED,NFFT,Rg)

c      CALL RND(5,NFFT,Re)
c      CALL RND(5,NFFT,Ra)
c      CALL RND(5,NFFT,Rc)
c      CALL RND(5,NFFT,Rg)
   
C	GENERATE IRREGULARITY PROFILES FOR ALL SIMULATIONS   	
      DO 10 I=1,NFFT
        X=POSI+DT*SPEED*(I-1)
        GAM=0.0D0
        DO 10 J=1,NFFT
          GAM=GAML+DGAM*(DFLOAT(J)-1./2.)
C          PHASE=RR(J)*PI2

C		FRA PSD
C	    Elevation:
          PSDe=Ae*(O2e**2)*(GAM**2+O1e**2)/((GAM**4)*(GAM**2+O2e**2)) !PSDe 7Nov07
          CR=DSQRT(4.*PSDe*DGAM)										!PSDe 7Nov07
		PHASE=Re(J)*PI2
          GEOe(I)=GEOe(I)+CR*DCOS(PI2*GAM*X-PHASE)

C	    Alignment:
          PSDa=Aa*(O2a**2)*(GAM**2+O1a**2)/((GAM**4)*(GAM**2+O2a**2)) !PSDa 7Nov07
          CR=DSQRT(4.*PSDa*DGAM)										!PSDa 7Nov07
		PHASE=Ra(J)*PI2
          GEOa(I)=GEOa(I)+CR*DCOS(PI2*GAM*X-PHASE)

C	    Cross:
          PSDc=Ac*(O2c**2)/((GAM**2+O1c**2)*(GAM**2+O2c**2))			!PSDc 7Nov07
          CR=DSQRT(4.*PSDc*DGAM)										!PSDc 7Nov07
		PHASE=Rc(J)*PI2
          GEOc(I)=GEOc(I)+CR*DCOS(PI2*GAM*X-PHASE)

C	    Gauge:
          PSDg=Ag*(O2g**2)/((GAM**2+O1g**2)*(GAM**2+O2g**2))			!PSDg 7Nov07
          CR=DSQRT(4.*PSDg*DGAM)										!PSDg 7Nov07
		PHASE=Rg(J)*PI2
          GEOg(I)=GEOg(I)+CR*DCOS(PI2*GAM*X-PHASE)

C		Print values of PSD functions for each irre. type (at only I=1) (added 7Nov07)
		IF(I.EQ.1)WRITE(100,1400) GAM,PSDe,PSDa,PSDc,PSDg
C		WRITE(100,1400) GAM,PSDe,PSDa,PSDc,PSDg

   10 CONTINUE
	
  100 CONTINUE

C	FIND THE AVERAGE IRREGULARITY PROFILES
	RSIM = DFLOAT(NOSIM)
      DO 20 I=1,NFFT
        GEOx(I)=POSI+DT*SPEED*(I-1)	!x-ordinate array for rail profiles
        GEOe(I)=GEOe(I)/RSIM
        GEOa(I)=GEOa(I)/RSIM
        GEOc(I)=GEOc(I)/RSIM
        GEOg(I)=GEOg(I)/RSIM

C	  GEOc(I)=GEOc(I)/(0.5*GAUGE-GEOg(I))  !Rotational irregularities
C	  Previous line moved down 12Jul07
   20 CONTINUE

c	Next added 07Nov07
	WRITE(100,1150)
C	--------------------------------------------------------------------
C	Asign maximum tolerable deviations for each FRA class (4-6), unit: m
C	Class (1-3) not available, Y.B.Yang's book, p 359
C	added 12July07
c	---------------------------------------------------------------------
	SELECT CASE(IRCLA)
	  CASE(1,2,3)
		WRITE(ITO,*)'WARNING: Max. deviations of FRA classes 1-3 not 
	1available. Rail irregularities generated are not normalized'
          WRITE(10,*)'WARNING: Max. deviations of FRA classes 1-3 not 
	1available. Rail irregularities generated are not normalized'
		GOTO 500
	  CASE(4)
	  	Te = 0.00405
		Ta = 0.00510
		Tc = 0.00150
		Tg = 0.00180
	  CASE(5)
	  	Te = 0.00338
		Ta = 0.00425
		Tc = 0.00125
		Tg = 0.00165			
	  CASE(6)
	  	Te = 0.00270
		Ta = 0.00340
		Tc = 0.00100
		Tg = 0.00130
	END SELECT

C	NORMALIZE THE IRREGULARITY PROFILES TO FRA MAX. TOLERABLE DEVIATIONS
C	Find max. absolute value of each profile
	CALL FMAXA(GEOe,NFFT,EMAX)
	CALL FMAXA(GEOa,NFFT,AMAX)
	CALL FMAXA(GEOc,NFFT,CMAX)
	CALL FMAXA(GEOg,NFFT,GMAX)

C	Normalize to FRA deviations
	DO 200 I=1,NFFT
		GEOe(I)=GEOe(I)*(Te/EMAX)
		GEOa(I)=GEOa(I)*(Ta/AMAX)
		GEOc(I)=GEOc(I)*(Tc/CMAX)
		GEOg(I)=GEOg(I)*(Tg/GMAX)
		GEOc(I)=GEOc(I)/(0.5*GAUGE-GEOg(I))  !Rotational irregularities
  200	CONTINUE

  500	DO 300 I=1,NFFT
  300 WRITE(100,1200) GEOx(I),GEOe(I),GEOa(I),GEOc(I)

c 1000	FORMAT (//1H#,15X,38(1H+)/
c     + 1H#,15X,36H+    3D RAIL IRREGULARITY GEOMETRY     1X,1H+/
c     + 1H#,15X,38(1H+)/)
C	Changed next 07Nov07
 1000	FORMAT (//1H#,15X,38(1H+)/
     + 1H#,15X,36H+   3D RAIL IRREGULARITY GENERATION    1X,1H+/
     + 1H#,15X,38(1H+)/
     +5X,'GENERATING METHOD: BY RANDOM PROCESS',/
     +5X,'PSD FUNCTIONS: by Federal Railway Administration (FRA),USA',/
     +5X,'RAIL IRREGULARITY CLASS (by FRA) . . . . . = ',I3/
     +5X,'Number of spectrums for total rail length. = ',I10/
     +5X,'Seed value to generate random numbers . . .= ',F10.4/
     +5X,'Shortest wave length considered (m). . . . = ',F10.4/
     +5X,'Longest wave length considered (m) . . . . = ',F10.4/
     +5X,'Number of simulations . . . . . . . . . . .= ',I10/
     +5X,'Nominal (design) gauge of rails (m)        = ',F10.4/)

c 1100 FORMAT(5X,'GENERATING METHOD: BY RANDOM PROCESS',/
c     +5X,'PSD FUNCTIONS: by Federal Railway Administration (FRA),USA',/
c     +5X,'RAIL IRREGULARITY CLASS (by FRA) . . . . . = ',I3//
c     +5X,'Number of spectrums for total rail length. = ',I10/
c     +5X,'Seed value to generate random numbers . . .= ',F10.4/
c     +5X,'Upper limit of wave number (1/m). . . . . .= ',F10.7/
c     +5X,'Lower limit of wave number (1/m). . . . . .= ',F10.7/
c     +5X,'Number of simulations . . . . . . . . . . .= ',I10/
c     +5X,'Nominal (design) gauge of rails (m)        = ',F10.4//
c     +4X,'X-RAIL (m)',3X,'ELEVATION(m)',3X,'ALIGNMENT (m)',3X,'
c     +ROTATION (Rad)')
C	Changed 07Nov07

 1100	FORMAT (/1H#,15X,'POWER SPECTRAL DENSITY FUNCTIONS'/
     +1H#,3X,70(1H-)/
     +4X,'FREQ (Hz)',4X,'PSD-elev(m**3)',2X,'PSD-alig(m**3)',2X,'PSD-cros
     +(m**3)',2X,'PSD-gaus(m**3)')

 1150	FORMAT (/1H#,15X,'3D RAIL IRREGULARITY GEOMETRY'/
     +1H#,3X,60(1H-)/
     +4X,'X-RAIL (m)',3X,'ELEVATION (m)',3X,'ALIGNMENT (m)',3X,'
     +ROTATION (Rad)')

 1200 FORMAT(2X,F10.5,1X,3E15.5)
 1400 FORMAT(2X,F10.5,1X,4E15.5) !Added 07Nov07 for PSD print
     
	RETURN
      END
C
C	=======================================================================
      SUBROUTINE RND(DSEED,NR,RR)
C
C	PROGRAM FOR RANDOM NUMBER GENERATION
C	------------------------------------
      INTEGER NR
      REAL*8 RR(NR)
      DOUBLE PRECISION DSEED
C	--------------------------------------
C     SPECIFICATIONS FOR LOCAL VARIABLES
C	--------------------------------------
C     INTEGER 1
      DOUBLE PRECISION D2P31M, D2P31
C	---------------------------------------------
C     D2P31M=(2**31)-1
C     D2P31=(2**31) (OR AN ADJUSTED VALUE)
C	---------------------------------------------
      DATA D2P31M/2147483647.D0/
      DATA D2P31/2147483648.D0/

C	FIRST EXECUTABLE STATEMENT
      DO 10 I=1,NR
      DSEED=DMOD(16807.D0*DSEED,D2P31M)
   10 RR(I)=DSEED/D2P31

      RETURN
      END
C
C	=======================================================================
      SUBROUTINE FMAXA(A,N,AMAX) 

      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	Program to find maximum of absolute values of an array A(N)	
c	July2007
	DIMENSION A(*)

	AMAX = 0.0
	DO I=1,N
		IF(ABS(A(I)).GT.AMAX) AMAX = ABS(A(I))
	ENDDO		

	RETURN
	END
C	=======================================================================
      SUBROUTINE CALIRR(GEOx,GEOe,GEOa,GEOr,WHPOS,JTE,NWHEEL,TIRAIL)		!Added JTE 11Feb07
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------
C     PROGRAM TO GENERATE 3D RAIL/TRACK IRREGULARITIES:
C		+ Interpolate track profiles at wheel position (TIRLAT,TIRROT,TIRVER) 
C		  from profiles generated at rail sections
C	-----------------------------------------------------
C	INPUT
C	POSI =		Wheel position
C	GEOx(NFFT): 
C	GEOe(NFFT)
C	GEOa(NFFT)
C	GEOr(NFFT)
C	JTE	 = Iteration number, if ITE>1, not calculate but use previous at JTE=1

C	OUPUT:
C	TIRLAT = Lateral (or Alignment) Irregularity
C	TIRROT = Vertical (or Elevation) Irregularity
C	TIRVER = Resulted Rotational Irregularities

C	INDEX:
C	IRIN   = GENERATING METHOD OF IRREGULARITY 
C		    =0: No Irregularity considered
C			=1: Text file (Measurement data)
C			=2: Generate by semi-periodic functions
C			=3: Generate by random functions, where:
C				IRAN = 1: by FRA Power Spectral Density (Federal Railway Administration, USA)
C				IRAN = 2: by SBCF Power Spectral Density (French National Railway)
C			=4: Rail joints
C	NFFT  = NUMBER OF FAST FOURIER TRANSFORM USED IN SUBROUTINE GEOMET WHICH IS
C        TO CALCULATE TRACK GEOMETRY(IF NOT SPECIFIED, NFFT=1024)
c	------------------------------------------------------ 
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM
	COMMON /RIR3/ RLENG,SLSTEP,DVMIN,DHMIN,GAUGE,DSEED,GAMU,GAML
	COMMON /IR3O/ TIRLAT,TIRROT,TIRVER
C	
      DIMENSION GEOx(NFFT),GEOe(NFFT),GEOa(NFFT),GEOr(NFFT)
      
	DIMENSION WHPOS(NWHEEL),TIRAIL(3,NWHEEL)
c	------------------------------------------------------
	IF(JTE.GT.1) GOTO 100  !not calculate but use previous IRRE at JTE=1, Added JTE 11Feb07
	
	DO 800 IWH = 1,NWHEEL
	
	POSI = WHPOS(IWH)
	
      IF(IRIN.EQ.0) GOTO 35  !No Irregularity considered
      IF(IRIN.EQ.3) GOTO 25

C	THE FOLLOWING STATEMENTS ARE TO CALCULATE GEOMETRY BASED ON PSD
   25 CONTINUE

      DO 10 I=1,NFFT-1  !Searching location of wheel along the rail
        IF(POSI.GE.GEOx(I).AND.POSI.LT.GEOx(I+1)) THEN

		 X1 = (POSI-GEOx(I))/(GEOx(I+1)-GEOx(I)) !added 0.1 temperarily 28Feb07
c		 TIRLAT = 0.1*(GEOa(I) + X1*(GEOa(I+1)-GEOa(I)))
c		 TIRROT = 0.1*(GEOr(I) + X1*(GEOr(I+1)-GEOr(I)))
c		 TIRVER = 0.1*(GEOe(I) + X1*(GEOe(I+1)-GEOe(I)))
c												!Removed 0.1 from previous 24June07
		 TIRLAT = GEOa(I) + X1*(GEOa(I+1)-GEOa(I))
		 TIRROT = GEOr(I) + X1*(GEOr(I+1)-GEOr(I))
		 TIRVER = GEOe(I) + X1*(GEOe(I+1)-GEOe(I))

           GOTO 45
        ENDIF
   10 CONTINUE

      TIRLAT=0.0D0
      TIRROT=0.0D0
      TIRVER=0.0D0

   45 GOTO 50

   35 CONTINUE
      TIRLAT=0.0D0
      TIRROT=0.0D0
      TIRVER=0.0D0

   50 CONTINUE
   
      TIRAIL(1,IWH) = TIRVER
      TIRAIL(2,IWH) = TIRROT
      TIRAIL(3,IWH) = TIRLAT
800   CONTINUE

      
  100 RETURN
      END
C
C	=======================================================================
      SUBROUTINE GENIRRE3D(WHPOS,TLENGTH,JTE,NWHEEL,TIRAIL) 

      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------
C     TO GENERATE 3D RAIL/TRACK IRREGULARITIES USING:
C		+ Summation of periodic profiles with large frequency range
C		+ Contain high occurent frequencies
C		+ Account for sleeper step, rail joint distance, and deriorated spots
C		+ Possible to account for all 4 type of rail irregularities 
C		  (elevation, alignment, cross, and gauge)
C	Programmed 6Sep2006 by NguyenDV
C	-----------------------------------------------------
C	INPUT:
C	POSI	= Position of current wheel
C	TLENGTH = Total length of path where wheel runs
C	JTE		= Iteration number, if ITE>1, not calculate but use previous at JTE=1

C	OUPUT:
C	TIRLAT = Lateral (or Alignment) Irregularity
C	TIRROT = Vertical (or Elevation) Irregularity
C	TIRVER = Resulted Rotational Irregularities
C	-----------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM
	COMMON /RIR3/ RLENG,SLSTEP,DVMIN,DHMIN,GAUGE,DSEED,GAMU,GAML
	COMMON /IR3O/ TIRLAT,TIRROT,TIRVER						  !Output block
	
	DIMENSION WHPOS(NWHEEL),TIRAIL(3,NWHEEL)
C	-----------------------------------------------------
	IF(JTE.GT.1) GOTO 900

      DO 800 IWH = 1,NWHEEL
      
      POSI = WHPOS(IWH)
      
      PI=3.141592654
	IF(IRCLA.EQ.0) THEN 
		TIRLAT = 0.
		TIRROT = 0.	
		TIRVER = 0.		
		GOTO 900
	ENDIF

C	Assign number of Irrgularity wave length samples (later, to be inputed)
	NV1=1
	NV2=1
	NH1=1
	NH2=1
	NC1=1
	NC2=1

C	Asign maximum tolerable deviations for each FRA class
c	-----------------------------------------------------
	SELECT CASE(IRCLA)
	  CASE(1,2,3)
		GOTO 900

	  CASE(4)
	  	Av = 0.00405
		Ah = 0.00510
		Ac = 0.00150
		Ag = 0.00180

	  CASE(5)
	  	Av = 0.00338
		Ah = 0.00425
		Ac = 0.00125
		Ag = 0.00165			

	  CASE(6)
	  	Av = 0.00270
		Ah = 0.00340
		Ac = 0.00100
		Ag = 0.00130

	END SELECT	
					   	
C	Generate Frequencies (wave) profiles and Normalized to amplitude of FRA deviations
c	----------------------------------------------------------------------------------
	XSL = PI*POSI/SLSTEP
	XDV = PI*POSI/DVMIN
	XHV = PI*POSI/DHMIN

	STL = DSIN(PI*POSI/TLENGTH)
	SRL = DSIN(PI*POSI/RLENG)
	SSL = DSIN(XSL)
	SDV = DSIN(XDV)
	SDH = DSIN(XHV)

c	Vertical (or Elevation) Irregularity:
	NV = NV1+NV2+2     
	SV1 = 0.
	IF(NV1.EQ.0) GOTO 100		!Not include effect of sleeper step
	DO I=1,NV1
	  XI = DFLOAT(I)*XSL
	  SV1 = SV1 + DSIN(XI)
	ENDDO

  100	SV2 = 0.
  	IF(NV2.EQ.0) GOTO 150
	DO I=1,NV2
	  XI = DFLOAT(I)*XDV		!Not include effect of vertical spots
	  SV2 = SV2 + DSIN(XI)
	ENDDO

  150	TIRVER = Av*(STL + SRL + SV1 + SV2)/NV
     
c	Lateral (or Alignment) Irregularity:
	NH = NH1+NH2+2     
	SH1 = 0.
	IF(NH1.EQ.0) GOTO 200		!Not include effect of sleeper step
	DO I=1,NH1
	  XI = DFLOAT(I)*XSL
	  SH1 = SH1 + DSIN(XI)
	ENDDO

  200	SH2 = 0.
  	IF(NH2.EQ.0) GOTO 250		!Not include effect of lateral spots
	DO I=1,NH2
	  XI = DFLOAT(I)*XHV
	  SH2 = SH2 + DSIN(XI)
	ENDDO

  250	TIRLAT = Ah*(STL + SRL + SH1 + SH2)/NH

C	Cross irregularities (Vertical direction)
	NC = NC1+NC2+2     
	SC1 = 0.
	IF(NC1.EQ.0) GOTO 300		!Not include effect of sleeper step
	DO I=1,NC1
	  XI = DFLOAT(I)*XSL
	  SC1 = SC1 + DSIN(XI)
	ENDDO

  300	SC2 = 0.
	IF(NC2.EQ.0) GOTO 350		!Not include effect of vertical spots
	DO I=1,NC2
	  XI = DFLOAT(I)*XDV
	  SC2 = SC2 + DSIN(XI)
	ENDDO
  350	TIRC = Ac*(STL + SRL + SC1 + SC2)/NC

C	Gauge irregularities (Lateral direction)
	TIRG = Ag*(STL + SRL + SDV + SDH)/4

C	Resultant Rotational Irregularities:
	TIRROT = TIRC/(GAUGE - TIRG)	

      TIRAIL(1,IWH) = TIRVER
      TIRAIL(2,IWH) = TIRROT
      TIRAIL(3,IWH) = TIRLAT
800   CONTINUE
	    
  900 RETURN
      END

C	=======================================================================
      SUBROUTINE GENIRRE3D1(POSI,TLENGTH,JTE) !Currently not use

      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------
C     TO GENERATE 3D RAIL/TRACK IRREGULARITIES USING:
C		+ Summation of periodic profiles with large frequency range
C		+ Contain high occurent frequencies
C		+ Account for sleeper step, rail joint distance, and deriorated spots
C		+ Possible to account for all 4 type of rail irregularities 
C		  (elevation, alignment, cross, and gauge)
C	Programmed 6Sep2006 by NguyenDV
C	-----------------------------------------------------
C	INPUT:
C	POSI	= Position of current wheel
C	TLENGTH = Total length of path where wheel runs

C	OUPUT:
C	TIRLAT = Lateral (or Alignment) Irregularity
C	TIRROT = Vertical (or Elevation) Irregularity
C	TIRVER = Resulted Rotational Irregularities
C	-----------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM
	COMMON /RIR3/ RLENG,SLSTEP,DVMIN,DHMIN,GAUGE,DSEED,GAMU,GAML
	COMMON /IR3O/ TIRLAT,TIRROT,TIRVER						  !Output block
C	-----------------------------------------------------
	IF(JTE.GT.1) GOTO 900

	WRITE(ISO,*)POSI

      PI=3.141592654
	IF(IRCLA.EQ.0) THEN 
		TIRLAT = 0.
		TIRROT = 0.	
		TIRVER = 0.		
		GOTO 900
	ENDIF

C	Asign maximum tolerable deviations for each FRA class
c	-----------------------------------------------------
	SELECT CASE(IRCLA)
	  CASE(1,2,3)
		GOTO 900

	  CASE(4)
	  	Av = 0.00405
		Ah = 0.00510
		Ac = 0.00150
		Ag = 0.00180

	  CASE(5)
	  	Av = 0.00338
		Ah = 0.00425
		Ac = 0.00125
		Ag = 0.00165			

	  CASE(6)
	  	Av = 0.00270
		Ah = 0.00340
		Ac = 0.00100
		Ag = 0.00130

	END SELECT					   	

C	Calculate  number of sample/frequencies contained
c	-------------------------------------------------
c	Number of sample less than sleeper step & larger than spot distance in vertical
	NV1 = IDNINT(SLSTEP/DVMIN)

c	Number of sample less than sleeper step & larger than spot distance in lateral
	NH1 = IDNINT(SLSTEP/DVMIN)

c	Number of rail samples which is larger than sleeper step	
	N2 = IDNINT(RLENG/SLSTEP)

c	Number of sample less than total length & larger than rail bar length
	N3 = IDNINT(TLENGTH/RLENG)
C	
C	Generate Frequencies (wave) profiles
c	--------------------------------------
	SV1 = 0.
	SV2 = 0.
	SV3 = 0.

	STL = DSIN(PI*POSI/TLENGTH)
	SRL = DSIN(PI*POSI/RLENG)
	SSL = DSIN(PI*POSI/SLSTEP)
	SDV = DSIN(PI*POSI/DVMIN)
	SDH = DSIN(PI*POSI/DHMIN)

c	Vertical (or Elevation) Irregularity:
	NV = NV1+N2+2
	X1 = PI*POSI/DVMIN
	DO I=1,NV1
	  XI = DFLOAT(I)*X1
	  SV1 = SV1 + DSIN(XI)
	ENDDO

	X2 = PI*POSI/SLSTEP
	DO I=1,N2
	  XI = DFLOAT(I)*X2
	  SV2 = SV2 + DSIN(XI)
	ENDDO	
c	TIRVER=Av*(DSIN(PI*POSI/TLENGTH)
c	1	     + DSIN(PI*POSI/RLENG)
c     2		 + SV1 + SV2)/NV 
C	TIRVER = Av*(DSIN(PI*POSI/TLENGTH) + DSIN(PI*POSI/RLENG) 
C	1	       + DSIN(PI*POSI/SLSTEP)  + DSIN(PI*POSI/DVMIN))/4  
	TIRVER = Av*(DSIN(PI*POSI/SLSTEP)  + DSIN(PI*POSI/DVMIN))/2
     
c	Lateral (or Alignment) Irregularity:
	NH = N2+N3+2     
	X3 = PI*POSI/RLENG
	DO I=1,N3
	  XI = DFLOAT(I)*X3
	  SV3 = SV3 + DSIN(XI)
	ENDDO		
c	TIRLAT = Ah*(DSIN(PI*POSI/TLENGTH)+SV3+SV2+DSIN(PI*POSI/DHMIN))/NH
	TIRLAT = Ah*(DSIN(PI*POSI/TLENGTH) + DSIN(PI*POSI/RLENG) 
	1	       + DSIN(PI*POSI/SLSTEP)  + DSIN(PI*POSI/DHMIN))/4

C	Cross irregularities (Vertical direction)
	TIRC = Ac*(DSIN(PI*POSI/TLENGTH) + DSIN(PI*POSI/RLENG) 
	1	       + DSIN(PI*POSI/SLSTEP)  + DSIN(PI*POSI/DVMIN))/4 

C	Gauge irregularities (Lateraldirection)
	TIRG = Ag*(DSIN(PI*POSI/TLENGTH) + DSIN(PI*POSI/RLENG) 
	1	     + DSIN(PI*POSI/SLSTEP)  + DSIN(PI*POSI/DHMIN))/4

C	Resulted Rotational Irregularities:
	TIRROT = TIRC/(GAUGE - TIRG)	
	    
  900 RETURN
      END

C	=======================================================================
      SUBROUTINE GENIRRE3D2(POSI,TLENGTH) !Currently not use

      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------
C     TO GENERATE 3D RAIL/TRACK IRREGULARITIES USING:
C		+ Summation of periodic profiles with large frequency range
C		+ Contain high occurent frequencies
C		+ Account for sleeper step, rail joint distance, and deriorated spots
C		+ Possible to account for all 4 type of rail irregularities 
C		  (elevation, alignment, cross, and gauge)
C	Programmed 6Sep2006 by NguyenDV
C	-----------------------------------------------------
C	INPUT:
C	POSI	= Position of current wheel
C	TLENGTH = Total length of path where wheel runs

C	OUPUT:
C	TIRLAT = Lateral (or Alignment) Irregularity
C	TIRROT = Vertical (or Elevation) Irregularity
C	TIRVER = Resulted Rotational Irregularities
C	-----------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM
	COMMON /RIR3/ RLENG,SLSTEP,DVMIN,DHMIN,GAUGE,DSEED,GAMU,GAML
	COMMON /IR3O/ TIRLAT,TIRROT,TIRVER						  !Output block
C	-----------------------------------------------------
      PI=3.141592654

	IF(IRCLA.EQ.0) THEN 
		TIRLAT = 0.
		TIRROT = 0.	
		TIRVER = 0.		
		GOTO 900
	ENDIF

C	Asign maximum tolerable deviations for each FRA class
c	-----------------------------------------------------
	SELECT CASE(IRCLA)
	  CASE(1,2,3)
		GOTO 900

	  CASE(4)
	  	Av = 0.00405
		Ah = 0.00510
		Ac = 0.00150
		Ag = 0.00180

	  CASE(5)
	  	Av = 0.00338
		Ah = 0.00425
		Ac = 0.00125
		Ag = 0.00165			

	  CASE(6)
	  	Av = 0.00270
		Ah = 0.00340
		Ac = 0.00100
		Ag = 0.00130

	END SELECT					   	

C	Calculate  number of sample/frequencies contained
c	-------------------------------------------------
c	Number of sample larger than sleeper step in vertical		
	NV = IDNINT(TLENGTH/SLSTEP)
c	Number of sample less than sleeper step & larger than spot distance in vertical	 
	MV = IDNINT(RLENG/DVMIN)

c	Number of sample larger than sleeper step in lateral
	NH = IDNINT(RLENG/SLSTEP)

c	Number of sample less than sleeper step & larger than spot distance in lateral
	MH = IDNINT(RLENG/DHMIN)
C	
C	Generate Frequencies (wave) profiles
c	--------------------------------------
	SNV = 0.
	SMV = 0.
	SNH = 0.
	SMH = 0.

	X0 = PI*POSI/SLSTEP

	DO I=1,NV
	  XI = DFLOAT(I)*X0
	  SNV = SNV + DSIN(XI) + DCOS(XI)
	ENDDO
	  SNV = SNV/NV			!Normalized to 1

	DO I=2,MV
	  YI = X0/DFLOAT(I)
	  SMV = SMV + DSIN(YI) + DCOS(YI)
	ENDDO
	  SMV = SMV/MV			!Normalized to 1

	DO I=1,NH
	  XI = DFLOAT(I)*X0
	  SNH = SNH + DSIN(XI) + DCOS(XI)
	ENDDO
	  SNH = SNH/NH			!Normalized

	DO I=2,MH
	  YI = X0/DFLOAT(I)
	  SMH = SMH + DSIN(YI) + DCOS(YI)
	ENDDO
	  SMH = SMH/MH			!Normalized
C

C	Generate Irregularities:
c	------------------------
c	Lateral (or Alignment) Irregularity:		
	TIRLAT = Ah*(SNH + SMH)

c	Vertical (or Elevation) Irregularity:	
	TIRVER = Av*(SNV + SMV)

C	Cross irregularities (Vertical direction)
	TIRC   = Ac*(SNV + SMV)

C	Gauge irregularities (Lateraldirection)
	TIRG   = Ag*(SNH + SMH)

C	Resulted Rotational Irregularities:
	TIRROT = TIRC/(GAUGE - TIRG)
     
  900 RETURN
      END

C	==========================================================================
      SUBROUTINE RDAMCO1(AK,AM,IDIAG,NEQ,NWM,NWK)     
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------------------------------
C     PROGRAM TO CALCULATE RELAIGH DAMPING COEFFICIENT USING MODAL DAMPING RATIOS:
C	Temperarily removed, later activate to solve for part of structure (more advanced)
C	---------------------------------------------------------------------------
C	INPUT:
C	AM(NWM)  = MASS MATRIX
C	AK(NWK)  = STIFFNESS MATRIX
C	NWK      = TOTAL NUMBER OF ELEMENTS OF STIFFNESS MATRIX
C	NWM      = TOTAL NUMBER OF ELEMENTS OF MASS MATRIX

C	OUTPUT:
C	RDM		 = RAYLEIGH DAMPING COEFFCIENT PROPORTIONAL TO MASS (ALPHA)
C	RDK		 = RAYLEIGH DAMPING COEFFCIENT PROPORTIONAL TO STIFFNESS (BETA)

C	INDEX:
C	ISOLV     = EIGEN VALUE SOLVER 
C				(1-SUBSPACE ITERATION,2-LANCZOS,3-INVERSE ITERATION)
C     IEIG      = TYPE OF EIGENVALUE ANALYSIS
C				(0-NO EIGEN SOLUTION,1-STANDARD EIGEN SOLUTION
C				 2-BUCKLING EIGEN SOLUTION,	3-NATUTAL FREQUANCEY)
C     IMASS     = OPTIONS FOR FORMATION OF MASS MATRIX
C				(0-NO MASS MATRIX,     1-LUMPED MASS MATRIX
C				 2-CONSISTETLY LUMPED, 3-CONSISTENT MASS MATRIX)
C     NSEIG     = NUMBER OF STEPS BETWEEN EIGENVALUE ANALYSIS
C     NROOT     = NUMBER OF EIGENVALUES REQUIRED
C	IVPRT     = NUMBER OF EIGENVECTOR REQUIRED
C     NC        = NUMBER OF ITERATION VECTORS USED		= NCROOT
C     NNC       = SIZE OF MATRICES (UPPER TRIANG.)		= NNCROOT
C     NITEM     = NUMBER OF SUBSPACE ITERATIONS PERMITTED
C				IF NITEM.EQ.-1, TWO EIGENVALUES SHOULD BE INPUTTED
C     IFSS      = FLAG FOR STURM SEQUENCE CHECK
C     EPS       = CONVERGENCE TOLERANCE ON EIGENVALUES    = RTOL
C	-------------------------------------------------------------------------
      COMMON /IOUNIT/IN,IOUT,NSTIF
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
C
C      DIMENSION AM(*),AK(*),IDIAG(*),ROOT(NROOT) 
      DIMENSION AM(NWM),AK(NWK),IDIAG(NEQ+1),ROOT(NROOT)	  
      DIMENSION R(NEQ,1),TT(1),W(1),EIGV(3),D(1),VEC(NC,1),AR(1),BR(1)
      DIMENSION RTOLV(1),BUP(1),BLO(1),BUPC(1) 	 
C	--------------------------------------------------------------------------
      IF(NROOT.EQ.0) NROOT=2
      IF(NITEM.EQ.0) NITEM=100
      IFSS1=1
      IFPR1=0
	NEQ1=NEQ+1

      IF(NITEM.GT.0) THEN
	  WRITE(ITO,*)'SOLVE FOR THE SMALLEST EIGENVALUES BY SUBSPACE'
        WRITE(10,*)'SOLVE FOR THE SMALLEST EIGENVALUES BY SUBSPACE'
	  WRITE(ISO,135)
	  NSTIF = 2
	  CALL SSPACE1(AK,AM,IDIAG,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,
     1			   BLO,BUPC,NEQ,NEQ1,NWK,NWM,NROOT,EPS,NC,NNC,NITEM,
     2			   IFSS1,IFPR1,NSTIF,ITO,ISO)

      ELSEIF(NITEM.EQ.-1) THEN
	  READ(ITI,*)
        READ(ITI,1000) (ROOT(I),I=1,NROOT)
        WRITE(ISO,2000) (ROOT(I),I=1,NROOT)
      ENDIF

      PI=3.141592654
      DO 40 I=1,NROOT
        TEMP=DSQRT(ROOT(I))
        TEMQ=TEMP/(2.*PI)
        TEMR=1./TEMQ
        WRITE(ISO,2020) I,TEMR,TEMQ,TEMP,ROOT(I)
   40 CONTINUE

      DO 20 I=1,NROOT
        ROOT(I)=DSQRT(ROOT(I))
   20 CONTINUE

      TEMP=2.*(ROOT(1)*ROOT(2))/(ROOT(2)**2-ROOT(1)**2)
      RDM =TEMP*(ROOT(2)*ZET1-ROOT(1)*ZET2)
      RDK =TEMP*(-ZET1/ROOT(2)+ZET2/ROOT(1))

 1000 FORMAT(2F10.0)
  135	FORMAT (//1H#,18X,49(1H*)/1H#,18X,1H*,47X,1H*/
     +1H#,18X,1H*,1X,45H SOLVE THE EIGENVALUES OF BRIDGE BY SUSBSPACE,
     +1X,1H*/1H#,18X,1H*,47X,1H*/1H#,18X,49(1H*))

 2000 FORMAT(//,1X,'EIGENVALUES(SQUARED NATURAL FREQUENCIES)',/
     &1X,'FIRST EIGENVALUE    = ',E10.3/
     &1X,'SECOND EIGENVALUE   = ',E10.3)

      WRITE(ISO,2010)
 2010 FORMAT(//,' MODE   PERIOD         FREQUENCY      FREQUENCY',
     &  '      EIGENVALUE',/
     &  '        (TIME)         (CYC/TIME)     (RAD/TIME)     ',
     &  '(RAD/TIME)**2')
 2020 FORMAT(I5,4(2X,E13.6))
       
	RETURN
      END
C
C=====================================================================   
C      SUBROUTINE TRGPARA3D(IT1,IT2,IT3)
      SUBROUTINE TRGPARA3D

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     TO CALCULATE TRAIN GENERAL PARAMETERS (NUMBERS OF COMPONENTS)
C	-------------------------------------------------------------
C	INPUT:
C	------
C	IT1	  = NUMBER OF FRONT POWER CARS
C	IT2	  = NUMBER OF INTERMEDIATE & PASSENGER CARS
C	IT3	  = NUMBER OF REAR POWER CARS
C	-------
C	OUTPUT:
C	-------
C	NACAB = Number of articulated car bodies = IT2
C	NCARB = Total number of car bodies = IT1+IT2+IT3
C	NBOGI = Total number of bogies     = IT1*2+(IT2+1)+IT3*2
C	NWHEEL= Total number of wheels     = NBOGI*2

C	NEQTR  = Number of equations of train's stiffness/mass /damping matrix = NWMTR
C	NWKTR  = Number of elements below skyline of upper tria. stiff/damping matrix of train
C     -----------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL

	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH

	COMMON /ITRM/ ITRK	!Added here 25Mar06
C	-----------------------------------------------------
C	Previous lines changed to next for case of classical train (or only power cars)
	IF(IT2.GE.3)THEN	
		NACAB = IT2
	ELSE IF(IT2.EQ.0)THEN
		NACAB = 0
	ENDIF
c	The whole train:
c	----------------
	NCARB = IT1 + IT2 + IT3

c	Next IF block added 27Nov2005
	IF(IT2.GE.3)THEN	
	  NBOGI = IT1*2 + (IT2+1 ) + IT3*2
	  NEQTR = (5+2*(5+3+3))*IT1 + IT2*(2*3+5+5)+11 + IT3*(5+2*(5+3+3))

	  IF(ITRK.EQ.1.OR.ITRK.EQ.2)THEN					 !Added here 25Mar06	
	    NWKTR = IT1*167 + (83+(IT2-1)*(78+29)+72) + IT3*167   !Without couple terms of car bodies	  
	  ELSE IF(ITRK.EQ.3)THEN	
	    NWKTR = IT1*167 + (83+(IT2-1)*(138+29)+72) + IT3*167 
	  ENDIF

	ELSE IF(IT2.EQ.0)THEN
	  NBOGI = IT1*2 + IT3*2
	  NEQTR = (5+2*(5+3+3))*IT1 + (5+2*(5+3+3))*IT3	
	  NWKTR = IT1*167 + IT3*167
	ENDIF

	NWMTR = NEQTR	!Lumped mass matrix
	NWHEEL= NBOGI*2

	RETURN	
	END
C	
C	=======================================================================
      SUBROUTINE RETRN3D(ITRAIN,PCAB,CADI,ADDI,ALSU,PRSS,PBOG,PRPS,PRWH,
	1				   WRCO,WPOS,IDATM)    
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	IDATM added 30Oct04 by NguyenDV
C	-----------------------------------------------------------------------
C	PROGRAM TO READ & WRITE PROPERTIES OF TRAIN (CAR BODY,SUSPENSION,BOGIE)
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	IT1	  = NUMBER OF FRONT POWER CARS
C	IT2	  = NUMBER OF INTERMEDIATE & PASSENGER CARS
C	IT3	  = NUMBER OF REAR POWER CARS

C	NACAB = Number of articulated car bodies = IT2
C	NCARB = Total number of car bodies = IT1+IT2+IT3
C	NBOGI = Total number of bogies     = IT1*2+(IT2+1)+IT3*2
C	NWHEEL= Total number of wheels     = NBOGI*2

C	IDATM: Data type
C		  = 0:Input data by text format file (Original & Long)
C		  = 1:Input data by GiD (Optimized by Input only typical data)
C	-------
C	OUTPUT:
C	-------
C	ITRAIN(3)	   = Train composition stored in an array
C	PCAB(NCARB,4)  = Car body mechanical properties (MASS,Jx,Jy,Jz)		!updated 19Jan08
C	CADI(NCARB,7)  = Dimension parameters of cars (sA,sB,qA,qB,hC,h1A,h1B) 
C	ADDI(NACAB,5)  = Additional Dimension Parameters for Articulated Cars (b3,h5A,h5B,h6A,h6B)
C	ALSU(NACAB-1,6)= Longitudinal Suspensions of Articulated Cars  =(kTH,kTV,cX)
C	PRSS(NBOGI,9)  = Secondary Suspensions (k2H,c2H,k2V,c2V,k2X,c2X,c2RX,b2,h2) 
C	PBOG(NBOGI,4)  = Bogie of whole train (MASS,Jx,Jy,Jz)				!updated 19Jan08
C	PRPS(NWHEEL,9) = Primary Suspensions (k1H,c1H,k1V,c1V,k1X,c1X,b1,h3,t)
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jy,Jz,rW,b0)	!updated 19Jan08
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwY,cwY)
C	WPOS(NWHEEL)   = Wheel position			
C	-------------------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS	
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH

	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL
	COMMON /WRCT/ Ewr,WMU,RDR,RDW,CANG,SH,D11,D22,D23,D33,DHV,DV		!Created here 6Mar08
C
      DIMENSION ITRAIN(*),PCAB(NCARB,4),CADI(NCARB,7),ADDI(NACAB,5)
	DIMENSION ALSU(NACAB-1,6),PRSS(NBOGI,9),PBOG(NBOGI,4)
	DIMENSION PRPS(NWHEEL,9),PRWH(NWHEEL,6),WRCO(NWHEEL,6)
	DIMENSION WPOS(NWHEEL),NWHNUM(NWHEEL)
C	----------------------------------------------------------------------  
      ITRAIN(1)=IT1
      ITRAIN(2)=IT2
      ITRAIN(3)=IT3

C	CHECK NUMBER OF TRAIN COMPONENTS :Number of articulated cars must =0 or >=3
      IF(ITRAIN(2).GT.0.AND.ITRAIN(2).LT.3) THEN
        WRITE(*,9000) ITRAIN(2)
      ENDIF
      
      WRITE(ITO,1005)
      WRITE(10,1005)
      WRITE(ISO,900)
 
 1005 FORMAT(/,1X,'READ TRAIN PROPERTIES (CAR,SUSPENSION,BOGIE,WHEEL)')
  900	FORMAT (//1H#,18X,49(1H*)/1H#,18X,1H*,47X,1H*/
     +1H#,18X,1H*,1X,45HTRAIN PROPERTIES (CAR,SUSPENSION,BOGIE,WHEEL),
     +1X,1H*/1H#,18X,1H*,47X,1H*/1H#,18X,49(1H*))
C	-------------------------------------------------
C	A.FULL TEXT FORMAT DATA IS USED: 
C	The long data, where all components are listed 
C	-------------------------------------------------	
	IF (IDATM.EQ.0) THEN

C	Read Title for Train Properties
	READ(ITI,*)
C	---------------------
C	1a.CAR BODY PROPERTIES 
C	---------------------
	READ(ITI,*)
	WRITE(ISO,1100)
      WRITE(ITO,*) ('READ & WRITE CAR BODY PROPERTIES')
      WRITE(10,*) ('READ & WRITE CAR BODY PROPERTIES')

      DO 10 I=1,NCARB
        READ(ITI,1110) (PCAB(I,J),J=1,4)
        WRITE(ISO,1120) I,(PCAB(I,J),J=1,4)
   10 CONTINUE

 1100 FORMAT(//,1X,'PROPERTIES OF CAR BODY',/
     +1X,'Car No.      Mass     Rolling moment   Yawing moment  Nodding 
     +moment')
 1110 FORMAT(4F10.0)
 1120 FORMAT(I3,4X,4E15.6)
C	-----------------------------------
C	2a.DIMMENSION PARAMETERS OF CAR BODY:
C	-----------------------------------
	READ(ITI,*)
      WRITE(ISO,1210)

      DO 20 I=1,NCARB
        READ(ITI,*) (CADI(I,J),J=1,7)
        WRITE(ISO,1220) I,(CADI(I,J),J=1,7)
   20 CONTINUE

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
C	3a.ADDITIONAL DIMMENSION PARAMETERS OF ARTICULATED CAR BODY:
C	-----------------------------------------------------------
	IF(IT2.GE.3) THEN

	READ(ITI,*)
      WRITE(ISO,1310)

      DO 30 I=1,NACAB
        READ(ITI,*) (ADDI(I,J),J=1,5)
        WRITE(ISO,1320) I,(ADDI(I,J),J=1,5)
   30 CONTINUE

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
C	4a.LONGITUDINAL SUSPENSIONS OF ARTICULATED CARS:
C	----------------------------------------------
	READ(ITI,*)
      WRITE(ISO,1410)

      DO 40 I=1,NACAB-1
        READ(ITI,*) (ALSU(I,J),J=1,3)
        WRITE(ISO,1420) I,(ALSU(I,J),J=1,3)
	  DO 40 J=4,6
		ALSU(I,J) = 0.0  ! Assign the nominal zero values
   40 CONTINUE

 1410 FORMAT(//,1X,'LONGITUDINAL SUSPENSIONS OF ARTICULATED CAR BODY ',/
     +1X,'Transverse stiffness of central hinge . . . . . kTH '/
     +1X,'Vertical stiffness of central hinge. . . .  . . kTV '/
     +1X,'Damping coefficient of longitudinal damper. . . .cX '/
     +1X,'ACarNo.   kTH	     kTV		cX')

 1420 FORMAT(I3,4X,3E15.4)

	ENDIF	! <---IF(IT2.GE.3)
C	-----------------------
C	5a.SECONDARY SUSPENSIONS:
C	-----------------------
	READ(ITI,*)
      WRITE(ISO,1500)

      DO 50 I=1,NBOGI
        READ(ITI,1510) (PRSS(I,J),J=1,9)
        WRITE(ISO,1520) I,(PRSS(I,J),J=1,9)
   50 CONTINUE

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
C	----------------------
C	6a.PROPERTIES OF BOGIES:
C	----------------------  
      READ(ITI,*)
	WRITE(ISO,1600)

      DO 60 I=1,NBOGI
        READ(ITI,1610) (PBOG(I,J),J=1,4)
        WRITE(ISO,1620) I,(PBOG(I,J),J=1,4)
   60 CONTINUE

 1600 FORMAT(//,1X,'PROPERTIES OF BOGIES',/
     +1X,'BogieNo.',2x,'Mass',10x,'RollingMoment',2x,'YawingMoment',3x,'
     +NoddingMoment')
 1610 FORMAT(4F10.0)
 1620 FORMAT(I3,4X,4E15.6)

C	----------------------
C	7a.PRIMARY SUSPENSIONS:
C	----------------------
	READ(ITI,*)
      WRITE(ISO,1700)

      DO 70 I=1,NWHEEL
        READ(ITI,1710) (PRPS(I,J),J=1,9)
        WRITE(ISO,1720)I,(PRPS(I,J),J=1,9)
   70 CONTINUE

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
C	8a.WHEEL AXLE SET PROPERTIES: 
C	---------------------------
      READ(ITI,*)
	WRITE(ISO,1800)

      DO 80 I=1,NWHEEL
        READ(ITI,1810) (PRWH(I,J),J=1,6)
        WRITE(ISO,1820) I,(PRWH(I,J),J=1,6)
   80 CONTINUE

 1800 FORMAT(//,1X,'PROPERTIES OF WHEEL-AXLE SET',/
     +1X,'WheelNo.',2x,'Mass',10x,'RollingMoment',2x,'YawingMoment', 
     +3x,'NoddingMoment',2x,'WheelRadious',3x,'WheelDistance')
 1810 FORMAT(6F10.0)
 1820 FORMAT(I3,4X,6E15.6)

C	---------------------------------------------------
C	9a.WHEEL-RAIL CONTACT MECHANISM (Three-dimmensional): 
C	---------------------------------------------------
	READ(ITI,*)
      WRITE(ISO,1900)

      DO 90 I=1,NWHEEL
        READ(ITI,1910) (WRCO(I,J),J=1,6)
        WRITE(ISO,1920) I,(WRCO(I,J),J=1,6)
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

C	---------------------------------
C	10a.WHEEL POSITION:
C	   Modified 30Oct2007 (add BIDIS)
C	---------------------------------   
      READ(ITI,*)
	WRITE(ISO,1950)
	 
c	Read the relative wheel positions
      READ(ITI,1960) (WPOS(I),I=1,NWHEEL)

	DO 95 I=1,NWHEEL	 
C		Write the relative wheel positions		 
		WRITE(ISO,1970) (I,WPOS(I))
c		Calculate absolute wheel position	!added 30Oct07
		WPOS(I) = WPOS(I)-BIDIS
   95	CONTINUE

c 1950 FORMAT(//,1X,'WHEEL POSITION',/ 
c     +1X,'WheelNo.	   Wheel Position')		!changed next 30Oct07
 1950 FORMAT(//,1X,'RELATIVE WHEEL POSITION',/
     +1X,'WheelNo.	   Position relative to 1st wheel')
 1960 FORMAT(2F10.0)
 1970 FORMAT(I4,6X,E15.4)

C	---------------------------------------------------------
C	B.CONDITION IF GID FORMAT DATA (NEW) IS USED:
C	To simplify the GiD input by input only the properties of 
C	typical cars/bodies
C	---------------------------------------------------------	
	ELSEIF (IDATM.EQ.1) THEN
C
C	Read Title for Train Properties
	READ(ITI,*)

C	---------------------
C	1b.CAR BODY PROPERTIES 
C	---------------------
	READ(ITI,*)
	WRITE(ISO,1100)
      WRITE(ITO,*) ('Read & write car body properties') 
      WRITE(10,*) ('Read & write car body properties') 

C	Read Number of Typical Car Bodies
	READ(ITI,*)NTCARB

	CALL GENTRAIN(NCARB,4,NTCARB,PCAB) 

      DO 110 I=1,NCARB
        WRITE(ISO,1120) I,(PCAB(I,J),J=1,4)
  110 CONTINUE


C	-----------------------------------
C	2b.DIMMENSION PARAMETERS OF CAR BODY:
C	-----------------------------------
	READ(ITI,*)
      WRITE(ISO,1210)
      WRITE(ITO,*) ('Read & write dimension parameters of car bodies')

C	Read Number of Typical Sets of Dimmension Parameters 
	READ(ITI,*)NTDIM

	CALL GENTRAIN(NCARB,7,NTDIM,CADI)

      DO 120 I=1,NCARB
        WRITE(ISO,1220) I,(CADI(I,J),J=1,7)
  120 CONTINUE
C	-----------------------------------------------------------
C	3b.ADDITIONAL DIMMENSION PARAMETERS OF ARTICULATED CAR BODY:
C	-----------------------------------------------------------
	IF(IT2.GE.3) THEN	!added 10Feb2006

	READ(ITI,*)
      WRITE(ISO,1310)
      WRITE(ITO,*) ('Read & write additional dimension parameters for 
	1Articulated Cars')
      WRITE(10,*) ('Read & write additional dimension parameters for 
	1Articulated Cars')

C	Read Number of Typical Sets of Additional Dimmension Parameters 
	READ(ITI,*)NTADIM

	CALL GENTRAIN(NACAB,5,NTADIM,ADDI)

      DO 130 I=1,NACAB
        WRITE(ISO,1320) I,(ADDI(I,J),J=1,5)
  130 CONTINUE
C	----------------------------------------------
C	4b.LONGITUDINAL SUSPENSIONS OF ARTICULATED CARS:
C	----------------------------------------------
	READ(ITI,*)
      WRITE(ISO,1410)
      WRITE(ITO,*) ('Read & write longitudinal suspensions of Articulate
	1d Cars')

C	Read Number of Typical Sets of Longitudinal Suspensions
	READ(ITI,*)NTLSUS

c	CALL GENTRAIN(NACAB-1,3,NTADIM,ALSU)	!Found error 8Feb2006
	CALL GENTRAIN(NACAB-1,3,NTLSUS,ALSU)

      DO 140 I=1,NACAB-1
        WRITE(ISO,1420) I,(ALSU(I,J),J=1,3)
	  DO 140 J=4,6
		ALSU(I,J) = 0.0  ! Assign the nominal zero values
  140 CONTINUE

	ENDIF !	IF(IT2.GE.3)
C	-----------------------
C	5b.SECONDARY SUSPENSIONS:
C	-----------------------
	READ(ITI,*)
      WRITE(ISO,1500)
      WRITE(ITO,*) ('Read & write secondary suspensions')
      WRITE(10,*) ('Read & write secondary suspensions')

C	Read Number of Typical Secondary Suspensions
	READ(ITI,*)NTSSUS

	CALL GENTRAIN(NBOGI,9,NTSSUS,PRSS)
      DO 150 I=1,NBOGI
        WRITE(ISO,1520) I,(PRSS(I,J),J=1,9)
  150 CONTINUE
C	-----------------------
C	6b.PROPERTIES OF BOGIES:
C	-----------------------  
      READ(ITI,*)
	WRITE(ISO,1600)
      WRITE(ITO,*) ('Read & write properties of bogie')
      WRITE(10,*) ('Read & write properties of bogie')

C	Read Number of Typical Properties of Bogie
	READ(ITI,*)NTPBOG

	CALL GENTRAIN(NBOGI,4,NTPBOG,PBOG)

      DO 160 I=1,NBOGI
        WRITE(ISO,1620) I,(PBOG(I,J),J=1,4)
  160 CONTINUE
C	----------------------
C	7b.PRIMARY SUSPENSIONS:
C	----------------------
	READ(ITI,*)
      WRITE(ISO,1700)
      WRITE(ITO,*) ('Read & write primary suspensions')
      WRITE(10,*) ('Read & write primary suspensions')

C	Read Number of Typical Primary Suspensions
	READ(ITI,*)NTPSUS

	CALL GENTRAIN(NWHEEL,9,NTPSUS,PRPS)

      DO 170 I=1,NWHEEL
        WRITE(ISO,1720)I,(PRPS(I,J),J=1,9)
  170 CONTINUE
C	---------------------------
C	8b.WHEEL-AXLE SET PROPERTIES: 
C	---------------------------
      READ(ITI,*)
	WRITE(ISO,1800)
      WRITE(ITO,*) ('Read & write wheel-axle set properties')
      WRITE(10,*) ('Read & write wheel-axle set properties')

C	Read Number of Typical Wheel-Axle Set Properties
	READ(ITI,*)NTWHPR

	CALL GENTRAIN(NWHEEL,6,NTWHPR,PRWH)

      DO 180 I=1,NWHEEL
        WRITE(ISO,1820) I,(PRWH(I,J),J=1,6)
  180 CONTINUE
C	---------------------------------------------------
C	9b.WHEEL-RAIL CONTACT MECHANISM (Three-dimmensional): WRCO(kwX,cwX,kwH,cwH,kwV,cwV)
C	---------------------------------------------------
	READ(ITI,*)
c      WRITE(ISO,1900)  !Moved down 13Apr08
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
	  WRITE(ISO,1900)
	  READ(ITI,*)NTWRCO		!added 13Apr08
	CALL GENTRAIN(NWHEEL,6,NTWRCO,WRCO)

      DO 190 I=1,NWHEEL
	    WRCO(I,4) = WRCO(I,4)/VEL   !cwH = {1.022284*E**(1/3)*P**(2/3)}/V(m/sec)
        WRITE(ISO,1920) I,(WRCO(I,J),J=1,6)
  190 CONTINUE

	ELSEIF(IWRIN.EQ.1) THEN !Added Dec2008,new
c	    -Input the wheel-rail properties (E,v,Rr, static/dynamic wheel force P)
c	    -Calculate contact stiffness & damping kwH,cwH,kwV
		READ(ITI,*)Ewr,WMU,RDR,CANG
C			Ewr = Elastic modulus of rail and wheel material	
C			WMU = Poisson ratio 
C			RDR = Radius of rail head
C			CANG= Contact angle between the wheel and rail
	    WRITE(ISO,2010) Ewr,WMU,RDR,CANG  !Added Dec2008,new
	ENDIF

c	Next added Dec2008,new
 2010 FORMAT(//,1X,'MECHANICAL PROPERTIES OF WHEEL-RAIL CONTACTS',/
     +1X,'Elastic modulus of rail and wheel material (Ewr) = ',E11.4/
     +1X,'Poisson ratio . . . . . . . . . . . . . . .(WMU) = ',F11.4/
     +1X,'Radius of rail head . . . . . . . . . . . .(RDR) = ',F11.4/
     +1X,'Contact angle between the wheel and rail. (CANG) = ',F11.4)
C	------------------
C	10b.WHEEL POSITION:
C	   Modified 30Oct2007 (add BIDIS)
C	------------------   
      READ(ITI,*)
	WRITE(ISO,1950)
      WRITE(ITO,*) ('Read & write wheel positions')
      WRITE(10,*) ('Read & write wheel positions')
	 
	DO 195 I=1,NWHEEL
c		Read the relative wheel positions
		READ(ITI,*) (NWHNUM(I),WPOS(I))	 
C		Write the relative wheel positions		 
		WRITE(ISO,1970) (I,WPOS(I))

c		Calculate absolute wheel position	!added 30Oct07
		WPOS(I) = WPOS(I)-BIDIS
  195	CONTINUE

C	<---------------End IF (IDATM.EQ.1) 
	ENDIF
C  
  600 RETURN

 9000 FORMAT(1X,'NUMBER OF INTERMEDIATE AND PASSENGER CARS SHOULD BE',/,
     &1X,'EQUAL TO 0 OR GREATER THAN AND EQUAL TO 3',/,
     &1X,'IT2 = ',I5)
      END
C
C	=======================================================================
      SUBROUTINE GENWHPRO(PRWH,FAXL,WRCO,VEL)             
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO GENERATE WHEEL-RAIL CONTACT STIFFNESS & DAMPING FROM MECHANICAL PROP
C	December 2008, new by Nguyen
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
	COMMON /NIAX/ NWHEEL
	COMMON /WRCT/ Ewr,WMU,RDR,RDW,CANG,SH,D11,D22,D23,D33,DHV,DV
C			
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
C	  C1=((Ewr/(1-WMU**2))**(2/3))*((PRWH(IW,5)*RDR)**(1/6))
C	  WRCO(IW,5) = ((1.5*FAXL(IW))**(1/3))*C1
	  WRCO(IW,5) = CKV*(FAXL(IW)**(1./3.))*(PRWH(IW,5)**(1./6.))
	  WRCO(IW,6) = 0.0
	  WRCO(IW,3) = 0.8134*WRCO(IW,5)
C	  WRCO(IW,4) = 1.02284*(Ewr**(1/3))*(FAXL(IW)**(2/3))/VEL
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
c      SUBROUTINE REPRCONT(IPRTTR,IPRTST,NPRTTR,NPRTST,IWHPR,IWDOFPR) !changed to next 14Nov06
C      SUBROUTINE REPRCONT(IPRTTR,IPRTST,NPRTTR,NPRTST)  !changed next 03Mar07
      SUBROUTINE REPRCO(IPRTTR,IPRTST,IPRWH,NPRTTR,NPRTST,NPRWH)

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO READ:
C	  - DEGREES OF FREEDOM OF TRAIN TO PRINT
C	  - WHEEL NUMBER & DOFs (TO PRINT RELATVE DISPLACEMENT)
C	  - NODES NUMBER OF STRUCTURE IPRTST(NPRTST,2)[NODE #, DOF #]
C	----------------------------------------------------------------------
C	INPUT:
C	NPRTTR : NUMBER OF DEGREES OF FREEDOM TO PRINT OF TRAIN
C	NPRTST : NUMBER OF NODES TO PRINT OF STRUCTURE
C	NPRWH  : Number of wheels to print relative displacements

C	OUTPUT:
C	------
C	IPRTTR(NPRTTR)   : DEGREES OF FREEDOM OF TRAIN TO PRINT
C	IPRTST(NPRTST,2) : NODES NUMBER OF STRUCTURE[NODE #, DOF #]
C	IPRWH(NPRWH,2)   : Wheel number & DoFs to print relative displacements
C				       DoF (Direction) of wheel: 
C						= 1: Lateral, = 2: Rolling, = 3: Vertical
C	IWHPR			 : WHEEL NUMBER TO PRINT AXLE LOADS
C	IWDOFPR			 : DoF (Direction) of wheel to print responses !added 6Nov06
C						= 1: Lateral
C						= 2: Rolling
C						= 3: Vertical
C	---------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION IPRTTR(NPRTTR),IPRTST(NPRTST,2),IPRWH(NPRWH,2) !IPRWH added 03Mar07 
C	---------------------------------------------------------------------
C	READ DEGREES OF FREEDOM OF TRAIN TO PRINT:
C	------------------------------------------
      READ(ITI,*)
      IF(NPRTTR.GT.0) THEN
	  DO 100 II = 1,NPRTTR  
  100		READ(ITI,*) IDT,IPRTTR(II)
  
c	  Next added 6Nov06 to read Wheel No. & Direction (DoF number) to print responses
C	  READ(ITI,*)
c       READ(ITI,*) IWHPR,IWDOFPR  !changed to Do-loop 150, 03Mar07 (read number of DOFs)

        WRITE(ISO,1000)
        WRITE(ISO,1100) (IPRTTR(I),I=1,NPRTTR)
	  WRITE(ITO,*)('READ TRAIN DEGREES OF FREEDOM TO PRINT RESPONSE')
        WRITE(10,*)('READ TRAIN DEGREES OF FREEDOM TO PRINT RESPONSE')
      ENDIF

C	READ WHEEL NUMBER & DOFs (TO PRINT RELATVE DISPLACEMENT), added 03Mar07 
C	---------------------------------------------------------
      READ(ITI,*)
      IF(NPRWH.GT.0) THEN
        DO 150 I=1,NPRWH
		READ(ITI,*) (IPRWH(I,J),J=1,2)
  150   CONTINUE

        WRITE(ISO,1300)
        DO 180 I=1,NPRTST
          WRITE(ISO,1100) (IPRWH(I,J),J=1,2)
  180   CONTINUE
	  WRITE(ITO,*)('READ WHEEL NUMBER AND ITS DOFs TO PRINT')
        WRITE(10,*)('READ WHEEL NUMBER AND ITS DOFs TO PRINT')
      ENDIF

C	READ NODES NUMBER OF STRUCTURE IPRTST(NPRTST,2)[NODE #, DOF #]:
C	------------------------------------------------------------
      READ(ITI,*)
      IF(NPRTST.GT.0) THEN
        DO 200 I=1,NPRTST
		READ(ITI,*) (IPRTST(I,J),J=1,2)
  200   CONTINUE

        WRITE(ISO,1500)
        DO 250 I=1,NPRTST
          WRITE(ISO,1100) (IPRTST(I,J),J=1,2)
  250   CONTINUE
	  WRITE(ITO,*)('READ NODE NUMBER AND ITS DOFs OF BRIDGE TO PRINT')
        WRITE(10,*)('READ NODE NUMBER AND ITS DOFs OF BRIDGE TO PRINT')
      ENDIF
C     
      RETURN

 1000 FORMAT(//,1X,'PRINTED DEGREES OF FREEDOM OF TRAIN')
 1100 FORMAT(10I5)
 1300 FORMAT(//,1X,'PRINTED WHEEL NUMBER AND ITS DOF',/
     *'  WHEEL NO.  DOF')  !added 03Mar07
 1500 FORMAT(//,1X,'PRINTED NODE NUMBER AND ITS DOF BRIDGE',/
     *'  NUM  DOF')
      END
C
C	=======================================================================
      SUBROUTINE PCAMAT(PCAB,CADI,PRSS,PBOG,PRPS,PRWH,WRCO,AMTR,ACTR,
	1				  AKTR,IPM,IPK,I)     
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
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL
C	Next added 21Mar06
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH
C
      DIMENSION PCAB(NCARB,*),CADI(NCARB,*),PRSS(NBOGI,*),PBOG(NBOGI,*)
	DIMENSION PRPS(NWHEEL,*),PRWH(NWHEEL,*),WRCO(NWHEEL,*)	
	DIMENSION AMTR(*),ACTR(*),AKTR(*)
C	--------------------------------------------------------------------------
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
C      SUBROUTINE ARMAT2(ITRAIN,PCAB,CADI,ADDI,ALSU,PRSS,PBOG,PRPS,
C	1				  PRWH,WRCO,ACTR,AKTR,IPK) 
C	---->changed next 19Jan08 (removed [PCAB], [PBOG])
      SUBROUTINE ARMAT2(ITRAIN,CADI,ADDI,ALSU,PRSS,PRPS,PRWH,WRCO,ACTR,
	1				  AKTR,IPK)            
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO FORM DAMPING & STIFFNESS MATRICES OF THE ARTICULATED GROUP
C	----------------------------------------------------------------------
C	INPUT:
C	------
C	ITRAIN(3)	   = Train composition stored in an array
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
C      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
C     1              IFPR(10),IFPL(10)
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL
C
C     DIMENSION ITRAIN(3),PCAB(NCARB,*),CADI(NCARB,*),ADDI(NACAB,*)
C	DIMENSION ALSU(NACAB-1,*),PRSS(NBOGI,*),PBOG(NBOGI,*)
C	DIMENSION PRPS(NWHEEL,*),PRWH(NWHEEL,*),WRCO(NWHEEL,*)
c	----->Previous changed next 19Jan08 (removed [PCAB], [PBOG])
      DIMENSION ITRAIN(3),CADI(NCARB,*),ADDI(NACAB,*),ALSU(NACAB-1,*)
	DIMENSION PRSS(NBOGI,*),PRPS(NWHEEL,*),PRWH(NWHEEL,*)		
	DIMENSION WRCO(NWHEEL,*)			
	DIMENSION ACTR(*),AKTR(*)  !activated 19Jan08		
c	DIMENSION ACTR(*),AKTR(NWKTR)
C	---------------------------------------------------------------------------
C	FORM THE STIFF & DAMPING MATRICES:
C	-------------------------------------------------------
C	1.The first transition car:
C	 Matrices of the first pair (car body & bogie) (K11,C11):
C	-------------------------------------------------------
	  IA = ITRAIN(1)+1		
	  JB = ITRAIN(1)+IA	
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
	  AKTR(IPK+83)=-2.*PRSS(JB,3)*CADI(IA,1)

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
	  ACTR(IPK+83)=-2.*PRSS(JB,4)*CADI(IA,1)
	  
	  IPK=IPK+83

C	---------------------------------------------------------------------------
C	2.The intermidate passenger cars (1 <n <N) and the last transition car (n=N):	
C	 Form the tri-diagonal matrices by looping over each pair (bogie & car body)
C		- The last transition car (n = N), several terms of matrices must be  
C		  formulated seperately since the chacteristics of the transition car
C	---------------------------------------------------------------------------
	DO 30 IA = ITRAIN(1)+2,ITRAIN(1)+ITRAIN(2)
	  JB = ITRAIN(1)+IA		
	  JW1=2*JB-1            
	  JW2=2*JB              
	  IAC= IA-ITRAIN(1)		
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
	  

	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars
		AKTR(IPK+79)= ALSU(IAC-1,1) + 2.*PRSS(JB,1) + ALSU(IAC,1)
	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car
		AKTR(IPK+79)= ALSU(IAC-1,1) + 2.*PRSS(JB,1) + 2.*PRSS(JB+1,1)
	  AKTR(IPK+79)=	2.*PRSS(JB,1) + ALSU(IAC,1)
	  ENDIF

	  AKTR(IPK+80)= 0.
	  AKTR(IPK+81)= 0.
	  AKTR(IPK+82)= 0.

	  AKTR(IPK+83)=-2.*PRSS(JB,1)*PRSS(JB,9)
	  AKTR(IPK+84)=-2.*PRSS(JB,1)
	  
  
	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 				
		AKTR(IPK+85)=   ALSU(IAC-1,1)*CADI(IA,6)**2 +
	1				 2.*PRSS(JB,1)*CADI(IA,6)**2 + 
	2				 2.*PRSS(JB,3)*PRSS(JB,8)**2 + 
	3					ALSU(IAC,1)*CADI(IA,7)**2

	    AKTR(IPK+86)=-ALSU(IAC-1,1)*CADI(IA,6)
     1				 - 2.*PRSS(JB,1)*CADI(IA,6)
	2			     -ALSU(IAC,1)*CADI(IA,7) 

	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car
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

	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 	
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

	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car   
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

	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 
	    AKTR(IPK+95)= ALSU(IAC-1,2) + 2.*PRSS(JB,3) + ALSU(IAC,2)
	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car 
	    AKTR(IPK+95)= ALSU(IAC-1,2) + 2.*PRSS(JB,3) + 2.*PRSS(JB+1,3)	 
	  AKTR(IPK+95)=	2.*PRSS(JB,3) + ALSU(IAC,2) 
	  ENDIF	  		
			
	  AKTR(IPK+96)= 0.
	  AKTR(IPK+97)= 0.
	  AKTR(IPK+98)= 0.
	  AKTR(IPK+99)= 0.

	  AKTR(IPK+100)=-2.*PRSS(JB,3)

	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 
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

        ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car 
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


	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars
		ACTR(IPK+79)= ALSU(IAC-1,4)	+ 2.*PRSS(JB,2) + ALSU(IAC,4)	
	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car
		ACTR(IPK+79)= ALSU(IAC-1,4) + 2.*PRSS(JB,2) + 2.*PRSS(JB+1,2)
	  ACTR(IPK+79)=	2.*PRSS(JB,2) + ALSU(IAC,4)	
	  ENDIF

	  ACTR(IPK+80)= 0.
	  ACTR(IPK+81)= 0.
	  ACTR(IPK+82)= 0.

	  ACTR(IPK+83)=-2.*PRSS(JB,2)*PRSS(JB,9)
	  ACTR(IPK+84)=-2.*PRSS(JB,2)

	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 				
	    ACTR(IPK+85)=    ALSU(IAC-1,4)*CADI(IA,6)**2 +
	1			      2.*PRSS(JB,2)*CADI(IA,6)**2 + 
	2				  2.*PRSS(JB,4)*PRSS(JB,8)**2 + 
	3				     ALSU(IAC,4)*CADI(IA,7)**2

	    ACTR(IPK+86)=    -ALSU(IAC-1,4)*CADI(IA,6)
	1				 - 2.*PRSS(JB,2)*CADI(IA,6)
	2			         -ALSU(IAC,4)*CADI(IA,7)  

	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car
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
	  
	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 	
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

	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car   
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

	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 
	    ACTR(IPK+95)= ALSU(IAC-1,5) + 2.*PRSS(JB,4) + ALSU(IAC,5)
	  ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car 
	    ACTR(IPK+95)= ALSU(IAC-1,5) + 2.*PRSS(JB,4) + 2.*PRSS(JB+1,4)	 
	  ACTR(IPK+95)=	2.*PRSS(JB,4) + ALSU(IAC,5)	 
	  ENDIF
	  	
	  ACTR(IPK+96)= 0.
	  ACTR(IPK+97)= 0.
	  ACTR(IPK+98)= 0.
	  ACTR(IPK+99)= 0.

	  ACTR(IPK+100)=-2.*PRSS(JB,4)			

	  IF(IAC.LT.ITRAIN(2)) THEN				!Intermidiate cars 
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

        ELSE IF (IAC.EQ.ITRAIN(2)) THEN		!Last transition car 
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

      SUBROUTINE TRNMAT3D(ITRAIN,PCAB,CADI,ADDI,ALSU,PRSS,PBOG,PRPS,
	1				       PRWH,WRCO,AMTR,ACTR,AKTR)     
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO FORM MASS, DAMPING & STIFFNESS MATRICES OF THE TRAIN SYSTEM
C	----------------------------------------------------------------------
C	INPUT:
C	------
C	ITRAIN(3)	   = Train composition stored in an array
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
C	NACAB = Number of articulated car bodies     = ITRAIN(2)
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
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL

	COMMON /ITRM/ ITRK	!Added here 25Mar06
C
      DIMENSION ITRAIN(3),PCAB(NCARB,*),CADI(NCARB,*),ADDI(NACAB,*)
	DIMENSION ALSU(NACAB-1,*),PRSS(NBOGI,*),PBOG(NBOGI,*)
	DIMENSION PRPS(NWHEEL,*),PRWH(NWHEEL,*),WRCO(NWHEEL,*)	
			
	DIMENSION AMTR(*),ACTR(*),AKTR(*)
c	DIMENSION AMTR(*),ACTR(*),AKTR(NWKTR)
C	---------------------------------------------------------------------------
C	CHECK TRAIN COMPONENTS:
      IF(ITRAIN(2).GT.0.AND.ITRAIN(2).LT.3) THEN
        WRITE(ITO,1000) ITRAIN(2)
        WRITE(10,1000) ITRAIN(2)
        STOP
      ENDIF
	
	CALL CLEARA(AKTR,NWKTR)      
C	------------------------------------------------------------------------------
C 	A.EVALUATION OF SYSTEM MATRICES(MASS, DAMPING, STIFFNESS) OF FRONT POWER CARS:                                                    
C	------------------------------------------------------------------------------
      IPM=0  ! Position pointer for element of mass matrix stored in array
	IPK=0  ! Position pointer for element of stiffness matrix stored in array

      IF(ITRAIN(1).GT.0) THEN
	  DO 10 I=1,ITRAIN(1)
	    CALL PCAMAT(PCAB,CADI,PRSS,PBOG,PRPS,PRWH,WRCO,AMTR,ACTR,AKTR,
	1		       IPM,IPK,I)   
	    IPM = IPM + 27
	    IPK = IPK + 167
   10	  CONTINUE
      ENDIF
C	------------------------------------------------------------
C     B.EVALUATION OF SYSTEM MATRICES(MASS, DAMPING, STIFFNESS) OF               
C	ARTICULATED GROUP (INTERMEDIATE AND PASSENGER CARS):                                    
C	------------------------------------------------------------       
      IF(ITRAIN(2).GE.3) THEN

C	FORM THE MASS MATRIX:
C	---------------------
C	Form the lumped mass matrix of each pair of articulated car (bogie & car body):
C	------------------------------------------------------------------------------	
	DO 20 IA = ITRAIN(1)+1,ITRAIN(1)+ITRAIN(2)	!loop over number of articulated cars
	      JB = ITRAIN(1)+IA		! Corresponding bogie number

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
	  JB = ITRAIN(1)*2+ITRAIN(2)+1

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
	CALL ARMAT2(ITRAIN,CADI,ADDI,ALSU,PRSS,PRPS,PRWH,WRCO,ACTR,AKTR,
	1			  IPK)
C	------------------------------------------------------------------------
C    3.For the last bogie: (n = N+1) and its coupling with the last (Nth) car
C	The matrices must be formulated seperately since its chacteristics
C	------------------------------------------------------------------------
	IA = ITRAIN(1)+ITRAIN(2)
	JB = ITRAIN(1)+IA+1	    
	JW1= 2*JB-1			
	JW2= 2*JB				
	IAC= IA-ITRAIN(1)		

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

      ENDIF   !(ITRAIN(2).GE.3)
C	-----------------------------------------------------------------------------
C	C.EVALUATION OF SYSTEM MATRICES(MASS, DAMPING, STIFFNESS) OF REAR POWER CARS                                                    .
C	-----------------------------------------------------------------------------
	IF(ITRAIN(3).GT.0) THEN 
       DO 50 I = ITRAIN(1)+ITRAIN(2)+1,ITRAIN(1)+ITRAIN(2)+ITRAIN(3)    !Loop over number of rear power cars
	   CALL PCAMAT(PCAB,CADI,PRSS,PBOG,PRPS,PRWH,WRCO,AMTR,ACTR,AKTR,
	1		       IPM,IPK,I)   
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
      SUBROUTINE MAXATRN(ITRAIN,NEQTR,NWKTR,NCHTR,NDTR)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------------------------
C	PROGRAM TO FORM VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS OF TRAIN
C	--------------------------------------------------------------------
C	INPUT:
C	------
C	ITRAIN(3) = TRAIN COMPOSITION
C	NEQTR  = Number of equations of stiffness matrix = 3*(ITRAIN(1)*5) + (2*ITRAIN(2)+1)*5 + 3*(ITRAIN(3)*5)
C	NWKTR  = NUMBER OF ELEMENTS BELOW SKYLINE OF UPPER TRIA. STIFFNESS MATRIX(=ITRAIN(1)*29+(12+ITRAIN(2)*17)+ITRAIN(3)*29)
C
C	OUTPUT:
C	------
C	NCHTR(NEQTR) = VECTOR CONTANING COLUMN HEIGHT OF TRAIN'S STIFF/DAMPING MATRIX
C	NDTR(NEQTR+1)= VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS OF TRAIN STIFFNESS/DAMPING MATRIX
C	---------------------------------------------------------------------	
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /ITRM/ ITRK	!Added here 25Mar06
C
	DIMENSION ITRAIN(*),NCHTR(NEQTR),NDTR(NEQTR+1)
C	--------------------------------------------------------------------------------------
C	MAKING VECTOR CONTANING NUMBER OF ELEMENTS IN COLUMNS OF TRAIN'S STIFF/DAMPING MATRIX
C	--------------------------------------------------------------------------------------
	CALL CLEARI(NCHTR,NEQTR)
	ICOL = 0

C	Front power cars:
C	-----------------
	IF(ITRAIN(1).GT.0) THEN
	DO 60 IP=1,ITRAIN(1)   
		NCHTR(ICOL+1 ) = 1
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 1
		NCHTR(ICOL+6 ) = 6
		NCHTR(ICOL+7 ) = 7
		NCHTR(ICOL+8 ) = 1
		NCHTR(ICOL+9 ) = 6
		NCHTR(ICOL+10) = 1
		NCHTR(ICOL+11) = 6
		NCHTR(ICOL+12) = 6
		NCHTR(ICOL+13) = 5
		NCHTR(ICOL+14) = 9	
		NCHTR(ICOL+15) = 9
		NCHTR(ICOL+16) = 8
		NCHTR(ICOL+17) = 17
		NCHTR(ICOL+18) = 18
		NCHTR(ICOL+19) = 1
		NCHTR(ICOL+20) = 17
		NCHTR(ICOL+21) = 1
		NCHTR(ICOL+22) = 6
		NCHTR(ICOL+23) = 6
		NCHTR(ICOL+24) = 5
		NCHTR(ICOL+25) = 9
		NCHTR(ICOL+26) = 9
		NCHTR(ICOL+27) = 8

		ICOL = ICOL+27
   60 CONTINUE
	ENDIF
C	Articulated group (Transition and Passenger Cars)
C	-----------------
	IF(ITRAIN(2).GE.3) THEN
C   	 1.The first transition car:
		NCHTR(ICOL+1 ) = 1
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 2	!Found error 18Jan
		NCHTR(ICOL+6 ) = 1
		
		NCHTR(ICOL+7 ) = 7
		NCHTR(ICOL+8 ) = 8
		NCHTR(ICOL+9 ) = 9
		NCHTR(ICOL+10) = 8
		NCHTR(ICOL+11) = 9
		
		NCHTR(ICOL+12) = 6				
		NCHTR(ICOL+13) = 7
		NCHTR(ICOL+14) = 8 !3
		NCHTR(ICOL+15) = 6
		NCHTR(ICOL+16) = 7
						 
		ICOL = ICOL+16	

C	 2.The passenger cars (1 <n <N),& last transition car (n = N)
C	   looping over each pair (bogie & car body):
	DO 70 IA = 2,ITRAIN(2)
		NCHTR(ICOL+1 ) = 1
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 2
		NCHTR(ICOL+6 ) = 1
		
		
		NCHTR(ICOL+7 ) = 7 + 0 + 5
		NCHTR(ICOL+8 ) = 8 + 0 + 5
		NCHTR(ICOL+9 ) = 9 + 0 + 5
		NCHTR(ICOL+10) = 8 + 2 + 5
		NCHTR(ICOL+11) = 9 + 2 + 5	

		
		IF(ITRK.EQ.1.OR.ITRK.EQ.2)THEN	!Added here 25Mar06				
			NCHTR(ICOL+12) = 6	 !here			
			NCHTR(ICOL+13) = 7
			NCHTR(ICOL+14) = 3
			NCHTR(ICOL+15) = 6
			NCHTR(ICOL+16) = 7

		ELSE IF(ITRK.EQ.3)THEN
			NCHTR(ICOL+12) = 17  !here				
			NCHTR(ICOL+13) = 18
			NCHTR(ICOL+14) = 19
			NCHTR(ICOL+15) = 17
			NCHTR(ICOL+16) = 18
		ENDIF
			
		ICOL = ICOL+16

   70 CONTINUE

C      3.The last bogie: (n = N+1) 
		NCHTR(ICOL+1 ) = 6
		NCHTR(ICOL+2 ) = 7
		NCHTR(ICOL+3 ) = 3 !1
		NCHTR(ICOL+4 ) = 6
		NCHTR(ICOL+5 ) = 7 !1 !HERE

		NCHTR(ICOL+6 ) = 6
		NCHTR(ICOL+7 ) = 6
		NCHTR(ICOL+8 ) = 5
		NCHTR(ICOL+9 ) = 9
		NCHTR(ICOL+10) = 9
		NCHTR(ICOL+11) = 8

		ICOL = ICOL+11 

	ENDIF
C	Rear power cars:
C	----------------
	IF(ITRAIN(3).GT.0)THEN
	DO 80 IP=1,ITRAIN(3)   
		NCHTR(ICOL+1 ) = 1
		NCHTR(ICOL+2 ) = 2
		NCHTR(ICOL+3 ) = 1
		NCHTR(ICOL+4 ) = 1
		NCHTR(ICOL+5 ) = 1
		NCHTR(ICOL+6 ) = 6
		NCHTR(ICOL+7 ) = 7
		NCHTR(ICOL+8 ) = 1
		NCHTR(ICOL+9 ) = 6
		NCHTR(ICOL+10) = 1
		NCHTR(ICOL+11) = 6
		NCHTR(ICOL+12) = 6
		NCHTR(ICOL+13) = 5
		NCHTR(ICOL+14) = 9	
		NCHTR(ICOL+15) = 9
		NCHTR(ICOL+16) = 8
		NCHTR(ICOL+17) = 17
		NCHTR(ICOL+18) = 18
		NCHTR(ICOL+19) = 1
		NCHTR(ICOL+20) = 17
		NCHTR(ICOL+21) = 1
		NCHTR(ICOL+22) = 6
		NCHTR(ICOL+23) = 6
		NCHTR(ICOL+24) = 5
		NCHTR(ICOL+25) = 9
		NCHTR(ICOL+26) = 9
		NCHTR(ICOL+27) = 8

		ICOL = ICOL+27
   80 CONTINUE
	ENDIF
C	-------------------------------------------------
C	MAKING VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS
C	-------------------------------------------------
	CALL CLEARI(NDTR,NEQTR+1)

      NDTR(1) = 1
      MK      = 0
      DO 100  IEQ=1,NEQTR
		NCHI = NCHTR(IEQ)-1			!COLUMN HEIGHT I (!added 27Nov05)
C		IF (NCHTR(IEQ).GT.MK) MK = NCHTR(IEQ)
C		IF (NCHI.GT.MK) MK = NCHI			  !Removed 14Feb06
C  100 NDTR(IEQ+1) = NDTR(IEQ) + NCHTR(IEQ) + 1 !(Changed 22Oct2005)
  100 NDTR(IEQ+1) = NDTR(IEQ) + NCHI + 1 !added 27Nov05
C  100 NDTR(IEQ+1) = NDTR(IEQ) + NCHTR(IEQ)
      MK  = MK+1

C	Check by number of elements below skyline of upper tria. stiff/damping matrix ??
      NWKTR1 = NDTR(NEQTR+1) - NDTR(1)
	IF(NWKTR1.NE.NWKTR) THEN
	 WRITE (*,1050)
	 STOP
	ENDIF

C	Mean bandwidth:
      BM  = FLOAT(NWKTR)/FLOAT(NEQTR)
C
      WRITE (ISO,1000)
      WRITE (ISO,1100) (IEQ,NCHTR(IEQ),NDTR(IEQ), IEQ=1,NEQTR)
      WRITE (ISO,1200) NDTR(NEQTR+1)
C
 1000 FORMAT(////18X,43(1H*)/18X,1H*,41X,1H*/
     1      18X,43H* COLUMN  HEIGHTS  AND  BANDWIDTH (NCHTR) */
     2      18X,43H* ADDRESSES OF DIAGONAL ELEMENTS (NDTR)   */
     3      18X,1H*,41X,1H*/18X,43(1H*)////
     4      22X,11HCOLUMN NR.I,5X,8HNCHTR(I),5X,7HNDTR(I)/22X,36(1H-)/)

 1050 FORMAT(5X,'*******  WARNING:INCORRECT NUMBER OF TRAIN STIFFNESS 
     1ELEMENTS NWKTR *******')

 1100 FORMAT ((25X,3(I5,7X)))
 1200 FORMAT (49X,I5)

      RETURN

	END
C
C	=====================================================================
      SUBROUTINE CALCS1(PRWH,WRCO,NWHEEL,CS1) 
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO CALCULATE CONSTANTS FOR FORCES APPLED INTO TRAIN/STRUCTURE
c	 11Feb2007
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	NBOGI = Total number of bogies     = IT1*2+(IT2+1)+IT3*2
C	NWHEEL= Total number of wheels     = NBOGI*2
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jz,Jy,rW,b0)
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwZ,cwZ)

C	OUTPUT:
C	-------
C	CS1(NWHEEL,2) = Constants for forces applied into wheel
C					 CS1(I,1) = 2*(kwV*bo**2 + kwH*rw**2)
C					 CS1(I,2) = 2*(cwV*bo**2 + cwH*rw**2)
C	-------------------------------------------------------------------------------------------
      DIMENSION PRWH(NWHEEL,6),WRCO(NWHEEL,6),CS1(NWHEEL,2)
C	-------------------------------------------------------
	CALL CLEARMAT(CS1,NWHEEL,2)

	DO 10 I = 1,NWHEEL	
	CS1(I,1)=2.*(WRCO(I,5)*(PRWH(I,6)**2)+WRCO(I,3)*(PRWH(I,5)**2))    
	CS1(I,2)=2.*(WRCO(I,6)*(PRWH(I,6)**2)+WRCO(I,4)*(PRWH(I,5)**2))
   10 CONTINUE	 
 
      RETURN
	END
C
C	=======================================================================  
      SUBROUTINE CALSAXL(PCAB,CADI,PBOG,PRWH,SCALE,FAXL)	
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	------------------------------------------------------------------------------------------------
C	PROGRAM TO CALCULATE STATIC AXLE FORCE (VERTICAL FORCE FROM WHEEL POSITION DUE TO TRAIN WEIGHT)
C	Programmed 19June2007
C	Later, should modify to use in Moving Axle Loads analysis (the input of PCAB,CADI,PBOG,PRWH)
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
C	NACAB		  = Number of articulated car bodies = IT2
C	NCARB		  = Total number of car bodies = IT1+IT2+IT3
C	NBOGI		  = Total number of bogies     = IT1*2+(IT2+1)+IT3*2
C	NWHEEL		  = Total number of wheels     = NBOGI*2
C	SCALE		  = Gravity acceleration (g)
C	-------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
c	COMMON /DMSO/ IPATH,ALPHA,BITA,DELT,VEL,EXTDIS,NMAX,TOLER,SCALE,
c     +			  NPRTTR,NPRTST,NPRWH,KPRTR,TLENGTH		
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	  !added IWRIN 07Mar08
	COMMON /NIAX/ NWHEEL
C
	DIMENSION PCAB(NCARB,4),CADI(NCARB,7),PBOG(NBOGI,4),PRWH(NWHEEL,6)
      DIMENSION FAXL(NWHEEL)
C	-----------------------------------------------------------------------------------------
C	CALCULATE TRAIN/RAIL INTERACTION LUMPED FORCE APPLIED TO STRUCTURE AT WHEEL POSITION{FSS}
C	-----------------------------------------------------------------------------------------
	CALL CLEARA(FAXL,NWHEEL)

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
	WRITE(ISO,1100)
	 
	DO 500 I=1,NWHEEL	 
		WRITE(ISO,1200) (I,FAXL(I))
  500	CONTINUE

 1100 FORMAT(//,1X,'CALCULATED STATIC AXLE FORCES',/
     +1X,'Wheel No.	   Vertical Axle Force')
 1200 FORMAT(I4,6X,E15.4)

	RETURN
      END
C	=======================================================================
      SUBROUTINE CALFSEC3D(FSS,FS,IDWHPOS,NEID,XI,NRL)
      
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO ADD BRIDGE ECCENTRICITY EFFECT INTO NODAL FORCE VECTOR 
C	Created 03Sep2006 by NguyenDV
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	FSS(NWHEEL*3) = Lumped force applied to structure at wheel position {Lat,FRxi,Ver}
C	FS(NEQ) = Equivalent nodal force vector aplied to structure (NO ECCENTRICITY)
C	NSF    = Maximum number of DOFs per node (=6 or 7, depends on the type of element)
C	NEQ    = NUMBER OF EQUATIONS OF STRUCTURE
C	NNO    = CURRENT NUMBER OF ELEMENT NODES
C	H4    = Verical distance from mass center of the wheel set to shear center of bridge
C	TREC  = Eccentricity of the track: lateral distance from track center to bridge center

C	NRL(NSF,NSN) = EQUATION NUMBER
C	IDWHPOS(NWHEEL) = ID NUMBER OF ELEMENT CORRESPONDING TO WHEEL POSITION
C	NEID(2,NELW) = NODE NUMBER OF SEGMENT WHERE WHEEL RUNS
C	XI(NWHEEL)      = Natural coordinate of wheel position at the current element
C					      XI(I)= DABS(WHPOS(I)-XYZI)/DABS(XYZI-XYZJ)
C	NWHEEL= Total number of wheels     = NBOGI*2
C	NEQTR  = Number of equations of train
C	NELW   = NUMBER OF SEGMENT WHERE WHEEL RUNS

C	OUTPUT:
C	-------
C	FS(NEQ) = Equivalent nodal force vector aplied to structure (WITH ECCENTRICITY)

C	INDEX:
C	------
C	SFBC0(2)= shape functions of C0-beam element in narural coordinate {1/2(1-r) 1/2(1+r)}
C	FSE(14) = Equivalent nodal forces vector of the current element (from lumped forces at wheel position)
C	-----------------------------------------------------------------------
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /NIAX/ NWHEEL
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
C
      DIMENSION FSS(NWHEEL*3),NEID(2,1),XI(1),NRL(NSF,1)
	DIMENSION IDWHPOS(NWHEEL),SFBC0(2),FSE(14),FS(NEQ)  !Found error FSE(2) 21Mar08
C	-----------------------------------------------------------------------
	CALL CLEARA(FSE,14)
C
      DO 40 I=1,NWHEEL
        IF(IDWHPOS(I).GT.0) THEN
	    ISS = (I-1)*3   

C	CALCULATE EQUIVALENT NODAL FORCES VECTOR OF THE CURRENT BEAM ELEMENT
          IENUM=IDWHPOS(I)		
		CALL CLEARA(SFBC0,2)
		SFBC0(1) = 1.-XI(I)
		SFBC0(2) =    XI(I)

C		Torsional moment at wheel position due to eccentricity
c	    FRXI   = FSS(ISS+1)*ECC + FSS(ISS+3)*H4 !changed to next 28Feb07
c	    FRXI   = FSS(ISS+1)*H4 + FSS(ISS+3)*ECC !changed to next 02Mar07
	    FRXI   = FSS(ISS+1)*H4 - FSS(ISS+3)*ECC

C		Torsional moment at nodes due to eccentricity
		FSE(1) = SFBC0(1)*FRXI	!at 1st node
		FSE(2) = SFBC0(2)*FRXI  !at 2nd node (last node of element)

C	ADD EQUIVALENT NODAL FORCE OF ELEMENT DUE TO ECC. INTO GLOBAL NODAL FORCE
          DO 50 J=1,2
            NDNUM=NEID(IENUM,J)
            IP=NRL(4,NDNUM)			!add to torsional DOF (Rx) only
            IF(IP.GT.0) THEN
              FS(IP)=FS(IP)+FSE(J)
            ENDIF
   50     CONTINUE
        ENDIF
   40 CONTINUE

      RETURN
      END
C
C	=======================================================================

      SUBROUTINE CALRNW(IDWHPOS,NEID,PELW,XI,ID,RS,RSW)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO:
C	  - INTERPOLATE NODAL RESPONSES OF STRUCTURE RS[D,V,A](NEQ) 
C		 INTO RESPONSES AT WHEEL POSITIONS RSW[D,V,A](NWHEEL) (in bridge coordinate, currently
C		 different from train' coordinates)
C	  - APPLCABLE FOR VARIOUS ELEMENT TYPE OF STRUCTURES
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	IDWHPOS(NWHEEL) = ID NUMBER OF ELEMENT CORRESPONDING TO WHEEL POSITION
C	NEID(2,NELW)  = NODE NUMBER OF SEGMENT WHERE WHEEL RUNS
C	XI(NWHEEL)		= NATURAL COORDINATE OF WHEEL POSITION
C	NRL(NSF,NSN)	= EQUATION NUMBER OF STRUCTURE
C	RS(NEQ,3)		= RESPONSE OF STRUCTURE[DIS,VEL,ACCE]
C	NWHEEL			= NUMBER OF WHEELS
C	NELW			= NUMBER OF SEGMENTS WHERE WHEEL RUNS
C	NNO			    = NUMBER OF NODES MODELLING A BEAM ELEMENT (=2 OR 3 FOR BEAM)
C	NSF			    = NUMBER OF DEGREES OF FREEDOM
C	NEQ				= NUMBER OF EQUATIONS OF STRUCTURE

C	OUTPUT:
C	-------
C	RSW(NWHEEL*3,2) = RESPONSES OF STRUCTURE AT WHEEL POSITION[DIS,VEL] (in bridge coordinate,
C					  currently different from train' coordinates)
C					  Where DIS = {Vertical,Rolling,Lateral}T displacements of each wheel i
C	INDEX:
C	-------
C	SBC0(3)	 = Shape functions of C0-beam element (2 or 3 node) in narural coordinate 
C	SBC1(6)	 = Shape functions of C1-beam element
C	SBC11(6)	 = First derivatives of shape functions of C1-beam element
C	RSE(NSF*NNO,3)= ELEMENT RESPONSES WHERE WHEEL RUNS[DIS,VEL,ACCE]
C	-------------------------------------------------------------------------
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
C
      DIMENSION IDWHPOS(NWHEEL),NEID(2,1),PELW(4,1),XI(1),ID(NSF,1)
      DIMENSION VR(3),VS(3),VT(3),TRANS(14,14),RL(14,3)
      DIMENSION RSE(14,3),RS(NEQ,3)
	DIMENSION RSW(NWHEEL*3,1)
C	------------------------------------------------------------------------

	ISS=0  

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
c         DO 60 J=1,3 !			[DIS,VEL,ACCE]  
          DO 60 J=1,2 !			[DIS,VEL]

            CI = 1.0-XI(I) ; CJ = XI(I)
            
			RSW(ISS+1,J)=CI*RL(2  ,J)+		    !Ywi !SONGSAK CHANGED NSF TO 7 OCT2007  vertical
     1					 CJ*RL(7+2,J)   
			RSW(ISS+2,J)=CI*RL(4  ,J)+         !Rxwi(Twist) !SONGSAK CHANGED NSF TO 7 OCT2007  rolling
	1                   CJ*RL(7+4,J)          
			RSW(ISS+3,J)=CI*RL(3  ,J)+		    !Zwi !SONGSAK CHANGED NSF TO 7 OCT2007  lateral
     1					 CJ*RL(7+3,J)


            RSW(ISS+1,J) = RSW(ISS+1,J) - RSW(ISS+2,J)*ECC
            
C     Form C0-shape function for beam
c            CALL CLEARA(SBC1,6)
c            SBC1(1)=1.0-XI(I)					    !Use C0-shape function
c            SBC1(2)=1.0-XI(I)                       !N1 (for shear, node 1)
c            SBC1(3)=0.0D0	                        !N2 (for moment,node 1)
c            SBC1(4)=XI(I)					        !Use C0-shape function
c            SBC1(5)=XI(I)                         !N3 (for shear, node 1)
c            SBC1(6)=0.0D0                         !N4 (for moment,node 1)     
c			RSW(ISS+1,J)=SBC1(2)*RSE(2  ,J)+SBC1(3)*RSE(6  ,J)+		    !Ywi !SONGSAK CHANGED NSF TO 7 OCT2007
c     1					 SBC1(5)*RSE(7+2,J)+SBC1(6)*RSE(7+6,J)   
c			RSW(ISS+2,J)=SBC1(1)*RSE(4  ,J)+SBC1(4)*RSE(7+4,J)          !Rxwi(Twist) !SONGSAK CHANGED NSF TO 7 OCT2007
c			RSW(ISS+3,J)=SBC1(2)*RSE(3  ,J)+SBC1(3)*RSE(5  ,J)+		    !Zwi !SONGSAK CHANGED NSF TO 7 OCT2007
c     1					 SBC1(5)*RSE(7+3,J)+SBC1(6)*RSE(7+5,J)
        
   60     CONTINUE

        ELSE
        
c         DO 70 J=1,3       ! [DIS,VEL,ACCE]  !changed to next 24Feb07
          DO 70 J=1,2       ! [DIS,VEL]
            RSW(ISS+1,J)=0.0D0
            RSW(ISS+2,J)=0.0D0
            RSW(ISS+3,J)=0.0D0
   70     CONTINUE
   
        ENDIF

	  ISS = ISS+3 !Increase of wheel DOFs
   10 CONTINUE

  100 RETURN
      END
      
C	=======================================================================
      SUBROUTINE SHFC0B(XI,SFBC0)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO: 
C	 - CALCULATE INTERPOLATION (SHAPE) FUNCTIONS FOR C0-BEAM ELEMENT
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	XI = XI(I)      = Current natural coordinate of wheel position at the current element
C					      XI(I)= DABS(WHPOS(I)-XYZI)/DABS(XYZI-XYZJ)
C	Note: If the whole array IDWHPOS(NWHEEL) and XI(NWHEEL), the transfer values are larger !
C
C	OUTPUT:
C	-------
C	SFBC0(3)= shape functions of C0-beam element in narural coordinate
C			  Two node beam: {1/2(1-r) 1/2(1+r)}
C	-------------------------------------------------------------------------------
C
      DIMENSION SFBC0(3)
C	--------------------------------------------------------------------
C	IF(ITYPE.EQ.5) THEN
C		SELECT CASE(ISTYP)
C	    CASE(1)	! BEAM21Q (LOW ORDER BEAM ELEMENT BY VIVEK GUPTA)
C		  CALL CLEARA(SFBC0,3)
C		  SFBC0(1) = 1.-XI
C		  SFBC0(2) =    XI
C		CASE(2)	! BEAM21F (COMMON (OLD) BEAM ELEMENT)
C		 SELECT CASE(NNO) !Number of nodes for element (=NNBE)
C	      CASE(2)		  !Two-node beam
C		   CALL CLEARA(SFBC0,3)
C		   SFBC0(1) = 1.-XI
C		   SFBC0(2) =    XI
C	      CASE(3)		  !Three-node beam (not available now)
C		   GOTO 100
C		  END SELECT
C		END SELECT
C	ENDIF
C	Previous changed next 21Mar08
		  CALL CLEARA(SFBC0,3)
		  SFBC0(1) = 1.-XI
		  SFBC0(2) =    XI

  100 RETURN
      END

C	=======================================================================
      SUBROUTINE CWHRES(RSW,WHPOS) 

      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO:
C	Calculate the motions of wheelsets by adding to bridge displacement at wheel positions:
C	 - Effect of eccentricities from bridge center to deck:
C	 - Wheel hunting and,
C	 - Track irregularities
C	----------------------------------------------------------------------
C	INPUT:
C	------
C	RSW(NWHEEL*3,3) = Bridge responses at wheel positions:
C	WHPOS(NWHEEL)   = Wheel positions array

C	OUTPUT:
C	-------
C	RSW(NWHEEL*3,3) = Motions of wheelsets [D,V,A]
C	----------------------------------------------------------------------
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
	COMMON /NIAX/ NWHEEL

	DIMENSION WHPOS(*),RSW(NWHEEL*3,*)
C	-----------------------------------------------------------------------

	ISW = 0
	DO 100 I = 1,NWHEEL

C       Add effects of eccentricities from bridge center to deck:
	  RSW(ISW+1,1) = RSW(ISW+1,1) + H4*RSW(ISW+2,1)
	  RSW(ISW+3,1) = RSW(ISW+3,1) + ECC*RSW(ISW+2,1)

C	  Add Wheel hunting (to only lateral displacement)
	  WHHUNT = 0.
	  RSW(ISW+1,1) = RSW(ISW+1,1) + WHHUNT

C	  Add Track irregularities:
	  CALL TIRRE3D(WHPOS(I),TLENGTH,TIRLAT,TIRROT,TIRVER)

	  RSW(ISW+1,1) = RSW(ISW+1,1) + TIRLAT
	  RSW(ISW+2,1) = RSW(ISW+2,1) + TIRROT
	  RSW(ISW+3,1) = RSW(ISW+3,1) + TIRVER 

	  ISW = ISW+3
  100 CONTINUE
     
      RETURN
      END

C	=======================================================================
      SUBROUTINE TIRRE3D(POSI,TLENGTH,TIRLAT,TIRROT,TIRVER) 

      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	PROGRAM TO COMPUTE/ GENERATE 3D TRACK IRREGULARITIES
C
      COMMON /IOUNIT/IN,IOUT,NSTIF

C	LATERAL IRREGULARITIES (Y-IRRE)
	TIRLAT = 0.

C	ROTATIONAL IRREGULARITIES (Rx-IRRE)	
	TIRROT = 0.
      
	PI = 3.14156  
C	VERITCAL IRREGULARITIES (Z-IRRE)
      TIRVER =SIN(0.5*2.*PI*POSI/TLENGTH)*0.0127
	TIRVER =0.0D0

      RETURN
      END

C	=======================================================================
      SUBROUTINE EFFSTIF3(AM,AC,AK,ND,DELT,ALPHA,BETA,NEQ,NWM,NWK)

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO CALCULATE EFFECTIVE STIFFNESS MATRIX
C			(WITH BOTH LUMPED AND CONSISTENT MASS MATRIX)
C	Programmed 28Oct05 and Revised 13Nov05
C	----------------------------------------------------------------------- 
C	INPUT:
C	------
C	AM(NWK) : MASS MATRIX
C	AC(NWK) : DAMPING MATRIX
C	AK(NWK) : STIFFNESS MATRIX
C	ND(NEQ+1) : VECTOR CONTAINING ADRESSES OF DIAGONAL ELEMENTS OF STIFFNESS MATRIX
C	DELT : INCREMENTAL TIME
C	ALPHA : NEWMARK COEFFICIENT
C	BETA
C	NEQ : NUMBER OF EQUATIONS
C	NWM : NUMBER OF ELEMENTS OF MASS MATRIX
C	NWK : NUMBER OF ELEMENTS OF STIFFNESS MATRIX

C	OUTPUT
C	------
C	AK(NWK) : [L]*[D]*[L]**T OF EFFECTIVE STIFFNESS MATRIX

C	AVERAGE ACCELERATION METHOD : UNCONDITIONALLY STABLE(RECOMMENDED)
C	ALPHA=1/4, BETA=1/2

C	LINEAR ACCELERATION METHOD : CONDITIONALLY STABLE
C	ALPHA=1/6, BETA=1/2
C	-----------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION AM(1),AC(1),AK(1),ND(NEQ+1)
C	-----------------------------------------------------------------------
      IF(ALPHA.GT.0.25.AND.BETA.NE.0.5) THEN
        WRITE(*,1000) ALPHA,BETA
        STOP
      ENDIF
C     ------------------------------------
C	CALCULATE EFFECTIVE STIFFNESS MATRIX
C	------------------------------------
      A1=1./(ALPHA*DELT**2)
      A3=BETA/(ALPHA*DELT)

C	ADD MASS MATRIX:
C	---------------
C	When mass matrix is consistent
 	IF(NWM.EQ.NEQ) GOTO 50
      DO 10 ITE=1,NWK
        AK(ITE)=AK(ITE) + A1*AM(ITE)
   10 CONTINUE
	GOTO 70						!Found error when added lumped mass 16Feb06

C	When mass matrix is lumped (diagonal)
   50	DO 60 I=1,NEQ
	 JDIA = ND(I)
	 AK(JDIA)=AK(JDIA)+ A1*AM(I)
   60 CONTINUE


C	ADD DAMPING MATRIX:
C	-------------------
   70 DO 80 ITE=1,NWK
        AK(ITE)=AK(ITE)+A3*AC(ITE)
   80 CONTINUE

C	-------------------------------------------------------
C	DECOMPOSE EFFECTIVE STIFFNESS MATRIX INTO [L]*[D]*[L]**T
C	-------------------------------------------------------
C	ISH              : EIGENVALUE SHIFT INDEX(ISH=0:NOT USED,ISH=1:USED)
      IF(NEQ.EQ.1) THEN
        AK(1) = 1.0D0/AK(1)
      ELSE
        ISH=0
        CALL DECOMP(AK,ND,NEQ,ISH,IOUT)	! Can be Moved outside 
      ENDIF  
      

      RETURN
 1000 FORMAT(//,1X,'INPUT ERROR IN NEWMARK COEFFICIENT',/
     &1X,'RECOMMENDED VALUES ALPHA = 0.15, BETA=0.5 BUT INPUT VALUES',/
     &1X,'ARE AS FOLLOWING',/
     &1X,'  ALPHA = ',F10.3/
     &1X,'  BETA  = ',F10.3)
      END
C
C	=======================================================================
c      SUBROUTINE TITBRTR(IPRTST,IPRTTR,SOLT,IDMA,IWHPR) Changed to next 12Feb07
C      SUBROUTINE TITBRTR(IPRTST,IPRTTR,SOLT,IDMA) !Changed to next 03Mar07
      SUBROUTINE TITBRTR(IPRTST,IPRTTR,IPRWH,SOLT,TRLENG,IDMA) !Added TRLENG 31Oct07
      
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO:
C		- PRINT OUT TITLE FOR BRIDGE/TRAIN IN DYNAMIC INTERACTION ANALYSIS
c	Input:
C	IPRTTR(NPRTTR)  : DEGREES OF FREEDOM OF TRAIN TO PRINT
C	IPRWH(NPRWH,2)  : Wheel number & DoFs to print relative displacements
C				       DoF (Direction) of wheel: 
C						= 1: Lateral, = 2: Rolling, = 3: Vertical
C	IPRTST(NPRTST,2): NODES NUMBER OF STRUCTURE[NODE #, DOF #]
C	-----------------------------------------------------------------------
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	DIMENSION IPRTST(NPRTST,2),IPRTTR(NPRTTR),IPRWH(NPRWH,2) !IPRWH added 03Mar07
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

	WRITE(ITO,200)	!Added 9Feb2006
      WRITE(10,200)	!Added 9Feb2006

C	PRINT TITLE FOR BRIDGE OUTPUT
	WRITE(100,500)
	WRITE(100,600)
c	WRITE(100,1000) VELKMH,TLENGTH,SOLT,TEXTRA !changed next 31Oct07
	WRITE(100,1000) VELKMH,VEL,TALENG,SOLT,BIDIS,BITIME,TLENGTH,
     +				TRLENG,BRTIME,EXTDIS,TEXTRA

	WRITE(100,1500)(IPRTST(I,1),I=1,NPRTST)
	WRITE(100,2000)(IPRTST(I,2),I=1,NPRTST)
C	WRITE(100,2200)		!added 6Sep06, Removed 17July07
	WRITE(100,2500)

C	PRINT TITLE FOR TRAIN OUTPUT
	IF (IDMA.EQ.3) THEN
		WRITE(500,500)
		WRITE(500,700)
c		WRITE(500,1000) VELKMH,TLENGTH,SOLT,TEXTRA !changed next 31Oct07
		WRITE(500,1000) VELKMH,VEL,TALENG,SOLT,BIDIS,BITIME,TLENGTH,
     +					TRLENG,BRTIME,EXTDIS,TEXTRA

		WRITE(500,2000)(IPRTTR(I),I=1,NPRTTR)
c		WRITE(500,2100)IWHPR,IWDOFPR !Changed to next 03Mar07
		WRITE(500,2100)	
		WRITE(500,2150)(IPRWH(I,1),I=1,NPRWH)
		WRITE(500,2000)(IPRWH(I,2),I=1,NPRWH)		
			
c		WRITE(500,2500) !changed to IF loop 04Mar07
	    IF(KPRTR.EQ.0) THEN
		   WRITE(500,2550)
	    ELSE IF(KPRTR.EQ.1) THEN
		WRITE(500,2500)
		ENDIF
	ENDIF
C
	RETURN

  200	FORMAT (/,1X,'>>> SOLVING TRAIN-TRACK-BRIDGE DYNAMIC INTERACTION')

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

 1500 FORMAT (7X,'CONTROL NODE . . ',13(2X,I5))
c 2000 FORMAT (7X,'DOFs SPECIFIED. .',13(2X,I5))		!changed next 08Dec07
 2000 FORMAT (7X,'DOFs SPECIFIED. .',30(2X,I5))

c 2100 FORMAT (5X,'WHEEL NO.=',I5,' TO PRINT AXLE LOADS APPLIED ON BRIDGE
c     + (the last three columns {Lateral, Rolling, Vertical}'/) !changed to next 12Feb07
 2100 FORMAT (/5X,'LAST COLUMNS ARE RELATIVE DISPLACEMENT OF WHEEL
     + AT SPECIFIED DOF (=1:Lateral, 2:Rolling, 3:Vertical)')
c 2150 FORMAT (7X,'WHEEL NUMBER . . ',13(2X,I5))		!changed next 08Dec07
 2150 FORMAT (7X,'WHEEL NUMBER . . ',30(2X,I5))

 2200 FORMAT (5X,'The last three columns are 3D Rail Irregularites 
     +{Lateral, Rolling, Vertical}'/)		!added 6Sep06

C 2550 FORMAT (/5X,'TIME',6X,'POSW1(m)',5X,'ACCE')			!added 04Mar07
C	Changed next 26Apr08 to print max number of iterations JTEM
 2500 FORMAT (/5X,'TIME',6X,'POSW1(m)',2X,'JTEM',4X,'DISP',10X,'ACCE')
 2550 FORMAT (/5X,'TIME',6X,'POSW1(m)',2X,'JTEM',4X,'ACCE')			!added 04Mar07
	END
C	=======================================================================	  
c      SUBROUTINE TRNDOF(RSW,RT,IPRWH,DISW)	!IPRWH added 03Mar07, changed next 11Jul07
      SUBROUTINE TRNDOF(WHPOS,IDWHPOS,RSW,RT,IPRWH,Gx,Ge,Ga,Gr,DISW,TIRAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM 
C		+ FIND THE SELECTED DEGREE OF FREEDOM OF TRAIN
C		+ CALCULATE RELATIVE DISPLACEMENT OF SELECTED WHEEL 
C	 14Nov2006, Modified 03Mar07
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	WHPOS(NWHEEL)  = Wheel position at current step (=WHPOSN)	
C	IDWHPOS(NWHEEL)= ID number of elements corresponding to wheel position
C	RSW(NWHEEL*3,2)= Response of structure part at wheel position [DWI,VWI] in bridge coordinate
C					 Where DWI = {Vertical,Rolling,Lateral} displacements of each wheel i.
C	RT(NEQTR,3)    = Response of Train part at train's DOFs [DTI,VTI,ATI]
C					  Where DTI = {Yi,Rxi,Rzi,Zi,Ryi}T ={Lateral,Rolling,Yawing,Vertical,Pitching} 
C					  displacements of each rigid body i (car body or bogie)
C	IPRWH(NPRWH,2) = Wheel number & DoFs to print relative displacements
C	                 DoF of wheel(= 1: Lateral, = 2: Rolling, = 3: Vertical)
C	-------
C	OUTPUT:
C	-------
C	DISW(NPRWH)    = Relative displacement of specified wheels and DOFs to the rail		
C	-------------------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
c	COMMON /TRIW/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR !changed next & added I_PRWH 03Mar07
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	  !added IWRIN 07Mar08
	COMMON /NIAX/ NWHEEL
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM		!Added 11Jul07
	COMMON /IR3O/ TIRLAT,TIRROT,TIRVER				!Added 11Jul07
C	!TLENGTH added 12Aug06, NPRWH & KPRTR added 03Mar07

	DIMENSION RSW(NWHEEL*3,2),RT(NEQTR,3)		!changed size of RSW, 28 Feb07
	DIMENSION IPRWH(NPRWH,2),DISW(NPRWH),DSTRW(NPRWH) !Added 03Mar07
	DIMENSION Gx(NFFT),Ge(NFFT),Ga(NFFT),Gr(NFFT),WHPOS(1),IDWHPOS(1)	!Added 11Jul07
	
	DIMENSION TIRAIL(3,1)
C	-------------------------------------------------------------------------------------------
C	Initialization
C	DSTRW = 0.
C	DISW  = 0.
	CALL CLEARA(DISW,NPRWH)
	CALL CLEARA(DSTRW,NPRWH)
C	----------------------------------------------------------------------
C	CALCULATE DISPLACEMENT OF STRUCTURE (RAIL) AT SELECTED WHEEL POSITION
C	----------------------------------------------------------------------
C	The displacement of structure at position of the specified wheels & DoFs at current time
	DO 50 NW = 1,NPRWH
c	  IP = IPRWH(NW,1)*3 - IPRWH(NW,2)  !changed next 11Jul07
	  IP = IPRWH(NW,1)*3-IPRWH(NW,2)+1	
	  DSTRW(NW) = RSW(IP,1) 
   50 CONTINUE
C	---------------------------------------------------------
C 	CALCULATE ABSOLUTE DISPLACEMENT OF SELECTED WHEEL AND DOF
C	---------------------------------------------------------
C	FRONT POWER CARS:     
                                               
      IPF=0  

	IF(IT1.GT.0)THEN	!Case of independent-car train
      DO 100 I=1,IT1		

C	  Power Car body:
	  IPF = IPF + 5	 
C	  Two Bogies:
	  DO 100 KB=1,2
		JB = 2*(I-1) + KB	
	    IPF = IPF + 5

C	    Two wheels of each bogie:
	    DO 100 KW=1,2
		  JW = 2*(JB-1) + KW	

C		  IF(JW.EQ.IWHPR) THEN
C			DISW = RT(IPF+IWDOFPR,1)
C			GOTO 500	!Exist loop since current wheel matches with selected wheel
C	      ENDIF
C		  Changed to next 03Mar07		  		  		  
		  DO 80 NW=1,NPRWH
		    IF(JW.EQ.IPRWH(NW,1)) THEN	  
			   IWDOF = IPRWH(NW,2)
			   DISW(NW) = RT(IPF+IWDOF,1)
	        ENDIF
   80		  CONTINUE
	      IPF = IPF + 3
  100 CONTINUE
	ENDIF

C 	ARTICULATED CARS:                                                    
C	------------------
	IF(IT2.GE.3)THEN	!for case of articulated train

	DO 200 I = IT1+1,IT1+IT2
	    JB = IT1 + I		
C	    Two wheels of front bogie:
	    DO 150 KW=1,2
		  JW = 2*(JB-1) + KW	

c		  IF(JW.EQ.IWHPR) THEN
c			DISW = RT(IPF+IWDOFPR,1)
c			GOTO 500	!Exist loop since current wheel matches with selected wheel
c	      ENDIF
C		  Changed to next 03Mar07		  		  		  
		  DO 180 NW=1,NPRWH
		    IF(JW.EQ.IPRWH(NW,1)) THEN	  
			   IWDOF = IPRWH(NW,2)
			   DISW(NW) = RT(IPF+IWDOF,1)
	        ENDIF
  180		  CONTINUE
	      IPF = IPF + 3
  150     CONTINUE

C		Front bogie 
	    IPF = IPF + 5

C		Car body:
	    IPF = IPF + 5

  200 CONTINUE

C	The last bogie of rear transition car:
C	--------------------------------------
	    JB = 2*IT1 + IT2 + 1			
C		Bogie: 
	    IPF = IPF + 5
C	    Two wheels of last bogie:
	    DO 250 KW=1,2
		  JW = 2*(JB-1) + KW	

c		  IF(JW.EQ.IWHPR) THEN
c			DISW = RT(IPF+IWDOFPR,1)
c			GOTO 500	!Exist loop since current wheel matches with selected wheel
c	      ENDIF
C		  Changed to next 03Mar07		  		  		  
		  DO 220 NW=1,NPRWH
		    IF(JW.EQ.IPRWH(NW,1)) THEN	   
			   IWDOF = IPRWH(NW,2)
			   DISW(NW) = RT(IPF+IWDOF,1)
	        ENDIF
  220		  CONTINUE

	      IPF = IPF + 3
  250     CONTINUE
	ENDIF	!IF(IT2.GE.3), for case of articulated train

C 	REAR POWER CARS:                                                    
C	------------------
	IF(IT3.GT.0)THEN
      DO 300 I=IT1+IT2+1,IT1+IT2+IT3	

C	  Power Car body:
	  IPF = IPF + 5	 
C	  Two Bogies:
	  DO 300 KB=1,2
		JB = 2*I-IT2-1+KB 	
	    IPF = IPF + 5
C	    Two wheels of each bogie:
	    DO 300 KW=1,2
		  JW = 2*(JB-1) + KW	

c		  IF(JW.EQ.IWHPR) THEN
c			DISW = RT(IPF+IWDOFPR,1)
c			GOTO 500	!Exist loop since current wheel matches with selected wheel
c	      ENDIF
C		  Changed to next 03Mar07		  		  		  
		  DO 280 NW=1,NPRWH
		    IF(JW.EQ.IPRWH(NW,1)) THEN	   
			   IWDOF = IPRWH(NW,2)
			   DISW(NW) = RT(IPF+IWDOF,1)
	        ENDIF
  280		  CONTINUE

	      IPF = IPF + 3
  300 CONTINUE
	ENDIF
C	-----------------------------------------------------------------	
C	RE-CALCULATE RAIL IRREGULARITIES AT POSITIONS OF SPECIFIED WHEELS
C	This part added 11July07
C	Currently used for only measured or randomly generated irregularity
C	-----------------------------------------------------------------
	IF(IRIN.NE.1.AND.IRIN.NE.3) GOTO 500

C     Generate rail irregularity
	CALL CALIRR(Gx,Ge,Ga,Gr,WHPOS,1,NWHEEL,TIRAIL) 
	
	DO 400 NW = 1,NPRWH
	  IWP = IPRWH(NW,1)
	  IWDOF = IPRWH(NW,2)

	  TIRVER = TIRAIL(1,NW)
	  TIRROT = TIRAIL(2,NW)
	  TIRLAT = TIRAIL(3,NW)
	  	
	  IF(IDWHPOS(IWP).GT.0) THEN			
		 IF(IWDOF.EQ.1) THEN
			DISW(NW) = DISW(NW)-DSTRW(NW)-TIRLAT
		 ELSEIF(IWDOF.EQ.2) THEN
			DISW(NW) = DISW(NW)-DSTRW(NW)-TIRROT
		 ELSEIF(IWDOF.EQ.3) THEN
			DISW(NW) = DISW(NW)-DSTRW(NW)-TIRVER
	     ENDIF

	  ELSE					
		 IF(IWDOF.EQ.1) THEN
			DISW(NW) = DISW(NW)-TIRLAT
		 ELSEIF(IWDOF.EQ.2) THEN
			DISW(NW) = DISW(NW)-TIRROT
		 ELSEIF(IWDOF.EQ.3) THEN
			DISW(NW) = DISW(NW)-TIRVER
	     ENDIF
	  ENDIF
  400 CONTINUE
	GOTO 1000

C	NO IRREGULARITY CASE
C	---------------------
  500	DO NW=1,NPRWH
		DISW(NW) = DISW(NW)-DSTRW(NW)
	ENDDO	
C
 1000 RETURN
      END
C
C	=======================================================================
C	Added coordinate of 1st wheel WCOD1 11Feb08, JTE added 26Apr08
      SUBROUTINE TROUT3D(RT,TT,IPRTTR,NEQTR,ITE,DDT,FSS,DISW,WCOD1,JTE,FWCON)
	      
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO:
C		- PRINT OUT THE TRAIN RESPONSE (DISPLACEMENT, ACCELERATION) & FIND MAX.
C	    - PRINT OUT RELATIVE DISPLACEMENT OF WHEELS TO RAIL & FIND MAX.
C		- FIND MAX. AXLE LOAD OF EVERY WHEEL UP TO CURRENT TIME
C	Totally modified 04Mar07 for various optional printings.
c	-----------------------------------------------------------------------
C
C	INPUT:
C	------
C	RT(NEQTR,3)		= (I)TH RESPONSE OF TRAIN[DIS, VEL, ACCE] 
C	IPRTTR(NPRTTR)  = DOFs OF TRAIN TO PRINT
C	NPRTTR			= NUMBER OF DOFs OF TRAIN TO PRINT
C	NPRWH		   = NUMBER OF DEGREES OF FREEDOM OF WHEEL TO PRINT
C	KPRTR		   = TRAIN PRINT FLAG (in order to control size of train output file)
C						= 0: displacements of train body is not print, only accelerations are printed
C						= 1: print
C	ITE				= I-TH SOLUTION STEP 
C	DDT				= TIME STEP INCREMENT
C	TT(NEQTR)       = DISPLACEMENT OF TRAIN AT ZERO STEP
C	FSS(NWHEEL*3)   = Lumped force applied to structure at wheel position {Fyi,FRxi,Fzi}
C	NWHEEL
C	DISW(NPRWH)    = Relative displacement of specified wheels and DOFs to the rail
C	JTE			   = Max. number of iterations needed at each time step

C	INDEX
C	-----
C	PRTTR(NPRTTR*2) = TEMPORARY STORAGE TO PRINT TRAIN RESPONSES
C	TIMMAX(4){MAXDIS,MINDIS,MAXACCE,MINACCE} = MAXIMUM TRAIN RESPONSE [D,A] AMONG ALL DOFs AT CURRENT TIME STEP
C	WHPRT(3)		= Temperary storage to print axle loads at selected wheel IWHPR
C	--------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
C
      COMMON /MAXDM/AMXTR(100,3),AMXST(100,6)
C	!Added 16April for Max. axial loads at wheels, may change the size of AMXAXL(NWHEEL,3)
	COMMON /MAX1/AMXAXL(100,4),AMXWH(50,2)
	COMMON /NIAX/ NWHEEL	!Added 16April	
C
CN	DIMENSION RT(NEQTR,*),IPRTTR(NPRTTR),PRTTR(NPRTTR*2),TT(NEQTR)		
	DIMENSION RT(NEQTR,3),IPRTTR(NPRTTR),PRTTR(NPRTTR*2),TT(NEQTR)
c	DIMENSION TIMMAX(4),FSS(NWHEEL*3),WHPRT(3),RSW(1,3)
	DIMENSION TIMMAX(4),FSS(NWHEEL*3),WHPRT(3),RSW(NWHEEL*3,2)  !changed size of RSW, 28 Feb07
	DIMENSION DISW(NPRWH)  !added 03Mar07
	DIMENSION FWCON(NWHEEL,3)
C	-----------------------------------------------------------------------
      IF(NPRTTR.LE.0) GOTO 1000
	CALL CLEARA(PRTTR,NPRTTR*2)
      IPOS=1

C	next IF added 04Mar07 for optional print of displacements in tran body
	IF(KPRTR.EQ.0)THEN
	  NTR = NPRTTR 
C	  FIND ONLY ACCELERATION [A] AT CURRENT DOFs
	  DO 50 I=1,NPRTTR
		IPRT=IPRTTR(I)
c		PRTTR(IPOS)=RT(IPRT,3)/SCALE	!Next line removed (/SCALE) 29July 2007, used 1Nov07
		PRTTR(IPOS)=RT(IPRT,3)
		IPOS=IPOS+1
   50	  CONTINUE
	ELSEIF(KPRTR.EQ.1)THEN
	  NTR = NPRTTR*2 
C	  FIND DISPLACEMENT & ACCELERATION [D,A] AT CURRENT DOFs
      DO 100 I=1,NPRTTR
        IPRT=IPRTTR(I)
        PRTTR(IPOS)=RT(IPRT,1)-TT(IPRT)
        IPOS=IPOS+1
c		PRTTR(IPOS)=RT(IPRT,3)/SCALE	!Next line removed (/SCALE) 29July 2007
		PRTTR(IPOS)=RT(IPRT,3)
        IPOS=IPOS+1
  100 CONTINUE
	ENDIF

C	FIND & STORE MAXIMUM TRAIN RESPONSE [D,A] AT CONSIDERED DOFs 
	IF(ITE.EQ.0)THEN
        DO 120 I=1,NTR
          AMXTR(I,1)=PRTTR(I)
          AMXTR(I,2)=PRTTR(I)
  120 CONTINUE
	ELSEIF(ITE.GT.0)THEN
        DO 150 K=1,NTR
          IF(PRTTR(K).GT.AMXTR(K,1)) AMXTR(K,1)=PRTTR(K)
          IF(PRTTR(K).LT.AMXTR(K,2)) AMXTR(K,2)=PRTTR(K)
  150   CONTINUE  
	ENDIF

C	FIND MAXIMUM TRAIN RESPONSE [D,A] AMONG ALL DOFs AT CURRENT TIME STEP
c	(currently not use)
        DO 180 I=1,NEQTR
	    RTI = RT(I,1) - TT(I)
          IF(RTI.GT.TIMMAX(1)) TIMMAX(1)=RTI
          IF(RTI.LT.TIMMAX(2)) TIMMAX(2)=RTI

          IF(RT(I,3).GT.TIMMAX(3)) TIMMAX(3)=RT(I,3)
          IF(RT(I,3).LT.TIMMAX(4)) TIMMAX(4)=RT(I,3)
  180   CONTINUE

C	FIND MAX. RELATIVE DISPLACEMENT OF WHEELS:
C	------------------------------------------
 1000 IF(NPRWH.LE.0) GOTO 2000

	IF(ITE.EQ.0)THEN	
        DO 200 I=1,NPRWH
          AMXWH(I,1)=DISW(I)
          AMXWH(I,2)=DISW(I)
  200   CONTINUE	
	ELSEIF(ITE.GT.0)THEN
        DO 250 K=1,NPRWH
          IF(DISW(K).GT.AMXWH(K,1)) AMXWH(K,1)=DISW(K)  !Find max
          IF(DISW(K).LT.AMXWH(K,2)) AMXWH(K,2)=DISW(K)  !Find min
  250   CONTINUE
      ENDIF

 2000	IF(NPRTTR.LE.0.AND.NPRWH.LE.0) GOTO 3000
  	
	IF(ITE.EQ.0)THEN
		TIME=0.0D0
	ELSEIF(ITE.GT.0) THEN
		TIME=ITE*DDT
	ENDIF
C	XW1 = VEL*TIME !added 03Mar07, changed next 11Feb08
	XW1 = WCOD1 + VEL*TIME 

	IF(NPRTTR.GT.0.AND.NPRWH.GT.0) THEN
		WRITE(500,8410) TIME,XW1,JTE,(PRTTR(I),I=1,NTR),
     1					(DISW(I),I=1,NPRWH)
	ELSEIF(NPRTTR.GT.0.AND.NPRWH.LE.0) THEN    
		WRITE(500,8410) TIME,XW1,JTE,(PRTTR(I),I=1,NTR)
	ELSEIF(NPRTTR.LE.0.AND.NPRWH.GT.0) THEN  
		WRITE(500,8410) TIME,XW1,JTE,(DISW(I),I=1,NPRWH)
	ENDIF    

C	FIND THE AXLE LOADS AT SELECTED WHEEL (added 16April)
 3000	DO	350 I=1,3
	  IP=(IWHPR-1)*3+I 
	  WHPRT(I)=FSS(IP)
  350 CONTINUE
  		
C	FIND Max.AXLE LOADS UPTO CURRENT TIME STEP:
C	------------------------------------------
	IF(ITE.EQ.0)THEN
	  DO 400 IP=1,3
		DO 400 IW=1,NWHEEL	
			AMXAXL(IW,IP)=FWCON(IW,IP)
  400	  CONTINUE
		DO 410 IW=1,NWHEEL	
			AMXAXL(IW,4 )=FWCON(IW,3 )
  410	  CONTINUE
	ELSEIF(ITE.GT.0)THEN
	  DO 450 IP=1,2
		DO 450 IW=1,NWHEEL	
		   IF(AMXAXL(IW,IP).LT.ABS(FWCON(IW,IP))) AMXAXL(IW,IP)=ABS(FWCON(IW,IP))
  450	  CONTINUE
		DO 460 IW=1,NWHEEL	
		   IF(AMXAXL(IW,3).GT.FWCON(IW,3)) AMXAXL(IW,3)=FWCON(IW,3)
  460	  CONTINUE
		DO 470 IW=1,NWHEEL	
		   IF(AMXAXL(IW,4).LT.FWCON(IW,3)) AMXAXL(IW,4)=FWCON(IW,3)
  470	  CONTINUE
	ENDIF
C
      RETURN
c 8410 FORMAT(F10.5,23(1X,E15.7)) !changed to next 03Mar07
c 8410 FORMAT(F10.5,3X,E10.4,30(1X,E15.7))  !changed to next 14July07
c 8410 FORMAT(F10.5,3X,F10.5,50(1X,E15.7))  !changed next 26Apr08 (added I3)
 8410 FORMAT(F10.5,3X,F10.5,1X,I3,50(1X,E15.7))  !may change if size AMXTR(10,2)
      END
C
C	=======================================================================
      SUBROUTINE DIAMAT(A,NDA,NEQ,NWA,DIA)

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------------
C	DIAGONAL ELEMENTS OF A COMPACTED TRIAGULARIZED MATRIX (STIFF,MASS)
C	----------------------------------------------------------------------
C	INPUT:
C	------
C	A(NWA)    = MATRIX STORED IN A COMPACTED TRIAGULARIZED FORM
C	NDA(NEQ+1)= ADDRESS OF DIAGONAL ELEMENTS
C	NEQ
C	NWA
C
C	OUTPUT:
C	-------
C	DIA(NEQ) = Diagonal elements
C	---------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION A(NWA),NDA(NEQ+1),DIA(NEQ)
C	---------------------------------------------------------------------
C	DIAGONAL ELEMENTS OF MATRIX
	DO 10 II=1,NEQ
		IDIA = NDA(II)
		DIA(II) = A(IDIA)
   10	CONTINUE
C     
      RETURN
      END
C	=======================================================================
      SUBROUTINE MOVMAT(A,B,NR,NC)

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------
C	MOVE MATRIA [A] INTO MATRIX [B] 
C	-------------------------------
C	INPUT:
C	A[NR,NC] = Matrix stored in full form
C	NR		 = Number of rows
C	NC		 = Number of columns
C	OUTPUT:
C	DIA(NEQ) = Diagonal elements
C	------------------------------
      DIMENSION A(NR,1),B(NR,1)

      DO 100 J=1,NC
	  DO 100 I=1,NR
		B(I,J)=A(I,J)
  100 CONTINUE

      RETURN
      END

C	=======================================================================
      SUBROUTINE CHOLUP(A,MAXA,NN)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------------------------------------------------
C	PROGRAM FOR MATRIX DECOMPOSION BY CHOLESKI METHOD OF
C		+ UPPER TRIANGULAR MATRIX STORED IN 1D ARRAY
C		+ VARIABLE BANDWIDTH
C		+ OUTPUT IS [A]=[U]T*[U], WHERE [U] IS STORED IN THE SAME ARAY [A]
C	By NguyenDV April2006
C	Input:
C	------
C	A:		= 1D array of upper triangular part of the symmetric matrix [A]
C	MAXA:   = Array contains the addresses of diagonal elements of [A]
C	NN		= Dimensions of symmetric matrix [A]
C
C	Output
C	------
C	A:		= the upper triangular matrix [U] stored in 1D array

C	Arguments:
C	----------
C	L	   = Row number of the first stored element in j-th column (last in 1D array) of [A] or [U]
C	KI,KJ  = Fictitious address of element Ai0, Aj0
C	X	   = Element Aij
C	---------------------------------------------------------------------   
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
      DIMENSION A(NWKTR),MAXA(NN+1)						!TEMP
c      DIMENSION A(*),MAXA(*)
C
	IF(A(1).LE.0.) THEN
            WRITE(*,2000) 1,A(1)
            STOP
	ENDIF
	A(1) = SQRT(A(1))
	DO 10	J=2,NN
	  KJ =  MAXA(J)+J
	  L  = -MAXA(J+1)+KJ+1
	  DO 20 I=L,J
		X = A(KJ-I)
		KI= MAXA(I)+I
		IF(I.EQ.1) GOTO 20
		LBAR = -MAXA(I+1)+KI+1
		LBAR =  MAX0(L,LBAR)

		IF(LBAR.EQ.J) GOTO 20
		DO 30 K=LBAR,I-1
   30		   X = X-A(KI-K)*A(KJ-K)

   20	  A(KJ-I)=X/A(KI-I)
      
	IF(X.LE.0.)THEN
            WRITE(*,2000) J,X
            STOP
	ENDIF
   10 A(KJ-J) = SQRT(X)
C
      RETURN
 2000 FORMAT(//1X,'*** STOP IN "CHOLUP"',
     1       //,1X,'STOP-STIFFNESS MATRIX NOT POSITIVE DEFINITE',
     2       //,1X,'NONPOSITIVE DIAGONAL CHOLESKY TERM',
     3       I4,//1X,'VALUE=',E20.12)
      END
C
C	=======================================================================
      SUBROUTINE TRSOL3D(ND,WHPOS,IPRTTR,IPRTST,IPRWH,FAXL,
	1				   AK,AM,AC,STIF,MASS,DAMP) 

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 STIF,MASS,DAMP
C	------------------------------------------------------------------------
C	INPUT:
C	------
C	NEID(2,NELW)= NODE NUMBERS CORRESPONDING TO ELEMENT PATH WHERE WHEELS RUN
C	NELW           = NUMBER OF ELEMENTS WHERE WHEEL RUNS
C	WHPOS(NWHEEL)  = WHEEL POSITION AT CURRENT TIME STEP
C	XI(NWHEEL)     = NATURAL COORDINATE OF WHEEL POSITION
C	FT(NEQTR)      = FORCE VECTOR APPLIED TO TRAIN
C	FS(NEQ)        = Equivalent nodal force vector applied to structure (at current time)
C	FSS(NWHEEL*3)  = Lumped force applied to structure at wheel position {Fyi,FRxi,Fzi}

C	RSW(NWHEEL*3,2)= RESPONSES OF TRAIN PATH OF STRUCTURE AT WHEEL POSITION[DIS, VEL]
C					  Where DIS = {Lateral, Rolling, Vertical}T displacements of each wheel i		
C
C	RT(NEQTR,3)    = Response of Train part at train's DOFs [DTI,VTI,ATI]
C					  Where DTI = {Yi,Rxi,Rzi,Zi,Ryi}T ={Lateral,Rolling,Yawing,Vertical,Nodding} 
C								   displacements of each rigid body i (car body or bogie)
C	RTN(NEQTR,3)   = (I+1)TH RESPONSE OF TRAIN

C	ORS(NEQ,3)  = PREVIOUS RESPONSE OF STRUCTURE AT NODE DOFs [DIS, VEL, ACCE] 
C	RS(NEQ,3)   = (I)TH RESPONSE OF STRUCTURE AT NODE DOFs [DIS, VEL, ACCE] 
C	RSN(NEQ,3)  = (I+1)TH RESPONSE OF STRUCTURE

C	AKTR(NWKTR) = STIFFNESS MATRIX OF TRAIN
C	ACTR(NWKTR) = DAMPING MATRIX OF TRAIN
C	AMTR(NWKTR) = MASS MATRIX OF TRAIN
C	AMOTR(NWKTR)= ORIGINAL MASS MATRIX OF TRAIN
C	AKOTR(NWKTR)= ORIGINAL STIFFNESS MATRIX 
C	NDTR(NEQTR+1)= VECTOR CONTAINING ADDRESS OF DIAGONAL ELEMENTS OF TRAIN

C	AK(NWK)  = STIFFNESS MATRIX OF STRUCTURE
C	AC(NWK)  = DAMPING MATRIX OF STRUCTURE
C	AM(NWK)  = MASS MATRIX OF STRUCTURE
C	DK(NEQ)  = Diagonal terms of structure stiffness matrix
C	AMO(NWK) = ORIGINAL MASS MATRIX OF STRUCTURE
C	AKO(NWK) = ORIGINAL STIFFNESS MATRIX OF STRUCTURE
C	DKO(NEQ) = Original Diagonal terms of structure stiffness matrix
C	ND(NEQ+1)= VECTOR CONTAINING ADDRESS OF DIAGONAL ELEMENTS OF STRUCTURE
C	IA(LID)  = NRL(NSF,NSN) : EQUATION NUMBER OF STRUCTURE

C	XX(NEQTR) : WORKING VECTOR
C	SS(NEQ) : WORKING VECTOR

C	ITRAIN(3)     : TRAIN COMPOSITION
C	ITRAIN(1)=IT1 : NUMBER OF FRONT POWER CARS
C	ITRAIN(2)=IT2 : NUMBER OF INTERMEDIATE & PASSENGER CARS
C	ITRAIN(3)=IT3 : NUMBER OF REAR POWER CARS

C	NELW : NUMBER OF ELEMENTS OF TRAIN PATH IN GENERAL
C	NACAB = Number of articulated car bodies = IT2
C	NCARB = Total number of car bodies = IT1+IT2+IT3
C	NBOGI = Total number of bogies     = IT1*2+(IT2+1)+IT3*2
C	NWHEEL= Total number of wheels     = NBOGI*2

C	NWKTR : NUMBER OF ELEMENTS OF STIFFNESS MATRIX OF TRAIN
C	NEQTR : NUMBER OF EQUATIONS OF TRAIN
C	NWK : NUMBER OF ELEMENTS OF STIFFNESS MATRIX OF STRUCTURE
C	NEQ : NUMBER OF EQUATIONS OF STRUCTURE
C	NDOF : NUMBER OF DEGREES OF FREEDOM OF NODE
C	NNBE : NUMBER OF NODES MAKING A BEAM ELEMENT(=2)

C	TLENGTH : TOTAL LENGTH OF PATH WHERE WHEEL RUNS
C	SCALE : GRAVITY ACCELERATION
C	TOLER : TOLERANCE FOR THE ITERATION
C	VEL : VELOCITY(M/SEC) OF TRAIN
C	DELT : DELTA T(TIME INCREMENT)
C	NMAX : NUMBER OF MAXIMUM ITERATIONS
C	EXTDIS : EXTRA DISTANCE

C	ALPHA : NEWMARK BETA COEFFICIENT
C	BETA

C	IPRTTR(NPRTTR) : DEGREES OF FREEDOM OF TRAIN TO PRINT
C	IPRWH(NPRWH,2) : Wheel number & DoFs to print relative displacements
C				       DoF (Direction) of wheel: 
C						= 1: Lateral, = 2: Rolling, = 3: Vertical
C	IPRTST(NPRTST,2) : NODES NUMBER OF STRUCTURE[NODE #, DOF #]
C	NPRTTR : NUMBER OF DEGREES OF FREEDOM TO PRINT OF TRAIN
C	KPRTR	     = TRAIN PRINT FLAG (in order to control size of train output file)
C				     = 0: displacements of train body is not print, only accelerations are printed
C					 = 1: print
C	NPRTST : NUMBER OF NODES TO PRINT OF STRUCTURE

C	INDEX:
C	------
C	WHPOSN(NWHEEL) : POSITION AT NEXT TIME STEP
C	IDWHPOS(NWHEEL) : ID NUMBER OF ELEMENT CORRESPONDING TO WHEEL POSITION
C	TT(NEQTR) : DISPLACEMENT OF TRAIN AT ZERO STEP
C	NSOLT : NUMBER OF TOTAL SOLUTION TIMES
C	------------------------------------------------------------------
      COMMON /IOUNIT/IN,IOUT,NSTIF	!Removed 1Dec2005
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      COMMON /MAXDM/AMXTR(100,3),AMXST(100,6)
	COMMON /MAX1/AMXAXL(100,4),AMXWH(50,2)
C
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
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
	COMMON /SPBC/ NSS,NLSS

      COMMON /DMCO/ IDMA,IANA,IDSO,NNBR,NNRA,NEBEBR,NPBEBR,NLRB,
     +			  NEWKTK,NPWKTK,NEBERA,NPBERA,NDOF,NNBE,NNWK

	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS

	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW

	COMMON /BRIW/ I_PELW,I_NEID,I_PRST

	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL

	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH

	COMMON /TRW3/ I_PCAB,I_CADI,I_ADDI,I_ALSU,I_PRSS,I_PBOG,I_PRPS,
	1			  I_PRWH,I_WRCO,I_WPOS,I_AMTR,I_ACTR,I_AKTR,I_AXLD
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM !added 31Jan2007 by Nguyen
	COMMON /IRW3/ IRRx,IRRe,IRRa,IRRr	 !added 31Jan2007 by Nguyen

      COMMON A(9000000),IA(9000000)
	COMMON /MEMW/ W(7000000),IW(7000000)
C	-------------------------------------------------------------------
      DIMENSION WHPOS(NWHEEL),WHPOSN(NWHEEL),
	1IDWHPOS(NWHEEL),NRL(NSF,NSN)
	DIMENSION XI(NWHEEL),FT(NEQTR),RT(NEQTR,3),RTN(NEQTR,3)
      DIMENSION FS(NEQ),FSS(NWHEEL*3),RS(NEQ,3),ORS(NEQ,2),RSN(NEQ,3)
      DIMENSION RSW(NWHEEL*3,2),AKOTR(NWKTR),AMOTR(NEQTR)
      DIMENSION ND(NEQ+1)
	DIMENSION XX(NEQTR),SS(NEQ),IPRTTR(NPRTTR),IPRTST(NPRTST,2)
	DIMENSION PRTTR(NPRTTR*2),PRTST(NPRTST*2),TT(NEQTR),FWR(NWHEEL,3)
	DIMENSION CS1(NWHEEL,2),ORSW(NWHEEL*3,2),IPRWH(NPRWH,2)
C	! CS1 added 11Feb07, ORSW added 24Feb07, IPRWH added 03Mar07
	DIMENSION DISW(NPRWH)	!added 03Mar07
	DIMENSION FAXL(NWHEEL)	!added 20June07
	
	DIMENSION TIRAIL(3,NWHEEL),RTO(NEQTR,3),FWCON(NWHEEL,3)
C	-------------------------------------------------------------------
	DIMENSION AK(1),AM(1),AC(1),DDK(NEQ),DDM(NEQ)

C	CALCULATE NUMBER OF SOLUTION STEPS NSOLT
c	Calculate the train length
	TRLENG = DABS(WHPOS(NWHEEL)-WHPOS(1))		!Added 30Oct07
c	BIDIS can be (+): a distance before entering
c				 (-): the train starts within bridge spans

c      SOLT=(DABS(WHPOS(NWHEEL))+TLENGTH+EXTDIS)/VEL	!Changed next 30Oct2007
      SOLT=(BIDIS + TLENGTH + EXTDIS + TRLENG)/VEL
      NSOLT=IDNINT(SOLT/DELT)
 
C	Total time considered for rail irregularities
c	TIRR=(DABS(WHPOS(NWHEEL))+TLENGTH+DABS(WHPOS(NWHEEL)-WHPOS(1)) !Changed next 30Oct2007
c    +       +EXTDIS)/VEL				
C	TIRR=(TRLENG + BIDIS + TLENGTH + EXTDIS + BIDIS)/VEL !Found error 23Dec07 (BIDIS), changed next
	TIRR=(TRLENG + BIDIS + TLENGTH + EXTDIS + TRLENG)/VEL
c	Generate 3-D rail irregularities by random functions (using FRA PSD),added 31Jan2007
	IF(IRIN.EQ.3.AND.IRAN.EQ.1)THEN
		CALL RAILFRA(WHPOS(NWHEEL),TIRR,VEL,W(IRRx),W(IRRe),W(IRRa),
     +			    W(IRRr))
	ENDIF
 
C	BACK UP TRAIN MATRICES
	CALL MOVE(W(I_AMTR),AMOTR,NEQTR)
	CALL MOVE(W(I_AKTR),AKOTR,NWKTR)


C	PRINT TITLE OF OUTJOB & OUTRAIN
c	CALL TITBRTR(IW(I_PRST),IW(I_PRTR),SOLT,IDMA) !Removed IWHPR 12Feb07, changed next 03Mar07
	CALL TITBRTR(IW(I_PRST),IW(I_PRTR),IW(I_PTWH),SOLT,TRLENG,IDMA) !Added TRLENG 31Oct07

C	Calculate constants for forces applied into wheels & structure, added 11Feb07
	CALL CALCS1(W(I_PRWH),W(I_WRCO),NWHEEL,CS1)

C	Calculate static axle loads due to train weights, added 20June2007
c	CALL CALSAXL(W(I_PCAB),W(I_CADI),W(I_PBOG),W(I_PRWH),SCALE,FAXL) !Moved outside Dec2008

C	***********************************************
C 	*		      ZERO STEP                       *                            
C     ***********************************************
C	Asumming Responses of structure at nodal DOFs: RS(NEQ,3)[DIS,VEL,ACCE] 
	CALL CLEARMAT(RS,NEQ,3)

C	Initialize Responses of structure at wheel position
	CALL CLEARMAT(RSW,NWHEEL*3,2)  

C	Identify element number correspoding to wheel position at zero step
      CALL IDWHEEL(NELW,IW(I_NEID),W(I_PELW),W(I_WPOS),IDWHPOS,XI,NWHEEL)

C     INITIALIZE THE INTERACTION FORCE ACTING ON STRUCTURE
      FSS = 0.0D0
      
C	  Initialize Responses of train at train's DOFs at zero step
      CALL CLEARMAT(RT,NEQTR,3)	
      CALL CLEARMAT(RTO,NEQTR,3)	     
C     ---------------------------------------------------------
C 	ITERATION PROCESS OF ZERO STEP: 
C	In order to cal. Responses of structure at nodal DOFs: RSN (later, check whether RS)                                   
C     ---------------------------------------------------------
      DO 20 JTE=1,5
C        WRITE(*,2000) JTE

	  IF(IRIN.EQ.2) THEN !IRREGULARITY GENERATION
		    CALL GENIRRE3D(W(I_WPOS),TLENGTH,JTE,NWHEEL,TIRAIL)			
	  ELSEIF(IRIN.EQ.0.OR.IRIN.EQ.1.OR.IRIN.EQ.3) THEN
		    CALL CALIRR(W(IRRx),W(IRRe),W(IRRa),W(IRRr),W(I_WPOS),JTE,NWHEEL,TIRAIL) 
	  ENDIF
 
C	  CALCULATE INTERACTION FORCE APPLIED TO TRAIN'DOFs AT ZERO STEP
	  CALL CALFT3D(RTO,W(I_PCAB),W(I_PBOG),W(I_PRWH),W(I_WRCO),	
	1		       W(I_WPOS),IDWHPOS,FWR,FT,JTE,CS1,FSS) 

C	  CALCULATE DISPLACEMENT THE TRAIN SYSTEM AT TRAIN'DOF:
C	  -----------------------------------------------------	
C	  Calculate equivalent inertial forces {XX]=[Mtr]*{ATr} 
        CALL MULTDM(XX,W(I_AMTR),RT(1,3),IW(I_NDTR),NEQTR,NEQTR)
C	  and add to RHS forces {FT} = {FT} - {XX}
	  CALL VECSUB(FT,XX,FT,FACT,VNORM,NEQTR)

C	  Calculate equivalent damping forces {XX]=[Ctr]*{VTr} and add to RHS forces
        CALL MULTDM(XX,W(I_ACTR),RT(1,2),IW(I_NDTR),NEQTR,NWKTR)
C	  and add to RHS forces {FT} = {FT} - {XX}
	  CALL VECSUB(FT,XX,FT,FACT,VNORM,NEQTR)

C	  Decompose stiffness matrix of train at first iteration JTE=1
C	  BEFORE DECOMP AKTR=ORIGINAL STIFFNESS MATRIX
C	  AFTER DECOMP AKTR=DECOMPOSED STIFFNESS MATRIX
        IF(JTE.EQ.1) THEN
          ISH=0				
          CALL DECOMP(W(I_AKTR),IW(I_NDTR),NEQTR,ISH,IOUT)
        ENDIF

C	  Reduce & back-substitute for train displacement	            
C	  CALL COLSOL (IW(I_NDTR),W(I_AKTR),A(LBB),DIKTR,FT,2,INDPD)	            
        CALL REDBAK(W(I_AKTR),FT,IW(I_NDTR),NEQTR)
	  CALL MOVE(FT,RT(1,1),NEQTR)
	  
C	  CALCULATE INTERACTION FORCE APPLIED TO TRAIN'DOFs AT ZERO STEP
	  CALL CALFT3D(RT,W(I_PCAB),W(I_PBOG),W(I_PRWH),W(I_WRCO),	
	1		       W(I_WPOS),IDWHPOS,FWR,FT,JTE,CS1,FSS) 	  
C
C	  CALCULATE INTERACTION FORCE APPLIED TO STRUCTURE'DOFs AT ZERO STEP
	  CALL CALFS3D(FSS,FS,RSW,FWR,IDWHPOS,IW(I_NEID),W(I_PELW),XI,
     1               IA(LID),W(I_PRWH),W(I_WRCO),CS1,FAXL,TIRAIL,FWCON)			

C	  CALCULATE DISPLACEMENT THE STRUCTURE SYSTEM AT DOFs:
C	  ----------------------------------------------------
C	  Calculate equivalent inertial forces {SS]=[Mst]*{Ast} 
C        CALL MULTDM(SS,AM,RS(1,3),ND,NEQ,NWK)
	  CALL MAMULT(ND,AM,RS(1,3),SS,MASS,'STD')

C	  and add to RHS forces	{FS} = {FS} - {SS}
	  CALL VECSUB(FS,SS,FS,FACT,VNORM,NEQ)

C	  Calculate equivalent damping forces {SS]=[Ctr]*{VTr}
C	  CALL MULTDM(SS,AC,RS(1,2),ND,NEQ,NWK)
	  CALL MAMULT(ND,AC,RS(1,2),SS,DAMP,'STD')

	
C	  and add to RHS forces	{FS} = {FS} - {SS}
	  CALL VECSUB(FS,SS,FS,FACT,VNORM,NEQ)

C	  Decompose stiffness matrix of structure at first iteration JTE=1
C	  BEFORE DECOMP AK=ORIGINAL STIFFNESS MATRIX
C	  AFTER DECOMP AK=DECOMPOSED STIFFNESS MATRIX
        IF(JTE.EQ.1) THEN

C		ISH=0
C		CALL DECOMP(AK,ND,NEQ,ISH,IOUT)
C		INDPD = 0
C		CALL COLSOL(ND,AK,DDK,DDK,1,INDPD,STIF,'TEMP')

        ENDIF


C	  Reduce & back-substitute for structure displacement
C	  CALL REDBAK(AK,FS,ND,NEQ)
	  CALL COLSOL(ND,AK,DDK,FS,1,INDPD,'STIF','TEMP')
	  CALL COLSOL(ND,AK,DDK,FS,2,INDPD,'STIF','TEMP')

	  CALL MOVE(FS,RS(1,1),NEQ)

C	  Store responses of bridge at wheel position RSW as of previous iteration ORSW, added 28Feb07
        CALL CLEARMAT(ORSW,NWHEEL*3,2)
        CALL MOVMAT(RSW,ORSW,NWHEEL*3,2)

C	  Interpolate bridge response at wheel position from current responses of bridge at nodal DOFs
	  IF (JTE.GT.1) THEN !at JTE=1, [RS] = [0]
	  CALL CALRNW(IDWHPOS,IW(I_NEID),W(I_PELW),XI,IA(LID),RS,RSW) 
	  ENDIF


c        WRITE(*,*) FWR(1,3),FSS(3),FWR(1,3)-FSS(3)-FAXL(1)
C	  Calculate the motions of wheelsets by adding to bridge displacement at wheel positions:
C	  - Effect of eccentricities from bridge center to deck,
C	  - Wheel hunting and,
cn	  CALL CWHRES(RSW,W(I_WHPOS)) 

C	  CALCULATE THE DIFFERENCE Btw. CURRENT & PREVIOUS DISPLACEMENTS
C	  Previous changed to next 28Feb07 (Use RSW instead)
C	  Calculate the difference btw current & previous displacements of struc at wheel positions       
        DO 100 K=1,2
          DO 100 L=1,NWHEEL
            ORSW(L,K)=RSW(L,K)-ORSW(L,K)
  100   CONTINUE
C	  Calculate SRSS of displacement differences in struc at wheel positions 
        TOL=0.0
        DO 110 K=1,NWHEEL			
          TOL=TOL+ORSW(K,1)**2
  110   CONTINUE

C	  CHECK FOR NEXT ITERATION
        TOL=DSQRT(TOL)
        IF(TOL.GT.TOLER) GOTO 20
        GOTO 130

   20 CONTINUE

  130 CONTINUE
C     ---------------------
C	OUTPUT FOR ZERO STEP
C	---------------------
      DO 145 I=1,NEQTR
        TT(I)=RT(I,1)
  145 CONTINUE       

C	Find the selected DOF of train and calculate relative wheel displacements
C	CALL TRNDOF(RSW,RT,IPRWH,DISW)		!IPRWH added 03Mar07 , changed next 11July07
	CALL TRNDOF(W(I_WPOS),IDWHPOS,RSW,RT,IPRWH,W(IRRx),W(IRRe),
	1			W(IRRa),W(IRRr),DISW,TIRAIL)

C	Print out train response (D,A) at zero step
C     CALL TROUT3D(RT,TT,IPRTTR,NEQTR,ITE,DELT,FSS,IWHPR)
C	Changed to the next 14Nov06 (moved IWHPR to common block, added DISW)
C	Added coordinate of 1st wheel WHPOS(1) 11Feb08, JTE added 26Apr08
      CALL TROUT3D(RT,TT,IPRTTR,NEQTR,ITE,DELT,FSS,DISW,WHPOS(1),JTE,FWCON)

C	Print out structures response (D,A) at zero step
C	Added original coordinate of 1st wheel WHPOS(1) 11Feb08
	CALL BROUT3D(IA(LID),RS,IW(I_PRST),0,DELT,WHPOS(1))

      IF(JTE.GT.NMAX) THEN
        WRITE(ITO,8476) TIME,NMAX
        WRITE(10,8476) TIME,NMAX
        STOP
      ENDIF
C	----------------------------
C	PREPARE FOR FOLLOWING STEPS
C	----------------------------
C	RESTORE ORIGINAL STIFFNESS MATRIX TO AKTR & AK 
C	(Contain stiffness matrices decomposed in previous soubroutine SUBROUTINE DECOMP)
	CALL MOVE(AKOTR,W(I_AKTR),NWKTR)


C	DECOMPOSED MASS MATRIX OF STRUCTURE
C	 -Before DECOMP,  AM   = Original mass matrix of structure
C	 -After  DECOMP,  AM   = Decomposed mass matrix of structure
C	ISH=0
C	CALL DECOMP(AM,ND,NEQ,ISH,IOUT)
C	INDPD = 0
C	CALL COLSOL(ND,AM,DDM,DDM,1,INDPD,MASS,'MAS0')


C	CALCULATE EFFECTIVE STIFFNESS MATRIX OF TRAIN
C	 - Before EFFSTIF, AKTR = Original stiffness matrix
C	 - After  EFFSTIF  AKTR = Decomposed Effective stiffness matrix
      CALL EFFSTIF3(AMOTR,W(I_ACTR),W(I_AKTR),IW(I_NDTR),DELT,ALPHA,
     1	     	  BITA,NEQTR,NWMTR,NWKTR)	

C	CALCULATE EFFECTIVE STIFFNESS MATRIX OF STRUCTURE
C	 - Before EFFSTIF, AK = Original stiffness matrix
C	 - After  EFFSTIF  AK = Decomposed Effective stiffness matrix
C		Check stiffness matrix of BRIDGE by dia-elements
cn	CALL PRDIAMAT(AK,ND,NEQ,NWK)


C	CALL EFFSTIF3(AMO,AC,AK,ND,DELT,ALPHA,BITA,NEQ,NWK,NWK)
      CALL EFFSTIFB(AM,AC,AK,ND,DELT,ALPHA,BITA,NEQ,NWK,DDK,
	1			 'STIF','MASS','DAMP','EFTF','TEMP')
C
C     ******************************************************************
C     *							FOLLOWING STEP                       *                                       
C     ******************************************************************
      DO 140 ITE=1,NSOLT
       WRITE(*,3000) ITE,DFLOAT(ITE)/NSOLT*100.

C	  CALCULATE (ITE+1)TH WHEEL POSITION
	  CALL CLEARA(WHPOSN,NWHEEL)	 
        DO 175 J=1,NWHEEL
          WHPOSN(J)=WHPOS(J)+ITE*DELT*VEL
  175   CONTINUE

C	  IDENTIFY ELEMENT NUMBER CORRESPONDING TO WHEEL POSITION AT (ITE+1)TH STEP
	  CALL CLEARI(IDWHPOS,NWHEEL)
        CALL IDWHEEL(NELW,IW(I_NEID),W(I_PELW),WHPOSN,IDWHPOS,XI,NWHEEL)

C       GATHER IRREGULARITY AT THe WHEEL POSITION
	  IF(IRIN.EQ.2) THEN !IRREGULARITY GENERATION
		    CALL GENIRRE3D(WHPOSN,TLENGTH,1,NWHEEL,TIRAIL)			
	  ELSEIF(IRIN.EQ.0.OR.IRIN.EQ.1.OR.IRIN.EQ.3) THEN
		    CALL CALIRR(W(IRRx),W(IRRe),W(IRRa),W(IRRr),WHPOSN,1,NWHEEL,TIRAIL) 
	  ENDIF
	  
        
      CALL MOVMAT(RT,RTN,NEQTR,3)
C     -----------------------------------
C 	ITERATION PROCESS OF FOLLOWING STEP                                
C     -----------------------------------
        DO 150 JTE=1,NMAX

C	    Store responses of bridge at wheel position RSW as of previous iteration ORSW, added 28Feb07
          CALL MOVMAT(RSW,ORSW,NWHEEL*3,2)
	  
C	    CALCULATE INTERACTION FORCE APPLIED TO TRAIN'DOFs AT (JTE)th STEP {FT}
	    CALL CALFT3D(RTN,W(I_PCAB),W(I_PBOG),W(I_PRWH),W(I_WRCO),
	1		         WHPOSN,IDWHPOS,FWR,FT,JTE,CS1,FSS) 

C		ANALYZE TRAIN SYSTEM: 
C	    ---------------------
C			AMOTR = Original mass matrix of train
C			ACTR  = Original damping matrix of train
C			AKOTR = Original stiffness matrix of train
C			AKTR  = Decomposed Effective stiffness matrix
C			AMTR  = Decomposed mass matrix of train
C		USING NEWMARK BETA METHOD TO CAL. DISP,VELO & ACCE AT (ITE+1)STEP [RTN(NEQTR,3)]
          CALL NMSOL(AMOTR,W(I_ACTR),W(I_AKTR),RT,RTN,FT,IW(I_NDTR),
     1               DELT,ALPHA,BITA,NEQTR,NWMTR,NWKTR)

C	    CALCULATE INTERACTION FORCE APPLIED TO TRAIN'DOFs AT (JTE)th STEP {FT}
	    CALL CALFT3D(RTN,W(I_PCAB),W(I_PBOG),W(I_PRWH),W(I_WRCO),
	1		         WHPOSN,IDWHPOS,FWR,FT,JTE,CS1,FSS) 
	
C		RE-CALCULATE THE ACCELERATION FROM THE EQUILIBRIUM EQUATION
C          CALL CALACC(FT,RTN,AMOTR,W(I_ACTR),AKOTR,XX,IW(I_NDTR),NEQTR,
C	1				NWMTR,NWKTR)

C	    CALCULATE INTERACTION FORCE APPLIED TO STRUCTURE'DOFs AT (ITE+1)STEP {FS}
	    CALL CALFS3D(FSS,FS,RSW,FWR,IDWHPOS,IW(I_NEID),W(I_PELW),XI,
     1                 IA(LID),W(I_PRWH),W(I_WRCO),CS1,FAXL,TIRAIL,FWCON)			
	
c		write(100,*)fs
C		ANALYZE STRUCTURE:
C		------------------ 
C			AMO = Original mass matrix of structure
C			AC  = Original damping matrix of structure
C			AK  = Decomposed Effective stiffness matrix
C			AM  = Decomposed mass matrix of structure
C			AKO = Original stiffness matrix of structure
C		Check stiffness matrix of train by dia-elements
C		CALL PRDIAMAT(AK,ND,NEQ,NWK)
C		WRITE(100,*)RS
C		USING NEWMARK BETA METHOD TO CAL. DISP,VELO & ACCE OF (ITE+1)STEP [RSN(NEQ,3)]

	CALL NMSOLB(AM,AC,AK,RS,RSN,FS,ND,DELT,ALPHA,BITA,NEQ,NWK,NWK,
	1			  DDK,'EFTF','MASS','DAMP')

C		RE-CALCULATE THE ACCELERATION FROM THE EQUILIBRIUM EQUATION [RSN(NEQ,1)]
CC	CALL CALACCB(FS,RSN,AM,AC,AK,SS,ND,NEQ,NWK,NWK,
CC	1			   DDM,'STIF','MASS','DAMP')    

          
C	    Interpolate bridge response at wheel position from responses of bridge at nodal DOFs [RSW]
		CALL CALRNW(IDWHPOS,IW(I_NEID),W(I_PELW),XI,IA(LID),RSN,RSW) 

C	    CALCULATE INTERACTION FORCE APPLIED TO STRUCTURE'DOFs AT (ITE+1)STEP {FS}
	    CALL CALFS3D(FSS,FS,RSW,FWR,IDWHPOS,IW(I_NEID),W(I_PELW),XI,
     1                 IA(LID),W(I_PRWH),W(I_WRCO),CS1,FAXL,TIRAIL,FWCON)	
     
C		Calculate the difference btw current & previous displacements of struc at wheel positions       
          DO 250 K=1,2
            DO 250 L=1,NWHEEL
              ORSW(L,K)=RSW(L,K)-ORSW(L,K)
  250     CONTINUE
C	    Calculate SRSS of displacement differences in struc at wheel positions 
          TOL=0.0
          DO 300 K=1,NWHEEL			
            TOL=TOL+ORSW(K,1)**2
  300     CONTINUE


C        WRITE(*,*) FWR(1,3),FSS(3),FWR(1,3)-FSS(3)
        
C	    CHECK FOR NEXT ITERATION  
          TOL=DSQRT(TOL)
          IF(TOL.GT.TOLER) GOTO 150
          GOTO 141

  150   CONTINUE
C	----------------------------> end of iteration for time step 
C	  STORE RESPONSES OF TRAIN & STRUCTURES FOR NEXT TIME STEP
  141   CALL CLEARMAT(RT,NEQTR,3)
        CALL MOVMAT(RTN,RT,NEQTR,3)
        
c        WRITE(110,*) FWR(1,3),FSS(3),FWR(1,3)-FSS(3)-FAXL(1)
c        WRITE(*,*) FWR(1,3),FSS(3),FWR(1,3)-FSS(3)-FAXL(1)
c      PAUSE
c        WRITE(100,2000) JTE

	  CALL CLEARMAT(RS,NEQ,3)
        CALL MOVMAT(RSN,RS,NEQ,3) 
C     ----------------------
C	OUTPUT FOR (ITE+1)STEP
C	----------------------
C	Find the selected DOF of train and calculate relative wheel displacements
c	CALL TRNDOF(RSW,RT,IPRWH,DISW)		!IPRWH added 03Mar07 , changed next 11July07
	CALL TRNDOF(WHPOSN,IDWHPOS,RSW,RT,IPRWH,W(IRRx),W(IRRe),
	1			W(IRRa),W(IRRr),DISW,TIRAIL)

C	Print out train response (D,A)
C	Added coordinate of 1st wheel WHPOS(1) 11Feb08, JTE added 26Apr08
      CALL TROUT3D(RT,TT,IPRTTR,NEQTR,ITE,DELT,FSS,DISW,WHPOS(1),JTE,FWCON)

C	Print out structure response (D,A)
	CALL BROUT3D(IA(LID),RS,IW(I_PRST),ITE,DELT,WHPOS(1))

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
	CALL PRNFLAG('NONE','NONE','NONE','NONE','DISP','NONE','NONE')
	CALL  PRNOUT('STND','PONE','NONE',ITE)
C	--------------------------------------------------------------

	TIME=ITE*DELT
      IF(JTE.GT.NMAX) THEN
          WRITE(*,8476) TIME,NMAX,TOL
          STOP
      ENDIF

  140 CONTINUE !Continue next time step <<------------------>>
C	-----------------------------------------------
C	OUTPUT OF EXTREMA RESPONSE FOR OVERALL ANALYSIS
C	-----------------------------------------------
C----	PRINT EXTREMA RESPONSE [D,A] AT CONSIDERED TRAIN'S DOFs
	IF(KPRTR.EQ.0) THEN
		NTRPR = NPRTTR
		WRITE(500,2550)
	ELSEIF (KPRTR.EQ.1) THEN 
		NTRPR = NPRTTR*2
c		WRITE(500,2500) !changed next 26Apr08
		WRITE(500,2520)
	ENDIF

c	Changed label 8420 to 8410 26Apr08 in Train'max print		
      DO 700 I=1,2
	 IF(NPRTTR.GT.0.AND.NPRWH.GT.0) THEN
        WRITE(500,8410) (AMXTR(J,I),J=1,NTRPR),(AMXWH(K,I),K=1,NPRWH)
       ELSEIF(NPRTTR.GT.0.AND.NPRWH.LE.0) THEN
           WRITE(500,8410) (AMXTR(J,I),J=1,NTRPR)
	 ELSEIF(NPRTTR.LE.0.AND.NPRWH.GT.0) THEN  
           WRITE(500,8410) (AMXWH(K,I),K=1,NPRWH)
	 ENDIF 
  700 CONTINUE
	IF (NPRTTR.GT.0)
     +  WRITE(500,8410)(MAX(ABS(AMXTR(J,1)),ABS(AMXTR(J,2))),J=1,NTRPR)

C----	PRINT EXTREMA RESPONSE [D,A] AT CONSIDERED STRUCTURE'S DOFs
	WRITE(100,2500)	
      DO 800 I=1,2
        WRITE(100,8420) (AMXST(J,I),J=1,NPRTST*2)
  800 CONTINUE
      WRITE(100,8420)(MAX(ABS(AMXST(J,1)),ABS(AMXST(J,2))),J=1,NPRTST*2)

C-----PRINT MAX.AXLE LOAD APPLIED ON BRIDGE UNDER EVERY WHEELS (added 16Apr06)
	WRITE(500,3200)	
	DO 900 I=1,NWHEEL
		WRITE(500,3500)I,(AMXAXL(I,IP),IP=1,4)
  900	CONTINUE
C  
      RETURN
 1000 FORMAT(1X,'Current Time Step =        0     0.00%')
 2000 FORMAT(1X,'Current Iteration Number = ',I8)
 3000 FORMAT(1X,'Current Time Step = ',I8,3X,F6.2,'%')

C 2500 FORMAT (//15X,'EXTREMA VALUES (MAX & MIN) OVER WHOLE HISTORY AT CU
C     +RRENT DOFs'/
C    +/15X,'DISP.EXTREMA',4X,'ACCE.EXTREMA')
C	Previous changed to next 03Mar07
 2500 FORMAT (//26X,'EXTREMA VALUES (MAX & MIN) OVER WHOLE HISTORY AT CU
     +RRENT DOFs'//26X,'DISP.EXTREMA',4X,'ACCE.EXTREMA')
C	Next label 2520 added 26Apr08	
 2520 FORMAT (//30X,'EXTREMA VALUES (MAX & MIN) OVER WHOLE HISTORY AT CU
     +RRENT DOFs'//30X,'DISP.EXTREMA',4X,'ACCE.EXTREMA')

c 2550 FORMAT (//26X,'EXTREMA VALUES (MAX & MIN) OVER WHOLE HISTORY AT CU
c    +RRENT DOFs'//26X,'ACCE.EXTREMA')
c	Previous changed to next 26Apr08
 2550 FORMAT (//30X,'EXTREMA VALUES (MAX & MIN) OVER WHOLE HISTORY AT CU
     +RRENT DOFs'//30X,'ACCE.EXTREMA')

c 3200 FORMAT (//26X,'MAXIMUM VALUES OF AXLE LOAD UNDER EVERY WHEEL APPLI
c     +ED ON BRIDGE '/
c     +/21X,'WHEEL NO.',4X,'LATERAL',8X,'ROLLING',8X,'VERTICAL')
c	Previous changed to next 26Apr08
 3200 FORMAT (//22X,'MAXIMUM VALUES OF AXLE LOAD UNDER EVERY WHEEL APPLI
     +ED ON BRIDGE '/
     +/17X,'WHEEL NO.',4X,'LATERAL',8X,'ROLLING',8X,'VERTICAL MAX',8X,'VERTICAL MIN')

c 3500 FORMAT(23X,I5,3X,3(1X,E15.7))	!changed to next 26Apr08
 3500 FORMAT(19X,I5,3X,4(1X,E15.7))

 8410 FORMAT(27X,50(1X,E15.7)) !added 26Apr08

c 8420 FORMAT(10X,20(1X,E15.7)) !changed to next 03Mar07
 8420 FORMAT(23X,50(1X,E15.7)) 

 8476 FORMAT(//,
     &1X,'DO NOT CONVERGE AT THIS TIME = ',F10.4,5X,'DURING NMAX = ',I5/
     &1X,'CURRENT DISP DIFFERENCE =',E11.4,5X,'--->PLEASE INCRESE NMAX')
      END
C
C	=======================================================================

      SUBROUTINE CALFT3D(RT,PCAB,PBOG,PRWH,WRCO,WHPOS,IDWHPOS,FWRO,FT,
	1				  JTE,CS1,FSS) 	!added JTE 7Sep06, CS1 added 11Feb07
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO CALCULATE FORCE VECTOR APPLIED TO TRAIN FOR 3D-MODELLING
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	RSW(NWHEEL*3,2)= Response of structure part at wheel position [DWI,VWI] in bridge coordinate
C					  Where DWI = {Vertical,Rolling,Lateral} displacements of each wheel i.
C	PCAB(NCARB,4)  = Car body mechanical properties (MASS,Jx,Jz,Jy)
C	PBOG(NBOGI,4)  = Bogie of whole train (MASS,Jx,Jz,Jy)
C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jz,Jy,rW,b0)
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwZ,cwZ)
C	WHPOS(NWHEEL)  = Wheel position at current step (=WHPOSN)	
C	IDWHPOS(NWHEEL)= ID NUMBER OF ELEMENT CORRESPONDING TO WHEEL POSITION
C	NACAB   = Number of articulated car bodies = IT2
C	NCARB   = Total number of car bodies = IT1+IT2+IT3
C	NBOGI   = Total number of bogies     = IT1*2+(IT2+1)+IT3*2
C	NWHEEL  = Total number of wheels     = NBOGI*2
C	TLENGTH = Total length of path where wheel runs
C	SCALE   = Gravity acceleration (g)
C	NEQTR   = Number of equations of train
C	JTE	    = Iteration number, if ITE>1, not calculate rail Irre. but use previous at JTE=1
C	CS1(NWHEEL,2) = Constants for forces applied into wheel
C					 CS1(I,1) = 2*(kwV*bo**2 + kwH*rw**2)
C					 CS1(I,2) = 2*(cwV*bo**2 + cwH*rw**2)
C	OUTPUT:
C	-------
C	FWR(NWHEEL,3) = Wheel-rail contact forces {Lateral, Rolling, Vertical}
C	FT(NEQTR)     = Force vector applied to train (at current time)		
C	-------------------------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL
	COMMON /IR3O/ TIRLAT,TIRROT,TIRVER		!3D rail irre, added 6Sep06
	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM !added 31Jan2007 by Nguyen
	COMMON /IRW3/ IRRx,IRRe,IRRa,IRRr		 !added 31Jan2007 by Nguyen
	COMMON /MEMW/ W(7000000),IW(7000000)
C
      DIMENSION PCAB(NCARB,4),PBOG(NBOGI,4),PRWH(NWHEEL,6)
	DIMENSION WRCO(NWHEEL,6),WHPOS(NWHEEL),IDWHPOS(*) !!changed size of RSW 28Feb07
	DIMENSION FWR(NWHEEL,3),FT(NEQTR),CS1(NWHEEL,2)			! CS1 added 11Feb07
	
	DIMENSION RT(NEQTR,1),FWRO(NWHEEL,3),FSS(1)
C	-------------------------------------------------------------------------------------------

      
C	FORM WHEEL-RAIL CONTACT FORCES:
C	-------------------------------
	CALL CLEARMAT(FWR,NWHEEL,3)

      IPF=0   !Position pointer for each train DOFs

C	Front Power Cars:
C	-----------------
	IF(IT1.GT.0)THEN
       DO 10 I=1,IT1		 !Loop over number of power cars
C	  Power Car body:
	  IPF = IPF + 5
	   
C	  Two Bogies:
	  DO 10 KB=1,2
	    JB = 2*(I-1) + KB	 
	    IPF = IPF + 5

C	    Two wheels of each bogie:
	    DO 10 KW=1,2
		  JW = 2*(JB-1) + KW	

		  FWRO(JW,1) = 2.0*WRCO(JW,3)*(RT(IPF+1,1)-PRWH(JW,5)*RT(IPF+2,1)) 
	1				  +2.0*WRCO(JW,4)*(RT(IPF+1,2)-PRWH(JW,5)*RT(IPF+2,2))      !Lateral
		  FWRO(JW,2) = CS1(JW,1)*RT(IPF+2,1)
     1				  -2.0*WRCO(JW,3)*PRWH(JW,5)*RT(IPF+1,1)
     2				  +CS1(JW,2)*RT(IPF+2,2)
     3			      -2.0*WRCO(JW,4)*PRWH(JW,5)*RT(IPF+1,2)                    !Rolling			 
		  FWRO(JW,3) = 2.0*WRCO(JW,5)*RT(IPF+3,1)+2.0*WRCO(JW,6)*RT(IPF+3,2)    !Vertical

          IPF = IPF + 3
  10  CONTINUE
	ENDIF
C	Articulated Cars:
C	-----------------
	IF(IT2.GE.3)THEN	!for case of Articulated train
	DO 20 I = IT1+1,IT1+IT2
	    JB = IT1 + I		
C	    Two wheels of front bogie:
	    DO 15 KW=1,2
		  JW = 2*(JB-1) + KW	

		  FWRO(JW,1) = 2.0*WRCO(JW,3)*(RT(IPF+1,1)-PRWH(JW,5)*RT(IPF+2,1)) 
	1				  +2.0*WRCO(JW,4)*(RT(IPF+1,2)-PRWH(JW,5)*RT(IPF+2,2))      !Lateral
		  FWRO(JW,2) = 2.0*(WRCO(JW,5)*(PRWH(JW,6)**2.0)
	1				  +WRCO(JW,3)*(PRWH(JW,5)**2))*RT(IPF+2,1)
	2				  -2.0*WRCO(JW,3)*PRWH(JW,5)*RT(IPF+1,1)
	3				  +2.0*(WRCO(JW,6)*(PRWH(JW,6)**2.0)
	4				  +WRCO(JW,4)*(PRWH(JW,5)**2))*RT(IPF+2,2)
	5				  -2.0*WRCO(JW,4)*PRWH(JW,5)*RT(IPF+1,2)                    !Rolling		
		  FWRO(JW,3) = 2.0*WRCO(JW,5)*RT(IPF+3,1)+2.0*WRCO(JW,6)*RT(IPF+3,2)    !Vertical

	    IPF = IPF + 3
  15     CONTINUE
C		Front bogie 
	    IPF = IPF + 5
C		Car body:
	    IPF = IPF + 5
  20  CONTINUE

C	The last bogie of rear transition car:
C	--------------------------------------
	    JB = 2*IT1 + IT2 + 1			
C		Bogie: 
	    IPF = IPF + 5
C	    Two wheels of last bogie:
	    DO 25 KW=1,2
		  JW = 2*(JB-1) + KW	

		  FWRO(JW,1) = 2.0*WRCO(JW,3)*(RT(IPF+1,1)-PRWH(JW,5)*RT(IPF+2,1)) 
	1				  +2.0*WRCO(JW,4)*(RT(IPF+1,2)-PRWH(JW,5)*RT(IPF+2,2))      !Lateral
		  FWRO(JW,2) = 2.0*(WRCO(JW,5)*(PRWH(JW,6)**2.0)
	1				  +WRCO(JW,3)*(PRWH(JW,5)**2))*RT(IPF+2,1)
	2				  -2.0*WRCO(JW,3)*PRWH(JW,5)*RT(IPF+1,1)
	3				  +2.0*(WRCO(JW,6)*(PRWH(JW,6)**2.0)
	4				  +WRCO(JW,4)*(PRWH(JW,5)**2))*RT(IPF+2,2)
	5				  -2.0*WRCO(JW,4)*PRWH(JW,5)*RT(IPF+1,2)                    !Rolling		
		  FWRO(JW,3) = 2.0*WRCO(JW,5)*RT(IPF+3,1)+2.0*WRCO(JW,6)*RT(IPF+3,2)    !Vertical

	    IPF = IPF + 3
  25     CONTINUE
	ENDIF	!IF(IT2.GE.3), for case of classical train

C	Rear Power Cars:
C	-----------------
	IF(IT3.GT.0)THEN
       DO 30 I=IT1+IT2+1,IT1+IT2+IT3		!Loop over number of rear power cars
C	  Power Car body:
	  IPF = IPF + 5 
C	  Two Bogies:
	  DO 30 KB=1,2
		JB = 2*I-IT2-1+KB 	
	    IPF = IPF + 5
C	    Two wheels of each bogie:
	    DO 30 KW=1,2
		  JW = 2*(JB-1) + KW	

		  FWRO(JW,1) = 2.0*WRCO(JW,3)*(RT(IPF+1,1)-PRWH(JW,5)*RT(IPF+2,1)) 
	1				  +2.0*WRCO(JW,4)*(RT(IPF+1,2)-PRWH(JW,5)*RT(IPF+2,2))  !Lateral
		  FWRO(JW,2) = CS1(JW,1)*RT(IPF+2,1)
	2				  -2.0*WRCO(JW,3)*PRWH(JW,5)*RT(IPF+1,1)
	3				  +CS1(JW,2)*RT(IPF+2,2)
	5			      -2.0*WRCO(JW,4)*PRWH(JW,5)*RT(IPF+1,2)                !Rolling		
		  FWRO(JW,3) = 2.0*WRCO(JW,5)*RT(IPF+3,1)+2*WRCO(JW,6)*RT(IPF+3,2)  !Vertical										

	    IPF = IPF + 3
  30  CONTINUE	
	ENDIF
	
	
      ISS = 0
	DO JW = 1,NWHEEL
	  FWR(JW,1) = FSS(ISS+1)
	  FWR(JW,2) = FSS(ISS+2)
	  FWR(JW,3) = FSS(ISS+3)
	  ISS = ISS + 3
      ENDDO
      
      
C	---------------------------------------------------
C 	CALCULATE FORCE VECTOR APPLIED TO FRONT POWER CARS:                                                    
C	---------------------------------------------------
	CALL CLEARA(FT,NEQTR)
      IPF=0  ! Position pointer for force vector

	IF(IT1.GT.0)THEN !Case of No power car or all cars are articulated
      DO 100 I=1,IT1		!Loop over number of power cars
C	  Power Car body:	
C	  ---------------
	  FT(IPF+1) = 0.0
	  FT(IPF+2) = 0.0
	  FT(IPF+3) = 0.0
c	  FT(IPF+4) =-PCAB(I,1)*SCALE	!Changed to next 17June07 (Gravity force applied to bridge) 
	  FT(IPF+4) = 0.0
	  FT(IPF+5) = 0.0
	  IPF = IPF + 5
	 
C	  Two Bogies:	
C	  -----------
	  DO 100 KB=1,2
		  JB = 2*(I-1) + KB	!Pointer for bogie of the power car  
	    FT(IPF+1) = 0.0
	    FT(IPF+2) = 0.0
	    FT(IPF+3) = 0.0
C	    FT(IPF+4) =-PBOG(JB,1)*SCALE  !Changed to next 17June07 (Gravity force applied to bridge) 
	    FT(IPF+4) = 0.0
	    FT(IPF+5) = 0.0

	    IPF = IPF + 5

C	    Two wheels of each bogie:
C	    -------------------------
	    DO 100 KW=1,2
		    JW = 2*(JB-1) + KW	!Pointer for wheel of current bogie 
	      FT(IPF+1) = FWR(JW,1)
	      FT(IPF+2) = FWR(JW,2)
c	      FT(IPF+3) =-PRWH(JW,1)*SCALE + FWR(JW,3)	 !Changed to next 17June07 (Gravity force applied to bridge)
	      FT(IPF+3) = FWR(JW,3)

	      IPF = IPF + 3

  100 CONTINUE
	ENDIF
C	---------------------------------------------------
C 	CALCULATE FORCE VECTOR APPLIED TO ARTICULATED CARS:                                                    
C	---------------------------------------------------
	IF(IT2.GE.3)THEN	!for case of articulated train
	DO 200 I = IT1+1,IT1+IT2
	    JB = IT1 + I			!Pointer for bogie   
C	    Two wheels of front bogie:
C	    -------------------------
	    DO 150 KW=1,2
		  JW = 2*(JB-1) + KW	!Pointer for wheel of current bogie 
	      FT(IPF+1) = FWR(JW,1)
	      FT(IPF+2) = FWR(JW,2)
c	      FT(IPF+3) =-PRWH(JW,1)*SCALE + FWR(JW,3)	!Changed to next 17June07 (Gravity force applied to bridge) 
	      FT(IPF+3) = FWR(JW,3)
	      IPF = IPF + 3
  150     CONTINUE

C		Front bogie 
C	    ------------
	    FT(IPF+1) = 0.0
	    FT(IPF+2) = 0.0
	    FT(IPF+3) = 0.0
c	    FT(IPF+4) =-PBOG(JB,1)*SCALE	!Changed to next 17June07 (Gravity force applied to bridge)
	    FT(IPF+4) = 0.0
	    FT(IPF+5) = 0.0
	    IPF = IPF + 5

C		Car body:
C		---------
	    FT(IPF+1) = 0.0
	    FT(IPF+2) = 0.0
	    FT(IPF+3) = 0.0
c	    FT(IPF+4) =-PCAB(I,1)*SCALE	   !Changed to next 17June07 (Gravity force applied to bridge)
	    FT(IPF+4) = 0.0
	    FT(IPF+5) = 0.0
	    IPF = IPF + 5

  200 CONTINUE

C	The last bogie of rear transition car:
C	--------------------------------------
	    JB = 2*IT1 + IT2 + 1			!Pointer for bogie  
C		Bogie: 
C	    ------
	    FT(IPF+1) = 0.0
	    FT(IPF+2) = 0.0
	    FT(IPF+3) = 0.0
c	    FT(IPF+4) =-PBOG(JB,1)*SCALE	!Changed to next 17June07 (Gravity force applied to bridge)
	    FT(IPF+4) = 0.0
	    FT(IPF+5) = 0.0
	    IPF = IPF + 5

C	    Two wheels of last bogie:
C	    -------------------------
	    DO 250 KW=1,2
		  JW = 2*(JB-1) + KW	!Pointer for wheel of current bogie 
	      FT(IPF+1) = FWR(JW,1)
	      FT(IPF+2) = FWR(JW,2)
c	      FT(IPF+3) =-PRWH(JW,1)*SCALE + FWR(JW,3)	   !Changed to next 17June07 (Gravity force applied to bridge)
	      FT(IPF+3) = FWR(JW,3)
	      IPF = IPF + 3
  250     CONTINUE

	ENDIF	!IF(IT2.GE.3), for case of articulated train
C	--------------------------------------------------
C 	CALCULATE FORCE VECTOR APPLIED TO REAR POWER CARS:                                                    
C	--------------------------------------------------
	IF(IT3.GT.0)THEN
      DO 300 I=IT1+IT2+1,IT1+IT2+IT3		!Loop over number of rear power cars
C	  Power Car body:
C	  ---------------
	  FT(IPF+1) = 0.0
	  FT(IPF+2) = 0.0
	  FT(IPF+3) = 0.0
c	  FT(IPF+4) =-PCAB(I,1)*SCALE	!Changed to next 17June07
	  FT(IPF+4) = 0.0
	  FT(IPF+5) = 0.0
	  IPF = IPF + 5
	 
C	  Two Bogies:
C	  -----------
	  DO 300 KB=1,2
		JB = 2*I-IT2-1+KB 	 
	    FT(IPF+1) = 0.0
	    FT(IPF+2) = 0.0
	    FT(IPF+3) = 0.0
c	    FT(IPF+4) =-PBOG(JB,1)*SCALE	!Changed next 17June07
	    FT(IPF+4) = 0.0
	    FT(IPF+5) = 0.0
	    IPF = IPF + 5

C	    Two wheels of each bogie:
C	    -------------------------
	    DO 300 KW=1,2
		  JW = 2*(JB-1) + KW	!Pointer for wheel of current bogie 
	      FT(IPF+1) = FWR(JW,1)
	      FT(IPF+2) = FWR(JW,2)
c	      FT(IPF+3) =-PRWH(JW,1)*SCALE + FWR(JW,3)	!Changed next 17June07
	      FT(IPF+3) = FWR(JW,3)
	      IPF = IPF + 3

  300 CONTINUE
	ENDIF
C
      RETURN
      END
C
C	=======================================================================
      SUBROUTINE CALFS3D(FSSO,FS,RSW,FWR,IDWHPOS,NEID,PELW,XI,ID,
     +                   PRWH,WRCO,CS1,FAXL,TIRAIL,FWCON)	
C			  !CS1 added 11Feb07, FAXL added 20June07
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO: 
C	 - LUMPED FORCE APPLIED TO STRUCTURE AT WHEEL POSITION {Fyi,FRxi,Fzi}
C	 - CALCULATE EQUIVALENT NODAL FORCE VECTOR APPLIED TO 3D-STRUCTURE {FS(NEQ)}
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	RT(NEQTR,3)     = Response of Train part at train's DOFs [DTI,VTI,ATI]
C					Where DTI = {Yi,Rxi,Rzi,Zi,Ryi}T ={Lateral,Rolling,Yawing,Vertical,Pitching} 
C								   displacements of each rigid body i (car body or bogie)
C								={Yi,Rxi,Zi}T ={Lateral,Rolling,Vertical} for each wheel-axle set
C	FWR(NWHEEL,3) = Wheel-rail contact forces {Lateral, Rolling, Vertical}

C	NSF    = Maximum number of DOFs per node (=6 or 7, depends on the type of element)
C	NNF    = NUMBER OF NODAL DEGREES OF FREEDOM
C	NEQ    = NUMBER OF EQUATIONS OF STRUCTURE
C	NNO    = CURRENT NUMBER OF ELEMENT NODES

C	H4    = Verical distance from mass center of the wheel set to shear center of bridge
C	ECC  = Eccentricity of the track: lateral distance from track center to bridge center

C	ID(NSF,NSN)  = BOUNDARY CODES AND EQUATION NUMBERSNRL = NRL(NSF,NSN) = EQUATION NUMBER
C	IDWHPOS(NWHEEL) = ID NUMBER OF ELEMENT CORRESPONDING TO WHEEL POSITION
C	NEID(2,NELW) = NODE NUMBER OF BEAM ELEMENT WHERE WHEEL RUNS
C	XI(NWHEEL)      = Natural coordinate of wheel position at the current element
C					      XI(I)= DABS(WHPOS(I)-XYZI)/DABS(XYZI-XYZJ)
C	NWHEEL= Total number of wheels     = NBOGI*2
C	NELW   = NUMBER OF SEGMENT WHERE WHEEL RUNS

C	PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jz,Jy,rW,b0)
C	WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwZ,cwZ)
C	WPOS(NWHEEL)   = Wheel position
C	CS1(NWHEEL,2)  = Constants for forces applied into wheel
C					 CS1(I,1) = kwV*bo**2 + kwH*rw**2	
C					 CS1(I,2) = cwV*bo**2 + cwH*rw**2	
C	FAXL(NWHEEL)   = Static axle force (vertical force due to train weights acting onto bridge)
c	-------
C	OUTPUT:
C	-------
C	FS(NEQ)       = Equivalent nodal force vector applied to structure (at current time)
C	FSS(NWHEEL*3) = Lumped force applied to structure at wheel position {Fyi,FRxi,Fzi}

C	INDEX:
C	------
C	SFBC0(2)= shape functions of C0-beam element in narural coordinate {1/2(1-r) 1/2(1+r)}
C	SFBC1(6)= shape functions of C1-beam element
C	FSE(14) = Equivalent nodal forces vector of the current element (from lumped forces at wheel position)
C	-------------------------------------------------------------------------------
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /DMCO/ IDMA,IANA,IDSO,NNBR,NNRA,NEBEBR,NPBEBR,NLRB,
     +			  NEWKTK,NPWKTK,NEBERA,NPBERA,NDOF,NNBE,NNWK
	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	COMMON /NIAX/ NWHEEL
	
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
C
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV	
	
      DIMENSION FS(NEQ),FSS(NWHEEL*3),FWR(NWHEEL,3)
	DIMENSION NEID(2,1),PELW(4,1),XI(1),ID(NSF,1),IDWHPOS(NWHEEL)
      DIMENSION PRWH(NWHEEL,6),WRCO(NWHEEL,6)
      DIMENSION CS1(NWHEEL,2)		!CS1 added 11Feb07
	DIMENSION FAXL(NWHEEL)		!added 20June07
	
	DIMENSION FSSO(NWHEEL*3),RSW(NWHEEL*3,2),TIRAIL(3,NWHEEL)
	
	
	DIMENSION VR(3),VS(3),VT(3),TRANS(14,14),FIXD(14),FIXG(14),FSE(7,2)
	DIMENSION FWCON(NWHEEL,3)
C	-----------------------------------------------------------------------------------------
C	CALCULATE TRAIN/RAIL INTERACTION LUMPED FORCE APPLIED TO STRUCTURE AT WHEEL POSITION{FSS}
C	-----------------------------------------------------------------------------------------

	ISS=0	!Pointer for Wheel disp & forces components applied on bridge at wheel position
	DO 100 IWH = 1,NWHEEL

        IF(IDWHPOS(IWH).GT.0) THEN
        
            FSSO(ISS+1) = 2.0*WRCO(IWH,3)*( RSW(IWH*3,1) + TIRAIL(3,IWH)
     1                                    - PRWH(IWH,5)*(RSW(IWH*3-1,1)+TIRAIL(2,IWH)) )
     2                   +2.0*WRCO(IWH,4)*( RSW(IWH*3,2) - PRWH(IWH,5)*RSW(IWH*3-1,2) )   !Lateral contact force

	      FSSO(ISS+2) = CS1(IWH,1)*( RSW(IWH*3-1,1) + TIRAIL(2,IWH) )					
	1			         -2.0*WRCO(IWH,3)*PRWH(IWH,5)*( RSW(IWH*3,1) + TIRAIL(3,IWH) )
     2			         +CS1(IWH,2)*RSW(IWH*3-1,2)
	3		             -2.0*WRCO(IWH,4)*PRWH(IWH,5)*RSW(IWH*3,2)                        !Rolling contact moment

	      FSSO(ISS+3) = 2.0*WRCO(IWH,5)*( RSW(IWH*3-2,1) + TIRAIL(1,IWH) )			    
	1		 	         +2.0*WRCO(IWH,6)*RSW(IWH*3-2,2)                                  !Vertical contact force
	
	      FSS(ISS+1) = FWR(IWH,1) - FSSO(ISS+1) 
	      FSS(ISS+2) = FWR(IWH,2) - FSSO(ISS+2) 
	      FSS(ISS+3) = FWR(IWH,3) - FSSO(ISS+3) - FAXL(IWH)
	      
	      IF(FSS(ISS+3).GT.0.0D0) FSS(ISS+3) = 0.0D0 !NO TENSION FORCE ON RAIL
	      
	  ELSE		!(the wheel is running on rail without bridge)
	  
C	    since without bridge, rail Irre still creates interaction forces		
	      FSSO(ISS+1) = 2.0*WRCO(IWH,3)*( TIRAIL(3,IWH) - PRWH(IWH,5)*TIRAIL(2,IWH) )

	      FSSO(ISS+2) = CS1(IWH,1)*TIRAIL(2,IWH)
     1                   -2.0*WRCO(IWH,3)*PRWH(IWH,5)*TIRAIL(3,IWH) 
          		   	
	      FSSO(ISS+3) = 2.0*WRCO(IWH,5)*TIRAIL(1,IWH)
		
	      FSS(ISS+1) = 0.0D0 
	      FSS(ISS+2) = 0.0D0 
	      FSS(ISS+3) = 0.0D0
	      
	  ENDIF

C     FOR PRINTING MAXIMUM CONTACT FORCES SEE ALSO "TROUT3D"	  
	  FWCON(IWH,1) = FWR(IWH,1) - FSSO(ISS+1)
	  FWCON(IWH,2) = FWR(IWH,2) - FSSO(ISS+2)
	  FWCON(IWH,3) = FWR(IWH,3) - FSSO(ISS+3) - FAXL(IWH)
	  IF(FWCON(IWH,3).GT.0.0D0) FWCON(IWH,3) = 0.0D0 !NO TENSION FORCE ON RAIL
	  
		ISS = ISS + 3
			
100   CONTINUE
				
C	-----------------------------------------------------------------------------------
C	EQUIVALENT NODAL FORCE CORRESPONDING TO INTERACTION LUMPED FORCE AT WHEEL POSITION
C	-----------------------------------------------------------------------------------
	CALL CLEARA(FS,NEQ)

      DO 40 I=1,NWHEEL
        IF(IDWHPOS(I).GT.0) THEN

	    ISS = (I-1)*3   

C	CALCULATE EQUIVALENT NODAL FORCES VECTOR OF THE CURRENT BEAM ELEMENT
            IENUM=IDWHPOS(I)		!ID of current element
            
C     LOCAL VECTOR & SEGMENT LENGTH		
            VR(1:3) = PELW(1:3,IENUM)
               SLEN = PELW(  4,IENUM)
               
            CALL FMVEVR(VR,VS,VT)
            RANG = 0.0D0
	      CALL ROMBAC(VR,VS,VT,RANG)  
	      CALL TRANLG(VR,VS,VT,TRANS)
    	
	      WR = VR(NGRAV)*FSS(ISS+3)
	      WS = VS(NGRAV)*FSS(ISS+3)
	      WT = VT(NGRAV)*FSS(ISS+3)+FSS(ISS+1)
	      MR = FSS(ISS+2) + (-WS)*ECC
    	        
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

 1000 RETURN
      END

C	=======================================================================	
      SUBROUTINE SSPACE1(A,B,IDIAG,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,
     &  BLO,BUPC,NN,NNM,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,
     &  ITO,ISO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------------
C	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
C     1              IFPR(10),IFPL(10)
CMKIM 2000. 3. 7
      DIMENSION A(NWK),B(*),R(NN,*),TT(*),W(*),EIGV(3),D(*),VEC(NC,*),
     &          AR(*),BR(*),RTOLV(*),BUP(*),BLO(1),BUPC(1)
CN      INTEGER   IDIAG(*)
	DIMENSION IDIAG(*)
C	-------------------------------------------------------------------
CN      write(90,8720) (idiag(i),i=1,nnm)
CN	WRITE(ISO,*)'VECTOR CONTAIN ADD. OF DIAG. ELEMS OF STIFF MATRIX A'
CN      write(ISO,8720) (idiag(i),i=1,nnm)
CN 8720 format(20i5)

CN	Added 18 Nov04
CN	EIGV(1)= 0.67217433832162
CN	EIGV(2)= 0.10576561553736E+01 
CN	EIGV(3)= 0.13192847760328E+04
CN	return
c
C9Dec  	WRITE (ISO,*) ('DIA STIFF MATRIX OF BRIDGE AFTER REARANGE AKBO')
C9Dec	DO I = 1,NN
C9Dec	II = IDIAG(I)
C9Dec	WRITE (ISO,*) A(II)
C9Dec	ENDDO

C	SET TOLERANCE FOR JACOBI ITERATION
C	----------------------------------
      TOLJ=1.D-10
C
C	INITIALIZATION
C	---------------
      ICONV=0
      NSCH=0
      NSMAX=12
      N1=NC+1
      NC1=NC-1
      REWIND NSTIF
      WRITE(NSTIF) A
      DO 60 I=1,NC
   60 D(I)=0.
C
C	ESTABLISH STARTING ITERATION VECTORS
C	-------------------------------------
      ND = NN/NC
      IF(NWM.GT.NN) GO TO 4
      J=0
      DO 2 I=1,NN
      II=IDIAG(I)
      R(I,1) = B(I)
      IF(B(I).GT.0.D0) J=J+1
    2 W(I)=B(I)/A(II)
      IF(NC.LE.J) GO TO 16
CN      WRITE(IOUT,1007)
CN	WRITE(ITO,1007)
      WRITE(ISO,1007)
      STOP
    4 DO 10 I=1,NN
      II = IDIAG(I)
      R(I,1) = B(II)
   10 W(I) = B(II)/A(II)
   16 DO 20 J=2,NC
      DO 20 I=1,NN
   20 R(I,J)=0.
C
      L=NN-ND
      DO 30 J=2,NC
      RT=0.
      DO 40 I=1,L
      IF(W(I).LT.RT) GO TO 40
      RT=W(I)
      IJ=I
   40 CONTINUE
      DO 50 I=L,NN
      IF(W(I).LE.RT) GO TO 50
      RT=W(I)
      IJ=I
   50 CONTINUE
      TT(J)=FLOAT(IJ)
      W(IJ)=0.
      L=L-ND
   30 R(IJ,J)=1.
CN      WRITE(IOUT,1008)
CN      WRITE(IOUT,1002) (TT(J),J=2,NC)
      WRITE(ISO,1008)
      WRITE(ISO,1002) (TT(J),J=2,NC)
C
C	FACTORIZE MATRIX A INTO (L)*(D)*(L(T)):
C	----------------------------------------
      ISH=0
CN      CALL DECOMP(A,IDIAG,NN,ISH,IOUT)
      CALL DECOMP(A,IDIAG,NN,ISH,ISO)
C	-----------------------------------------------
C	S T A R T    O F   I T E R A T I O N   L O O P:
C	-----------------------------------------------
      NITE=0
  100 NITE=NITE+1
      IF(IFPR.EQ.0) GO TO 90
C      WRITE(IOUT,1010) NITE
      WRITE(ISO,1010) NITE
C
C	CALCULATE THE PROJECTION OF A AND B:
C	------------------------------------
   90 IJ=0
      DO 110 J=1,NC
      DO 120 K=1,NN
  120 TT(K)=R(K,J)
      CALL REDBAK(A,TT,IDIAG,NN)
      DO 130 I=J,NC
      ART=0.
      DO 140 K=1,NN
  140 ART=ART+R(K,I)*TT(K)
      IJ=IJ+1
  130 AR(IJ)=ART
      DO 150 K=1,NN
  150 R(K,J)=TT(K)
  110 CONTINUE
      IJ=0
      DO 160 J=1,NC
      CALL MULTDM(TT,B,R(1,J),IDIAG,NN,NWM)
      DO 180 I=J,NC
      BRT=0.
      DO 190 K=1,NN
  190 BRT=BRT+R(K,I)*TT(K)
      IJ=IJ+1
  180 BR(IJ)=BRT
      IF(ICONV.GT.0) GO TO 160 
      DO 200 K=1,NN
  200 R(K,J)=TT(K)
  160 CONTINUE
C
C	SOLVE FOR EIGENSYSTEM OF SUBSPACE OPERATORS:
C	--------------------------------------------
      IF(IFPR.EQ.0) GO TO 320
      IND=1
CN  210 WRITE(IOUT,1020)
  210 WRITE(ISO,1020)
      II=1
      DO 300 I=1,NC
      ITEMP=II+NC-I

CN      WRITE(IOUT,1005) (AR(J),J=II,ITEMP)
	WRITE(ISO,1005) (AR(J),J=II,ITEMP)
  300 II=II+N1-I
CN       WRITE(IOUT,1030)
	 WRITE(ISO,1030)
      II=1
      DO 310 I=1,NC
      ITEMP=II+NC-I

CN      WRITE(IOUT,1005) (BR(J),J=II,ITEMP)
	WRITE(ISO,1005) (BR(J),J=II,ITEMP)
  310 II=II+N1-I
      IF(IND.EQ.2) GO TO 350
C
CN  320 CALL JACOBI1(AR,BR,VEC,EIGV,W,NC,NNC,TOLJ,NSMAX,IFPR,IOUT)
  320 CALL JACOBI1(AR,BR,VEC,EIGV,W,NC,NNC,TOLJ,NSMAX,IFPR,ISO)
C
      IF(IFPR.EQ.0) GO TO 350
CN      WRITE(IOUT,1040)
      WRITE(ISO,1040)
      IND=2
      GO TO 210
C
C	ARRANGE EIGENVALUE IN ASCENDING ORDER:
C	--------------------------------------
  350 IS=0
      II=1
      DO 360 I=1,NC1
      ITEMP=II+N1-I
      IF(EIGV(I+1).GE.EIGV(I)) GO TO 360
      IS=IS+1
      EIGVT=EIGV(I+1)
      EIGV(I+1)=EIGV(I)
      EIGV(I)=EIGVT
      BT=BR(ITEMP)
      BR(ITEMP)=BR(II)
      BR(II)=BT
      DO 370 K=1,NC
      RT=VEC(K,I+1)
      VEC(K,I+1)=VEC(K,I)
  370 VEC(K,I)=RT
  360 II=ITEMP
      IF(IS.GT.0) GO TO 350
      IF(IFPR.EQ.0) GO TO 375

C      WRITE(IOUT,1035)
      WRITE(ISO,1035)
CN      WRITE(IOUT,1006) (EIGV(I),I=1,NC)
      WRITE(ISO,1006) (EIGV(I),I=1,NC)
C
C	CALCULATE B TIMES APROXIMATE EIGENVECTORS (ICONV.EQ.0)
C	OR FINAL EIGENVECTORS APROXIMATIONS (ICONV.GT.0)
C	------------------------------------------------------
  375 DO 420 I=1,NN
      DO 422 J=1,NC
  422 TT(J)=R(I,J)
      DO 424 K=1,NC
      RT=0.
      DO 430 L=1,NC
  430 RT=RT+TT(L)*VEC(L,K)
  424 R(I,K)=RT
  420 CONTINUE
      IF(ICONV.GT.0) GO TO 500
C
C	CHECK CONVERGENCE OF EIGENVALUES:
C	---------------------------------
      DO 380 I=1,NC
      DIF=DABS(EIGV(I)-D(I))
  380 RTOLV(I)=DIF/EIGV(I)
      IF(IFPR.EQ.0) GO TO 385

CN      WRITE(IOUT,1050)
      WRITE(ISO,1050)

CN      WRITE(IOUT,1005) (RTOLV(I),I=1,NC)
	WRITE(ISO,1005) (RTOLV(I),I=1,NC)
C
  385 DO 390 I=1,NROOT
      IF(RTOLV(I).GT.RTOL) GO TO 400
  390 CONTINUE
CN      WRITE(IOUT,1060) RTOL
      WRITE(ISO,1060) RTOL
      ICONV=1
      GO TO 100
  400 IF(NITE.LT.NITEM) GO TO 410
CN      WRITE(IOUT,1070)
      WRITE(ISO,1070)
      ICONV=2
      IFSS=0
      GO TO 100
C
  410 DO 440 I=1,NC
  440 D(I)=EIGV(I)
      GO TO 100
C
C	<====E N D    O F   I T E R A T I O N   L O O P
C
CN  500 WRITE(IOUT,1100)
CN      WRITE(IOUT,1006) (EIGV(I),I=1,NROOT)
  500 WRITE(ISO,1100)
      WRITE(ISO,1006) (EIGV(I),I=1,NROOT)

CMKIM 2000. 2. 17.
C	CALCULATE AND PRINT ERROR NORMS
C	--------------------------------
      REWIND NSTIF
      READ(NSTIF) A
      DO 580 L=1,NROOT
      RT=EIGV(L)
      CALL MULTDM(TT,A,R(1,L),IDIAG,NN,NWK)
      VNORM=0.
      DO 590 I=1,NN
  590 VNORM=VNORM+TT(I)*TT(I)
      CALL MULTDM(W,B,R(1,L),IDIAG,NN,NWM)
      WNORM=0.
      DO 600 I=1,NN
      TT(I) = TT(I) - RT*W(I)
  600 WNORM=WNORM+TT(I)*TT(I)
      VNORM=DSQRT(VNORM)
      WNORM=DSQRT(WNORM)
      D(L)=WNORM/VNORM
  580 CONTINUE
C      WRITE(IOUT,1115)
C      WRITE(IOUT,1006) (D(I),I=1,NROOT)
      WRITE(ISO,1115)
      WRITE(ISO,1006) (D(I),I=1,NROOT)
C
C	APPLY STURM SEQUENCE CHECK
C	----------------------------
      IF(IFSS.EQ.0) GO TO 700
CN      CALL SCHECK1(EIGV,RTOLV,BUP,BLO,BUPC,D,NC,NEI,RTOL,SHIFT,IOUT)
      CALL SCHECK1(EIGV,RTOLV,BUP,BLO,BUPC,D,NC,NEI,RTOL,SHIFT,ISO)
C
CN      WRITE(IOUT,1120) SHIFT
      WRITE(ISO,1120) SHIFT
C
C	SHIFT MATRIX A
C	---------------
      REWIND NSTIF
      READ(NSTIF) A
      IF(NWM.GT.NN) GO TO 645
      DO 640 I=1,NN
      II=IDIAG(I)
  640 A(II)=A(II)-B(I)*SHIFT
      GO TO 660
  645 DO 650 I=1,NWK
  650 A(I)=A(I)-B(I)*SHIFT
C
C	FACTORIZE SHIFT MATRIX
C	----------------------
  660 ISH=1
CN      CALL DECOMP(A,IDIAG,NN,ISH,IOUT)
      CALL DECOMP(A,IDIAG,NN,ISH,ISO)

C	COUNT OF NEGATIVE DIGONAL ELEMENTS
      NSCH=0
      DO 664 I=1,NN
      II=IDIAG(I)
      IF(A(II).LT.0.D0) NSCH=NSCH+1
  664 CONTINUE
      IF(NSCH.EQ.NEI) GO TO 670
      NMIS=NSCH-NEI
CN      WRITE(IOUT,1130) NMIS
      WRITE(ISO,1130) NMIS
      GO TO 700
CN  670 WRITE(IOUT,1140) NSCH
  670 WRITE(ISO,1140) NSCH
  700 RETURN
C
CMKIM 2001. 8. 15
 1002 FORMAT(7F10.0)
C1002 FORMAT(1X,10F10.0)
 1005 FORMAT(1X,15F7.4)
CMKIM 2001. 8. 15.
C1006 FORMAT(5(2X,E13.6))
 1006 FORMAT(1X,6E22.14)
 1007 FORMAT(/1X,'STOP, NC IS LARGER THAN THE NUMBER OF MASS ',
     &  'DEGREES OF FREDOM')
 1008 FORMAT(//1X,'DREGREES OF FREDOM EXCITED BY UNIT STARTING ',
     &  'ITERATION VECTORS')
 1010 FORMAT(//,1X,'I T E R A T I O N    N U M B E R ',I4/)
 1020 FORMAT(/1X,'PROJECTION OF A (MATRIX AR)')
 1030 FORMAT(/1X,'PROJECTION OF B (MATRIX BR)')
 1035 FORMAT(/1X,'EIGENVALUE OF AR-LAMBDA*BR')
 1040 FORMAT(/1X,'AR AND BR AFTER JACOBI DIGONALIZATION')
 1050 FORMAT(/1X,'RELATIVE TOLERANCE REACHED EIGENVALUES')
 1060 FORMAT(/1X,'CONVERGENCE REACHED FOR RTOL',E11.4/)
 1070 FORMAT(//,'NO CONVERGENCE IN MAXIMUM NUMBER OF ',
     &  'ITERATIONS PERMITED'/1X,'WE ACCEPT CURRENT ',
     &  'ITERATION VLAUE'/1X,'THE STURM SEQUENCE CHECK ',
     &  'IS NOT PERFOMED')
CMKIM 2000. 2. 17
C1100 FORMAT(//1X,'THE SQUARE ROOTS OF CALCULATD EIGENVALUES ARE',/)
 1100 FORMAT(//1X,'THE CALCULATD EIGENVALUES ARE',/)
 1115 FORMAT(//1X,'PRINT ERROR NORMS ON THE EIGENVALUES')
 1110 FORMAT(//1X,'THE CALCULATED EIGENVECTORS ARE'//)
 1120 FORMAT(//1X,'CHECK APPLIED AT SHIFT',E22.14)
 1130 FORMAT(//1X,'THERE ARE',I4,3X,'EIGENVALUES MISSING')
 1140 FORMAT(//1X,'WE FOUND THE LOWEST',I4,3X,'EIGENVALUES')
      END
C
C	=======================================================================
      SUBROUTINE JACOBI1(A,B,X,EIGV,D,N,NWA,RTOL,NSMAX,IFPR,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------                                                            
C     PROGRAM  TO SOLVE THE GENERALIZED EIGENPROGRAM USING THE GENERALIZED 
C			 JACOBI ITERATION                            
C	-----------------------------------------------------------------------    
C     DIMENSION A(NWA),B(NWA),X(N,N),EIGV(N),D(N)
      DIMENSION A(*),B(*),X(N,*),EIGV(*),D(*)
C
C     INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES
C
      N1=N+1
      II=1
      DO 10 I=1,N
      IF(A(II).GT.0.D0.AND.B(II).GT.0.D0) GO TO 4
CN      IF(A(II).GT.0.D0.OR.B(II).GT.0.D0) GO TO 4
      WRITE(IOUT,2020) II,A(II),B(II)
      STOP
    4 D(I)=A(II)/B(II)
      EIGV(I)=D(I)
   10 II=II+N1-I
      DO 30 I=1,N
      DO 20 J=1,N
   20 X(I,J)=0.
   30 X(I,I)=1.
      IF(N.EQ.1) RETURN
C
C     INITIALIZE SWEEP COUNTER AND BEGIN ITERATION
C
      NSWEEP=0
      NR=N-1
   40 NSWEEP=NSWEEP+1
      IF(IFPR.EQ.1) WRITE(IOUT,2000) NSWEEP
C
C     CHECK IF PRESENT OFF-DIAGONAL ELEMENT IS LARGE ENOUGH TO
C     REQUIRE ZEROING 
C
      EPS=(.01**NSWEEP)**2
      DO 210 J=1,NR
      JP1=J+1
      JM1=J-1
      LJK=JM1*N-JM1*J/2
      JJ=LJK+J
      DO 210 K=JP1,N
      KP1=K+1
      KM1=K-1
      JK=LJK+K
      KK=KM1*N-KM1*K/2+K
      EPTOLA=(A(JK)*A(JK))/(A(JJ)*A(KK))
      EPTOLB=(B(JK)*B(JK))/(B(JJ)*B(KK))
      IF((EPTOLA.LT.EPS).AND.(EPTOLB.LT.EPS)) GO TO 210
C
C     IF ZEROING IS REQUIRED , CALCULATE THE ROTATION MATRIX
C
C     ELEMENTS CA AND CG.
      AKK=A(KK)*B(JK)-B(KK)*A(JK)
      AJJ=A(JJ)*B(JK)-B(JJ)*A(JK)
      AB=A(JJ)*B(KK)-A(KK)*B(JJ)
      CHECK=(AB*AB+4.*AKK*AJJ)/4.
      IF(CHECK) 50,60,60
   50 WRITE(IOUT,2020)
      STOP
   60 SQCH=DSQRT(CHECK)
      D1=AB/2.+SQCH
      D2=AB/2.-SQCH
      DEN=D1
      IF(DABS(D2).GT.DABS(D1)) DEN=D2
      IF(DEN) 80,70,80
   70 CA=0.
      CG=-A(JK)/A(KK)
      GO TO 90
   80 CA=AKK/DEN
      CG=-AJJ/DEN
C
C     PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT
C     OFF-DIAGONAL ELEMENT
   90 IF(N-2)100,190,100
  100 IF(JM1-1)130,110,110
  110 DO 120 I=1,JM1
      IM1=I-1
      IJ=IM1*N-IM1*I/2+J
      IK=IM1*N-IM1*I/2+K
      AJ=A(IJ)
      BJ=B(IJ)
      AK=A(IK)
      BK=B(IK)
      A(IJ)=AJ+CG*AK
      B(IJ)=BJ+CG*BK
      A(IK)=AK+CA*AJ
  120 B(IK)=BK+CA*BJ
  130 IF(KP1-N)140,140,160
  140 LJI=JM1*N-JM1*J/2
      LKI=KM1*N-KM1*K/2
      DO 150 I=KP1,N
      JI=LJI+I
      KI=LKI+I
      AJ=A(JI)
      BJ=B(JI)
      AK=A(KI)
      BK=B(KI)
      A(JI)=AJ+CG*AK
      B(JI)=BJ+CG*BK
      A(KI)=AK+CA*AJ
  150 B(KI)=BK+CA*BJ
  160 IF(JP1-KM1) 170,170,190
  170 LJI=JM1*N-JM1*J/2
      DO 180 I=JP1,KM1
      JI=LJI+I
      IM1=I-1
      IK=IM1*N-IM1*I/2+K
      AJ=A(JI)
      BJ=B(JI)
      AK=A(IK)
      BK=B(IK)
      A(JI)=AJ+CG*AK
      B(JI)=BJ+CG*BK
      A(IK)=AK+CA*AJ
  180 B(IK)=BK+CA*BJ
  190 AK=A(KK)
      BK=B(KK)
      A(KK)=AK+2.*CA*A(JK)+CA*CA*A(JJ)
      B(KK)=BK+2.*CA*B(JK)+CA*CA*B(JJ)
      A(JJ)=A(JJ)+2.*CG*A(JK)+CG*CG*AK
      B(JJ)=B(JJ)+2.*CG*B(JK)+CG*CG*BK
      A(JK)=0.
      B(JK)=0.
C
C	UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION
C
      DO 200 I=1,N
      XJ=X(I,J)
      XK=X(I,K)
      X(I,J)=XJ+CG*XK
  200 X(I,K)=XK+CA*XJ
  210 CONTINUE
C
C	UPDATE THE EIGENVALUES AFTER EACH SWEEP
C
      II=1
      DO 220 I=1,N
      IF(A(II).GT.0.D0.AND.B(II).GT.0.D0) GO TO 215
CN      IF(A(II).GT.0.D0.OR.B(II).GT.0.D0) GO TO 215
      WRITE(IOUT,2020)II,A(II),B(II)
      STOP
  215 EIGV(I)=A(II)/B(II)
  220 II=II+N1-I
      IF(IFPR.EQ.0) GO TO 230
      WRITE(IOUT,2030)
      WRITE(IOUT,2010)(EIGV(I),I=1,N)
C	CHECK FOR CONVERGENCE
  230 DO 240 I=1,N
      TOL=RTOL*D(I)
      DIF=DABS(EIGV(I)-D(I))
      IF(DIF.GT.TOL) GO TO 280
  240 CONTINUE
C
C     CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP IS
C        REQUIRED
      EPS=RTOL**2
      DO 250 J=1,NR
      JM1=J-1
      JP1=J+1
      LJK=JM1*N-JM1*J/2
      JJ=LJK+J
      DO 250 K=JP1,N
      KM1=K-1
      JK=LJK+K
      KK=KM1*N-KM1*K/2+K
      EPSA=(A(JK)*A(JK))/(A(JJ)*A(KK))
      EPSB=(B(JK)*B(JK))/(B(JJ)*B(KK))
      IF((EPSA.LT.EPS).AND.(EPSB.LT.EPS)) GO TO 250
      GO TO 280
  250 CONTINUE
C
C	FILL OUT BOTTOM TRIANGLE OF RESULTANT MATRICES AND
C         SCALE EIGENVECTORS
  255 II=1
      DO 275 I=1,N
      BB=DSQRT(B(II))
      DO 270 K=1,N
  270 X(K,I)=X(K,I)/BB
  275 II=II+N1-I
      RETURN
C
C	UPDATE D MATRIX AND START NEW SWEEP,IF ALLOWED
C
  280 DO 290 I=1,N
  290 D(I)=EIGV(I)
      IF(NSWEEP.LT.NSMAX) GO TO 40
      GO TO 255
 2000 FORMAT(/1X,'SWEEP NUMBER IN  *JACOBI* =',I4)
 2010 FORMAT(1X,6E20.12)
 2020 FORMAT(/1X,'*** ERROR SOLUTION STOP'/1X,'MATRICES NOT  
     &  POSITIVE DEFINITE'/1X,'II=',I4,5X,'A(II)=',
     &  E20.12,5X,'B(II)=',E20.12)
 2030 FORMAT(/1X,'CURRENT EIGENVALUES IN *JACOBI* ARE',/)
      END
C
C	=======================================================================
C ****************************************************************
C *                                                              *
C *   P R O G R A M       TO 
C *           EVALUATE SHIFT FOR STURM SEQUENCE CHECK            *
C *                                                              *
C ****************************************************************
      SUBROUTINE SCHECK1 (EIGV,RTOLV,BUP,BLO,BUPC,NEIV,NC,NEI,
     &  RTOL,SHIFT,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      
CMKIM 2000. 3. 7
C     DIMENSION EIGV(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC),NEIV(NC)
      DIMENSION EIGV(*),RTOLV(*),BUP(*),BLO(*),BUPC(*),NEIV(*)
C
      FTOL=0.01
C
      DO 100 I=1,NC
      BUP(I)=EIGV(I)*(1.+FTOL)
  100 BLO(I)=EIGV(I)*(1.-FTOL)
      NROOT=0
      DO 120 I=1,NC
  120 IF(RTOLV(I).LT.RTOL) NROOT=NROOT+1
      IF(NROOT.GE.1) GOTO 200
      WRITE(IOUT,1010)
      STOP
C
C     FIND UPPER BOUNDS ON EIGENVALUE CLUSTERS
C
  200 DO 240 I=1,NROOT
  240 NEIV(I)=1
      IF(NROOT.NE.1) GOTO 260
      BUPC(1)=BUP(1)
      LM=1
      L=1
      I=2
      GOTO 295
  260 L=1
      I=2
  270 IF(BUP(I-1).LE.BLO(I)) GOTO 280
      NEIV(L)=NEIV(L)+1
      I=I+1
      IF(I.LE.NROOT) GOTO 270
  280 BUPC(L)=BUP(I-1)
      IF(I.GT.NROOT) GOTO 290
      L=L+1
      I=I+1
      IF(I.LE.NROOT) GOTO 270
      BUPC(L)=BUP(I-1)
  290 LM=L
      IF(NROOT.EQ.NC) GOTO 300
  295 IF(BUP(I-1).LE.BLO(I)) GOTO 300
      IF(RTOLV(I).GT.RTOL) GOTO 300
      BUPC(L)=BUP(I)
      NEIV(L)=NEIV(L)+1
      NROOT=NROOT+1
      IF(NROOT.EQ.NC) GOTO 300
      I=I+1
      GO TO 295
C
C     FIND SHIFT
C
  300 WRITE(IOUT,1020)
      WRITE(IOUT,1005)(BUPC(I),I=1,LM)
      WRITE(IOUT,1030)
      WRITE(IOUT,1006)(NEIV(I),I=1,LM)
      LL=LM-1
      IF(LM.EQ.1) GO TO 310
  330 DO 320 I=1,LL
  320 NEIV(L)=NEIV(L)+NEIV(I)
      L=L-1
      LL=LL-1
      IF(L.NE.1) GO TO 330 
  310 WRITE(IOUT,1040)
      WRITE(IOUT,1006)(NEIV(I),I=1,1,LM)
      L=0
      DO 340 I=1,LM
      L=L+1
      IF(NEIV(I).GE.NROOT) GO TO 350
  340 CONTINUE
  350 SHIFT=BUPC(L)
      NEI=NEIV(L)
C
      RETURN
C
 1005 FORMAT(1X,6E22.14)
 1006 FORMAT(1X,6I22)
 1010 FORMAT(/1X,'****ERROR SOLUTION STOP IN *SCHECK*'/12X,
     & 'NO EIGENVALUES FOUND',/1X)
 1020 FORMAT(/,1X,'UPPER BOUNDS EIGENVALUE CLUSTERS')
 1030 FORMAT(/1X,'NO OF EIGENVALUES IN EACH CLUSTER')
 1040 FORMAT(/1X,'NO OF EIGENVALUES LESS THAN UPPER BOUNDS')
      END
C
C	=======================================================================

      SUBROUTINE EFFSTIFB(AM,AC,AK,ND,DELT,ALPHA,BETA,NEQ,NWK,DDK,
	1				   TYP1,TYP2,TYP3,EFTF,TEMP)

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP1,TYP2,TYP3,EFTF,TEMP
C	-----------------------------------------------------------------------
C	INPUT:
C	------
C	AM(NWK) : MASS MATRIX
C	AC(NWK) : DAMPING MATRIX
C	AK(NWK) : STIFFNESS MATRIX
C	ND(NEQ+1) : VECTOR CONTAINING ADRESSES OF DIAGONAL ELEMENTS OF STIFFNESS MATRIX
C	DELT : INCREMENTAL TIME
C	ALPHA : NEWMARK COEFFICIENT
C	BETA
C	NEQ : NUMBER OF EQUATIONS
C	NWK : NUMBER OF ELEMENTS OF STIFFNESS MATRIX

C	OUTPUT
C	AK(NWK) : [L]*[D]*[L]**T OF EFFECTIVE STIFFNESS MATRIX

C	AVERAGE ACCELERATION METHOD : UNCONDITIONALLY STABLE(RECOMMENDED)
C	ALPHA=1/4, BETA=1/2

C	LINEAR ACCELERATION METHOD : CONDITIONALLY STABLE
C	ALPHA=1/6, BETA=1/2
C	-----------------------------------------------------------------------
      COMMON /IOUNIT/IN,IOUT,NSTIF
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION AM(1),AC(1),AK(1),ND(1),DDK(1)
c      DIMENSION AM(NWK),AC(NWK),AK(NWK),ND(NEQ+1) !7Nov2005
C	-----------------------------------------------------------------------
      IF(ALPHA.GT.0.25.AND.BETA.NE.0.5) THEN
        WRITE(*,1000) ALPHA,BETA
        STOP
      ENDIF

      A1=1./(ALPHA*DELT**2)
      A3=BETA/(ALPHA*DELT)

C	CALCULATE EFFECTIVE STIFFNESS MATRIX
C	DO 10 ITE=1,NWK
C	    AK(ITE)=A1*AM(ITE)+A3*AC(ITE)+AK(ITE)
C  10 CONTINUE
	CALL MDOPER(EFTF,TYP1,TYP2,1.0D0,A1,'ADD')
	CALL MDOPER(EFTF,EFTF,TYP3,1.0D0,A3,'ADD')


C	DECOMPOSE EFFECTIVE STIFFNESS MATRIX INTO [L]*[D]*[L]**T
C	ISH              : EIGENVALUE SHIFT INDEX(ISH=0:NOT USED,ISH=1:USED)

C	ISH=0
C	CALL DECOMP(AK,ND,NEQ,ISH,IOUT)
C	INDPD = 0
C	CALL COLSOL(ND,AK,DDK,DDK,1,INDPD,EFTF,TEMP)



      RETURN
 1000 FORMAT(//,1X,'INPUT ERROR IN NEWMARK COEFFICIENT',/
     &1X,'RECOMMENDED VALUES ALPHA = 0.15, BETA=0.5 BUT INPUT VALUES',/
     &1X,'ARE AS FOLLOWING',/
     &1X,'  ALPHA = ',F10.3/
     &1X,'  BETA  = ',F10.3)
      END
C
C	=======================================================================
      SUBROUTINE NMSOLB(AM,AC,AKE,R,RN,FN,ND,DELT,ALPHA,BETA,NEQ,NWM,
	1				  NWK,DDK,EFTF,MASS,DAMP)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 EFTF,MASS,DAMP
C	----------------------------------------------------------------------
C	INPUT
C	------
C	AM(NWK)          : MASS MATRIX
C	AC(NWK)          : DAMPING MATRIX
C	AKE(NWK)         : [L]*[D]*[L]**T OF EFFECTIVE STIFFNESS MATRIX
C	R(NEQ,3)         : (I)TH RESPONCE[DIS,VEL,ACCE]
C	FN(NEQ)          : (I+1)TH EXTERNAL FORCE
C	ND(NEQ+1)        : VECTOR CONTAINING ADDRESSES OF DIAGONAL ELEMENTS
C	DELT             : TIME INCREMENT
C	ALPHA            : NEWMARK COEFFICIENT
C	BETA     
C	NEQ              : NUMBER OF EQUATIONS
C	NWK              : NUMBER OF ELEMENTS OF STIFFNESS MATRIX
C	NWM              : NUMBER OF ELEMENTS OF MASS MATRIX
C					   = NWK FOR CONSISTENT MASS MATRIX
C					   = NEQ FOR LUMPED MASS MATRIX

C	OUTPUT
C	-------
C	RN(NEQ,3)        : (I+1)TH RESPONCE[DIS,VEL,ACCE]
C	-----------------------------------------------------------------------
CN      COMMON /IOUNIT/IN,IOUT,NSTIF
C
C     DIMENSION AM(NWK),AC(NWK),AKE(NWK),R(NEQ,3),RN(NEQ,3),FN(NEQ) !CHANED TO NEXT 13OCT2005
      DIMENSION AM(NWM),AC(NWK),AKE(NWK),R(NEQ,3),RN(NEQ,3),FN(NEQ)
      DIMENSION ND(NEQ+1),DDK(1)
C	-----------------------------------------------------------------------
      A1=1./(ALPHA*DELT**2)
      A2=1./(ALPHA*DELT)
      A3=BETA*A2
      A4=1./(2.*ALPHA)-1.
      A5=BETA/ALPHA-1.
      A6=BETA/ALPHA-2.
      A7=1.-BETA
      A8=1.-BETA/(2.*ALPHA)
      A8=A8*DELT

      DO 5 JTE=1,3
        DO 5 ITE=1,NEQ
          RN(ITE,JTE)=0.0D0          
    5 CONTINUE
      DO 10 ITE=1,NEQ
        RN(ITE,2)=A1*R(ITE,1)+A2*R(ITE,2)+A4*R(ITE,3)
   10 CONTINUE
C      CALL MULTDM(RN(1,3),AM,RN(1,2),ND,NEQ,NWK)
C	Previous line changed to the next by Nguyen DV 13Oct2005 for various mass matrix

C	CALL MULTDM(RN(1,3),AM,RN(1,2),ND,NEQ,NWM)
	CALL MAMULT(ND,AM,RN(1,2),RN(1,3),MASS,'STD')

      DO 15 ITE=1,NEQ
        RN(ITE,1)=FN(ITE)+RN(ITE,3)
   15 CONTINUE

      DO 20 ITE=1,NEQ
        RN(ITE,2)=A3*R(ITE,1)+A5*R(ITE,2)+A6*DELT/2.*R(ITE,3)
   20 CONTINUE

C	CALL MULTDM(RN(1,3),AC,RN(1,2),ND,NEQ,NWK)
	CALL MAMULT(ND,AC,RN(1,2),RN(1,3),DAMP,'STD')

      DO 25 ITE=1,NEQ
        RN(ITE,1)=RN(ITE,1)+RN(ITE,3)
   25 CONTINUE
      
C	CALL REDBAK(AKE,RN(1,1),ND,NEQ)
	CALL COLSOL(ND,AKE,DDK,RN(1,1),1,INDPD,EFTF,'TEMP')
	CALL COLSOL(ND,AKE,DDK,RN(1,1),2,INDPD,EFTF,'TEMP')


      DO 30 ITE=1,NEQ
c        RN(ITE,2)=A3*(RN(ITE,1)-R(ITE,1))-A5*R(ITE,2)+A8*R(ITE,3)
c        RN(ITE,3)=A1*(RN(ITE,1)-R(ITE,1))-A2*R(ITE,2)-A4*R(ITE,3)
C	  PREVIOUS CHANGED TO THE NEXT (use previous, the structure vibrate nagatively only???)
CMKIM 2001. 7. 13.
      RN(ITE,3)=A1*(RN(ITE,1)-R(ITE,1))-A2*R(ITE,2)-A4*R(ITE,3) !Re-activated for DMF1D
      RN(ITE,2)=R(ITE,2)+A7*DELT*R(ITE,3)+BETA*DELT*RN(ITE,3)

CMKIM THERE IS A SUBROUTINE TO CALCULATE THE ACCELERATION FROM THE EQUILIBRIUM EQUATION
CMKIM OUTSIDE OF THIS SUBROUTINE.
CMKIM THE NAME OF THAT SUBROUTINE IS CALACC.
   30 CONTINUE

      RETURN
      END
C
C	=======================================================================
      SUBROUTINE CALACCB(F,RN,AM,AC,AK,TT,ND,NEQ,NWM,NWK,
	1				  DDM,STIF,MASS,DAMP)
      IMPLICIT REAL*8(A-H,O-Z)
	CHARACTER*4 STIF,MASS,DAMP
c	---------------------------------------------------------------------
C	PROGRAM TO RE-CALCULATE THE ACCELERATION BY THE EQUILIBRIUM EQUATION:
C	---------------------------------------------------------------------
C	INPUT:
C	-------
C	F(NEQ)	  = (I+1)TH EXTERNAL FORCE
C	RN(NEQ,3) = (I+1)TH RESPONSES(DIS,VEL,ACC)
C                 THIS RESULTS ARE CALCULATED IN THE NEWMARK BETA METHOD
C	AM(NWK)   = DECOMPOSED MASS MATRIX
C	AC(NWK)   = DAMPING MATRIX
C	AK(NWK)   = STIFFNESS MATRIX
C	ND(NEQ+1) = VECTOR CONTAINING ADDRESSES OF DIAGONAL ELEMENTS
C	NEQ		  = NUMBER OF EQUATIONS
C	NWK		  = NUMBER OF ELEMENTS
C	NWM       = NUMBER OF ELEMENTS OF MASS MATRIX
C				= NWK FOR CONSISTENT MASS MATRIX
C				= NEQ FOR LUMPED MASS MATRIX

C	OUTPUT:
C	-------
C	RN(1,3)	  = ACCELERATION CALCULATED FROM EQUILIBRIUM EQUATION

C	INDEX:
C	TT(NEQ)   = AUXILIARY VECTOR
C	------------------------------------------------------------------
c      COMMON /IOUNIT/IN,IOUT,NSTIF
C
      DIMENSION F(NEQ),RN(NEQ,3),AM(1),AC(1),AK(1),TT(NEQ),
     1		  ND(NEQ+1),DDM(1)

C	CALCULATE TT=AC*RN(1,2)

C	CALL MULTDM(TT,AC,RN(1,2),ND,NEQ,NWK)
	CALL MAMULT(ND,AC,RN(1,2),TT,DAMP,'STD')

      DO 10 I=1,NEQ
        F(I)=F(I)-TT(I)
   10 CONTINUE

C	CALCULATE TT=AK*RN(1,1)
C	CALL MULTDM(TT,AK,RN(1,1),ND,NEQ,NWK)
	CALL MAMULT(ND,AK,RN(1,1),TT,STIF,'STD')


      DO 20 I=1,NEQ
        F(I)=F(I)-TT(I)
   20 CONTINUE

C	RE-CALCULATE ACCELERATION BY THE EQUILIBRIUM EQUATION
C	changed 13Oct0025 by NguyenDV to for lumped mass matrix

C	CALL REDBAK(AM,F,ND,NEQ)
	INDPD = 0
	CALL COLSOL(ND,AM,DDM,F,1,INDPD,MASS,'TEMP')
	CALL COLSOL(ND,AM,DDM,F,2,INDPD,MASS,'TEMP')


C	RESOTRE ACCELERATION TO RN(1,3)
      DO 30 I=1,NEQ
        RN(I,3)=F(I)
   30 CONTINUE

      RETURN
      END
C
C	=======================================================================





C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE PRINSTIFC(A,MAXA,NEQ)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------

      DIMENSION A(1),MAXA(1)
      
      DIMENSION ROW(NEQ),COL(NEQ),LROW(NEQ)
      
      
      DO IROW = 1,NEQ
          ROW(1:NEQ) = 0.0D0
          KCOL = 1
          DO ICOL = IROW,NEQ
              NCOH = MAXA(ICOL+1) - MAXA(ICOL+0)
              DO ICOH = 1,NCOH
              JROW = ICOL - ICOH + 1
                  IF(JROW.EQ.IROW) THEN
                      INDX = MAXA(ICOL) + ICOH - 1
                      ROW(ICOL) = A(INDX)
                      KCOL = ICOL
                      GOTO 50
                  ENDIF
              ENDDO
50            CONTINUE    
          ENDDO
          WRITE(110,100) ROW(1:KCOL)
      ENDDO          
         
      

      DO IROW = 1,NEQ
          LROW(1:NEQ) = 0.0D0
          KCOL = 1
          DO ICOL = IROW,NEQ
              NCOH = MAXA(ICOL+1) - MAXA(ICOL+0)
              DO ICOH = 1,NCOH
              JROW = ICOL - ICOH + 1
                  IF(JROW.EQ.IROW) THEN
                      INDX = MAXA(ICOL) + ICOH - 1
                      LROW(ICOL) = INDX
                      KCOL = ICOL
                      GOTO 80
                  ENDIF
              ENDDO
80            CONTINUE    
          ENDDO
          WRITE(110,200) LROW(1:KCOL)
      ENDDO   
      
100   FORMAT(500E13.5)  
200   FORMAT(500I6)  
         
         
      DO IROW = 1,NEQ
      IDIR = MAXA(IROW)
      WRITE(*,*) IROW,A(IDIR)
      ENDDO
      
     
      RETURN
      END
C	=======================================================================
C	=======================================================================
C	=======================================================================






