      SUBROUTINE OPENS
     
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /IOUNIT/IN,IOUT,NSTIF
      
      IN=10
      IOUT=20
      NSTIF=30

C      OPEN(UNIT=10,FILE='trtb.dat',STATUS='UNKNOWN',FORM='FORMATTED')
C	Next line blocked dueto Unit 20 is matching with Gid output 
CN      OPEN(UNIT=20,FILE='trtb.out',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=30,FILE='stif.txt',STATUS='UNKNOWN',FORM='UNFORMATTED')

      OPEN(UNIT=17,FILE='train.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=27,FILE='strut.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=47,FILE='trnmax.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=57,FILE='strmax.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=37,FILE='axlod.dat',STATUS='UNKNOWN',FORM='FORMATTED')


      OPEN(UNIT=90,FILE='maxa.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=91,FILE='stif.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=92,FILE='mass.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=93,FILE='damp.dat',STATUS='UNKNOWN',FORM='FORMATTED')
     
      RETURN
      END
C
C	=======================================================================
      SUBROUTINE EXITS
      
      COMMON /IOUNIT/IN,IOUT,NSTIF
      
      CLOSE(IN)
      CLOSE(IOUT)
      CLOSE(NSTIF)
      CLOSE(17)
      CLOSE(27)
      CLOSE(37)
      CLOSE(90)
      CLOSE(91)
      CLOSE(92)
      CLOSE(93)
      
      RETURN
      END
C
C	=======================================================================
      SUBROUTINE DMCONT
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	June2004 by NguyenDV
C     -------------------------------------------------
C     READS CONTROL VARIABLES FOR MOVING LOAD ANALYSIS
C	-------------------------------------------------
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
C
C	NNBR	= TOTAL NUMBER OF NODES OF BEAM FOR BRIDGE MODELLING
C	NNRA	= TOTAL NUMBER OF NODE OF BEAM FOR RAIL MODELLING
C	NEBEBR  = NUMBER OF BEAM ELEMENTS FOR BRIDGE MODELING
C	NPBEBR  = NUMBER OF MATERIAL PROPERTY CARDS OF BEAM ELEMENT FOR BRIDGE MODELLING
C	NLRB	= CARD NUMBER OF ELASTIC BEARINGS
C	NEWKTK  = NUMBER OF WINKLER ELEMENTS FOR TRACK MODELLING
C	NPWKTK  = NUMBER OF MATERIAL PROPERTY CARDS TO DEFINE WINKLER ELEMENT FOR TRACK MODELLING
C	NEBERA  = NUMBER OF BEAM ELEMENTS FOR RAIL MODELLING
C	NPBERA  = NUMBER OF MATERIAL PROPERTY CARDS OF BEAM ELEMENT MODELLING RAIL

C	NDOF	= NUMBER OF DEGREES OF FREEDOM OF NODE (=3)
C	NNBE	= NUMBER OF NODES MODELLING A BEAM (=2)
C	NNWK	= NUMBER OF NODES MODELLING A WINKLER ELEMENT(=4)
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /DMCO/ IDMA,IANA,IDSO,NNBR,NNRA,NEBEBR,NPBEBR,NLRB,
     +			  NEWKTK,NPWKTK,NEBERA,NPBERA,NDOF,NNBE,NNWK
C
C	Next added 30-10-2005
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW

C	----------------------------------------------------------------
	IF(ISOLOP.NE.11) GOTO 500
C	---------------------------------
C     READ DYNAMIC MOVING ANALYSIS TYPE:
C     ---------------------------------
	READ (ITI,*)

CN	READ (ITI,*,END=100) IDMA,IANA
	READ (ITI,*) IDMA,IANA,IDSO
   50 WRITE (ITO,1000)
      WRITE (10,1000)

C     -------------------------
C     READ SYSTEM CONTROL DATA:
C     -------------------------
	SELECT CASE (IDMA)

	  CASE(1)

        WRITE(ISO,1500)IDMA,IANA,IDSO !added IDSO 17Jun06

	  CASE(2)
	  NNBE = 2
		NNWK = 4
		NDOF = 3

		READ (ITI,*,END=200)NNBR,NNRA,NEBEBR,NPBEBR,NLRB,NEWKTK,
     +					    NPWKTK,NEBERA,NPBERA
     
  200		WRITE (ISO,2000) IDMA,IANA,IDSO,NNBR,NNRA,NEBEBR,NPBEBR,	!added IDSO 17Jun06
     +					 NLRB,NEWKTK,NPWKTK,NEBERA,NPBERA

	  CASE(3)
	  NNBE = 2
		NNWK = 4
		NDOF = 7

		NNRA   = 0
		NPBEBR = 0
		NLRB   = 0
		NEWKTK = 0
		NPWKTK = 0
		NEBERA = 0
		NPBERA = 0	

  300		WRITE (ISO,3000) IDMA,IANA,IDSO,NNRA,NLRB,		!added IDSO 17Jun06
     +					 NEWKTK,NPWKTK,NEBERA,NPBERA     			  	
	END SELECT
C     -------------------------
C     READ SYSTEM CONTROL DATA:
C     -------------------------
 1000 FORMAT (1X,'READ CONTROL VARIABLES FOR MOVING LOAD ANALYSIS ')

 1500 FORMAT (//18X,58(1H*)/18X,1H*,56X,1H*/
     +18X,'* CONTROLS FOR DYNAMICS OF INPUTED-AXLE LOAD OF VEHICLES *'/
     +18X,1H*,56X,1H*/18X,58(1H*)//
     +14X,'DYNAMIC MOVING MODELLING METHOD(1D,2D OR 3D). . .IDMA =',I9/
     +14X,'ANALYSIS METHOD OF MOVING LOAD . . . . . . . . . IANA =',I9/
     +14X,'DYNAMIC SOLUTION (Integ or Modal . . . . . . . . IDSO =',I9)

 2000 FORMAT (//18X,59(1H*)/18X,1H*,57X,1H*/
     +18X,'* CONTROLS FOR 2D DYNAMICS OF BRIDGE-TRAIN INTERACTION *'/
     +18X,1H*,57X,1H*/18X,59(1H*)//
     +14X,'DYNAMIC MOVING MODELLING METHOD(1D,2D OR 3D). . .IDMA =',I9/
     +14X,'ANALYSIS METHOD OF MOVING LOAD . . . . . . . . . IANA =',I9/
     +14X,'DYNAMIC SOLUTION (Integ or Modal . . . . . . . . IDSO =',I9/
     +14X,'TOTAL NUM. OF BEAM NODES FOR BRIDGE MODELLING . . .NBR =',I9/
     +14X,'TOTAL NUMBER OF BEAM NODES FOR RAIL MODELLING. . .NNRA =',I9/
     +14X,'NUM. OF BEAM ELEMENTS FOR BRIDGE MODELING. . . .NEBEBR =',I9/
     +14X,'NUM. OF MAT. PROP. CARDS FOR BRIDGE MODELLING. .NPBEBR =',I9/
     +14X,'CARD NUMBER OF ELASTIC BEARINGS . . . . . .  . .  NLRB =',I9/
     +14X,'NUM. OF WINKLER ELEMENTS FOR TRACK MODELLING. . NEWKTK =',I9/
     +14X,'NUM. OF MAT.PROP. ARDS TO DEFINE WINKLER ELEMENT FOR TRACK 
     +MODELLING .  NPWKTK =',I9/
     +14X,'NUM. OF BEAM ELEMENTS FOR RAIL MODELLING . . .  NEBERA =',I9/
     +14X,'NUM. OF MAT.PROP.CARDS OF BEAM ELEMENT MODELLING RAIL
     +. . . . . . .NPBERA = ',I9)

 3000 FORMAT (//18X,59(1H*)/18X,1H*,57X,1H*/
     +18X,'* CONTROLS FOR 3D DYNAMICS OF BRIDGE-VEHICLES INTERACTION *'/
     +18X,1H*,57X,1H*/18X,59(1H*)//
     +14X,'DYNAMIC MOVING MODELLING METHOD(1D,2D OR 3D). . .IDMA =',I9/
     +14X,'ANALYSIS METHOD OF MOVING LOAD . . . . . . . . . IANA =',I9/
     +14X,'DYNAMIC SOLUTION (Integ or Modal . . . . . . . . IDSO =',I9/
     +14X,'TOTAL NUMBER OF BEAM NODES FOR RAIL MODELLING. . .NNRA =',I9/
     +14X,'CARD NUMBER OF ELASTIC BEARINGS . . . . . .  . .  NLRB =',I9/
     +14X,'NUM. OF WINKLER ELEMENTS FOR TRACK MODELLING. . NEWKTK =',I9/
     +14X,'NUMBER OF MATERIAL PROPERTY SETS TO DEFINE '/ 
     +14X,'WINKLER ELEMENT FOR TRACK MODELLING . . . . . . NPWKTK =',I9/
     +14X,'NUM. OF BEAM ELEMENTS FOR RAIL MODELLING . . .  NEBERA =',I9/
     +14X,'NUMBER OF MATERIAL PROPERTY SETS OF BEAM ELEMENT '/ 
     +14X,'FOR RAIL MODELLING  . . . . . . . . . . . . . .NPBERA = ',I9)

  500 RETURN
      END
C
C=====================================================================
      SUBROUTINE LOCAW(IAW,IBW) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------
C     COMPUTES THE FIRST GLOBAL ADDRESSES OF ALL THE FLEXIBLE 
C	ARRAYS IN W (ADDITIONAL STORAGE ARRAY)
C	IAW = MAIN OPTION
C	IBW = SUB OPTION, if there is no sub-option, leave it =0
C	added IAW,IBW & changed 26/10/2005
C	---------------------------------------------
C     STORES THESE ADDRESSES IN COMMON BLOCK /LOCW/
C	---------------------------------------------
C     ARRAY IW      ADDRESS   CONTENTS
C	-------       -------   --------

C	 3D train-bridge
C	NEID(2,NELW) = I_NEID : NODE NUMBER OF WHEEL PATH ELEMENT
      
C	NCHTR(NEQTR) = I_NCHTR :VECTOR CONTANING COLUMN HEIGHT OF TRAIN'S STIFF/DAMPING MATRIX
C	NDTR(NEQTR+1)= I_NDTR  :VECTOR CONTAINING DIAGONAL ELEMENT ADDRESS OF TRAIN STIFFNESS/DAMPING MATRIX
C	IPRTR(NPRTR) = I_PRTR :DEGREES OF FREEDOM OF TRAIN TO PRINT
C	IPRST(NPRST,2)= I_PRST: NODES NUMBER OF STRUCTURE[NODE #, DOF #]

C     ARRAY W     ADDRESS  CONTENTS
C	-------     -------  --------



C     --------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C	Next added 22Oct2005
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C
	COMMON /SEIS/ ISEIDA,ISEIOP,IADIR,IAFM,ICOMB,IRES,NPICK,
	1			  NSTIM,NEAQ,NAV,NSAV

      COMMON /LOCW/ N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,NSTI,NAX,NAY,
     +			  NAZ,NSTIN,NAXI,NAYI,NAZI

      COMMON /DMCO/ IDMA,IANA,IDSO,NNBR,NNRA,NEBEBR,NPBEBR,NLRB,
     +			  NEWKTK,NPWKTK,NEBERA,NPBERA,NDOF,NNBE,NNWK

	COMMON /DMSOI/ IPATH,NMAX,NPRTTR,NPRTST,NPRWH,KPRTR
	COMMON /DMSOR/ ALPHA,BITA,DELT,VEL,EXTDIS,TOLER,SCALE,TLENGTH,BIDIS

C	Next common /LODW/ created Jul04 to store Dynamic Moving solution in W-array
	COMMON /LODW/ I_MPWK,I_LNWK,I_PRWK,I_MPBEBR,I_LNBEBR,I_PRBEBR

	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
	COMMON /BRIW/ I_PELW,I_NEID,I_PRST
	
	COMMON /TRN3/ NCARB,NACAB,NBOGI,NEQTR,NWKTR,NWMTR,IWRIN	
	
C	Next added 15Oct05
	COMMON /IWTR/ IT1,IT2,IT3,I_ITRAIN,I_NCHTR,I_NDTR,I_PRTR,I_PTWH
	COMMON /TRW3/ I_PCAB,I_CADI,I_ADDI,I_ALSU,I_PRSS,I_PBOG,I_PRPS,
	1			  I_PRWH,I_WRCO,I_WPOS,I_AMTR,I_ACTR,I_AKTR,I_AXLD
	COMMON /NIAX/ NWHEEL

	COMMON /IIR3/ IRIN,IRCLA,IRAN,NFFT,NOSIM !added 31Jan2007 by Nguyen
	COMMON /IRW3/ IRRx,IRRe,IRRa,IRRr	     !created 31Jan2007 by Nguyen


C	Next common /WLAS/ created Jul04 to store the last pointers in IW and W-arrays
	COMMON /WLAS/ MERW,MEIW,NLASI,NLASW

      COMMON A(9000000),IA(9000000)
C	COMMON /MEMW/ W(5000000)
C	Previous lined changed to the next, July04
	COMMON /MEMW/ W(7000000),IW(7000000)
C	----------------
C	W-ARRAY POINTERS
C	----------------
C	Next lines added Jul04 by NguyenDV
	GOTO(100,90,300,400,500),IAW
C	----------------------------------------------------
C	Recall the local pointers from "LANCZOS" subroutines
c	----------------------------------------------------
C	MR: Order of the reduced system (derived triangular matrix)
c	IF(NROOT.GT.3) MR = 4*NROOT
C	MR = 10
C	Label 300 added Jul04 by NguyenDV to separate range of pointers
C	W-ARRAY
  100	CONTINUE !NLASI = 1
	CONTINUE !NLASW = 1
	MR = 2*NROOT
      IF (NROOT.LE.20) MR = IFIX(20.0+(FLOAT(NROOT)-5.0)*50.0/45.0)+1
      IF (NROOT.LE.05) MR = 20

      IF (MR.GT.NEQ) MR = NEQ
C
	NN = MAX0(NEQ,2*MR)
C
C	N1  = 1  + NEQ
C	Previous line changed to the next Jul04 by NguyenDV
	N1  = NLASW  + NEQ
	N2  = N1 + NEQ	
	N3  = N2 + NEQ
	N4  = N3 + NEQ + 2
	N5  = N4 + MR
	N6  = N5 + MR
	N7  = N6 + MR
	N8  = N7 + NEQ
C	N8  = N7 + MR*(2*NEQ-MR+1)/2   SONGSAK SUPPRESS THIS TO REDUCE MEMORY USAGE OCT2019...WARNING LANCZOS WILL NOT WORK [LANC.FOR]...TRY TO USE SUBSPACE INSTEAD FOR EIGEN
	N9  = N8 + NROOT
	N10 = N9 + MR + ((NROOT - 1)*NN+1)/2
	N11 = N10+ NEQ
C	N11 = N10+ NROOT*NEQ           SONGSAK SUPPRESS THIS TO REDUCE MEMORY USAGE OCT2019...WARNING LANCZOS WILL NOT WORK [LANC.FOR]...TRU TO USE SUBSPACE INSTEAD FOR EIGEN

	N12 = N11+ MR         !CHANGE  NR  TO MR  SONGSAK JAN2007
c	NguyenDV,Oct11 added N12 to make N11 store eigenvalue pointer of original system

	IF (ISOLOP.EQ.9) GOTO 210
C	NLAST = N11
	NLASW = N12
	GOTO 90
C	--------------------------------------------
C	MORE GLOBAL POINTERS FOR W-ARRAY ARE CREATED
C	--------------------------------------------
C	Pointers for generated acce. components storage
C	------------------------
C   10 NSTI = N11 + NAV
  210 NSTI = N12 + NAV
	NAX  = NSTI+ NAV
	NAY  = NAX + NAV
	NAZ  = NAY + NAV

C	IF(NSTIM.GT.1) GOTO 20
C	Previous line changed to the next 7Apr04 by NguyenDV
C	to make NSTIN,NAXI,NAYI,NAZI store either Picked acce or Interpolated acce.
	IF(NSTIM.GT.1.OR.NPICK.GT.1) GOTO 220

	NLASW = NAZ 
	GOTO 90
C	WRITE(100,*)N1,NAV
C
C	Pointers for interpolated acce. components storage
C	-------------------------
  220	NSTIN = NAZ   + NSAV
	NAXI  = NSTIN + NSAV
	NAYI  = NAXI  + NSAV
	NAZI  = NAYI  + NSAV

	NLASW = NAZI
	GOTO 90
C	----------------------------------
C	LOCATI FOR 2D DYNAMIC MOVING LOADS:
C	----------------------------------
C	FOR WINKLER ELEMENT:
  300	CONTINUE !NLASI = 1
	CONTINUE !NLASW = 1

C	IW-ARRAY:
  	I_MPWK = NLASI
	I_LNWK = I_MPWK + NEWKTK
	NLASI  = I_LNWK + NEWKTK*NNWK

C	W-ARRAY:
	I_PRWK = NLASW
	NLASW  = I_PRWK + NPWKTK*2

C	TO STORE BRIDGE BEAM PROPERTIES
C	IW-ARRAY:
      I_MPBEBR = NLASI
      I_LNBEBR = I_MPBEBR+NEBEBR
      NLASI = I_LNBEBR+NEBEBR*NNBE 

C	W-ARRAY:
      I_PRBEBR = NLASW
C	Next IF to store total material set (Bridge & Rail)
	IF (IDATM.EQ.0) THEN
		NLASW=I_PRBEBR+NPBEBR*9
	ELSEIF (IDATM.EQ.1) THEN
		NLASW=I_PRBEBR+(NPBEBR+NPBERA)*9
	ENDIF

  	GOTO 90
C	----------------------------------
C	LOCATI FOR 3D DYNAMIC MOVING LOADS:
c	Added 6Sep2005
C	----------------------------------
C	STORE STRUCTURE PARAMETERS
  400	GOTO(410,420,430),IBW
C
  410	CONTINUE !NLASI = 1
	CONTINUE !NLASW = 1

C	IW-ARRAY:

C	W-ARRAY:

	GOTO 90

c	IW-ARRAY:
C	To store bridge's path element numbers (loadway elements)
  420	I_NEID = NLASI 
	NLASI  = I_NEID + 2*NELW      !2 Nodes per loadway segment

C	W-ARRAY:
      I_PELW = NLASW				!Local vector and length of segment
      NLASW  = I_PELW + 4*NELW	
	GOTO 90

C	STORE TRAIN PARAMETERS
C	IW-ARRAY:
  430	IF (IDMA.EQ.3.AND.IANA.EQ.1) THEN
	  I_ITRAIN=NLASI		!Train component poiter
	  NLASI=I_ITRAIN+3
C	  Next added 12Oct05
	  I_NCHTR= NLASI
	  I_NDTR = I_NCHTR + NEQTR	
	  NLASI  = I_NDTR  + NEQTR + 1
	ENDIF
C	Next added 15Oct05 to store DOFs of train and Bridge to print
  	IF (IDMA.EQ.3.AND.IANA.EQ.1) THEN
	  I_PRTR = NLASI
c	  I_PRST = I_PRTR + NPRTTR !changed next & added I_PRWH 03Mar07
c	  I_PRWH = I_PRTR + NPRTTR
c	  I_PRST = I_PRWH + NPRWH*2	 
c	  08May2008,found I_PRWH match with other pointer, changed to I_PTWH
	  I_PTWH = I_PRTR + NPRTTR
	  I_PRST = I_PTWH + NPRWH*2	   
	  	    
	  NLASI  = I_PRST + NPRTST*2
	ELSE
	  I_PRST = NLASI
	  NLASI  = I_PRST + NPRTST*2
	ENDIF

C	W-ARRAY:
	IF (IDMA.EQ.3.AND.IANA.EQ.1) THEN
C	Store the train properties
C	I_PCAB: PCAB(NCARB,4)  = Car body mechanical properties (MASS,Jx,Jz,Jy)
C	I_CADI: CADI(NCARB,7)  = Dimension parameters of cars (sA,sB,qA,qB,hC,h1A,h1B) 
C	I_ADDI: ADDI(NACAB,5)  = Additional Dimension Parameters for Articulated Cars (b3,h5A,h5B,h6A,h6B)
C	I_ALSU: ALSU(NACAB-1,3)= Longitudinal Suspensions of Articulated Cars  =(kTH,kTV,cX)
C	I_PRSS: PRSS(NBOGI,9)  = Secondary Suspensions (k2H,c2H,k2V,c2V,k2X,c2X,c2RX,b2,h2) 
C	I_PBOG: PBOG(NBOGI,4)  = Bogie of whole train (MASS,Jx,Jz,Jy)
C	I_PRPS: PRPS(NWHEEL,9)= Primary Suspensions (k1H,c1H,k1V,c1V,k1X,c1X,b1,h3,t)
C	I_PRWH: PRWH(NWHEEL,6) = Wheel-axle set properties (MASS,Jx,Jz,Jy,rW,b0)
C	I_WRCO: WRCO(NWHEEL,6) = Wheel-rail Contact Mechanism (Hertzian Spring) (kwX,cwX,kwY,cwY,kwY,cwY)
C	I_WPOS: WPOS(NWHEEL)   = Wheel position	
C    	
	  I_PCAB = NLASW
	  I_CADI = I_PCAB + NCARB*4
	  I_ADDI = I_CADI + NCARB*7

C	  I_ALSU = I_ADDI + NACAB*5
C	  I_PRSS = I_ALSU + (NACAB-1)*8  !Found wrong 13April2008
C	  Changed next 13April2008
	  IF(NACAB.GE.3) THEN
	  I_ALSU = I_ADDI + NACAB*5
	    I_PRSS = I_ALSU + (NACAB-1)*6
	  ELSEIF(NACAB.EQ.0) THEN
		I_ALSU = I_ADDI
	    I_PRSS = I_ALSU
	  ENDIF

	  I_PBOG = I_PRSS + NBOGI*9
	  I_PRPS = I_PBOG + NBOGI*4		

	  I_PRWH = I_PRPS + NWHEEL*9
c	  I_WRCO = I_PRWH + NWHEEL*5	!added b0 to the next 8Jan06
	  I_WRCO = I_PRWH + NWHEEL*6

	  I_WPOS = I_WRCO + NWHEEL*6	

C	  NLASW  = I_WPOS + NWHEEL !New, changed next two lines Dec2008
	  I_AXLD = I_WPOS + NWHEEL
	  NLASW  = I_AXLD + NWHEEL

C	  Store the train matrix in upper triangular column-wise
	  I_AMTR = NLASW
	  I_ACTR = I_AMTR + NEQTR
	  I_AKTR = I_ACTR + NWKTR
	  NLASW  = I_AKTR + NWKTR

C	  STORE 3-D IRREGULARITIES OF RAIL, 
c	  added 31Jan2007
	  IF(IRIN.EQ.3) THEN
		IRRx = NLASW
		IRRe = IRRx + NFFT
		IRRa = IRRe + NFFT
		IRRr = IRRa + NFFT
		NLASW = IRRr + NFFT
	  ENDIF

	ELSEIF (IDMA.EQ.1) THEN
	  I_AXLD = NLASW 
	  I_WPOS= I_AXLD + NWHEEL	!Found error WHPOS 8Mar06
	  NLASW  = I_WPOS+ NWHEEL 
	ENDIF

	GOTO 90
C	----------------------------------
C	INITIALIZE FIRST POSITION OF MEMORY
C	----------------------------------
500   NLASI = 1
	NLASW = 1
	
	GOTO 90
C	--------------------------------------------------
C	Memory requirement check and clear array W and IW:
C	Re-activated 13Aug06 by NguyenDV
C	--------------------------------------------------
	
   90	IF (NLASW.LT.MERW) GOTO 900 
	WRITE(6,2000)MERW,NLASW			
  
  900 IF (NLASI.LT.MEIW) GOTO 950	
	WRITE(6,2100)MEIW,NLASI		

C 950	CALL PZERO(W,NLASW)		

  950	RETURN
C  90 RETURN
 2000 FORMAT(/,5X,43H*******  WARNING : MEMORY SHORTAGE  *******,/,
     *       /,14X,22HMEMORY CAPACITY (IW) = ,I12,/,14X,
     *             22HMEMORY REQUIRED (IW) = ,I12)
 2100 FORMAT(/,5X,43H*******  WARNING : MEMORY SHORTAGE  *******,/,
     *       /,14X,22HMEMORY CAPACITY (W) = ,I12,/,14X,
     *             22HMEMORY REQUIRED (W) = ,I12)
C
	END
C
C=====================================================================
      SUBROUTINE GENTRAIN(NROW,NCOL,NTYP,PROP) 
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	PROGRAM TO READ & GENERATE PROPERTIES OF TRAIN (CAR BODY,SUSPENSION,BOGIE)
C		Programmed 30Oct04 by NguyenDV
C	-----------------------------------------------------------------------
C	INPUT:
C	------     
C	NROW	= NUMBER OF ROWS = NUMBER OF FULLY GENERATED INPUT PROPERTIES
C	NCOL	= NUMBER OF COLUMNS
C	NTYP	= NUMBER OF TYPICAL PROPERTIES
C	
C	OUTPUT:
C	-------
C	PROP(NROW,NCOL)	= INPUT PROPERTIES (FULLY GENERATED)	
C	----------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C
c	DIMENSION NTYPROP(NTYP),NTYP2(NTYP),PROP(NROW,*)
	DIMENSION NTYPROP(NTYP),NTYP2(NTYP),PROP(NROW,NCOL)
C	---------------------------------------------------------------------- 
C	Read Typical Properties Numbers 
	CALL CLEARI(NTYBOD,NTYP)

	READ(ITI,*) (NTYPROP(I),I=1,NTYP)
 
      DO 120 I=1,NROW
		DO 100 ITB = 1,NTYP
			IF (I.EQ.NTYPROP(ITB)) THEN
				READ(ITI,*) NTYP2(ITB),(PROP(I,J),J=1,NCOL)
			ENDIF

			IF(ITB.EQ.1) GOTO 100
C			!As the 1st propperty is always different from others

C			Generate all properties from typical properties 
			NT1 = NTYPROP(ITB-1)
			NT = NTYPROP(ITB)
			NNT = NT - NT1
			IF (NNT.LE.0) THEN
 				WRITE (ITO,*)'CHECK,INPUT TYPICAL BODY NUMBER ERROR'
                  WRITE (10,*)'CHECK,INPUT TYPICAL BODY NUMBER ERROR'
				GOTO 600
				STOP
			ELSEIF (NNT.EQ.1) THEN
			    GOTO 100	
			ELSEIF (NNT.GE.2) THEN
				DO 105 K=1,NNT - 1
					DO 105 J = 1,NCOL
  105					PROP(NT1+K,J) = PROP(NT1,J)

C			Check the Input of typical properties
			  AVPROCI = 0.
			  TOLI = 1.0E-12
			  DO 110 J = 1,NCOL
  110			  AVPROCI = 0.25*(AVPROCI + PROP(NT,J)- PROP(NT1,J))

			  IF (ABS(AVPROCI).EQ.TOLI) THEN  ! Need to check later ???
  115                 WRITE (ITO,*)'CHECK, INPUT PROPERTIES ERROR'
                      WRITE (10,*)'CHECK, INPUT PROPERTIES ERROR'
				GOTO 600
				STOP
			  ENDIF
			ENDIF
				
  100 CONTINUE
  120 CONTINUE
 
  600 RETURN
      END
C
C	=======================================================================