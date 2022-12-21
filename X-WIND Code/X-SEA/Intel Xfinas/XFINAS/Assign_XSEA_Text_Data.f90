      SUBROUTINE OPENXSEADATA (IOPENCASE)
	USE IFPORT
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      !OGICAL(4) MAKEDIR
      CHARACTER*200   PATH

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1,KFATL !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      COMMON /SOILDYNA/ NSOIL
      COMMON A(9000000),IA(9000000)
      
	COMMON /WLAS/ MERW,MEIW,NLASI,NLASW
	! ------------ COMMON FOR OFFSHORE ------------------------------
	COMMON /WARNING/ WARNING,RAMDA(100),RK
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION
	! ---------------------------------------------------------------
	
	COMMON /MEMW/ W(7000000),IW(7000000)

!      CALL TESTWRIT
C	----------------------------------
C	FILE NUMBER LIST FOR WHOLE PROGRAM
C	----------------------------------
	
C	XFINAS.FOR		DATA   INPUT  FILE DATIN1.DAT		ITI		= 1							SEQUENTIAL
C	XFINAS.FOR		BACKUP FILE FOR NONLINEAR		  NPLOT		= 4							SEQUENTIAL (SEE ALSO RETAKE,RETAKEC,RTKDYN)	
C	XFINAS.FOR		SCREEN OUTPUT						ITO		= 6							SCREEN
C	XFINAS.FOR		DATA   OUTPUT FILE OUTPUT.OUT		ISO		= 7							SEQUENTIAL

C	INPUT.FOR		    SUPPRESS THE FLYING DOF                 = 8							SEQUENTIAL (DELETE AFTER CALLING NNSUPPRS OPTION 'SUPP' -- AFTER ELEMIN)

C	XFINAS.FOR		GID    OUTPUT FILE GIDOUT.FLAVIA.RES		= 101						SEQUENTIAL											

C	INPUT.FOR		FILE FOR ELEMENT GROUP INTERGER DATA        = 10  + NEG					SEQUENTIAL (SEE ALSO ELEMIN)
C	INPUT.FOR		FILE FOR ELEMENT GROUP REAL     DATA        = 30  + NEG					SEQUENTIAL (SEE ALSO ELEMIN)


C	XFINAS.FOR		OUTPUT FILE					OUTJOB.OUT		= 100						SEQUENTIAL	

C	INPUT.FOR		    OUTPUT GAUSSPOINT POSITION              = 102						SEQUENTIAL (SEE ALSO GAUSSPRIN)	

C	CSALDCOM.FOR	    TENDON OUTPUT                           = 103						SEQUENTIAL	

C	XFPrinout.FOR    OUTPUT FOR MAKING EXCEL REPORT	(DISP.)     = 106						SEQUENTIAL	
C	XFPrinout.FOR    OUTPUT FOR MAKING EXCEL REPORT	(ELEMENT)	= 107						SEQUENTIAL	
C	XFPrinout.FOR    OUTPUT FOR MAKING EXCEL REPORT	(SUPPORT)	= 108						SEQUENTIAL	
C	XFPrinout.FOR    OUTPUT FOR MAKING EXCEL REPORT	(LINK)	    = 109						SEQUENTIAL	


C	XFINAS.FOR		TEMPORARY FILE				OUTSUP.OUT		= 110						SEQUENTIAL	

C	GENFIBER.FOR		GID    OUTPUT FILE FRAMESEC.SEC         = 111						SEQUENTIAL			

C	PASSWORD.FOR	FOR SECURITY CHECKING                       = 202						SEQUENTIAL
C	PASSWORD.FOR	FOR SECURITY CHECKING                       = 203						SEQUENTIAL
C	PASSWORD.FOR	FOR SECURITY CHECKING                       = 303						SEQUENTIAL
C	PASSWORD.FOR	FOR SECURITY CHECKING                       = 505						SEQUENTIAL	
	
C	XFINAS.FOR		BRIDGE-TRAIN INTERACTION	OUTRAIN.OUT		= 500						SEQUENTIAL	

C	NEWSOLVER.FOR	GLOBAL STIFFNESS                            = 301						SEQUENTIAL			
C	NEWSOLVER.FOR	GLOBAL MASS                                 = 302						SEQUENTIAL		
C	NEWSOLVER.FOR	GLOBAL DAMPING                              = 303						SEQUENTIAL		
C	NEWSOLVER.FOR	GLOBAL GEOMETIC STIFFNESS                   = 304						SEQUENTIAL	
C	NEWSOLVER.FOR	BACK UP GLOBAL STIFFNESS                    = 305						SEQUENTIAL			
C	NEWSOLVER.FOR	BACK UP GLOBAL MASS                         = 306						SEQUENTIAL		
C	NEWSOLVER.FOR	BACK UP GLOBAL DAMPING                      = 307						SEQUENTIAL		
C	NEWSOLVER.FOR	BACK UP GLOBAL GEOMETIC STIFFNESS           = 308						SEQUENTIAL
C	NEWSOLVER.FOR	ELEMENT STIFFNESS FOR ASSEMBLY              = 309						SEQUENTIAL	
C	NEWSOLVER.FOR	BACK UP GLOBAL STIFFNESS FOR RETAKE         = 310						SEQUENTIAL (SEE ALSO RETAKE,RETAKEC,RTKDYN)	
C	NEWSOLVER.FOR	GLOBAL EFFECTIVE STIFFNESS FOR DYNAMIC      = 311						SEQUENTIAL			
C	NEWSOLVER.FOR	TEMPORARY FILE                              = 312						SEQUENTIAL			
C	NEWSOLVER.FOR	COLUMN INDEX FOR EACH VALUE IN STIFFNESS    = 313						SEQUENTIAL			
C	NEWSOLVER.FOR	IMCOMPLETE CHOLESKI FACTOR                  = 314						SEQUENTIAL		
C	NEWSOLVER.FOR	TEMPORARY FILE 1                            = 315						SEQUENTIAL (SEE ALSO MDMOVE )		
C	NEWSOLVER.FOR	TEMPORARY FILE 2                            = 316						SEQUENTIAL (SEE ALSO IMC_LDL) 	
C	FrameElemLoad.FOR	FILE FOR STORING OFFSHORE LOAD          = 317						DIRECT ACCESS (SEE ALSO SUBROUTINE INPUT & OFFSHFILE) VARY LAOD
C	FrameElemLoad.FOR	FILE FOR STORING OFFSHORE LOAD          = 318						DIRECT ACCESS (SEE ALSO SUBROUTINE INPUT & OFFSHFILE) CONSTANT LOAD

C	MOVINGLOAD.FOR	LANEDATA									= 600  + NLANE				SEQUENTIAL
C	MOVINGLOAD.FOR	LANE FIXEND FORCE DATA						= 650  + NLANE				SEQUENTIAL
C	MOVINGLOAD.FOR	RESPONSE DUE TO UNIT LOAD AT EACH STATION	= 700						DIRECT ACCESS
C	MOVINGLOAD.FOR	MAX RESPONSE DUE TO VEHICLE					= 701						DIRECT ACCESS
C	MOVINGLOAD.FOR	MAX RESPONSE DUE TO UNIT UNITFORM LOAD		= 702						DIRECT ACCESS
C	MOVINGLOAD.FOR	MAX RESPONSE DUE TO UNIT LOAD				= 703						DIRECT ACCESS
C	MOVINGLOAD.FOR	MAX RESPONSE DUE TO VEHICLE CLASS			= 704						DIRECT ACCESS
C	MOVINGLOAD.FOR	VAHICLE CLASS DATA							= 705						SEQUENTIAL
C	MOVINGLOAD.FOR	DUMMY FILE 1								= 706						DIRECT ACCESS
C	MOVINGLOAD.FOR	DUMMY FILE 2								= 707						DIRECT ACCESS
C	MOVINGLOAD.FOR	DUMMY FILE 3								= 708						DIRECT ACCESS
C	MOVINGLOAD.FOR	DUMMY FILE 4								= 709						DIRECT ACCESS

C	XFOUTPUT.FOR	LOAD CASE OUTPUT FILE                       = 800						DIRECT ACCESS
C	XFOUTPUT.FOR	LOAD COMBINATION INPUT                      = 801						SEQUENTIAL
C	XFOUTPUT.FOR	ELEMENT GROUP STORAGE FILE                  = 803  + 6*NEG				DIRECT ACCESS
C	XFOUTPUT.FOR	LOAD COMBINATION OUTPUT                     = 850						DIRECT ACCESS
C	SEISMA.FOR	RESPONSE SPECTRUM ANALYSIS DATA                 = 851						DIRECT ACCESS
C	SEISMA.FOR	DUMMY FOR MODE COMB FILE                        = 852						DIRECT ACCESS
C	INPUT.FOR		EIGENVALUE SOLVER                           = 853						DIRECT ACCESS  ... SEE ALSO IN SUBROUTINE EIGMOD
C	GENFIBER.FOR	XFINAS FRAME SECTION                        = 860  + NEG				DIRECT ACCESS  ... SEE ALSO SUB. FRAMSEC
C	GENFIBER.FOR	XFINAS FRAME SECTION                        = 880  + NEG				DIRECT ACCESS  ... STORE ALL INTEGRATED SECTION PROPERTY ALONG FRAME GAUSS POINT

C	CSASECTION.FOR    CONSTRUCTION STAGE MODULE                 = 3000			 			DIRECT ACCESS (FRAME SECTION)
C	CSAMESHIN.FOR	CONSTRUCTION STAGE MODULE                   = 3000 + 10*NEG				DIRECT ACCESS (ELEMENT GROUP + TENDON)
C	NEWSOLVER.FOR	BACK UP WORKING ARRAY FILE                  = 3200 + NEG                DIRECT ACCESS	
C	NEWSOLVER.FOR	NORMAL  WORKING ARRAY FILE                  = 3300 + NEG                DIRECT ACCESS
C	WARRAY.FOR	EAS ARRAY                                       = 3400 + NEG                DIRECT ACCESS

C	TENDONSEARCH.FOR  TEMPORARY DATA USE FOR SEARCHING          = 3500			 			SEQUENTIAL
C	TENDONSEARCH.FOR  TENDON ELEMENT (NODAL CONNECTIVITY)       = 3501			 			DIRECT ACCESS 
C	TENDONSEARCH.FOR  TENDON ELEMENT (JACKING FORCE)            = 3502			 			DIRECT ACCESS 

C	CSASPRING.FOR	CONSTRUCTION SPRING SUPPORT                 = 5000                      DIRECT ACCESS
C	INPUT.FOR 	DYNAMIC LOAD CURVE                              = 5001                      DIRECT ACCESS (SEE ALSO LOADCRV2)
      IF (IOPENCASE.EQ.1) THEN
C	---------------------------------------------------------
C     OPEN FILES
C	---------------------------------------------------------

      OPEN(UNIT=4,STATUS='SCRATCH',FORM='UNFORMATTED')

      OPEN(UNIT=1  ,FILE='datin1.dat'       ,STATUS='OLD'    )
      OPEN(UNIT=7  ,FILE='output.out'       ,STATUS='UNKNOWN')
      OPEN(UNIT=10,FILE='screenjob.dat'       ,STATUS='UNKNOWN')
      OPEN(UNIT=100,FILE='outjob.out'       ,STATUS='UNKNOWN')
      OPEN(UNIT=101,FILE='gidout.flavia.res',STATUS='UNKNOWN')
	OPEN(UNIT=102,FILE='gaussp.flavia.res',STATUS='UNKNOWN')
	OPEN(UNIT=103,FILE='tedout.flavia.res',STATUS='UNKNOWN')
	OPEN(UNIT=103,FILE='tedout.flavia.res',STATUS='UNKNOWN')
	OPEN(UNIT=110,FILE='outsup.out'       ,STATUS='UNKNOWN')
	OPEN(UNIT=111,FILE='framesec.sec'     ,STATUS='UNKNOWN')
	OPEN(UNIT=500,FILE='outrain.out'      ,STATUS='UNKNOWN')
	
C	---------------------------------------------------------
C     OPEN OUTPUT FILE FOR EXCEL REPORT
      !OPEN(UNIT=3  ,FILE='Output_Design_Code/Out Designcode.csv'       ,STATUS='UNKNOWN'    )
      OPEN(UNIT=106,FILE='Out Displacement.csv'                        ,STATUS='UNKNOWN')
      OPEN(UNIT=107,FILE='Out Element.csv'                             ,STATUS='UNKNOWN')
      OPEN(UNIT=108,FILE='Out Support.csv'                             ,STATUS='UNKNOWN')
      OPEN(UNIT=109,FILE='Out Link.csv'                                ,STATUS='UNKNOWN')
      OPEN(UNIT=120,FILE='Out Stress.csv'                              ,STATUS='UNKNOWN')    
C	---------------------------------------------------------

C     OPEN INOUT FILE FOR FAST      
      OPEN(UNIT=551,FILE='NCOUPLE.dat'         ,STATUS='UNKNOWN',SHARE='DENYNONE')
      OPEN(UNIT=201,FILE='dynamic_mass.dat'      ,STATUS='UNKNOWN')
      GOTO 100
      ENDIF
      
      IF (KOFFL.NE.0.0D0) IOPENCASE_OCEA = 1
      IF (NSOIL.NE.0.0D0) IOPENCASE_OCEA_SOIL = 1
      IF (NSUPER.EQ.1)    IOPENCASE_SUPER = 1
      IF (KFTTD.EQ.1.OR.KFTTD.EQ.2) IOPENCASE_FATIGUE_TIME = 1
      
      IF (IOPENCASE_FREQ.EQ.1) THEN  ! FREQUENCY DOMAIN
      MAKEDIR = MAKEDIRQQ('Result of Internal force')
      MAKEDIR = MAKEDIRQQ('FREQ DOMAIN')
      OPEN(UNIT=251,FILE='FREQ DOMAIN/Transfer Function Step.dat'               ,STATUS='UNKNOWN')
      OPEN(UNIT=252,FILE='FREQ DOMAIN/Displacment Step.dat'                     ,STATUS='UNKNOWN')
      OPEN(UNIT=254,FILE='FREQ DOMAIN/Temporary_file_displacement.dat'          ,STATUS='UNKNOWN')
      OPEN(UNIT=255,FILE='FREQ DOMAIN/Temporary_file_force.dat'                 ,STATUS='UNKNOWN')
      OPEN(UNIT=256,FILE='FREQ DOMAIN/Testing_freq.dat'                         ,STATUS='UNKNOWN')
      OPEN(UNIT=5200,FILE='Result of Internal force/Internal Force (Freq.).dat'              ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5201,FILE='Result of Internal force/Internal Force (Time.).dat'              ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5202,FILE='Result of Internal force/Pre-Internal Force (Time.).dat'          ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5203,FILE='Result of Internal force/Reaction Force  (Freq.).dat'             ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5204,FILE='Result of Internal force/Reaction Moment (Freq.).dat'             ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5205,FILE='Result of Internal force/Reaction Force  (Time.).dat'             ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5206,FILE='Result of Internal force/Reaction Moment (Time.).dat'             ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5207,FILE='Result of Internal force/Pre-Internal Reaction Force (Time.).dat' ,STATUS='UNKNOWN'  ,SHARE='DENYNONE')
      OPEN(UNIT=5208,FILE='Result of Internal force/Pre-Internal Reaction Moment (Time.).dat',STATUS='UNKNOWN' ,SHARE='DENYNONE')
      ENDIF

      IF (IOPENCASE_OCEA.EQ.1) THEN ! OCEAN ANALYSIS
      MAKEDIR = MAKEDIRQQ('Result of ocean analysis') 
      OPEN(UNIT=67,  FILE='Result of ocean analysis/Water Surface Elevation.out'                       ,STATUS='UNKNOWN')
      ENDIF

      IF (IOPENCASE_OCEA_SOIL.EQ.1)THEN
      MAKEDIR = MAKEDIRQQ('PILE-SOIL INTERACTION OUTPUT')
      OPEN(UNIT=171,FILE='PILE-SOIL INTERACTION OUTPUT/Pullout and End Bearing Capacity.dat'       ,STATUS='UNKNOWN')
      OPEN(UNIT=7120,FILE='PILE-SOIL INTERACTION OUTPUT/Summary Pile-Soil Interaction.dat'         ,STATUS='UNKNOWN')
      OPEN(UNIT=13004,FILE='PILE-SOIL INTERACTION OUTPUT/Superelement stiffness.dat'         ,STATUS='UNKNOWN') 
      ENDIF
      
      IF (IOPENCASE_FATIGUE_TIME.EQ.1)THEN
      MAKEDIR = MAKEDIRQQ('Result of Fatigue') 
      OPEN(UNIT=5100,FILE='Result of Fatigue/Internal Force of X-SEA.out'              ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5101,FILE='Result of Fatigue/Internal Force of X-SEA (Input data).txt' ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      OPEN(UNIT=5102,FILE='Result of Fatigue/Fatigue_Model_Input.dat'                  ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      ENDIF
      
      IF(IOPENCASE_SUPER.EQ.1)THEN
      MAKEDIR = MAKEDIRQQ('Superelement') 
      OPEN(UNIT=13005,FILE='Superelement/Superelement stiffness.dat'              ,STATUS='UNKNOWN' ,SHARE='DENYNONE')
      ENDIF
      
      

100   END SUBROUTINE
C     =======================================================================================================================================================
C     ==============================================
C     ==============================================   
C     ==============================================     
      SUBROUTINE TESTWRIT 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      CHARACTER*8 CVAL 
      
      DIMENSION A(100),B(100)
 
      RETURN
      
      LUN = 8111
      OPEN (LUN,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='STREAM')
      REWIND(LUN)
      
      LUN = 8112
      OPEN (LUN,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='STREAM')
      REWIND(LUN)
      
      RETURN
      
      
      LUN = 488
      OPEN (LUN,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='STREAM')
      
       DO I = 1,100
          A(I) = FLOAT(I)
      ENDDO
      
      WRITE(LUN,POS=1) A(1:10)
      N = 2*8+1
      REWIND(LUN)
      WRITE(LUN,POS=N) A(1)   
      READ(LUN,POS=1) B(1:10)
      
      WRITE(*,*) B(1:10)
      STOP
      
      
      
      
      
      
      LUN = 488
      OPEN (LUN,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='STREAM')
      
      DO I = 1,100
          A(I) = FLOAT(I)
      ENDDO
      
      REWIND(LUN)
      WRITE(LUN) A(1:10)
      N = 5*8+1
      WRITE(LUN,POS=N) A(21:30)
      
      N = 5*8+1
      READ(LUN,POS=N) B(1:10)
      
      WRITE(*,*) B(1:3)
      
      STOP
      RETURN
      
      
      DO I = 1,100
          A(I) = FLOAT(I)
      ENDDO
      
      WRITE(488) A(1:5)
      WRITE(488) A(1:4)
      
      REWIND(488)
      READ(488) B(1:9)
      
      WRITE(*,*) 'KAK1'
      WRITE(*,*) B(1:9)
      
      
      REWIND(488)
      A(1:4) = A(1:4)*10.0
      WRITE(488) A(1:4)

      REWIND(488)
      READ(488) B(1:9)
      
      WRITE(*,*) 'KAK1'
      WRITE(*,*) B(1:9)
      
      REWIND(488)
      READ(488) B(1:3)
      READ(488) B(1:3)
      
      WRITE(*,*) 'KAK2'
      WRITE(*,*) B(1:3)
      
      STOP
      END
C     ==============================================
C     ==============================================   
C     ==============================================  
      SUBROUTINE PRNVTEMP(IGM,EDATA,NOUT,LENGTH,MAPP,NV,METHOD,
     1                      LFPRIN,NPR1,NPR2)  !ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD
	
C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRN1,LFPRN2,LFPRN3,LFPRN4	
	
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(20),MAPP(1)
C     ----------------------------------------------------------------


	NUM = NOUT/LENGTH

      NLINE = (NPR2-NPR1)+1
      
	DO IPT = NPR1,NPR2

	DO J = 1,NV
	PDATA(J) = 0.0D0
	M = MAPP(J)

	IF(M.NE.0) THEN
	SELECTCASE(METHOD)
	CASE('MAX')
	PDATA(J) = EDATA(M+LENGTH*(IPT-1))
	CASE('MIN')
	PDATA(J) = EDATA(M+LENGTH*(IPT-1)+NOUT)
	ENDSELECT
	IF(ABS(PDATA(J)).LT.1.0E-12) PDATA(J) = 0.0D0
	ENDIF

	ENDDO


      IF(IPT.EQ.1) THEN
	  IF(KEXCEL.EQ.0) WRITE(LFPRIN) FLOAT(IGM),PDATA(1:NV)
        IF(KEXCEL.EQ.1) EXIT    
      ENDIF 
             
      IF(IPT.NE.1) THEN
        IF(KEXCEL.EQ.0) WRITE(LFPRIN) PDATA(1:NV)
        IF(KEXCEL.EQ.1) EXIT
      ENDIF
        
      ENDDO



      RETURN
      END
C

C	=====================================================================
C	=====================================================================	     
      SUBROUTINE PRNVALL(LFPRIN,LFTEMP,NUMV,NV)  !ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD
	
      ALLOCATABLE DATOUT(:)
      
      
      ALLOCATE(DATOUT(NUMV))

      REWIND(LFTEMP)
      READ(LFTEMP) DATOUT(1:NUMV)
      
      
      IF(NV.EQ.3) THEN
      WRITE(LFPRIN,100) DATOUT(1:NUMV)
      ENDIF
      
      IF(NV.EQ.4) THEN 
      WRITE(LFPRIN,110) DATOUT(1:NUMV)
      ENDIF
      
      IF(NV.EQ.6) THEN 
      WRITE(LFPRIN,120) DATOUT(1:NUMV)
      ENDIF
      
100	FORMAT(2X,F12.2,3E15.6/,2X,12X,3E15.6/,2X,12X,3E15.6/,2X,12X,3E15.6)
110	FORMAT(2X,F12.2,4E15.6/,2X,12X,4E15.6/,2X,12X,4E15.6/,2X,12X,4E15.6)
120	FORMAT(2X,F12.2,6E15.6/,2X,12X,6E15.6/,2X,12X,6E15.6/,2X,12X,6E15.6)
   
      DEALLOCATE(DATOUT)
      
      RETURN
      END
C

C	=====================================================================
C	=====================================================================	 

	SUBROUTINE PRNTYP9_O(ISTYP,ISTEP,METHOD)  !SHELL ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(8)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
     
      
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(11)
	HSTP = 'XShell3Q'
	CASE(1)
	HSTP = ' XShell4'
	CASE(4)
	HSTP = 'XShell4Q'
	CASE(12)
	HSTP = ' XShell9'
	CASE(3)
	HSTP = ' XShell8'
      CASE(13) ! FOR SHEAR WALL BY BJ
	HSTP = '  XWALL4'	
	CASE DEFAULT
	RETURN
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


      
      NUMV1 = 0
      NUMV2 = 0
      NUMV3 = 0
      NUMV4 = 0
      NUMV5 = 0
      NUMV6 = 0
      NUMV7 = 0
      NUMV8 = 0
      REWIND(8111)
      REWIND(8112)
      REWIND(8113)
      REWIND(8114)
      REWIND(8115)
      REWIND(8116)
      REWIND(8117)
      REWIND(8118)
      
C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
      
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

      NV = 3
	MAPP(1:3) = [1,2,0]
      LFTEMP = 8111
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV1 = NUMV1 + (1+NV*NPT)

      NV = 3
	MAPP(1:3) = [3,4,5]
      LFTEMP = 8112
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV2 = NUMV2 + (1+NV*NPT)
      
      NV = 3
	MAPP(1:3) = [6,7,8]
      LFTEMP = 8113
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV3 = NUMV3 + (1+NV*NPT)
      
      GOTO 1001
      
	NV = 6
	MAPP(1:6) = [9,10,0,11,0,0]
      LFTEMP = 8114
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV4 = NUMV4 + (1+NV*NPT)
      
	NV = 6
	MAPP(1:6) = [12,13,0,14,0,0]
      LFTEMP = 8115
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV5 = NUMV5 + (1+NV*NPT)
      
      NV = 6
	MAPP(1:6) = [9,10,11,12,13,14]
      LFTEMP = 8116
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV6 = NUMV6 + (1+NV*NPT)
      
      NV = 6
	MAPP(1:6) = [15,16,17,18,19,20]
      LFTEMP = 8117
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV7 = NUMV7 + (1+NV*NPT)
      
      NV = 6
	MAPP(1:6) = [21,22,23,24,25,26]
      LFTEMP = 8118
	CALL PRNVTEMP(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFTEMP,1,NPT)
      NUMV8 = NUMV8 + (1+NV*NPT)
      
1001	CONTINUE

      ENDDO
      
      
      
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3010) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL 
      NV = 3
      LFTEMP = 8111
	CALL PRNVALL(LFPRIN,LFTEMP,NUMV1,NV)
	CALL PRNEND(METHOD)

      
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3011) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
      NV = 3
      LFTEMP = 8112
	CALL PRNVALL(LFPRIN,LFTEMP,NUMV2,NV)
	CALL PRNEND(METHOD)


	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3002) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3012) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
      NV = 3
      LFTEMP = 8113
	CALL PRNVALL(LFPRIN,LFTEMP,NUMV3,NV)
	CALL PRNEND(METHOD)



      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A8,'-Membrane-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"'/,
	4		'Values')
3010	FORMAT(/'Result "',A8,'-Membrane-Stress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T"'/,
	4		   'Values')	

3001	FORMAT(/'Result "',A8,'-Bending-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-RR" "S-SS" "S-RS"'/,
	4		'Values')
3011	FORMAT(/'Result "',A8,'-Bending-Stress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-RR", "S-SS", "S-RS"'/,
	4		'Values')
	
3002	FORMAT(/'Result "',A8,'-Shear-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST" "S-RS"'/,
	4		'Values')
3012	FORMAT(/'Result "',A8,'-Shear-Stress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-RT", "S-ST", "S-RS"'/,
	4		'Values')

3003	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
3013	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S" ,"S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

3004	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
3014	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')


4001	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4011	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

4002	FORMAT(/'Result "',A8,'-Fiber-Stress-Mid ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4012	FORMAT(/'Result "',A8,'-Fiber-Stress-Mid ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

4003	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4013	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================