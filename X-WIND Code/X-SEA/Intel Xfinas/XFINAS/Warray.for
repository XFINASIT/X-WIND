C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOKN (IEG,ITYPE,ISTYP,MTMOD,NPT,NGT,NGPS,ISECT,
	1				   NWA,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)



	SELECTCASE(ITYPE)

	CASE(2)
	CALL INWOK2 (ISTYP,MTMOD,NWG)
	CASE(5)
	CALL INWOK5 (ISTYP,MTMOD,NWG,IEG,NGPS,ISECT)
	CASE(6)
	CALL INWOK6 (ISTYP,MTMOD,NWG)
	CASE(8)
	CALL INWOK8 (ISTYP,MTMOD,NWG)
	CASE(9)
	CALL INWOK9 (ISTYP,MTMOD,NWG,NGT)
	CASE(10)
	CALL INWOK10(ISTYP,MTMOD,NWG,NGT)
	CASE(11)
	CALL INWOK11(ISTYP,MTMOD,NWG)
	CASE(13)
	CALL INWOK13(ISTYP,MTMOD,NWG)
	CASE(14)
	CALL INWOK14(ISTYP,MTMOD,NWG)
	CASE(15)
	CALL INWOK15(ISTYP,MTMOD,NWG)
	CASE(16)
	CALL INWOK16(ISTYP,MTMOD,NWG)
	CASE(20)
	CALL INWOK20(ISTYP,MTMOD,NWG)
	CASE(17)
	CALL INWOK17(ISTYP,MTMOD,NWG)

	ENDSELECT


	NWA = NWG*NPT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK2 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	NWG = 1

	SELECTCASE(ISTYP)

	CASE(3) !TRUSS
	IF(MTMOD.EQ.3) NWG = 4

	CASE(4) !PARABOLIC CABLE
	NWG = 4

	CASE(5) !GAP
	NWG = 2

	CASE(6) !BONDLINK
	NWG = 6

	CASE(7) !CATENARY CABLE 
C	-- 1 FOR RESULTANT FORCE 
C	-- 6 STORE INITIAL NODAL FORCES 
C	-- 1 FOR FLAG TO STORE THE INITIAL FORCE
C	-- 1 FOR INITIAL LENGTH
	NWG = 1 + 6 + 1 + 1 
      
      CASE(8) !CATENARY CABLE 
C	-- 1 FOR RESULTANT FORCE 
C	-- 6 STORE INITIAL NODAL FORCES 
C	-- 1 FOR FLAG TO STORE THE INITIAL FORCE
C	-- 1 FOR INITIAL LENGTH
	NWG = 1 + 6 + 1 + 1 

	ENDSELECT



      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK5 (ISTYP,MTMOD,NWG,IEG,NGPS,ISECT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)


	NWG = 9

C	---------------------------------------------
C	MAXIMUM NUMBER OF FIBER FOR SECTION
C	---------------------------------------------
	IF(ISECT.EQ.3) THEN
	NFMAX   = 0
	DO I = 1,NGPS
	CALL XFSECTION(IEG,I,0)
	INAME(1:4) = [5,0,1,IEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FNFIB,1,4 ,0) !
	NFIB = INT(FNFIB)
	IF(NFIB.GE.NFMAX) NFMAX = NFIB
	ENDDO
C	==============================
	ENDIF
C	---------------------------------------------



C	------------------------------
C	------------------------------
	SELECTCASE(ISTYP)


	CASE(5) !FRAME ELEMENT

C	------------------------------
	SELECTCASE(MTMOD)

	CASE(1)
	NWG = 7

	CASE(3)
	IF(ISECT.EQ.3) NWG  = 3*NFMAX

	CASE(5)
	IF(ISECT.EQ.3) NWG  = 8				!FIBER STRAIN
	IF(ISECT.EQ.3) NWG  = NWG + 3*NFMAX	!FIBER STRESS, MAXIMUM STRESS, FLAG

	ENDSELECT
C	------------------------------

	ENDSELECT
C	------------------------------
C	------------------------------



      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK6 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	
	NWG = 4

	SELECTCASE(ISTYP)

	CASE(0,1,2) !8 NODE
	IF(MTMOD.EQ.1) NWG =  4          !4 STS 
	IF(MTMOD.EQ.3) NWG =  4 + 4 + 1  !4 STS 4 STN 1 PLS
	IF(MTMOD.EQ.6) NWG =  4 + 4 + 1 + 2  !4 STS 4 STN 1 PLS 2 COULOMB (THICKNESS 1 STS 1 STN)

	CASE(3,4,5) !4 EAS
	IF(MTMOD.EQ.1) NWG =  4          !4 STS 
	IF(MTMOD.EQ.3) NWG =  4 + 4 + 1  !4 STS 4 STN 1 PLS
	IF(MTMOD.EQ.4) NWG =  4 + 4 + 1 + 2  !4 STS 4 STN 1 PLS 2 COULOMB (THICKNESS 1 STS 1 STN)
	IF(MTMOD.EQ.5) NWG =  4 + 4 + 1 + 2  !4 STS 4 STN 1 PLS 2 DRUCKER (THICKNESS 1 STS 1 STN)

	CASE(6) !ZERO THICKNESS INTERFACE MEMBRANE
	NWG = 5

	ENDSELECT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK8 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	NWG = 4

	SELECTCASE(ISTYP)

	CASE(3,4,5) !4 EAS
	IF(MTMOD.EQ.1) NWG =  4          !4 STS 
	IF(MTMOD.EQ.3) NWG =  4 + 4 + 1  !4 STS 4 STN 1 PLS
	IF(MTMOD.EQ.4) NWG =  4 + 4 + 1 + 2  !4 STS 4 STN 1 PLS 2 COULOMB (THICKNESS 1 STS 1 STN)
	IF(MTMOD.EQ.5) NWG =  4 + 4 + 1 + 2  !4 STS 4 STN 1 PLS 2 DRUCKER (THICKNESS 1 STS 1 STN)

	ENDSELECT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK9 (ISTYP,MTMOD,NWG,NGT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	NWG = 9

	SELECTCASE(ISTYP)

	CASE(1) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)
	IF(MTMOD.EQ.5) NWG = 20*NGT
	IF(MTMOD.EQ.6) NWG = 43*NGT

	CASE(2) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8  + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)

	CASE(3) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8  + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)

	CASE(4) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8  + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)
	IF(MTMOD.EQ.5) NWG = 20*NGT
	IF(MTMOD.EQ.6) NWG = 43*NGT

	CASE(8) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)
	IF(MTMOD.EQ.5) NWG = 20*NGT
	IF(MTMOD.EQ.6) NWG = 43*NGT
	
	CASE(9) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)
	IF(MTMOD.EQ.5) NWG = 20*NGT
	IF(MTMOD.EQ.6) NWG = 43*NGT
	
	CASE(10) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8  + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)
	
	CASE(11) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8  + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)

	CASE(12) 
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8  + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)
      
      CASE(13) !SHEAR WALL SHELL 4 EAS BY BJ
	IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG =  8 + 8 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.3) NWG =  8 + 8 + 1 + 1 !8 STRESS 8 STRAIN 1 PLASTIC STRAIN 1 YIELD
	IF(MTMOD.EQ.4) NWG =  9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 PLS)
	IF(MTMOD.EQ.5) NWG = 20*NGT
	IF(MTMOD.EQ.6) NWG = 43*NGT

	ENDSELECT

	NWG = NWG + 18  !18 FOR STORE 2 SET OF LOCAL AXIS FOR EACH GAUSS POINT JAN09  SHELL LAX

      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK10 (ISTYP,MTMOD,NWG,NGT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


	NWG = 6

	SELECTCASE(ISTYP)

      CASE(1) 
          
      NWG = 9    
          
      CASE(4) !FOR SOLID-SHELL(LIU) BY BJ
	IF(MTMOD.EQ.1) NWG = 9
      !IF(MTMOD.EQ.1) NWG =  9 + (6*3)     !6-STRESS FOR 3 LAYER (TOP MID BOT)
	IF(MTMOD.EQ.2) NWG = 9 + 9 + (9*NGT)   !9 RESULTANT + NGT*(4 STS 4 STN 1 DAMAGE FLAG)
	IF(MTMOD.EQ.4) NWG = 9 + 7*NGT
	IF(MTMOD.EQ.5) NWG = 20*NGT
	IF(MTMOD.EQ.6) NWG = 43*NGT

      NWG = NWG + 18  !18 FOR STORE 2 SET OF LOCAL AXIS FOR EACH GAUSS POINT JAN09  SHELL LAX
      
	CASE(6) 
	IF(MTMOD.EQ.1) NWG =  6          !6 STS 
	IF(MTMOD.EQ.3) NWG =  6 + 6 + 1  !6 STS 6 STN 1 PLS
	IF(MTMOD.EQ.4) NWG =  6 + 6 + 1  !6 STS 6 STN 1 MOHR
	IF(MTMOD.EQ.5) NWG =  6 + 6 + 1  !6 STS 6 STN 1 DRUK
	IF(MTMOD.EQ.6) NWG =  32         !DAMAGE

	CASE(8) 
	IF(MTMOD.EQ.1) NWG =  6          !6 STS 
	IF(MTMOD.EQ.3) NWG =  6 + 6 + 1  !6 STS 6 STN 1 PLS
	IF(MTMOD.EQ.6) NWG =  32         !DAMAGE


	CASE(9) 
	NWG = 7


	CASE(13) 
	IF(MTMOD.EQ.1) NWG =  6          !6 STS 
	IF(MTMOD.EQ.3) NWG =  6 + 6 + 1  !6 STS 6 STN 1 PLS
	IF(MTMOD.EQ.4) NWG =  6 + 6 + 1  !6 STS 6 STN 1 MOHR COEF
	IF(MTMOD.EQ.5) NWG =  6 + 6 + 1  !6 STS 6 STN 1 DRUK
	IF(MTMOD.EQ.6) NWG =  32         !DAMAGE
	IF(MTMOD.EQ.7) NWG =  6 + 6 + 1 + 6 !6 STS 6 STN 1 PLS 6 VISCO VECTOR


	ENDSELECT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK11 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


	NWG = 6

	SELECTCASE(ISTYP)

	CASE(2) 
	IF(MTMOD.EQ.1) NWG = 6
	IF(MTMOD.EQ.3) NWG = 6 + 6 + 1      !6 STS 6 STN 1 MISE COEF
	IF(MTMOD.EQ.4) NWG = 6 + 6 + 1		!6 STS 6 STN 1 MOHR COEF
C	IF(MTMOD.EQ.4) NWG = 6 + 6 + 1 + 6  !6 STS 6 STN 1 PLS 6 VISCO VECTOR

	ENDSELECT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK13 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	NWG = 2     !2 VELOCITY COMPONENT FOR SEEPAGE

	SELECTCASE(ISTYP)

	CASE(1) 
	NWG = 2

	ENDSELECT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK14 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


	NWG = 3     !3 VELOCITY COMPONENT FOR SEEPAGE

	SELECTCASE(ISTYP)

	CASE(1) 
	NWG = 3

	ENDSELECT


      RETURN

      END


C	====================================================================
C	====================================================================
C	====================================================================
      SUBROUTINE INWOK15 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	NWG = 2     

	SELECTCASE(ISTYP)

	CASE(1) 
	NWG = 2

	ENDSELECT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK16 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


	NWG = 3     

	SELECTCASE(ISTYP)

	CASE(1) 
	NWG = 3

	ENDSELECT


      RETURN

      END


C	====================================================================
C	====================================================================
C	====================================================================

      SUBROUTINE INWOK20 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	
	NWG = 4

	SELECTCASE(ISTYP)

	CASE(0,1,2) !4 EAS
	NWG = 4          

	ENDSELECT


      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INWOK17 (ISTYP,MTMOD,NWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	TENDON WORKING ARRAY
	NWG = 20

      RETURN

      END

C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE EASCOM(IEG,IFLAG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	====================================================================
C	LOAD MM AND MM1 FOR ENHANCE ANALYSIS  SONGSAK MAR2006
C	====================================================================
	ALLOCATABLE A(:)
C	=============================================================

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /MMENH/ MM,MM1,MM2,NDIMC


	MM  = 1
	MM1 = 1

	MM2   = 0  !FOR CONSOLIDATION SACHARUCK DEC2006
	NDIMC = 3  !SET PROBLEM DIMENSION = 3D FOR CONSOLIDATION

	IF(ITYPE.EQ.6) THEN ! MEMBRANE
	SELECT CASE(ISTYP)
	CASE(3)   !AXIS SYMMETRIC
	MM  = 5
	MM1 = 8
	CASE(4,5) !PLANE ELEMENT
	MM  = 7
	MM1 = 8
	END SELECT
	NDIMC = 2  !SET PROBLEM DIMENSION = 2D
	ENDIF


C	ITYPE = 8 MEMBRANE EAS CONSOLIDATION	
	IF(ITYPE.EQ.8) THEN 
	SELECT CASE(ISTYP)
	CASE(3)   !AXIS SYMMETRIC
	MM  = 5
	MM1 = 8
	MM2 = 4
	CASE(4,5) !PLANE ELEMENT
	MM  = 7
	MM1 = 8
	MM2 = 4
	END SELECT
	NDIMC = 2  !SET PROBLEM DIMENSION = 2D
	ENDIF

	IF(ITYPE.EQ.9) THEN ! SHELL 4 NODE 
	SELECT CASE(ISTYP)
	CASE(1) !SHELL 4 EAS - SONG 7 TERM
	MM  = 7
	MM1 = 24
	CASE(8) !SHELL 4 EAS - BO
	MM  = 7
	MM1 = 24
	CASE(9) !SHELL 4 EAS - SONG 5 TERM
	MM  = 5
	MM1 = 24
      CASE(13) !SHEAR WALL SHELL 4 EAS - SONG 7 TERM BY BJ
	MM  = 7
	MM1 = 24
	END SELECT	
	NDIMC = 3  !SET PROBLEM DIMENSION = 2D
	ENDIF

	IF(ITYPE.EQ.10) THEN ! SOLID
	IF(ISTYP.EQ.2) THEN
	MM  = 5
	MM1 = 24
	ELSEIF(ISTYP.EQ.3) THEN
	MM  = 21
	MM1 = 24	
	ELSEIF(ISTYP.EQ.4) THEN !SOLID-SHELL(LIU) BY BJ
	MM  = 7
	MM1 = 24	
	ELSEIF(ISTYP.EQ.1.OR.ISTYP.EQ.6.OR.ISTYP.EQ.7) THEN !SOLID-SHELL(P) BY BJ
	MM  = 24
	IF(NLOPT.GT.0) MM = 24 !FOR NONLINEAR REDUCE EAS TERM SONGSAK NOV2005
	MM1 = 24
	ELSEIF(ISTYP.EQ.10) THEN
	MM  = 6
	MM1 = 24
	ELSEIF(ISTYP.EQ.12) THEN !CONSOLIDATION
	MM  = 24
	MM1 = 24  
	MM2 = 8	
	ENDIF	
	NDIMC = 3  !SET PROBLEM DIMENSION = 3D	
	ENDIF

C	ADDED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
C	SOLID CONTINUUM EAS CONSOLIDATION
	IF(ITYPE.EQ.11.AND.ISTYP.EQ.2)THEN
	MM   = 24
	MM1  = 24  
	MM2  = 8
	NDIMC = 3  !SET PROBLEM DIMENSION = 3D
	END IF


	IF(ITYPE.EQ.20) THEN ! MEMBRANE FLUID EAS
	SELECT CASE(ISTYP)
	CASE(0)   !AXIS SYMMETRIC
	MM  = 2
	MM1 = 8
	CASE(1,2) !PLANE ELEMENT
	MM  = 3
	MM1 = 8
	END SELECT
	NDIMC = 2  !SET PROBLEM DIMENSION = 2D
	ENDIF


	IF(IFLAG.EQ.2) RETURN

	LNALI = 1
	LNSEL = LNALI  + MM
	LNALP = LNSEL  + MM*MM1
	LNINF = LNALP  + MM
	LNDIS = LNINF  + MM
	LNCON = LNDIS  + (MM1+MM2)   !CONSO  + MM2 (FLUID DOF TERMS)
	LNINP = LNCON  + MM*MM2		 !CONSO
	LATES = LNINP  + MM
	NEAS  = LATES


	NBIT  = 2*NEAS
	NFL0  = 3400
	KWREC = NFL0 + IEG
	CALL DIROPEN(KWREC,NBIT) 

	CALL MINTFIL('GEAS',NELE ,IEG,1 ,1)
	CALL MINTFIL('GEAS',NEAS ,IEG,2 ,1)
	CALL MINTFIL('GEAS',LNALI,IEG,3 ,1)
	CALL MINTFIL('GEAS',LNSEL,IEG,4 ,1)
	CALL MINTFIL('GEAS',LNALP,IEG,5 ,1)
	CALL MINTFIL('GEAS',LNINF,IEG,6 ,1)
	CALL MINTFIL('GEAS',LNDIS,IEG,7 ,1)
	CALL MINTFIL('GEAS',LNCON,IEG,8 ,1)
	CALL MINTFIL('GEAS',LNINP,IEG,9 ,1)
	CALL MINTFIL('GEAS',LATES,IEG,10,1)

	ALLOCATE(A(NEAS))

	A(1:NEAS) = 0.0D0
	DO IELE = 1,NELE
	IRC = IELE 
	WRITE(KWREC,REC=IRC) A(1:NEAS)
	ENDDO 

	DEALLOCATE(A)

      RETURN

      END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE EASCALL(ASEL,ASEQ,AINF,AINP,APHA,APHI,ADIS,EDIS,ALPHA,
	1				   ALPHI,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	====================================================================
C	MODULE FOR CALLING THE EAS PARAMETER
C	OUTPUT VARIABLE (ALPHA,ALPHI)
C	====================================================================
	COMMON /MMENH/ MM,MM1,MM2,NDIMC
C	====================================================================
	DIMENSION DL(MM,MM1),ALPHA(MM),ALPHH(MM)
	DIMENSION DQ(MM,MM2),ALPHI(MM),ALPHP(MM)

	DIMENSION EDIS(MM1+MM2),EADSI(MM1+MM2)
	DIMENSION EDIS1(MM1),EDIS2(MM2)
	DIMENSION EDIST1(MM1),EDIST2(MM2)

	DIMENSION ASEL(MM,MM1),ASEQ(MM,MM2),APHA(MM),
	1		  AINF(MM),AINP(MM),ADIS(MM1+MM2),
	2		  APHI(MM)

	DIMENSION ALPHA1(MM),ALPHA2(MM),ALPHA3(MM)
	DIMENSION ALPHI1(MM),ALPHI2(MM),ALPHI3(MM)
C	====================================================================
	DO MMI = 1,MM
	DO MMJ = 1,MM1
	DL(MMI,MMJ)  = ASEL(MMI,MMJ)
	ENDDO
	ALPHH(MMI)   = AINF(MMI)
	ENDDO

	DO MMI = 1,MM1+MM2                         !INCREMENTAL
	EADSI(MMI)   = EDIS(MMI) - ADIS(MMI)
	ENDDO

	IF(IND.EQ.2) GOTO 100
C	----------------------------------------------------
C	STANDARD EAS
C	----------------------------------------------------
	DO MMI = 1,MM
	ALPHA(MMI) = ALPHH(MMI) + APHA(MMI)
	DO MMJ = 1,MM1
	ALPHA(MMI) = ALPHA(MMI) + DL(MMI,MMJ)*EADSI(MMJ)
	ENDDO
	ENDDO 
C	----------------------------------------------------
	IF(IND.EQ.1) GOTO 200
C	----------------------------------------------------

100	CONTINUE

C	----------------------------------------------------
C	COUPLING EAS (DISPLACEMENT AND PRESSURE)
C	----------------------------------------------------
	DO MMI = 1,MM  
	DO MMJ = 1,MM2
	DQ(MMI,MMJ)  = ASEQ(MMI,MMJ)
	ENDDO
	ALPHP(MMI)   = AINP(MMI)
	ENDDO
C	----------------------------------------------------
	EDIS1  = 0.0
	EDIS2  = 0.0
	EDIST1 = 0.0
	EDIST2 = 0.0
C	----------------------------------------------------
	IF(NDIMC.EQ.3) THEN               !SOLID 3D
	CALL EXTRACTDISP3D(EADSI,EDIS1)   !INCREMENTAL
	CALL EXTRACTPRES3D(EADSI,EDIS2)
	CALL EXTRACTDISP3D(EDIS,EDIST1)   !TOTAL
	CALL EXTRACTPRES3D(EDIS,EDIST2)
	ELSEIF(NDIMC.EQ.2) THEN           !MEMBRANE 2D
	CALL EXTRACTDISP  (EADSI,EDIS1)   !INCREMENTAL
	CALL EXTRACTPRES  (EADSI,EDIS2)
	CALL EXTRACTDISP  (EDIS,EDIST1)   !TOTAL
	CALL EXTRACTPRES  (EDIS,EDIST2)
	ENDIF
C	----------------------------------------------------
	ALPHA1 = 0.0
	ALPHA2 = 0.0
	ALPHA3 = 0.0
	ALPHI1 = 0.0
	ALPHI2 = 0.0
	ALPHI3 = 0.0
C	----------------------------------------------------
	DO MMI = 1,MM
	ALPHA1(MMI) = APHA(MMI)
	ALPHI1(MMI) = APHI(MMI)
	DO MMJ = 1,MM1
	ALPHA2(MMI) =  ALPHA2(MMI) + DL(MMI,MMJ)*EDIS1(MMJ)  !INCREMENTAL
	ALPHI2(MMI) =  ALPHI2(MMI) + DL(MMI,MMJ)*EDIS1(MMJ)  !INCREMENTAL
C	ALPHI2(MMI) =  ALPHI2(MMI) + DL(MMI,MMJ)*EDIST1(MMJ) !TOTAL
	ENDDO
	DO MMJ = 1,MM2
	ALPHA3(MMI) =  ALPHA3(MMI) + DQ(MMI,MMJ)*EDIS2(MMJ)  !INCREMENTAL
	ALPHI3(MMI) =  ALPHI3(MMI) + DQ(MMI,MMJ)*EDIS2(MMJ)  !INCREMENTAL
C	ALPHI3(MMI) =  ALPHI3(MMI) + DQ(MMI,MMJ)*EDIST2(MMJ) !TOTAL
	ENDDO
	ENDDO 
C	----------------------------------------------------
	ALPHA = ALPHA1 + (ALPHH+ALPHP+ALPHA2+ALPHA3)  !INCREMENTAL
	ALPHI = ALPHI1 + (ALPHH+ALPHP+ALPHI2+ALPHI3)  !INCREMENTAL
C	ALPHI =        + (ALPHI2+ALPHI3)              !TOTAL  *REMARK* THIS INTRODUCED BY SONGSAK REFER TO THE EAS FORMULATION OF HINFC AND HINFP
C	----------------------------------------------------
200	CONTINUE

C
C	STORE THE UPDATED TOTAL ENHANCED STRAIN PARAMETERS
	DO MMI = 1,MM
	APHA(MMI) = ALPHA(MMI)
	APHI(MMI) = ALPHI(MMI)
	ENDDO

      RETURN

      END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE EASSTOR(SEDI,SEL,SEQ,HINFC,HINFP,ASEL,ASEQ,AINF,AINP,
	1				   ADIS,EDIS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	====================================================================
C	MODULE FOR STORING THE EAS PARAMETER
C	====================================================================
C	====================================================================
	COMMON /MMENH/ MM,MM1,MM2,NDIMC
C	====================================================================
	DIMENSION DL(MM,MM1),HINFC(MM),HINFP(MM)
	DIMENSION EDIS(MM1+MM2)
	DIMENSION DQ(MM,MM2)
	DIMENSION SEDI(MM,MM),SEL(MM,MM1),SEQ(MM,MM2)

	DIMENSION ASEL(MM,MM1),ASEQ(MM,MM2),
	1		  AINF(MM),AINP(MM),ADIS(MM1+MM2)

C	====================================================================
	DL = MATMUL(SEDI,SEL)
      DQ = MATMUL(SEDI,SEQ)

	DO MMI = 1,MM       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	AINF(MMI) = 0.0D0
	AINP(MMI) = 0.0D0
	DO MMJ = 1,MM
	AINF(MMI) = AINF(MMI) + SEDI(MMI,MMJ)*HINFC(MMJ)
	AINP(MMI) = AINP(MMI) + SEDI(MMI,MMJ)*HINFP(MMJ)
	ENDDO

	DO MMJ = 1,MM1
	ASEL(MMI,MMJ) = DL(MMI,MMJ)
	ENDDO

	DO MMJ = 1,MM2
	ASEQ(MMI,MMJ) = DQ(MMI,MMJ)
	ENDDO

	ENDDO              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	ELEMENT TOTAL DISPLACMENT
	DO MMI = 1,MM1+MM2
	ADIS(MMI) = EDIS(MMI)
	ENDDO



      RETURN

      END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE EASTYP(ITYPE,ISTYP,IESPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	====================================================================
C	MODULE FOR CALL THE EAS FLAG 
C	====================================================================
C	====================================================================

	IESPT = 0

	IF(ITYPE.EQ.10)    THEN !SOLID EAS

	SELECT CASE (ISTYP)
	CASE(1,2,3,4,6,7,10,12) !BY BJ -- ADD SOLID-SHELL(P) 
	IESPT = 1
	END SELECT ! FOR SOLID EAS

	ELSEIF(ITYPE.EQ.6) THEN !MEMBRANE EAS

	SELECT CASE (ISTYP)
	CASE(3,4,5)  
	IESPT = 1
	END SELECT

	ELSEIF(ITYPE.EQ.8) THEN !MEMBRANE EAS CONSOLIDATION

	SELECT CASE (ISTYP)
	CASE(3,4,5)  
	IESPT = 2
	END SELECT

	ELSEIF(ITYPE.EQ.9) THEN !SHELL EAS

	SELECT CASE (ISTYP)
	CASE(1,8,9,10,13)  
	IESPT = 1
	END SELECT

	ELSEIF(ITYPE.EQ.11) THEN !SOLID EAS CONSOLIDATION

	SELECT CASE (ISTYP)
	CASE(2)  
	IESPT = 2
	END SELECT

	ELSEIF(ITYPE.EQ.20) THEN !MEMBRANE EAS FLUID

	SELECT CASE (ISTYP)
	CASE(0,1,2)  
	IESPT = 1
	END SELECT

	ENDIF

      RETURN

      END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE ADREAS(IEG,IELE,APHI,ASEL,APHA,AINF,ADIS,ASEQ,AINP,
	1				  OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	CHARACTER*3 OPER
	DIMENSION APHI(1),ASEL(1),APHA(1),AINF(1),ADIS(1),ASEQ(1),AINP(1)
	ALLOCATABLE A(:)
C	=============================================================
	IF(OPER.EQ.'RED') IND = 1
	IF(OPER.EQ.'WRT') IND = 2

	NFL0  = 3400
	KWREC = NFL0 + IEG

	CALL MINTFIL('GEAS',NELE ,IEG,1 ,0)
	CALL MINTFIL('GEAS',NEAS ,IEG,2 ,0)
	CALL MINTFIL('GEAS',LNALI,IEG,3 ,0)
	CALL MINTFIL('GEAS',LNSEL,IEG,4 ,0)
	CALL MINTFIL('GEAS',LNALP,IEG,5 ,0)
	CALL MINTFIL('GEAS',LNINF,IEG,6 ,0)
	CALL MINTFIL('GEAS',LNDIS,IEG,7 ,0)
	CALL MINTFIL('GEAS',LNCON,IEG,8 ,0)
	CALL MINTFIL('GEAS',LNINP,IEG,9 ,0)
	CALL MINTFIL('GEAS',LATES,IEG,10,0)

	ALLOCATE(A(NEAS))

	IRC   = IELE 

	SELECTCASE(IND)

	CASE(1)
	 READ(KWREC,REC=IRC) A(1:NEAS)

      N = (LNSEL-1) - LNALI + 1
      IF(N.GT.0) APHI(1:N) = A(LNALI:LNSEL-1)
      
      N = (LNALP-1) - LNSEL + 1
      IF(N.GT.0) ASEL(1:N) = A(LNSEL:LNALP-1)
      
      N = (LNINF-1) - LNALP + 1
      IF(N.GT.0) APHA(1:N) = A(LNALP:LNINF-1)
      
      N = (LNDIS-1) - LNINF + 1
      IF(N.GT.0) AINF(1:N) = A(LNINF:LNDIS-1)
      
      N = (LNCON-1) - LNDIS + 1
      IF(N.GT.0) ADIS(1:N) = A(LNDIS:LNCON-1)
      
      N = (LNINP-1) - LNCON + 1
      IF(N.GT.0) ASEQ(1:N) = A(LNCON:LNINP-1)
      
      N = (LATES-1) - LNINP + 1
      IF(N.GT.0) AINP(1:N) = A(LNINP:LATES-1)

	CASE(2)

      N = (LNSEL-1) - LNALI + 1 
      IF(N.GT.0) A(LNALI:LNSEL-1) = APHI(1:N)
      
      N = (LNALP-1) - LNSEL + 1
      IF(N.GT.0) A(LNSEL:LNALP-1) = ASEL(1:N)
      
      N = (LNINF-1) - LNALP + 1
      IF(N.GT.0) A(LNALP:LNINF-1) = APHA(1:N)
      
      N = (LNDIS-1) - LNINF + 1
      IF(N.GT.0) A(LNINF:LNDIS-1) = AINF(1:N)
      
      N = (LNCON-1) - LNDIS + 1
      IF(N.GT.0) A(LNDIS:LNCON-1) = ADIS(1:N)
      
      N = (LNINP-1) - LNCON + 1
      IF(N.GT.0) A(LNCON:LNINP-1) = ASEQ(1:N)
      
      N = (LATES-1) - LNINP + 1
      IF(N.GT.0) A(LNINP:LATES-1) = AINP(1:N)
      
	WRITE(KWREC,REC=IRC) A(1:NEAS)

	ENDSELECT

	DEALLOCATE(A)

	RETURN
      
      
C     BACKUP OLD CODE
	SELECTCASE(IND)

	CASE(1)
	 READ(KWREC,REC=IRC) A(1:NEAS)

	K = 0
	DO I = LNALI,LNSEL-1
	K = K + 1
	APHI(K) = A(I)
	ENDDO
	K = 0
	DO I = LNSEL,LNALP-1
	K = K + 1
	ASEL(K) = A(I)
	ENDDO
	K = 0
	DO I = LNALP,LNINF-1
	K = K + 1
	APHA(K) = A(I)
	ENDDO
	K = 0
	DO I = LNINF,LNDIS-1
	K = K + 1
	AINF(K) = A(I)
	ENDDO
	K = 0
	DO I = LNDIS,LNCON-1
	K = K + 1
	ADIS(K) = A(I)
	ENDDO
	K = 0
	DO I = LNCON,LNINP-1
	K = K + 1
	ASEQ(K) = A(I)
	ENDDO
	K = 0
	DO I = LNINP,LATES-1
	K = K + 1
	AINP(K) = A(I)
	ENDDO

	CASE(2)

	K = 0
	DO I = LNALI,LNSEL-1
	K = K + 1
	A(I) = APHI(K)
	ENDDO
	K = 0
	DO I = LNSEL,LNALP-1
	K = K + 1
	A(I) = ASEL(K)
	ENDDO
	K = 0
	DO I = LNALP,LNINF-1
	K = K + 1
	A(I) = APHA(K)
	ENDDO
	K = 0
	DO I = LNINF,LNDIS-1
	K = K + 1
	A(I) = AINF(K)
	ENDDO
	K = 0
	DO I = LNDIS,LNCON-1
	K = K + 1
	A(I) = ADIS(K)
	ENDDO
	K = 0
	DO I = LNCON,LNINP-1
	K = K + 1
	A(I) = ASEQ(K)
	ENDDO
	K = 0
	DO I = LNINP,LATES-1
	K = K + 1
	A(I) = AINP(K)
      ENDDO

	WRITE(KWREC,REC=IRC) A(1:NEAS)

      ENDSELECT

      
	END

C	=============================================================
C	=============================================================
C	=============================================================