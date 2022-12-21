

C	===================================
C	SOLID CONSOLIDATION ELEMENT LIBRARY
C	ISTYP = 11
C	==========

      SUBROUTINE SOLIDCON3(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,
     +		   MWG,ALPHA,SEL,SEDI,FIN,HINFC,MSET,SEQ,HINFP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     ===================================
C     CONSOLIDATION SOLID ELEMENT LIBRARY
C	===================================

C	  ISTYP	  	  ELEMENT DESCRIPTION
C	  =====		  ===================
 
C	ISTYP = 1	COMPATIBLE SOLID CONSOLIDATION ELEMENT
C	ISTYP = 2	ENHANCED ASSUMED STRAIN CONSOLIDATION SOLID ELEMENT
C	ISTYP = 3   -
C	ISTYP = 4	-

C     ======================================================================

      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


C	EAS ARRAY SONGSAK JAN2007
	COMMON /MMENH/ MM,MM1,MM2,NDIMC

C	=================================================
C	ENHANCED ASSUMED STRAIN SOLID ELEMENT DESCRIPTION
C	=================================================

C	  VARIABLE       DESCRIPTION
C	  ========		 ===========

C	    MM1			 24/DOF 
C		MM2			  8/DOF
C	  	 MM			SETTING 24 PARAMETERS
C	=====================================

      DIMENSION PROPM(1),PROPG(1),NODEX(1),WA(1),S(1),COORD(1),EDIS(1)
      DIMENSION EDISI(1),ELOD(1),FIN(1)

C	=====================================
C	ENHANCED ASSUMED STRAIN VECTOR/MATRIX
C	=====================================

	DIMENSION ALPHA(MM,1),SEL(MM,24),SEDI(MM,MM),RH(MM,1),
	1		  HINFC(MM)
	DIMENSION SEQ(MM,MM2),HINFP(MM) !CONSO


	GOTO (10,20,30,40),ISTYP

 
10	CONTINUE

C	======================================
C	COMPATIBLE SOLID CONSOLIDATION ELEMENT
C	======================================

C	CALL SOLIDCONSO(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,
C	1				RE,NWG,SDIS,FLVEC)

	RETURN

20	CONTINUE

C	===================================================
C	ENHANCED ASSUMED STRAIN CONSOLIDATION SOLID ELEMENT
C	===================================================
	CALL CONSOLEAS3D(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,
	1            RE,NWG,ALPHA,SEL,SEDI,FIN,SEQ,HINFC,HINFP)


	RETURN

30	CONTINUE


	RETURN

40	CONTINUE



	RETURN


	END
