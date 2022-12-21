C=========================================================================
C	START BEAM ELEMENT ROUTINES
C=====================================================================
      SUBROUTINE BEAM21 (PROPM,PROPG,PROPO,WA,S,COORD,
	1			 	EDIS,EDISI,RE,MWG,FIN,IPIN,
     2				MCF,ISET,FIXEN,FIXLR,FIXEO,FIXLO,PMATRL,
     3                FIXEN_OFF,FIXLR_OFF,FIXEO_OFF,FIXLO_OFF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------------
C     A 2 & 3 NODED SPACE BEAM ELEMENT WITH 7 D.O.F. PER NODE FOR LARGE
C     DEFLECTION ELASTO-PLASTIC ANALYSIS OF SPACE FRAMES, BEAM COLUMNS
C     OR PLATE/SHELL STIFFENERS HAVING THIN-WALLED OPEN CROSS-SECTIONS
C     OF ARBITRARY SHAPE  (2003). 
C     ------------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      DIMENSION PROPM(*),PROPG(*),WA(1),S(105),COORD(3,2),EDIS(14)
      DIMENSION EDISI(*),RE(14),DR(49),PROPO(*)
	DIMENSION FIN(1)
	DIMENSION FIXEN(NEF),FIXLR(NEF),FIXEO(NEF),FIXLO(NEF),IPIN(14)

	DIMENSION PMATRL(1)




	IPI = 0 !NO END RELEASE IN 104,107
C	---------------------------------------------------------
      GO TO (101,102,103,104,105,105,107) ISTYP

C	---------------------------------------------------------
C	DIRECT STIFFNESS METHORDS IN FRAME ANALYSIS (BY DE SILVA)
C	---------------------------------------------------------
 104  CONTINUE
	RETURN
C	-----------------------------------
C	HIGH ORDER BEAM ELEMENT (BY GILSON)
C	-----------------------------------
 103  CONTINUE							
	RETURN
C	--------------------------
C	COMMON (OLD) BEAM ELEMENT
C	--------------------------
 102  CONTINUE					
	RETURN
C	----------------------------------------
C	LOW ORDER BEAM ELEMENT (BY VIVEK GUPTA)
C	----------------------------------------
 101  CONTINUE
	RETURN

 105  CONTINUE
	CALL BFRAMEW (PROPM,PROPG,S,COORD,EDIS,EDISI,RE,PROPO,
	1			  FIN,ISET,WA,FIXEN,FIXLR,FIXEO,FIXLO,IPIN,PMATRL,
     1              FIXEN_OFF,FIXLR_OFF,FIXEO_OFF,FIXLO_OFF)
C	DIRECT STIFFNESS FRAME (SPC CONCRETE) SONGSAK FEB2006
C	CALL BFRASPC (PROPM,PROPG,S,COORD,EDIS,EDISI,RE,PROPO,
C	1			  FIN,ISET,WA,FIXEN,FIXLR,FIXEO,FIXLO,IPIN,PMATRL)
	RETURN
C	---------------------------------------------------------
C	FRAME ELEMENT FOR CONSTRUCTION STAGE ANALYSIS
C	---------------------------------------------------------
107	CONTINUE
	RETURN

	END
C=====================================================================

C======================================================================
C	RANG1 ADDED BY DE SILVA TO INTRODUCE MEMBER ROTATION ABOUT BEAM AXIS-OCT2003
      SUBROUTINE BMRIGD (COORD,COORDI,EDIS,REDIS,VSO,VR,VS,VT,NNO,
     1                   NLOPT,EC,R,ELN,RANG1)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C     EVALUATES LOCAL DIRECTION COSINES AND MODIFIES TOTAL DISPLACEMENT
C     BY DEDUCTING RIGID BODY TRANSLATIONS AND ROTATION 
C	-----------------------------------------------------------------
C     COORD(3,NNO)  = CURRENT NODAL COORDINATES (X,Y,Z)
C     COORDI(3,NNO) = ORIGINAL NODAL COORDINATES (X,Y,Z)
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     REDIS(14)     = CO-ROTATIONAL FORM OF DISPLACEMENT VECTOR
C     VR,VS,VT      = DIRECTION COSINE VECTORS OF CURRENT R,S,T AXES
C                     (EACH 3)
C     ARC           = 2*PI
C     NNO           = NUMBER OF NODES
C     NLOPT         = NONLINEAR OPTION CODE
C     IWRP          = WARPING OPTION CODE
C     ----------------------------------------------------------------
      DIMENSION COORD(3,2),COORDI(3,2),EDIS(14),ROT(3)
      DIMENSION ROV(3),VRO(3),VSO(3),VTO(3),VR(3),VS(3),VT(3)
      DIMENSION VRM(3),VSM(3),VTM(3),VV(3)
      DIMENSION REDIS(14),TM(3,3),CDO(3,3),XYZ(3)
      DIMENSION EC(3),R(1,2),S(2),T(2),DUMMY(3)
C	FOLLOWING LINES ADDED BY DESILVA-OCT2003
	DIMENSION VVS(3)

	DIMENSION XYZI(3),CDOI(3,3)
C
C	PASS ANGLE OF ROTATION ABOUT r DIRECTION AND POSITIVE IF ROTATION FOLLOWS
C	RIGHT HAND RULE( RANG- ANGLE OF ROTAION(DEGREE) )
	
	RANG = RANG1

	RANG = RANG*0.01745329252
c      CSR=COS(RANG) !CHANGED MARCH2005 (INTEGRATE)
c      SNR=SIN(RANG)
      CSR = COSD(RANG1)
      SNR = SIND(RANG1)

      ARC=6.2831853071796

	
C     ------------------------------------
C     FIND ELEMENT CENTER COORDINATES - EC
C     ------------------------------------
      DO 30 I=1,3
   30	EC(I)=(COORDI(I,1)+COORDI(I,2))/2.0
C     ----------------------------------------
C     FIND DIRECTION COSINE VECTORS VR, VS, VT
C     ----------------------------------------
      DO 40 I=1,3
      VR(I)  = (COORD(I,2)-COORD(I,1))/2.
   40	VRO(I) = (COORDI(I,2)-COORDI(I,1))/2.


      CALL SCALEN (VR,VR,DUM,3)
      CALL SCALEN (VRO,VRO,DUM,3)
      CALL VECPRD (VRO,VSO,VTO)
      CALL VECPRD (VTO,VRO,VSO)
C	------------------------------------------------------------
C	ROTATION ABOUT THE AXIS OF THE BEAM-ADDED BY DESILVA OCT2003
C	------------------------------------------------------------
      DO 65 I=1,3
   65 VVS(I)=VSO(I)*CSR+VTO(I)*SNR

      CALL VECPRD (VRO,VVS,VTO)
      CALL VECPRD (VTO,VRO,VVS)

	DO 2 I=1,3
	VS(I)=VVS(I)
     	VT(I)=VTO(I)
    2 CONTINUE


      DO 50 I=1,2
      DO 60 J=1,3
      DUMMY(J)=COORDI(J,I)-EC(J)
   60	CONTINUE
	R(1,I)= DUMMY(1)*VRO(1)+DUMMY(2)*VRO(2)+DUMMY(3)*VRO(3)
	S(I)  = DUMMY(1)*VVS(1)+DUMMY(2)*VVS(2)+DUMMY(3)*VVS(3)
      T(I)  = DUMMY(1)*VTO(1)+DUMMY(2)*VTO(2)+DUMMY(3)*VTO(3)
   50 CONTINUE

C     --------------------------------------
C     FINDS ELN    = ELEMENT BOUNDARY LENGTH
C     --------------------------------------
      ELN=R(1,2)-R(1,1)

	IF (NLOPT.LE.1) RETURN
C	RETURN
C     -------------------------------------
C     FIND GLOBAL ROTATION COMPONENTS (ROT)
C     -------------------------------------
C     FIND ROTATIONAL COMPONENTS AT THE ELEMENT CENTER/REF POINT
      CALL CLEARA (ROT,3)
	DO 1 I=1,3
      ROT(I)=(EDIS(I+3)+EDIS(I+10))/2.0
   1	CONTINUE

C     ---------------------------------------------------------------
C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
C     AND HENCE DETERMINE MATERIALLY ATTACHED FRAME (VRM,VSM,VTM)
C     ---------------------------------------------------------------
      CALL SCALEN (ROT,ROV,RO,3)
      IF (RO.EQ.0.) RETURN
      CALL ROTVEC (VRO,ROV,RO,VRM)
      CALL ROTVEC (VVS,ROV,RO,VSM)

C     --------------------------------------------------------------
C     NOW DEFINE LOCAL FRAME (VR,VS,VT) BY DEDUCTING SHEAR ROTATIONS
C     --------------------------------------------------------------
      CALL SCAPRD (VRM,VR,CS,3)
      IF (CS.LT.0.999) GOTO 110
      CALL VECPRD (VRM,VSM,VTM)
      CALL VECPRD (VTM,VR,VS)
      CALL SCALEN (VS,VS,DUM,3)
      GOTO 130
  110 CALL VECPRD (VRM,VR,VV)
      CALL SCALEN (VV,VV,DUM,3)
      RSH=ACOS(CS)
      CALL ROTVEC (VSM,VV,RSH,VS)
  130 CALL VECPRD (VR,VS,VT)

C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
C     ELEMENT INITIAL CENTER CORD/INITIAL REF. POINT COORD - XYZ
      CALL CLEARA (XYZ,3)
	DO 11 I=1,3
       XYZ(I) = (COORDI(I,1) + COORDI(I,2))/2.
	XYZI(I) = ( COORD(I,1) +  COORD(I,2))/2.
   11	CONTINUE

C     ----------------------------------------------------------
C     INITIAL NODAL DISTANCE RELATIVE TO INITIAL REF POINT - CDO     
C     ----------------------------------------------------------
      DO 150 I=1,NNO
      DO 150 J=1,3
       CDO(J,I) = COORDI(J,I) -  XYZ(J)
  150 CDOI(J,I) =  COORD(J,I) - XYZI(J)

C     ------------------------------------
C     FIND ORTHOGONAL ROTATION MATRIX - TM
C     ------------------------------------
C	R = T*ToT
      DO 160 I=1,3
      VV(I)=ROV(I)
      DO 160 J=1,3
  160 TM(I,J)=VR(I)*VRO(J)+VS(I)*VVS(J)+VT(I)*VTO(J)


C     ---------------------------
C     FIND ANGLE OF ROTATION - RR
C     ---------------------------
C	========================================================
C	CHANGED BY SONGSAK TO CONTROL THE ROUND-OFF ERROR MAR2006
	RZETA = (.5*(TM(1,1)+TM(2,2)+TM(3,3)-1.0))
	IF(ABS(RZETA).GT.1.0) RZETA = RZETA/ABS(RZETA)
      RR = ACOS(RZETA)
C	========================================================

C     RR = ACOS(.5*(TM(1,1)+TM(2,2)+TM(3,3)-1.))

      SN=SIN(RR)
      IF (SN.EQ.0.) GOTO 190
C     ------------------------------
C     FIND ROTATION AXIS VECTOR - VV
C     ------------------------------
C
C			    [r32 - r23]
C	e =   1     [r13 - r31]
C		2sin(o) [r21 - r12]
C
C     ------------------------------
      F =0.5/SN
      VV(1) = F*(TM(3,2)-TM(2,3))
      VV(2) = F*(TM(1,3)-TM(3,1))
      VV(3) = F*(TM(2,1)-TM(1,2))
C     ------------------------------
C     ------------------------------

      CALL SCAPRD (ROV,VV,CS,3)
      IF (CS.GE.0.) GOTO 190

      RR=-RR
      DO 180 I=1,3
  180 VV(I)=-VV(I)

  190 RD=RO-RR
      N=(RD+SIGN(ARC/2.01,RD))/ARC
      RR=RR+N*ARC
C     --------------------------------
C     FIND CO-ROTATIONAL DISPLACEMENTS
C     --------------------------------
      K=1
      DO 220 I=1,NNO
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)

C	=================================
C	CHANGED BY SONGSAK MAR2006
      REDIS(K) = CDOI(J,I)  - TCDO
C	REDIS(K) = COORD(J,I) - TCDO
C	=================================

      REDIS(K+3) = EDIS(K+3) - RR*VV(J)
  210 K=K+1
  220 K=K+4

      RETURN

      END
C
C	=====================================================================

