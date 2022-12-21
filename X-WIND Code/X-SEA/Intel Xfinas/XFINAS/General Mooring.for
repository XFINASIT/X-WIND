	SUBROUTINE CABLE_LENGTH_ALL_TYPE (ALB,COOR_H,COOR_V,ALENGTH,IEL,OPT,COOR_REF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)   
      COMMON /STORE_MOOR/ MOORING_LINE(500,10),MOORISET(5000),MOORIELE(5000),INDEX_MOORING,INDEX_MOORING_ISET
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /MGRAV/ NGRAV

      ALLOCATABLE AX(:),AY(:),AZ(:),NODE(:)
      CHARACTER*4 OPT
      DIMENSION COOR_REF(3)
      
      ! --- INPUT ---
      ! IEL = LOOP NUMBER
      ! --- OUTPUT --
      ! COOR_H1
      ! COOR_V
      ! ALC = 
      NLINE = MOORISET(IEL)
      CALL INTFILL('#GMP',NEGEL,1,2,0)  ! TOTAL NUMBER OF ELEMENT
      
         ALLOCATE (AX(500),AY(500),AZ(500),NODE(500))  
         AX = 0.0D0
         AY = 0.0D0
         AZ = 0.0D0
         INDEX  = 0.0D0
         
         DO 100 I = 1,NELE
         IF (MOORING_LINE(I,NLINE).EQ.0) GOTO 100
         CALL CALLNUMNODE_F (MOORING_LINE(I,NLINE),NODEX,NODEY) 

         CALL RELFILL('@XYZ',AX1,1,NODEX,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',AY1,2,NODEX,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',AZ1,3,NODEX,0)  !GETTING NODAL COORDINATE	 
         
         CALL RELFILL('@XYZ',AX2,1,NODEY,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',AY2,2,NODEY,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',AZ2,3,NODEY,0)  !GETTING NODAL COORDINATE	
         
         INDEX = INDEX + 1 
         AX(INDEX)   = AX1
         AY(INDEX)   = AY1
         AZ(INDEX)   = AZ1
         NODE(INDEX) = NODEX
         INDEX = INDEX + 1 
         AX(INDEX)   = AX2
         AY(INDEX)   = AY2
         AZ(INDEX)   = AZ2
         NODE(INDEX) = NODEY
100      CONTINUE
         
         SELECTCASE (NGRAV)
         CASE (1) ! X-DIRECTION
            AMAX = MAXVAL(AX(1:INDEX))
            AMIN = MINVAL(AX(1:INDEX))
            NMAX_LOC = NODE(MAXLOC(AX(1:INDEX),DIM=1))
            NMIN_LOC = NODE(MINLOC(AX(1:INDEX),DIM=1))
         CASE (2) ! Y-DIRECTION
            AMAX = MAXVAL(AY(1:INDEX))
            AMIN = MINVAL(AY(1:INDEX))
            NMAX_LOC = NODE(MAXLOC(AY(1:INDEX),DIM=1))
            NMIN_LOC = NODE(MINLOC(AY(1:INDEX),DIM=1))
         CASE (3) ! Z-DIRECTION
            AMAX = MAXVAL(AZ(1:INDEX))
            AMIN = MINVAL(AZ(1:INDEX))
            NMAX_LOC = NODE(MAXLOC(AZ(1:INDEX),DIM=1))
            NMIN_LOC = NODE(MINLOC(AZ(1:INDEX),DIM=1))
         ENDSELECT
         DEALLOCATE (AX,AY,AZ,NODE)
         
         ! CALCULATE AS A LOCAL COORDINATE
         ! VESSEL LOCATION
         CALL RELFILL('@XYZ',COOR_X_N2,1,NMAX_LOC,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',COOR_Y_N2,2,NMAX_LOC,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',COOR_Z_N2,3,NMAX_LOC,0)  !GETTING NODAL COORDINATE
         
         CALL RELFILL('@XYZ',COOR_X_N1,1,NMIN_LOC,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',COOR_Y_N1,2,NMIN_LOC,0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',COOR_Z_N1,3,NMIN_LOC,0)  !GETTING NODAL COORDINATE
         
          ALENGTH = SQRT((COOR_X_N2-COOR_X_N1)**2+(COOR_Y_N2-COOR_Y_N1)**2+(COOR_Z_N2-COOR_Z_N1)**2)

          IF (NGRAV.EQ.1)     THEN ! X-DIRECTION
          COOR_V = ABS(COOR_X_N2 - COOR_X_N1)
          COOR_H = SQRT((COOR_Y_N2 - COOR_Y_N1)**2+(COOR_Z_N2 - COOR_Z_N1)**2)
          ELSEIF (NGRAV.EQ.2) THEN ! Y-DIRECTION
          COOR_V = ABS(COOR_Y_N2 - COOR_Y_N1)
          COOR_H = SQRT((COOR_X_N2 - COOR_X_N1)**2+(COOR_Z_N2 - COOR_Z_N1)**2)    
          ELSEIF (NGRAV.EQ.3) THEN ! Z-DIRECTION
          COOR_V = ABS(COOR_Z_N2 - COOR_Z_N1)
          COOR_H = SQRT((COOR_X_N2 - COOR_X_N1)**2+(COOR_Y_N2 - COOR_Y_N1)**2)
          ENDIF
          
          IF (OPT.EQ."CALC") THEN
          ALENGTH_ALB = ALENGTH - ALB
          RATIO =  1- ALENGTH_ALB/ALENGTH
          
          COOR_REF(1)  = COOR_X_N1 + (COOR_X_N2 - COOR_X_N1) * RATIO
          COOR_REF(2)  = COOR_Y_N1 + (COOR_Y_N2 - COOR_Y_N1) * RATIO
          COOR_REF(3)  = COOR_Z_N1 + (COOR_Z_N2 - COOR_Z_N1) * RATIO
          ENDIF
          
          
          


      
      END
      
! =================================================================     
      SUBROUTINE MOORING_FORCE (FAIRLEADH,FAIRLEADV,IEL,W,CB,ALB,COOR_H_INPUT,COOR_V_INPUT,HA,VA,COOR_H,OPT_CAL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)   
      CHARACTER*3 OPT_CAL
      COMMON /STORE_MOOR/ MOORING_LINE(500,10),MOORISET(5000),MOORIELE(5000),INDEX_MOORING,INDEX_MOORING_ISET
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /MGRAV/ NGRAV
      DIMENSION LINE_MOORING(500)
      DIMENSION AX(500),AY(500),AZ(500)
      DIMENSION AX_RE(500),AY_RE(500),AZ_RE(500)
      DIMENSION NMOORING_NODE(1000),NMOORING_NODE_RE(500)
      ! 
      ! IEL = NUMBER ELEMENT LOOP
      ! ALB = SEABED CONTRACT LENGTH
      
      NELEMENT_MOORING = MOORIELE(IEL)
      ISET             = MOORISET(IEL)
      INDEX            = 0.0D0
      
      DO I =1,500
          IF (MOORING_LINE(I,ISET).NE.0.0D0) INDEX = INDEX + 1
      ENDDO
      
      CALL CALLNUMNODE_F (NELEMENT_MOORING,NODE1,NODE2) 
      
      CALL RELFILL('@XYZ',COOR_X_N2,1,NODE2,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',COOR_Y_N2,2,NODE2,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',COOR_Z_N2,3,NODE2,0)  !GETTING NODAL COORDINATE
      
      CALL RELFILL('@XYZ',COOR_X_N1,1,NODE1,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',COOR_Y_N1,2,NODE1,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',COOR_Z_N1,3,NODE1,0)  !GETTING NODAL COORDINATE
      
      ! CHECK LOCATION OF THE LINE
      SELECTCASE (NGRAV)
      CASE (1)
      IF (COOR_X_N1.GT.COOR_X_N2) THEN
      AMIN = COOR_Z_N2  
      NMIN = NODE2
      AYY_MAX  = COOR_Y_N1
      AZZ_MAX  = COOR_Z_N1
      AYY_MIN  = COOR_Y_N2
      AZZ_MIN  = COOR_Z_N2
      ELSEIF (COOR_X_N2.GE.COOR_X_N1) THEN
      AMAX = COOR_X_N1  
      NMIN = NODE1 
      AYY_MAX  = COOR_Y_N2
      AZZ_MAX  = COOR_Z_N2
      AYY_MIN  = COOR_Y_N1
      AZZ_MIN  = COOR_Z_N1
      ENDIF
      CASE (2)
      IF (COOR_Y_N1.GT.COOR_Y_N2) THEN
      AMIN = COOR_Y_N2  
      NMIN = NODE2
      AXX_MAX  = COOR_X_N1
      AZZ_MAX  = COOR_Z_N1
      AXX_MIN  = COOR_X_N2
      AZZ_MIN  = COOR_Z_N2
      ELSEIF (COOR_Y_N2.GE.COOR_Y_N1) THEN
      AMAX = COOR_Y_N1  
      NMIN = NODE1    
      AXX_MAX  = COOR_X_N2
      AZZ_MAX  = COOR_Z_N2
      AXX_MIN  = COOR_X_N1
      AZZ_MIN  = COOR_Z_N1
      ENDIF
      CASE (3)
      IF (COOR_Z_N1.GT.COOR_Z_N2) THEN
      AMAX = COOR_Z_N1
      AMIN = COOR_Z_N2  
      NMIN = NODE2
      NMAX = NODE1
      AXX  = COOR_X_N1
      AYY  = COOR_Y_N1
      AXX_MIN  = COOR_X_N2
      AYY_MIN  = COOR_Y_N2
      ELSEIF (COOR_Z_N2.GE.COOR_Z_N1) THEN
      AMAX = COOR_Z_N2
      AMIN = COOR_Z_N1  
      NMIN = NODE1
      NMAX = NODE2
      AXX_MAX  = COOR_X_N2
      AYY_MAX  = COOR_Y_N2
      AXX_MIN  = COOR_X_N1
      AYY_MIN  = COOR_Y_N1
      ENDIF
      ENDSELECT
      
      ! CHECK COODINATE ON THE LINE
      INDEX_NODE = 0.0D0
      DO I = 1 ,INDEX
      CALL CALLNUMNODE_F (MOORING_LINE(I,ISET),N1,N2)
      INDEX_NODE               = INDEX_NODE + 1
      NMOORING_NODE(INDEX_NODE) = N1
      INDEX_NODE               = INDEX_NODE + 1
      NMOORING_NODE(INDEX_NODE) = N2
      ENDDO
      
      INDEX_NODE_CAL = 0.0D0
      DO 100 I =1,INDEX_NODE
          NODE = NMOORING_NODE(I)
          IF (NODE.EQ.0D0) GOTO 100
          OPT  = 0.0D0
          DO J = 1,INDEX_NODE
             IF (NODE.EQ.NMOORING_NODE(J)) THEN
             NMOORING_NODE(J) = 0.0D0
              IF (OPT.EQ.0D0) THEN
              INDEX_NODE_CAL                   = INDEX_NODE_CAL + 1
              NMOORING_NODE_RE(INDEX_NODE_CAL) = NODE
              OPT                              = 1D0
              ENDIF
             ENDIF
          ENDDO
100   CONTINUE
      
      ! CREATE A SEQ. NUMBER
      DO I = 1,INDEX_NODE_CAL
      CALL RELFILL('@XYZ',AX(I),1,NMOORING_NODE_RE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY(I),2,NMOORING_NODE_RE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ(I),3,NMOORING_NODE_RE(I),0)  !GETTING NODAL COORDINATE
      ENDDO
      
      MMAX_VALU_LINE = MAXLOC(AZ(1:INDEX_NODE_CAL),DIM=1)
      MMIN_VALU_LINE = MINLOC(AZ(1:INDEX_NODE_CAL),DIM=1)
      
      ! BEGIN POSITION LINE
      IF (OPT_CAL.EQ."BEG") THEN
      SELECTCASE (NGRAV)
      CASE (1)
      COOR_V = ABS(AX(MMIN_VALU_LINE)   - AMIN)
      COOR_H = SQRT((AY(MMIN_VALU_LINE) - AYY_MIN)**2+(AZ(MMIN_VALU_LINE) - AZZ_MIN)**2)
      CASE (2)
      COOR_V = ABS(AY(MMIN_VALU_LINE)   - AMIN)
      COOR_H = SQRT((AX(MMIN_VALU_LINE) - AXX_MIN)**2+(AZ(MMIN_VALU_LINE) - AZZ_MIN)**2)
      CASE (3)
      COOR_V = ABS(AZ(MMIN_VALU_LINE)   - AMIN)
      COOR_H = SQRT((AX(MMIN_VALU_LINE) - AXX_MIN)**2+(AY(MMIN_VALU_LINE) - AYY_MIN)**2)
      ENDSELECT
      
      ! END POSITION OF THE LINE
      ELSEIF (OPT_CAL.EQ."END")THEN
      SELECTCASE (NGRAV)
      CASE (1)
      COOR_V = ABS(AX(MMIN_VALU_LINE)   - AMAX)
      COOR_H = SQRT((AY(MMIN_VALU_LINE) - AYY_MAX)**2+(AZ(MMIN_VALU_LINE) - AZZ_MAX)**2)
      CASE (2)
      COOR_V = ABS(AY(MMIN_VALU_LINE)   - AMAX)
      COOR_H = SQRT((AX(MMIN_VALU_LINE) - AXX_MAX)**2+(AZ(MMIN_VALU_LINE) - AZZ_MAX)**2)
      CASE (3)
      COOR_V = ABS(AZ(MMIN_VALU_LINE)   - AMAX)
      COOR_H = SQRT((AX(MMIN_VALU_LINE) - AXX_MAX)**2+(AY(MMIN_VALU_LINE) - AYY_MAX)**2)
      ENDSELECT
      ENDIF
      
      IF (COOR_H.LE.ALB)THEN ! SEABED CONTRACT
      HA = FAIRLEADH - CB*W*(ALB - COOR_H)  
      VA = 0.0D0
      ELSEIF (COOR_H.GT.ALB) THEN ! SUSPENDED CABLE
      HA = FAIRLEADH - CB*W*(COOR_H_INPUT-COOR_H)  
      VA = FAIRLEADV - W*(COOR_V_INPUT-COOR_V)
      ENDIF
      
      END