      SUBROUTINE WIND_ASCE7_02    
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /WIND_LATERAL/ WIND_D(3,10000),WIND_CASE(500),IN_FLOOR
      
      ! ILOAD         = TOTAL NUMBER OF WIND STANDARD
      ! IOPT_WIND     = OPTION OF WIND STANDARD
      ! IOPT_ADD_WIND = ADDITONAL WIND OPTION
      ! ITO_ADD       = TOTAL WIND ADDITIONAL WIND LOAD
      ! ILCN          = LOAD CASE NUMBER
      READ (ITI,*) ! READ HEAD
      READ (ITI,*) ILOAD
      
      DO III = 1,ILOAD
      READ (ITI,*) IOPT_WIND,IOPT_ADD_WIND,ITO_ADD,ILCN
      SELECTCASE (IOPT_WIND)
      CASE (1)
      ! SIMPLIFIRED PROCEDURE
      ! V_WIND         = BASIC WIND SPEED
      ! E              = EXPOSURE FACTOR
      ! IMPORT_FACTOR  = IMPORTANCE FACTOR
      ! OPT_EXPOSURE   = EXPOSURE CATEGORY
      READ (ITI,*) V_WIND,IOPT_EXPOSURE,IMPORT_FACTOR,DIREC_COFF1,DIREC_COFF2
      WIND_CASE(III) = ILCN
      CALL IMPORT_FACTOR_ASCE (IMPORT_FACTOR,AIMPORT)
      CALL EXPOSURE_CATEGOTY_98 (IOPT_EXPOSURE,E) ! CALLING EXPOSURE FACTOR
      CALL DESIGN_WIND_PRESSURE (V_WIND,PB)       ! CALLING WIND PRESSURE 
      !CALL EXPOSURE_CATEGOTY_02 (OPT_EXPOSURE,HEIGHT,E)
      P = PB*AIMPORT*E
      
      CALL STORY_DATA_SIMPLIFIRED (P,DIREC_COFF1,DIREC_COFF2) ! CALCULATE WIND FORCE 
      
      ! ADDITIONAL WIND FORCE
      IF (IOPT_ADD_WIND.EQ.1) CALL WIND_ADD_FORCE (ITO_ADD)
      
      ! PRINTING CALCULATION SHEET
      
      
      CASE (2)
      ! ANALYTICAL PROCEDURE
      ! V_WIND        = BASIC WIND SPEED
      ! IOPT_EXPOSURE = OPTION OF EXPOSURE (B,C,D)
      ! IMPORT_FACTOR = IMPORTANCE FACTOR
      ! WIND_DIREC_FACTOR = WIND DIRECTION FACTOR
      ! AMEAN_HEIGHT  = MEAN ROOF HEIGHT
      ! AGX           = GUST EFFECT FACTOR
      ! AGY           = GUST EFFECT FACTOR
      ! IOPT_TOPO     = OPTION OF TOPOGRAHIC FACTOR
      ! IOPT_HILL     = OPTION OF HILL SHAPE (2D-RID, 2D-ESC, 3D-AXIS)   
      ! IOPT_BUILD    = OPTIOIN OF BUILDING LOCATION (UPWIND, DOWNWIND)
      ! AHILL_HEIGHT  = HILL HEIGHT
      ! AHILL_LENGTH  = HILL LENGTH
      ! CREAT_DISTAN  = CREST-NUILDING DISTANCE
      ! IOPT_FORCE_COEFF = OPTION OF LOAD EVALUATION
      ! FOCE_COFF     = FORCE COEFFICIENT
      ! DIREC_COFF1   = DIRECTION FACTOR
      ! DIREC_COFF2   = DIRECTION FACTOR
      ! ILCN          = LOADCASE NUMBER
      READ (ITI,*) V_WIND,IOPT_EXPOSURE,IMPORT_FACTOR,WIND_DIREC_FACTOR,AMEAN_HEIGHT,AGX,AGY,IOPT_TOPO,IOPT_HILL,IOPT_BUILD
     1             ,AHILL_HEIGHT,AHILL_LENGTH,CREAT_DISTAN,IOPT_FORCE_COEFF,FOCE_COFF,DIREC_COFF1,DIREC_COFF2
    
      WIND_CASE(III) = ILCN
      CALL IMPORT_FACTOR_ASCE (IMPORT_FACTOR,AIMPORT)
      IF (IOPT_TOPO.EQ.1) CALL TOPOGRAPHIC_FACTOR_KZT (AHILL_HEIGHT,AHILL_LENGTH,IOPT_TOPO,AKZT)
      IF (IOPT_TOPO.NE.1) AKZT = 1.0D0
      
      ! MODIFY INSIDE STORY_DATA_PROCEDURE 
      AKZ = 1.0D0
      AKH = 1.0D0
      
      ! VELOCITY PRESSURE
       !QZ = 0.00256D0*AKZ*AKZT*WIND_DIREC_FACTOR*(V_WIND**2D0)*AIMPORT ! LB/FT^2
        QZ = 0.613D0*AKZ*AKZT*WIND_DIREC_FACTOR*(V_WIND**2D0)*AIMPORT  ! N/M^2
        QH = 0.613D0*AKH*AKZT*WIND_DIREC_FACTOR*(V_WIND**2D0)*AIMPORT  ! N/M^2
      ! DESIGN WIND PRESSURE
        IF (IOPT_FORCE_COEFF.EQ.0) FOCE_COFF = 1.0D0
        CALL STORY_DATA_PROCEDURE (QZ,QH,AGX,AGY,AMEAN_HEIGHT,IOPT_EXPOSURE,FOCE_COFF,DIREC_COFF1,DIREC_COFF2,ILCN) ! CALCULATE WIND FORCE 
        
      ! ADDITIONAL WIND FORCE
        IF (IOPT_ADD_WIND.EQ.1) CALL WIND_ADD_FORCE (ITO_ADD,ILCN)
        
      ! PRINTING CALCULATION SHEET
        
      CASE (3) ! PROJECTED AREA METHOD (UBC1997)
      WIND_CASE(III) = ILCN    
      CALL UBC1997 (IOPT_WIND,IOPT_ADD_WIND)   
      
      CASE (4) ! NORMAL FORCE METHOD (UBC1997)
      WIND_CASE(III) = ILCN    
      CALL UBC1997 (IOPT_WIND,IOPT_ADD_WIND)    
      
      CASE (5) ! SIMPLIPIED PROCEDURE (NBC1995)
      WIND_CASE(III) = ILCN    
      CALL WIND_NBC1995 (IOPT_WIND,IOPT_ADD_WIND)
      
      CASE (6) ! DETAILED PROCEDURE (NBC1995)
      WIND_CASE(III) = ILCN    
      CALL WIND_NBC1995 (IOPT_WIND,IOPT_ADD_WIND)
      
      CASE (7) ! SIMPLIFIED PROCEDURE (EUROCODE1-1997)
      WIND_CASE(III) = ILCN    
      CALL WIND_EUROCODE_1997 (IOPT_WIND)    
      CASE (8) ! DETAILED PROCEDURE (EUROCODE1-1997)
      WIND_CASE(III) = ILCN    
      CALL WIND_EUROCODE_1997 (IOPT_WIND)
      
      CASE (9) ! STANDATD METHOD (BS6399-1997)
      WIND_CASE(III) = ILCN    
      CALL WIND_BS6399_1997 (IOPT_WIND,IOPT_ADD_WIND)
      CASE (10) ! DIRECTIONAL METHOD (BS6399-1997)    
      WIND_CASE(III) = ILCN    
      CALL WIND_BS6399_1997 (IOPT_WIND,IOPT_ADD_WIND)
      
      
      CASE (11) ! ANSI (1982)
      WIND_CASE(III) = ILCN    
      ENDSELECT
      
          
      
      ENDDO ! END OF LOAD 
      
      END
C	==================================================================
      SUBROUTINE WIND_ADD_FORCE (ITO_ADD,ILCN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /MGRAV/ NGRAV
      COMMON /WIND_LATERAL/ WIND_D(3,10000),WIND_CASE(500),IN_FLOOR
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      
      IN_FLOOR_ADD = IN_FLOOR-IFLOOR + 1
      DO I = 1,ITO_ADD
       READ (ITI,*) INDEX_STORY,F1,F2,ROT
        SELECTCASE(NGRAV)
        CASE(1)
        WIND_D(1,IN_FLOOR_ADD+I) = 0.0D0
        WIND_D(2,IN_FLOOR_ADD+I) = WIND_D(2,I) + F1
        WIND_D(3,IN_FLOOR_ADD+I) = WIND_D(3,I) + F2
        CASE(2)
        WIND_D(1,IN_FLOOR_ADD+I) = WIND_D(1,I) + F1
        WIND_D(2,IN_FLOOR_ADD+I) = 0.0D0
        WIND_D(3,IN_FLOOR_ADD+I) = WIND_D(3,I) + F2
        CASE(3)
        WIND_D(1,IN_FLOOR_ADD+I) = WIND_D(1,I) + F1
        WIND_D(2,IN_FLOOR_ADD+I) = WIND_D(2,I) + F2
        WIND_D(3,IN_FLOOR_ADD+I) = 0.0D0
        ENDSELECT
      ENDDO
      
      END
C	==================================================================
      SUBROUTINE IMPORT_FACTOR_ASCE (IMPORT_FACTOR,AIMPORT) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      SELECTCASE (IMPORT_FACTOR)
      CASE (1) 
      AIMPORT = 0.77D0
      CASE (2) 
      AIMPORT = 0.87D0
      CASE (3) 
      AIMPORT = 1.00D0
      CASE (4) 
      AIMPORT = 1.15D0
      ENDSELECT
      END
C	==================================================================
      SUBROUTINE DESIGN_WIND_PRESSURE (V_WIND,PB) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION PB_DATA(20)
      DATA PB_DATA/12.0D0,14.0D0,17.0D0,20.0D0,24.0D0,29.0D0,33.0D0,38.0D0,43.0D0,49.0D0
     1              ,85.0D0,90.0D0,100.0D0,110.0D0,120.0D0,130.0D0,140.0D0,150.0D0,160.0D0,170.0D0/
      
        DO 100 I = 1,10
            REAL_SPEED = PB_DATA(10+I)*0.44704D0 ! MILE/H > M/SEC
            IF (V_WIND.LE.PB_DATA(11)*0.44704D0)      THEN
            PB = 12.00D0
            EXIT
            ELSEIF (V_WIND.LE.REAL_SPEED) THEN
            IF(I.NE.1) DELTA_PRESSURE = PB_DATA(I)    -  PB_DATA(I-1)
            IF(I.NE.1) DELTA_SPEED    = PB_DATA(I+10) -  PB_DATA(I+10-1)
            ERR_CHECK = REAL_SPEED-V_WIND  
            IF (ERR_CHECK.NE.0.) PB = PB_DATA(I-1) + ((DELTA_PRESSURE/DELTA_SPEED)*(REAL_SPEED-V_WIND))
            IF (ERR_CHECK.EQ.0.) PB = PB_DATA(I)
            EXIT 
            ELSEIF (V_WIND.GT.PB_DATA(20)*0.44704D0)  THEN
            PB = 49.0D0
            EXIT
            ENDIF
100     CONTINUE
        
      PB = PB * 47.880258888889D0 ! UNIT CONVENTION
      END
C	==================================================================
      SUBROUTINE EXPOSURE_CATEGOTY_98 (IOPT_EXPOSURE,E) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      SELECTCASE (IOPT_WIND)
      CASE (1) ! B
      E = 1.00D0
      CASE (2) ! C
      E = 1.40D0
      CASE (3) ! D
      E = 1.66D0  
      ENDSELECT
      END
C	==================================================================
      SUBROUTINE EXPOSURE_CATEGOTY_02 (OPT_EXPOSURE,HEIGHT,E) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION B_DATA(20),C_DATA(20),D_DATA(20)
      DATA B_DATA/1.00D0,1.00D0,1.00D0,1.00D0,1.05D0,1.09D0,1.12D0,1.16D0,1.19D0,1.22D0
     1           ,15D0,20D0,25D0,30D0,35D0,40D0,45D0,50D0,55D0,60D0/
      DATA C_DATA/1.21D0,1.29D0,1.35D0,1.40D0,1.45D0,1.49D0,1.53D0,1.56D0,1.59D0,1.62D0
     1           ,15D0,20D0,25D0,30D0,35D0,40D0,45D0,50D0,55D0,60D0/
      DATA D_DATA/1.00D0,1.00D0,1.00D0,1.00D0,1.05D0,1.09D0,1.12D0,1.16D0,1.19D0,1.22D0
     1           ,15D0,20D0,25D0,30D0,35D0,40D0,45D0,50D0,55D0,60D0/
      
      END
C	==================================================================
      SUBROUTINE STORY_DATA_SIMPLIFIRED (PRESSURE,DIREC_COFF1,DIREC_COFF2) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /MGRAV/ NGRAV
      COMMON /WIND_LATERAL/ WIND_D(3,10000),WIND_CASE(500),IN_FLOOR
      DIMENSION AREA_OUT1(IFLOOR),AREA_OUT2(IFLOOR) 

      AREA1 = 0.
      AREA2 = 0.
      CALL FLOOR_AREA (AREA_OUT1,AREA_OUT2)
      
      DO I = 1,IFLOOR 
      SELECTCASE(NGRAV)
      CASE(1)
      IN_FLOOR           = IN_FLOOR + 1
      WIND_D(1,IN_FLOOR) = 0.0D0
      WIND_D(2,IN_FLOOR) = PRESSURE*AREA_OUT2(I)*DIREC_COFF1
      WIND_D(3,IN_FLOOR) = PRESSURE*AREA_OUT1(I)*DIREC_COFF2
      CASE(2)
      IN_FLOOR           = IN_FLOOR + 1
      WIND_D(1,IN_FLOOR) = PRESSURE*AREA_OUT2(I)*DIREC_COFF1
      WIND_D(2,IN_FLOOR) = 0.0D0
      WIND_D(3,IN_FLOOR) = PRESSURE*AREA_OUT1(I)*DIREC_COFF2
      CASE(3)
      IN_FLOOR           = IN_FLOOR + 1
      WIND_D(1,IN_FLOOR) = PRESSURE*AREA_OUT2(I)*DIREC_COFF1
      WIND_D(2,IN_FLOOR) = PRESSURE*AREA_OUT1(I)*DIREC_COFF2
      WIND_D(3,IN_FLOOR) = 0.0D0
      ENDSELECT
      ENDDO
      
      END
C	==================================================================
      SUBROUTINE STORY_DATA_PROCEDURE (QZ,QH,AGX,AGY,AMEAN_HEIGHT,IOPT_EXPOSURE,FOCE_COFF,DIREC_COFF1,DIREC_COFF2,ILCN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /MGRAV/ NGRAV
      COMMON /WIND_LATERAL/ WIND_D(3,10000),WIND_CASE(500),IN_FLOOR
      DIMENSION AREA_OUT1(IFLOOR),AREA_OUT2(IFLOOR) 
      DIMENSION ELEVATION(IFLOOR)
      DIMENSION FORCE1(IFLOOR),FORCE2(IFLOOR)
      
      AREA_OUT1 = 0.
      AREA_OUT2 = 0.
      CALL FLOOR_AREA_SPAN (AREA_OUT1,AREA_OUT2,ELEVATION)
      CALL L_B_LENGTH (ALENGTH_L,ALENGTH_B,ELEVATION)
      CALL WALL_PRESSURE_COEFF (ALENGTH_L,ALENGTH_B,CP_WIND,CP_LEE1,CP_LEE2)
      
      
      DO I = 1,IFLOOR
       IF (I.NE.IFLOOR-1) CALL VELOCITY_PRESSURE_KH_CALL (ELEVATION(I),IOPT_EXPOSURE,AKZ1)      ! WINDWARD
       IF (I.EQ.IFLOOR-1) CALL VELOCITY_PRESSURE_KH_CALL (ELEVATION(IFLOOR),IOPT_EXPOSURE,AKZ1) ! WINDWARD
       CALL VELOCITY_PRESSURE_KH_CALL (AMEAN_HEIGHT,IOPT_EXPOSURE,AKZ2)                         ! LEEWARD
      QZ_CAL = QZ*AKZ1 ! WINDWARD
      QH_CAL = QH*AKZ2 ! LEEWARD
      ! QH = 0.613D0*KH*KZT*KD*(V**2D0)*I  ! N/M^2
      PRESSURE1 = QZ_CAL*(AGX*CP_WIND) - QH_CAL*(AGX*CP_LEE1)
      PRESSURE2 = QZ_CAL*(AGY*CP_WIND) - QH_CAL*(AGY*CP_LEE2)
      
      FORCE1(I) = PRESSURE1*AREA_OUT2(I)*DIREC_COFF1*FOCE_COFF
      FORCE2(I) = PRESSURE2*AREA_OUT1(I)*DIREC_COFF2*FOCE_COFF
      ENDDO
      
      DO I = 1,IFLOOR
      SELECTCASE(NGRAV)
      CASE(1)
        IN_FLOOR           = IN_FLOOR + 1
        IF (I.EQ.1.OR.I.EQ.IFLOOR) THEN
        WIND_D(1,IN_FLOOR) = 0.0D0
        WIND_D(2,IN_FLOOR) = FORCE1(I)/2D0
        WIND_D(3,IN_FLOOR) = FORCE2(I)/2D0
        ELSEIF (I.NE.1) THEN
        WIND_D(1,IN_FLOOR) = 0.0D0
        WIND_D(2,IN_FLOOR) = (FORCE1(I) + FORCE1(I-1))/2D0
        WIND_D(3,IN_FLOOR) = (FORCE2(I) + FORCE2(I-1))/2D0
        ENDIF
      CASE(2)
        IN_FLOOR           = IN_FLOOR + 1
        IF (I.EQ.1.OR.I.EQ.IFLOOR) THEN
        WIND_D(1,IN_FLOOR) = FORCE1(I)/2D0
        WIND_D(2,IN_FLOOR) = 0.0D0
        WIND_D(3,IN_FLOOR) = FORCE2(I)/2D0
        ELSEIF (I.NE.1) THEN
        WIND_D(1,IN_FLOOR) = (FORCE1(I) + FORCE1(I-1))/2D0
        WIND_D(2,IN_FLOOR) = 0.0D0
        WIND_D(3,IN_FLOOR) = (FORCE2(I) + FORCE2(I-1))/2D0
        ENDIF
      CASE(3)
        IN_FLOOR           = IN_FLOOR + 1
        IF (I.EQ.1.OR.I.EQ.IFLOOR) THEN
        WIND_D(1,IN_FLOOR) = FORCE1(I)/2D0
        WIND_D(2,IN_FLOOR) = 0.0D0
        WIND_D(3,IN_FLOOR) = FORCE2(I)/2D0
        ELSEIF (I.NE.1) THEN
        WIND_D(1,IN_FLOOR) = (FORCE1(I) + FORCE1(I-1))/2D0
        WIND_D(2,IN_FLOOR) = (FORCE2(I) + FORCE2(I-1))/2D0
        WIND_D(3,IN_FLOOR) = 0.0D0
        ENDIF
      ENDSELECT
      ENDDO
      
      END
C	==================================================================
      SUBROUTINE WALL_PRESSURE_COEFF (RATIO_L,RATIO_B,QZ,QH1,QH2)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      QZ = 0.8D0
      
      IF (RATIO_L.LE.1.0D0)THEN
      QH1 = -0.5D0
      ELSEIF (RATIO_L.GT.1.0D0.AND.RATIO_L.LE.2.0)THEN
      DELTA_RATIO = 2.0D0 - 1.0D0
      DELTA_CP    = -0.3D0
      QH1         = -0.5D0 + (DELTA_CP/DELTA_RATIO)*(RATIO_L-1.0D0)
         IF (RATIO_L.EQ.2D0) QH1 = -0.3D0
         IF (RATIO_L.EQ.1D0) QH1 = -0.5
      ELSEIF (RATIO_L.GE.4.0D0)THEN
      QH1         = -0.2D0
      ENDIF
      
      IF (RATIO_B.LE.1.0D0)THEN
      QH2 = -0.5D0
      ELSEIF (RATIO_B.GT.1.0D0.AND.RATIO_B.LE.2.0)THEN
      DELTA_RATIO = 2.0D0 - 1.0D0
      DELTA_CP    = -0.3D0
      QH2         = -0.5D0 + (DELTA_CP/DELTA_RATIO)*(RATIO_B-1.0D0)
         IF (RATIO_B.EQ.2D0) QH2 = -0.3D0
         IF (RATIO_B.EQ.1D0) QH2 = -0.5
      ELSEIF (RATIO_B.GE.4.0D0)THEN
      QH2         = -0.2D0
      ENDIF
      END
C	==================================================================
      SUBROUTINE L_B_LENGTH (RATIO_L,RATIO_B,ELEVATION_LENGTH_FLOOR)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /DIAPH_DATA/ NSLAVE_NODE(500,100),NSLAVE(500)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /MGRAV/ NGRAV
      DIMENSION AREA_OUT1(IFLOOR),AREA_OUT2(IFLOOR)
      DIMENSION FLOOR_MID_ELEVATION(IFLOOR)
      ALLOCATABLE NSLAVE_NODE_FLOOR(:),A1(:),A2(:),ELEVATION(:)
      ALLOCATABLE A1_LENGTH_FLOOR(:),A2_LENGTH_FLOOR(:)
      DIMENSION ELEVATION_LENGTH_FLOOR(IFLOOR)
      
      ! FLOOR PROJECTED AREA
      ALLOCATE (A1_LENGTH_FLOOR(IFLOOR),A2_LENGTH_FLOOR(IFLOOR))
      DO I = 1,IFLOOR
         ALLOCATE (NSLAVE_NODE_FLOOR(NSLAVE(I)),A1(NSLAVE(I)),A2(NSLAVE(I)),ELEVATION(NSLAVE(I)))
         A1 = 0.
         A2 = 0.
         ELEVATION = 0.
         NSLAVE_NODE_FLOOR = 0.
         NSLAVE_NODE_FLOOR(1:NSLAVE(I)) = NSLAVE_NODE(I,1:NSLAVE(I))
         
         DO J = 1,NSLAVE(I)
         SELECTCASE (NGRAV)
         CASE (1) ! X-DIRECTION
         CALL RELFILL('@XYZ',ELEVATION(I),1,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A1(J),2,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A2(J),3,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
         CASE (2) ! Y-DIRECTION  
         CALL RELFILL('@XYZ',A1(J),1,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',ELEVATION(I),2,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A2(J),3,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
         CASE (3) ! Z-DIRECTION
         CALL RELFILL('@XYZ',A1(J),1,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A2(J),2,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',ELEVATION(J),3,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
         ENDSELECT 
         ENDDO
         
         ! FLOOR PROJECTED AREA
         COOR_MAX_A1_FLOOR        = MAXVAL(A1)
         COOR_MAX_A2_FLOOR        = MAXVAL(A2)
         COOR_MAX_ELEVATION_FLOOR = MAXVAL(ELEVATION)
         
         COOR_MIN_A1_FLOOR        = MAXVAL(A1)
         COOR_MIN_A2_FLOOR        = MAXVAL(A2)
         COOR_MIN_ELEVATION_FLOOR = MAXVAL(ELEVATION) 
         
         A1_LENGTH_FLOOR(I) = COOR_MAX_A1_FLOOR - COOR_MIN_A1_FLOOR
         A2_LENGTH_FLOOR(I) = COOR_MAX_A2_FLOOR - COOR_MIN_A2_FLOOR
         ELEVATION_LENGTH_FLOOR(I) = COOR_MAX_ELEVATION_FLOOR - COOR_MIN_ELEVATION_FLOOR
         IF (COOR_MAX_A1_FLOOR.EQ.COOR_MIN_A1_FLOOR) A1_LENGTH_FLOOR(I) = COOR_MAX_A1_FLOOR
         IF (COOR_MAX_A2_FLOOR.EQ.COOR_MIN_A2_FLOOR) A2_LENGTH_FLOOR(I) = COOR_MAX_A2_FLOOR
         IF (COOR_MAX_ELEVATION_FLOOR.EQ.COOR_MIN_ELEVATION_FLOOR) ELEVATION_LENGTH_FLOOR(I) = COOR_MAX_ELEVATION_FLOOR
      
         
         DEALLOCATE (NSLAVE_NODE_FLOOR,A1,A2,ELEVATION)
      ENDDO
      
      ! CALCULATE PROJECTED AREA EACH FLOOR
      AREA_OUT1 = 0.
      AREA_OUT2 = 0.
      DO I =1,IFLOOR
         IF (I.EQ.1)         THEN ! FIRST FLOOR
         ELEVATION_FLOOR        = ELEVATION_LENGTH_FLOOR(I+1) - ELEVATION_LENGTH_FLOOR(I)
         FLOOR_MID_ELEVATION(I) = (ELEVATION_FLOOR)/2D0
         AREA_OUT1(I)           = A1_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         AREA_OUT2(I)           = A2_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         ELSEIF (I.EQ.IFLOOR)THEN ! LAST FLOOR
         ELEVATION_FLOOR        = ELEVATION_LENGTH_FLOOR(I) - ELEVATION_LENGTH_FLOOR(I-1)
         FLOOR_MID_ELEVATION(I) = (ELEVATION_FLOOR)/2D0
         AREA_OUT1(I)           = A1_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         AREA_OUT2(I)           = A2_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         ELSEIF (I.NE.1)     THEN
         ELEVATION_FLOOR1        = ELEVATION_LENGTH_FLOOR(I+1) - ELEVATION_LENGTH_FLOOR(I)    
         ELEVATION_FLOOR2        = ELEVATION_LENGTH_FLOOR(I) - ELEVATION_LENGTH_FLOOR(I-1)  
         FLOOR_MID_ELEVATION(I) = (ELEVATION_FLOOR1+ELEVATION_FLOOR2)/2D0
         AREA_OUT1(I)           = A1_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         AREA_OUT2(I)           = A2_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         ENDIF
      ENDDO
      
      SELECTCASE (NGRAV)
      CASE (1) ! X-DIRECTION
      RATIO_L = MAXVAL(A1_LENGTH_FLOOR)/MAXVAL(A2_LENGTH_FLOOR) ! Y/Z
      RATIO_B = MAXVAL(A2_LENGTH_FLOOR)/MAXVAL(A1_LENGTH_FLOOR) ! Z/Y
      CASE (2) ! Y-DIRECTION
      RATIO_L = MAXVAL(A1_LENGTH_FLOOR)/MAXVAL(A2_LENGTH_FLOOR) ! X/Z
      RATIO_B = MAXVAL(A2_LENGTH_FLOOR)/MAXVAL(A1_LENGTH_FLOOR) ! Z/X  
      CASE (3) ! Z-DIRECTION
      RATIO_L = MAXVAL(A1_LENGTH_FLOOR)/MAXVAL(A2_LENGTH_FLOOR) ! X/Y
      RATIO_B = MAXVAL(A2_LENGTH_FLOOR)/MAXVAL(A1_LENGTH_FLOOR) ! Y/X
      ENDSELECT
      
      DEALLOCATE (A1_LENGTH_FLOOR,A2_LENGTH_FLOOR)
      END
C	==================================================================
      SUBROUTINE FLOOR_AREA (AREA_OUT1,AREA_OUT2)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /DIAPH_DATA/ NSLAVE_NODE(500,100),NSLAVE(500)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /MGRAV/ NGRAV
      DIMENSION AREA_OUT1(IFLOOR),AREA_OUT2(IFLOOR)
      DIMENSION FLOOR_MID_ELEVATION(IFLOOR)
      ALLOCATABLE NSLAVE_NODE_FLOOR(:),A1(:),A2(:),ELEVATION(:)
      ALLOCATABLE A1_LENGTH_FLOOR(:),A2_LENGTH_FLOOR(:),ELEVATION_LENGTH_FLOOR(:)
      
      ! FLOOR PROJECTED AREA
      ALLOCATE (A1_LENGTH_FLOOR(IFLOOR),A2_LENGTH_FLOOR(IFLOOR),ELEVATION_LENGTH_FLOOR(IFLOOR))
      DO I = 1,IFLOOR
         ALLOCATE (NSLAVE_NODE_FLOOR(NSLAVE(I)),A1(NSLAVE(I)),A2(NSLAVE(I)),ELEVATION(NSLAVE(I)))
         A1 = 0.
         A2 = 0.
         ELEVATION = 0.
         NSLAVE_NODE_FLOOR = 0.
         NSLAVE_NODE_FLOOR(1:NSLAVE(I)) = NSLAVE_NODE(I,1:NSLAVE(I))
         
         DO J = 1,NSLAVE(I)
         SELECTCASE (NGRAV)
         CASE (1) ! X-DIRECTION
         CALL RELFILL('@XYZ',ELEVATION(I),1,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A1(J),2,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A2(J),3,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
         CASE (2) ! Y-DIRECTION  
         CALL RELFILL('@XYZ',A1(J),1,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',ELEVATION(I),2,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A2(J),3,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
         CASE (3) ! Z-DIRECTION
         CALL RELFILL('@XYZ',A1(J),1,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',A2(J),2,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	   CALL RELFILL('@XYZ',ELEVATION(J),3,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
         ENDSELECT 
         ENDDO
         
         ! FLOOR PROJECTED AREA
         COOR_MAX_A1_FLOOR        = MAXVAL(A1)
         COOR_MAX_A2_FLOOR        = MAXVAL(A2)
         COOR_MAX_ELEVATION_FLOOR = MAXVAL(ELEVATION)
         
         COOR_MIN_A1_FLOOR        = MAXVAL(A1)
         COOR_MIN_A2_FLOOR        = MAXVAL(A2)
         COOR_MIN_ELEVATION_FLOOR = MAXVAL(ELEVATION) 
         
         A1_LENGTH_FLOOR(I) = COOR_MAX_A1_FLOOR - COOR_MIN_A1_FLOOR
         A2_LENGTH_FLOOR(I) = COOR_MAX_A2_FLOOR - COOR_MIN_A2_FLOOR
         ELEVATION_LENGTH_FLOOR(I) = COOR_MAX_ELEVATION_FLOOR - COOR_MIN_ELEVATION_FLOOR
         IF (COOR_MAX_A1_FLOOR.EQ.COOR_MIN_A1_FLOOR) A1_LENGTH_FLOOR(I) = COOR_MAX_A1_FLOOR
         IF (COOR_MAX_A2_FLOOR.EQ.COOR_MIN_A2_FLOOR) A2_LENGTH_FLOOR(I) = COOR_MAX_A2_FLOOR
         IF (COOR_MAX_ELEVATION_FLOOR.EQ.COOR_MIN_ELEVATION_FLOOR) ELEVATION_LENGTH_FLOOR(I) = COOR_MAX_ELEVATION_FLOOR
      
         
         DEALLOCATE (NSLAVE_NODE_FLOOR,A1,A2,ELEVATION)
      ENDDO
      
      ! CALCULATE PROJECTED AREA EACH FLOOR
      AREA_OUT1 = 0.
      AREA_OUT2 = 0.
      DO I =1,IFLOOR
         IF (I.EQ.1)         THEN ! FIRST FLOOR
         ELEVATION_FLOOR        = ELEVATION_LENGTH_FLOOR(I+1) - ELEVATION_LENGTH_FLOOR(I)
         FLOOR_MID_ELEVATION(I) = (ELEVATION_FLOOR)/2D0
         AREA_OUT1(I)           = A1_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         AREA_OUT2(I)           = A2_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         ELSEIF (I.EQ.IFLOOR)THEN ! LAST FLOOR
         ELEVATION_FLOOR        = ELEVATION_LENGTH_FLOOR(I) - ELEVATION_LENGTH_FLOOR(I-1)
         FLOOR_MID_ELEVATION(I) = (ELEVATION_FLOOR)/2D0
         AREA_OUT1(I)           = A1_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         AREA_OUT2(I)           = A2_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         ELSEIF (I.NE.1)     THEN
         ELEVATION_FLOOR1        = ELEVATION_LENGTH_FLOOR(I+1) - ELEVATION_LENGTH_FLOOR(I)    
         ELEVATION_FLOOR2        = ELEVATION_LENGTH_FLOOR(I) - ELEVATION_LENGTH_FLOOR(I-1)  
         FLOOR_MID_ELEVATION(I) = (ELEVATION_FLOOR1+ELEVATION_FLOOR2)/2D0
         AREA_OUT1(I)           = A1_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         AREA_OUT2(I)           = A2_LENGTH_FLOOR(I)*FLOOR_MID_ELEVATION(I)
         ENDIF
      ENDDO
      
      DEALLOCATE (A1_LENGTH_FLOOR,A2_LENGTH_FLOOR,ELEVATION_LENGTH_FLOOR)
      
      END
C	==================================================================
      SUBROUTINE VELOCITY_PRESSURE_KH_CALL (ELEVATION_IN1,IOPT_EXPOSURE,AKZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      ! THE KZ IS CALCULATED IN UNIT FT
      ELEVATION = ELEVATION_IN1/0.3048D0 ! M > FT
      
      SELECTCASE (IOPT_EXPOSURE)
      CASE (1)
      ALPHA = 7.0D0
      ZG    = 1200D0
      AL    = 320D0
      CASE (2)
      ALPHA = 9.5D0
      ZG    = 900D0
      AL    = 500D0
      CASE (3)
      ALPHA = 11.5D0
      ZG    = 70D0
      AL    = 650D0
      ENDSELECT
      ALENGTH_MAXIMUM = 15D0
      
      IF (ALENGTH_MAXIMUM.LE.ELEVATION.AND.ALENGTH_MAXIMUM.LE.ZG)THEN
         AKZ = 2.01D0*(ELEVATION/ZG)**(2D0/ALPHA) 
      ELSEIF (ELEVATION.LT.ALENGTH_MAXIMUM)THEN
         AKZ = 2.01D0*((15D0/ZG)**(2D0/ALPHA)) 
      ENDIF
      
      
      
      END
C	==================================================================
      SUBROUTINE VELOCITY_PRESSURE_KH_TABLE6_3 (ELEVATION,IOPT_EXPOSURE,AKZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION BKZ_DATA(44),CKZ_DATA(44),DKZ_DATA(44)
      DIMENSION DATA_AKZ(44)
      DATA BKZ_DATA/0.57D0,0.62D0,0.66D0,0.70D0,0.760D0,0.81D0,0.85D0,0.89D0,0.93D0,0.96D0
     1             ,0.99D0,1.04D0,1.09D0,1.13D0,1.17D0,1.20D0,1.28D0,1.35D0,1.41D0,1.47D0
     1             ,1.52D0,1.56D0
     1             ,4.6D0,6.1D0,7.6D0,9.1D0,12.2D0,15.2D0,18.0D0,21.3D0,24.4D0,27.4D0
     1             ,30.5D0,36.6D0,42.7D0,48.8D0,54.9D0,61.0D0,76.2D0,91.4D0,106.7D0,121.9D0
     1             ,137.2D0,152.4D0/
      DATA CKZ_DATA/0.85D0,0.90D0,0.94D0,0.98D0,1.04D0,1.09D0,1.13D0,1.17D0,1.21D0,1.24D0
     1             ,1.26D0,1.31D0,1.36D0,1.39D0,1.43D0,1.46D0,1.53D0,1.59D0,1.64D0,1.69D0
     1             ,1.73D0,1.77D0
     1             ,4.6D0,6.1D0,7.6D0,9.1D0,12.2D0,15.2D0,18.0D0,21.3D0,24.4D0,27.4D0
     1             ,30.5D0,36.6D0,42.7D0,48.8D0,54.9D0,61.0D0,76.2D0,91.4D0,106.7D0,121.9D0
     1             ,137.2D0,152.4D0/
      DATA DKZ_DATA/1.03D0,1.08D0,1.12D0,1.16D0,1.22D0,1.27D0,1.31D0,1.34D0,1.38D0,1.40D0
     1             ,1.43D0,1.48D0,1.52D0,1.55D0,1.58D0,1.61D0,1.68D0,1.73D0,1.78D0,1.82D0
     1             ,1.86D0,1.89D0
     1             ,4.6D0,6.1D0,7.6D0,9.1D0,12.2D0,15.2D0,18.0D0,21.3D0,24.4D0,27.4D0
     1             ,30.5D0,36.6D0,42.7D0,48.8D0,54.9D0,61.0D0,76.2D0,91.4D0,106.7D0,121.9D0
     1             ,137.2D0,152.4D0/
      
      SELECTCASE (IOPT_EXPOSURE)
      CASE (1)
      DATA_AKZ(1:44) =  BKZ_DATA(1:44)   
      CASE (2)
      DATA_AKZ(1:44) =  CKZ_DATA(1:44)       
      CASE (3)
      DATA_AKZ(1:44) =  DKZ_DATA(1:44)       
      ENDSELECT
      
        DO 100 I = 1,22
            IF (ELEVATION.LE.DATA_AKZ(23))      THEN
            SELECTCASE (IOPT_EXPOSURE)
            CASE (1)
            AKZ = 0.57D0 
            CASE (2)
            AKZ = 0.85D0    
            CASE (3)
            AKZ = 1.03D0    
            ENDSELECT
            EXIT
            ELSEIF (ELEVATION.LE.DATA_AKZ(22+I)) THEN
            IF(I.NE.1) DELTA_PRESSURE  = DATA_AKZ(I)    -  DATA_AKZ(I-1)
            IF(I.NE.1) DELTA_ELEVATION = DATA_AKZ(I+22) -  DATA_AKZ(I+22-1)
            ERR_CHECK = ELEVATION-DATA_AKZ(I+22-1)
            IF (ERR_CHECK.NE.0.) AKZ = DATA_AKZ(I-1) + ((DELTA_PRESSURE/DELTA_ELEVATION)*(ELEVATION-DATA_AKZ(I+22-1)))
            IF (ERR_CHECK.EQ.0.) AKZ = DATA_AKZ(I)
            EXIT 
            ELSEIF (ELEVATION.GT.DATA_AKZ(44))  THEN
            SELECTCASE (IOPT_EXPOSURE)
            CASE (1)
            AKZ = 1.56D0 
            CASE (2)
            AKZ = 1.77D0    
            CASE (3)
            AKZ = 1.89D0    
            ENDSELECT
            EXIT
            ENDIF
100     CONTINUE

      END
C	==================================================================
      SUBROUTINE TOPOGRAPHIC_FACTOR_KZT (HILL_HEIGHT,HILL_LENGTH,IOPT_TOPOGRAPHIC,AKZT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION K1_RID_DATA(14),K1_ESC_DATA(14),K1_AXI_DATA(14)  
      DIMENSION K2_RID_DATA(18),K2_ALL_DATA(18)
      DIMENSION K3_RID_DATA(26),K3_ESC_DATA(26),K3_AXI_DATA(26)  
      DIMENSION DATA_AK1(14),DATA_AK2(18),DATA_AK3(26)
      DATA K1_RID_DATA/0.29D0,0.36D0,0.43D0,0.51D0,0.58D0,0.65D0,0.72D0
     1                 ,0.20D0,0.25D0,0.30D0,0.35D0,0.40D0,0.45D0,0.50D0/
      DATA K1_ESC_DATA/0.17D0,0.21D0,0.26D0,0.30D0,0.34D0,0.38D0,0.43D0
     1                 ,0.20D0,0.25D0,0.30D0,0.35D0,0.40D0,0.45D0,0.50D0/
      DATA K1_AXI_DATA/0.21D0,0.26D0,0.32D0,0.37D0,0.42D0,0.47D0,0.53D0
     1                 ,0.20D0,0.25D0,0.30D0,0.35D0,0.40D0,0.45D0,0.50D0/
      
      DATA K2_RID_DATA/1.00D0,0.88D0,0.75D0,0.63D0,0.50D0,0.38D0,0.25D0,0.13D0,0.00D0
     1                 ,0.00D0,0.5D0,1.00D0,1.50D0,2.00D0,2.50D0,3.00D0,3.50D0,4.00D0/
      DATA K2_ALL_DATA/1.00D0,0.67D0,0.33D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0
     1                 ,0.00D0,0.5D0,1.00D0,1.50D0,2.00D0,2.50D0,3.00D0,3.50D0,4.00D0/
      
      DATA K3_RID_DATA/1.00D0,0.74D0,0.55D0,0.41D0,0.30D0,0.22D0,0.17D0,0.12D0,0.09D0,0.07D0,0.05D0,0.01D0,0.00D0
     1                 ,0.00D0,0.10D0,0.20D0,0.30D0,0.40D0,0.50D0,0.60D0,0.70D0,0.80D0,0.90D0,1.00D0,1.50D0,2.00D0/
      DATA K3_ESC_DATA/1.00D0,0.78D0,0.61D0,0.47D0,0.37D0,0.29D0,0.22D0,0.17D0,0.14D0,0.11D0,0.08D0,0.02D0,0.00D0
     1                 ,0.00D0,0.10D0,0.20D0,0.30D0,0.40D0,0.50D0,0.60D0,0.70D0,0.80D0,0.90D0,1.00D0,1.50D0,2.00D0/
      DATA K3_AXI_DATA/1.00D0,0.67D0,0.5D0,0.30D0,0.20D0,0.14D0,0.09D0,0.06D0,0.04D0,0.03D0,0.02D0,0.00D0,0.00D0
     1                 ,0.00D0,0.10D0,0.20D0,0.30D0,0.40D0,0.50D0,0.60D0,0.70D0,0.80D0,0.90D0,1.00D0,1.50D0,2.00D0/
      ! DEFINED K1
      RATIO = HILL_HEIGHT/HILL_LENGTH
      SELECTCASE(IOPT_TOPOGRAPHIC)
      CASE (1) ! 2-D RIDGE
      DATA_AK1(1:14) =  K1_RID_DATA(1:14)       
      CASE (2) ! 2-D ESCARP
      DATA_AK1(1:14) =  K1_ESC_DATA(1:14)   
      CASE (3) ! 3-D AXISYM HILL
      DATA_AK1(1:14) =  K1_AXI_DATA(1:14)  
      ENDSELECT
          
      DO 100 I = 1,7
            IF (RATIO.LE.DATA_AK1(8))      THEN
              SELECTCASE(IOPT_TOPOGRAPHIC)
              CASE (1) ! 2-D RIDGE
              AK1 = 0.29D0
              CASE (2) ! 2-D ESCARP
              AK1 = 0.17D0
              CASE (3) ! 3-D AXISYM HILL
              AK1 = 0.21D0
              ENDSELECT
              EXIT
            ELSEIF (RATIO.LE.DATA_AK1(7+I)) THEN
            IF(I.NE.1) DELTA_PRESSURE  = DATA_AK1(I)   -  DATA_AK1(I-1)
            IF(I.NE.1) DELTA_RATIO     = DATA_AK1(I+7) -  DATA_AK1(I+7-1)
            ERR_CHECK = RATIO-DATA_AK1(I+7-1)
            IF (ERR_CHECK.NE.0) AK1 = DATA_AK1(I-1) + ((DELTA_PRESSURE/DELTA_RATIO)*(RATIO-DATA_AK1(I+7-1)))
            IF (ERR_CHECK.EQ.0) AK1 = DATA_AK1(I)
            EXIT 
            ELSEIF (RATIO.GT.DATA_AK1(14))  THEN
              SELECTCASE(IOPT_TOPOGRAPHIC)
              CASE (1) ! 2-D RIDGE
              AK1 = 0.72D0
              CASE (2) ! 2-D ESCARP
              AK1 = 0.43D0
              CASE (3) ! 3-D AXISYM HILL
              AK1 = 0.53D0
              ENDSELECT
            EXIT
            ENDIF
100   CONTINUE
      
      ! DEFINED K2
      RATIO = HILL_HEIGHT/HILL_LENGTH
      SELECTCASE(IOPT_TOPOGRAPHIC)
      CASE (1) ! 2-D RIDGE
      DATA_AK2(1:18) =  K2_ALL_DATA(1:18)       
      CASE (2) ! 2-D ESCARP
      DATA_AK2(1:18) =  K2_RID_DATA(1:18)   
      CASE (3) ! 3-D AXISYM HILL
      DATA_AK2(1:18) =  K2_ALL_DATA(1:18)  
      ENDSELECT
      
      DO 200 I = 1,9
            IF (RATIO.LE.DATA_AK2(10))      THEN
              SELECTCASE(IOPT_TOPOGRAPHIC)
              CASE (1) ! 2-D RIDGE
              AK2 = 1.00D0
              CASE (2) ! 2-D ESCARP
              AK2 = 1.00D0
              CASE (3) ! 3-D AXISYM HILL
              AK2 = 1.00D0
              ENDSELECT
              EXIT
            ELSEIF (RATIO.LE.DATA_AK2(9+I)) THEN
            IF(I.NE.1) DELTA_PRESSURE  = DATA_AK2(I)   -  DATA_AK2(I-1)
            IF(I.NE.1) DELTA_RATIO     = DATA_AK2(I+9) -  DATA_AK2(I+9-1)
            ERR_CHECK = RATIO-DATA_AK2(I+9-1)
            IF (ERR_CHECK.NE.0) AK2 = DATA_AK2(I-1) + ((DELTA_PRESSURE/DELTA_RATIO)*(RATIO-DATA_AK2(I+9-1)))
            IF (ERR_CHECK.EQ.0) AK2 = DATA_AK2(IS)
            EXIT 
            ELSEIF (RATIO.GT.DATA_AK2(18))  THEN
              SELECTCASE(IOPT_TOPOGRAPHIC)
              CASE (1) ! 2-D RIDGE
              AK2 = 0.00D0
              CASE (2) ! 2-D ESCARP
              AK2 = 0.00D0
              CASE (3) ! 3-D AXISYM HILL
              AK2 = 0.00D0
              ENDSELECT
            EXIT
            ENDIF
200   CONTINUE
      
      
      ! DEFINED K3
      RATIO = HILL_HEIGHT/HILL_LENGTH
      SELECTCASE(IOPT_TOPOGRAPHIC)
      CASE (1) ! 2-D RIDGE
      DATA_AK3(1:26) =  K3_RID_DATA(1:26)       
      CASE (2) ! 2-D ESCARP
      DATA_AK3(1:26) =  K3_ESC_DATA(1:26)   
      CASE (3) ! 3-D AXISYM HILL
      DATA_AK3(1:26) =  K3_AXI_DATA(1:26)  
      ENDSELECT
      
      DO 300 I = 1,13
            IF (RATIO.LE.DATA_AK3(14))      THEN
              SELECTCASE(IOPT_TOPOGRAPHIC)
              CASE (1) ! 2-D RIDGE
              AK3 = 1.00D0
              CASE (2) ! 2-D ESCARP
              AK3 = 1.00D0
              CASE (3) ! 3-D AXISYM HILL
              AK3 = 1.00D0
              ENDSELECT
              EXIT
            ELSEIF (RATIO.LE.DATA_AK3(13+I)) THEN
            IF(I.NE.1) DELTA_PRESSURE  = DATA_AK3(I)   -  DATA_AK3(I-1)
            IF(I.NE.1) DELTA_RATIO     = DATA_AK3(I+13) -  DATA_AK3(I+13-1)
            ERR_CHECK = RATIO-DATA_AK3(I+13-1)
            IF (ERR_CHECK.NE.0) AK3 = DATA_AK3(I-1) + ((DELTA_PRESSURE/DELTA_RATIO)*(RATIO-DATA_AK3(I+13-1)))
            IF (ERR_CHECK.EQ.0) AK3 = DATA_AK3(I)
            EXIT 
            ELSEIF (RATIO.GT.DATA_AK3(26))  THEN
              SELECTCASE(IOPT_TOPOGRAPHIC)
              CASE (1) ! 2-D RIDGE
              AK3 = 0.00D0
              CASE (2) ! 2-D ESCARP
              AK3 = 0.00D0
              CASE (3) ! 3-D AXISYM HILL
              AK3 = 0.00D0
              ENDSELECT
            EXIT
            ENDIF
300   CONTINUE
      
      AKZT = (1D0+AK1*AK2*AK3)**2D0
      END
      