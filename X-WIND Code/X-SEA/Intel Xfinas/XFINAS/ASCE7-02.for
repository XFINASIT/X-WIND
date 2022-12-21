      SUBROUTINE ASCE7_02    
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      
      READ (ITI,*) ! READ HEAD
      READ (ITI,*) IOPT_WIND,IOPT_ADD_WIND
      SELECTCASE (IOPT_WIND)
      CASE (1)
      ! SIMPLIFIRED PROCEDURE
      ! V_WIND         = BASIC WIND SPEED
      ! E              = EXPOSURE FACTOR
      ! IMPORT_FACTOR  = IMPORTANCE FACTOR
      ! OPT_EXPOSURE   = EXPOSURE CATEGORY
      
      READ (ITI,*) V_WIND,IOPT_EXPOSURE,IMPORT_FACTOR,E,SFX,SFY,ILCN
      CALL IMPORT_FACTOR_ASCE (IMPORT_FACTOR,AIMPORT)
      CALL EXPOSURE_CATEGOTY_98 (IOPT_EXPOSURE,E) ! CALLING EXPOSURE FACTOR
      CALL DESIGN_WIND_PRESSURE (V_WIND,PB)       ! CALLING WIND PRESSURE 
      !CALL EXPOSURE_CATEGOTY_02 (OPT_EXPOSURE,HEIGHT,E)
      P = PB*AIMPORT*E
      
      CALL STORY_DATA_SIMPLIFIRED (P) ! CALCULATE WIND FORCE 
      
      
      CASE (2)
      ! ANALYTICAL PROCEDURE
      READ (ITI,*) V_WIND,IOPT_EXPOSURE,AIMPORT,WIND_DIREC_FACTOR,AMEAN_HEIGHT,AGX,AGY,IOPT_TOPO,IOPT_HILL,IOPT_BUILD
     1             AHILL_HEIGHT,AHILL_LENGTH,CREAT_DISTAN,FOCE_COFF
      
      CALL IMPORT_FACTOR_ASCE (IMPORT_FACTOR,AIMPORT)
      CALL VELOCITY_PRESSURE_KH (ELEVATION,OPT_EXPOSURE,AKZ)
      CALL TOPOGRAPHIC_FACTOR_KZT (HILL_HEIGHT,HILL_LENGTH,OPT_TOPOGRAPHIC,AKZT)
      
      ! VELOCITY PRESSURE
       !QZ = 0.00256D0*AKZ*AKZT*WIND_DIREC_FACTOR*(V_WIND**2D0)*AIMPORT ! LB/FT^2
        QZ = 0.0613D0*AKZ*AKZT*WIND_DIREC_FACTOR*(V_WIND**2D0)*AIMPORT  ! N/M^2
        QH = 0.0613D0*AKZ*AKZT*WIND_DIREC_FACTOR*(V_WIND**2D0)*AIMPORT  ! N/M^2
      ! DESIGN WIND PRESSURE
        CALL STORY_DATA_PROCEDURE (QZ,QH,AGX,AGY) ! CALCULATE WIND FORCE 
          
      ENDSELECT
      
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
            !IF(I.EQ.1) DELTA_PRESSURE = PB_DATA(I+1)    -  0.0D0
            PB = PB_DATA(I) + ((DELTA_PRESSURE/DELTA_SPEED)*(REAL_SPEED-V_WIND))
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
      SUBROUTINE STORY_DATA_SIMPLIFIRED (PRESSURE) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /MGRAV/ NGRAV
      COMMON /WIND_LATERAL/ WIND_D(3,500)
      DIMENSION AREA_OUT1(IFLOOR),AREA_OUT2(IFLOOR) ! 500 FLOOR

      AREA1 = 0.
      AREA2 = 0.
      CALL FLOOR_AREA (AREA_OUT1,AREA_OUT2)
      
      DO I = 1,IFLOOR 
      SELECTCASE(NGRAV)
      CASE(1)
      WIND_D(1,I) = 0.0D0
      WIND_D(2,I) = PRESSURE*AREA_OUT1(I)
      WIND_D(3,I) = PRESSURE*AREA_OUT2(I)
      CASE(2)
      WIND_D(1,I) = PRESSURE*AREA_OUT1(I)
      WIND_D(2,I) = 0.0D0
      WIND_D(3,I) = PRESSURE*AREA_OUT2(I)
      CASE(3)
      WIND_D(1,I) = PRESSURE*AREA_OUT1(I)
      WIND_D(2,I) = PRESSURE*AREA_OUT2(I)
      WIND_D(3,I) = 0.0D0
      ENDSELECT
      ENDDO
      
      END
C	==================================================================
      SUBROUTINE STORY_DATA_PROCEDURE (QZ,QH,AGX,AGY) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /MGRAV/ NGRAV
      COMMON /WIND_LATERAL/ WIND_D(3,500)
      DIMENSION AREA_OUT1(IFLOOR),AREA_OUT2(IFLOOR) ! 500 FLOOR

      AREA1 = 0.
      AREA2 = 0.
      CALL FLOOR_AREA (AREA_OUT1,AREA_OUT2)
      CALL L_B_LENGTH (ALENGTH_L,ALENGTH_B)
      CALL WALL_PRESSURE_COEFF (RATIO_L,RATIO_B,CP1,CP2)

      PRESSURE = QZ*(AGX*CP1) - QH*(AGY*CP2)
      
      DO I = 1,IFLOOR 
      SELECTCASE(NGRAV)
      CASE(1)
      WIND_D(1,I) = 0.0D0
      WIND_D(2,I) = PRESSURE*AREA_OUT1(I)
      WIND_D(3,I) = PRESSURE*AREA_OUT2(I)
      CASE(2)
      WIND_D(1,I) = PRESSURE*AREA_OUT1(I)
      WIND_D(2,I) = 0.0D0
      WIND_D(3,I) = PRESSURE*AREA_OUT2(I)
      CASE(3)
      WIND_D(1,I) = PRESSURE*AREA_OUT1(I)
      WIND_D(2,I) = PRESSURE*AREA_OUT2(I)
      WIND_D(3,I) = 0.0D0
      ENDSELECT
      ENDDO
      
      END
C	==================================================================
      SUBROUTINE WALL_PRESSURE_COEFF (RATIO_L,RATIO_B,QZ,QH)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      QZ = 1.0D0
      
      IF (RATIO_L.LE.1.0D0)THEN
      QH = -0.5D0
      ELSEIF (RATIO_L.GT.1.0D0.AND.RATIO_L.LE.2.0)THEN
      DELTA_RATIO = 2.0D0 - 1.0D0
      DELTA_CP    = -0.3D0
      QH          = -0.5D0 + (DELTA_CP/DELTA_RATIO)*(RATIO_L-1.0D0)
         IF (RATIO_L.EQ.2D0) QH = -0.3D0
         IF (RATIO_L.EQ.1D0) QH = -0.5
      ENDIF
      
      IF (RATIO_B.LE.1.0D0)THEN
      QH = -0.5D0
      ELSEIF (RATIO_B.GT.1.0D0.AND.RATIO_B.LE.2.0)THEN
      DELTA_RATIO = 2.0D0 - 1.0D0
      DELTA_CP    = -0.3D0
      QH          = -0.5D0 + (DELTA_CP/DELTA_RATIO)*(RATIO_B-1.0D0)
         IF (RATIO_B.EQ.2D0) QH = -0.3D0
         IF (RATIO_B.EQ.1D0) QH = -0.5
      ENDIF
      END
C	==================================================================
      SUBROUTINE L_B_LENGTH (RATIO_L,RATIO_B)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0
      COMMON /DIAPH_DATA/ NSLAVE_NODE(500,100),NSLAVE(500)
      COMMON /MGRAV/ NGRAV
      ALLOCATABLE NSLAVE_NODE_FLOOR(:)
       
      DO I =1,IFLOOR
          ALLOCATE (NSLAVE_NODE_FLOOR(NSLAVE(I)))
          IF (I.EQ.IFLOOR)THEN ! TOP FLOOR
             DO J = 1,NSLAVE(I)
             SELECTCASE (NGRAV)
             CASE (1) ! X-DIRECTION
             CALL RELFILL('@XYZ',ELEVATION,1,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
	       CALL RELFILL('@XYZ',A1,2,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	       CALL RELFILL('@XYZ',A2,3,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
             CASE (2) ! Y-DIRECTION  
             CALL RELFILL('@XYZ',A1,1,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	       CALL RELFILL('@XYZ',ELEVATION,2,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
	       CALL RELFILL('@XYZ',A2,3,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
             CASE (3) ! Z-DIRECTION
             CALL RELFILL('@XYZ',A1,1,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	       CALL RELFILL('@XYZ',A2,2,NSLAVE_NODE_FLOOR(J),0)         !GETTING NODAL COORDINATE
	       CALL RELFILL('@XYZ',ELEVATION,3,NSLAVE_NODE_FLOOR(J),0)  !GETTING NODAL COORDINATE
             ENDSELECT 
             ENDDO
         ENDIF
         DEALLOCATE (NSLAVE_NODE_FLOOR)
      ENDDO
      
      SELECTCASE (NGRAV)
      CASE (1) ! X-DIRECTION
      RATIO_L = A1/A2 ! Y/Z
      RATIO_B = A2/A1 ! Z/Y
      CASE (2) ! Y-DIRECTION
      RATIO_L = A1/A2 ! X/Z
      RATIO_B = A2/A1 ! Z/X  
      CASE (3) ! Z-DIRECTION
      RATIO_L = A1/A2 ! X/Y
      RATIO_B = A2/A1 ! Y/X
      ENDSELECT
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
      SUBROUTINE VELOCITY_PRESSURE_KH (ELEVATION,IOPT_EXPOSURE,AKZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION BKZ_DATA(44),CKZ_DATA(44),DKZ_DATA(44)
      DIMENSION DATA_AKZ(44)
      DATA BKZ_DATA/0.70D0,0.70D0,0.70D0,0.70D0,0.76D0,0.81D0,0.85D0,0.89D0,0.93D0,0.96D0
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
            AKZ = 0.70D0
            EXIT
            ELSEIF (ELEVATION.LE.DATA_AKZ(22+I)) THEN
            IF(I.NE.1) DELTA_PRESSURE  = DATA_AKZ(I)    -  DATA_AKZ(I-1)
            IF(I.NE.1) DELTA_ELEVATION = DATA_AKZ(I+22) -  DATA_AKZ(I+22-1)
            AKZ = DATA_AKZ(I) + ((DELTA_PRESSURE/DELTA_ELEVATION)*(ELEVATION-DATA_AKZ(I+22-1)))
            EXIT 
            ELSEIF (ELEVATION.GT.DATA_AKZ(44))  THEN
            AKZ = 1.56D0
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
            AK1 = DATA_AK1(I) + ((DELTA_PRESSURE/DELTA_RATIO)*(RATIO-DATA_AK1(I+7-1)))
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
            AK2 = DATA_AK2(I) + ((DELTA_PRESSURE/DELTA_RATIO)*(RATIO-DATA_AK2(I+9-1)))
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
            AK3 = DATA_AK3(I) + ((DELTA_PRESSURE/DELTA_RATIO)*(RATIO-DATA_AK3(I+13-1)))
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
      