!      ORIGINAL OFFSHORE WAVE FORCE
!      DO NOT USED IN THIS PRESENT TIME      

!      ==============================================================================
!      GENERATE WAVE LOADING
!      ==============================================================================
!      RAMDA   = WAVE LENGTH
!      HIGHT   = WAVE HIGHTH
!      HEIGHT  = TOTAL HEIGHT OF STURCTURE FROM STILL WATER LEVEL
!      TIME    = TIME
!      RK      = WAVE NUMBER = 2*Pi/(RAMDA)
!      C       = WAVE SPEED 
!      DEPTH   = WATER DEPTH FROM STILL WATER LEVEL ( TAKE AS POSITIVE )
!      DEPTH0  = REFERENCE DEPTH FPT WIND-DENERATE CURRENT; DEPTH0=50M
!      W       = ANGULAR ACCELERATION
!      G       = GAVITY ACCLERATION
!      U       = THE HORIZONTAL WAVE-INDUCE VELOCITY
!      V       = VERTICAL WATER VELOCITY
!      AX      = THE HORIZONTAL WAVE-INDUCE ACCRLERATION
!      D       = DIAMETER OF STRUCTURE
!      ROLL    = THE MASS DENSITY OF THE FLUID ( WATER )
!      ROLLA   = THE DENSITY AND VISCOSITY OF THE AIR 
!      VWIND0  = WIND-GENERATE CURRENT AT STILL WATER LEVEL
!      VTIDE0  = TIDAL CURRENT AT SILL WATER LEVEL
!      VZ      = TOTAL CURRENT VELOCITY AT LEVEL Z
!      UM      = THE MAXIMUM VALUE OF THE ORBITAL VELOCITY AT THE BED ( LINEAR THEROY )
!      PERIOD  = WAVE PERIOD
!      HF      = ARM MEASURED FORM THE SEABED
!      CD      = DRAG COEFFICIENTS
!      CM      = INERTIA COEFFICIENTS
!      MASS    = MASS PER UNIT LENGTH
!      MASSP   = ADD MASS PER UNIT LENGTH OF PILE
!      A       = AREA ON THE STRUCTURE EXPOSED TO TH SLAMMING FORCE
!      ICASE   = SELECT CASE
!      H0      = REFERENCE DEPTH FOR WIND-GENERATED CURRENT; H0=50M
!      LEVEL   = WATER LEVEL ( ASSUME )
!      ALPHA   = EXPONET IN POWER-LAW MODEL FOR WIND SPEED PROFILE
!      UH      = 10 MINUTE MEAN WIND SPEEDV
!      PK      = SHAPE PARAMETER
!      HIGH    = REFERENCE HEIGHT
!      Z0      = THE ROUGHNESS LENGTH OF THE SURFACE
!      S       = THE PENETRATION DISTANCE S FOR A SECTION IN QUESTION IS THE HORIXONTAL DISTANCE
!                FROM THE PEROPHERY ON THE WET SIDE
!      X,Y,ZP  = COORDINATE WE INTERESTING
!      =============================================================================
      
       SUBROUTINE WAVE_LOADING(HIGHT,DEPTH,THIGHT,X,ZP,RAMDA,G,ROLL,ROLLA,ZETA,VTIDE0,VWIND0,H0,AP,S,CS,HIGH,UH,ALPHA,Z0,TIME,CD,CM
     1                         ,D,ICASE,Y,UT,VT,WT,TF)
	   IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*4 (I-N)
       DIMENSION U(5),V(5),AX(5),AY(5),FF(5)
     
!      (TRY) INPUT VALUE ***********************************************************
!       RAMDA=375D0                                                        
!       HIGHT=35D0
!       DEPTH=75D0
!       TIME=0.0D0
!       G=32.2D0
!       X=0.0D0
!       Y=0.0D0
!       ZP=0.0D0
!       ZETA=90D0
!       CD=1.0D0
!       CM=2.0D0
       
!       D=4.0D0
!       ROLL=1.99D0
 !     **********************************************************************************
 !     INPUT FOR CURRENT   
!       VTIDE0=20D0
!       VWIND0=20D0
 !     H0 REFERENCE DEPTH FOR WIND-GENERATED CURRENT; 50 M    
!       H0=50D0
!       ROLLA=0.00238D0
!       HEIGHT=120.0D0
 !     **********************************************************************************
 !     INPUT FOR PLUNGING   
!       AP=1
 !     **********************************************************************************
 !     INPUT FOR SURGING
!       S=1  
 !     ********************************************************************************** 
 !     INPUT FOR WATER LEVEL  
!       LEVEL=10
 !     **********************************************************************************
 !     INPUT FOR WIND LOAD
!       HIGH=10
 !       PK=2
 !       UH=210
 !       TR=1
 !       ALPHA=0.1D0
 !       ROLLA=0.00238D0
 !     **********************************************************************************
 !     INPUT FOR PLUNGING AND SURING
 !       CS=3.0
 !     **********************************************************************************
 !     INPUT FOR SELECT CASE
 !     CASE 1 WAVE LOAD
 !     CASE 2 CURRENT
 !     CASE 3 WATER LEVEL
 !     CASE 4 WAVE BREAK
 !     CASE 5 WIND LOAD
 !     -----------------
 !     SELECTCASE 1,2,3,4,5 
 !     **********************************************************************************       
 !     GENERATE VELOCITY ****************************************************************
 !     ********************************************************************************** 
       Pi=3.141592654
       RK=2*Pi/RAMDA
       RATIO=DEPTH/RAMDA
       MORRI=0.2*RAMDA
       
       SELECT CASE (ICASE)
       
       CASE(1)
!     ------------------------------------------------------------------------------
!     WHEN D>0.2*RAMDA MORISON'S EQUATION IS NOT VALID
!     -----------------------------------------------------------------------------
       IF (D.GT.MORRI)THEN
       DKR=D/2
       CALL NVMORISON (FX,ROLL,AW,RK,DEPTH,G,HF,D)
       GOTO 11
       ENDIF
!     -----------------------------------------------------------------------------
!     WHEN D<=0.2*RAMDA MORISON'S EQUATION IS VAID
!     -----------------------------------------------------------------------------   
       IF (D.LE.MORRI)THEN
       CALL VELOCITY (U,V,DEPTH,RK,G,TIME,RATIO,ZETA,X,Y,HIGHT,UT,VT,WT,ZP,CD,PI,ROLL,D,UTT,VTT,WTT,TY)
       CALL ACCELERATION (G,RK,DEPTH,X,Y,TIME,AX,AY,RATIO,AT,HIGHT,AXT,AYT,AZT,ZETA,CM,ROLL,Pi,D)
       UT=UTT+AXT
       VT=VTT+AYT
       WT=WTT+AZT
       GOTO 12
       ENDIF
       
!      WIND CURRENT LOAD    
       CASE(2)
       CALL CURRENT (VWIND0,DEPTH,H0,ROLLA,Y,TF,D)
       GOTO 12
       
!      WATER TIDE LOAD 
       CASE(3)
       CALL TIDE_LOAD (VTIDE0,DEPTH,ROLL,Y,THIGHT,TF,D)
       GOTO 12
     
 !     WATER LEVEL ( ASSUME IN KOREA 10 M )
       CASE(4)
 !     **********************************************************************************
 !     INPUT FOR SELECT CASE
 !     CASE 1 ASSUME IN KOREA 10 M
 !     CASE 2 USER DEFINE
 !     -----------------
 !     SELECTCASE 1,2 
       KCASE=1
       SELECT CASE (KCASE)
       CASE(1)
       LEVEL=10
       CASE(2)
       LEVEL=5
       ENDSELECT
!     **********************************************************************************
       DEPTH=LEVEL+DEPTH
!     ------------------------------------------------------------------------------
!     WHEN D>0.2*RAMDA MORISON'S EQUATION IS NOT VALID
!     -----------------------------------------------------------------------------
       IF (D.GT.MORRI)THEN
       DKR=D/2
       CALL NVMORISON (FX,ROLL,AW,RK,DEPTH,G,HF,D)
       GOTO 11
       ENDIF
!     -----------------------------------------------------------------------------
!     WHEN D<=0.2*RAMDA MORISON'S EQUATION IS VALID
!     -----------------------------------------------------------------------------
       IF (D.LE.MORRI)THEN
       CALL VELOCITY (U,V,DEPTH,RK,G,TIME,RATIO,ZETA,X,Y,HIGHT,UT,VT,WT,ZP,CD,PI,ROLL,D,UTT,VTT,WTT,TY)
       CALL ACCELERATION (G,RK,DEPTH,X,Y,TIME,AX,AY,RATIO,AT,HIGHT,AXT,AYT,AZT,ZETA,CM,ROLL,Pi,D)
       UT=UTT+AXT
       VT=VTT+AYT
       WT=WTT+AZT
       GOTO 12
       ENDIF
       
       
 !     GENERRATE FORCE FOR PLUNGING *************************************
       CASE(5)
       CALL VELOCITY (U,V,DEPTH,RK,G,TIME,RATIO,ZETA,X,Y,HIGHT,UT,VT,WT,ZP,CD,PI,ROLL,D,UTT,VTT,WTT,TY)
       CALL PLUNGING (CS,PI,ROLL,UT,AP,F2)
       TF=F2
       GOTO 12
       
 !     GENERRATE FORCE FOR SURING *************************************
       CASE(6)
       CALL VELOCITY (U,V,DEPTH,RK,G,TIME,RATIO,ZETA,X,Y,HIGHT,UT,VT,WT,ZP,CD,PI,ROLL,D,UTT,VTT,WTT,TY)
       CALL SURING (FSP,HIGHT,DEPTH,ROLL,PI,CM,CD,FF,S,D)
       TF=FSP
       GOTO 12
       
       CASE(7)
 !     **********************************************************************************
 !     CASE 1 OPEN SEA WITHOUT WAVE Z0=0.0001m
 !     CASE 2 COASTAL AREAS WITH ON-SHORE WIND
 !     CASE 3 THE ROUGHNESS PARAMETER MAY BE SOLVE IMPLICITLY FROM THE FOLLOWING EQUATION Z0=0.003m
 !            CASE 3.1 AC IS CHARNOCK'S CONSTANT =0.011 FOR OPEN SEA   AC=0.011
 !            CASE 3.2 AC IS CHARNOCK'S CONSTANT =0.034 FOR NEAR-COASTAL LOCATION   AC=0.034
 !             Z0=(AC/G)*((0.4*UH)/(LOG(YZ/Z0)))
 !     CASE 4 USER DEFINE  Z0
       CALL VWIND (HIGH,Z0,UH,ALPHA,TF,UR,ROLLA,Y,DEPTH,D) 
       
       
       ENDSELECT    
       
11     CONTINUE
12     CONTINUE
     
       RETURN
       END
!     ==============================================================================
!     SUBROUTINE INTERPOLATION_F ( VALUE OF WAVE PROFILE PARAMETHERS ) 
!     SOURCE:VALUE BASE ON SKJELBREIA AND HENDRICKSON (1961)
!     R10F22 MEAN R=REAL 10 = H/RAMDA=0.10 F22=F(22)
!     ==============================================================================
      SUBROUTINE INTERPOLATION_F (F22,F24,F33,F35,F44,F55,RATIO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!     F(22)-------------------------------------------------------------------------
      R10F22=3.892D0
      R15F22=1.539D0
      R20F22=0.927D0
      R25F22=0.699D0
      R30F22=0.599D0
      R35F22=0.551D0
      R40F22=0.527D0
      R50F22=0.507D0
      R60F22=0.502D0
!     F(24)-------------------------------------------------------------------------
      R10F24=-28.61D0
      R15F24=1.344D0
      R20F24=1.398D0
      R25F24=1.064D0
      R30F24=0.893D0
      R35F24=0.804D0
      R40F24=0.759D0
      R50F24=0.722D0
      R60F24=0.712D0
!     F(33)-------------------------------------------------------------------------
      R10F33=13.09D0
      R15F33=2.381D0
      R20F33=0.996D0
      R25F33=0.630D0
      R30F33=0.495D0
      R35F33=0.435D0
      R40F33=0.410D0
      R50F33=0.384D0
      R60F33=0.377D0
!     F(35)--------------------------------------------------------------------------
      R10F35=-138.6D0
      R15F35=6.935D0
      R20F35=3.679D0
      R25F35=2.244D0
      R30F35=1.685D0
      R35F35=1.438D0
      R40F35=1.330D0
      R50F35=1.230D0
      R60F35=1.205D0
!     F(44)-------------------------------------------------------------------------
      R10F44=44.99D0
      R15F44=4.147D0
      R20F44=1.259D0
      R25F44=0.676D0
      R30F44=0.484D0
      R35F44=0.407D0
      R40F44=0.371D0
      R50F44=0.344D0
      R60F44=0.337D0
!     F(55)------------------------------------------------------------------------
      R10F55=163.8D0
      R15F55=7.935D0
      R20F55=1.734D0
      R25F55=0.797D0
      R30F55=0.525D0
      R35F55=0.420D0
      R40F55=0.373D0
      R50F55=0.339D0
      R60F55=0.329D0
!     ------------------------------------------------------------------------------
!     INTERPOLATION BY ETC
      IF (RATIO.LT.0.1D0)THEN
      !WRITE (*,13)
      !WRITE (*,9)
      !WRITE (*,10) 
9     FORMAT ('---------------- OFFSHORE ERROR MESSAGE ( STOKE FIFTH ORDER )----------------')   
10    FORMAT ('- WAVE HIGHT / RAMDA LESS THAN 0.10')
13    FORMAT('')
      F22=0.0D0
      F24=0.0D0
      F33=0.0D0
      F35=0.0D0
      F44=0.0D0
      F55=0.0D0
      STOP
      ELSEIF (RATIO.GE.0.1D0.AND.RATIO.LE.0.15D0)THEN
      F22=(((R15F22-R10F22)/0.05D0)*(RATIO-0.1D0))+R10F22
      F24=(((R15F24-R10F24)/0.05D0)*(RATIO-0.1D0))+R10F24
      F33=(((R15F33-R10F33)/0.05D0)*(RATIO-0.1D0))+R10F33 
      F35=(((R15F35-R10F35)/0.05D0)*(RATIO-0.1D0))+R10F35
      F44=(((R15F44-R10F44)/0.05D0)*(RATIO-0.1D0))+R10F44
      F55=(((R15F55-R10F55)/0.05D0)*(RATIO-0.1D0))+R10F55
      ELSEIF (RATIO.GE.0.15D0.AND.RATIO.LE.0.20D0)THEN
      F22=(((R20F22-R15F22)/0.05D0)*(RATIO-0.15D0))+R15F22
      F24=(((R20F24-R15F24)/0.05D0)*(RATIO-0.15D0))+R15F24    
      F33=(((R20F33-R15F33)/0.05D0)*(RATIO-0.15D0))+R15F33
      F35=(((R20F35-R15F35)/0.05D0)*(RATIO-0.15D0))+R15F35
      F44=(((R20F44-R15F44)/0.05D0)*(RATIO-0.15D0))+R15F44
      F55=(((R20F55-R15F55)/0.05D0)*(RATIO-0.15D0))+R15F55
      ELSEIF (RATIO.GE.0.2D0.AND.RATIO.LE.0.25D0)THEN
      F22=(((R25F22-R20F22)/0.05D0)*(RATIO-0.2D0))+R20F22
      F24=(((R25F24-R20F24)/0.05D0)*(RATIO-0.2D0))+R20F24
      F33=(((R25F33-R20F33)/0.05D0)*(RATIO-0.2D0))+R20F33
      F35=(((R25F35-R20F35)/0.05D0)*(RATIO-0.2D0))+R20F35
      F44=(((R25F44-R20F44)/0.05D0)*(RATIO-0.2D0))+R20F44
      F55=(((R25F55-R20F55)/0.05D0)*(RATIO-0.2D0))+R20F55
      ELSEIF (RATIO.GE.0.25D0.AND.RATIO.LE.0.30D0)THEN
      F22=(((R30F22-R25F22)/0.05D0)*(RATIO-0.25D0))+R25F22
      F24=(((R30F24-R25F24)/0.05D0)*(RATIO-0.25D0))+R25F24
      F33=(((R30F33-R25F33)/0.05D0)*(RATIO-0.25D0))+R25F33
      F35=(((R30F35-R25F35)/0.05D0)*(RATIO-0.25D0))+R25F35
      F44=(((R30F44-R25F44)/0.05D0)*(RATIO-0.25D0))+R25F44
      F55=(((R30F55-R25F55)/0.05D0)*(RATIO-0.25D0))+R25F55
      ELSEIF (RATIO.GE.0.30D0.AND.RATIO.LE.0.35D0)THEN
      F22=(((R35F22-R30F22)/0.05D0)*(RATIO-0.30D0))+R30F22   
      F24=(((R35F24-R30F24)/0.05D0)*(RATIO-0.30D0))+R30F24
      F33=(((R35F33-R30F33)/0.05D0)*(RATIO-0.30D0))+R30F33
      F35=(((R35F35-R30F35)/0.05D0)*(RATIO-0.30D0))+R30F35
      F44=(((R35F44-R30F44)/0.05D0)*(RATIO-0.30D0))+R30F44
      F55=(((R35F55-R30F55)/0.05D0)*(RATIO-0.30D0))+R30F55
      ELSEIF (RATIO.GE.0.35D0.AND.RATIO.LE.0.40D0)THEN
      F22=(((R40F22-R35F22)/0.05D0)*(RATIO-0.35D0))+R35F22
      F24=(((R40F24-R35F24)/0.05D0)*(RATIO-0.35D0))+R35F24
      F33=(((R40F33-R35F33)/0.05D0)*(RATIO-0.35D0))+R35F33
      F35=(((R40F35-R35F35)/0.05D0)*(RATIO-0.35D0))+R35F35
      F44=(((R40F44-R35F44)/0.05D0)*(RATIO-0.35D0))+R35F44
      F55=(((R40F55-R35F55)/0.05D0)*(RATIO-0.35D0))+R35F55
      ELSEIF (RATIO.GE.0.40D0.AND.RATIO.LE.0.50D0)THEN
      F22=(((R50F22-R40F22)/0.1D0)*(RATIO-0.40D0))+R40F22
      F24=(((R50F24-R40F24)/0.1D0)*(RATIO-0.40D0))+R40F24
      F33=(((R50F33-R40F33)/0.1D0)*(RATIO-0.40D0))+R40F33
      F35=(((R50F35-R40F35)/0.1D0)*(RATIO-0.40D0))+R40F35
      F44=(((R50F44-R40F44)/0.1D0)*(RATIO-0.40D0))+R40F44
      F55=(((R50F55-R40F55)/0.1D0)*(RATIO-0.40D0))+R40F55
      ELSEIF (RATIO.GE.0.50D0.AND.RATIO.LE.0.60D0)THEN
      F22=(((R60F22-R50F22)/0.1D0)*(RATIO-0.50D0))+R50F22
      F24=(((R60F24-R50F24)/0.1D0)*(RATIO-0.50D0))+R50F24
      F33=(((R60F33-R50F33)/0.1D0)*(RATIO-0.50D0))+R50F33
      F35=(((R60F35-R50F35)/0.1D0)*(RATIO-0.50D0))+R50F35
      F44=(((R60F44-R50F44)/0.1D0)*(RATIO-0.50D0))+R50F44
      F55=(((R60F55-R50F55)/0.1D0)*(RATIO-0.50D0))+R50F55
      ELSEIF (RATIO.GT.0.60D0)THEN
      !WRITE (*,16)
      !WRITE (*,14)
      !WRITE (*,15) 
14     FORMAT ('---------------- OFFSHORE ERROR MESSAGE ( STOKE FIFTH ORDER )----------------')   
15    FORMAT ('- WAVE HIGHT / RAMDA GREATER THAN 0.60')
16    FORMAT('')
      F22=0.0D0
      F24=0.0D0
      F33=0.0D0
      F35=0.0D0
      F44=0.0D0
      F55=0.0D0
      STOP 
      ENDIF
      END
!     ==============================================================================
!     SUBROUTINE INTERPOLATION_G ( VALUE OF VELOCITY PARAMETHERS ) 
!     SOURCE:VALUE BASE ON SKJELBREIA AND HENDRICKSON (1961)
!     ==============================================================================
      SUBROUTINE INTERPOLATION_G (G11,G13,G15,G22,G24,G33,G35,G44,G55,RATIO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!     G(11)-------------------------------------------------------------------------
      G11=1.000D0
!     G(13)-------------------------------------------------------------------------
      R10G13=-7.394D0
      R15G13=-2.320D0
      R20G13=-1.263D0
      R25G13=-0.911D0
      R30G13=-0.765D0
      R35G13=-0.696D0
      R40G13=-0.662D0
      R50G13=-0.635D0
      R60G13=-0.628D0
!     G(15)-------------------------------------------------------------------------
      R10G15=-12.73D0
      R15G15=-4.864D0
      R20G15=-2.266D0
      R25G15=-1.415D0
      R30G15=-1.077D0
      R35G15=-0.925D0
      R40G15=-0.850D0
      R50G15=-0.790D0
      R60G15=-0.777D0
!     G(22)-------------------------------------------------------------------------
      R10G22=2.996D0
      R15G22=0.860D0
      R20G22=0.326D0
      R25G22=0.154D0
      R30G22=0.076D0
      R35G22=0.038D0
      R40G22=0.020D0
      R50G22=0.006D0
      R60G22=0.002D0
!     G(24)-------------------------------------------------------------------------
      R10G24=-48.14D0
      R15G24=-0.907D0
      R20G24=0.680D0
      R25G24=0.673D0
      R30G24=0.601D0
      R35G24=0.556D0
      R40G24=0.528D0
      R50G24=0.503D0
      R60G24=0.502D0
!     G(33)-------------------------------------------------------------------------
      R10G33=5.942D0
      R15G33=0.310D0
      R20G33=-0.017D0
      R25G33=-0.030D0
      R30G33=-0.020D0
      R35G33=-0.012D0
      R40G33=-0.006D0
      R50G33=-0.002D0
      R60G33=-0.001D0
!     G(35)-------------------------------------------------------------------------
      R10G35=-121.7D0
      R15G35=2.843D0
      R20G35=1.093D0
      R25G35=0.440D0
      R30G35=0.231D0
      R35G35=0.152D0
      R40G35=0.117D0
      R50G35=0.092D0
      R60G35=0.086D0
!     G(44)-------------------------------------------------------------------------
      R10G44=7.671D0
      R15G44=-0.167D0
      R20G44=-0.044D0
      R25G44=-0.005D0
      R30G44=0.002D0
      R35G44=0.002D0
      R40G44=0.001D0
      R50G44=0.000D0
      R60G44=0.000D0
!     G(55)-------------------------------------------------------------------------
      R10G55=0.892D0
      R15G55=-0.257D0
      R20G55=0.006D0
      R25G55=0.005D0
      R30G55=0.001D0
      R35G55=0.000D0
      R40G55=0.000D0
      R50G55=0.000D0
      R60G55=0.000D0
!     ------------------------------------------------------------------------------
!     INTERPOLATION BY ETC
      IF (RATIO.LT.0.1D0)THEN
!      WRITE (6,30) 
30    FORMAT (2X,'H/L < 0.10 ( NOT INCLUDE IN PROGRAMS )')
      ELSEIF (RATIO.GE.0.1D0.AND.RATIO.LE.0.15D0)THEN
      G13=(((R15G13-R10G13)/0.05D0)*(RATIO-0.1D0))+R10G13
      G15=(((R15G15-R10G15)/0.05D0)*(RATIO-0.1D0))+R10G15
      G22=(((R15G22-R10G22)/0.05D0)*(RATIO-0.1D0))+R10G22
      G24=(((R15G24-R10G24)/0.05D0)*(RATIO-0.1D0))+R10G24
      G33=(((R15G33-R10G33)/0.05D0)*(RATIO-0.1D0))+R10G33
      G35=(((R15G35-R10G35)/0.05D0)*(RATIO-0.1D0))+R10G35
      G44=(((R15G44-R10G44)/0.05D0)*(RATIO-0.1D0))+R10G44
      G55=(((R15G55-R10G55)/0.05D0)*(RATIO-0.1D0))+R10G55
      ELSEIF (RATIO.GE.0.15D0.AND.RATIO.LE.0.20D0)THEN
      G13=(((R20G13-R15G13)/0.05D0)*(RATIO-0.15D0))+R15G13
      G15=(((R20G15-R15G15)/0.05D0)*(RATIO-0.15D0))+R15G15
      G22=(((R20G22-R15G22)/0.05D0)*(RATIO-0.15D0))+R15G22
      G24=(((R20G24-R15G24)/0.05D0)*(RATIO-0.15D0))+R15G24
      G33=(((R20G33-R15G33)/0.05D0)*(RATIO-0.15D0))+R15G33
      G35=(((R20G35-R15G35)/0.05D0)*(RATIO-0.15D0))+R15G35
      G44=(((R20G44-R15G44)/0.05D0)*(RATIO-0.15D0))+R15G44
      G55=(((R20G55-R15G55)/0.05D0)*(RATIO-0.15D0))+R15G55
      ELSEIF (RATIO.GE.0.20D0.AND.RATIO.LE.0.25D0)THEN
      G13=(((R25G13-R20G13)/0.05D0)*(RATIO-0.2D0))+R20G13
      G15=(((R25G15-R20G15)/0.05D0)*(RATIO-0.2D0))+R20G15
      G22=(((R25G22-R20G22)/0.05D0)*(RATIO-0.2D0))+R20G22
      G24=(((R25G24-R20G24)/0.05D0)*(RATIO-0.2D0))+R20G24
      G33=(((R25G33-R20G33)/0.05D0)*(RATIO-0.2D0))+R20G33
      G35=(((R25G35-R20G35)/0.05D0)*(RATIO-0.2D0))+R20G35
      G44=(((R25G44-R20G44)/0.05D0)*(RATIO-0.2D0))+R20G44
      G55=(((R25G55-R20G55)/0.05D0)*(RATIO-0.2D0))+R20G55
      ELSEIF (RATIO.GE.0.25D0.AND.RATIO.LE.0.30D0)THEN
      G13=(((R30G13-R25G13)/0.05D0)*(RATIO-0.25D0))+R25G13
      G15=(((R30G15-R25G15)/0.05D0)*(RATIO-0.25D0))+R25G15
      G22=(((R30G22-R25G22)/0.05D0)*(RATIO-0.25D0))+R25G22
      G24=(((R30G24-R25G24)/0.05D0)*(RATIO-0.25D0))+R25G24
      G33=(((R30G33-R25G33)/0.05D0)*(RATIO-0.25D0))+R25G33
      G35=(((R30G35-R25G35)/0.05D0)*(RATIO-0.25D0))+R25G35
      G44=(((R30G44-R25G44)/0.05D0)*(RATIO-0.25D0))+R25G44
      G55=(((R30G55-R25G55)/0.05D0)*(RATIO-0.25D0))+R25G55
      ELSEIF (RATIO.GE.0.30D0.AND.RATIO.LE.0.35D0)THEN
      G13=(((R35G13-R30G13)/0.05D0)*(RATIO-0.30D0))+R30G13
      G15=(((R35G15-R30G15)/0.05D0)*(RATIO-0.30D0))+R30G15
      G22=(((R35G22-R30G22)/0.05D0)*(RATIO-0.30D0))+R30G22
      G24=(((R35G24-R30G24)/0.05D0)*(RATIO-0.30D0))+R30G24
      G33=(((R35G33-R30G33)/0.05D0)*(RATIO-0.30D0))+R30G33
      G35=(((R35G35-R30G35)/0.05D0)*(RATIO-0.30D0))+R30G35
      G44=(((R35G44-R30G44)/0.05D0)*(RATIO-0.30D0))+R30G44
      G55=(((R35G55-R30G55)/0.05D0)*(RATIO-0.30D0))+R30G55
      ELSEIF (RATIO.GE.0.35D0.AND.RATIO.LE.0.40D0)THEN
      G13=(((R40G13-R35G13)/0.05D0)*(RATIO-0.35D0))+R35G13
      G15=(((R40G15-R35G15)/0.05D0)*(RATIO-0.35D0))+R35G15
      G22=(((R40G22-R35G22)/0.05D0)*(RATIO-0.35D0))+R35G22
      G24=(((R40G24-R35G24)/0.05D0)*(RATIO-0.35D0))+R35G24
      G33=(((R40G33-R35G33)/0.05D0)*(RATIO-0.35D0))+R35G33
      G35=(((R40G35-R35G35)/0.05D0)*(RATIO-0.35D0))+R35G35
      G44=(((R40G44-R35G44)/0.05D0)*(RATIO-0.35D0))+R35G44
      G55=(((R40G55-R35G55)/0.05D0)*(RATIO-0.35D0))+R35G55
      ELSEIF (RATIO.GE.0.40D0.AND.RATIO.LE.0.50D0)THEN
      G13=(((R50G13-R40G13)/0.1D0)*(RATIO-0.40D0))+R40G13
      G15=(((R50G15-R40G15)/0.1D0)*(RATIO-0.40D0))+R40G15
      G22=(((R50G22-R40G22)/0.1D0)*(RATIO-0.40D0))+R40G22
      G24=(((R50G24-R40G24)/0.1D0)*(RATIO-0.40D0))+R40G24
      G33=(((R50G33-R40G33)/0.1D0)*(RATIO-0.40D0))+R40G33
      G35=(((R50G35-R40G35)/0.1D0)*(RATIO-0.40D0))+R40G35
      G44=(((R50G44-R40G44)/0.1D0)*(RATIO-0.40D0))+R40G44
      G55=(((R50G55-R40G55)/0.1D0)*(RATIO-0.40D0))+R40G55
      ELSEIF (RATIO.GE.0.50D0.AND.RATIO.LE.0.60D0)THEN
      G13=(((R60G13-R50G13)/0.1D0)*(RATIO-0.50D0))+R50G13
      G15=(((R60G15-R50G15)/0.1D0)*(RATIO-0.50D0))+R50G15
      G22=(((R60G22-R50G22)/0.1D0)*(RATIO-0.50D0))+R50G22
      G24=(((R60G24-R50G24)/0.1D0)*(RATIO-0.50D0))+R50G24
      G33=(((R60G33-R50G33)/0.1D0)*(RATIO-0.50D0))+R50G33
      G35=(((R60G35-R50G35)/0.1D0)*(RATIO-0.50D0))+R50G35
      G44=(((R60G44-R50G44)/0.1D0)*(RATIO-0.50D0))+R50G44
      G55=(((R60G55-R50G55)/0.1D0)*(RATIO-0.50D0))+R50G55
      ELSEIF (RATIO.GT.0.60D0)THEN
!      WRITE (6,40)
!40    FORMAT (2X,'H/L >0.6 ( NOT INCLUDE IN PROGRAMS )') 
      ENDIF
      END
!     ==============================================================================
!     SUBROUTINE INTERPOLATION_C ( FREQUENCY AND PRESSURE PARAMETHERS ) 
!     SOURCE:VALUE BASE ON SKJELBREIA AND HENDRICKSON (1961)
!     ==============================================================================
      SUBROUTINE INTERPOLATION_C (C1,C2,C3,C4,RATIO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!     C(1)-------------------------------------------------------------------------
      R10C1=8.791D0
      R15C1=2.646D0
      R20C1=1.549D0
      R25C1=1.229D0
      R30C1=1.107D0
      R35C1=1.055D0
      R40C1=1.027D0
      R50C1=1.008D0
      R60C1=1.002D0
!     C(2)-------------------------------------------------------------------------
      R10C2=383.7D0
      R15C2=19.82D0
      R20C2=5.044D0
      R25C2=2.568D0
      R30C2=1.833D0
      R35C2=1.532D0
      R40C2=1.393D0
      R50C2=1.283D0
      R60C2=1.240D0
!     C(3)-------------------------------------------------------------------------
      R10C3=-0.310D0
      R15C3=-0.155D0
      R20C3=-0.082D0
      R25C3=-0.043D0
      R30C3=-0.023D0
      R35C3=-0.012D0
      R40C3=-0.007D0
      R50C3=-0.001D0
      R60C3=-0.001D0
!     C(4)-------------------------------------------------------------------------
      R10C4=-0.060D0
      R15C4=0.257D0
      R20C4=0.077D0
      R25C4=0.028D0
      R30C4=0.010D0
      R35C4=0.004D0
      R40C4=0.002D0
      R50C4=0.000D0
      R60C4=0.000D0
!     ------------------------------------------------------------------------------
!     INTERPOLATION BY ETC
      IF (RATIO.LT.0.1D0)THEN
!      WRITE (6,50) 
50    FORMAT (2X,'H/L < 0.10 ( NOT INCLUDE IN PROGRAMS )')
      ELSEIF (RATIO.GE.0.1D0.AND.RATIO.LE.0.15D0)THEN
      C1=(((R15C1-R10C1)/0.05D0)*(RATIO-0.1D0))+R10C1
      C2=(((R15C2-R10C2)/0.05D0)*(RATIO-0.1D0))+R10C2
      C3=(((R15C3-R10C3)/0.05D0)*(RATIO-0.1D0))+R10C3
      C4=(((R15C4-R10C4)/0.05D0)*(RATIO-0.1D0))+R10C4
      ELSEIF (RATIO.GE.0.15D0.AND.RATIO.LE.0.20D0)THEN
      C1=(((R20C1-R15C1)/0.05D0)*(RATIO-0.15D0))+R15C1
      C2=(((R20C2-R15C2)/0.05D0)*(RATIO-0.15D0))+R15C2
      C3=(((R20C3-R15C3)/0.05D0)*(RATIO-0.15D0))+R15C3
      C4=(((R20C4-R15C4)/0.05D0)*(RATIO-0.15D0))+R15C4
      ELSEIF (RATIO.GE.0.2D0.AND.RATIO.LE.0.25D0)THEN
      C1=(((R25C1-R20C1)/0.05D0)*(RATIO-0.2D0))+R20C1
      C2=(((R25C2-R20C2)/0.05D0)*(RATIO-0.2D0))+R20C2
      C3=(((R25C3-R20C3)/0.05D0)*(RATIO-0.2D0))+R20C3
      C4=(((R25C4-R20C4)/0.05D0)*(RATIO-0.2D0))+R20C4
      ELSEIF (RATIO.GE.0.25D0.AND.RATIO.LE.0.30D0)THEN
      C1=(((R30C1-R25C1)/0.05D0)*(RATIO-0.25D0))+R25C1
      C2=(((R30C2-R25C2)/0.05D0)*(RATIO-0.25D0))+R25C2
      C3=(((R30C3-R25C3)/0.05D0)*(RATIO-0.25D0))+R25C3
      C4=(((R30C4-R25C4)/0.05D0)*(RATIO-0.25D0))+R25C4
      ELSEIF (RATIO.GE.0.30D0.AND.RATIO.LE.0.35D0)THEN
      C1=(((R35C1-R30C1)/0.05D0)*(RATIO-0.30D0))+R30C1
      C2=(((R35C2-R30C2)/0.05D0)*(RATIO-0.30D0))+R30C2
      C3=(((R35C3-R30C3)/0.05D0)*(RATIO-0.30D0))+R30C3
      C4=(((R35C4-R30C4)/0.05D0)*(RATIO-0.30D0))+R30C4
      ELSEIF (RATIO.GE.0.35D0.AND.RATIO.LE.0.40D0)THEN
      C1=(((R40C1-R35C1)/0.05D0)*(RATIO-0.35D0))+R35C1
      C2=(((R40C2-R35C2)/0.05D0)*(RATIO-0.35D0))+R35C2
      C3=(((R40C3-R35C3)/0.05D0)*(RATIO-0.35D0))+R35C3
      C4=(((R40C4-R35C4)/0.05D0)*(RATIO-0.35D0))+R35C40
      ELSEIF (RATIO.GE.0.40D0.AND.RATIO.LE.0.50D0)THEN
      C1=(((R50C1-R40C1)/0.1D0)*(RATIO-0.40D0))+R40C1
      C2=(((R50C2-R40C2)/0.1D0)*(RATIO-0.40D0))+R40C2
      C3=(((R50C3-R40C3)/0.1D0)*(RATIO-0.40D0))+R40C3
      C4=(((R50C4-R40C4)/0.1D0)*(RATIO-0.40D0))+R40C4
      ELSEIF (RATIO.GE.0.50D0.AND.RATIO.LE.0.60D0)THEN
      C1=(((R60C1-R50C1)/0.1D0)*(RATIO-0.50D0))+R50C1
      C2=(((R60C2-R50C2)/0.1D0)*(RATIO-0.50D0))+R50C2
      C3=(((R60C3-R50C3)/0.1D0)*(RATIO-0.50D0))+R50C3
      C4=(((R60C4-R50C4)/0.1D0)*(RATIO-0.50D0))+R50C4
      ELSEIF (RATIO.GT.0.60D0)THEN
!      WRITE (6,60)
!60    FORMAT (2X,'H/L >0.6 ( NOT INCLUDE IN PROGRAMS )') 
      ENDIF
      END
!     ==============================================================================
!     SUBROUTINE WAVE_HIGHT ( WAVE HIGHT PARAMETER ) 
!     ==============================================================================
      SUBROUTINE WAVE_HIGHT (RK,HIGHT,A,RATIO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION AI(6),X(3)
      REAL*8 FA,K,H,DEL1,DEL0,A
      CALL INTERPOLATION_F (F22,F24,F33,F35,F44,F55,RATIO)
       X0=-99
       X1=0
       X2=100
300    FA=F35+F55
       AI(1)=-RK*HIGHT
       AI(2)=2D0
       AI(3)=0D0
       AI(4)=2*F33
       AI(5)=0D0
       AI(6)=2*FA
       F0=(AI(6)*(X0**5))+(AI(4)*(X0**3))+(AI(2)*X0)+AI(1)
       F1=(AI(6)*(X1**5))+(AI(4)*(X1**3))+(AI(2)*X1)+AI(1)
       F2=(AI(6)*(X2**5))+(AI(4)*(X2**3))+(AI(2)*X2)+AI(1)
       H0=X1-X0
       H1=X2-X1
       DEL0=(F1-F0)/(X1-X0)
       DEL1=(F2-F1)/(X2-X1)
       AA=(DEL1-DEL0)/(H1+H0)
       BB=AA*H1+DEL1
       CC=F2
       DISC=SQRT(BB**2-(4*AA*CC))
       X3=X2+((-2*CC)/(BB+DISC))
       ERROR=ABS((X3-X2)/X3)*100
       X(1)=X1
       X(2)=X2
       X(3)=X3
      IF (ERROR.GT.0.001)THEN
         X0=X(1)
         X1=X(2)
         X2=X(3)
      GO TO 300
      ENDIF
      A=X2
      END         
!     ==============================================================================
!     SUBROUTINE VELOCITY ( HORIZONTAL AND VERTICAL VELOCITY,VELOCITY COEFFICIENTS ) 
!     ==============================================================================
      SUBROUTINE VELOCITY (U,V,DEPTH,RK,G,TIME,RATIO,ZETA,X,Y,HIGHT,UT,VT,WT,ZP,CD,PI,ROLL,D,UTT,VTT,WTT,TY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION U(5),V(5),UU(5),FN(10),VV(5),WZ(5),F(5),FU(5),FD(5)
      CALL INTERPOLATION_F (F22,F24,F33,F35,F44,F55,RATIO)
      CALL INTERPOLATION_G (G11,G13,G15,G22,G24,G33,G35,G44,G55,RATIO)
      CALL INTERPOLATION_C (C1,C2,C3,C4,RATIO)
      CALL WAVE_HIGHT (RK,HIGHT,A,RATIO)
      WS=(G*RK)
      WSS=(1+((A**2)*C1)+((A**4)*C2))
      WSSS=TANH(RK*DEPTH)
      W=SQRT(WS*WSS*WSSS)
      !============================== 
      W = 0.813039803031561
      TIME = 6.0D0/W
      !============================== 
      G1=(A*G11)+((A**3)*G13)+((A**5)*G15)
      G2=2.0*(((A**2)*G22)+((A**4)*G24))
      G3=3.0*((A**3)*G33)+((A**5)*G35)
      G4=4.0*(A**4)*G44
      G5=5.0*A**5*G55
!     THE FREE-SURFACE WATER DEFLECTION ---------------------------------------------
      F(1)=A
      F(2)=((A**2)*F22)+((A**4)*F24)
      F(3)=((A**3)*F33)+((A**5)*F35)
      F(4)=((A**4)*F44)
      F(5)=((A**5)*F55)
      AA=((RK*Z*COS(ZETA))+(RK*X*COS(ZETA))-(W*TIME))
      DO I=1,5
      FU(I)=(1/RK)*F(I)*COS(I*AA)
      ENDDO
!     DEFLECTION
      TFNU=FU(1)+FU(2)+FU(3)+FU(4)+FU(5)
      TY=DEPTH+TFNU
!     FOR X HORIZONTAL VELOCITY -----------------------------------------------------
      U(1)=(W/RK)*G1*(COS(ZETA))*(COSH(RK*Y)/SINH(RK*DEPTH))*COS((RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      U(2)=(W/RK)*G2*(COS(ZETA))*(COSH(2.0*RK*Y)/SINH(2.0*RK*DEPTH))*COS(2.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      U(3)=(W/RK)*G3*(COS(ZETA))*(COSH(3.0*RK*Y)/SINH(3.0*RK*DEPTH))*COS(3.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      U(4)=(W/RK)*G4*(COS(ZETA))*(COSH(4.0*RK*Y)/SINH(4.0*RK*DEPTH))*COS(4.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      U(5)=(W/RK)*G5*(COS(ZETA))*(COSH(5.0*RK*Y)/SINH(5.0*RK*DEPTH))*COS(5.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
!     FOR Y VERTICAL VELOCITY --------------------------------------------------------
      V(1)=(W/RK)*G1*(SINH(RK*Y)/SINH(RK*DEPTH))*SIN((RK*ZP*SIN(ZETA))+(RK*X*COS(ZETA))-(W*TIME))
      V(2)=(W/RK)*G2*(SINH(2.0*RK*Y)/SINH(2.0*RK*DEPTH))*SIN(2.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      V(3)=(W/RK)*G3*(SINH(3.0*RK*Y)/SINH(3.0*RK*DEPTH))*SIN(3.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      V(4)=(W/RK)*G4*(SINH(4.0*RK*Y)/SINH(4.0*RK*DEPTH))*SIN(4.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      V(5)=(W/RK)*G5*(SINH(5.0*RK*Y)/SINH(5.0*RK*DEPTH))*SIN(5.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
!     FOR Z HORIZONTAL VELOCITY ------------------------------------------------------
      WZ(1)=(W/RK)*G1*(SIN(ZETA))*(COSH(RK*Y)/SINH(RK*DEPTH))*COS((RK*ZP*SIN(ZETA))+(RK*X*COS(ZETA))-(W*TIME))
      WZ(2)=(W/RK)*G2*(SIN(ZETA))*(COSH(2.0*RK*Y)/SINH(2.0*RK*DEPTH))*COS(2.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      WZ(3)=(W/RK)*G3*(SIN(ZETA))*(COSH(3.0*RK*Y)/SINH(3.0*RK*DEPTH))*COS(3.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      WZ(4)=(W/RK)*G4*(SIN(ZETA))*(COSH(4.0*RK*Y)/SINH(4.0*RK*DEPTH))*COS(4.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      WZ(5)=(W/RK)*G5*(SIN(ZETA))*(COSH(5.0*RK*Y)/SINH(5.0*RK*DEPTH))*COS(5.0*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      UT=U(1)+U(2)+U(3)+U(4)+U(5)
      VT=V(1)+V(2)+V(3)+V(4)+V(5)
      WT=WZ(1)+WZ(2)+WZ(3)+WZ(4)+WZ(5)
      UTT=(0.50*ROLL*CD*D)*(U(1)+U(2)+U(3)+U(4)+U(5))*ABS(U(1)+U(2)+U(3)+U(4)+U(5))
      VTT=(0.50*ROLL*CD*D)*(V(1)+V(2)+V(3)+V(4)+V(5))*ABS(V(1)+V(2)+V(3)+V(4)+V(5))
      WTT=(0.50*ROLL*CD*D)*(WZ(1)+WZ(2)+WZ(3)+WZ(4)+WZ(5))*ABS(WZ(1)+WZ(2)+WZ(3)+WZ(4)+WZ(5))
      END
!     ===============================================================================
!     SUBROUTINE ACCELERATION( HORIZONTAL AND VERTICAL ACCELERATION,VELOCITY COEFFICIENTS ) 
!     ===============================================================================
      SUBROUTINE ACCELERATION (G,RK,DEPTH,X,Y,TIME,AX,AY,RATIO,AT,HIGHT,AXT,AYT,AZT,ZETA,CM,ROLL,Pi,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION AX(5),AY(5),AU(5),AV(5),GG(5),R(5),S(5),AZ(5)
      CALL INTERPOLATION_G (G11,G13,G15,G22,G24,G33,G35,G44,G55,RATIO)
      CALL INTERPOLATION_C (C1,C2,C3,C4,RATIO)
      CALL WAVE_HIGHT (RK,HIGHT,A,RATIO)
      WS=(G*RK)
      WSS=(1+((A**2)*C1)+((A**4)*C2))
      WSSS=TANH(RK*DEPTH)
      W=SQRT(WS*WSS*WSSS)
      !============================== 
      W = 0.813039803031561
      TIME = 6.0D0/W
      !============================== 
      C=SQRT((G/RK)*(1.0+(A**2)*C1+(A**4)*C2)*TANH(RK*DEPTH))
      GG(1)=(A*G11)+((A**3)*G13)+((A**5)*G15)
      GG(2)=2.0*(((A**2)*G22)+((A**4)*G24))
      GG(3)=3.0*((A**3)*G33)+((A**5)*G35)
      GG(4)=4.0*(A**4)*G44
      GG(5)=5.0*A**5*G55
!     FOR VELOCITY COEFFICIENTS -------------------------------------------
      DO I=1,5
      AU(I)=GG(I)*(COSH(I*RK*Y)/SINH(I*RK*DEPTH))
      AV(I)=GG(I)*(SINH(I*RK*Y)/SINH(I*RK*DEPTH))
      ENDDO
!     CALCULATE R(N) S(N) -------------------------------------------------   
      R(1)=2.0*AU(1)-AU(1)*AU(2)-AV(1)*AV(2)-AU(2)*AU(3)-AV(2)*AV(3)
      R(2)=4.0*AU(2)-AU(1)**2+AV(1)**2-2.0*AU(1)*AU(3)-2.0*AV(1)*AV(3)
      R(3)=6.0*AU(3)-3.0*AU(1)*AU(2)+3.0*AV(1)*AV(2)-3.0*AU(1)*AU(4)-3.0*AV(1)*AV(4)
      R(4)=8.0*AU(4)-2.0*AU(2)**2+2.0*AV(2)**2-4.0*AU(1)*AU(3)+4.0*AV(1)*AV(3)
      R(5)=10.0*AU(5)-5.0*AU(1)*AU(4)-5.0*AU(2)*AU(3)+5.0*AV(1)*AV(4)+5.0*AV(2)*AV(3)
      S(1)=-2.0*AU(1)*AV(1)+4*AU(2)*AV(1)
      S(2)=2.0*AV(1)-3.0*AU(1)*AV(2)-3.0*AU(2)*AV(1)-5.0*AU(2)*AV(3)-5.0*AU(5)*AV(2)
      S(3)=6.0*AV(3)-AU(1)*AV(2)+AU(2)*AV(1)-5.0*AU(1)*AV(4)-5.0*AU(4)*AV(1)
      S(4)=8.0*AV(4)-2.0*AU(1)*AV(3)+2.0*AU(3)*AV(1)+4.0*AU(2)*AV(2)
      S(5)=10.0*AV(5)-3.0*AU(1)*AV(4)+3.0*AU(4)*AV(1)-AU(2)*AV(3)+AU(3)*AV(2)
!     CALCULATE HORIZONTAL ACCELERATION -----------------------------------
      DO I=1,5
      AX(I)=(1.0/2.0)*(RK*(C**2))*R(I)*SIN(I*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))*COS(ZETA)
!     CALCULATE HORIZONTAL ACCELERATION -----------------------------------
      AY(I)=(-1.0/2.0)*(RK*(C**2))*S(I)*COS(I*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AZ(I)=(1.0/2.0)*SIN(ZETA)*(RK*(C**2))*R(I)*SIN(I*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      ENDDO
      AXT=(ROLL*CM*Pi*(D**2)/4)*(AX(1)+AX(2)+AX(3)+AX(4)+AX(5))
      AYT=(ROLL*CM*Pi*(D**2)/4)*(AY(1)+AY(2)+AY(3)+AY(4)+AY(5))
      AZT=(ROLL*CM*Pi*(D**2)/4)*(AZ(1)+AZ(2)+AZ(3)+AZ(4)+AZ(5))
      END
!     ===============================================================================
!     SUBROUTINE MORISON ( RELATIVE MOTION ) 
!     =============================================================================== 
      SUBROUTINE MORISON (RATIO,RK,HIGHT,Y,X,ZETA,TIME,DEPTH,DF,ROLL,D,G,CD,FF,ZP,Pi,CM)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION F(5),FF(5),DX(3),GG(5),XX(3),WF(3),DDX(5),DXX(3),DA(5),DAT(3)
     1           ,D1(3),D2(3),D3(3),D4(3),D5(3),DDXT(3),AU1(3),AU2(3),AU3(3)
     1           ,AU4(3),AU5(3),AV1(3),AV2(3),AV3(3),AV4(3),AV5(3),R1(3),R2(3)
     1           ,R3(3),R4(3),R5(3),S1(3),S2(3),S3(3),S4(3),S5(3),AX1(3),AX2(3)
     1           ,AX3(3),AX4(3),AX5(3),AY1(3),AY2(3),AY3(3),AY4(3),AY5(3),AXT(3)
     1           ,AYT(3),AXIT(3),AYIT(3),ALX(3),AUX(3),ALY(3),AUY(3),AZ1(3),AZ2(3)
     1           ,AZ3(3),AZ4(3),AZ5(3),AZIT(3),ALZ(3),AUZ(3),AZT(3)
      CALL INTERPOLATION_G (G11,G13,G15,G22,G24,G33,G35,G44,G55,RATIO)
      CALL INTERPOLATION_C (C1,C2,C3,C4,RATIO)
      CALL WAVE_HIGHT (RK,HIGHT,A,RATIO)
      WS=(G*RK)
      WSS=(1+((A**2)*C1)+((A**4)*C2))
      WSSS=TANH(RK*DEPTH)
      W=SQRT(WS*WSS*WSSS)
      C=SQRT((G/RK)*(1.0+(A**2)*C1+(A**4)*C2)*TANH(RK*DEPTH))
      GG(1)=(A*G11)+((A**3)*G13)+((A**5)*G15)
      GG(2)=2.0*(((A**2)*G22)+((A**4)*G24))
      GG(3)=3.0*((A**3)*G33)+((A**5)*G35)
      GG(4)=4.0*(A**4)*G44
      GG(5)=5.0*A**5*G55
!     GENERATE DRAG FORCE USING GAUSS QUADRATURE INTERGRATION 
!     DRAG FORCE X AXIAL  -----------------------------------------------------------------
!     XX = SAMPLING POINT LOCATION IN 3 NODE 
      XX(1)=-SQRT(0.6)
      XX(2)=0D0
      XX(3)=SQRT(0.6)
!     WF = WEIGHT FACTORS
      WF(1)=5.0/9.0
      WF(2)=8.0/9.0
      WF(3)=5.0/9.0
!     Y1 = LOWER FORM COORDINATE
!     Y2 = UPPER FORM COORDINATE
      Y1=0
      Y2=97.2
!     MAPPING
      A1=(Y2-Y1)/2
      A2=(Y2+Y1)/2
      DO I=1,3
      Y=A1*XX(I)+A2
!     CALCULATE DRAG FORCR VELOCITY HORIZONTAL DIRECTION AT ANY HIGHT(Y) 
      D1(I)=(W/RK)*GG(1)*(COS(ZETA))*(COSH(RK*Y)/SINH(RK*DEPTH))*COS((RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D2(I)=(W/RK)*GG(2)*(COS(ZETA))*(COSH(2*RK*Y)/SINH(2*RK*DEPTH))*COS(2*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D3(I)=(W/RK)*GG(3)*(COS(ZETA))*(COSH(3*RK*Y)/SINH(3*RK*DEPTH))*COS(3*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D4(I)=(W/RK)*GG(4)*(COS(ZETA))*(COSH(4*RK*Y)/SINH(4*RK*DEPTH))*COS(4*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D5(I)=(W/RK)*GG(5)*(COS(ZETA))*(COSH(5*RK*Y)/SINH(5*RK*DEPTH))*COS(5*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
!     DAT = TOTAL VELOCITY IN HORIZONTAL DIRECTION AT ANY HIGHT(Y) 
      DAT(I)=(D1(I)+D2(I)+D3(I)+D3(I)+D4(I)+D5(I))
      DDXT(I)=0.5*ROLL*CD*D*(DAT(I)*ABS(DAT(I)))*A1
!     DX  = FORCE AT LOWER NODE AT ANY HIGHT(Y) 
!     DXX = FORCE AT UPPER NODE AT ANY HIGHT(Y) 
      DX(I)=(1-(Y/(Y2-Y1)))*DDXT(I)*(WF(I))
      DXX(I)=(Y/(Y2-Y1))*DDXT(I)*(WF(I))
      ENDDO
!     DTYU= TOTAL FORCE UPPER NODE AT ANY HIGHT(Y)
!     DFYL= TOTAL FORCE AT LOWER NODE AT ANY HIGHT(Y)
      TRY=DDXT(1)*WF(1)+DDXT(2)*WF(2)+DDXT(3)*WF(3) 
      DTYU=DXX(1)+DXX(2)+DXX(3)
      DTYL=DX(1)+DX(2)+DX(3)
!     -------------------------------------------------------------------------------------
!     DRAG FORCE Y AXIAL  -----------------------------------------------------------------
      DO I=1,3
      Y=A1*XX(I)+A2
!     CALCULATE DRAG FORCR VELOCITY VERTICAL DIRECTION AT ANY HIGHT(Y) 
      D1(I)=(W/RK)*GG(1)*(SINH(RK*Y)/SINH(RK*DEPTH))*SIN((RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D2(I)=(W/RK)*GG(2)*(SINH(2*RK*Y)/SINH(2*RK*DEPTH))*SIN(2*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D3(I)=(W/RK)*GG(3)*(SINH(3*RK*Y)/SINH(3*RK*DEPTH))*SIN(3*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D4(I)=(W/RK)*GG(4)*(SINH(4*RK*Y)/SINH(4*RK*DEPTH))*SIN(4*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D5(I)=(W/RK)*GG(5)*(SINH(5*RK*Y)/SINH(5*RK*DEPTH))*SIN(5*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
!     DAT = TOTAL VELOCITY IN VERTICAL DIRECTION AT ANY HIGHT(Y) 
      DAT(I)=(D1(I)+D2(I)+D3(I)+D3(I)+D4(I)+D5(I))
      DDXT(I)=0.5*ROLL*CD*D*DAT(I)*ABS(DAT(I))*A1
!     DX  = FORCE AT LOWER NODE AT ANY HIGHT(Y) 
!     DXX = FORCE AT UPPER NODE AT ANY HIGHT(Y) 
      DX(I)=(1-(Y/(Y2-Y1)))*DDXT(I)*(WF(I))
      DXX(I)=(Y/(Y2-Y1))*DDXT(I)*(WF(I))
      ENDDO
!     DTYU= TOTAL FORCE UPPER NODE AT ANY HIGHT(Y)
!     DFYL= TOTAL FORCE AT LOWER NODE AT ANY HIGHT(Y) 
      DTXU=DXX(1)+DXX(2)+DXX(3)
      DTXL=DX(1)+DX(2)+DX(3)
!     -------------------------------------------------------------------------------------
!     DRAG FORCE Z AXIAL  -----------------------------------------------------------------   
      DO I=1,3
      Y=A1*XX(I)+A2
!     CALCULATE DRAG FORCR VELOCITY HORIZONTAL DIRECTION AT ANY HIGHT(Y) 
      D1(I)=(W/RK)*GG(1)*(SIN(ZETA))*(COSH(RK*Y)/SINH(RK*DEPTH))*COS((RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D2(I)=(W/RK)*GG(2)*(SIN(ZETA))*(COSH(2*RK*Y)/SINH(2*RK*DEPTH))*COS(2*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D3(I)=(W/RK)*GG(3)*(SIN(ZETA))*(COSH(3*RK*Y)/SINH(3*RK*DEPTH))*COS(3*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D4(I)=(W/RK)*GG(4)*(SIN(ZETA))*(COSH(4*RK*Y)/SINH(4*RK*DEPTH))*COS(4*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      D5(I)=(W/RK)*GG(5)*(SIN(ZETA))*(COSH(5*RK*Y)/SINH(5*RK*DEPTH))*COS(5*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
!     DAT = TOTAL VELOCITY IN HORIZONTAL DIRECTION AT ANY HIGHT(Y) 
      DAT(I)=(D1(I)+D2(I)+D3(I)+D3(I)+D4(I)+D5(I))
      DDXT(I)=0.5*ROLL*CD*D*DAT(I)*ABS(DAT(I))*A1
!     DX  = FORCE AT LOWER NODE AT ANY HIGHT(Y) 
!     DXX = FORCE AT UPPER NODE AT ANY HIGHT(Y) 
      DX(I)=(1-(Y/(Y2-Y1)))*DDXT(I)*(WF(I))
      DXX(I)=(Y/(Y2-Y1))*DDXT(I)*(WF(I))
      ENDDO
!     DTYU= TOTAL FORCE UPPER NODE AT ANY HIGHT(Y)
!     DFYL= TOTAL FORCE AT LOWER NODE AT ANY HIGHT(Y) 
      DTZU=DXX(1)+DXX(2)+DXX(3)
      DTZL=DX(1)+DX(2)+DX(3)
!     FOR VELOCITY COEFFICIENTS -------------------------------------------
      DO I=1,3
      Y=A1*XX(I)+A2
      AU1(I)=GG(1)*(COSH(1*RK*Y)/SINH(1*RK*DEPTH))
      AU2(I)=GG(2)*(COSH(2*RK*Y)/SINH(2*RK*DEPTH))
      AU3(I)=GG(3)*(COSH(3*RK*Y)/SINH(3*RK*DEPTH))
      AU4(I)=GG(4)*(COSH(4*RK*Y)/SINH(4*RK*DEPTH))
      AU5(I)=GG(5)*(COSH(5*RK*Y)/SINH(5*RK*DEPTH))
      AV1(I)=GG(1)*(SINH(1*RK*Y)/SINH(1*RK*DEPTH))
      AU2(I)=GG(2)*(COSH(2*RK*Y)/SINH(2*RK*DEPTH))
      AU3(I)=GG(3)*(COSH(3*RK*Y)/SINH(3*RK*DEPTH))
      AU4(I)=GG(4)*(COSH(4*RK*Y)/SINH(4*RK*DEPTH))
      AU5(I)=GG(5)*(COSH(5*RK*Y)/SINH(5*RK*DEPTH))
      ENDDO
!     CALCULATE R(N) S(N) -------------------------------------------------
      DO I=1,3   
      R1(I)=2.0*AU1(I)-AU1(I)*AU2(I)-AV1(I)*AV2(I)-AU2(I)*AU3(I)-AV2(I)*AV3(I)
      R2(I)=4.0*AU2(I)-AU1(I)**2+AV1(I)**2-2.0*AU1(I)*AU3(I)-2.0*AV1(I)*AV3(I)
      R3(I)=6.0*AU3(I)-3.0*AU1(I)*AU2(I)+3.0*AV1(I)*AV2(I)-3.0*AU1(I)*AU4(I)-3.0*AV1(I)*AV4(I)
      R4(I)=8.0*AU4(I)-2.0*AU2(I)**2+2.0*AV2(I)**2-4.0*AU1(I)*AU3(I)+4.0*AV1(I)*AV3(I)
      R5(I)=10.0*AU5(I)-5.0*AU1(I)*AU4(I)-5.0*AU2(I)*AU3(I)+5.0*AV1(I)*AV4(I)+5.0*AV2(I)*AV3(I)
      S1(I)=-2.0*AU1(I)*AV1(I)+4*AU2(I)*AV1(1)
      S2(I)=2.0*AV1(I)-3.0*AU1(I)*AV2(I)-3.0*AU2(I)*AV1(I)-5.0*AU2(I)*AV3(I)-5.0*AU5(I)*AV2(I)
      S3(I)=6.0*AV3(I)-AU1(I)*AV2(I)+AU2(I)*AV1(I)-5.0*AU1(I)*AV4(I)-5.0*AU4(I)*AV1(I)
      S4(I)=8.0*AV4(I)-2.0*AU1(I)*AV3(I)+2.0*AU3(I)*AV1(I)+4.0*AU2(I)*AV2(I)
      S5(I)=10.0*AV5(I)-3.0*AU1(I)*AV4(I)+3.0*AU4(I)*AV1(I)-AU2(I)*AV3(I)+AU3(I)*AV2(I)
      ENDDO
!     CALCULATE HORIZONTAL VERTICAL DIRECTION AT ANY HIGHT(Y)--------- 
      DO I=1,3
      Y=A1*XX(I)+A2
      AX1(I)=(1.0/2.0)*COS(ZETA)*(RK*(C**2))*R1(I)*SIN(1*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AX2(I)=(1.0/2.0)*COS(ZETA)*(RK*(C**2))*R2(I)*SIN(2*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AX3(I)=(1.0/2.0)*COS(ZETA)*(RK*(C**2))*R3(I)*SIN(3*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AX4(I)=(1.0/2.0)*COS(ZETA)*(RK*(C**2))*R4(I)*SIN(4*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AX5(I)=(1.0/2.0)*COS(ZETA)*(RK*(C**2))*R5(I)*SIN(5*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
!     CALCULATE VELOCITY VERTICAL DIRECTION AT ANY HIGHT(Y)----------- 
      AY1(I)=(-1.0/2.0)*(RK*(C**2))*S1(I)*COS(1*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AY2(I)=(-1.0/2.0)*(RK*(C**2))*S2(I)*COS(2*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AY3(I)=(-1.0/2.0)*(RK*(C**2))*S3(I)*COS(3*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AY4(I)=(-1.0/2.0)*(RK*(C**2))*S4(I)*COS(4*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AY5(I)=(-1.0/2.0)*(RK*(C**2))*S5(I)*COS(5*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
!     CALCULATE HORIZONTAL VERTICAL DIRECTION AT ANY HIGHT(Y)---------  
      AZ1(I)=(1.0/2.0)*SIN(ZETA)*(RK*(C**2))*R1(I)*SIN(1*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AZ2(I)=(1.0/2.0)*SIN(ZETA)*(RK*(C**2))*R2(I)*SIN(2*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AZ3(I)=(1.0/2.0)*SIN(ZETA)*(RK*(C**2))*R3(I)*SIN(3*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AZ4(I)=(1.0/2.0)*SIN(ZETA)*(RK*(C**2))*R4(I)*SIN(4*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME))
      AZ5(I)=(1.0/2.0)*SIN(ZETA)*(RK*(C**2))*R5(I)*SIN(5*(RK*ZP*SIN(ZETA)+RK*X*COS(ZETA)-W*TIME)) 
!     AXT = TOTAL ACCELERATION IN HORIZONTAL DIRECTION AT ANY HIGHT(Y)   
      AXT(I)=AX1(I)+AX2(I)+AX3(I)+AX4(I)+AX5(I)
!     AYT = TOTAL VERTICAL IN HORIZONTAL DIRECTION AT ANY HIGHT(Y)     
      AYT(I)=AY1(I)+AY2(I)+AY3(I)+AY4(I)+AY5(I)
!     AZT = TOTAL ACCELERATION IN HORIZONTAL DIRECTION AT ANY HIGHT(Y)
      AZT(I)=AZ1(I)+AZ2(I)+AZ3(I)+AZ4(I)+AZ5(I)
!     INERIA FORCR AT ANY HIGHT(Y) --------- 
      AXIT(I)=(ROLL*CM*Pi*(D**2)/4)*AXT(I)*A1
      AYIT(I)=(ROLL*CM*Pi*(D**2)/4)*AYT(I)*A1
      AZIT(I)=(ROLL*CM*Pi*(D**2)/4)*AZT(I)*A1
!     ALX  = FORCE AT LOWER NODE AT ANY HIGHT(Y) 
!     AUX = FORCE AT UPPER NODE AT ANY HIGHT(Y) 
      ALX(I)=(1-(Y/(Y2-Y1)))*AXIT(I)*(WF(I))
      AUX(I)=(Y/(Y2-Y1))*AXIT(I)*(WF(I))
!     ALX  = FORCE AT LOWER NODE AT ANY HIGHT(Y) 
!     AUX = FORCE AT UPPER NODE AT ANY HIGHT(Y) 
      ALY(I)=(1-(Y/(Y2-Y1)))*AYIT(I)*(WF(I))
      AUY(I)=(Y/(Y2-Y1))*AYIT(I)*(WF(I))
!     ALZ  = FORCE AT LOWER NODE AT ANY HIGHT(Y) 
!     AUZ = FORCE AT UPPER NODE AT ANY HIGHT(Y) 
      ALZ(I)=(1-(Y/(Y2-Y1)))*AZIT(I)*(WF(I))
      AUZ(I)=(Y/(Y2-Y1))*AZIT(I)*(WF(I))
      ENDDO
!     TIXFU= TOTAL FORCE UPPER NODE AT ANY HIGHT(Y)
!     TIXFL TOTAL FORCE AT LOWER NODE AT ANY HIGHT(Y) 
      TIXFU=ALX(1)+ALX(2)+ALX(3)
      TIXFL=AUX(1)+AUX(2)+AUX(3)
!     TIYFU= TOTAL FORCE UPPER NODE AT ANY HIGHT(Y)
!     TIYFL TOTAL FORCE AT LOWER NODE AT ANY HIGHT(Y)
      TIYFU=ALY(1)+ALY(2)+ALY(3)
      TIYFL=AUY(1)+AUY(2)+AUY(3)
!     TIXFU= TOTAL FORCE UPPER NODE AT ANY HIGHT(Y)
!     TIXFL TOTAL FORCE AT LOWER NODE AT ANY HIGHT(Y) 
      TIZFU=ALZ(1)+ALZ(2)+ALZ(3)
      TIZFL=AUZ(1)+AUZ(2)+AUZ(3)  
      END
!     ===============================================================================
!     SUBROUTINE COEFFICIENTS ( DATA COEFFICIENTS ) 
!     ( DNV-OS-J101 JUNE 2004 PAGE 42 - SEC.4 )
!     ===============================================================================
      SUBROUTINE COEFFICIENTS (BB,ANS1,ANS2)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(95,3)
!    COEFFICIENTS FROM DNV-OS-J101 JUNE 2004 (TABLE E5)========================================================= 
      DATA A/0.02D0,0.04D0,0.06D0,0.08D0,0.1D0,0.12D0,0.14D0,0.16D0,0.18D0,0.2D0,0.22D0,0.24D0,0.26D0,0.28D0,0.3D0
     1     ,0.32D0,0.34D0,0.36D0,0.38D0,0.4D0,0.42D0,0.44D0,0.46D0,0.48D0,0.5D0,0.52D0,0.54D0,0.56D0,0.58D0,0.6D0
     1      ,0.62D0,0.64D0,0.66D0,0.68D0,0.7D0,0.72D0,0.74D0,0.76D0,0.78D0,0.8D0,0.82D0,0.84D0,0.86D0,0.88D0,0.9D0
     1      ,0.92D0,0.94D0,0.96D0,0.98D0,1D0,1.2D0,1.4D0,1.6D0,1.8D0,2D0,2.2D0,2.4D0,2.6D0,2.8D0,3D0,3.2D0,3.4D0
     1      ,3.6D0,3.8D0,4D0,4.2D0,4.4D0,4.6D0,4.8D0,5D0,5.2D0,5.4D0,5.6D0,5.8D0,6D0,6.2D0,6.4D0,6.6D0,6.8D0,7D0
     1      ,7.2D0,7.4D0,7.6D0,7.8D0,8D0,8.2D0,8.4D0,8.6D0,8.8D0,9D0,9.2D0,9.4D0,9.6D0,9.8D0,10D0
     1      ,0.00063D0,0.00252D0,0.00568D0,0.01012D0,0.01586D0,0.0229D0,0.03126D0,0.04095D0,0.05197D0,0.06433D0
     1      ,0.07802D0,0.09304D0,0.10937D0,0.12701D0,0.14591D0,0.16606D0,0.1874D0,0.20989D0,0.23347D0,0.25808D0
     1      ,0.28364D0,0.31008D0,0.33732D0,0.36526D0,0.39381D0,0.42287D0,0.45234D0,0.48214D0,0.51216D0,0.5423D0
     1      ,0.57249D0,0.60262D0,0.63236D0,0.66244D0,0.69197D0,0.72118D0,0.75D0,0.77839D0,0.80631D0,0.83373D0
     1      ,0.86062D0,0.88697D0,0.91276D0,0.93797D0,0.96261D0,0.98667D0,1.01016D0,1.03308D0,1.05545D0,1.07726D0
     1      ,1.26842D0,1.42148D0,1.54958D0,1.66133D0,1.76191D0,1.85448D0,1.94099D0,2.02275D0,2.10063D0,2.17524D0
     1      ,2.24705D0,2.31641D0,2.3836D0,2.44882D0,2.51228D0,2.57411D0,2.63444D0,2.6934D0,2.75107D0,2.80754D0
     1      ,2.86288D0,2.91717D0,2.97045D0,3.0228D0,3.07425D0,3.12485D0,3.17465D0,3.22368D0,3.27197D0,3.31956D0
     1      ,3.36648D0,3.41275D0,3.45841D0,3.50348D0,3.54797D0,3.59192D0,3.63533D0,3.67824D0,3.72064D0,3.76258D0
     1      ,3.80405D0,3.84508D0,3.88567D0,3.92585D0,3.96562D0,0.018D0,0.072D0,0.162D0,0.289D0,0.453D0,0.653D0
     1      ,0.889D0,1.162D0,1.471D0,1.816D0,2.195D0,2.609D0,3.056D0,3.534D0,4.043D0,4.581D0,5.145D0,5.733D0
     1      ,6.343D0,6.972D0,7.617D0,8.276D0,8.944D0,9.619D0,10.298D0,10.976D0,11.65D0,12.318D0,12.975D0,13.618D0
     1      ,14.245D0,14.852D0,15.438D0,15.998D0,16.532D0,17.036D0,17.51D0,17.952D0,18.36D0,18.734D0,19.073D0
     1      ,19.376D0,19.643D0,19.874D0,20.068D0,20.226D0,20.349D0,20.435D0,20.487D0,20.504D0,18.94D0,14.804D0
     1      ,8.844D0,1.611D0,-6.522D0,-15.306D0,-24.572D0,-34.201D0,-44.112D0,-54.245D0,-64.555D0,-75.009D0,-85.581D0
     1      ,-96.253D0,-107.007D0,-117.832D0,-128.717D0,-139.654D0,-150.637D0,-161.659D0,-172.716D0,176.196D0
     1      ,165.081D0,153.941D0,142.779D0,131.598D0,120.399D0,109.183D0,97.953D0,86.71D0,75.454D0,64.188D0,52.91D0
     1      ,41.624D0,30.328D0,19.025D0,7.714D0,-3.604D0,-14.929D0,-26.26D0,-37.596D0,-48.938D0,-60.284D0,-71.635D0
     1      ,-82.991D0/
!     =========================================================================================================
!     FORM DATA THEN INTERPOLATION COEFFICIENTS AT ANY KR MORISON EQUATION IS NOT VALID D > 0.2*RAMDA
        DO I=1,94
          IF(BB.GE.A(I,1).AND.BB.LE.A(I+1,1))THEN
          X1=A(I,1)
          X2=A(I+1,1)
          Y1=A(I,2)
          Y2=A(I+1,2)
          Z1=A(I,3)
          Z2=A(I+1,3)
          ANS1=(((Y2-Y1)/(X2-X1))*(BB-X1))+Y1
          ANS2=(((Z2-Z1)/(X2-X1))*(BB-X1))+Z1
          ENDIF
        ENDDO
       END
!     ===============================================================================
!     SUBROUTINE NVMORISON ( WHEN D>0.2*RAMDA MORISON EQUATION IS NOT VALID )
!     THE DIFFRACTION THEROY DNV-OS-J101 JUNE 2004 SEC 4 PAGE 41
!     ===============================================================================
      SUBROUTINE NVMORISON (FX,ROLL,AW,RK,DEPTH,G,HF,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      BB=RK*(D/2)
      CALL COEFFICIENTS (BB,ANS1,ANS2)
      XI=ANS1
      ALPHA=ANS2
      A=(4*G*ROLL)/(RK**2)
      AA=(SINH(RK*(DEPTH+AW*SIN(ALPHA))))/(TANH(RK*DEPTH))
      FX=A*AA*XI
      AAA=RK*DEPTH*SINH(RK*DEPTH)-COSH(RK*DEPTH)+1
      AAAA=RK*DEPTH*SINH(RK*DEPTH)
      HF=AAA/AAAA
      END
!     ===============================================================================
!     SUBROUTINE CURRENT ( CURRENT DATA )
!     - THE CRRENT IS REPRESENT BY WIND-GENERATESD CURRENT LOAD
!     ( DNV-OS-J101 JUNE 2004 PAGE 29 - SEC.3 )
!     ===============================================================================
      SUBROUTINE CURRENT (VWIND0,DEPTH,H0,ROLLA,HEIGHT,TF,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
     
      VWIND=VWIND0*((DEPTH+HEIGHT)/(H0))
      FWT=0.5*ROLLA*VWIND*VWIND*D
      
      TF=FWT
      END
!     ===============================================================================
!     SUBROUTINE CURRENT ( CURRENT DATA )
!     - THE CRRENT IS REPRESENT BY WIND-GENERATESD CURRENT ,TIDE CURRENT MAY EXIST
!     ( DNV-OS-J101 JUNE 2004 PAGE 29 - SEC.3 )
!     ===============================================================================
      SUBROUTINE TIDE_LOAD (VTIDE0,DEPTH,ROLL,HEIGHT,THIGHT,TF,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      IF(HEIGHT.GT.THIGHT) THEN
      FTT = 0.0D0
      ELSE
      VTIDE=VTIDE0*(((DEPTH+HEIGHT)/DEPTH)**(1.0/7.0))
      FTT=0.5*ROLL*VTIDE*VTIDE*D
      ENDIF
      
      TF=FTT
      END
!     ===============================================================================
!     SUBROUTINE PLUNGING ( PLUNGING WAVE,AN IMPACT MODEL )
!     - FOR PLUNGING WAVE, AN IMPACT MODEL MAY BE USED TO CALCULATE THE WAVE FORCE ON
!     A STUCYURE. THE IMPACT FORCE FORM PLUNGING WAVE ( DNV-OS-J101 JUNE 2004 PAGE 42 - SEC.4 )
!     ===============================================================================
      SUBROUTINE PLUNGING (CS,PI,ROOL,UT,AP,F2)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      AA=2.0D0*PI
      IF (CS.LE.3.0D0.AND.CS.GT.AA)THEN
!      FOR PLUNGING WAVE FOR SMOOTH CIRCULAR CYCINDER 
 90    WRITE (6,950)
950    FORMAT ('THE SLAMMING COEFFICIENT SHOULD NOT BE TAKEN LESS THAN 3.0 (CS<3.0) THE UPPER LIMIT 2*PI(CD<2*PI)')
      ENDIF
      F2=0.5D0*ROOL*CS*AP*(UT**2)
      END
!     ===============================================================================
!     SUBROUTINE SURING ( SURING AND SPILLING )
!     - ( DNV-OS-J101 JUNE 2004 PAGE 42 - SEC.4 )
!     ===============================================================================
      SUBROUTINE SURING (FSP,HIGHT,DEPTH,ROLL,PI,CM,CD,FF,S,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION FF(5)
!     CALL MORISON FOR GENERATE STORK FORCE BY ORDER --------------------------------
!     GET (FF) BY EVERY ORDER THEN ----------------------------------------------------
      CS=5.15*((D/(D+19*S))+(0.107*S/D))
      IF (DEPTH.NE.HIGHT)THEN
      FSP=0.5*ROLL*CS*D*(FF(1)+FF(2)+FF(3)+FF(4)+FF(5))
      ENDIF
      IF (DEPTH.EQ.HIGHT)THEN
      FSP=0.5*CD*ROLL*D*(FF(1)+FF(2)+FF(3)+FF(4)+FF(5)) 
      ENDIF 
      END
!     ===============================================================================
!     SUBROUTINE VWIND ( WIND VELOCITY )
!     INPUT  >>> HIGH,Z,Z0,UH,ALPHA,TR,ROLLA,HEIGHT
!     OUTPUT >>> FUM( THE LONG TERM DISTRIBUTION OF 10 MINUTE MEAN WIND SPEED ),TR(RETURN PERIOD)
!     ===============================================================================
      SUBROUTINE VWIND (HIGH,Z0,UH,ALPHA,FUM,UR,ROLLA,HEIGHT,DEPTH,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION UZ(3),FUT(3),FUB(3),FUA(3),XX(3),WF(3)
!     THE SCALE PARAMETER A(HEIGHT)
!     HEIGHT   = DISTANCE FROM STILL WATER LEVEL,POSITIVE UPWARDS
!     Z0       = THE ROUGHNESS LENGTH OF THE SURFACE
!     HIGH     = THE REFERNCE HEIGHT
      UZT=UH*((HEIGHT/HIGH)**ALPHA)
      FUM=0.5*ROLLA*D*UZT*UZT
      END
!     ==============================================================================
!     SUBROUTINE VELOCITY ( HORIZONTAL AND VERTICAL VELOCITY,VELOCITY COEFFICIENTS ) 
!     ==============================================================================
      SUBROUTINE PEAKLEV (DEPTH,RAMDA,G,TIME,ZETA,HIGHT,TY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION F(5),FU(5)
      
       Pi=3.141592654
       RK=2*Pi/RAMDA
       RATIO=DEPTH/RAMDA
       MORRI=0.2*RAMDA
          
      CALL INTERPOLATION_F (F22,F24,F33,F35,F44,F55,RATIO)
      CALL INTERPOLATION_G (G11,G13,G15,G22,G24,G33,G35,G44,G55,RATIO)
      CALL INTERPOLATION_C (C1,C2,C3,C4,RATIO)
      CALL WAVE_HIGHT (RK,HIGHT,A,RATIO)
      WS=(G*RK)
      WSS=(1.0+((A**2.0)*C1)+((A**4)*C2))
      WSSS=TANH(RK*DEPTH)
      W=SQRT(WS*WSS*WSSS)
      G1=(A*G11)+((A**3)*G13)+((A**5)*G15)
      G2=2.0*(((A**2)*G22)+((A**4)*G24))
      G3=3.0*((A**3)*G33)+((A**5)*G35)
      G4=4.0*(A**4)*G44
      G5=5.0*A**5*G55
!     THE FREE-SURFACE WATER DEFLECTION ---------------------------------------------
      F(1)=A
      F(2)=((A**2)*F22)+((A**4)*F24)
      F(3)=((A**3)*F33)+((A**5)*F35)
      F(4)=((A**4)*F44)
      F(5)=((A**5)*F55)
      AA=((RK*Z*COS(ZETA))+(RK*X*COS(ZETA))-(W*TIME))
      DO I=1,5
      FU(I)=(1/RK)*F(I)*COS(I*AA)
      ENDDO
!     DEFLECTION
      TFNU=FU(1)+FU(2)+FU(3)+FU(4)+FU(5)
      TY=DEPTH+TFNU
      END
!     ===============================================================================
!     SUBROUTINE ACCELERATION( HORIZONTAL AND VERTICAL ACCELERATION,VELOCITY COEFFICIENTS ) 
!     ===============================================================================
   
