C     =========================================================================        
       SUBROUTINE  BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                               WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                               VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                               VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,NBI,STOPFUNCTION,
     1                               NGROWTH,GROWTH)
       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*4 (I-N)
       DIMENSION LWCASE(7)
       DIMENSION VCURRENTL(5)
       
       COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER
       
       IF (NBI.EQ.1.0)THEN
      ! ----------------------------------------------- WARNING MESSAGE --------------------------------------------------
      IF (LWCASE(1).EQ.0.AND.LWCASE(2).EQ.0.AND.LWCASE(3).EQ.0.AND.LWCASE(4).EQ.0.AND.LWCASE(5).EQ.0.AND.LWCASE(6).EQ.0
     1    .AND.LWCASE(7).EQ.0.0D0.AND.LWCASE(8).EQ.0.0D0)THEN
      WRITE (*,8)
      WRITE (*,34)
      WRITE (*,35)
34    FORMAT ('---------------- OFFSHORE ERROR MESSAGE ( FRAME ELEMENT )----------------')
35    FORMAT ('- PLEASE CHECK OFFSHORE LOAD TYPES')
      STOP
      ENDIF
      
      !  --------- ERROR MESSAGE FOR AUTOMATIC GENERATE DRAG ANS INERTIA COEFFICIENTS ---------
      IF (NFUNCTION.EQ.3)THEN
         IF (CD.LE.0.0D0.OR.CM.LE.0.0D0.OR.DIAM2.LE.0.0D0)THEN
         WRITE (*,8)
         WRITE (*,36)
36       FORMAT ('---------------- OFFSHORE ERROR MESSAGE ( FRAME ELEMENT )----------------')
            IF (CD.LE.0.0D0)THEN
            WRITE (*,37)
37          FORMAT ('- DRAG COEFFICIENT MUST BE GRATER THAN ZERO')
            ENDIF
            IF (CM.LE.0.0D0)THEN
            WRITE (*,38)
38          FORMAT ('- INERTIA COEFFICIENT MUST BE GRATER THAN ZERO')
            ENDIF
             IF (NDI.EQ.2.0D0)THEN
                IF (DIAM2.LE.0.0D0)THEN
                WRITE (*,39)
39              FORMAT ('- NORMINAL DIMAETER MUST BE GRATER THAN ZERO')
                ENDIF
             ENDIF
          STOPFUNCTION=1D0    
         ENDIF
      ENDIF
      ! ----------------------------------------------------------------------------------------
      
       ELSEIF (NBI.EQ.2.0D0)THEN
       ! --- ERROR MESSAGE  --- 20 / JULY / 2012
       ! -------------------------------------- WAVE PARAMETER ERROR MESSAGE ------------------------------------ 
       IF (LWCASE(1).EQ.1.0D0)THEN ! WAVE LOAD
           
           IF ( KFTTD .NE.1 .AND. KFATJ.EQ.0 )THEN
               
          IF (WVHIGHT.EQ.0.0D0.OR.WDEPTH.EQ.0.0D0.OR.PERIOD.EQ.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,4) OFFSELECT
4         FORMAT ('----- OFFSHORE ERROR MESSAGE FOR FRAME ( WAVE PARAMETER )-- LOAD PARAMETER',F4.0)

             IF (WVHIGHT.EQ.0.0D0)THEN
             WRITE (*,1)
1            FORMAT ('- WAVE HEIGHT MUST BE GREATER THAN ZERO')
             STOPFUNCTION=1D0
             ENDIF
             
             IF (KFATJ.EQ.0)THEN
             IF (PERIOD.EQ.0.0D0)THEN
             WRITE (*,3)
3            FORMAT ('- PERIOD MUST BE GREATER THAN ZERO')
             STOPFUNCTION=1D0
             ENDIF
             ENDIF
             
             IF (WDEPTH.EQ.0.0D0)THEN
             WRITE (*,2)
2            FORMAT ('- WATER DEPTH MUST BE GREATER THAN ZERO')
             STOPFUNCTION=1D0
             ENDIF
             
          ENDIF
          
      
          
          ELSEIF ( KFTTD .EQ.1 )THEN
               
          IF (WDEPTH.EQ.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,4) OFFSELECT
!4         FORMAT ('----- OFFSHORE ERROR MESSAGE FOR FRAME ( WAVE PARAMETER )-- LOAD PARAMETER',F4.0)          
          IF (WDEPTH.EQ.0.0D0)THEN
          WRITE (*,2)
!2            FORMAT ('- WATER DEPTH MUST BE GREATER THAN ZERO')          
          STOPFUNCTION=1D0
          ENDIF
          ENDIF
          ENDIF
      ENDIF
      
8     FORMAT ('') ! BLANK

      ! -------------------------------------- CURRENT PARAMETER ERROR MESSAGE ------------------------------------ 
      IF (LWCASE(2).EQ.1.0D0)THEN ! CURRENT LOAD
      ! NCURRENT = 1 ; DNV
      ! NCURRENT = 2 ; API
      ! NCURRENT = 3 ; IEC
      ! NCURRENT = 4 ; POWER LAW
      ! NCURRENT = 5 ; LINEAR CURRENT 5 VALUE
      
         IF (NCURRENT.EQ.1.OR.NCURRENT.EQ.3)THEN
            IF (VTIDE.EQ.0.0D0.OR.VWIND0.EQ.0.0D0)THEN
               WRITE (*,8)
               WRITE (*,7) OFFSELECT
7              FORMAT ('----- OFFSHORE ERROR MESSAGE FOR FRAME ( CURRENT PAMETERS FOR FRAME )-- LOAD PARAMETER',F4.0)
               IF (VWIND0.EQ.0.0D0)THEN
               WRITE (*,5)
5              FORMAT ('- WIND-GENERATE CURRENT VELOCITY ( DNV,IEC ) MUST BE GREATER THAN ZERO')
               ENDIF
               IF (VTIDE.EQ.0.0D0)THEN
               WRITE (*,6)
6              FORMAT ('- WIND-GENERATE TIDAL VELOCITY ( DNV,IEC ) MUST BE GREATER THAN ZERO')
               ENDIF
               STOPFUNCTION=1D0
            ENDIF
            
          ELSEIF (NCURRENT.EQ.2)THEN
             IF (MAXVAL(VCURRENTL).EQ.0.0D0)THEN
                 WRITE (*,8)
                 WRITE (*,25) OFFSELECT
25               FORMAT ('----- OFFSHORE ERROR MESSAGE FOR FRAME ( CURRENT PAMETERS )-- LOAD PARAMETER',F4.0)
                 WRITE (*,26)
26               FORMAT ('- CURRENT VELOCITY (API) MUST BE GREATER THAN ZERO')
             STOPFUNCTION=1D0
             ENDIF
             
          ELSEIF (NCURRENT.EQ.4)THEN
             IF (POWERLAW.EQ.0.0D0.OR.VCURRENTP.EQ.0.0D0)THEN
                 WRITE (*,8)
                 WRITE (*,27) OFFSELECT
27               FORMAT ('----- OFFSHORE ERROR MESSAGE FOR FRAME ( CURRENT PAMETERS )-- LOAD PARAMETER',F4.0)
                 WRITE (*,28)
28               FORMAT ('- CURRENT VELOCITY ( POWER LAW ) MUST BE GREATER THAN ZERO')
             STOPFUNCTION=1D0
             ENDIF
          ELSEIF (NCURRENT.EQ.5)THEN
             IF (VCURRENTL(1).EQ.0.0D0)THEN
                WRITE (*,8)
                WRITE (*,29) OFFSELECT
29              FORMAT ('----- OFFSHORE ERROR MESSAGE FOR FRAME ( CURRENT PAMETERS )-- LOAD PARAMETER',F4.0)
                WRITE (*,40)
40              FORMAT ('- CURRENT VELOCITY ( LINEAR CURRENT ) MUST BE GREATER THAN ZERO')
             STOPFUNCTION=1D0
             ENDIF
          ENDIF
      ENDIF
      
      ! -------------------------------------- WIND PARAMETER ERROR MESSAGE ------------------------------------ 
      IF (LWCASE(7).EQ.1.0D0)THEN ! WIND LOAD
          ! NWINDO = 1 ; LOGARITHMIC PROFILE   (DNV) 
          ! NWINDO = 2 ; POWER LAW PROFILE     (DNV)
          ! NWINDO = 3 ; WIND PROFILE AND GUST (API)
          ! NWINDO = 4 ; IEC 614000-1
          ! NWINDO = 5 ; USER DEFINED
       IF(NWINDO.EQ.1.OR.NWINDO.EQ.2)THEN
         IF (ALPHA.EQ.0.0D0.OR.UH.EQ.0.0D0.OR.RHIGH.EQ.0.0D0)THEN
         WRITE (*,8)
         WRITE (*,9) OFFSELECT
9        FORMAT ('----- OFFSHORE ERROR MESSAGE ( WIND LOAD FOR FRAME )-- LOAD PARAMETER',F4.0)
            IF(UH.EQ.0.0D0)THEN !
            WRITE (*,10)
10          FORMAT ('- MEAN WIND SPEED (DNV) MUST BE GREATER THAN ZERO')
            ENDIF
            IF (ALPHA.EQ.0.0)THEN
            WRITE (*,11)
11          FORMAT ('- EXPONENT IN POWER-LAW MODEL FOR WIND SPEED PROFILE MUST BE GREATER THAN ZERO')
            ENDIF
            IF (RHIGH.EQ.0.0D0)THEN
            WRITE (*,15)
15          FORMAT ('- REFERENCE HIGHT MUST BE GREATER THAN ZERO')  
            ENDIF
         STOPFUNCTION=1D0
         ENDIF
        ELSEIF (NWINDO.EQ.3)THEN
           IF (UHAPI.EQ.0.0D0.OR.AVERAGE.EQ.0.0D0.OR.RHIGH.EQ.0.0D0)THEN
           WRITE (*,8)
           WRITE (*,16) OFFSELECT
16         FORMAT ('----- OFFSHORE ERROR MESSAGE ( WIND LOAD FOR FRAME )-- LOAD PARAMETER',F4.0)
           IF (UHAPI.EQ.0.0D0)THEN
           WRITE (*,17)
17         FORMAT ('- MEAN WIND SPEED (API) MUST BE GREATER THAN ZERO')             
           ENDIF
           IF (AVERAGE.EQ.0.0D0)THEN
           WRITE (*,18)
18         FORMAT ('- AVERAGING TIME MUST BE GREATER THAN ZERO AND LESS THAN 3600 SEC')             
           ENDIF
           IF (RHIGH.EQ.0.0D0)THEN
           WRITE (*,19)
19         FORMAT ('- REFERENCE HIGHT MUST BE GREATER THAN ZERO')  
           ENDIF
           STOPFUNCTION=1D0
           ENDIF
        ELSEIF (NWINDO.EQ.4)THEN
         IF (ALPHA.EQ.0.0D0.OR.UH.EQ.0.0D0.OR.RHIGH.EQ.0.0D0)THEN
         WRITE (*,8)
         WRITE (*,20) OFFSELECT
20        FORMAT ('----- OFFSHORE ERROR MESSAGE ( WIND LOAD FOR FRAME )-- LOAD PARAMETER',F4.0)
            IF(UH.EQ.0.0D0)THEN !
            WRITE (*,21)
21          FORMAT ('- MEAN WIND SPEED MUST BE GREATER THAN ZERO')
            ENDIF
            IF (ALPHA.EQ.0.0)THEN
            WRITE (*,21)
22          FORMAT ('- EXPONENT IN POWER-LAW MODEL FOR WIND SPEED PROFILE MUST BE GREATER THAN ZERO')
            ENDIF
            IF (RHIGH.EQ.0.0D0)THEN
            WRITE (*,23)
23          FORMAT ('- REFERENCE HIGHT MUST BE GREATER THAN ZERO')  
            ENDIF   
            STOPFUNCTION=1D0      
         ENDIF
        ELSEIF (NWINDO.EQ.5)THEN ! UN
          IF (ALPHA.EQ.0.0D0.OR.UHD.EQ.0.0D0.OR.RHIGH.EQ.0.0D0)THEN
            WRITE (*,8)
            WRITE (*,20) OFFSELECT
            IF(UHD.EQ.0.0D0)THEN
            WRITE (*,21)
            ENDIF
            IF (ALPHA.EQ.0.0)THEN
            WRITE (*,21)
            ENDIF
            IF (RHIGH.EQ.0.0D0)THEN
            WRITE (*,23)
            ENDIF   
            STOPFUNCTION=1D0      
          ENDIF
        ENDIF
          ENDIF
          
      ELSEIF (NDI.EQ.3)THEN ! FOR MARINE GROWTH
           IF (NGROWTH.EQ.0.0D0)THEN
                IF (GROWTH.NE.0.0D0)THEN
                   WRITE (*,59)
                   WRITE (*,60)
                   WRITE (*,61)
                   WRITE (*,59)
                   WRITE (*,62)
                   WRITE (*,63)
                   WRITE (*,64)
59                 FORMAT ('')                   
60                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( FRAME ELEMENT )----------------')
61                 FORMAT ('- INSIDE MARINE GROWTH FUNCTION HAVE SOME VALUE. II WILL EFFECT ON THE RESULT')
62                 FORMAT ('**** COMMENT *****')
63                 FORMAT (' CALCULATE     = SELECT MARINE GROWTH BLOCK')
64                 FORMAT (' NOT CALAULATE = UNSELECT MARINE GROWTH BLOCK')
                  STOPFUNCTION=1D0 
                  RETURN
                ELSEIF (GROWTH.EQ.0.0D0)THEN
                ! CONTINUE 
                ENDIF
             ELSEIF (NGROWTH.NE.0.0D0)THEN
                IF (GROWTH.GT.0.0D0)THEN
                ! CONTINUE
                ELSEIF (GROWTH.LE.0.0D0)THEN
                   WRITE (*,65)
                   WRITE (*,66)
                   WRITE (*,67)
                   WRITE (*,65)
                   WRITE (*,68)
                   WRITE (*,69)
                   WRITE (*,70)
65                 FORMAT ('')
66                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( FRAME ELEMENT )----------------')
67                 FORMAT ('- MARINE GROWTH MUST BE GRATER THAN ZERO ')
68                 FORMAT ('**** COMMENT *****')
69                 FORMAT (' CALCULATE     = SELECT MARINE GROWTH BLOCK')
70                 FORMAT (' NOT CALAULATE = UNSELECT MARINE GROWTH BLOCK')
                   ! STOPFUNCTION USE FOR STOP PROGRAM WHEN ERROR MESSAGE APPEAR
                   STOPFUNCTION=1D0 
                   ! GOTO 101 : OUT READ INPUT FOR STOP THE PROGRAM
                   RETURN
                ENDIF
             ENDIF
      ENDIF
      ! ----------------------------------------------------------------------------------------------
       
          END
C	=======================================================================
