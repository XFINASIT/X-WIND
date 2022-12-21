
C  ========================================================================================================        
      SUBROUTINE SECTIONPROPERTIES (ISET,MGPS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /DISEC/ DIMENPROP(5000,11)
      DIMENSION PROPSEC(10)
      
      READ(ITI,*)
      READ(ITI,*) NUMSECTION
      READ(ITI,*) PROPSEC(1:10)
      DIMENPROP(ISET,1)    = NUMSECTION
      DIMENPROP(ISET,2:11) = PROPSEC(1:10)
      END
C  ========================================================================================================    
      SUBROUTINE DESIGNREAD
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI 
      COMMON /DESIGNREADATA/ N_TOTAL_GROUP,N_TOTAL_ELEMENT,NELEMENT(500000),NAGROUP(500000),N_OPTION_CODE,NSTEEL
     1                       ,NGROUP(1000),NUNBRACED(1000),NBUCKLING(1000),ANOMINALT(1000),ANOMINALS(1000),UNBRACEDT(1000)
     1                       ,UNBRACEDS(1000),AFACTOR(1000),NCM(1000),ALENGTHRING(1000),NFFMEMBER(1000)          
      ! READ HEADER
      READ(ITI,*) 
      READ(ITI,*) N_OPTION_CODE
      READ(ITI,*) N_TOTAL_GROUP,N_TOTAL_ELEMENT
            
      ! READ HEADER
      READ(ITI,*) 
      DO I=1,N_TOTAL_GROUP
      READ(ITI,*) NGROUP(I),NUNBRACED(I),NBUCKLING(I),ANOMINALT(I),ANOMINALS(I),UNBRACEDT(I),UNBRACEDS(I),AFACTOR(I)
     1            ,NCM(I),ALENGTHRING(I),NFFMEMBER(I)
      ENDDO
      
      ! READ HEADER
      READ(ITI,*)
      DO J=1,N_TOTAL_ELEMENT
      READ(ITI,*) NELEMENT(J),NAGROUP(J)
      ENDDO
      
      END
      
C  ========================================================================================================         
      SUBROUTINE DESIGNFACTOR (NE,FACTOR,NS,NOP,NUNBRA,NBUCK,ANOMT,ANOMS,UNBRACET,UNBRACES,NNCM,AALENGTHRING,NFMEMBER)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /DESIGNREADATA/ N_TOTAL_GROUP,N_TOTAL_ELEMENT,NELEMENT(500000),NAGROUP(500000),N_OPTION_CODE,NSTEEL
     1                       ,NGROUP(1000),NUNBRACED(1000),NBUCKLING(1000),ANOMINALT(1000),ANOMINALS(1000),UNBRACEDT(1000)
     1                       ,UNBRACEDS(1000),AFACTOR(1000),NCM(1000),ALENGTHRING(1000),NFFMEMBER(1000)     
      ! ---- INITIAL SETTING ----
        BNOMINAL = 0.0D0
        FACTOR   = 0.0D0
        NS       = 1
        NOP      = N_OPTION_CODE 
      ! -------------------------
      DO I=1,N_TOTAL_ELEMENT
         IF (NE.EQ.NELEMENT(I))THEN
         MGROUP  =  NAGROUP(I)
         NS      =  0.0
            DO J=1,N_TOTAL_GROUP
               IF (MGROUP.EQ.NGROUP(J))THEN
               FACTOR    =  AFACTOR(J)
               NUNBRA    =  NUNBRACED(J)
               NBUCK     =  NBUCKLING(J)
               ANOMT     =  ANOMINALT(J)
               ANOMS     =  ANOMINALS(J)
               UNBRACET  =  UNBRACEDT(J)
               UNBRACES  =  UNBRACEDS(J)
               NS        =  0.0
               NNCM      =  NCM(J)   
               AALENGTHRING = ALENGTHRING(J)
               NFMEMBER   = NFFMEMBER(J)
               RETURN
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      END 

C  ========================================================================================================
      SUBROUTINE CALLENGTH (NELEMENT,VGV,WATERDEPTH,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                      ,ATMIDX,ATMIDY,ATMIDZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
      COMMON /MGRAV/ NGRAV   
      
      DIMENSION VGV(3),XL(100),YL(100),ZL(100)
      DIMENSION ATMIDX(100),ATMIDY(100),ATMIDZ(100)
      
      DIMENSION LMNM(100) !VECTOR STORE ELEMENT CONNECTIVITY NUMBER
      
C     CALLING NUMBER OF NODE NSN AND DIMENSION NSC ... SONGSAK OCT2019      
      CALL LOCATN('@XYZ',KXYZ,NSC,NSN,2) 
      
CC         NODEX=N1(NELEMENT) 
CC         NODEY=N2(NELEMENT) 
      CALL CALLNUMNODE_F (NELEMENT,NODEX,NODEY) !SONGSAK OCT2019
      
      
      ! CALCULATE LENGTH
CC      DO 100 J=1,NN
CC         IF (NND(J).EQ.NODEX)THEN
CC         AX1 = AX(J)
CC         AY1 = AY(J)
CC         AZ1 = AZ(J)
CC         EXIT
CC         ENDIF
CC100   CONTINUE    
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX1,1,NODEX,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY1,2,NODEX,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ1,3,NODEX,0)  !GETTING NODAL COORDINATE
      
      
CC      DO 101 J=1,NN
CC         IF (NND(J).EQ.NODEY)THEN
CC         AX2 = AX(J)
CC         AY2 = AY(J)
CC         AZ2 = AZ(J)
CC         EXIT
CC         ENDIF
CC101   CONTINUE
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019  
	CALL RELFILL('@XYZ',AX2,1,NODEY,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY2,2,NODEY,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ2,3,NODEY,0)  !GETTING NODAL COORDINATE
      
        IF (AX1.GT.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        ELSEIF (AX1.LT.AX2)THEN
        AMAXX = AX2
        AMINX = AX1
        ELSEIF (AX1.EQ.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        ENDIF
        
        IF (AY1.GT.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        ELSEIF (AY1.LT.AY2)THEN
        AMAXY = AY2
        AMINY = AY1
        ELSEIF (AY1.EQ.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        ENDIF
        
        IF (AZ1.GT.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        ELSEIF (AZ1.LT.AZ2)THEN
        AMAXZ = AZ2
        AMINZ = AZ1
        ELSEIF (AZ1.EQ.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        ENDIF
      
      ! OUTPUT ELEMENT LENGTH
        DX = ABS(AMAXX-AMINX)
        DY = ABS(AMAXY-AMINY)
        DZ = ABS(AMAXZ-AMINZ)
        ALENGTH  = SQRT((DX**2)+(DY**2)+(DZ**2))
        
      IF (NOP.EQ.2.0D0)THEN  
      ! OUTPUT THE COORDINATE OF THE CENTER OF GRAVITY OF EACH COMPOSITE PART OF THE BODY
        AMIDX = ((AX1+AX2)/2D0)
        AMIDY = ((AY1+AY2)/2D0)
        AMIDZ = ((AZ1+AZ2)/2D0)
   
      ! -------------------------- OUTPUT FOR BUONCY  --------------------------
      ELSEIF (NOP.EQ.3.0D0)THEN
        NBUX = 0.0
        NBUY = 0.0
        NBUZ = 0.0
        NEX  = 0.0
        NEY  = 0.0
        NEZ  = 0.0

        IF (NGRAV.EQ.1)THEN
        IF (AMAXX.GT.WATERDEPTH.AND.AMINX.LT.WATERDEPTH)THEN ! 
        NEX = 1    ! BETWEEN SUMMERGED AND IN AIR
        ELSEIF (AMAXX.LE.WATERDEPTH.AND.AMINX.LE.WATERDEPTH)THEN ! 
        NEX = 2    ! SUBMERGED
        ELSEIF (AMAXX.GE.WATERDEPTH.AND.AMINX.GE.WATERDEPTH)THEN ! 
        NEX = 3    ! IN AIR
        AX1 = 0.0D0
        AX2 = 0.0D0
        AY1 = 0.0D0
        AY2 = 0.0D0
        AZ1 = 0.0D0
        AZ2 = 0.0D0
        ENDIF
           IF (NEX.EQ.1.0)THEN
              !IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Z
                 IF (AX1.GE.AX2) AX1 = WATERDEPTH
                 IF (AX1.LT.AX2) AX2 = WATERDEPTH
              !  SLOPE 3D LINEAR
                 ! FINDING VECTOR
                 VEC_Y = AMAXY - AMINY
                 VEC_X = AMAXX - AMINX
                 VEC_Z = AMAXZ - AMINZ
                 ! LINEAR 3D (UNKNOWN X,Y)
                 ANEW_Y = (((WATERDEPTH-AMINX)/VEC_X)*VEC_Y) + AMINY
                 ANEW_Z = (((WATERDEPTH-AMINX)/VEC_X)*VEC_Z) + AMINZ
                 AX1    = WATERDEPTH
                 AX2    = AMINX
                 AY1    = ANEW_Y
                 AY2    = AMINY
                 AZ1    = ANEW_Z
                 AZ2    = AMINZ    
              !ENDIF
           ELSEIF (NEX.EQ.2.0)THEN
           ! SUBMERGED CASE
           ELSEIF (NEX.EQ.3.0)THEN
               IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! X
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0   
              ENDIF           
           ENDIF
        ELSEIF (NGRAV.EQ.2)THEN
        IF (AMAXY.GT.WATERDEPTH.AND.AMINY.LT.WATERDEPTH)THEN ! 
        NEY = 1    ! BETWEEN SUMMERGED AND IN AIR
        ELSEIF (AMAXY.LE.WATERDEPTH.AND.AMINY.LE.WATERDEPTH)THEN ! 
        NEY = 2    ! SUBMERGED
        ELSEIF (AMAXY.GE.WATERDEPTH.AND.AMINY.GE.WATERDEPTH)THEN ! 
        NEY = 3    ! IN AIR
        AX1 = 0.0D0
        AX2 = 0.0D0
        AY1 = 0.0D0
        AY2 = 0.0D0
        AZ1 = 0.0D0
        AZ2 = 0.0D0
        ENDIF
           IF (NEY.EQ.1.0)THEN
              !IF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y  
              !  CALCULATE NEW LOCATION OF ELEMENET
                 IF (AY1.GE.AY2) AY1 = WATERDEPTH
                 IF (AY1.LT.AY2) AY2 = WATERDEPTH
              !  SLOPE 3D LINEAR
                 ! FINDING VECTOR
                 VEC_Y = AMAXY - AMINY
                 VEC_X = AMAXX - AMINX
                 VEC_Z = AMAXZ - AMINZ
                 ! LINEAR 3D (UNKNOWN X,Y)
                 ANEW_X = (((WATERDEPTH-AMINY)/VEC_Y)*VEC_X) + AMINX
                 ANEW_Z = (((WATERDEPTH-AMINY)/VEC_Y)*VEC_Z) + AMINZ
                 AX1    = ANEW_X
                 AX2    = AMINX
                 AY1    = WATERDEPTH
                 AY2    = AMINY
                 AZ1    = ANEW_Z
                 AZ2    = AMINZ    
                    !IF (VEC_Z.EQ.0) THEN
                    !AZ1    = AMAXZ
                    !AZ2    = AMINZ
                    !ENDIF
              !ENDIF
           ELSEIF (NEY.EQ.2.0)THEN
           ! SUBMERGED CASE          
           ELSEIF (NEY.EQ.3.0)THEN
               IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! X
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0   
              ENDIF           
           ENDIF
        ELSEIF (NGRAV.EQ.3)THEN
        IF (AMAXZ.GT.WATERDEPTH.AND.AMINZ.LT.WATERDEPTH)THEN ! 
        NEZ = 1    ! BETWEEN SUMMERGED AND IN AIR
        ELSEIF (AMAXZ.LE.WATERDEPTH.AND.AMINZ.LE.WATERDEPTH)THEN ! 
        NEZ = 2    ! SUBMERGED
        ELSEIF (AMAXZ.GE.WATERDEPTH.AND.AMINZ.GE.WATERDEPTH)THEN ! 
        NEZ = 3    ! IN AIR
        AX1 = 0.0D0    ! UPDATE 2018
        AX2 = 0.0D0    ! UPDATE 2018
        AY1 = 0.0D0    ! UPDATE 2018
        AY2 = 0.0D0    ! UPDATE 2018
        AZ1 = 0.0D0    ! UPDATE 2018
        AZ2 = 0.0D0    ! UPDATE 2018
        ENDIF
           IF (NEZ.EQ.1.0)THEN
              !  IF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z  ! UPDATE 2018
              !  CALCULATE NEW LOCATION OF ELEMENET
                 IF (AZ1.GE.AZ2) AZ1 = WATERDEPTH
                 IF (AZ1.LT.AZ2) AZ2 = WATERDEPTH
              !  SLOPE 3D LINEAR
                 ! FINDING VECTOR
                 VEC_Y = AMAXY - AMINY
                 VEC_X = AMAXX - AMINX
                 VEC_Z = AMAXZ - AMINZ
                 ! LINEAR 3D (UNKNOWN X,Y)
                 ANEW_X = (((WATERDEPTH-AMINZ)/VEC_Z)*VEC_X) + AMINX
                 ANEW_Y = (((WATERDEPTH-AMINZ)/VEC_Z)*VEC_Y) + AMINY
                 AX1    = ANEW_X
                 AX2    = AMINX
                 AY1    = ANEW_Y
                 AY2    = AMINY
                 AZ1    = WATERDEPTH
                 AZ2    = AMINZ
              !ENDIF
           ELSEIF (NEZ.EQ.2.0)THEN
           ! SUBMERGED CASE 
           ELSEIF (NEZ.EQ.3.0)THEN
               IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! X
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0   
              ENDIF           
           ENDIF
        ENDIF
            
      ! CALCULATE NEW LENGTH
        DX = ABS(AX2-AX1)
        DY = ABS(AY2-AY1)
        DZ = ABS(AZ2-AZ1)
        ALENGTH  = SQRT((DX**2)+(DY**2)+(DZ**2))
        
      ! OUTPUT THE COORDINATE OF THE CENTER OF GRAVITY OF EACH COMPOSITE PART OF THE BODY
        BUOCYX = ((AX1+AX2)/2D0)
        BUOCYY = ((AY1+AY2)/2D0)
        BUOCYZ = ((AZ1+AZ2)/2D0)
      
      ! ---------------------------------------------------------------------------
      
      ELSEIF (NOP.EQ.4.0D0)THEN
      ! FOR TAPPER SECTION
        NBUX = 0.0
        NBUY = 0.0
        NBUZ = 0.0
        NEX  = 0.0
        NEY  = 0.0
        NEZ  = 0.0
        IF (AX1.GT.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        ELSEIF (AX1.LT.AX2)THEN
        AMAXX = AX2
        AMINX = AX1
        ELSEIF (AX1.EQ.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        NBUX  = 1.0
        ENDIF
        
        IF (AY1.GT.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        ELSEIF (AY1.LT.AY2)THEN
        AMAXY = AY2
        AMINY = AY1
        ELSEIF (AY1.EQ.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        NBUY  = 1.0
        ENDIF
        
        IF (AZ1.GT.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        ELSEIF (AZ1.LT.AZ2)THEN
        AMAXZ = AZ2
        AMINZ = AZ1
        ELSEIF (AZ1.EQ.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        NBUZ  = 1.0
        ENDIF
        
           DO I = 1,100
          
           IF (AX1.LE.AX2)THEN
           XL(0)  = AX1
           ELSEIF (AX2.GT.AX1)THEN
           XL(0)  = AX2
           ENDIF
           
           IF (AY1.LE.AY2)THEN
           YL(0)  = AY1
           ELSEIF (AY2.GT.AY1)THEN
           YL(0)  = AY2
           ENDIF
           
           IF (AZ1.LE.AZ2)THEN
           ZL(0)  = AZ1
           ELSEIF (AZ2.GT.AZ1)THEN
           ZL(0)  = AZ2
           ENDIF
           
           
           ALX    = (AMAXX-AMINX)/100D0
           XL(I)  = XL(I-1)+ALX
           
           ALY    = (AMAXY-AMINY)/100D0
           YL(I)  = YL(I-1)+ALY
           
           ALZ    = (AMAXZ-AMINZ)/100D0
           ZL(I)  = ZL(I-1)+ALZ
           
          ! OUTPUT THE COORDINATE OF THE CENTER OF GRAVITY OF EACH COMPOSITE PART OF THE BODY ( TAPPER )
           ATMIDX(I) = (XL(I-1)+XL(I))/2D0
           ATMIDY(I) = (YL(I-1)+YL(I))/2D0
           ATMIDZ(I) = (ZL(I-1)+ZL(I))/2D0
           
           ENDDO
           
           
      ELSEIF (NOP.EQ.5.0D0)THEN
        NBUX = 0.0
        NBUY = 0.0
        NBUZ = 0.0
        NEX  = 0.0
        NEY  = 0.0
        NEZ  = 0.0
        IF (AX1.GT.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        ELSEIF (AX1.LT.AX2)THEN
        AMAXX = AX2
        AMINX = AX1
        ELSEIF (AX1.EQ.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        NBUX  = 1.0
        ENDIF
        
        IF (AY1.GT.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        ELSEIF (AY1.LT.AY2)THEN
        AMAXY = AY2
        AMINY = AY1
        ELSEIF (AY1.EQ.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        NBUY  = 1.0
        ENDIF
        
        IF (AZ1.GT.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        ELSEIF (AZ1.LT.AZ2)THEN
        AMAXZ = AZ2
        AMINZ = AZ1
        ELSEIF (AZ1.EQ.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        NBUZ  = 1.0
        ENDIF
        
        IF (NGRAV.EQ.1)THEN
        IF (AMAXX.GT.WATERDEPTH.AND.AMINX.LT.WATERDEPTH)THEN ! 
        NEX = 1    ! BETWEEN SUMMERGED AND IN AIR
        ELSEIF (AMAXX.LE.WATERDEPTH.AND.AMINX.LE.WATERDEPTH)THEN ! 
        NEX = 2    ! SUBMERGED
        ELSEIF (AMAXX.GE.WATERDEPTH.AND.AMINX.GE.WATERDEPTH)THEN ! 
        NEX = 3    ! IN AIR
        ENDIF
           IF (NEX.EQ.1.0)THEN
              IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Z
                 IF (AX1.GE.AX2)THEN
                 AX1 = WATERDEPTH
                 ELSEIF (AY1.LT.AY2)THEN
                 AX2 = WATERDEPTH
                 ENDIF        
              ENDIF
           ELSEIF (NEX.EQ.2.0)THEN
           ! SUBMERGED CASE
           ELSEIF (NEX.EQ.3.0)THEN
               IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! X
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0   
              ENDIF           
           ENDIF
        ELSEIF (NGRAV.EQ.2)THEN
        IF (AMAXY.GT.WATERDEPTH.AND.AMINY.LT.WATERDEPTH)THEN ! 
        NEY = 1    ! BETWEEN SUMMERGED AND IN AIR
        ELSEIF (AMAXY.LE.WATERDEPTH.AND.AMINY.LE.WATERDEPTH)THEN ! 
        NEY = 2    ! SUBMERGED
        ELSEIF (AMAXY.GE.WATERDEPTH.AND.AMINY.GE.WATERDEPTH)THEN ! 
        NEY = 3    ! IN AIR
        ENDIF
           IF (NEY.EQ.1.0)THEN
              IF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y  
                 IF (AY1.GE.AY2)THEN
                 AY1 = WATERDEPTH
                 ELSEIF (AY1.LT.AY2)THEN
                 AY2 = WATERDEPTH
                 ENDIF        
              ENDIF
           ELSEIF (NEY.EQ.2.0)THEN
           ! SUBMERGED CASE          
           ELSEIF (NEY.EQ.3.0)THEN
               IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! X
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0   
              ENDIF           
           ENDIF
        ELSEIF (NGRAV.EQ.3)THEN
        IF (AMAXZ.GT.WATERDEPTH.AND.AMINZ.LT.WATERDEPTH)THEN ! 
        NEZ = 1    ! BETWEEN SUMMERGED AND IN AIR
        ELSEIF (AMAXZ.LE.WATERDEPTH.AND.AMINZ.LE.WATERDEPTH)THEN ! 
        NEZ = 2    ! SUBMERGED
        ELSEIF (AMAXZ.GE.WATERDEPTH.AND.AMINZ.GE.WATERDEPTH)THEN ! 
        NEZ = 3    ! IN AIR
        ENDIF
           IF (NEZ.EQ.1.0)THEN
              IF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z
                 IF (AZ1.GE.AZ2)THEN
                 AZ1 = WATERDEPTH
                 ELSEIF (AY1.LT.AY2)THEN
                 AZ2 = WATERDEPTH
                 ENDIF        
              ENDIF
           ELSEIF (NEZ.EQ.2.0)THEN
           ! SUBMERGED CASE 
           ELSEIF (NEZ.EQ.3.0)THEN
               IF (NBUY.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! X
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUZ.EQ.1.0)THEN ! Y
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0
              ELSEIF (NBUX.EQ.1.0.AND.NBUY.EQ.1.0)THEN ! Z
                 AX1 = 0.0D0
                 AX2 = 0.0D0
                 AY1 = 0.0D0
                 AY2 = 0.0D0
                 AZ1 = 0.0D0
                 AZ2 = 0.0D0   
              ENDIF           
           ENDIF
        ENDIF
          
           DO I = 1,100
          
           IF (AX1.LE.AX2)THEN
           XL(0)  = AX1
           ELSEIF (AX2.GT.AX1)THEN
           XL(0)  = AX2
           ENDIF
           
           IF (AY1.LE.AY2)THEN
           YL(0)  = AY1
           ELSEIF (AY2.GT.AY1)THEN
           YL(0)  = AY2
           ENDIF
           
           IF (AZ1.LE.AZ2)THEN
           ZL(0)  = AZ1
           ELSEIF (AZ2.GT.AZ1)THEN
           ZL(0)  = AZ2
           ENDIF
           
        IF (AX1.GT.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        ELSEIF (AX1.LT.AX2)THEN
        AMAXX = AX2
        AMINX = AX1
        ELSEIF (AX1.EQ.AX2)THEN
        AMAXX = AX1
        AMINX = AX2
        NBUX  = 1.0
        ENDIF
        
        IF (AY1.GT.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        ELSEIF (AY1.LT.AY2)THEN
        AMAXY = AY2
        AMINY = AY1
        ELSEIF (AY1.EQ.AY2)THEN
        AMAXY = AY1
        AMINY = AY2
        NBUY  = 1.0
        ENDIF
        
        IF (AZ1.GT.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        ELSEIF (AZ1.LT.AZ2)THEN
        AMAXZ = AZ2
        AMINZ = AZ1
        ELSEIF (AZ1.EQ.AZ2)THEN
        AMAXZ = AZ1
        AMINZ = AZ2
        NBUZ  = 1.0
        ENDIF
           
           
           ALX    = (AMAXX-AMINX)/100D0
           XL(I)  = XL(I-1)+ALX
              
           ALY    = (AMAXY-AMINY)/100D0
           YL(I)  = YL(I-1)+ALY
           
           ALZ    = (AMAXZ-AMINZ)/100D0
           ZL(I)  = ZL(I-1)+ALZ
           
          ! OUTPUT THE COORDINATE OF THE CENTER OF GRAVITY OF EACH COMPOSITE PART OF THE BODY ( TAPPER )
           ATMIDX(I) = (XL(I-1)+XL(I))/2D0
           ATMIDY(I) = (YL(I-1)+YL(I))/2D0
           ATMIDZ(I) = (ZL(I-1)+ZL(I))/2D0
           
              IF (AX1.EQ.0.0D0.AND.AX2.EQ.0.0D0)THEN
              ATMIDX(I)  = 0.0D0
              ENDIF
              IF (AY1.EQ.0.0D0.AND.AY2.EQ.0.0D0)THEN
              ATMIDY(I)  = 0.0D0
              ENDIF
              IF (AZ1.EQ.0.0D0.AND.AZ2.EQ.0.0D0)THEN
              ATMIDZ(I)  = 0.0D0
              ENDIF
           ENDDO
       
        ! CALCULATE NEW LENGTH
        DX = ABS(AX1-AX2)
        DY = ABS(AY1-AY2)
        DZ = ABS(AZ1-AZ2)
        ALENGTH  = SQRT((DX**2)+(DY**2)+(DZ**2))
           
      ENDIF 
      
      END

C  ========================================================================================================
   ! SET AND SELECT SECTION ELEMENT PROPERTIES 
      SUBROUTINE ARRAYELEMENT (NELEMENT,NSECTION)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /DISEC/ DIMENPROP(5000,11)
     	COMMON /SECTIONPROP/ SECTIONPROP(5000,32)
     	
     	! OUTPUT FOR CALCULATION
     	COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
     

CC      DO 10 J=1,NELE
CC           IF (NELEMENT.EQ.NM(J))THEN
CC              I = IGEO(J)
CC              EXIT
CC           ENDIF
CC10    CONTINUE 

C     SONGSAK REMOVE ABOVE LOOP AND REPLACE HERE OCT2019      
      IGM = NELEMENT
	CALL INTFILL('OGDM',IGST,3,IGM,0) !CALL GEOMETRIC SET INDEX
      I = IGST
      
         ! DIMENSION PROPERTIES
         NSECTION = DIMENPROP(I,1)
         
         IF (NSECTION.EQ.1)THEN  ! RECTANGULAR
         HREC = DIMENPROP(I,2)
         BREC = DIMENPROP(I,3)
         SECTIONT = SECTIONPROP(I,15)/(HREC/2D0)
         SECTIONS = SECTIONPROP(I,14)/(BREC/2D0)
         PLASTICT = BREC*(HREC**2D0)/4D0
         PLASTICS = HREC*(BREC**2D0)/4D0
         
         ELSEIF (NSECTION.EQ.2)THEN ! I SHAPE
         BFISHAPE1  = DIMENPROP(I,5)+DIMENPROP(I,6)+DIMENPROP(I,7)
         BFISHAPE2  = DIMENPROP(I,3)+DIMENPROP(I,4)+DIMENPROP(I,7)
          IF (BFISHAPE1.GT.BFISHAPE2) BFISHAPE = BFISHAPE2
          IF (BFISHAPE1.LE.BFISHAPE2) BFISHAPE = BFISHAPE1
         TFISHAPE   = MINVAL(DIMENPROP(I,8:11))
         DISHAPE    = DIMENPROP(I,2)-(TFISHAPE*2D0)
         BWISHAPE   = DIMENPROP(I,2)
         TWISHAPE   = DIMENPROP(I,7)
         ROOTISHAPE = 0.0D0
         SECTIONT   = SECTIONPROP(I,15)/(BWISHAPE/2D0)
         SECTIONS   = SECTIONPROP(I,14)/(BFISHAPE/2D0)
         PLASTICT   = ((BFISHAPE*TFISHAPE)*(BWISHAPE-TFISHAPE))+(0.25D0*TWISHAPE*(BWISHAPE-2D0*TFISHAPE)**2D0)
         PLASTICS   = ((BFISHAPE**2D0)*TFISHAPE/2D0)+((0.25D0*(TWISHAPE**2D0))*(BWISHAPE-2D0*TFISHAPE))
         
         ELSEIF (NSECTION.EQ.3)THEN ! H SHAPE
         BFHSHAPE1  = DIMENPROP(I,5)+DIMENPROP(I,6)+DIMENPROP(I,7)
         BFHSHAPE2  = DIMENPROP(I,3)+DIMENPROP(I,4)+DIMENPROP(I,7)
          IF (BFHSHAPE1.GT.BFHSHAPE2) BFHSHAPE = BFHSHAPE2
          IF (BFHSHAPE1.LE.BFHSHAPE2) BFHSHAPE = BFHSHAPE1
         TFHSHAPE   = MINVAL(DIMENPROP(I,8:11))
         DHSHAPE    = DIMENPROP(I,2)-(TFHSHAPE*2D0)
         BWHSHAPE   = DIMENPROP(I,2)
         TWHSHAPE   = DIMENPROP(I,7)
         ROOTHSHAPE = 0.0D0
         SECTIONT   = SECTIONPROP(I,15)/(BWHSHAPE/2D0)
         SECTIONS   = SECTIONPROP(I,14)/(BFHSHAPE/2D0)
         PLASTICT   = ((BFHSHAPE*TFHSHAPE)*(BWHSHAPE-TFHSHAPE))+(0.25D0*TWHSHAPE*(BWHSHAPE-2D0*TFHSHAPE)**2D0)
         PLASTICS   = ((BFHSHAPE**2D0)*TFHSHAPE/2D0)+((0.25D0*(TWHSHAPE**2D0))*(BWHSHAPE-2D0*TFHSHAPE))
         
         ELSEIF (NSECTION.EQ.4)THEN ! ANGLE
         BANGLE   =  DIMENPROP(I,2)
         TANGLE   =  DIMENPROP(I,4)
         HANGLE   =  DIMENPROP(I,3)
         THANGLE  =  DIMENPROP(I,5) 
         AXBAR    =  SECTIONPROP(I,5)
         SECTIONT = SECTIONPROP(I,15)/(SECTIONPROP(I,5))
         SECTIONS = SECTIONPROP(I,14)/(SECTIONPROP(I,4))
         PLASTICT =  ((HANGLE-TANGLE/2D0)*THANGLE)*ABS(SECTIONPROP(I,4)-THANGLE/2D0)
     1               +((BANGLE-THANGLE/2D0)*TANGLE)*ABS(SECTIONPROP(I,4)-BANGLE/2D0)
         PLASTICS =  ((HANGLE-TANGLE/2D0)*THANGLE)*ABS(SECTIONPROP(I,5)-THANGLE/2D0)
     1               +((BANGLE-THANGLE/2D0)*TANGLE)*ABS(SECTIONPROP(I,5)-BANGLE/2D0)  
         
         ELSEIF (NSECTION.EQ.5)THEN ! IL ANGLE
         BANGLE   =  DIMENPROP(I,2)
         TANGLE   =  DIMENPROP(I,4)
         HANGLE   =  DIMENPROP(I,3)
         THANGLE  =  DIMENPROP(I,5)
         AXBAR    =  SECTIONPROP(I,5)
         SECTIONT = SECTIONPROP(I,15)/(SECTIONPROP(I,5))
         SECTIONS = SECTIONPROP(I,14)/(SECTIONPROP(I,4))
         PLASTICT =  ((HANGLE-TANGLE/2D0)*THANGLE)*ABS(SECTIONPROP(I,4)-THANGLE/2D0)
     1               +((BANGLE-THANGLE/2D0)*TANGLE)*ABS(SECTIONPROP(I,4)-BANGLE/2D0)
         PLASTICS =  ((HANGLE-TANGLE/2D0)*THANGLE)*ABS(SECTIONPROP(I,5)-THANGLE/2D0)
     1               +((BANGLE-THANGLE/2D0)*TANGLE)*ABS(SECTIONPROP(I,5)-BANGLE/2D0)  
         
         ELSEIF (NSECTION.EQ.6)THEN ! C1 SHAPE
         BFCHANNEL  = DIMENPROP(I,2)
         TFCHANNEL  = MAXVAL(DIMENPROP(I,5:6))
         DCHANNEL   = DIMENPROP(I,3)-(TFCHANNEL*2D0)
         TWCHANNEL  = DIMENPROP(I,6)
         HCHANNEL   = DIMENPROP(I,3)
         ROOTCSHAPE = 0.0D0
         CXBAR      = SECTIONPROP(I,5)
         SECTIONT   = ((BFCHANNEL*HCHANNEL**2D0)/6D0)-(((BFCHANNEL-TWCHANNEL)*DCHANNEL**2D0)/(6D0*HCHANNEL))
         SECTIONS   = SECTIONPROP(I,14)/SECTIONPROP(I,4)
         PLASTICT   = (HCHANNEL*TWCHANNEL)*ABS(SECTIONPROP(I,4)-TWCHANNEL/2D0)+
     1                (BFCHANNEL*TFCHANNEL)*ABS(SECTIONPROP(I,4)-BFCHANNEL/2D0)*2D0
         PLASTICS   = (HCHANNEL*TWCHANNEL)*ABS(SECTIONPROP(I,5)-TWCHANNEL/2D0)+
     1                (BFCHANNEL*TFCHANNEL)*ABS(SECTIONPROP(I,5)-BFCHANNEL/2D0)*2D0       
         
         ELSEIF (NSECTION.EQ.7)THEN ! C2 SHAPE
         BFCHANNEL  = MAXVAL(DIMENPROP(I,2:3))
         TFCHANNEL1 = MAXVAL(DIMENPROP(I,7:8))
         TFCHANNEL2 = MAXVAL(DIMENPROP(I,10:11))
            IF (TFCHANNEL1.GT.TFCHANNEL2) TFCHANNEL = TFCHANNEL2
            IF (TFCHANNEL1.LE.TFCHANNEL2) TFCHANNEL = TFCHANNEL1
         DCHANNEL   = DIMENPROP(I,4)-(TFCHANNEL*2D0)
         TWCHANNEL  = DIMENPROP(I,9)
         HCHANNEL   = DIMENPROP(I,4)
         ROOTCSHAPE = 0.0D0
         CXBAR      = SECTIONPROP(I,5)
         SECTIONT   = ((BFCHANNEL*HCHANNEL**2D0)/6D0)-(((BFCHANNEL-TWCHANNEL)*DCHANNEL**2D0)/(6D0*HCHANNEL))
         SECTIONS   = SECTIONPROP(I,14)/SECTIONPROP(I,4)
         PLASTICT   = (HCHANNEL*TWCHANNEL)*ABS(SECTIONPROP(I,4)-TWCHANNEL/2D0)+
     1                (BFCHANNEL*TFCHANNEL)*ABS(SECTIONPROP(I,4)-BFCHANNEL/2D0)*2D0
         PLASTICS   = (HCHANNEL*TWCHANNEL)*ABS(SECTIONPROP(I,5)-TWCHANNEL/2D0)+
     1                (BFCHANNEL*TFCHANNEL)*ABS(SECTIONPROP(I,5)-BFCHANNEL/2D0)*2D0 
         
         ELSEIF (NSECTION.EQ.10)THEN ! T SHAPE
         BFTSHAPE  = DIMENPROP(I,3)+DIMENPROP(I,4)+DIMENPROP(I,5)
         TFTSHAPE  = MAXVAL(DIMENPROP(I,6:7))
         DTSHAPE   = DIMENPROP(I,2)- TFTSHAPE 
         BWTSHAPE  = DIMENPROP(I,2)
         TWTSHAPE  = DIMENPROP(I,5)
         TYBAR     = SECTIONPROP(I,4)            
         SECTIONT  = SECTIONPROP(I,15)/SECTIONPROP(I,5)
         SECTIONS  = SECTIONPROP(I,14)/SECTIONPROP(I,4)
         PLASTICT  = SECTIONT*1.14D0 ! AISC
         PLASTICS  = SECTIONS*1.14D0 ! AISC
         
         ELSEIF (NSECTION.EQ.11)THEN ! DOUBLE T SHAPE
         BFTSHAPE  = DIMENPROP(I,3)+DIMENPROP(I,4)+DIMENPROP(I,5)+DIMENPROP(I,9)+DIMENPROP(I,10)
         TFTSHAPE  = MAXVAL(DIMENPROP(I,6:8))
         DTSHAPE   = DIMENPROP(I,2)- TFTSHAPE 
         BWTSHAPE  = DIMENPROP(I,2)
         TWTSHAPE  = DIMENPROP(I,9)+DIMENPROP(I,10)
         TYBAR     = SECTIONPROP(I,4)
         SECTIONT  = SECTIONPROP(I,15)/SECTIONPROP(I,5)
         SECTIONS  = SECTIONPROP(I,14)/SECTIONPROP(I,4)
         PLASTICT  = SECTIONT*1.14D0 ! AISC
         PLASTICS  = SECTIONS*1.14D0 ! AISC
         
         ELSEIF (NSECTION.EQ.14)THEN ! CIRCULAR SOILD SHAPE
         DROUND    = DIMENPROP(I,2)
         SECTIONT  = SECTIONPROP(I,15)/(DROUND/2D0) 
         SECTIONS  = SECTIONPROP(I,14)/(DROUND/2D0)
         PLASTICT  = (DROUND**3D0)/6D0
         PLASTICS  = (DROUND**3D0)/6D0
         
         
         ELSEIF (NSECTION.EQ.15)THEN ! SOILD ROUND
         SECTIONT = SECTIONPROP(I,15)/SECTIONPROP(I,5)
         SECTIONS = SECTIONPROP(I,14)/SECTIONPROP(I,4)
         PLASTICT  = SECTIONT*1.14D0 ! AISC
         PLASTICS  = SECTIONS*1.14D0 ! AISC
         ELSEIF (NSECTION.EQ.16)THEN ! DOUBLE C SHAPE
         BFCHANNEL  = DIMENPROP(I,3)*2D0
         TFCHANNEL  = DIMENPROP(I,5)
         DCHANNEL   = DIMENPROP(I,2)-(TFCHANNEL*2D0)
         TWCHANNEL  = DIMENPROP(I,4)*2D0
         HCHANNEL   = DIMENPROP(I,2)
         ROOTCSHAPE = 0.0D0
         CXBAR      = SECTIONPROP(I,5)
         SECTIONT   = SECTIONPROP(I,15)/(HCHANNEL/2D0) 
         SECTIONS   = SECTIONPROP(I,14)/(BFCHANNEL/2D0)
         PLASTICT   = SECTIONT*1.14D0 ! AISC
         PLASTICS   = SECTIONS*1.14D0 ! AISC
         
         ELSEIF (NSECTION.EQ.17)THEN ! DOUBLE L SHAPE
         BANGLE    =  DIMENPROP(I,3)*2D0
         TANGLE    =  DIMENPROP(I,5)
         HANGLE    =  DIMENPROP(I,2)
         THANGLE   =  DIMENPROP(I,4)*2D0
         AXBAR     =  SECTIONPROP(I,5)
         SECTIONT  = SECTIONPROP(I,15)/(HANGLE)
         SECTIONS  = SECTIONPROP(I,14)/(BANGLE)
         PLASTICT  = SECTIONT*1.14D0 ! AISC
         PLASTICS  = SECTIONS*1.14D0 ! AISC
         
         ELSEIF (NSECTION.EQ.18)THEN ! DOUBLE IL SHAPE
         BANGLE  =  DIMENPROP(I,3)*2D0
         TANGLE  =  DIMENPROP(I,5)
         HANGLE  =  DIMENPROP(I,2)
         THANGLE =  DIMENPROP(I,4)*2D0
         AXBAR   =  SECTIONPROP(I,5)
         SECTIONT  = SECTIONPROP(I,15)/(HANGLE)
         SECTIONS  = SECTIONPROP(I,14)/(BANGLE)
         PLASTICT  = SECTIONT*1.14D0 ! AISC
         PLASTICS  = SECTIONS*1.14D0 ! AISC
         
         ELSEIF (NSECTION.EQ.20)THEN ! BOX
         BBOX      = DIMENPROP(I,3)
         TFBOX     = MINVAL(DIMENPROP(I,5:6))
         DBOX      = DIMENPROP(I,2)-(TFBOX*2D0)
         HBOX      = DIMENPROP(I,2)
         TWBOX     = DIMENPROP(I,4)
         ROOTBOX   = 0.0D0
         SECTIONT  = SECTIONPROP(I,15)/(HBOX/2D0)
         SECTIONS  = SECTIONPROP(I,14)/(BBOX/2D0)
         PLASTICT  = (BBOX*(HBOX**2D0)/4D0)-((BBOX-2D0*TFBOX)*(((HBOX/2D0)-TFBOX)**2D0))
         PLASTICS  = (HBOX*(BBOX**2D0)/4D0)-((HBOX-2D0*TFBOX)*(((BBOX/2D0)-TFBOX)**2D0))
         
         ELSEIF (NSECTION.EQ.21)THEN ! PIPE ( DIMENSION )
         DPIPE    = DIMENPROP(I,2)*2D0
         TPIPE    = DIMENPROP(I,2)-DIMENPROP(I,3)
         SECTIONT = 3.141592654D0*(DPIPE**4D0-(DPIPE-2D0*TPIPE)**4D0)/(32D0*DPIPE)
         SECTIONS = 3.141592654D0*(DPIPE**4D0-(DPIPE-2D0*TPIPE)**4D0)/(32D0*DPIPE)
         PLASTICT = ((DPIPE**3D0)/6D0)-(((DPIPE-(2D0*TPIPE))**3D0)/6D0)
         PLASTICS = ((DPIPE**3D0)/6D0)-(((DPIPE-(2D0*TPIPE))**3D0)/6D0)
         ELSEIF (NSECTION.NE.0)THEN
         SECTIONT  = SECTIONPROP(I,15)/(DPIPE/2D0)
         SECTIONS  = SECTIONPROP(I,14)/(DPIPE/2D0)
         PLASTICT  = SECTIONT*1.14D0 ! AISC
         PLASTICS  = SECTIONS*1.14D0 ! AISC
         ELSEIF (NSECTION.EQ.0)THEN
         ! CALAULTE AS RECTANGULAR
         NSECTION = 1
         AREA     = SECTIONPROP(I,1)
         HREC     = SQRT(AREA)
         BREC     = SQRT(AREA)
         SECTIONT = SECTIONPROP(I,15)/(HREC/2D0)
         SECTIONS = SECTIONPROP(I,14)/(BREC/2D0)
         PLASTICT = BREC*(HREC**2D0)/4D0
         PLASTICS = HREC*(BREC**2D0)/4D0
         ENDIF
         ! SECTION PROPERTIES
         AJ       = SECTIONPROP(I,29)
         CW       = SECTIONPROP(I,31)
         AREA     = SECTIONPROP(I,1)
         AIS      = SECTIONPROP(I,14)
         AIT      = SECTIONPROP(I,15)
         ARGS     = SQRT(AREA/AIS)
         ARGT     = SQRT(AREA/AIT)
         !NSECTION = 1
C         ENDIF      
C      ENDDO
      END
C  ========================================================================================================
      ! THIS SUBRUTINE FOR GET ONLY CIRCULAR SECTION
      SUBROUTINE MVALUE (SS,TT,MGPS,ISET,ND)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /STORESECTION/ D1(10000),D2(10000),D3(10000),D4(10000),XBAR(10000)
     1                      ,YBAR(10000),MGEO
      
      DIMENSION SS(2,40),TT(40),FIRST(40),SECOND(40),A(MGPS),B(MGPS),C(MGPS)
      
      IF (ND.EQ.1.0D0)THEN
      DO I=1,40,2
      FIRST(I)=SS(1,I)
      ENDDO
      ELSEIF (ND.EQ.2.0D0)THEN
      DO I=1,40,2
      FIRST(I)=SS(2,I)
      ENDDO
      ENDIF
      
      IF (SS(1,1).NE.SS(1,3))THEN ! CIRCULAR SECTION
      
      DO J=2,40,2
        FIRST(J)=0.0D0
      ENDDO
      
      A(ISET)=MAXVAL(FIRST)
      
      DO I=1,40
        IF(A(ISET).EQ.FIRST(I))THEN
        FIRST(I)=0.0D0
        ENDIF
      ENDDO
      
      IF (ND.EQ.2.0D0)THEN
      C(ISET) = FIRST(7)
      D3(ISET)= A(ISET)
      D4(ISET)= C(ISET)
      MGEO=MGPS
      ELSEIF (ND.EQ.1.0D0)THEN
      B(ISET) = FIRST(7)
      D1(ISET)= A(ISET)
      D2(ISET)= B(ISET)
      MGEO=MGPS
      ENDIF
      
      ELSEIF (SS(1,1).EQ.SS(1,3))THEN ! RECTANGURA SECTION
        IF (ND.EQ.2)THEN
        XBAR(ISET) = FIRST(1)
        ELSEIF (ND.EQ.1)THEN
        YBAR(ISET) = FIRST(1)
        ENDIF
      
      
      ENDIF
      
      END
C  ========================================================================================================
      ! SELECT SECTION DIAMETER
      SUBROUTINE MAXSECTION (NELEMENT,DOUT,DIN,DOUT1,DIN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /STORESECTION/ D1(10000),D2(10000),D3(10000),D4(10000),XBAR(10000)
     1                      ,YBAR(10000),MGEO
      
CC        DO I=1,NELE
CC           IF (NELEMENT.EQ.NM(I))THEN
CC              IGE=IGEO(I)
CC              DOUT  = D1(IGE) ! R NOT DIAMETER
CC              DIN   = D2(IGE) ! R NOT DIAMETER
CC              DOUT1 = D3(IGE) ! R NOT DIAMETER
CC              DIN2  = D4(IGE) ! R NOT DIAMETER
CC           ENDIF
CC        ENDDO 
        
C     SONGSAK REMOVE ABOVE LOOP AND REPLACE HERE OCT2019      
      IGM = NELEMENT
	CALL INTFILL('OGDM',IGST,3,IGM,0) !CALL GEOMETRIC SET INDEX
      
      IF(IGST.GT.0) THEN
          IGE = IGST
          DOUT  = D1(IGE) ! R NOT DIAMETER
          DIN   = D2(IGE) ! R NOT DIAMETER
          DOUT1 = D3(IGE) ! R NOT DIAMETER
          DIN2  = D4(IGE) ! R NOT DIAMETER
      ELSE
          WRITE(*,*) '*********************************************' 
          WRITE(*,*) 'WARNING WRONG IGSET NUMBER IN SUB. MAXSECTION' 
          WRITE(*,*) '*********************************************' 
          STOP
      ENDIF
             
      END
C  ========================================================================================================
      ! STOREAGE SECTION PROPERTIES
      SUBROUTINE PROPERTIES (PROP,MEMGID,ISET,IMP,ITYPE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /STOREPROP/ AM(80000),AP(80000),AD(80000),NPROP(80000),NSTEEL(80000),AFU(80000),AFY(80000)
      COMMON /STORE_MOOR/ MOORING_LINE(500,10),MOORISET(5000),MOORIELE(5000),INDEX_MOORING,INDEX_MOORING_ISET
      COMMON /INDEX_MOORING/ ISET_OLD
      DIMENSION PROP(IMP)
      NPROP(MEMGID)=ISET
      AM(MEMGID) = PROP(1)
      AP(MEMGID) = PROP(2)
      AD(MEMGID) = PROP(5)
      NSTEEL(MEMGID) = PROP(14)
      AFY(MEMGID)     = PROP(15)
      AFU(MEMGID)     = PROP(16)
      
      IF (ITYPE.EQ.2)THEN ! MOORING LINE
      INDEX_MOORING_ISET = INDEX_MOORING_ISET + 1
      INDEX_MOORING      = INDEX_MOORING + 1
      MOORISET(INDEX_MOORING_ISET) = ISET
      MOORIELE(INDEX_MOORING_ISET) = MEMGID
      IF (ISET.NE.ISET_OLD) INDEX_MOORING = 1
      ISET_OLD      = ISET
      MOORING_LINE(INDEX_MOORING,ISET) = MEMGID
      ENDIF
      
      END      
C  ========================================================================================================
      ! SELECT SECTION PROPERTIES
      SUBROUTINE SELECTPROP (NELEMENT,AMODULUS,POSSION,DENSITY,NSTEELGRAD,FU,FY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /STOREPROP/ AM(80000),AP(80000),AD(80000),NPROP(80000),NSTEEL(80000),AFU(80000),AFY(80000)
      
      AMODULUS    = AM(NELEMENT)
      POSSION     = AP(NELEMENT)
      DENSITY     = AD(NELEMENT)
      NSTEELGRAD  = NSTEEL(NELEMENT)
      FU          = AFU(NELEMENT)
      FY          = AFY(NELEMENT)
      
      END

! =======================================================================================  
      SUBROUTINE CALLNUMNODE_F (NELEMENT,N1,N2)
C     PRODUCE BY SONGSAK TO CALL N1,N2 FOR FRAME ELEMENT
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      DIMENSION LMNM(100)
      
      IGM = NELEMENT
	CALL INTFILL('OGDM',IEG,1,IGM,0) !CALL ELEMENT GROUP NUMBER
      CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !CALL FILE STORING THE CONNECTIVITY FOR SPECIFIC ELEM. GROUP IEG
      LMNM = 0 !VECTOR STORE ELEMENT CONNECTIVITY NUMBER
      CALL READCON(IGM,LMNM,NFL4)
      N1 = LMNM(1)
      N2 = LMNM(2)
         
      END
C  ========================================================================================================
      SUBROUTINE CALLNUMFRAME (NFRAM)
C     PRODUCE BY SONGSAK TO NUMBER OF ALL FRAME ELEMENT
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
C     ---------------------------------------      
C     SONGSAK MODIFY TO PRINT FASTER  OCT2019   
      
C	CALLING NEG
	CALL INTFILL('%NUB',NEG,1,5 ,0)
      
      NFRAM = 0
      DO 100 IEG = 1,NEG
          
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
      
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !NELE
          
      IF(ITYP.EQ.5) THEN
          NFRAM = NFRAM+NELE
      ENDIF
      
100   ENDDO
          
C     ---------------------------------------   
         
      END
C  ========================================================================================================
      SUBROUTINE WRITE_FRAME_CON (NFILE)
C     PRODUCE BY SONGSAK TO WRITE FRAME CONNECTIVITY TO SPECIFIC FILE ... NFILE
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      ALLOCATABLE LVAL(:,:)
      
C     ---------------------------------------      
C     SONGSAK MODIFY TO PRINT FASTER  OCT2019   
      
C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
      
      ALLOCATE(LVAL(3,NGIDM))
      
      II = 0
      DO 100 IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)  !IELE
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP

      IF(ITYP.EQ.5) THEN
          CALL CALLNUMNODE_F (IGM,N1,N2)
          II = II+1
          LVAL(1,II) = IGM
          LVAL(2,II) = N1
          LVAL(3,II) = N2
      ENDIF
      
100   ENDDO
      
      IF(II.GT.0) WRITE(NFILE,200) LVAL(1:3,1:II)
      
      DEALLOCATE(LVAL)
C     ---------------------------------------   
         
200   FORMAT (I8,2X,I8,2X,I8)

      END
C  ========================================================================================================