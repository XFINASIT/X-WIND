C     ----------------------------------------------------------------
      SUBROUTINE SOIL_PILE_INTERACTION (SD,ALSS2,ALSS3,XSTIFF,LID,LMA,DISP,OPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPT
      CHARACTER(100) NUMBERFILE,DIRPRIMARY
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /STCAS/ ILC
      
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
      COMMON / SOILTEST/ KTOEY
      COMMON / ITERATION / APSTIFF(3,50000),ATSTIFF(3,50000),AQSTIFF(3,50000),AUSTIFF(3,50000),ATSTIFF_RO(3,50000)
      COMMON / SOILRE/ NSOILRETURN
      
      COMMON / SOILSTIFFNESS / SSTIFF(6,1000)
      
      COMMON /NSOILT/ NNODEA(10000)
      
      COMMON / DISSOIL / DATADISP(10000)
      
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV  
      COMMON/COUNT/ICOUNT,NCOUNT
      
      COMMON /TZ/ PREDATATZSAND(6),PREDATATZCLAY(16)
      
      COMMON/PLOTSOIL/ KPCOUNT
      
      COMMON/OVER/POTOL
      
      COMMON/LARGE_PILE/ NSUC
      
C      COMMON /TEST_TOEY/ KTOEY
      
      DIMENSION PDATA(20)
      DIMENSION MAPP(6)
      DIMENSION NPM(10),NPI(10)
      DIMENSION SD(6,NSN)
      
      DIMENSION SBC(37)
      
      DIMENSION DISP(NSN*7)
      
!      DIMENSION NODE(2),NODEOUT(2)
      
      DIMENSION ALSS2(1),ALL3(1),LID(1),LMA(1)
      DIMENSION EDAT(20)
      
      DIMENSION STIFF(3),STIFF_RO(3)
      
      DIMENSION ERROR(15,NSN),DUMEDAT(3)
      
      DIMENSION X(16),NODETOPMAX(10000),NODEBOTTOMAX(10000),NTOP(10000)
      
      DIMENSION DATATZSAND(6),DATATZCLAY(16),DATAPYSAND(30)
      DIMENSION DATAPYCLAY(14),PYCURVE(14,200)
      
      DIMENSION PLOTP(15,1000),PLOTY(15,1000)
      
      DIMENSION CURVEDATA(100,NSN)
      
      DIMENSION NSTIFF_SOIL(1000)
      

C     ==================== COMMENT ====================
C     1) IN PSI PROGRAM Q-Z CURVE OR END-BEARING STIFFNESS WAS CALCULATE RELATE WITH THE T-Z CURVE BY AUTOMATICALY
C     2) PROGRAM HAS BEEN SEPERATE 2 PART 1/USER DEFIND FUNCTION 2/API FUNCITON
C     3) ALL RESULT AND VALIDATION BY FOLLOWING SACS PROGRAM
C     4) PILE SHOULD NOT BE CONNECT ONE ELEMENT TWO NODE
C     5) VERTICAL OR INCLINE POSSIBLE TO CALCULATED
C     =================================================

C     CALL TOTAL GID ELEMENT NUMBER NGIDM ... SONGSAK OCT2019
      CALL LOCATN('OGDM',KGDM,NGDAT,NGIDM,1) 

      
      IF (KTOEY.EQ.0) THEN
         CALL MODIFY_INPUT_SOIL_SPRING
         KTOEY = 1
      RETURN
      ENDIF

      IF (OPT.EQ."") THEN
          
      ELSEIF (OPT.EQ."ITER") THEN
          

      ENDIF
      
C     INTAIL SETTING 
      ICOUNT  = ICOUNT +1
      POTOTAL = 0.0D0
      KPCOUNT  = 1
      WRITE(7101,'(I5)') ICOUNT
      EDAT  = 0.
      ERROR = 0.
      SBC   = 0.
      NCODE = 0.
      NT    = 0.
      NOPT  = 2D0
      ISP   = 0.D0
      POTOL = 0.0D0
      CURVEDATA = 0.D0
      PLOTP  = 0.
      PLOTY  = 0.
      ! READ SOIL DATA FROM STOREAGE
      CALL READSOILDATA ("CALT",NTSOIL,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
       
      PREDATATZSAND = 0.0D0
      PREDATATZCLAY = 0.0D0
      
      DO NISN = 1,NTSOIL
          NCOUNT = NCOUNT + 1

         CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
         CALL READSOILDATA ("CALN",NTSOIL,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
         
         IF (NCODE.EQ.2) THEN
         CALL DIRECT_STIFFNESS_SOIL_PILE (NISN,NSTIFF_SOIL,NNODE_DIREC)
            DO IJK = 1,NNODE_DIREC
            ISN = NSTIFF_SOIL(IJK)
            IF (ISN.EQ.0)THEN
            SD(1,ISN) = SD(1,ISN) + SBC(2) 
            SD(2,ISN) = SD(2,ISN) + SBC(3) 
            SD(3,ISN) = SD(3,ISN) + SBC(4)
            ENDIF
            ENDDO
         GOTO 100
         ENDIF
	   DO I = 1,20
	   CALL RELFILL('ONDS',EDAT(I),I,ISN,0)  
         ENDDO 
         
         IF (ISOLOP.EQ.21.OR.ISOLOP.EQ.1)THEN
          EDAT(1) =  DISP((ISN*6)-5)
          EDAT(2) =  DISP((ISN*6)-4)
          EDAT(3) =  DISP((ISN*6)-3)
          EDAT(4) =  DISP((ISN*6)-2)
          EDAT(5) =  DISP((ISN*6)-1)
          EDAT(6) =  DISP((ISN*6))
         ENDIF
          
          ! DPIPE  = 2.0D0
          GRMMA  = SBC(6)  !20D0
          STRAIN = SBC(16) !0.35D0
          C      = SBC(7)  !100D0
          FREE   = SBC(6)
          ISAND  = SBC(5)     
          NSUC   = SBC(24)
          IF (NCODE.EQ.3) CALL NODETOELEMENT (ISN,NELEMENT,NISN,DEPTH,NNOUT,"LENGT")
          IF (NCODE.NE.0) CALL NODETOELEMENT (ISN,NELEMENT,NISN,DEPTH,NNOUT,"ABOVE")
          IF (NCODE.NE.0) CALL ARRAYELEMENT (NELEMENT,NSECTION)
          IF (NELEMENT.EQ.0.AND.NCODE.NE.3) GOTO 100
          IF (NCODE.NE.3) CALL CALLENGTH (NELEMENT,VREW,0,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                  ,ATMIDX,ATMIDY,ATMIDZ)
      
         IF (NCODE.EQ.3)THEN
            A1 = ABS(MOD(NISN,2))
            IF (A1.EQ.1) THEN
            ISP = ISN
            GOTO 100
            ENDIF
            
            CALL NODETOELEMENT (ISN,NELEMENT,ISP,DEPTH,NNOUT,"LENGG")
            
            CALL ELEMENTLENGTH (ISP,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,"PY")

            IF (DEPTH.EQ.1)THEN
               CALL NODETOELEMENT (ISN,NELEMENT,ISP,DEPTH,NNOUT,"LENGT")
               N_OLD_ISN = ISN
               ISN = ISP
	         DO I = 1,20
	         CALL RELFILL('ONDS',EDAT(I),I,ISN,0) 
               ISP = N_OLD_ISN
               ENDDO 
         
              IF (ISOLOP.EQ.21.OR.ISOLOP.EQ.1)THEN
              EDAT(1) =  DISP((ISN*6)-5)
              EDAT(2) =  DISP((ISN*6)-4)
              EDAT(3) =  DISP((ISN*6)-3)
              IF (NSUC.EQ.1) EDAT(4) =  DISP((ISN*6)-2)
              IF (NSUC.EQ.1) EDAT(5) =  DISP((ISN*6)-1)
              IF (NSUC.EQ.1) EDAT(6) =  DISP((ISN*6))
              ENDIF
              CALL NODETOELEMENT (ISN,NELEMENT,NISN,DEPTH,NNOUT,"ABOVE")
              CALL ARRAYELEMENT (NELEMENT,NSECTION)
            ELSEIF (DEPTH.EQ.2)THEN
              CALL NODETOELEMENT (ISN,NELEMENT,ISP,DEPTH,NNOUT,"LENGT")
            ENDIF
            
             CALL ELEMENTLENGTH (ISP,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,"TZ")
         ! USER DEFINE FUNCTION CURVE NONLINEAR SOIL SPRING SUPPORT
         CALL USERFUNCTIONCURVE(ISN,ISP,EDAT(1:3),NT,SBC,STIFF,ALENGTHPY,DEPTH,NTSOIL)   
         WRITE(7101,'(I5,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5)') ISN,STIFF(1),STIFF(2),STIFF(3)
     1                                                          ,EDAT(1),EDAT(2),EDAT(3)
C         ---- ADD AND UPDATE STIFFNESS ----
C        ERROR (10:12,NUMBER OF NODE POINT) = % ERROR COMPARED WITH OLD STIFFNESS
C        SD(1,NUMBER OF NODE POINT)         = STIFFNESS ADDITIONAL ON X DIRCTION
C        SD(2,NUMBER OF NODE POINT)         = STIFFNESS ADDITIONAL ON Y DIRCTION
C        SD(3,NUMBER OF NODE POINT)         = STIFFNESS ADDITIONAL ON Z DIRCTION
C        AUSTIFF(1:3,NUMBER OF NODE POINT)    = STORAGE OLD STIFFNESS FOR FIND % ERROR
         IF (SBC(1).EQ.1)THEN ! P-Y CURVE
              IF (NGRAV.EQ.1)THEN ! GRAVITY X-DIRECTION
              ERROR(10,ISN)     =  0.0D0
              ERROR(11,ISN)     =  ABS(AUSTIFF(2,ISN) - STIFF(2))/ABS(AUSTIFF(2,ISN))*100D0
              ERROR(12,ISN)     =  ABS(AUSTIFF(3,ISN) - STIFF(3))/ABS(AUSTIFF(3,ISN))*100D0
                 IF (STIFF(2).EQ.0D0) ERROR(11,ISN) = 0.0D0
                 IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(2,ISN)      = SD(2,ISN) - AUSTIFF(2,ISN) 
                 SD(3,ISN)      = SD(3,ISN) - AUSTIFF(3,ISN) 
                 AUSTIFF(2,ISN) = STIFF(2)
                 AUSTIFF(3,ISN) = STIFF(3)
                 SD(2,ISN)      = SD(2,ISN) + STIFF(2)
                 SD(3,ISN)      = SD(3,ISN) + STIFF(3)
              ELSEIF (NGRAV.EQ.2)THEN ! GRAVITY Y-DIRECTION
              ERROR(10,ISN)     =  ABS(AUSTIFF(1,ISN) - STIFF(1))/ABS(AUSTIFF(1,ISN))*100D0
              ERROR(11,ISN)     =  0.0D0
              ERROR(12,ISN)     =  ABS(AUSTIFF(3,ISN) - STIFF(3))/ABS(AUSTIFF(3,ISN))*100D0
                 IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
                 IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(1,ISN)      = SD(1,ISN) - AUSTIFF(1,ISN)
                 SD(3,ISN)      = SD(3,ISN) - AUSTIFF(3,ISN)
                 AUSTIFF(1,ISN) = STIFF(1)
                 AUSTIFF(3,ISN) = STIFF(3)
                 SD(1,ISN)      = SD(1,ISN) + STIFF(1)
                 SD(3,ISN)      = SD(3,ISN) + STIFF(3)
              ELSEIF (NGRAV.EQ.3)THEN  !GRAVITY Z-DIRECTION
              ERROR(10,ISN)     =  ABS(AUSTIFF(1,ISN) - STIFF(1))/ABS(AUSTIFF(1,ISN))*100D0
              ERROR(11,ISN)     =  ABS(AUSTIFF(2,ISN) - STIFF(2))/ABS(AUSTIFF(2,ISN))*100D0
              ERROR(12,ISN)     =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
                 IF (STIFF(2).EQ.0D0) ERROR(11,ISN) = 0.0D0
                 SD(1,ISN)      = SD(1,ISN) - AUSTIFF(1,ISN)
                 SD(2,ISN)      = SD(2,ISN) - AUSTIFF(2,ISN) 
                 AUSTIFF(1,ISN) = STIFF(1)
                 AUSTIFF(2,ISN) = STIFF(2)
                 SD(1,ISN)      = SD(1,ISN) + STIFF(1)
                 SD(2,ISN)      = SD(2,ISN) + STIFF(2)
              ENDIF
              
          ELSEIF (SBC(2).EQ.1)THEN ! T-Z CURVE
              IF (NGRAV.EQ.1)THEN
              ERROR(10,ISN)      =  ABS(ATSTIFF(1,ISN) - STIFF(1))/ABS(ATSTIFF(1,ISN))*100D0
              ERROR(11,ISN)      =  0.0D0
              ERROR(12,ISN)      =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
                 SD(1,ISN) = SD(1,ISN) - ATSTIFF(1,ISN)
                 ATSTIFF(1,ISN)  = STIFF(1)
                 SD(1,ISN)       = SD(1,ISN) + STIFF(1)
              ELSEIF (NGRAV.EQ.2)THEN
              ERROR(10,ISN)      =  0.0D0
              ERROR(12,ISN)      =  0.0D0
              ERROR(11,ISN)      =  ABS(ATSTIFF(2,ISN) - STIFF(2))/ABS(ATSTIFF(2,ISN))*100D0
                 IF (STIFF(2).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(2,ISN)       = SD(2,ISN) - ATSTIFF(2,ISN) 
                 ATSTIFF(2,ISN)  = STIFF(2)
                 SD(2,ISN)       = SD(2,ISN) + STIFF(2)
              ELSEIF (NGRAV.EQ.3)THEN
              ERROR(10,ISN)      =  0.0D0
              ERROR(11,ISN)      =  0.0D0
              ERROR(12,ISN)      =  ABS(ATSTIFF(3,ISN) - STIFF(3))/ABS(ATSTIFF(3,ISN))*100D0
                 IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(3,ISN)       = SD(3,ISN) - ATSTIFF(3,ISN) 
                 ATSTIFF(3,ISN)  = STIFF(3)
                 SD(3,ISN)       = SD(3,ISN) + STIFF(3)
              ENDIF
              
          ELSEIF (SBC(3).EQ.1)THEN ! Q-Z CURVE
              IF (NGRAV.EQ.1)THEN
              ERROR(11,ISN)      =  ABS(AQSTIFF(1,ISN) - STIFF(1))/ABS(AQSTIFF(1,ISN))*100D0
              ERROR(12,ISN)      =  0.0D0
              ERROR(10,ISN)      =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
                 SD(1,ISN) = SD(1,ISN) - AQSTIFF(1,ISN)
                 AQSTIFF(1,ISN)  = STIFF(1)
                 SD(1,ISN)       = SD(1,ISN) + STIFF(1)
              ELSEIF (NGRAV.EQ.2)THEN
              ERROR(10,ISN)      =  0.0D0
              ERROR(12,ISN)      =  0.0D0
              ERROR(11,ISN)      =  ABS(AQSTIFF(2,ISN) - STIFF(2))/ABS(AQSTIFF(2,ISN))*100D0
                 IF (STIFF(2).EQ.0D0) ERROR(11,ISN) = 0.0D0
                 SD(2,ISN)       = SD(2,ISN) - AQSTIFF(2,ISN) 
                 AQSTIFF(2,ISN)  = STIFF(2)
                 SD(2,ISN)       = SD(2,ISN) + STIFF(2)
              ELSEIF (NGRAV.EQ.3)THEN
              ERROR(10,ISN)      =  0.0D0
              ERROR(11,ISN)      =  0.0D0
              ERROR(12,ISN)      =  ABS(AQSTIFF(3,ISN) - STIFF(3))/ABS(AQSTIFF(3,ISN))*100D0
                 IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(3,ISN)       = SD(3,ISN) - AQSTIFF(3,ISN) 
                 AQSTIFF(3,ISN)  = STIFF(3)
                 SD(3,ISN)       = SD(3,ISN) + STIFF(3)
              ENDIF
          ENDIF
          
C          IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
C          IF (STIFF(2).EQ.0D0) ERROR(11,ISN) = 0.0D0
C          IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
          
C          SD(1,ISN) = SD(1,ISN) - AUSTIFF(1,ISN)
C          SD(2,ISN) = SD(2,ISN) - AUSTIFF(2,ISN) 
C          SD(3,ISN) = SD(3,ISN) - AUSTIFF(3,ISN) 
          
C          AUSTIFF(1,ISN) = STIFF(1)
C          AUSTIFF(2,ISN) = STIFF(2)
C          AUSTIFF(3,ISN) = STIFF(3)
          
C          SD(1,ISN) = SD(1,ISN) + STIFF(1)
C          SD(2,ISN) = SD(2,ISN) + STIFF(2)
C          SD(3,ISN) = SD(3,ISN) + STIFF(3)
         GOTO 100
      ENDIF
      
            ! MOD = CLASSIFY ODD AND EVEN
            A1 = ABS(MOD(NISN,2))
            IF (A1.EQ.1) THEN
            ISP = ISN
            GOTO 100
            ENDIF
            
            NTOPOPT = 0.0D0
            
            CALL NODETOELEMENT (ISN,NELEMENT,ISP,DEPTH,NNOUT,"LENGG")
            ! CALCULATE LAYER LENGTH FOR P-Y, T,Z
            ! ------------------------------------------ 
            ! -                                        - |
            ! -                                        - | L1/2
            ! -                                        - |
            ! ------------------------------------------
            ! -                                        - | L2/2
            ! -                                        - |
            ! ------------------------------------------
            ! -                                        - |
            ! -                                        - | L3/2
            ! -                                        - |
            ! -                                        - |
            ! ------------------------------------------
            CALL ELEMENTLENGTH (ISP,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,"PY")
            
            IF (DEPTH.EQ.1)THEN
               CALL NODETOELEMENT (ISN,NELEMENT,ISP,DEPTH,NNOUT,"LENGT")
               N_OLD_ISN = ISN
               ISN = ISP
	         DO I = 1,20
	         CALL RELFILL('ONDS',EDAT(I),I,ISN,0)  
               ISP = N_OLD_ISN
               ENDDO 
               
         
              IF (ISOLOP.EQ.21.OR.ISOLOP.EQ.1)THEN
              EDAT(1) =  DISP((ISN*6)-5)
              EDAT(2) =  DISP((ISN*6)-4)
              EDAT(3) =  DISP((ISN*6)-3)
              IF (NSUC.EQ.1) EDAT(4) =  DISP((ISN*6)-2)
              IF (NSUC.EQ.1) EDAT(5) =  DISP((ISN*6)-1)
              IF (NSUC.EQ.1) EDAT(6) =  DISP((ISN*6))
              ENDIF
              CALL NODETOELEMENT (ISN,NELEMENT,NISN,DEPTH,NNOUT,"ABOVE")
              CALL ARRAYELEMENT (NELEMENT,NSECTION)
            ELSEIF (DEPTH.EQ.2)THEN
              CALL NODETOELEMENT (ISN,NELEMENT,ISP,DEPTH,NNOUT,"LENGT")
            ENDIF
            
            ! FIND THE NUMBER OF TOP AND BOTTOM POINT 
            CALL ADDSOILSTIFF (SD,NTSOIL,NODETOPMAX,NODEBOTTOMAX,KCOUNBOTTOM,KCOUNTOP)
            
            DO JI = 1 ,KCOUNTOP
C              CALL NODETOELEMENT (NODETOPMAX(JI),NELEMENT,NISN,DEPTH,NODESECOND,"ABOVE")
                IF (ISN.EQ.NODETOPMAX(JI)) NTOPOPT = 1D0
                IF (ISP.EQ.NODETOPMAX(JI)) NTOPOPT = 1D0
            ENDDO

            CALL ELEMENTLENGTH (ISP,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,"TZ")
            
            ! BY FOLLOW THE SACS PROGRAM
            IF (SBC(1).EQ.1)THEN ! P-Y CURVE 
              ! GENERATE DATA OF LAYER SOIL LAYER 
              IF (SBC(4).EQ.1) CALL GENERGAPH (ISP,ISN,SBC,ALENGTHPP,ALENGTHQZ,POTOTAL,ALENGTHPY,DATATZSAND,DATATZCLAY
     1                                        ,DATAPYCLAY,NTOPOPT,NODETOPMAX,"P-Y")
              ! FOR P-Y CURVE CLAY 
              ! CALCULATE STIFFNESS BY USING ACTUAL DISPLACEMENT TO GET STIFFNESS 
              ! EDAT(1:3)  = ACTUAL DISPLACEMENT 
              ! EDAT(4:6)  = ACTUAL ROTATION
              IF (SBC(4).EQ.1) CALL APICURVE(EDAT(1:3),EDAT(4:6),SBC,DATATZSAND,DATATZCLAY,DATAPYCLAY,STIFF,STIFF_RO
     1                              ,DEPTH,ALENGTH,NTSOIL,"P-Y")
             
              ! FOR P-Y CURVE SAND 
              IF (SBC(4).NE.1) THEN
              CALL P_Y_CURVE_SAND  (ISN,ISP,EDAT(1:3),SBC(10),DPIPE,GRMMA,ALENGTHPY,ISAND,NCODE,NT,SBC(1:37)
     1                                              ,STIFF,PLOTP,PLOTY)
              CURVEDATA(51:65,ISN)   = PLOTP(1:15,KPCOUNT-1)
              CURVEDATA(66:80,ISN)   = PLOTY(1:15,KPCOUNT-1)
              ENDIF
              
              
              ! FOR PLOTTING GRAPH DATA 
              CURVEDATA(1:16,ISN)    = DATATZCLAY(1:16)
              CURVEDATA(21:26,ISN)   = DATATZSAND(1:6)
              CURVEDATA(31:44,ISN)   = DATAPYCLAY(1:14)
 
              ! SD(1:6,ISN)  = ADD STIFFNESS VALUE RELATE WITH NODAL(ISN)
              IF (NGRAV.EQ.1)THEN
              ! ERROR  = PERCENT ERROR BETWEEN OLD AND NEW STIFFNESS 
              ERROR(10,ISN)     =  0.0D0
              ERROR(11,ISN)     =  ABS(AUSTIFF(2,ISN) - STIFF(2))/ABS(AUSTIFF(2,ISN))*100D0
              ERROR(12,ISN)     =  ABS(AUSTIFF(3,ISN) - STIFF(3))/ABS(AUSTIFF(3,ISN))*100D0
                 IF (STIFF(2).EQ.0D0) ERROR(11,ISN) = 0.0D0
                 IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(2,ISN)      = SD(2,ISN) - AUSTIFF(2,ISN) 
                 SD(3,ISN)      = SD(3,ISN) - AUSTIFF(3,ISN) 
                 AUSTIFF(2,ISN) = STIFF(2)
                 AUSTIFF(3,ISN) = STIFF(3)
                 SD(2,ISN)      = SD(2,ISN) + STIFF(2)
                 SD(3,ISN)      = SD(3,ISN) + STIFF(3)
              IF (NSUC.EQ.1)THEN
              ! ROTATION
              ERROR(13,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(1))/ABS(ATSTIFF_RO(1,ISN))*100D0
              ERROR(14,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(2))/ABS(ATSTIFF_RO(2,ISN))*100D0
              ERROR(15,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(3))/ABS(ATSTIFF_RO(3,ISN))*100D0
                 IF (STIFF_RO(1).EQ.0D0) ERROR(13,ISN) = 0.0D0
                 IF (STIFF_RO(2).EQ.0D0) ERROR(14,ISN) = 0.0D0
                 IF (STIFF_RO(3).EQ.0D0) ERROR(15,ISN) = 0.0D0
                 SD(1,ISN)          = SD(1,ISN) - ATSTIFF_RO(1,ISN)
                 SD(2,ISN)          = SD(2,ISN) - ATSTIFF_RO(2,ISN)
                 SD(3,ISN)          = SD(3,ISN) - ATSTIFF_RO(3,ISN)
                 ATSTIFF_RO(1,ISN)  = STIFF_RO(1)
                 ATSTIFF_RO(2,ISN)  = STIFF_RO(2)
                 ATSTIFF_RO(3,ISN)  = STIFF_RO(3)
                 SD(1,ISN)          = SD(1,ISN) + STIFF_RO(1)
                 SD(2,ISN)          = SD(2,ISN) + STIFF_RO(2)
                 SD(3,ISN)          = SD(3,ISN) + STIFF_RO(3)
                 ENDIF
              ELSEIF (NGRAV.EQ.2)THEN
              ERROR(10,ISN)     =  ABS(AUSTIFF(1,ISN) - STIFF(1))/ABS(AUSTIFF(1,ISN))*100D0
              ERROR(11,ISN)     =  0.0D0
              ERROR(12,ISN)     =  ABS(AUSTIFF(3,ISN) - STIFF(3))/ABS(AUSTIFF(3,ISN))*100D0
                 IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
                 IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(1,ISN)      = SD(1,ISN) - AUSTIFF(1,ISN)
                 SD(3,ISN)      = SD(3,ISN) - AUSTIFF(3,ISN)
                 AUSTIFF(1,ISN) = STIFF(1)
                 AUSTIFF(3,ISN) = STIFF(3)
                 SD(1,ISN)      = SD(1,ISN) + STIFF(1)
                 SD(3,ISN)      = SD(3,ISN) + STIFF(3)
              IF (NSUC.EQ.1)THEN
              ! ROTATION
              ERROR(13,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(1))/ABS(ATSTIFF_RO(1,ISN))*100D0
              ERROR(14,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(2))/ABS(ATSTIFF_RO(2,ISN))*100D0
              ERROR(15,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(3))/ABS(ATSTIFF_RO(3,ISN))*100D0
                 IF (STIFF_RO(1).EQ.0D0) ERROR(13,ISN) = 0.0D0
                 IF (STIFF_RO(2).EQ.0D0) ERROR(14,ISN) = 0.0D0
                 IF (STIFF_RO(3).EQ.0D0) ERROR(15,ISN) = 0.0D0
                 SD(1,ISN)          = SD(1,ISN) - ATSTIFF_RO(1,ISN)
                 SD(2,ISN)          = SD(2,ISN) - ATSTIFF_RO(2,ISN)
                 SD(3,ISN)          = SD(3,ISN) - ATSTIFF_RO(3,ISN)
                 ATSTIFF_RO(1,ISN)  = STIFF_RO(1)
                 ATSTIFF_RO(2,ISN)  = STIFF_RO(2)
                 ATSTIFF_RO(3,ISN)  = STIFF_RO(3)
                 SD(1,ISN)          = SD(1,ISN) + STIFF_RO(1)
                 SD(2,ISN)          = SD(2,ISN) + STIFF_RO(2)
                 SD(3,ISN)          = SD(3,ISN) + STIFF_RO(3)
                 ENDIF
              ELSEIF (NGRAV.EQ.3)THEN
              ERROR(10,ISN)     =  ABS(AUSTIFF(1,ISN) - STIFF(1))/ABS(AUSTIFF(1,ISN))*100D0
              ERROR(11,ISN)     =  ABS(AUSTIFF(2,ISN) - STIFF(2))/ABS(AUSTIFF(2,ISN))*100D0
              ERROR(12,ISN)     =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
                 IF (STIFF(2).EQ.0D0) ERROR(11,ISN) = 0.0D0
                 SD(1,ISN)      = SD(1,ISN) - AUSTIFF(1,ISN)
                 SD(2,ISN)      = SD(2,ISN) - AUSTIFF(2,ISN) 
                 AUSTIFF(1,ISN) = STIFF(1)
                 AUSTIFF(2,ISN) = STIFF(2)
                 SD(1,ISN)      = SD(1,ISN) + STIFF(1)
                 SD(2,ISN)      = SD(2,ISN) + STIFF(2)
              IF (NSUC.EQ.1)THEN
              ! ROTATION
              ERROR(13,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(1))/ABS(ATSTIFF_RO(1,ISN))*100D0
              ERROR(14,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(2))/ABS(ATSTIFF_RO(2,ISN))*100D0
              ERROR(15,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(3))/ABS(ATSTIFF_RO(3,ISN))*100D0
                 IF (STIFF_RO(1).EQ.0D0) ERROR(13,ISN) = 0.0D0
                 IF (STIFF_RO(2).EQ.0D0) ERROR(14,ISN) = 0.0D0
                 IF (STIFF_RO(3).EQ.0D0) ERROR(15,ISN) = 0.0D0
                 SD(1,ISN)          = SD(1,ISN) - ATSTIFF_RO(1,ISN)
                 SD(2,ISN)          = SD(2,ISN) - ATSTIFF_RO(2,ISN)
                 SD(3,ISN)          = SD(3,ISN) - ATSTIFF_RO(3,ISN)
                 ATSTIFF_RO(1,ISN)  = STIFF_RO(1)
                 ATSTIFF_RO(2,ISN)  = STIFF_RO(2)
                 ATSTIFF_RO(3,ISN)  = STIFF_RO(3)
                 SD(1,ISN)          = SD(1,ISN) + STIFF_RO(1)
                 SD(2,ISN)          = SD(2,ISN) + STIFF_RO(2)
                 SD(3,ISN)          = SD(3,ISN) + STIFF_RO(3)
                 ENDIF
              ENDIF
            ENDIF
            
            IF (SBC(2).EQ.1)THEN ! T-Z CURVE
                
              CALL GENERGAPH (ISP,ISN,SBC,ALENGTH,ALENGTH,POTOTAL,ALENGTHPY,DATATZSAND,DATATZCLAY,DATAPYCLAY,NTOPOPT
     1                       ,NODETOPMAX,"T-Z")
              
              CALL APICURVE(EDAT(1:3),EDAT(4:6),SBC,DATATZSAND,DATATZCLAY,DATAPYCLAY,STIFF,STIFF_RO,DEPTH,ALENGTH,NTSOIL,"T-Z")
              
              CURVEDATA(1:16,ISN)    = DATATZCLAY(1:16)
              CURVEDATA(21:26,ISN)   = DATATZSAND(1:6)
              CURVEDATA(31:44,ISN)   = DATAPYCLAY(1:14)
              
              IF (NGRAV.EQ.1)THEN
              ! DISPLACEMENT
              ERROR(10,ISN)      =  ABS(ATSTIFF(1,ISN) - STIFF(1))/ABS(ATSTIFF(1,ISN))*100D0
              ERROR(11,ISN)      =  0.0D0
              ERROR(12,ISN)      =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(10,ISN) = 0.0D0
                 SD(1,ISN)       = SD(1,ISN) - ATSTIFF(1,ISN)
                 ATSTIFF(1,ISN)  = STIFF(1)
                 SD(1,ISN)       = SD(1,ISN) + STIFF(1)
              !IF (NSUC.EQ.1)THEN
              !! ROTATION
              !ERROR(13,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(1))/ABS(ATSTIFF_RO(1,ISN))*100D0
              !ERROR(14,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(2))/ABS(ATSTIFF_RO(2,ISN))*100D0
              !ERROR(15,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(3))/ABS(ATSTIFF_RO(3,ISN))*100D0
              !   IF (STIFF_RO(1).EQ.0D0) ERROR(13,ISN) = 0.0D0
              !   IF (STIFF_RO(2).EQ.0D0) ERROR(14,ISN) = 0.0D0
              !   IF (STIFF_RO(3).EQ.0D0) ERROR(15,ISN) = 0.0D0
              !   SD(1,ISN)          = SD(1,ISN) - ATSTIFF_RO(1,ISN)
              !   SD(2,ISN)          = SD(2,ISN) - ATSTIFF_RO(2,ISN)
              !   SD(3,ISN)          = SD(3,ISN) - ATSTIFF_RO(3,ISN)
              !   ATSTIFF_RO(1,ISN)  = STIFF_RO(1)
              !   ATSTIFF_RO(2,ISN)  = STIFF_RO(2)
              !   ATSTIFF_RO(3,ISN)  = STIFF_RO(3)
              !   SD(1,ISN)          = SD(1,ISN) + STIFF_RO(1)
              !   SD(2,ISN)          = SD(2,ISN) + STIFF_RO(2)
              !   SD(3,ISN)          = SD(3,ISN) + STIFF_RO(3)
              !ENDIF
              ELSEIF (NGRAV.EQ.2)THEN
              ! DISPLACEMENT
              ERROR(10,ISN)      =  0.0D0
              ERROR(12,ISN)      =  0.0D0
              ERROR(11,ISN)      =  ABS(ATSTIFF(2,ISN) - STIFF(2))/ABS(ATSTIFF(2,ISN))*100D0
                 IF (STIFF(2).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(2,ISN)       = SD(2,ISN) - ATSTIFF(2,ISN) 
                 ATSTIFF(2,ISN)  = STIFF(2)
                 SD(2,ISN)       = SD(2,ISN) + STIFF(2)
              !IF (NSUC.EQ.1)THEN
              !! ROTATION
              !ERROR(13,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(1))/ABS(ATSTIFF_RO(1,ISN))*100D0
              !ERROR(14,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(2))/ABS(ATSTIFF_RO(2,ISN))*100D0
              !ERROR(15,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(3))/ABS(ATSTIFF_RO(3,ISN))*100D0
              !   IF (STIFF_RO(1).EQ.0D0) ERROR(13,ISN) = 0.0D0
              !   IF (STIFF_RO(2).EQ.0D0) ERROR(14,ISN) = 0.0D0
              !   IF (STIFF_RO(3).EQ.0D0) ERROR(15,ISN) = 0.0D0
              !   SD(1,ISN)          = SD(1,ISN) - ATSTIFF_RO(1,ISN)
              !   SD(2,ISN)          = SD(2,ISN) - ATSTIFF_RO(2,ISN)
              !   SD(3,ISN)          = SD(3,ISN) - ATSTIFF_RO(3,ISN)
              !   ATSTIFF_RO(1,ISN)  = STIFF_RO(1)
              !   ATSTIFF_RO(2,ISN)  = STIFF_RO(2)
              !   ATSTIFF_RO(3,ISN)  = STIFF_RO(3)
              !   SD(1,ISN)          = SD(1,ISN) + STIFF_RO(1)
              !   SD(2,ISN)          = SD(2,ISN) + STIFF_RO(2)
              !   SD(3,ISN)          = SD(3,ISN) + STIFF_RO(3)
              !ENDIF
              ELSEIF (NGRAV.EQ.3)THEN
              ! DISPLACEMENT
              ERROR(10,ISN)      =  0.0D0
              ERROR(11,ISN)      =  0.0D0
              ERROR(12,ISN)      =  ABS(ATSTIFF(3,ISN) - STIFF(3))/ABS(ATSTIFF(3,ISN))*100D0
                 IF (STIFF(3).EQ.0D0) ERROR(12,ISN) = 0.0D0
                 SD(3,ISN)       = SD(3,ISN) - ATSTIFF(3,ISN) 
                 ATSTIFF(3,ISN)  = STIFF(3)
                 SD(3,ISN)       = SD(3,ISN) + STIFF(3)
                ! IF (ILC.EQ.1) SD(3,ISN)  = 10000000000D0
              !IF (NSUC.EQ.1)THEN
              !! ROTATION
              !ERROR(13,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(1))/ABS(ATSTIFF_RO(1,ISN))*100D0
              !ERROR(14,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(2))/ABS(ATSTIFF_RO(2,ISN))*100D0
              !ERROR(15,ISN)      =  ABS(ATSTIFF_RO(1,ISN) - STIFF_RO(3))/ABS(ATSTIFF_RO(3,ISN))*100D0
              !   IF (STIFF_RO(1).EQ.0D0) ERROR(13,ISN) = 0.0D0
              !   IF (STIFF_RO(2).EQ.0D0) ERROR(14,ISN) = 0.0D0
              !   IF (STIFF_RO(3).EQ.0D0) ERROR(15,ISN) = 0.0D0
              !   SD(1,ISN)          = SD(1,ISN) - ATSTIFF_RO(1,ISN)
              !   SD(2,ISN)          = SD(2,ISN) - ATSTIFF_RO(2,ISN)
              !   SD(3,ISN)          = SD(3,ISN) - ATSTIFF_RO(3,ISN)
              !   ATSTIFF_RO(1,ISN)  = STIFF_RO(1)
              !   ATSTIFF_RO(2,ISN)  = STIFF_RO(2)
              !   ATSTIFF_RO(3,ISN)  = STIFF_RO(3)
              !   SD(1,ISN)          = SD(1,ISN) + STIFF_RO(1)
              !   SD(2,ISN)          = SD(2,ISN) + STIFF_RO(2)
              !   SD(3,ISN)          = SD(3,ISN) + STIFF_RO(3)
              !ENDIF
              ENDIF
               ! START Q-Z CURVE
              
C              CALL ADDSOILSTIFF (SD,NTSOIL,NODETOPMAX,NODEBOTTOMAX,KCOUNBOTTOM,KCOUNTOP)
              
C              DO IJ = 1,KCOUNBOTTOM
              
C              CALL OVERBUND (POTOTAL)
                
C              CALL GENERGAPH (ISP,ISN,SBC,ALENGTH,POTOTAL,ALENGTHPY,DATATZSAND,DATATZCLAY,DATAPYCLAY,NTOPOPT,"Q-Z")
              
C              CALL APICURVE(EDAT(1:3),EDAT(4:6),SBC,DATATZSAND,DATATZCLAY,DATAPYCLAY,STIFF,DEPTH,ALENGTH,NTSOIL,NODETOPMAX,"Q-Z")

C              IF (NGRAV.EQ.1)THEN
C              ERROR(11,NODEBOTTOMAX(IJ))      =  ABS(AQSTIFF(1,NODEBOTTOMAX(IJ)) - STIFF(1))/ABS(AQSTIFF(1,NODEBOTTOMAX(IJ)))*100D0
C              ERROR(12,NODEBOTTOMAX(IJ))      =  0.0D0
C              ERROR(10,NODEBOTTOMAX(IJ))      =  0.0D0
C                 IF (STIFF(1).EQ.0D0) ERROR(10,NODEBOTTOMAX(IJ)) = 0.0D0
C                 SD(1,NODEBOTTOMAX(IJ)) = SD(1,NODEBOTTOMAX(IJ)) - AQSTIFF(1,NODEBOTTOMAX(IJ))
C                 AQSTIFF(1,NODEBOTTOMAX(IJ))  = STIFF(1)
C                 SD(1,NODEBOTTOMAX(IJ))       = SD(1,NODEBOTTOMAX(IJ)) + STIFF(1)
C              ELSEIF (NGRAV.EQ.2)THEN
C              ERROR(10,NODEBOTTOMAX(IJ))      =  0.0D0
C              ERROR(12,NODEBOTTOMAX(IJ))      =  0.0D0
C              ERROR(11,NODEBOTTOMAX(IJ))      =  ABS(AQSTIFF(2,NODEBOTTOMAX(IJ)) - STIFF(2))/ABS(AQSTIFF(2,NODEBOTTOMAX(IJ)))*100D0
C                 IF (STIFF(2).EQ.0D0) ERROR(11,NODEBOTTOMAX(IJ)) = 0.0D0
C                 SD(2,NODEBOTTOMAX(IJ))       = SD(2,NODEBOTTOMAX(IJ)) - AQSTIFF(2,NODEBOTTOMAX(IJ)) 
C                 AQSTIFF(2,NODEBOTTOMAX(IJ))  = STIFF(2)
C                 SD(2,NODEBOTTOMAX(IJ))       = SD(2,NODEBOTTOMAX(IJ)) + STIFF(2)
C              ELSEIF (NGRAVE.EQ.3)THEN
C             ERROR(10,NODEBOTTOMAX(IJ))      =  0.0D0
C              ERROR(11,NODEBOTTOMAX(IJ))      =  0.0D0
C              ERROR(12,NODEBOTTOMAX(IJ))      =  ABS(AQSTIFF(3,NODEBOTTOMAX(IJ)) - STIFF(3))/ABS(AQSTIFF(3,NODEBOTTOMAX(IJ)))*100D0
C                 IF (STIFF(3).EQ.0D0) ERROR(12,NODEBOTTOMAX(IJ)) = 0.0D0
C                 SD(3,NODEBOTTOMAX(IJ))       = SD(3,NODEBOTTOMAX(IJ)) - AQSTIFF(3,NODEBOTTOMAX(IJ)) 
C                 AQSTIFF(3,NODEBOTTOMAX(IJ))  = STIFF(3)
C                 SD(3,NODEBOTTOMAX(IJ))       = SD(3,ISN) + STIFF(3)
C              ENDIF
              
C              ENDDO
            ENDIF
            
            
C         ORIGINAL SOIL PILE INTERACTION ( **** DO NOT USED **** )
          ! ORIGINAL P-Y CURVE SOLVER WITHOUT DIVIDE LENGTH ***** SKIP *****
          IF(SBC(1).EQ.10.0D0)THEN
C         ------------- P - Y CURVE -------------
          ! 1 = SOFT CLAY
          ! 2 = STFF CLAY
          ! 3 = SAND
          STIFF = 0.
          IF (SBC(4).EQ.1) THEN
             AJ = 0.5D0
             CALL P_Y_CURVE_SOFT_CLAY (EDAT(1:3),C,DPIPE,GRMMA,AJ,ALENGTH,STRAIN,NCODE,NT,SBC(1:37),STIFF)

             WRITE(170,170) STIFF(1),EDAT(1)
170          FORMAT (E12.5,E12.5)    
             GOTO 500
          ENDIF
      
          IF (SBC(4).EQ.2) THEN
             AJ = 0.25D0
             CALL P_Y_CURVE_SOFT_CLAY (EDAT(1:3),C,DPIPE,GRMMA,AJ,ALENGTH,STRAIN,NCODE,NT,SBC(1:37),STIFF)
             GOTO 500
          ENDIF
          
          IF (SBC(4).EQ.3) THEN
             GRMMA  = SBC(8)  !20D0
             CALL P_Y_CURVE_SAND (ISN,ISP,EDAT(1:3),FREE,DPIPE,GRMMA,ALENGTH,ISAND,NCODE,NT,SBC(1:37),STIFF,PLOTP,PLOTY)
             GOTO 500
          ENDIF
          
          
C         UPDATE STIFFNESS
500        IF (SBC(1).EQ.1)THEN ! P-Y CURVE
              IF (NGRAV.EQ.1)THEN
              ERROR(1,ISN)     =  0.0D0
              ERROR(2,ISN)     =  ABS(AUSTIFF(2,ISN) - STIFF(2))/ABS(AUSTIFF(2,ISN))*100D0
              ERROR(3,ISN)     =  ABS(AUSTIFF(3,ISN) - STIFF(3))/ABS(AUSTIFF(3,ISN))*100D0
                 IF (STIFF(2).EQ.0D0) ERROR(2,ISN) = 0.0D0
                 IF (STIFF(3).EQ.0D0) ERROR(3,ISN) = 0.0D0
                 SD(2,ISN)      = SD(2,ISN) - AUSTIFF(2,ISN) 
                 SD(3,ISN)      = SD(3,ISN) - AUSTIFF(3,ISN) 
                 AUSTIFF(2,ISN) = STIFF(2)
                 AUSTIFF(3,ISN) = STIFF(3)
                 SD(2,ISN)      = SD(2,ISN) + STIFF(2)
                 SD(3,ISN)      = SD(3,ISN) + STIFF(3)
              ELSEIF (NGRAV.EQ.2)THEN
              ERROR(1,ISN)     =  ABS(AUSTIFF(1,ISN) - STIFF(1))/ABS(AUSTIFF(1,ISN))*100D0
              ERROR(2,ISN)     =  0.0D0
              ERROR(3,ISN)     =  ABS(AUSTIFF(3,ISN) - STIFF(3))/ABS(AUSTIFF(3,ISN))*100D0
                 IF (STIFF(1).EQ.0D0) ERROR(1,ISN) = 0.0D0
                 IF (STIFF(3).EQ.0D0) ERROR(3,ISN) = 0.0D0
                 SD(1,ISN)      = SD(1,ISN) - AUSTIFF(1,ISN)
                 SD(3,ISN)      = SD(3,ISN) - AUSTIFF(3,ISN)
                 AUSTIFF(1,ISN) = STIFF(1)
                 AUSTIFF(3,ISN) = STIFF(3)
                 SD(1,ISN)      = SD(1,ISN) + STIFF(1)
                 SD(3,ISN)      = SD(3,ISN) + STIFF(3)
              ELSEIF (NGRAV.EQ.3)THEN
              ERROR(1,ISN)     =  ABS(AUSTIFF(1,ISN) - STIFF(1))/ABS(AUSTIFF(1,ISN))*100D0
              ERROR(2,ISN)     =  ABS(AUSTIFF(2,ISN) - STIFF(2))/ABS(AUSTIFF(2,ISN))*100D0
              ERROR(3,ISN)     =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(1,ISN) = 0.0D0
                 IF (STIFF(2).EQ.0D0) ERROR(2,ISN) = 0.0D0
                 SD(1,ISN)      = SD(1,ISN) - AUSTIFF(1,ISN)
                 SD(2,ISN)      = SD(2,ISN) - AUSTIFF(2,ISN) 
                 AUSTIFF(1,ISN) = STIFF(1)
                 AUSTIFF(2,ISN) = STIFF(2)
                 SD(1,ISN)      = SD(1,ISN) + STIFF(1)
                 SD(2,ISN)      = SD(2,ISN) + STIFF(2)
              ENDIF
          ENDIF
          
C          IF (ISN.EQ.8)THEN
C          WRITE(7101,'(E12.5,E12.5,E12.5,E12.5)') EDAT(1),STIFF(1),EDAT(3),STIFF(3)
C          ENDIF
          
          IF(ICOUNT.EQ.1)THEN
          EDAT(1:3) = 1D0
          ENDIF
                    
          IF (NGRAV.EQ.1)THEN
          
          ELSEIF (NGRAV.EQ.2)THEN
              
          ELSEIF (NGRAV.EQ.3)THEN
              
          ENDIF
          
C         ------------- P - Y CURVE -------------
          ENDIF
          
          IF (SBC(2).EQ.10.0D0)THEN
C         ------------- T - Z CURVE -------------          
          ! 1 = CLAY
          ! 2 = SAND
C          C     = SBC(10)
C          AK    = SBC(11)
C          PO    = SBC(12)
C          DELTA = SBC(13)
C          STIFF = 0.
C          IF (SBC(9).EQ.1) CALL T_Z_CURVE_CLAY (ALENGTH,C,PO,DPIPE,EDAT(1:3),NCODE,NT,SBC(1:35),STIFF)
C          IF (SBC(9).EQ.2) CALL T_Z_CURVE_SAND (ALENGTH,PO,DPIPE,AK,DELTA,EDAT(1:3),NCODE,NT,SBC(1:35),STIFF)
          
C         UPDATE STIFFNESS
          IF (SBC(2).EQ.1)THEN ! T-Z CURVE
              IF (NGRAV.EQ.1)THEN
              ERROR(4,ISN)      =  ABS(ATSTIFF(1,ISN) - STIFF(1))/ABS(ATSTIFF(1,ISN))*100D0
              ERROR(5,ISN)      =  0.0D0
              ERROR(6,ISN)      =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(4,ISN) = 0.0D0
                 SD(1,ISN) = SD(1,ISN) - ATSTIFF(1,ISN)
                 ATSTIFF(1,ISN)  = STIFF(1)
                 SD(1,ISN)       = SD(1,ISN) + STIFF(1)
              ELSEIF (NGRAV.EQ.2)THEN
              ERROR(4,ISN)      =  0.0D0
              ERROR(5,ISN)      =  ABS(ATSTIFF(2,ISN) - STIFF(2))/ABS(ATSTIFF(2,ISN))*100D0
              ERROR(6,ISN)      =  0.0D0
                 IF (STIFF(2).EQ.0D0) ERROR(5,ISN) = 0.0D0
                 SD(2,ISN)       = SD(2,ISN) - ATSTIFF(2,ISN) 
                 ATSTIFF(2,ISN)  = STIFF(2)
                 SD(2,ISN)       = SD(2,ISN) + STIFF(2)
              ELSEIF (NGRAV.EQ.3)THEN
              ERROR(4,ISN)      =  0.0D0
              ERROR(5,ISN)      =  0.0D0
              ERROR(6,ISN)      =  ABS(ATSTIFF(3,ISN) - STIFF(3))/ABS(ATSTIFF(3,ISN))*100D0
                 IF (STIFF(3).EQ.0D0) ERROR(6,ISN) = 0.0D0
                 SD(3,ISN)       = SD(3,ISN) - ATSTIFF(3,ISN) 
                 ATSTIFF(3,ISN)  = STIFF(3)
                 SD(3,ISN)       = SD(3,ISN) + STIFF(3)
              ENDIF
          ENDIF
C         ------------- T - Z CURVE ------------- 
          ENDIF
          
          IF (SBC(3).EQ.10.0D0)THEN
C         ------------- Q - Z CURVE -------------       
          STIFF = 0.
          ANQ  =  SBC(14)
          C    =  SBC(15)
          PO   =  SBC(17)
          CALL Q_Z_CURVE (ANQ,PO,C,DPIPE,EDAT(1:3),NCODE,NT,SBC(1:37),STIFF,NOPT)
          
C         UPDATE STIFFNESS
          IF (SBC(3).EQ.1)THEN ! Q-Z CURVE
              IF (NGRAV.EQ.1)THEN
              ERROR(7,ISN)      =  ABS(AQSTIFF(1,ISN) - STIFF(1))/ABS(AQSTIFF(1,ISN))*100D0
              ERROR(8,ISN)      =  0.0D0
              ERROR(9,ISN)      =  0.0D0
                 IF (STIFF(1).EQ.0D0) ERROR(7,ISN) = 0.0D0
                 SD(1,ISN) = SD(1,ISN) - AQSTIFF(1,ISN)
                 AQSTIFF(1,ISN)  = STIFF(1)
                 SD(1,ISN)       = SD(1,ISN) + STIFF(1)
              ELSEIF (NGRAV.EQ.2)THEN
              ERROR(7,ISN)      =  0.0D0
              ERROR(8,ISN)      =  ABS(AQSTIFF(2,ISN) - STIFF(2))/ABS(AQSTIFF(2,ISN))*100D0
              ERROR(9,ISN)      =  0.0D0
                 IF (STIFF(2).EQ.0D0) ERROR(8,ISN) = 0.0D0
                 SD(2,ISN)       = SD(2,ISN) - AQSTIFF(2,ISN) 
                 AQSTIFF(2,ISN)  = STIFF(2)
                 SD(2,ISN)       = SD(2,ISN) + STIFF(2)
              ELSEIF (NGRAV.EQ.3)THEN
              ERROR(7,ISN)      =  0.0D0
              ERROR(8,ISN)      =  0.0D0
              ERROR(9,ISN)      =  ABS(AQSTIFF(3,ISN) - STIFF(3))/ABS(AQSTIFF(3,ISN))*100D0
                 IF (STIFF(3).EQ.0D0) ERROR(9,ISN) = 0.0D0
                 SD(3,ISN)       = SD(3,ISN) - AQSTIFF(3,ISN) 
                 AQSTIFF(3,ISN)  = STIFF(3)
                 SD(3,ISN)       = SD(3,ISN) + STIFF(3)
              ENDIF
          ENDIF
C         ------------- Q - Z CURVE -------------   
C          WRITE (7100,*) STIFF(2)
          ENDIF

          WRITE(7101,'(I5,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5)') ISN,STIFF(1),STIFF(2),STIFF(3)
     1                                                          ,EDAT(1),EDAT(2),EDAT(3)
         
100     ENDDO
               ! ---------------------------- DO NOT USED------------------------------
   
         
         
         
         
         
         
   
          ! FIND THE BOTTOM OF NODE AT EACH PILE SECTION
          CALL ADDSOILSTIFF (SD,NTSOIL,NODETOPMAX,NODEBOTTOMAX,KCOUNBOTTOM,KCOUNTOP)
          
          DO IJ = 1 ,KCOUNBOTTOM
          CALL NODETOELEMENT (NODEBOTTOMAX(IJ),NELEMENT,NISN,DEPTH,NODESECOND,"ABOVE")
          CALL READSOILDATA ("CALE",NODEBOTTOMAX(IJ),DFAC,NCODE,NT,SBC(1:37),NODESECOND,ISN,PROPOUT)

          IF (NCODE.EQ.3.AND.NCODE.EQ.4)THEN
          DO I = 1,16
          X(I) = SBC(I+3)    
          ENDDO
          ELSEIF (NCODE.EQ.1)THEN
           IF (SBC(1).EQ.1)THEN
               IF (SBC(4).EQ.1)THEN
                   ! DO NOT ADD STIFFNESS
               ELSEIF (SBC(4).NE.1)THEN
                   ! DO NOT ADD STIFFNESS
               ENDIF
           ELSEIF (SBC(2).EQ.1)THEN
               IF (SBC(15).EQ.1)THEN
               ! CLAY
               X(1:8) = CURVEDATA(1:8,NODESECOND)
               ELSEIF (SBC(15).EQ.2) THEN
               ! SAND   
               X(1:3) = CURVEDATA(21:23,NODESECOND)
               ENDIF
           ENDIF   
          ENDIF

C          DISP1 =  ABS(DISP((NODEBOTTOMAX(IJ)*6)-5))
C          DISP2 =  ABS(DISP((NODEBOTTOMAX(IJ)*6)-4))
C          DISP3 =  ABS(DISP((NODEBOTTOMAX(IJ)*6)-3))
           IF (NCODE.EQ.1)THEN
           EDAT(1) =  ABS(DISP((NODEBOTTOMAX(IJ)*6)-5))
           EDAT(2) =  ABS(DISP((NODEBOTTOMAX(IJ)*6)-4))
           EDAT(3) =  ABS(DISP((NODEBOTTOMAX(IJ)*6)-3))
          
C          IF (SBC(1).EQ.1)THEN
C          CALL APICURVE(EDAT(1:3),EDAT(4:6),SBC,CURVEDATA(21:23,NODESECOND),CURVEDATA(1:16,NODESECOND),DATAPYCLAY,
C     1                STIFF,DEPTH,ALENGTH,NTSOIL,"P-Y")    
C          ELSEIF (SBC(2).EQ.1)THEN
          CALL APICURVE(EDAT(1:3),EDAT(4:6),SBC,CURVEDATA(21:23,NODESECOND),CURVEDATA(1:16,NODESECOND),DATAPYCLAY,
     1                STIFF,STIFF_RO,DEPTH,ALENGTH,NTSOIL,"T-Z")
C          ENDIF
         
         SD(1,NODEBOTTOMAX(IJ)) = STIFF(1)/2.4D0
         SD(2,NODEBOTTOMAX(IJ)) = STIFF(2)/2.4D0
         SD(3,NODEBOTTOMAX(IJ)) = STIFF(3)/2.4D0
           ELSEIF (NCODE.EQ.3)THEN
               
               IF (DISP1.GE.MAXVAL(X)) SD(1,NODEBOTTOMAX(IJ)) = SD(1,NODESECOND)/2.45D0
               IF (DISP2.GE.MAXVAL(X)) SD(2,NODEBOTTOMAX(IJ)) = SD(2,NODESECOND)/2.45D0
               IF (DISP3.GE.MAXVAL(X)) SD(3,NODEBOTTOMAX(IJ)) = SD(3,NODESECOND)/2.45D0
               
               IF (DISP1.LT.MAXVAL(X)) SD(1,NODEBOTTOMAX(IJ)) = SD(1,NODESECOND)/2.0D0 
               IF (DISP2.LT.MAXVAL(X)) SD(2,NODEBOTTOMAX(IJ)) = SD(2,NODESECOND)/2.0D0 
               IF (DISP3.LT.MAXVAL(X)) SD(3,NODEBOTTOMAX(IJ)) = SD(3,NODESECOND)/2.0D0
           ENDIF 
         
          
          ENDDO
         ! FIND THE TOP OF NODE AT EACH PILE SECTION
          DO IK = 1,KCOUNTOP
          CALL NODETOELEMENT (NODETOPMAX(IK),NELEMENT,NISN,DEPTH,NODESECOND,"BOTTO") 
          CALL READSOILDATA ("CALE",NODETOPMAX(IK),DFAC,NCODE,NT,SBC(1:37),NODESECOND,ISN,PROPOUT)
          
          IF (NCODE.EQ.3)THEN
          DO I = 1,16
          X(I) = SBC(I+3)    
          ENDDO
          ELSEIF (NCODE.EQ.1)THEN
           IF (SBC(1).EQ.1)THEN
               IF (SBC(4).EQ.1)THEN
                   ! DO NOT ADD STIFFNESS
               ELSEIF (SBC(4).NE.1)THEN
                   ! DO NOT ADD STIFFNESS
               ENDIF
           ELSEIF (SBC(2).EQ.1)THEN
               IF (SBC(15).EQ.1)THEN
               ! CLAY
               X(1:8) = CURVEDATA(1:8,NODESECOND)
               ELSEIF (SBC(15).EQ.2) THEN
               ! SAND   
               X(1:3) = CURVEDATA(21:23,NODESECOND)
               ENDIF
           ENDIF   
          ENDIF
          
          IF (NCODE.EQ.1.AND.NCODE.EQ.4)THEN         
          EDAT(1) =  ABS(DISP((NODETOPMAX(IK)*6)-5))
          EDAT(2) =  ABS(DISP((NODETOPMAX(IK)*6)-4))
          EDAT(3) =  ABS(DISP((NODETOPMAX(IK)*6)-3))
          
C          IF (SBC(1).EQ.1)THEN
C          CALL APICURVE(EDAT(1:3),EDAT(4:6),SBC,CURVEDATA(21:23,NODESECOND),CURVEDATA(1:16,NODESECOND),DATAPYCLAY,
C     1                STIFF,DEPTH,ALENGTH,NTSOIL,"P-Y")    
C          ELSEIF (SBC(2).EQ.1)THEN
          CALL APICURVE(EDAT(1:3),EDAT(4:6),SBC,CURVEDATA(21:23,NODESECOND),CURVEDATA(1:16,NODESECOND),DATAPYCLAY,
     1                STIFF,STIFF_RO,DEPTH,ALENGTH,NTSOIL,"T-Z")
C          ENDIF
         
         !SD(1,NODETOPMAX(IK)) = STIFF(1)/2D0
         !SD(2,NODETOPMAX(IK)) = STIFF(2)/2D0
         !SD(3,NODETOPMAX(IK)) = STIFF(3)/2D0
         
         ELSEIF (NCODE.EQ.3)THEN
          DISP1 =  ABS(DISP((NODETOPMAX(IK)*6)-5))
          DISP2 =  ABS(DISP((NODETOPMAX(IK)*6)-4))
          DISP3 =  ABS(DISP((NODETOPMAX(IK)*6)-3))
          IF (DISP1.GE.MAXVAL(X)) SD(1,NODETOPMAX(IK)) = SD(1,NODESECOND)/2.45D0
          IF (DISP2.GE.MAXVAL(X)) SD(2,NODETOPMAX(IK)) = SD(2,NODESECOND)/2.45D0
          IF (DISP3.GE.MAXVAL(X)) SD(3,NODETOPMAX(IK)) = SD(3,NODESECOND)/2.45D0
          
          IF (DISP1.LT.MAXVAL(X)) SD(1,NODETOPMAX(IK)) = SD(1,NODESECOND)/2.0D0 
          IF (DISP2.LT.MAXVAL(X)) SD(2,NODETOPMAX(IK)) = SD(2,NODESECOND)/2.0D0 
          IF (DISP3.LT.MAXVAL(X)) SD(3,NODETOPMAX(IK)) = SD(3,NODESECOND)/2.0D0
         ENDIF
         
          
          ENDDO

       NCOUNT   = 0.D0
       ERRORMAX = MAXVAL(ERROR)
       !IF (ERRORMAX.LT.0.01) THEN ! 2018
       IF (ICOUNT.EQ.20D0) THEN ! ERROR
           ! STORAGE SOIL STIFFNESS FOR PRINTTING
           SSTIFF(1:6,1:NSN) = SD(1:6,1:NSN)

       ! P-Y CURVE FOR SAND    
       DO KB = 1,KPCOUNT
                       WRITE (7102,700)
                       !WRITE (7120,700)
                       WRITE (7102,701)
                       WRITE (7102,702)
                       WRITE (7102,703)   
                       !DATAPYSAND(1:15)  = PLOTP(I:15,KB)
                       !DATAPYSAND(16:30) = PLOTY(I:15,KB)
                       DO I =1,15
                       WRITE (7102,704) PLOTP(I,KB),PLOTY(I,KB)
                       !WRITE (7120,704) PLOTP(I,KB),PLOTY(I,KB)
                       ENDDO
                       WRITE (7102,705)
                       WRITE (7102,706)
                       WRITE (7102,706)
700                    FORMAT ('# Graph: "P-Y SAND"')
701                    FORMAT ("#")      
702                    FORMAT ('# X: "Z" Y: "T"')
703                    FORMAT ("# Units: - -")  
704                    FORMAT (F12.5,2X,F12.5)
705                    FORMAT ("#END")
706                    FORMAT ("")    
                       
       ENDDO
       
       NSOILRETURN = 1
       IF (NSOILRETURN.EQ.1) CALL SOILDATA_CURVE (ISN,ISP,DATAPYCLAY,SBC,"","PRIN") 
       ENDIF

      END
C     ----------------------------------------------------------------
      SUBROUTINE INTERNALFORCE (IGM,EDATA,NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,IPT,NUM,PDATA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 METHOD
      DIMENSION EDATA(1),PDATA(20),MAPP(1)
      
	   DO J = 1,NV
         PDATA(J) = 0.0D0
	   M        = MAPP(J)

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

      
      END
C     ----------------------------------------------------------------
      SUBROUTINE P_Y_CURVE_SOFT_CLAY (YT,C,D,GRMMA,AJ,X,STRAIN,NCODE,NT,SBC,STIFF)
C      (P,PU,C,D,GRMMA,AJ,X,STRAIN,STIFF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION STIFF(3),YT(3)
      DIMENSION PPU(7),YYC(7)
      DIMENSION SBC(1)
      DATA PPU/0.00D0,0.23D0,0.33D0,0.50D0,0.72D0,1.00D0,1.00D0/
      DATA YYC/0.00D0,0.10D0,0.30D0,1.00D0,3.00D0,8.00D0,8.00D0/
      
      COMMON /MGRAV/ NGRAV  
      
      
      COMMON/COUNT/ICOUNT,NCOUNT
      
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
     1                        ,SECTIONT,SECTIONS,BJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      

C     ---- INPUT DATA ----
C     C       = UNDRAINED SHEAR STRENGTH FOR UNDISTURBED CLAY SOIL SAMPLE 
C     D       = PILE DIAMETER 
C     GRAMMA  = EFFECTIVE UNIT WEIGHT OF SOIL 
C     AJ      = DIMENSIONLESS EMPIRICAL CONSTANT WITH VALUES RANGING FROM 0.25 TO 0.5
C     X       = DEPTH BELOW SOIL SURFACE
C     STRAIN  = ELASTIC STRAIN
      
      DO I = 1,3
          
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.1) Y = 0.0D0
         IF (I.EQ.2) Y = YT(2)
         IF (I.EQ.3) Y = YT(3)
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) Y = YT(1)
         IF (I.EQ.2) Y = 0.0D0
         IF (I.EQ.3) Y = YT(3)
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) Y = YT(1)
         IF (I.EQ.2) Y = YT(2)
         IF (I.EQ.3) Y = 0.0D0
      ENDIF
      
      
      IF (NCODE.EQ.1)THEN
C     --- API 2A  `-WSD ---
      XR = 6D0*D/((GRMMA*D/C)+AJ)
      IF (X.GT.XR)THEN
      PU = (9D0*C)*DPIPE
      ELSEIF (X.LE.XR)THEN
      PU = (3D0*C+GRMMA*X+(AJ*C*X/D))*DPIPE
      ENDIF
      
C     INTERPOLATION
      YC    = 2.5D0*STRAIN*D
      SYYC  = ABS(Y/YC)
C      SYYC  = 0.05D0

C      IF     (SYYC.GE.YYC(1).AND.SYYC.LT.YYC(2))THEN
C             AIYYC = YYC(2) - YYC(1)
C             AIPPU = PPU(2) - PPU(1)
C             UUPPU = PPU(1) + (AIPPU*(SYYC -  YYC(1))/AIYYC)
C      ELSEIF (SYYC.GE.YYC(2).AND.SYYC.LT.YYC(3))THEN
C             AIYYC = YYC(3) - YYC(2)
C             AIPPU = PPU(3) - PPU(2)
C             UUPPU = PPU(2) + (AIPPU*(SYYC -  YYC(2))/AIYYC)
C      ELSEIF (SYYC.GE.YYC(3).AND.SYYC.LT.YYC(4))THEN
C             AIYYC = YYC(4) - YYC(3)
C             AIPPU = PPU(4) - PPU(3)
C             UUPPU = PPU(3) + (AIPPU*(SYYC -  YYC(3))/AIYYC)
C      ELSEIF (SYYC.GE.YYC(4).AND.SYYC.LT.YYC(5))THEN
C             AIYYC = YYC(5) - YYC(4)
C             AIPPU = PPU(5) - PPU(4)
C             UUPPU = PPU(4) + (AIPPU*(SYYC -  YYC(4))/AIYYC)
C      ELSEIF (SYYC.GE.YYC(6).AND.SYYC.LT.YYC(7))THEN
C             AIYYC = YYC(6) - YYC(5)
C             AIPPU = PPU(6) - PPU(5)
C             UUPPU = PPU(5) + (AIPPU*(SYYC -  YYC(5))/AIYYC)
C      ELSEIF (SYYC.GE.YYC(7))THEN
C             UUPPU = PPU(7)
C      ENDIF
      UUPPU = 0.5D0*((ABS(Y/YC))**(1D0/3D0))
      IF (SYYC.GE.8.0D0) UUPPU = 1.0D0
      

      
      ENDIF
      AROUND = 2D0*3.141592654D0*(DPIPE/2D0)
C     SPRING STIFFNESS
      STIFF(I) = ABS(UUPPU*PU*X/Y)
      IF (Y.EQ.0) STIFF(I) = 0.0D0
      ENDDO
       
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE P_Y_CURVE_SAND (ISN,ISP,YT,AK,D,GRMMA,H,ISAND,NCODE,NT,SBC,STIFF,PLOTP,PLOTY)
C      (P,PU,C,D,GRMMA,AJ,X,STRAIN,STIFF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION C1(21),C2(21),C3(21)
      DIMENSION AK2(4),AK3(4),AK4(4)
      DIMENSION SBC(37)
      COMMON /MGRAV/ NGRAV  
      DATA C1/0.60D0,0.70D0,0.80D0,0.90D0,1.00D0,1.10D0,1.20D0,1.30D0,1.40D0,1.60D0,
     1        1.70D0,1.90D0,2.10D0,2.30D0,2.50D0,2.80D0,3.10D0,3.40D0,3.80D0,4.20D0,
     1        4.60D0/
      DATA C2/1.50D0,1.60D0,1.70D0,1.80D0,1.90D0,2.00D0,2.10D0,2.20D0,2.30D0,2.50D0,
     1        2.60D0,2.70D0,2.90D0,3.00D0,3.20D0,3.40D0,3.60D0,3.80D0,4.00D0,4.20D0,
     1        4.40D0/
      DATA C3/8.50D0,9.60D0,10.80D0,12.20D0,13.80D0,15.60D0,17.60D0,19.90D0,22.50D0,25.40D0,
     1        28.70D0,32.40D0,36.60D0,41.40D0,46.70D0,52.80D0,59.60D0,67.40D0,76.10D0,86.00D0, 
     1        101.50D0/
      DATA AK2/426.30D0,553.60D0,774.60D0,996.40D0/
      DATA AK3/1356.30D0,1716.10D0,2026.10D0,2491.10D0/
      DATA AK4/2850.90D0,3293.80D0,3792.00D0,4262.60D0/
      DIMENSION STIFF(3),YT(3)
      DIMENSION PLOTP(15,1000),PLOTY(15,1000)
      DIMENSION DATAPYSAND(30)
      
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,BJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      COMMON/COUNT/ICOUNT,NCOUNT
      
      COMMON/PLOTSOIL/ KPCOUNT

      IF (ICOUNT.EQ.1) YT(1:3) = 0.00001D0
      ISAND = SBC(5)
      GRMMA = SBC(6)
      AK    = SBC(10)
      XFAC  = SBC(12)
      YFAC  = SBC(13)
      
      DO I = 1,3
          
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.2) THEN
          Y = 0.0D0
          GOTO 100
         ENDIF
         IF (I.EQ.2) Y = YT(2)
         IF (I.EQ.3) THEN
          Y      = YT(3)
          PLOTYA = YT(3)
         ENDIF
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) THEN
          Y      = YT(1)
          PLOTYA = YT(1)
         ENDIF
         IF (I.EQ.2) THEN
          Y = 0.0D0
          GOTO 100
          ENDIF
         IF (I.EQ.3) Y = YT(3)
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) THEN
          Y      = YT(1)
          PLOTYA = YT(1)
         ENDIF
         IF (I.EQ.2) Y = YT(2)
         IF (I.EQ.3) THEN
          Y = 0.0D0
          GOTO 100
          ENDIF
      ENDIF

      
      
      IF (NCODE.EQ.1) THEN
      ! INTERPORATION
C      INDEX = (FREE-20)+1
      
C      AC11 = C1(INDEX)
C      AC22 = C2(INDEX)
C      AC33 = C3(INDEX)
      
C      AC12 = C1(INDEX+1)
C      AC23 = C2(INDEX+1)
C      AC34 = C3(INDEX+1)
      
C      AIV1 = AC12 - AC11
C      AIV2 = AC23 - AC22
C      AIV3 = AC34 - AC33
      
C      IVAL = FREE
      
C      AC1  = AC11 + ((FREE-IVAL))*AIV1
C      AC2  = AC22 + ((FREE-IVAL))*AIV2
C      AC3  = AC33 + ((FREE-IVAL))*AIV3

      IF (ISAND.EQ.1)THEN ! SILT 20
          AC1 = 0.6D0
          AC2 = 1.5D0
          AC3 = 8.5D0
          
      ELSEIF (ISAND.EQ.2)THEN ! SANDY SILT 25
          AC1 = 1.1D0
          AC2 = 2.0D0
          AC3 = 15.6D0
      ELSEIF (ISAND.EQ.3)THEN ! SILTY SAND  30
          AC1 = 1.7D0
          AC2 = 2.6D0
          AC3 = 28.7D0
      ELSEIF (ISAND.EQ.4)THEN ! CLEAN SAND 35
          AC1 = 3.0D0
          AC2 = 3.5D0
          AC3 = 52.8D0
          
      ELSEIF (ISAND.EQ.5)THEN ! GRAVEL 40
          AC1 = 4.6D0
          AC2 = 4.4D0
          AC3 = 101.5D0
      ENDIF
      
      
      ! LATERAL BEARING CAPACITY FOR SAND
      PUS = (AC1*H+AC2*D)*GRMMA*H
      PUD = AC3*D*GRMMA*H
      
      IF (PUS.GT.PUD) THEN
      PU  =  PUD
      A   = 3.0D0-(0.8D0*H/D) ! FOR STACTIC LOADING
      IF (A.GE.0.9D0) A = 0.92D0
      ENDIF
      IF (PUS.LE.PUD) THEN
      PU =  PUS!*H
      A   = 3.0D0-(0.8D0*H/D) ! FOR STACTIC LOADING
      IF (A.GE.0.9D0) A = 0.9D0
      ENDIF
      
      
C      ISAND = 4
C      SELECTCASE(ISAND)
C      CASE (1) ! VERY LOOSE 
C      AK  = 265.7D0     
C      CASE (2) ! LOOSE SAND
C      ! INTERPORATION
C      AINDEX = FREE-29D0
C          IF     (AINDEX.GE.0.0D0.AND.AINDEX.LT.0.334D0)THEN
C            AK = AK2(1) + (AK2(2) - AK2(1)) * AINDEX
C          ELSEIF (AINDEX.GE.0.334D0.AND.AINDEX.LT.0.667D0)THEN
C            AK = AK2(2) + (AK2(3) - AK2(2)) * AINDEX  
C          ELSEIF (AINDEX.GE.0.667D0.AND.AINDEX.LE.1.00D0)THEN
C            AK = AK2(3) + (AK2(4) - AK2(3)) * AINDEX    
C          ENDIF
C      CASE (3) ! MEDIUM DENSE
C         IF     (FREE.GE.30D0.AND.FREE.LT.32D0)THEN
C         AINDEX1 = FREE-30D0  
C         AK      = AK3(1) + (AK3(2) - AK3(1)) * AINDEX1
C         ELSEIF (FREE.GE.32D0.AND.FREE.LT.34D0)THEN
C         AINDEX2 = FREE-32D0
C        AK      = AK3(2) + (AK3(3) - AK3(2)) * AINDEX2  
C         ELSEIF (FREE.GE.34D0.AND.FREE.LE.36D0)THEN
C         AINDEX3 = FREE-34D0
C         AK      = AK3(3) + (AK3(4) - AK3(3)) * AINDEX3
C         ENDIF
C      CASE (4) ! DENSE
C         IF     (FREE.GE.36D0.AND.FREE.LT.37.334D0)THEN
C         AINDEX1 = FREE-36D0  
C         AK      = AK4(1) + (AK4(2) - AK4(1)) * AINDEX1
C         ELSEIF (FREE.GE.37.334D0.AND.FREE.LT.38.667D0)THEN
C         AINDEX2 = FREE-37.334D0
C         AK      = AK4(2) + (AK4(3) - AK4(2)) * AINDEX2  
C         ELSEIF (FREE.GE.38.667D0.AND.FREE.LE.40D0)THEN
C         AINDEX3 = FREE-38.667D0
C         AK      = AK4(3) + (AK4(4) - AK4(3)) * AINDEX3
C         ENDIF
          
C      ENDSELECT
      P   = (A*PU*TANH(AK*H*abs(Y)/(A*PU)))
      ENDIF

      AROUND = 2D0*3.141592654D0*(DPIPE/2D0)
C     SPRING STIFFNESS
      AROUND = 2D0*3.141592654D0*(DPIPE/2D0)
      STIFF(I) = ABS((P*AROUND*YFAC)/(Y*XFAC))*0.9d0
100   IF (Y.EQ.0) STIFF(I) = 0.0D0
      
      
      ENDDO
      
				   
      ! CALCULATE SAND GRAPH AND SENT OUT
      DO J = 2,16
      PLOTP(1,KPCOUNT)   = 0.0d0 
      PLOTY(1,KPCOUNT)   = 0.0d0 
      PLOTY(J,KPCOUNT)   = PLOTY(J-1,KPCOUNT) + abs(PLOTYA/30D0)!abs(PLOTY(J-1,KPCOUNT)+PLOTYA/15D0)
      PLOTP(J,KPCOUNT)   = (A*PU*TANH(AK*H*PLOTY(J,KPCOUNT)/(A*PU)))
								
      ENDDO
      
      DATAPYSAND(1:15)  = PLOTP(1:15,KPCOUNT)
      DATAPYSAND(16:30) = PLOTY(1:15,KPCOUNT)
      CALL SOILDATA_CURVE (ISN,ISP,DATAPYSAND,SBC,"PY_SAND","WRIT") 
      KPCOUNT  = KPCOUNT + 1
      PLOTYA   = 0.0d0
      END
C     ----------------------------------------------------------------   
      SUBROUTINE T_Z_CURVE_CLAY (ALENGTH,C,PO,D,YT,NCODE,NT,SBC,STIFF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION ZD(7),TTM(7)
      DIMENSION SBC(1)
      COMMON /MGRAV/ NGRAV  
      DATA ZD /0.0D0,0.30D0,0.50D0,0.75D0,0.90D0,1.0D0,0.90D0/
      DATA TTM/0.00D0,0.0016D0,0.0031D0,0.0057D0,0.0080D0,0.0010D0,0.0020D0/
      
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
      COMMON/COUNT/ICOUNT,NCOUNT
      
      
      DIMENSION STIFF(3),YT(3)
      
      IF (ICOUNT.EQ.1) YT(1:3) = 0.00001D0
      DO I = 1,3
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.1) Y = YT(1)
         IF (I.EQ.2) Y = 0.0D0
         IF (I.EQ.3) Y = 0.0D0
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) Y = 0.0D0
         IF (I.EQ.2) Y = YT(2)
         IF (I.EQ.3) Y = 0.0D0
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) Y = 0.0D0
         IF (I.EQ.2) Y = 0.0D0
         IF (I.EQ.3) Y = YT(3)
      ENDIF
      
C     --- INOPU
C     C  = UNDRAINED SHEAR STRENGTH
C     PO = EFFECTIVE OVERBURDEN PRESSURE
C     D  = PILE DIMETER
      
C     A  = DIMENSIONLESS FACTOR
      
      IF (NCODE.EQ.1)THEN
      DELTA = C/PO
      IF (DELTA.LE.1.0) A = 0.5D0*DELTA**(-0.5D0)
      IF (DELTA.GT.1.0) A = 0.5D0*DELTA**(-0.25D0)
      
      TMAX = A*C*AREA
      
C     INTERPOLATION
      STTM  = ABS(Y/D)

      IF     (STTM.GE.TTM(1).AND.STTM.LT.TTM(2))THEN
             AITTM = TTM(2) - TTM(1)
             AIZD  = ZD(2)  - ZD(1)
             UUZD  = ZD(1)  + AIZD*(STTM -  TTM(1))/AITTM
      ELSEIF  (STTM.GE.TTM(2).AND.STTM.LT.TTM(3))THEN
             AITTM = TTM(3) - TTM(2)
             AIZD  = ZD(3)  - ZD(2)
             UUZD  = ZD(2)  + AIZD*(STTM -  TTM(2))/AITTM
      ELSEIF  (STTM.GE.TTM(3).AND.STTM.LT.TTM(4))THEN  
             AITTM = TTM(4) - TTM(3)
             AIZD  = ZD(4)  - ZD(3)
             UUZD  = ZD(3)  + AIZD*(STTM -  TTM(3))/AITTM
      ELSEIF  (STTM.GE.TTM(4).AND.STTM.LT.TTM(5))THEN  
             AITTM = TTM(5) - TTM(4)
             AIZD  = ZD(5)  - ZD(4)
             UUZD  = ZD(4)  + AIZD*(STTM -  TTM(4))/AITTM
      ELSEIF  (STTM.GE.TTM(5).AND.STTM.LE.TTM(6))THEN
             AITTM = TTM(6) - TTM(5)
             AIZD  = ZD(6)  - ZD(5)
             UUZD  = ZD(5)  + AIZD*(STTM -  TTM(5))/AITTM
      ELSEIF  (STTM.GE.TTM(6).AND.STTM.LE.TTM(7))THEN
             AITTM = TTM(7) - TTM(6)
             AIZD  = ZD(7)  - ZD(6)
             UUZD  = ZD(6)  + AIZD*(STTM -  TTM(6))/AITTM
      ELSEIF  (STTM.GT.TTM(7))THEN
             UUZD  = 0.70D0
      ENDIF
      
      ENDIF
      AROUND = 2D0*3.141592654D0*(DPIPE/2D0)
C     SPRING STIFFNESS
      STIFF(I) = ABS(UUZD*TMAX*AROUND*ALENGTH/Y)
      IF (Y.EQ.0) STIFF(I) = 0.0D0
      ENDDO

      
      END
C     ----------------------------------------------------------------
      SUBROUTINE T_Z_CURVE_SAND (ALENGTH,PO,D,AK,DELTA,YT,NCODE,NT,SBC,STIFF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION ZD(2),TTM(2)
      DIMENSION SBC(1)
      DATA ZD /0.000D0,1.000D0/
      DATA TTM/0.000D0,0.100D0/
      DIMENSION STIFF(3),YT(3)
      COMMON /MGRAV/ NGRAV  
      
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
      COMMON/COUNT/ICOUNT,NCOUNT
      
      IF (ICOUNT.EQ.1) YT(1:3) = 0.0001D0
       
      DO I = 1,3
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.1) Y = YT(1)
         IF (I.EQ.2) Y = 0.0D0
         IF (I.EQ.3) Y = 0.0D0
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) Y = 0.0D0
         IF (I.EQ.2) Y = YT(2)
         IF (I.EQ.3) Y = 0.0D0
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) Y = 0.0D0
         IF (I.EQ.2) Y = 0.0D0
         IF (I.EQ.3) Y = YT(3)
      ENDIF
      
C     --- INPUT ---
C     PO    = EFFECTIVE OVERBURDEN PRESSURE
C     AK    = DIMENSION LESS COEFFICIENT OF LATERAL EARTH PRESSURE
C           = 1.0 PLUGGE
C           = 0.8 UNPLUGGE
C     DELTA = FRICTION ANGLE
      
      IF (NCODE.EQ.1) THEN
          
      TMAX  = AK*PO*TANH(DELTA) !*AREA
      
C      TMAX  = AK*PO*1680D0*TAN(0.5235987756) !*AREA
      
!      IF (TMAX.GE.0.0138889) TMAX = 0.0138889D0
      
      IF (ICOUNT.EQ.2)THEN
          A = 1
      ENDIF
C     INTERPOLATION
      STTM  = ABS(Y)
      IF     (STTM.GE.TTM(1).AND.STTM.LT.TTM(2))THEN
             AITTM = TTM(2) - TTM(1)
             AIZD  = ZD(2)  - ZD(1)
             UUZD  = ZD(1)  + AIZD*(STTM -  TTM(1))/AITTM
      ELSEIF  (STTM.GE.TTM(2))THEN
             UUZD  = 1.0D0
      ENDIF
           
      ENDIF
C     SPRING STIFFNESS
      AROUND = 2D0*3.141592654D0*(DPIPE/2D0)
      STIFF(I) = ABS(UUZD*TMAX*AROUND*ALENGTH/Y)
      IF (Y.EQ.0) STIFF(I) = 0.0D0
      ENDDO
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE Q_Z_CURVE (ANQ,PO,C,D,YT,NCODE,NT,SBC,STIFF,NOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION ZD(5),QQP(5)
      COMMON /MGRAV/ NGRAV  
      DATA ZD /0.250D0,0.500D0,0.750D0,0.900D0,1.000D0/
      DATA QQP/0.002D0,0.013D0,0.042D0,0.073D0,0.100D0/
      DIMENSION STIFF(3),YT(3)
      
      COMMON/COUNT/ICOUNT,NCOUNT
      
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
      
C     --- INPUT ---
C     ANQ  = END BEARING FACTOR
C     C    = UNDRAINED SHEAR STRENGTH
C     PO   = EFFECTIVE OVERBURDEN

      IF (ICOUNT.EQ.1) YT(1:3) = 0.001D0
      
      DO I = 1,3
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.1) Y = YT(1)
         IF (I.EQ.2) Y = 0.0D0
         IF (I.EQ.3) Y = 0.0D0
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) Y = 0.0D0
         IF (I.EQ.2) Y = YT(2)
         IF (I.EQ.3) Y = 0.0D0
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) Y = 0.0D0
         IF (I.EQ.2) Y = 0.0D0
         IF (I.EQ.3) Y = YT(3)
      ENDIF
      
      
      IF (NCODE.EQ.1)THEN
      PI = 3.141592654D0
      
      IF (NOPT.EQ.1)    Q = 9D0*C
      IF (NOPT.EQ.2)    Q = PO*ANQ

      
      QP  = Q*(PI*(D**2)/4D0)
      
C     INTERPOLATION
      QQPI = ABS(Y/D)
      IF (ICOUNT.EQ.1) QQPI = 0.002D0
      

      IF  (QQPI.GE.QQP(1).AND.QQPI.LT.QQP(2))THEN
             AIQQP = QQP(2) - QQP(1)
             AIZD  = ZD(2)  - ZD(1)
             UUZD  = ZD(1)  + AIZD*(QQPI -  QQP(1))/AIQQP
      ELSEIF (QQPI.GE.QQP(2).AND.QQPI.LT.QQP(3))THEN
             AIQQP = QQP(3) - QQP(2)
             AIZD  = ZD(3)  - ZD(2)
             UUZD  = ZD(2)  + AIZD*(QQPI -  QQP(2))/AIQQP
      ELSEIF (QQPI.GE.QQP(3).AND.QQPI.LT.QQP(4))THEN
             AIQQP = QQP(4) - QQP(3)
             AIZD  = ZD(4)  - ZD(3)
             UUZD  = ZD(3)  + AIZD*(QQPI -  QQP(3))/AIQQP
      ELSEIF (QQPI.GE.QQP(4).AND.QQPI.LE.QQP(5))THEN
             AIQQP = QQP(5) - QQP(4)
             AIZD  = ZD(5)  - ZD(4)
             UUZD  = ZD(4)  + AIZD*(QQPI -  QQP(4))/AIQQP    
      ELSEIF (QQPI.LE.QQP(1))THEN
             AIQQP = QQP(1) - 0.0D0
             AIZD  = ZD(1)  - 0.0D0
             UUZD  = 0.0D0  + AIZD*(QQPI -  0.0D0)/AIQQP
      ELSE
      UUZD = 1.0
      ENDIF
      
      ENDIF
C     SPRING STIFFNESS
      IF (ICOUNT.EQ.2)THEN
      A= 1
      ENDIF
      STIFF(I) = ABS(UUZD*QP/Y)
      IF (Y.EQ.0) STIFF(I) = 0.0D0
      
10    IF (I.EQ.1) STIFF(1) = 0.0D0
      IF (I.EQ.3) STIFF(3) = 0.0D0
      ENDDO
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE NODETOELEMENT (NODE,KELEMENT,NISN,DEPTH,NODEOUT,OPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      CHARACTER*5 OPT
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NA,NEF,NELTOT,NMV,MTYP,ISECT
CC      COMMON /STOREN/ IELEMENT(100000),N1(100000),N2(100000)
      
      COMMON /MGRAV/ NGRAV 
      COMMON/COUNT/ICOUNT,NCOUNT
      COMMON/STNO/ NODEOLD
      
      ALLOCATABLE NELEMENT(:) !SONGSAK OCT2019
CC      DIMENSION NELEMENT(NELE)
      
      DIMENSION NNODE(2),NSNODE(NSN)

C	CALL TOTAL GID ELEMENT NUMBER  NGIDM  ... SONGSAK2019
	CALL LOCATN('OGDM',KGDM,NGDAT,NGIDM,1) 
      ALLOCATE(NELEMENT(NGIDM))
      
      
      IF (OPT.EQ."ABOVE")THEN
      KELEMENT = 0
      NNODE    = 0
      NSNODE   = 0
      NELEMENT = 0
      INDEX    = 1
      
CC      DO I = 1, NELE
CC       IF (N1(I).EQ.NODE)THEN
CC          NELEMENT(INDEX) = IELEMENT(I)
CC          INDEX = INDEX + 1
CC       ENDIF
CC       IF (N2(I).EQ.NODE)THEN
CC          NELEMENT(INDEX) = IELEMENT(I)
CC          INDEX = INDEX + 1
CC       ENDIF
CC      ENDDO
      
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT LOOP ... OCT2019  

      DO 15 IGM = 1,NGIDM    
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
              CALL CALLNUMNODE_F (IGM,N1,N2) !CALL FRAME CONNECTIVITY NODE
              
               IF (N1.EQ.NODE)THEN
                  NELEMENT(INDEX) = IGM
                  INDEX = INDEX + 1
               ENDIF
               IF (N2.EQ.NODE)THEN
                  NELEMENT(INDEX) = IGM
                  INDEX = INDEX + 1
               ENDIF
          ENDIF
15    ENDDO
C     ----------------------------------------------------------------  
      
      
      DO  J = 1,INDEX-1
      CALL CALLNUMNODE_F (NELEMENT(J),NNODE(1),NNODE(2))
       IF (NODE.EQ.NNODE(1))THEN
          NSNODE(J) = NNODE(2)
       ELSEIF (NODE.EQ.NNODE(2))THEN
          NSNODE(J) = NNODE(1)
       ENDIF
      ENDDO
      
      ! NODE ABOVE
      CALL NODEABOVE (NSNODE,NODEOUT)
      
CC      DO 10 K = 1,NELE
CC          IF (NODEOUT.EQ.N1(K))THEN
CC              IF (NODE.EQ.N2(K))THEN
CC              KELEMENT = IELEMENT(K)
CC              EXIT
CC              ENDIF
CC          ENDIF          
CC           IF (NODEOUT.EQ.N2(K))THEN
CC              IF (NODE.EQ.N1(K))THEN
CC              KELEMENT =  IELEMENT(K) 
CC              EXIT
CC              ENDIF
CC          ENDIF
CC10    CONTINUE

C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT LOOP ... OCT2019  
 
      DO 20 IGM = 1,NGIDM    
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
              CALL CALLNUMNODE_F (IGM,N1,N2) !CALL FRAME CONNECTIVITY NODE
              
              IF (NODEOUT.EQ.N1)THEN
                  IF (NODE.EQ.N2)THEN
                  KELEMENT = IGM
                  EXIT
                  ENDIF
              ENDIF
          
               IF (NODEOUT.EQ.N2)THEN
                  IF (NODE.EQ.N1)THEN
                  KELEMENT =  IGM
                  EXIT
                  ENDIF
              ENDIF
          ENDIF
20    ENDDO
C     ----------------------------------------------------------------  
      
      ELSEIF (OPT.EQ."BOTTO")THEN
      KELEMENT = 0
      NNODE    = 0
      NSNODE   = 0
      NELEMENT = 0
      INDEX    = 1
      
CC      DO I = 1, NELE
CC       IF (N1(I).EQ.NODE)THEN
CC          NELEMENT(INDEX) = IELEMENT(I)
CC          INDEX = INDEX + 1
CC       ENDIF
CC       IF (N2(I).EQ.NODE)THEN
CC          NELEMENT(INDEX) = IELEMENT(I)
CC          INDEX = INDEX + 1
CC       ENDIF
CC      ENDDO

C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT LOOP ... OCT2019  
  
      DO 30 IGM = 1,NGIDM    
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
              CALL CALLNUMNODE_F (IGM,N1,N2) !CALL FRAME CONNECTIVITY NODE
              
               IF (N1.EQ.NODE)THEN
                  NELEMENT(INDEX) = IGM
                  INDEX = INDEX + 1
               ENDIF
               IF (N2.EQ.NODE)THEN
                  NELEMENT(INDEX) = IGM
                  INDEX = INDEX + 1
               ENDIF
          ENDIF
30    ENDDO
C     ----------------------------------------------------------------       
      DO 100 J = 1,INDEX-1
      CALL CALLNUMNODE_F (NELEMENT(J),NNODE(1),NNODE(2))
       IF (NODE.EQ.NNODE(1))THEN
          NSNODE(1) = NNODE(2)
          NSNODE(2) = NNODE(2)
          EXIT
       ELSEIF (NODE.EQ.NNODE(2))THEN
          NSNODE(1) = NNODE(1)
          NSNODE(2) = NNODE(2)
          EXIT
       ENDIF
100   CONTINUE
      
      ! NODE BOTTOM
      CALL NODEBOTTOM (NSNODE,NODEOUT)
      
      ELSEIF (OPT.EQ."LENGT")THEN
          
C          A1 = ABS(MOD(NISN,2))
C          IF (NISN.EQ.1) THEN
C          NODEOLD = NODE
C          RETURN
C          ENDIF
          
CC          AX1 = AX(NODE)
CC          AY1 = AY(NODE)
CC          AZ1 = AZ(NODE)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX1,1,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY1,2,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ1,3,NODE,0)  !GETTING NODAL COORDINATE
          
CC          AX2 = AX(NISN)
CC          AY2 = AY(NISN)
CC          AZ2 = AZ(NISN)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX2,1,NISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY2,2,NISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ2,3,NISN,0)  !GETTING NODAL COORDINATE
          
          IF (AX1.GT.AX2) AXT = AX1 - AX2
          IF (AX1.LE.AX2) AXT = AX2 - AX1
          IF (AY1.GT.AY2) AYT = AY1 - AY2
          IF (AY1.LE.AY2) AYT = AY2 - AY1
          IF (AZ1.GT.AZ2) AZT = AZ1 - AZ2
          IF (AZ1.LE.AZ2) AZT = AZ2 - AZ1
       
          IF (NGRAV.EQ.1)THEN
          DEPTH = AXT  
          ELSEIF (NGRAV.EQ.2)THEN
          DEPTH = AYT
          ELSEIF (NGRAV.EQ.3)THEN
          DEPTH = AZT  
          ENDIF
      ELSEIF (OPT.EQ."LENGG") THEN

CC          AX1 = AX(NISN) ! OLD
CC          AY1 = AY(NISN)
CC          AZ1 = AZ(NISN)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX1,1,NISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY1,2,NISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ1,3,NISN,0)  !GETTING NODAL COORDINATE
          
CC          AX2 = AX(NODE) ! NEW
CC          AY2 = AY(NODE)
CC          AZ2 = AZ(NODE)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX2,1,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY2,2,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ2,3,NODE,0)  !GETTING NODAL COORDINATE
          
          IF (NGRAV.EQ.1)THEN
          IF (AX1.GT.AX2)THEN
          DEPTH = 2.0D0 
          ELSEIF (AX1.LE.AX2)THEN
          DEPTH = 1.0D0     
          ENDIF
          ENDIF
          
          IF (NGRAV.EQ.2)THEN
          IF (AY1.GT.AY2)THEN
          DEPTH = 2.0D0  
          ELSEIF (AY1.LE.AY2)THEN
          DEPTH = 1.0D0     
          ENDIF
          ENDIF
          
          IF (NGRAV.EQ.3)THEN
          IF (AZ1.GT.AZ2)THEN
          DEPTH = 2.0D0  
          ELSEIF (AY1.LE.AXZ)THEN
          DEPTH = 1.0D0     
          ENDIF
          ENDIF
          
      ENDIF
      
      DEALLOCATE(NELEMENT) !SONGSAK OCT2019
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE NODEABOVE (NSNODE,NODEOUT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /MGRAV/ NGRAV
      DIMENSION NSNODE(NSN)
      DIMENSION ANNX(NSN),ANNY(NSN),ANNZ(NSN)
      
      ANNX   = 0.
      ANNY   = 0.
      ANNZ   = 0.
      NPOTIN = 0.
      COORMAX = 0.0D0
      NCCCO   = 0
      DO IJ = 1,NSN
          IF (NSNODE(IJ).NE.0)THEN
              NCCCO = NCCCO + 1
          ENDIF
      ENDDO
      
      DO 100 I =1,NSN
      IF (NSNODE(I).EQ.0) EXIT
      
CC      DO 10 J=1,NSN
CC         IF (NND(J).EQ.NSNODE(I))THEN
CC         ANNX(I) = AX(J)
CC         ANNY(I) = AY(J)
CC         ANNZ(I) = AZ(J)
CC         EXIT
CC         ENDIF
CC10    CONTINUE
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',ANNX(I),1,NSNODE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',ANNY(I),2,NSNODE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',ANNZ(I),3,NSNODE(I),0)  !GETTING NODAL COORDINATE
      
100   CONTINUE
      
      IF (NGRAV.EQ.1)THEN
      NPOINTMAX = MAXLOC(ANNX(1:NCCCO),DIM=1)
      NPOINTLOW = MINLOC(ANNX(1:NCCCO),DIM=1)
         IF (ANNX(NPOINTMAX).GE.ANNX(NPOINTLOW))THEN ! 2018 ADD ABS() FOR - COOR
             NPOINT = NPOINTMAX
         ELSEIF (ANNX(NPOINTMAX).LT.ANNX(NPOINTLOW))THEN
             NPOINT = NPOINTLOW
         ENDIF
      ELSEIF (NGRAV.EQ.2)THEN
      NPOINTMAX = MAXLOC(ANNY(1:NCCCO),DIM=1)
      NPOINTLOW = MINLOC(ANNY(1:NCCCO),DIM=1)
         IF ((ANNY(NPOINTMAX)).GE.(ANNY(NPOINTLOW)))THEN ! 2018 ADD ABS() FOR - COOR
             NPOINT = NPOINTMAX
         ELSEIF ((ANNY(NPOINTMAX)).LT.(ANNY(NPOINTLOW)))THEN
             NPOINT = NPOINTLOW
         ENDIF
      ELSEIF (NGRAV.EQ.3)THEN
      NPOINTMAX = MAXLOC(ANNZ(1:NCCCO),DIM=1)
      NPOINTLOW = MINLOC(ANNZ(1:NCCCO),DIM=1)
         IF (ANNZ(NPOINTMAX).GE.ANNZ(NPOINTLOW))THEN ! 2018 ADD ABS() FOR - COOR
             NPOINT = NPOINTMAX
         ELSEIF (ANNZ(NPOINTMAX).LT.ANNZ(NPOINTLOW))THEN
             NPOINT = NPOINTLOW
         ENDIF
      ENDIF
      
      
      NODEOUT = NSNODE(NPOINT)
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE NODEABOVE2 (NSNODE,NODEOUT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /MGRAV/ NGRAV
      DIMENSION NSNODE(NSN)
      DIMENSION ANNX(NSN),ANNY(NSN),ANNZ(NSN)
      
      ANNX   = 0.
      ANNY   = 0.
      ANNZ   = 0.
      NPOTIN = 0.
      COORMAX = 0.0D0
      NCCCO   = 0
      DO IJ = 1,2
          IF (NSNODE(IJ).NE.0)THEN
              NCCCO = NCCCO + 1
          ENDIF
      ENDDO
      
      DO 100 I =1,2
      IF (NSNODE(I).EQ.0) EXIT
      
CC      DO 10 J=1,NSN
CC         IF (NND(J).EQ.NSNODE(I))THEN
CC         ANNX(I) = AX(J)
CC         ANNY(I) = AY(J)
CC         ANNZ(I) = AZ(J)
CC         EXIT
CC         ENDIF
CC10    CONTINUE
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',ANNX(I),1,NSNODE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',ANNY(I),2,NSNODE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',ANNZ(I),3,NSNODE(I),0)  !GETTING NODAL COORDINATE
      
100   CONTINUE
      
      IF (NGRAV.EQ.1)THEN
      NPOINTMAX = MAXLOC(ANNX(1:NCCCO),DIM=1)
      NPOINTLOW = MINLOC(ANNX(1:NCCCO),DIM=1)
         IF (ANNX(NPOINTMAX).GE.ANNX(NPOINTLOW))THEN ! 2018 ADD ABS() FOR - COOR
             NPOINT = NPOINTMAX
         ELSEIF (ANNX(NPOINTMAX).LT.ANNX(NPOINTLOW))THEN
             NPOINT = NPOINTLOW
         ENDIF
      ELSEIF (NGRAV.EQ.2)THEN
      NPOINTMAX = MAXLOC(ANNY(1:NCCCO),DIM=1)
      NPOINTLOW = MINLOC(ANNY(1:NCCCO),DIM=1)
         IF ((ANNY(NPOINTMAX)).GE.(ANNY(NPOINTLOW)))THEN ! 2018 ADD ABS() FOR - COOR
             NPOINT = NPOINTMAX
         ELSEIF ((ANNY(NPOINTMAX)).LT.(ANNY(NPOINTLOW)))THEN
             NPOINT = NPOINTLOW
         ENDIF
      ELSEIF (NGRAV.EQ.3)THEN
      NPOINTMAX = MAXLOC(ANNZ(1:NCCCO),DIM=1)
      NPOINTLOW = MINLOC(ANNZ(1:NCCCO),DIM=1)
         IF (ANNZ(NPOINTMAX).GE.ANNZ(NPOINTLOW))THEN ! 2018 ADD ABS() FOR - COOR
             NPOINT = NPOINTMAX
         ELSEIF (ANNZ(NPOINTMAX).LT.ANNZ(NPOINTLOW))THEN
             NPOINT = NPOINTLOW
         ENDIF
      ENDIF
      
      
      NODEOUT = NSNODE(NPOINT)
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE NODEBOTTOM (NSNODE,NODEOUT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /MGRAV/ NGRAV
      DIMENSION NSNODE(NSN)
      DIMENSION ANNX(NSN),ANNY(NSN),ANNZ(NSN)
      
      ANNX   = 123456789.
      ANNY   = 123456789.
      ANNZ   = 123456798.
      NPOTIN = 0.
      
      DO 100 I =1,NSN
      IF (NSNODE(I).EQ.0) EXIT
      
CC      DO 10 J=1,NN
CC         IF (NND(J).EQ.NSNODE(I))THEN
CC         ANNX(I) = AX(J)
CC         ANNY(I) = AY(J)
CC         ANNZ(I) = AZ(J)
CC         EXIT
CC         ENDIF
CC10    CONTINUE
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',ANNX(I),1,NSNODE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',ANNY(I),2,NSNODE(I),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',ANNZ(I),3,NSNODE(I),0)  !GETTING NODAL COORDINATE
      
100   CONTINUE
      
      IF (NGRAV.EQ.1)THEN
      NPOINT = MINLOC(ANNX,DIM=1) 
      ELSEIF (NGRAV.EQ.2)THEN
      NPOINT = MINLOC(ANNY,DIM=1)
      ELSEIF (NGRAV.EQ.3)THEN
      NPOINT = MINLOC(ANNZ,DIM=1)   
      ENDIF
      
      NODEOUT = NSNODE(NPOINT)
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE READSOILDATA (NOPT,NODEOUT,DEFAC,NCODEO,NT,SEBC,NODEIN,ISN,PROPOUT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 NOPT
      DIMENSION SEBC(37),PROPOUT(3)
      !DIMENSION NPROP(10000)
      
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOILOPT/ NSOIL,NNSO
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
      COMMON /NSOILT/ NNODEA(10000),NPROP(10000)
       
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NA,NEF,NELTOT,NMV,MTYP,ISECT
      
      COMMON /CALENODE/ NODECALE1(1000),NODECALE2(1000)
      
      
C      COMMON /SOILDATA/ NODEA(2000),DFAC(2000),SBC(6,2000)
      
      COMMON /SOILDATA/ NODEA(1000),SBC(37,1000),ITER(1000),NCODE(1000),PROP(2,1000),START_ELV(1000),END_ELV(1000)
     1                 ,NNODE_DIRECT(1000,100),NTOTAL_NODE(100)
      
      DIMENSION NNODEAR(500)
      
      IF (ISOLOP.NE.5) THEN
          
      CALL SOILMODULE (ISOLOP,"READD")
      
C      IF (ISOLOP.EQ.1) RETURN
      
      IF (NOPT.EQ."READ")THEN
      ! BY FOLOOWING TCL/TK SBC(1:25)
      ! SBC(1)  = ANALYSIS P-Y CUVRE CONTROL (0,1)
      ! SBC(2)  = ANALYSIS T-Z CUVRE CONTROL (0,1)
      ! SBC(3)  = ANALYSIS Q-Z CUVRE CONTROL (0,1) RERATE WITH T-Z CURVE
      ! SBC(4)  = TYPE OF STIFFNESS (API,USER DEFINED)
      READ(ITI,*) 
      
      READ(ITI,*) NSOIL !,NNSO
      READ(ITI,*)
      DO I = 1,NSOIL
         READ(ITI,*) NCODE(I)
         IF(NCODE(I).EQ.1.OR.NCODE(I).EQ.4)THEN
         READ(ITI,*) NODEA(I),START_ELV(I),END_ELV(I),SBC(1:25,NODEA(I))
         ELSEIF (NCODE(I).EQ.2)THEN
         READ(ITI,*) NODEA(I),NTOTAL_NODE(I),NNODE_DIRECT(1:NTOTAL_NODE(I),I),SBC(1:3,NODEA(I))    
         ELSEIF (NCODE(I).EQ.3)THEN
         READ(ITI,*) NODEA(I),START_ELV(I),END_ELV(I),ITER(I),SBC(1:3,NODEA(I)),SBC(36,NODEA(I)),SBC(37,NODEA(I))
         READ(ITI,*) SBC(4:19,NODEA(I))   
         READ(ITI,*) SBC(20:35,NODEA(I))
         ENDIF
      ENDDO
      
      ! ORIGINAL 2 NODES DATA
      !READ(ITI,*) NSOIL ,NNSO
      !READ(ITI,*)
      ! DO I = 1,NSOIL
      !   READ(ITI,*) NCODE(I)
      !   IF(NCODE(I).EQ.1.OR.NCODE(I).EQ.4)THEN
      !   READ(ITI,*) NODEA(I),SBC(1:25,NODEA(I))
      !   ELSEIF (NCODE(I).EQ.2)THEN
      !   READ(ITI,*) NODEA(I),SBC(1:3,NODEA(I))    
      !   ELSEIF (NCODE(I).EQ.3)THEN
      !   READ(ITI,*) NODEA(I),ITER(I),SBC(1:3,NODEA(I)),SBC(36,NODEA(I)),SBC(37,NODEA(I))
      !   READ(ITI,*) SBC(4:19,NODEA(I))   
      !   READ(ITI,*) SBC(20:35,NODEA(I))
      !   ENDIF
      !ENDDO
       
      !READ(ITI,*)
      
      !! ---- SKIP 09/2021 ----
      !DO I =1,NNSO
      !READ(ITI,*)  NNODEA(I),NPROP(I)
      !ENDDO
      !
      !! ---- SKIP 09/2021 ----
      !TCOUNT = 1
      !DO I = 2,NNSO,2
      !    NODECALE1(TCOUNT) = NNODEA(I)
      !    NODECALE2(TCOUNT) = NNODEA(I-1)
      !    TCOUNT       = TCOUNT + 1
      !ENDDO
      
      
      ELSEIF (NOPT.EQ."CALL")THEN
      SEBC = 0.D0
      SEBC(1:37) = SBC(1:37,NPROP(NODEIN))
      NT         = ITER(NPROP(NODEIN))
      NCODEO     = NCODE(NPROP(NODEIN))
      ISN        = NNODEA(NODEIN)
      
      ELSEIF (NOPT.EQ."CALO")THEN ! FOR OUTPUT
      SEBC = 0.D0
      DO 10 I = 1,NNSO
         IF (NNODEA(I).EQ.NODEIN) THEN
         INDEX_NODE = I
         EXIT
         ENDIF
10    ENDDO
      
      SEBC(1:37) = SBC(1:37,NPROP(INDEX_NODE))
      NT         = ITER(NPROP(INDEX_NODE))
      NCODEO     = NCODE(NPROP(INDEX_NODE))
      ISN        = NNODEA(INDEX_NODE)
C          EXIT
C          ENDIF
C10    CONTINUE
      ELSEIF (NOPT.EQ."CALE")THEN
      DO 100 I =1,NNSO
          IF (NODEOUT.EQ.NODECALE1(I))THEN
              IF (NODEIN.EQ.NODECALE2(I))THEN
                  NNNC = I*2
                  EXIT
              ENDIF
          ELSEIF (NODEOUT.EQ.NODECALE2(I))THEN
              IF (NODEIN.EQ.NODECALE1(I))THEN
                  NNNC = I*2
                  EXIT
              ENDIF
          ENDIF
          
100   CONTINUE
      SEBC       = 0.D0
      SEBC(1:37) = SBC(1:37,NPROP(NNNC))
      NT         = ITER(NPROP(NNNC))
      NCODEO     = NCODE(NPROP(NNNC))
      ISN        = NNODEA(NNNC)
      
      ELSEIF (NOPT.EQ."CALN")THEN
      NCODEO     = NCODE(NPROP(NODEIN))
      
      ELSEIF (NOPT.EQ."CALT")THEN
      NODEOUT = NNSO
      
      ELSEIF (NOPT.EQ."CALR")THEN
      
          DO I = 1,NNSO
              INDEX = 1
              DO J =1,NSOIL
              IF (NPROP(J).EQ.I)THEN
                  NNODEAR(INDEX) = NNODEA(J)
                  INDEX = INDEX + 1
              ENDIF
              ENDDO
!              CALL REARRAYNODE (NNODEAR,INDEX)
          ENDDO
          
      ENDIF
      
      ELSEIF (ISOLOP.EQ.5) THEN ! DYNAMIC CASE
          
      IF (NOPT.EQ."READ")THEN
      !READ(ITI,*) 
      !READ(ITI,*) NSOIL,NNSO
      !READ(ITI,*)
      !DO I = 1,NSOIL
      !   READ(ITI,*) DUM,PROP(1:3,I)
      !ENDDO
      !
      !READ(ITI,*)
      !
      !DO I =1,NNSO
      !READ(ITI,*)  NNODEA(I),NPROP(I)
      !ENDDO
          
          
      ! BY FOLOOWING TCL/TK SBC(1:25)
      ! SBC(1)  = ANALYS5IS P-Y CUVRE CONTROL (0,1)
      ! SBC(2)  = ANALYSIS T-Z CUVRE CONTROL (0,1)
      ! SBC(3)  = ANALYSIS Q-Z CUVRE CONTROL (0,1) RERATE WITH T-Z CURVE
      ! SBC(4)  = TYPE OF STIFFNESS (API,USER DEFINED)
      READ(ITI,*) 
      READ(ITI,*) NSOIL,NNSO
      READ(ITI,*)
      DO I = 1,NSOIL
         READ(ITI,*) NCODE(I)
         IF(NCODE(I).EQ.1)THEN
         READ(ITI,*) NODEA(I),SBC(1:25,NODEA(I))
         ELSEIF (NCODE(I).EQ.2)THEN
         READ(ITI,*) NODEA(I),SBC(1:3,NODEA(I))    
         ELSEIF (NCODE(I).EQ.3)THEN
         READ(ITI,*) NODEA(I),ITER(I),SBC(1:3,NODEA(I)),SBC(36,NODEA(I)),SBC(37,NODEA(I))
         READ(ITI,*) SBC(4:19,NODEA(I))   
         READ(ITI,*) SBC(20:35,NODEA(I))
         ENDIF
      ENDDO
      
      READ(ITI,*)
      
      DO I =1,NNSO
      READ(ITI,*)  NNODEA(I),NPROP(I)
      ENDDO
      
      TCOUNT = 1
      DO I = 2,NNSO,2
          NODECALE1(TCOUNT) = NNODEA(I)
          NODECALE2(TCOUNT) = NNODEA(I-1)
          TCOUNT       = TCOUNT + 1
      ENDDO
      
      
      
      ELSEIF (NOPT.EQ."CALL")THEN
      PROPOUT = 0.D0
      
      PROPOUT(1:3) = PROP(1:3,NPROP(NODEIN))
      
      ELSEIF (NOPT.EQ."CALT")THEN
      NODEOUT = NNSO
      
       
      ENDIF
      
      ENDIF    
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE ADDSPISTIFF (ID,MAXA,SD,SD2,SDWOK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

C	NEW SPRING SUPPORT SONGSAK JUL2006
	COMMON /SPBC/ NSS,NLSS

	DIMENSION ID(NSF,NSN),SD(6,NSN),MAXA(1)
	DIMENSION SD2(42,NSS),SDWOK(30,NSS)
	DIMENSION ST2(7,6),SDW(5,6)

      DO I = 1,NSN
         DO J = 1,6
         YOUNG = SD(J,I)
         IEQ   = ID(J,I)  
	   CALL MESTIF(IEQ,YOUNG,1,1,'WRT')
         ENDDO
      ENDDO

	RETURN
      END
C	=====================================================================
C	=====================================================================
C	===================================================================== 
      SUBROUTINE SOILMODULE (ISOLOP,OPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*5 OPT
      COMMON /NSOIL/ NVSOIL
      IF (OPT.EQ."STORE")THEN
      NVSOIL = ISOLOP
      ISOLOP = 1
      ELSEIF (OPT.EQ."READD")THEN
      ISOLOP = NVSOIL
      IF (NVSOIL.EQ.0)  ISOLOP = 1
      ENDIF
      END
C	=====================================================================
C	=====================================================================
C	===================================================================== 
      SUBROUTINE USERFUNCTIONCURVE (ISN,ISP,YT,NT,SBC,STIFF,ALENGTH,DEPTH,NTSOIL)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION SBC(37)
      DIMENSION X(16),Y(16)
      COMMON /MGRAV/ NGRAV 
      DIMENSION YT(3),STIFF(3)
      
      COMMON/COUNT/ICOUNT,NCOUNT
      
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
      COMMON / DISSOIL / DATADISP(10000)
      STIFF = 0.
      ! DEFIND DATA ON X AXIS
      DO I = 1,16
      X(I) = SBC(I+3)    
      ENDDO
      ! DEFIND DATA ON Y AXIS
      DO I = 1,16
      Y(I) = SBC(I+19)!*1000D0   
      ENDDO
      
      ! MUTIPLY FACTOR
      X(1:16) = X(1:16) * SBC(36)
      Y(1:16) = Y(1:16) * SBC(37)
      
      CALL SOILDATA_CURVE (ISN,ISP,DATAPYCLAY,SBC,"USERDEF","WRIT") 
      
      
      IF (ICOUNT.EQ.1) YT(1:3) = (MAXVAL(X)-MINVAL(X))/10D0
      
      ! 3 STIFFNESS ( KX,KY,KZ ) 
      DO I= 1,3
      IF (SBC(1).EQ.1)THEN ! P-Y CURVE
          
      IF (NGRAV.EQ.1)THEN ! X-DIRECTION
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = ABS(YT(3))
      ELSEIF (NGRAV.EQ.2)THEN ! Y-DIRECTION
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = ABS(YT(3))
      ELSEIF (NGRAV.EQ.3)THEN ! Z-DIRECTION
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = 0.0D0
      ENDIF
      
      ELSEIF (SBC(2).EQ.1)THEN ! T-Z CURVE
          
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = ABS(YT(3))
      ENDIF
      
      ELSEIF (SBC(3).EQ.1)THEN ! Q-Z CURVE
          
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = ABS(YT(3))
      ENDIF
      
      ENDIF
      
      ! INTERPORATION OF DATA
      IF     (AOU.GE.X(1).AND.AOU.LT.X(2))THEN
             AIY = Y(2) - Y(1)
             AIX = X(2) - X(1)
             PP  = Y(1) + ABS((AIY/AIX)*(X(1)-AOU))
      ELSEIF (AOU.GE.X(2).AND.AOU.LT.X(3))THEN
             AIY = Y(3) - Y(2)
             AIX = X(3) - X(2)
             PP  = Y(2) + ABS((AIY/AIX)*(X(2)-AOU))
      ELSEIF (AOU.GE.X(3).AND.AOU.LT.X(4))THEN
             AIY = Y(4) - Y(3)
             AIX = X(4) - X(3)
             PP  = Y(3) + ABS((AIY/AIX)*(X(3)-AOU))
      ELSEIF (AOU.GE.X(4).AND.AOU.LT.X(5))THEN
             AIY = Y(5) - Y(4)
             AIX = X(5) - X(4)
             PP  = Y(4) + ABS((AIY/AIX)*(X(4)-AOU))
      ELSEIF (AOU.GE.X(5).AND.AOU.LT.X(6))THEN
             AIY = Y(6) - Y(5)
             AIX = X(6) - X(5)
             PP  = Y(5) + ABS((AIY/AIX)*(X(5)-AOU))
      ELSEIF (AOU.GE.X(6).AND.AOU.LT.X(7))THEN
             AIY = Y(7) - Y(6)
             AIX = X(7) - X(6)
             PP  = Y(6) + ABS((AIY/AIX)*(X(6)-AOU))
      ELSEIF (AOU.GE.X(7).AND.AOU.LT.X(8))THEN
             AIY = Y(8) - Y(7)
             AIX = X(8) - X(7)
             PP  = Y(7) + ABS((AIY/AIX)*(X(7)-AOU))
      ELSEIF (AOU.GE.X(8).AND.AOU.LT.X(9))THEN
             AIY = Y(9) - Y(8)
             AIX = X(9) - X(8)
             PP  = Y(8) + ABS((AIY/AIX)*(X(8)-AOU))
      ELSEIF (AOU.GE.X(9).AND.AOU.LT.X(10))THEN
             AIY = Y(10) - Y(9)
             AIX = X(10) - X(9)
             PP  = Y(9) + ABS((AIY/AIX)*(X(9)-AOU))
      ELSEIF (AOU.GE.X(10).AND.AOU.LT.X(11))THEN
             AIY = Y(11) - Y(10)
             AIX = X(11) - X(10)
             PP  = Y(10) + ABS((AIY/AIX)*(X(10)-AOU))
      ELSEIF (AOU.GE.X(11).AND.AOU.LT.X(12))THEN
             AIY = Y(12) - Y(11)
             AIX = X(12) - X(11)
             PP  = Y(11) + ABS((AIY/AIX)*(X(11)-AOU))
      ELSEIF (AOU.GE.X(12).AND.AOU.LT.X(13))THEN
             AIY = Y(13) - Y(12)
             AIX = X(13) - X(12)
             PP  = Y(12) + ABS((AIX*(AOU-Y(12))/AIY))
      ELSEIF (AOU.GE.X(13).AND.AOU.LT.X(14))THEN
             AIY = Y(14) - Y(13)
             AIX = X(14) - X(13)
             PP  = Y(13) + ABS((AIY/AIX)*(X(13)-AOU))
      ELSEIF (AOU.GE.X(14).AND.AOU.LT.X(15))THEN
             AIY = Y(15) - Y(14)
             AIX = X(15) - X(14)
             PP  = Y(14) + ABS((AIY/AIX)*(X(14)-AOU))
      ELSEIF (AOU.GE.X(15).AND.AOU.LT.X(16))THEN
             AIY = Y(16) - Y(15)
             AIX = X(16) - X(15)
             PP  = Y(15) + ABS((AIY/AIX)*(X(15)-AOU))
      ELSEIF (AOU.GE.16)THEN
             PP  = Y(16)
      ELSE
             PP  =MAXVAL(Y)*1.0
      ENDIF
      
      AROUND = 2D0*3.141592654D0*(DPIPE/2D0)
      IF (SBC(1).EQ.1)THEN ! P-Y CURVE
      PP = PP*DEPTH
      ELSEIF (SBC(2).EQ.1)THEN ! T-Z CURVE
      PP = PP*AROUND*DEPTH
      ELSEIF (SBC(3).EQ.1)THEN ! Q-Z CURVE
          
      ENDIF
      ! STIFFNESS
      STIFF(I) = ABS(PP)/AOU
      IF (AOU.EQ.0) STIFF(I) = 0.0D0
      
!      IF (I.EQ.2) WRITE(7101,'(E12.5)') PP
      ENDDO
      
      END
C     ----------------------------------------------------------------
      SUBROUTINE DYNAMIC_SOIL (S,NELEMENT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION S(105),PROPOUT(3)
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
!      DO 10 I=1,NELE
!       IF (IELEMENT(I).EQ.NELEMENT) THEN
!       NNODE1 = N1(I)
!       NNODE2 = N2(I)
!       EXIT
!       ENDIF
!10    ENDDO
      
CC       NNODE1 = N1(NELEMENT)
CC       NNODE2 = N2(NELEMENT)
      CALL CALLNUMNODE_F (NELEMENT,NNODE1,NNNODE2) !SONGSAK ADD THIS FOR N1 N2 ... OCT2019
       
      CALL READSOILDATA ("CALL",NODEOUT,DEFAC,NCODEO,NT,SEBC,NNODE1,ISN,PROPOUT)
      
      IF (PROPOUT(1).EQ.0.0D0.AND.PROPOUT(2).EQ.0D0.AND.PROPOUT(3).EQ.0D0) GOTO 30
      
      AKHH = PROPOUT(1)
      AKHM = PROPOUT(2)
      AKMM = PROPOUT(3)
      
      SELECTCASE(NGRAV)
          
      CASE(1)
      ! X DIRECTION
          
      CASE(2)
      ! Y DIRECTION
      S(1)   = S(1)  + AKHH
      S(28)  = S(28) + AKHH
      S(40)  = S(40) + AKMM
      S(61)  = S(61) + AKMM
      
      S(6)   = S(6)  + AKHM
      S(29)  = S(29) + AKHM
      
      CASE(3)
      ! Z DIRECTION
      S(1)   = S(1)  + AKHH
      S(15)  = S(15) + AKHH
      S(40)  = S(40) + AKMM
      S(51)  = S(51) + AKMM
      
      S(6)   = S(6)  + AKHM
      S(17)  = S(17) + AKHM
      ENDSELECT
      GOTO 40
      
30    CALL READSOILDATA (NOPT,NODEOUT,DEFAC,NCODEO,NT,SEBC,NNODE2,ISN,PROPOUT)
      
      IF (PROPOUT(1).EQ.0.0D0.AND.PROPOUT(2).EQ.0D0.AND.PROPOUT(3).EQ.0D0) GOTO 40
      
      AKHH = PROPOUT(1)
      AKHM = PROPOUT(2)
      AKMM = PROPOUT(3)
      
      SELECTCASE(NGRAV)
          
      CASE(1)
      ! X DIRECTION
          
      CASE(2)
      ! Y DIRECTION
      S(78)   = S(78)  + AKHH
      S(91)   = S(91)  + AKHH
      S(96)   = S(96)  + AKMM
      S(103)  = S(103) + AKMM
      
      S(83)   = S(83)  + AKHM
      S(92)   = S(92)  + AKHM
      
      CASE(3)
      ! Z DIRECTION
      S(78)   = S(78)  + AKHH
      S(85)   = S(85)  + AKHH
      S(96)   = S(96)  + AKMM
      S(100)  = S(100) + AKMM
      
      S(83)   = S(83)  + AKHM
      S(87)   = S(87)  + AKHM
      ENDSELECT
      
40    RETURN 
      END
      
C     ----------------------------------------------------------------
      SUBROUTINE PULLOUT (ILC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*200 LHDGRP
      CHARACTER*200 NAMEG
      CHARACTER(100) NUMBERFILE,DIRPRIMARY


      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
      COMMON /SOILOPT/ NSOIL,NNSO
CC      COMMON /STOREN/ IELEMENT(100000),N1(100000),N2(100000)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /NSOILT/ NNODEA(10000)
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                         ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      
      COMMON / LOCALFORCE / DESIGNAXIAL(100000),DESIGNSHEART(100000),DESIGNSHEARS(100000)
     1                     ,DESIGNTORSION(100000),DESIGNMOMENTT(100000),DESIGNMOMENTS(100000)
            
      COMMON /GHEADER/ LHDGRP(500)
      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      
      COMMON / SOILSTIFFNESS / SSTIFF(6,1000)

C     SONGSAK OCT2019 / TOEY 2020      
CC      ALLOCATABLE NCONNECT1(:,:),NCONNECT2(:,:)
CC      ALLOCATABLE NELEMENTGR(:,:),NELEMENTGRNEW(:,:)
CC      ALLOCATABLE N1GROUP(:,:),N2GROUP(:,:)
CC      ALLOCATABLE N1GROUPNEW(:,:),N2GROUPNEW(:,:)
CC      ALLOCATABLE ADATA(:,:),BDATA(:),NSOILELEMENT(:,:),NPRINT(:)
CC      ALLOCATABLE NSOILELEMENT1(:,:),NSOILELEMENT2(:,:)
CC      ALLOCATABLE NSOILELEMENT1NEW(:,:),NSOILELEMENT2NEW(:,:)
CC      ALLOCATABLE ADATANEW(:,:),BDATANEW(:)
CC      ALLOCATABLE BEARINGDATA(:,:),BEARINGDATANEW(:,:)
CC      ALLOCATABLE NSOILELEMENTPILE(:,:),NSOILELEMENTPILENEW(:,:)
CC      ALLOCATABLE ADDSTIFFNESS1(:,:),ADDSTIFFNESS2(:,:),ADDSTIFFNESS3(:,:)
CC      ALLOCATABLE NODEINDEX(:,:)
      
      DIMENSION NCONNECT1(NELE,500),NCONNECT2(NELE,500)
      DIMENSION NELEMENTGR(NELE,500),NELEMENTGRNEW(NELE,500)
      DIMENSION N1GROUP(NELE,500),N2GROUP(NELE,500)
      DIMENSION N1GROUPNEW(NELE,500),N2GROUPNEW(NELE,500)      
      DIMENSION ADATA(1000,500),BDATA(1000),NSOILELEMENT(1000,2),NPRINT(1000)
      DIMENSION NSOILELEMENT1(1000,500),NSOILELEMENT2(1000,500)
      DIMENSION NSOILELEMENT1NEW(1000,500),NSOILELEMENT2NEW(1000,500)
      DIMENSION ADATANEW(1000,500),BDATANEW(1000)
      DIMENSION BEARINGDATA(1000,500),BEARINGDATANEW(1000,500)
      DIMENSION NSOILELEMENTPILE(1000,500),NSOILELEMENTPILENEW(1000,500)
      DIMENSION ADDSTIFFNESS1(1000,500),ADDSTIFFNESS2(1000,500),ADDSTIFFNESS3(1000,500)
      DIMENSION NODEINDEX(1000,500)

      DIMENSION SBC(37)
      
      DIMENSION X(16),NODETOPMAX(10000),NODEBOTTOMAX(10000),NTOP(10000)
      
      DIMENSION IINDEX(500)
      
      DIMENSION SUMTOTAL(500)
      
      DIMENSION NELEMENTINDEX(500),NELEMENTINDEXBEARING(500)
      
      DIMENSION SD(6,NSN),NSNODE(NSN)
      
      DIMENSION POTOTALGROUP(500)
      
      DIMENSION INDEX_TOTOAL_OUT(500),NNCODE_OUT(1000,500),NOUT_NODE(1000,500)
      
      COMMON /SOIL_LOCATE/NOUT_NODE_SUPER(1000,500),INDEX_OUT_G(500),NGROUP
      
C     PULLOUT SUBROTINE HAS CAPACITY TO PRINTING THE RESULT BY FOLLOWING
C     - PULLOUT AND END BEARING AT EACH PILE SECTION ON FILE NAME "Pullout and End Bearing Capacity"
C     - SOIL STIFFNESS FOR USER TO INPUT IN REAL MODEL "/SOIL-STIFFNESS/SOIL_STIFFNESS1" (RELATE WITH THE LOCASE NUMBER)      
      
      ! INTIAL SETTING
      QF         = 0D0
      NCOUNT     = 0
      NPRINT     = 0
      NELEMENTGR = 0
      N1GROUP    = 0
      N2GROUP    = 0
      NELEMENTGRNEW = 0
      N1GROUPNEW    = 0
      N2GROUPNEW    = 0
      NSOILELEMENTPILE    = 0
      NSOILELEMENTPILENEW = 0
      NSOILELEMENT1NEW    = 0
      NSOILELEMENT2NEW    = 0
      IINDEX              = 0
      INDEX_TOTOAL_OUT = 0.
      INDEX_OUT_G      = 0.
      
      AMAXBEARING         = 0.0D0
      BEARINGDATANEW      = 0.0D0
      NOUT_NODE = 0.
      
C     SONGSAK CREATE NEXT BLOCK ... OCT2019  TWEAK MEMORY & SPEED    
	CALL LOCATN('OGDM',KGDM,NGDAT,NGIDM,1)  !CALLING TOTAL NUMBER OF GID ELEMENT   ... SONGSAK OCT2019
      CALL INTFILL('%NUB',NGGEO,1,13,0)       !CALLING TOTAL NUMBER OF GEOMETRIC SET ... SONGSAK OCT2019
      CALL CALLNUMFRAME(NFRAM)                !CALLING NUMBER OF FRAME ELEMENT FROM OVERALL GROUP
      
CC      ALLOCATE(NCONNECT1(NGIDM,NGGEO),NCONNECT2(NGIDM,NGGEO))
CC      ALLOCATE(NELEMENTGR(NGIDM,NGGEO),NELEMENTGRNEW(NGIDM,NGGEO))
CC      ALLOCATE(N1GROUP(NGIDM,NGGEO),N2GROUP(NGIDM,NGGEO))
CC      ALLOCATE(N1GROUPNEW(NGIDM,NGGEO),N2GROUPNEW(NGIDM,NGGEO))
CC      ALLOCATE(ADATA(NGIDM,NGGEO),BDATA(NGIDM),NSOILELEMENT(NGIDM,2),NPRINT(NGIDM))
CC      ALLOCATE(NSOILELEMENT1(NGIDM,NGGEO),NSOILELEMENT2(NGIDM,NGGEO))
CC      ALLOCATE(NSOILELEMENT1NEW(NGIDM,NGGEO),NSOILELEMENT2NEW(NGIDM,NGGEO))
CC      ALLOCATE(ADATANEW(NGIDM,NGGEO),BDATANEW(NGIDM))
CC      ALLOCATE(BEARINGDATA(NGIDM,NGGEO),BEARINGDATANEW(NGIDM,NGGEO))
CC      ALLOCATE(NSOILELEMENTPILE(NGIDM,NGGEO),NSOILELEMENTPILENEW(NGIDM,NGGEO))
CC      ALLOCATE(ADDSTIFFNESS1(NGIDM,NGGEO),ADDSTIFFNESS2(NGIDM,NGGEO),ADDSTIFFNESS3(NGIDM,NGGEO))
CC      ALLOCATE(NODEINDEX(NGIDM,NGGEO))      
      
      ! CLASS GROUP TO CONNECTIVITY
CC      DO I = 1,NELE
CC          NCGROUP               = NGEOMETRY(I)
CC          NELEMENTGR(I,NCGROUP) = IELEMENTGEOMETRY(I)
CC          N1GROUP(I,NCGROUP)    = N1(IELEMENTGEOMETRY(I))
CC          N2GROUP(I,NCGROUP)    = N2(IELEMENTGEOMETRY(I))
CC      ENDDO
      NCONNECT1            = 0.
      NCONNECT2            = 0.
      NELEMENTGR           = 0.
      NELEMENTGRNEW        = 0.
      N1GROUP              = 0.
      N2GROUP              = 0.
      N1GROUPNEW           = 0.
      N2GROUPNEW           = 0.
      ADATA                = 0.
      BDATA                = 0.
      NSOILELEMENT         = 0.
      NPRINT               = 0.
      NSOILELEMENT1        = 0.
      NSOILELEMENT2        = 0.
      NSOILELEMENT1NEW     = 0.
      NSOILELEMENT2NEW     = 0.
      ADATANEW             = 0.
      BDATANEW             = 0.
      BEARINGDATA          = 0.
      BEARINGDATANEW       = 0.
      NSOILELEMENTPILE     = 0.
      NSOILELEMENTPILENEW  = 0.
      ADDSTIFFNESS1        = 0.
      ADDSTIFFNESS2        = 0.
      ADDSTIFFNESS3        = 0.
      NODEINDEX            = 0.
      NSNODE               = 0.
C     ---------------------------------------------- 
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT LOOP ... OCT2019    
      IFRAM = 0
      DO IGM = 1,NGIDM    
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
              IFRAM = IFRAM + 1
	        CALL INTFILL('OGDM',IGST,3,IGM,0) !CALL GEOMETRIC SET INDEX ... IGST
              CALL CALLNUMNODE_F(IGM,N1,N2)     !CALL FRAME CONNECTIVITY NODE ... N1 N2
              NELEMENTGR(IFRAM,IGST) = IGM
              N1GROUP(IFRAM,IGST)    = N1
              N2GROUP(IFRAM,IGST)    = N2
          ENDIF
      ENDDO
C     ---------------------------------------------- 
      
C     ELEMENT GROUP PROVIDE
CC      NGROUP = MAXVAL(NGEOMETRY)
      NGROUP = NGGEO  !SONGSAK OCT2019
      
      DO I = 1,NGROUP
          INDEX = 1
          DO J = 1,NFRAM !CHANGE NELE TO NFRAM ... SONGSAK OCT2019
          IF (NELEMENTGR(J,I).NE.0) THEN
          NELEMENTGRNEW(INDEX,I) = NELEMENTGR(J,I)
          INDEX = INDEX + 1
          ENDIF
          ENDDO
      ENDDO
C     ONE ELEMENTS HAVE 2 CONNECTIVITY 
C     CONNECTIVITY 1 IN PILE GROUP SECTION   
      DO I = 1,NGROUP        
          INDEX = 1
          DO J = 1,NFRAM !CHANGE NELE TO NFRAM ... SONGSAK OCT2019
          IF (N1GROUP(J,I).NE.0) THEN
          N1GROUPNEW(INDEX,I) = N1GROUP(J,I)
          INDEX = INDEX + 1
          ENDIF
          ENDDO
      ENDDO
C     CONNECTIVITY 2 IN PILE GROUP SECTION     
      DO I = 1,NGROUP
          INDEX = 1
          DO J = 1,NFRAM !CHANGE NELE TO NFRAM ... SONGSAK OCT2019
          IF (N2GROUP(J,I).NE.0) THEN
          N2GROUPNEW(INDEX,I) = N2GROUP(J,I)
          INDEX = INDEX + 1
          ENDIF
          ENDDO
      ENDDO
C    FOR FINDING THE TOTAL NUMBER OF ELEMENT IN EACH PILE GROUP SECTION
      DO I = 1,NGROUP
          DO  J = 1,NFRAM !CHANGE NELE TO NFRAM ... SONGSAK OCT2019
          IF (N2GROUPNEW(J,I).NE.0)THEN
          IINDEX(I) = IINDEX(I) + 1
          ENDIF
          ENDDO
      ENDDO
      POTOTALGROUP = 0.0D0
      ! READ NTSOIL = TOTAL NUMBER OF SOIL INPUT DATA
      CALL READSOILDATA ("CALT",NTSOIL,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
     
      DO 100 NISN = 2,NTSOIL,2
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN-1,ISN1,PROPOUT)
           
           ! RE-CHECK ELEMENT AND NODE NUMBER FOR FINDING GROUP
           ! --------------------------------------------------------
           ! EXAMPLE. ELEMENT NO.1 ( CONNECTIVITY 7 8 ) IN GROUP 2
           ! ISN  = SOIL DATA 1  = 8
           ! ISN2 = SOIL DATA 2  = 7
           ! THIS WILL BE GIVE NUMBER OF GROUP IN 2 ( IGROUP = 2 )
           ! CONCLUDE : ELEMENT NO.1 HAS BEEN STORGE IN GROUP NO.2
           ! --------------------------------------------------------
           DO 7000 KGROUP = 1,NGROUP
              NGROUPINDEX = IINDEX(NGROUP) 
               NCORRECT1 = 0
               NCORRECT2 = 0
               DO JJ = 1,NGROUPINDEX
                   IF (N1GROUPNEW(JJ,KGROUP).EQ.ISN)THEN
                      NCORRECT1 = 1
                   ELSEIF (N1GROUPNEW(JJ,KGROUP).EQ.ISN1)THEN
                      NCORRECT2 = 1
                   ENDIF
                   IF (N2GROUPNEW(JJ,KGROUP).EQ.ISN1)THEN
                      NCORRECT2 = 1
                   ELSEIF (N2GROUPNEW(JJ,KGROUP).EQ.ISN)THEN
                      NCORRECT1 = 1
                   ENDIF
               ENDDO
                IF (NCORRECT1.EQ.1.AND.NCORRECT2.EQ.1)THEN
                   IGROUP = KGROUP
                    EXIT
                ENDIF
7000       CONTINUE
           
      ! FOR INPUT DATA
      NSNODE(1) = ISN
      NSNODE(2) = ISN1
      
      ! NODEABOVE FOR FINDING HIGHTEST NODE
      CALL NODEABOVE (NSNODE,NODEOUT)
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN,ISNT,PROPOUT)
      IF (ISNT.NE.NODEOUT) THEN
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN,ISNT,PROPOUT)
      ENDIF
      
CC      DO 20 I = 1,NELE
CC          IF (N1(I).EQ.ISN.AND.N2(I).EQ.ISN1)THEN
CC          NELEMENT = I
CC          EXIT
CC          ENDIF
CC          IF (N2(I).EQ.ISN.AND.N1(I).EQ.ISN1)THEN
CC          NELEMENT = I
CC          EXIT
CC          ENDIF
CC20    CONTINUE
      
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT LOOP ... OCT2019  

      DO 20 IGM = 1,NGIDM    
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
              CALL CALLNUMNODE_F (IGM,N1,N2) !CALL FRAME CONNECTIVITY NODE
              
              IF (N1.EQ.ISN.AND.N2.EQ.ISN1)THEN
                  NELEMENT = IGM
                  EXIT
              ENDIF
              IF (N2.EQ.ISN.AND.N1.EQ.ISN1)THEN
                  NELEMENT = IGM
                  EXIT
              ENDIF
          ENDIF
20    ENDDO
C     ----------------------------------------------------------------       
      
      ! CALL SECTION PROPERTIES
      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      ! CALL LENGTH OF THE ELEMENT
      CALL CALLENGTH (NELEMENT,VREW,0,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                ,ATMIDX,ATMIDY,ATMIDZ)
      
      ! START SOLUTION FOR CALCULATE SKIN FRICTION
      IF (SBC(2).NE.0.0D0) THEN 
          NUNIT  = SBC(14)
          NSOILY = SBC(15)
          C      = SBC(16)
          GRAMMA = SBC(17)
          DELTA  = SBC(18)
          ANQ    = SBC(19)
          BETA   = SBC(20)
          TFAC   = SBC(21)
          ZFAC   = SBC(22)
          AOPT   = SBC(23)
      
              
           ! CLAY
           IF (NSOILY.EQ.1) THEN
              FREE  = C/(GRAMMA*ALENGTH+POTOTALGROUP(IGROUP))
              IF (FREE.LE.1.0D0)THEN
                 ALPHA = 0.5D0*FREE**(-0.5D0)*FREE 
              ELSEIF (FREE.GT.1.0D0)THEN
                 ALPHA = 0.5D0*FREE**(-0.25D0)*FREE 
              ENDIF
              
              IF (ALPHA.GT.1.0D0) ALPHA = 1.0D0
              
              FRICTION  = ALPHA*C
              OUTAREA   = 2D0*3.141592654D0*(DPIPE/2D0)
              INAREA    = 2D0*3.141592654D0*((DPIPE-(TPIPE*2D0))/2D0)
              TOTALAREA = OUTAREA + INAREA
              
              QF        = TOTALAREA*FRICTION*ALENGTH
              
           ! SAND   
           ELSEIF (NSOILY.EQ.2) THEN
              CALL ELEMENTLENGTH (ISN1,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,"QZ")
              FRICTION  = BETA*(GRAMMA*ALENGTH+POTOTALGROUP(IGROUP))
              OUTAREA   = 2D0*3.141592654D0*(DPIPE/2D0)
              AINAREA    = 2D0*3.141592654D0*((DPIPE-(TPIPE*2D0))/2D0)
              TOTALAREA = OUTAREA + AINAREA
              
              QF        = TOTALAREA*FRICTION*ALENGTH
           ENDIF
           
           POTOTALGROUP(IGROUP) = GRAMMA*ALENGTH+POTOTALGROUP(IGROUP)
          
          
      ENDIF 
      
      IF (SBC(3).NE.0.0D0.OR.SBC(1).NE.0.0D0) THEN
C          ANQ  =  SBC(14)
C          C    =  SBC(15)
C          PO   =  SBC(17)
          
C          FRICTION  = ANQ*PO
          
C          OUTAREA   = 2D0*3.141592654D0*(DPIPE/2D0)
C          INAREA    = 2D0*3.141592654D0*((DPIPE-(TPIPE*2D0))/2D0)
C          TOTALAREA = OUTAREA + INAREA

C          QF        = QF + (TOTALAREA*FRICTION*ALENGTH)
          
C      ENDIF
      
      NCOUNT                            = NCOUNT +1
      ADATA(NCOUNT,IGROUP)              = QF
      NSOILELEMENT(NCOUNT,1)            = NELEMENT
      NSOILELEMENTPILE(NCOUNT,IGROUP)   = NELEMENT
      
      NPRINT(NELEMENT)                  = NELEMENT
      BDATA(NELEMENT)                   = ADATA(NCOUNT,IGROUP)
      
      ENDIF
C      ENDDO

C      ENDDO
      
100   CONTINUE
      
      ! FIND THE TOTAL GROUP NUMBER 
CC      NGROUP = MAXVAL(NGEOMETRY)
      NGROUP = NGGEO  !SONGSAK OCT2019
      
      DO I = 1,NGROUP
      INDEX = 1
      DO J = 1,NTSOIL
          IF (ADATA(J,I).NE.0.0)THEN
              ADATANEW(INDEX,I) = ADATA(J,I)
              INDEX = INDEX + 1 
          ENDIF
      ENDDO
      ENDDO
      
      DO I = 1,NGROUP
      INDEX = 1
      DO J = 1,NTSOIL
          IF (NSOILELEMENTPILE(J,I).NE.0.0)THEN
              NSOILELEMENTPILENEW(INDEX,I) = NSOILELEMENTPILE(J,I)
              NELEMENTINDEX(I) = INDEX
              INDEX = INDEX + 1 
          ENDIF
      ENDDO
      ENDDO
      
      
      DO I =1,NGROUP
      SUMTOTAL(I)             = SUM(ADATA(1:NCOUNT,I))
      ENDDO
      
      WRITE (171,*)
      WRITE (171,200) ILC
200   FORMAT ("*********************** LOAD CASE",1X,I5,"**********************") 
      WRITE (171,201) 
201   FORMAT ("--------- SKIN FRICTION CAPACITY API 2A-WSD ( PULLOUT LOAD ) ---------")
      
      WRITE (171,202)
202   FORMAT ("NUM. ELEMENT       SKIN FRICTION(KN)      EXTERNAL FORCE(UNIT)")
      
      DO J = 1,NGROUP
      IF (NELEMENTINDEX(J).EQ.0) GOTO 1000
      CALL PRNGPHD(NAMEG,LENGTH,J)
      WRITE (171,270) NAMEG(1:LENGTH)
270   FORMAT ("--",A,"--")
      DO I = 1,NELEMENTINDEX(J) 
          DO 1100 IGM = 1,NGIDM  !SONGSAK CHANGE HERE FROM [K =1,NSN] TO [IGM = 1,NGIDM] ... OCT2019
              
CC             IF (NSOILELEMENTPILENEW(I,J).EQ.IELEMENTGEOMETRY(K))THEN
CC             INPUTINDEX = K
C             SUPPRESS ABOVE STATEMENT AND REPLACE WITH NEXT STATEMENT ... SONGSAK OCT2019            
             IF (NSOILELEMENTPILENEW(I,J).EQ.IGM) THEN
             INPUTINDEX = IGM
             
             EXIT
             ENDIF
1100      CONTINUE
      CALL AXIALFORCE (INPUTINDEX,AXIALFORCEOUT)
      WRITE (171,203) NSOILELEMENTPILENEW(I,J),ADATANEW(I,J),AXIALFORCEOUT
203   FORMAT (I5,13X,E12.5,12X,E12.5)
      ENDDO
      WRITE (171,204) SUMTOTAL(J)
C      WRITE (171,210)
1000  ENDDO
      

204   FORMAT ("TOTAL SKIN FRICTION(KN) = ",1X,E12.5)   
      
      WRITE (LFPRIN,205)  
205   FORMAT ('Result "(Pile-Soil) Skin Friction"  "Pile Result"      1  Scalar  OnGaussPoints  "XFrame2"')
 
      WRITE (LFPRIN,206) 
206   FORMAT ('ComponentNames "Skin Friction"')  
      
      WRITE (LFPRIN,207) 
207   FORMAT ('Values')  
      
      DO I = 1,NFRAM !CHANGE NELE TO NFRAM ... SONGSAK OCT2019
          DUMMY = 0.0D0
          IF (NPRINT(I).EQ.0) WRITE (LFPRIN,208) I,DUMMY
          IF (NPRINT(I).NE.0) WRITE (LFPRIN,208) I,BDATA(I)
208       FORMAT (I5,5X,E12.5)
          IF (NPRINT(I).EQ.0) WRITE (LFPRIN,209) DUMMY
          IF (NPRINT(I).NE.0) WRITE (LFPRIN,209) BDATA(I)
209       FORMAT (5X,5X,E12.5)
      ENDDO
      
      WRITE (LFPRIN,250) 
250   FORMAT ("End Values")    
      
       DO J = 1,NGROUP
           INDEX_OUT = 1D0
           DO I = 1,NELEMENTINDEX(J) 
               NSNODE    = 0.
               
CC               NSNODE(1) = N1(NSOILELEMENTPILENEW(I,J)) 
CC               NSNODE(2) = N2(NSOILELEMENTPILENEW(I,J))
C             SONGSAK OCT2019
               CALL CALLNUMNODE_F (NSOILELEMENTPILENEW(I,J),N1,N2)
               NSNODE(1) = N1
               NSNODE(2) = N2
               
               CALL NODEABOVE (NSNODE,NODEOUT)
               NODEINDEX(INDEX_OUT,J)       =  NSNODE(1)
               INDEX_OUT = INDEX_OUT + 1
               NODEINDEX(INDEX_OUT,J)       =  NSNODE(2)
               INDEX_OUT = INDEX_OUT + 1
               INDEX_TOTOAL_OUT(J) = INDEX_OUT
               
               !ADDSTIFFNESS1(INDEX_OUT,J) = SSTIFF(1,NSNODE(1))
               !ADDSTIFFNESS2(INDEX_OUT,J) = SSTIFF(2,NSNODE(1)) 
               !ADDSTIFFNESS3(INDEX_OUT,J) = SSTIFF(3,NSNODE(1))
               !ADDSTIFFNESS1(INDEX_OUT+1,J) = SSTIFF(1,NSNODE(2))
               !ADDSTIFFNESS2(INDEX_OUT+1,J) = SSTIFF(2,NSNODE(2))
               !ADDSTIFFNESS3(INDEX_OUT+1,J) = SSTIFF(3,NSNODE(2))
               !NODEINDEX(INDEX_OUT,J)       =  NSNODE(1)
               !NODEINDEX(INDEX_OUT+1,J)     =  NSNODE(1)
               !INDEX_OUT = INDEX_OUT + 2
                   
               !IF (I.EQ.1D0)THEN
               !    IF (NODEOUT.EQ.NSNODE(2))THEN
               !    ADDSTIFFNESS1(I,J) = SSTIFF(1,NSNODE(1))
               !    ADDSTIFFNESS2(I,J) = SSTIFF(2,NSNODE(1))
               !    ADDSTIFFNESS3(I,J) = SSTIFF(3,NSNODE(1))
               !    ADDSTIFFNESS1(I+1,J) = SSTIFF(1,NSNODE(2))
               !    ADDSTIFFNESS2(I+1,J) = SSTIFF(2,NSNODE(2))
               !    ADDSTIFFNESS3(I+1,J) = SSTIFF(3,NSNODE(2))
               !    NODEINDEX(I,J)       =  NSNODE(1)
               !    NODEINDEX(I+1,J)     =  NSNODE(2)
               !    ELSEIF (NODEOUT.EQ.NSNODE(1))THEN
               !    ADDSTIFFNESS1(I,J) = SSTIFF(1,NSNODE(2))
               !    ADDSTIFFNESS2(I,J) = SSTIFF(2,NSNODE(2))
               !    ADDSTIFFNESS3(I,J) = SSTIFF(3,NSNODE(2))
               !    ADDSTIFFNESS1(I+1,J) = SSTIFF(1,NSNODE(1))
               !    ADDSTIFFNESS2(I+1,J) = SSTIFF(2,NSNODE(1))
               !    ADDSTIFFNESS3(I+1,J) = SSTIFF(3,NSNODE(1))
               !    NODEINDEX(I,J)       =  NSNODE(2)
               !    NODEINDEX(I+1,J)     =  NSNODE(1)
               !    ENDIF
               !ELSEIF (I.NE.1D0)THEN
               !    ADDSTIFFNESS1(I+1,J) = SSTIFF(1,NODEOUT)
               !    ADDSTIFFNESS2(I+1,J) = SSTIFF(2,NODEOUT)
               !    ADDSTIFFNESS3(I+1,J) = SSTIFF(3,NODEOUT)
               !    NODEINDEX(I+1,J)       =  NODEOUT
              !ENDIF
           ENDDO
           
       ENDDO
       
            
       DO J = 1,NGROUP
           DO I = 1,INDEX_TOTOAL_OUT(J) - 1
               NNODE_INDEX = NODEINDEX(I,J)
               OPTION_INDEX = 0
               DO II = 1,INDEX_TOTOAL_OUT(J) - 1
               IF (NNODE_INDEX.EQ.NODEINDEX(II,J).AND.OPTION_INDEX.EQ.0)THEN
                   OPTION_INDEX = 1
               ELSEIF(NNODE_INDEX.EQ.NODEINDEX(II,J).AND.OPTION_INDEX.EQ.1)THEN
                   NODEINDEX(I,J) = 0.0D0
               ENDIF
               ENDDO
           ENDDO
       ENDDO
       
       INDEX_OUT_G = 0.
       DO J = 1,NGROUP
           INDEX_OUT_G(J) = 1D0
           DO I = 1,INDEX_TOTOAL_OUT(J) - 1
               IF (NODEINDEX(I,J).NE.0)THEN
               NOUT_NODE(INDEX_OUT_G(J),J) = NODEINDEX(I,J)
               NOUT_NODE_SUPER(INDEX_OUT_G(J),J) = NODEINDEX(I,J)
               INDEX_OUT_G(J) = INDEX_OUT_G(J) + 1
               ENDIF
           ENDDO
       ENDDO
       
       DO J = 1,NGROUP
           INDEX_OUT = 1
           DO I = 1,INDEX_OUT_G(J)- 1
              NMIN_OUT =  MINVAL(NOUT_NODE(1:INDEX_OUT_G(J)-1,J))
              NMIN_OUT_LOC =  MINLOC(NOUT_NODE(1:INDEX_OUT_G(J)-1,J),DIM=1)
              NOUT_NODE(NMIN_OUT_LOC,J) = 10000D0
              NNCODE_OUT(INDEX_OUT,J) = NMIN_OUT
              INDEX_OUT = INDEX_OUT + 1
           ENDDO
       ENDDO
          
      
       OPEN(UNIT=7110  ,FILE="INTERGER_TO_STRING"      ,STATUS='UNKNOWN'    )   
       WRITE (7110,*) ILC
       REWIND (7110)
       READ (7110,*)  NUMBERFILE
       CLOSE(7110)
       !DIRPRIMARY     = "PILE-SOIL INTERACTION OUTPUT"//"/"//"SOIL_STIFFNESS"//TRIM(NUMBERFILE)//".DAT"
       !OPEN(UNIT=7100  ,FILE=DIRPRIMARY      ,STATUS='UNKNOWN'    )
       OPEN(UNIT=7050  ,FILE="PILE-SOIL INTERACTION OUTPUT/Normal spring stiffness.dat"     ,STATUS='UNKNOWN'    )
       WRITE (7050,253) ILC,LCS
253   FORMAT ("NODE     K-X           K-Y           K-Z --------- LOADCASE",I5,"/",I5,"---------")
       DO J = 1,NGROUP
       CALL PRNGPHD(NAMEG,LENGTH,J)
       WRITE (7050,260) NAMEG(1:LENGTH)
260    FORMAT ("--",A,"--")
C       WRITE (7050,251)
C251    FORMAT ("TOTAL OF NODE NUMBER") 
       WRITE (7050,252) NELEMENTINDEX(J)+1
252    FORMAT (I5)     
       !DO I =1,NELEMENTINDEX(J)+1
       !WRITE (7050,254) NODEINDEX(I,J),ADDSTIFFNESS1(I,J),ADDSTIFFNESS2(I,J),ADDSTIFFNESS3(I,J)
       DO I =1,INDEX_OUT_G(J) - 1
       WRITE (7050,254) NNCODE_OUT(I,J),SSTIFF(1,NNCODE_OUT(I,J)),SSTIFF(2,NNCODE_OUT(I,J)),SSTIFF(3,NNCODE_OUT(I,J))
254    FORMAT (I5,2X,E12.5,2X,E12.5,2X,E12.5)
       ENDDO
       ENDDO
       WRITE (7050,"")
       WRITE (7050,255)
255    FORMAT ("**********************")
       WRITE (7050,"")
        
      ADATA   = 0.0D0
      BEARING   = 0.D0
      NCOUNT  = 0D0
      POTOTALGROUP = 0.0D0
      DO 110 NISN = 2,NTSOIL,2
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT) 
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN-1,ISP,PROPOUT)
      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      
           DO 7001 KGROUP = 1,NGROUP
              NGROUPINDEX = IINDEX(NGROUP) 
               NCORRECT1 = 0
               NCORRECT2 = 0
               DO JJ = 1,NGROUPINDEX
                   IF (N1GROUPNEW(JJ,KGROUP).EQ.ISN)THEN
                      NCORRECT1 = 1
                   ELSEIF (N1GROUPNEW(JJ,KGROUP).EQ.ISP)THEN
                      NCORRECT2 = 1
                   ENDIF
                   IF (N2GROUPNEW(JJ,KGROUP).EQ.ISP)THEN
                      NCORRECT2 = 1
                   ELSEIF (N2GROUPNEW(JJ,KGROUP).EQ.ISN)THEN
                      NCORRECT1 = 1
                   ENDIF
               ENDDO
                IF (NCORRECT1.EQ.1.AND.NCORRECT2.EQ.1)THEN
                   IGROUP = KGROUP
                    EXIT
                ENDIF
7001           CONTINUE
C      IF (SBC(2).NE.0.0D0) THEN 
C           C      = SBC(10)
C           AK     = SBC(11)
C           PO     = SBC(12)
C           DELTA  = SBC(13)
C           ANQ    = SBC(18)
C           BETA   = SBC(19)
C           ! CLAY
C           IF (SBC(9).EQ.1) THEN
C              BEARING   = 9D0*C
C              QP        = AREA*BEARING
C           ! SAND   
C           ELSEIF (SBC(9).EQ.2) THEN
C              BEARING   = DELTA*PO
C              QP        = AREA*BEARING
C           ENDIF
C          
C      ENDIF
      
      IF (SBC(3).NE.0.0D0.OR.SBC(1).NE.0.0D0) THEN
          NUNIT  = SBC(14)
          NSOILY = SBC(15)
          C      = SBC(16)
          GRAMMA = SBC(17)
          DELTA  = SBC(18)
          ANQ    = SBC(19)
          BETA   = SBC(20)
          TFAC   = SBC(21)
          ZFAC   = SBC(22)
          AOPT   = SBC(23)
          
          
          CALL ELEMENTLENGTH (ISP,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,"QZ")
          IF (NSOILY.EQ.1) THEN ! CLAY
          BEARING   = 9D0*C
          QP        = AREA*BEARING
          ELSEIF (NSOILY.EQ.2)THEN ! SAND
          POTOTALGROUP(KGROUP)   = GRAMMA*ALENGTH + POTOTALGROUP(KGROUP)
          
          BEARING   = ANQ*POTOTALGROUP(KGROUP)
          QP        = AREA*BEARING
          ENDIF
      ENDIF
      
      NCOUNT                        = NCOUNT +1
      BEARINGDATA(NCOUNT,IGROUP)    = QP    ! TEXT 
      BDATA(NCOUNT)                 = QP    ! GID  
      NSOILELEMENT1(NCOUNT,IGROUP)  = ISN   ! TEXT 
      NSOILELEMENT2(NCOUNT,IGROUP)  = ISP   ! TEXT 
      NSOILELEMENT(NCOUNT,1)        = ISN   ! GID  
      NSOILELEMENT(NCOUNT,2)        = ISP   ! GID  
110   CONTINUE
      
CC      NGROUP = MAXVAL(NGEOMETRY)
      NGROUP = NGGEO  !SONGSAK OCT2019
      
      DO I = 1,NGROUP
      INDEX = 1
      DO J = 1,NTSOIL
          IF (BEARINGDATA(J,I).NE.0.0)THEN
              BEARINGDATANEW(INDEX,I) = BEARINGDATA(J,I)
              INDEX = INDEX + 1 
          ENDIF
      ENDDO
      ENDDO
      
      DO I = 1,NGROUP
      INDEX = 1
      DO J = 1,NTSOIL
          IF (NSOILELEMENT1(J,I).NE.0.0)THEN
              NSOILELEMENT1NEW(INDEX,I) = NSOILELEMENT1(J,I)
              INDEX = INDEX + 1 
          ENDIF
      ENDDO
      ENDDO
      
      DO I = 1,NGROUP
      INDEX = 1
      DO J = 1,NTSOIL
          IF (NSOILELEMENT2(J,I).NE.0.0)THEN
              NSOILELEMENT2NEW(INDEX,I) = NSOILELEMENT2(J,I)
              NELEMENTINDEXBEARING(I)   = INDEX
              INDEX = INDEX + 1 
          ENDIF
      ENDDO
      ENDDO
      
      NPRINT  = 0.
      CALL ADDSOILSTIFF (SD,NTSOIL,NODETOPMAX,NODEBOTTOMAX,KCOUNBOTTOM,KCOUNTOP)
      DO I = 1,NCOUNT
         IS1 =NSOILELEMENT(I,1)
         IS2 =NSOILELEMENT(I,2) 
         
         NTOP1   = 0.
         NTOP2   = 0.
         DO 500 J = 1,KCOUNTOP
             IF (NODETOPMAX(J).EQ.IS1) THEN
                NTOP1      = 1 
                NPRINT(I)  = IS1
                EXIT
             ENDIF
             IF (NODETOPMAX(J).EQ.IS2) THEN
                NTOP2       = 1 
                NPRINT(I)   = IS2
                EXIT
             ENDIF
500      CONTINUE
        
         IF (NTOP1.EQ.1)THEN
             ADATA(NPRINT(I),NGROUP) = BDATA(I)/2D0 
             ADATA(IS2,NGROUP)       = BDATA(I)
         ELSEIF (NTOP2.EQ.1)THEN
             ADATA(NPRINT(I),NGROUP) = BDATA(I)/2D0 
             ADATA(IS1,NGROUP)       = BDATA(I)
         ENDIF 
         
         IF (NTOP1.EQ.0.AND.NTOP2.EQ.0)THEN
             
         PREA    = 0.D0
         PREB    = 0.D0
         PREA    = ADATA(IS1,NGROUP)
         PREB    = ADATA(IS2,NGROUP)
         
         ANEWDATAA = BDATA(I)
         ANEWDATAB = BDATA(I)
         
         IF (PREA.EQ.0.0D0) THEN
            ADATA(IS1,NGROUP) = ANEWDATAA
         ELSEIF (PREA.NE.0.0D0) THEN
            IF (ANEWDATAA.GT.PREA)THEN
                ADATA(IS1,NGROUP) = PREA!ANEWDATAA
            ELSEIF (ANEWDATAA.LE.PREA)THEN
                ADATA(IS1,NGROUP) = ANEWDATAA
            ENDIF
         ENDIF 
         
         IF (PREB.EQ.0.0D0) THEN
            ADATA(IS2,NGROUP) = ANEWDATAB
         ELSEIF (PREA.NE.0.0D0) THEN
            IF (ANEWDATAB.GT.PREB)THEN
                ADATA(IS2,NGROUP) = PREB
            ELSEIF (ANEWDATAB.LE.PREB)THEN
                ADATA(IS2,NGROUP) = ANEWDATAB
            ENDIF
         ENDIF 
         
         
         ENDIF
          
        ENDDO
      
      
      WRITE (171,210)
      WRITE (171,210)
210   FORMAT ("")
      
      WRITE (171,211) 
211   FORMAT ("--------- END BEARING CAPACITY API 2A-WSD ---------")
      
      WRITE (171,212)
212   FORMAT ("NUM. NODE       END BEARING(KN)")
      
      DO J = 1,NGROUP
      IF (NELEMENTINDEXBEARING(J).EQ.0) GOTO 1001
      CALL PRNGPHD(NAMEG,LENGTH,J)
      WRITE (171,271) NAMEG(1:LENGTH)
271   FORMAT ("--",A,"--")
      DO I = 1,NELEMENTINDEXBEARING(J)
      WRITE (171,213) NSOILELEMENT1NEW(I,J),NSOILELEMENT2NEW(I,J),BEARINGDATANEW(I,J)
213   FORMAT (I3,"  -",I3,13X,E12.5)
      ENDDO
      AMAXBEARING = MAXVAL(BEARINGDATANEW(:,J))
      WRITE (171,214) AMAXBEARING
1001  CC = 1
      ENDDO
      
214   FORMAT ("MAXIMUM END BEARING CAPACITY =",E12.5)
      
      WRITE (LFPRIN,"")
      
      WRITE (LFPRIN,215)  
215   FORMAT ('Result "Pile-Soil End Bearing"  "Pile Result"      1  Scalar OnNodes')
 
      WRITE (LFPRIN,216) 
216   FORMAT ('ComponentNames "End Bearing"')  
      
      WRITE (LFPRIN,217) 
217   FORMAT ('Values')  
      
      DO I = 1,NSN
          WRITE (LFPRIN,218) I,ADATA(I,NGROUP)
218       FORMAT (I5,5X,E12.5)
      ENDDO
      
      
      WRITE (LFPRIN,219) 
219   FORMAT ("End Values")    
      

      
C     SONGSAK OCT2019      
CC      DEALLOCATE(NCONNECT1,NCONNECT2)
CC      DEALLOCATE(NELEMENTGR,NELEMENTGRNEW)
CC      DEALLOCATE(N1GROUP,N2GROUP)
CC      DEALLOCATE(N1GROUPNEW,N2GROUPNEW)
CC      DEALLOCATE(ADATA,BDATA,NSOILELEMENT,NPRINT)
CC      DEALLOCATE(NSOILELEMENT1,NSOILELEMENT2)
CC      DEALLOCATE(NSOILELEMENT1NEW,NSOILELEMENT2NEW)
CC      DEALLOCATE(ADATANEW,BDATANEW)
CC      DEALLOCATE(BEARINGDATA,BEARINGDATANEW)
CC      DEALLOCATE(NSOILELEMENTPILE,NSOILELEMENTPILENEW)
CC      DEALLOCATE(ADDSTIFFNESS1,ADDSTIFFNESS2,ADDSTIFFNESS3)
CC      DEALLOCATE(NODEINDEX)  
      
C      WRITE (171,205)
C      WRITE (171,205)

      
C      DO I = 1,NSOIL
C      WRITE (171,211) I
C      WRITE (171,212)
C211   FORMAT ("--------- SOIL LAYER",2I,"---------")   
C212   FORMAT ("SKIN FRICTION =")  
C      ENDDO
      
      END
C     ----------------------------------------
      SUBROUTINE ADDSOILSTIFF (SD,NTSOIL,NTOP,NBOTTOM,KCOUNBOTTOM,KCOUNTOP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

      COMMON /MGRAV/ NGRAV 
      DIMENSION SD(6,NSN),SBC(37)
      DIMENSION NODESOIL(NTSOIL)
      DIMENSION AXNODE(NTSOIL),AYNODE(NTSOIL),AZNODE(NTSOIL)
      DIMENSION NTOP(10000),NBOTTOM(10000)
      
      AXNODE = 0.0D0
      AYNODE = 0.0D0
      AZNODE = 0.0D0
      DO NISN = 1,NTSOIL
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
      
      NODESOIL(NISN)  = ISN
      
CC      AXNODE(NISN)    = AX(ISN)
CC      AYNODE(NISN)    = AY(ISN)
CC      AZNODE(NISN)    = AZ(ISN)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AXNODE(NISN),1,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AYNODE(NISN),2,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZNODE(NISN),3,ISN,0)  !GETTING NODAL COORDINATE
      
      ENDDO
      
      NAXMAX = MAXLOC(AXNODE,DIM=1)
      NAYMAX = MAXLOC(AYNODE,DIM=1)
      NAZMAX = MAXLOC(AZNODE,DIM=1)
      
      NAXMIN = MINLOC(AXNODE,DIM=1)
      NAYMIN = MINLOC(AYNODE,DIM=1)
      NAZMIN = MINLOC(AZNODE,DIM=1)
      
      IF (NGRAV.EQ.1)THEN
         NODETOP    = NODESOIL(NAXMAX)
         NODEBOTTOM = NODESOIL(NAXMIN)
      ELSEIF (NGRAV.EQ.2)THEN
         NODETOP    = NODESOIL(NAYMAX)
         NODEBOTTOM = NODESOIL(NAYMIN)
      ELSEIF (NGRAV.EQ.3)THEN
         NODETOP    = NODESOIL(NAZMAX)
         NODEBOTTOM = NODESOIL(NAZMIN)  
      ENDIF
      
      KCOUNTOP    = 0.0D0
      KCOUNBOTTOM = 0.0D0
      DO I =1,NTSOIL
          IF (NGRAV.EQ.1)THEN

 
          CALL RELFILL('@XYZ',AXN_T,1,NODETOP    ,0)  !GETTING NODAL COORDINATE  
          CALL RELFILL('@XYZ',AXN_S,1,NODESOIL(I),0)  !GETTING NODAL COORDINATE            
              
C          IF (AX(NODETOP).EQ.AX(NODESOIL(I)))THEN
          IF (AXN_T.EQ.AXN_S) THEN !SONGSAK OCT2019
              KCOUNTOP       = KCOUNTOP + 1
              NTOP(KCOUNTOP) = NODESOIL(I)
          ENDIF
          
          CALL RELFILL('@XYZ',AXN_B,1,NODEBOTTOM ,0)  !GETTING NODAL COORDINATE  
          CALL RELFILL('@XYZ',AXN_S,1,NODESOIL(I),0)  !GETTING NODAL COORDINATE
          
C          IF (AX(NODEBOTTOM).EQ.AX(NODESOIL(I)))THEN
          IF (AXN_B.EQ.AXN_S) THEN !SONGSAK OCT2019
               KCOUNBOTTOM          = KCOUNBOTTOM+ 1
               NBOTTOM(KCOUNBOTTOM) = NODESOIL(I)
          ENDIF    
          
          ELSEIF (NGRAV.EQ.2)THEN
              
          CALL RELFILL('@XYZ',AYN_T,2,NODETOP    ,0)  !GETTING NODAL COORDINATE  
          CALL RELFILL('@XYZ',AYN_S,2,NODESOIL(I),0)  !GETTING NODAL COORDINATE  
          
C          IF (AY(NODETOP).EQ.AY(NODESOIL(I)))THEN
          IF (AYN_T.EQ.AYN_S) THEN !SONGSAK OCT2019
              KCOUNTOP       = KCOUNTOP + 1
              NTOP(KCOUNTOP) = NODESOIL(I)
          ENDIF
          
          CALL RELFILL('@XYZ',AYN_B,2,NODEBOTTOM ,0)  !GETTING NODAL COORDINATE  
          CALL RELFILL('@XYZ',AYN_S,2,NODESOIL(I),0)  !GETTING NODAL COORDINATE
          
C          IF (AY(NODEBOTTOM).EQ.AY(NODESOIL(I)))THEN
          IF (AYN_B.EQ.AYN_S) THEN !SONGSAK OCT2019
               KCOUNBOTTOM          = KCOUNBOTTOM+ 1
               NBOTTOM(KCOUNBOTTOM) = NODESOIL(I)
          ENDIF
          
          ELSEIF (NGRAV.EQ.3)THEN
              
          CALL RELFILL('@XYZ',AZN_T,3,NODETOP    ,0)  !GETTING NODAL COORDINATE  
          CALL RELFILL('@XYZ',AZN_S,3,NODESOIL(I),0)  !GETTING NODAL COORDINATE  
          
C          IF (AZ(NODETOP).EQ.AZ(NODESOIL(I)))THEN
          IF (AZN_T.EQ.AZN_S) THEN !SONGSAK OCT2019
              KCOUNTOP       = KCOUNTOP + 1
              NTOP(KCOUNTOP) = NODESOIL(I)
          ENDIF
          
          CALL RELFILL('@XYZ',AZN_B,3,NODEBOTTOM ,0)  !GETTING NODAL COORDINATE  
          CALL RELFILL('@XYZ',AZN_S,3,NODESOIL(I),0)  !GETTING NODAL COORDINATE
          
C          IF (AZ(NODEBOTTOM).EQ.AZ(NODESOIL(I)))THEN
          IF (AZN_B.EQ.AZN_S) THEN !SONGSAK OCT2019
               KCOUNBOTTOM          = KCOUNBOTTOM+ 1
               NBOTTOM(KCOUNBOTTOM) = NODESOIL(I)
          ENDIF
          ENDIF
          
      ENDDO
      
      END
C     ----------------------------------------
      SUBROUTINE EACHPILE 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
CC      COMMON /STOREN/ IELEMENT(100000),N1(100000),N2(100000)

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NA,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /MGRAV/ NGRAV
      DIMENSION NCONNET(NSN,NSN),NTCONNET(NSN,NSN),ICONNECTION(NSN,NSN)
      DIMENSION NPOINT(NSN),NPLAN(NSN)
      
CC      DO K = 1,NSN
CC        NTCO = 1
CC        NODE = K
CC       DO J = 1,NSN
CC        DO 100 I = 1 ,NELE
CC          IF (NODE.EQ.N1(I))THEN
CC             NCONNET(NTCO,K)     = N2(I)  
CC             NODE                = NCONNET(NTCO,K)
CC             NTCONNET(NODE,K)    = NTCONNET(NODE,K) + 1
CC             EXIT
CC          ENDIF 
CC          IF (NODE.EQ.N2(I))THEN
CC             NCONNET(NTCO,K)    = N1(I) 
CC             NODE               = NCONNET(NTCO,K)
CC             NTCONNET(NODE,K)   = NTCONNET(NODE,K) + 1
CC             EXIT
CC          ENDIF     
CC100   CONTINUE
CC      NTCO = NTCO + 1
CC       ENDDO
CC      ENDDO
      
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT LOOP ... OCT2019  
      
C	CALL TOTAL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NGDAT,NGIDM,1) 
      
      DO K = 1,NSN
      NTCO = 1
      NODE = K
      DO J = 1,NSN
           
      DO 100 IGM = 1,NGIDM    
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
              CALL CALLNUMNODE_F (IGM,N1,N2) !CALL FRAME CONNECTIVITY NODE
              
              IF (NODE.EQ.N1)THEN
                 NCONNET(NTCO,K)     = N2
                 NODE                = NCONNET(NTCO,K)
                 NTCONNET(NODE,K)    = NTCONNET(NODE,K) + 1
                 EXIT
              ENDIF
      
              IF (NODE.EQ.N2)THEN
                 NCONNET(NTCO,K)    = N1
                 NODE               = NCONNET(NTCO,K)
                 NTCONNET(NODE,K)   = NTCONNET(NODE,K) + 1
                 EXIT
              ENDIF
          ENDIF
100    ENDDO
       
      NTCO = NTCO + 1
      ENDDO
      ENDDO  
C     ----------------------------------------------------------------  
      
      DO IK = 1,NSN
          NTTCO = 1
          DO IJ =1,NSN
              IF (NTCONNET(IJ,IK).NE.0)THEN
                  ICONNECTION(NTTCO,IK) = IJ
                  NTTCO = NTTCO + 1
              ENDIF
          ENDDO
      ENDDO

      DO IT = 1,NSN
       NPOINT(1:NSN) = ICONNECTION(1:NSN,IT)
       
          DO IP = 1,NSN
              DO 200 IC = 1,NSN
                 IF (NPOINT(IC).EQ.ICONNECTION(IC,IP))THEN
                 NPLAN(IP) = IT
                 EXIT
                 ENDIF
200           CONTINUE
          ENDDO
          
      ENDDO
      
      END
C     ----------------------------------------
      SUBROUTINE ELEMENTLENGTH (ISP,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,OPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*2 OPT

      COMMON /MGRAV/ NGRAV
      
      DIMENSION NTOP(1000),NBOTTOM(1000)
      
      
      IF (OPT.EQ."PY")THEN
          
CC      AX1 = AX(ISP)
CC      AY1 = AY(ISP)
CC      AZ1 = AZ(ISP)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX1,1,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY1,2,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ1,3,ISP,0)  !GETTING NODAL COORDINATE	
      
CC      AX2 = AX(ISN)
CC      AY2 = AY(ISN)
CC      AZ2 = AZ(ISN)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX2,1,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY2,2,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ2,3,ISN,0)  !GETTING NODAL COORDINATE	
      
      CALL ADDSOILSTIFF (SD,NTSOIL,NTOP,NBOTTOM,KCOUNBOTTOM,KCOUNTOP)
      
      CALL RELFILL('@XYZ',AXTOP,1,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AYTOP,2,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZTOP,3,NTOP(1),0)  !GETTING NODAL COORDINATE	
      
       IF (NGRAV.EQ.1)THEN
           ALENGTHPY = ABS(AXTOP - (AX1+AX2)/2D0) ! UPDATE 2022
       ELSEIF (NGRAV.EQ.2)THEN
           ALENGTHPY = ABS(AYTOP - (AY1+AY2)/2D0) ! UPDATE 2022
       ELSEIF (NGRAV.EQ.3)THEN
           ALENGTHPY = ABS(AZTOP - (AZ1+AZ2)/2D0) ! UPDATE 2022
       ENDIF
      
      CALL NODETOELEMENT (ISN,NELEMENT,NISN,DEPTH,NNOUT,"ABOVE") 
      
CC      AX1 = AX(ISP)
CC      AY1 = AY(ISP)
CC      AZ1 = AZ(ISP)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX1,1,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY1,2,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ1,3,ISP,0)  !GETTING NODAL COORDINATE	
      
CC      AX2 = AX(ISN)
CC      AY2 = AY(ISN)
CC      AZ2 = AZ(ISN)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX2,1,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY2,2,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ2,3,ISN,0)  !GETTING NODAL COORDINATE	

      
      IF (NGRAV.EQ.1)THEN
          IF (AX1.GT.AX2) NODE = ISN
          IF (AX1.LE.AX2) NODE = ISP
          
CC          AX3 = ABS(AX(NODE))
CC          AY3 = ABS(AY(NODE))
CC          AZ3 = ABS(AZ(NODE))
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX3,1,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY3,2,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ3,3,NODE,0)  !GETTING NODAL COORDINATE
          AX3 = ABS(AX3)
          AY3 = ABS(AY3)
          AZ3 = ABS(AZ3)	
          
CC          AX4 = ABS(AX(NTOP(1)))
CC          AY4 = ABS(AY(NTOP(1)))
CC          AZ4 = ABS(AZ(NTOP(1)))
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX4,1,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY4,2,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ4,3,NTOP(1),0)  !GETTING NODAL COORDINATE
          AX4 = ABS(AX4)
          AY4 = ABS(AY4)
          AZ4 = ABS(AZ4)	
          
          ALENGTHPP  = SQRT((AX3-AX4)**2+(AY3-AY4)**2+(AZ3-AZ4)**2)
      ELSEIF (NGRAV.EQ.2)THEN
          
          IF (AY1.GT.AY2) NODE = ISN
          IF (AY1.LE.AY2) NODE = ISP
          
CC          AX3 = ABS(AX(NODE))
CC          AY3 = ABS(AY(NODE))
CC          AZ3 = ABS(AZ(NODE))
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX3,1,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY3,2,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ3,3,NODE,0)  !GETTING NODAL COORDINATE
          AX3 = ABS(AX3)
          AY3 = ABS(AY3)
          AZ3 = ABS(AZ3)	
          
CC          AX4 = ABS(AX(NTOP(1)))
CC          AY4 = ABS(AY(NTOP(1)))
CC          AZ4 = ABS(AZ(NTOP(1)))
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX4,1,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY4,2,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ4,3,NTOP(1),0)  !GETTING NODAL COORDINATE
          AX4 = ABS(AX4)
          AY4 = ABS(AY4)
          AZ4 = ABS(AZ4)	
          
          ALENGTHPP  = SQRT((AX3-AX4)**2+(AY3-AY4)**2+(AZ3-AZ4)**2)
          
      ELSEIF (NGRAV.EQ.3)THEN
          
          IF (AZ1.GT.AZ2) NODE = ISN
          IF (AZ1.LE.AZ2) NODE = ISP
          
CC          AX3 = ABS(AX(NODE))
CC          AY3 = ABS(AY(NODE))
CC          AZ3 = ABS(AZ(NODE))
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX3,1,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY3,2,NODE,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ3,3,NODE,0)  !GETTING NODAL COORDINATE
          AX3 = ABS(AX3)
          AY3 = ABS(AY3)
          AZ3 = ABS(AZ3)	
          
CC          AX4 = ABS(AX(NTOP(1)))
CC          AY4 = ABS(AY(NTOP(1)))
CC          AZ4 = ABS(AZ(NTOP(1)))
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX4,1,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY4,2,NTOP(1),0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ4,3,NTOP(1),0)  !GETTING NODAL COORDINATE
          AX4 = ABS(AX4)
          AY4 = ABS(AY4)
          AZ4 = ABS(AZ4)
          
          ALENGTHPP  = SQRT((AX3-AX4)**2+(AY3-AY4)**2+(AZ3-AZ4)**2)
          
      ENDIF

      ELSEIF (OPT.EQ."TZ")THEN
          
CC      AX1 = AX(ISP)
CC      AY1 = AY(ISP)
CC      AZ1 = AZ(ISP)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX1,1,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY1,2,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ1,3,ISP,0)  !GETTING NODAL COORDINATE	
      
CC      AX2 = AX(ISN)
CC      AY2 = AY(ISN)
CC      AZ2 = AZ(ISN)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX2,1,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY2,2,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ2,3,ISN,0)  !GETTING NODAL COORDINATE	
       
       ALENGTH  = SQRT((AX1-AX2)**2+(AY1-AY2)**2+(AZ1-AZ2)**2)
      
      ELSEIF (OPT.EQ."QZ")THEN

CC      AX1 = AX(ISP)
CC      AY1 = AY(ISP)
CC      AZ1 = AZ(ISP)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX1,1,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY1,2,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ1,3,ISP,0)  !GETTING NODAL COORDINATE	
      
CC      AX2 = AX(ISN)
CC      AY2 = AY(ISN)
CC      AZ2 = AZ(ISN)
C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE BLOCK ... REPLACE BY NEXT BLOCK ... OCT2019    
	CALL RELFILL('@XYZ',AX2,1,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AY2,2,ISN,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AZ2,3,ISN,0)  !GETTING NODAL COORDINATE	
       
       ALENGTH  = SQRT((AX1-AX2)**2+(AY1-AY2)**2+(AZ1-AZ2)**2)
      ENDIF
      END
C     ----------------------------------------   
      SUBROUTINE GENERGAPH (ISP,ISN,SBC,ALENGTHPP,ALENGTH,POTOTAL,ALENGTHPY,DATATZSAND,DATATZCLAY,DATAPYCLAY,NTOPOPT,NODETOPMAX,OPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPT
      DIMENSION SBC(37)
      DIMENSION DATATZSAND(6),DATATZCLAY(16)
      DIMENSION DATAPYCLAY(14)
      DIMENSION DATAQZCLAY(10)
      DIMENSION ABC(37)
      COMMON/OVER/POTOL
      COMMON/COUNT/ICOUNT,NCOUNT
      COMMON /TZ/ PREDATATZSAND(6),PREDATATZCLAY(16),PRERESULT(17),NTZZAND
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

C     ---------------------- OUTPUT ------------------------
C     DATATZSAND    =  T-Z CURVE SAND DATA
C     DATATZCLAY    =  T-Z CURVE CLAY DATA
C     DATAPYCLAY    =  P-Y CURVE CLAY DATA
C     -------------------------------------------------------

C     ******** THIS SUBROUTINE IS FOR GENERATE SOIL PROFILE ********
C     THE METHOD OF GENERATE IS RELATE WITH THE SACS PROGRAM 

      ! INTAIL SETTING
      DATATZSAND = 0.D0
      DATATZCLAY = 0.D0
      DATAPYCLAY = 0.D0
      ABC        = 0.D0
      
200   FORMAT("#END")
201   FORMAT("")
      
      ! --- GENERATE CURVE OF SUCTION PILE ----
      NSUC  = SBC(24) ! CONTROL OPTION (0,1)
      
      IF (OPT.EQ."P-Y")THEN
      IF (SBC(1).EQ.1)THEN ! P-Y CUVRE
          ISAND  = SBC(5)   ! OPTION SELECTING OF THE SAND TYPE
          GRMMA  = SBC(6)   ! UNIT WEIGHT OF SOIL
          C      = SBC(7)   ! UNDRAINED SHEAR STRENGTH
          FREE   = SBC(8)   ! FRICTION ANGLE
          AJ     = SBC(9)   ! EMPRIRICAL PARAMETER
          STRAIN = SBC(11)  ! SOIL STRAIN
          XFAC   = SBC(12)  ! X AXIS FACTOR
          YFAC   = SBC(13)  ! Y AXIS FACTOR
          
           !NSUC = 1
          IF (NSUC.EQ.1)THEN  ! FOR SUCTION PILE (Jeanjean 2017)
             C     =  SBC(25) ! C = RATE OF INCREASE OF SHEAR STRENGTH WITH DEPTH 
             ! FINDING SHEAR STRENGTH INTERCEPTE AT SEAFLOOR USING UNDRAINED SHEAR STRENGTH AT THE TOP OF THE LAYER (NODETOPMAX)
             CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,ABC(1:37),NODETOPMAX,ISN,PROPOUT)
             SU1   =  ABC(25)
             RAMDA =  C/(SU1*DPIPE)
             IF (RAMDA.LT.6D0)THEN       ! SI UNIT
                 ETA = 0.25D0+0.05*RAMDA
             ELSEIF (RAMDA.GE.6D0)THEN
                 ETA = 0.55D0
             ENDIF
             ANP = 12D0-(4D0*EXP(-ETA*ALENGTHPY/DPIPE))
             PU  = ANP*C
             ! THE LIMIT OF SUCTION PILE IS 0.9999... BECASE OF ATANH(P/PU)
             ! GMAX = 1000*C 
             Y1  = 0.0D0 
             Y2  = (ATANH(0.23D0)*0.1D0)**2D0
             Y3  = (ATANH(0.33D0)*0.1D0)**2D0
             Y4  = (ATANH(0.50D0)*0.1D0)**2D0
             Y5  = (ATANH(0.72D0)*0.1D0)**2D0
             Y6  = (ATANH(0.998D0)*0.1D0)**2D0
             Y7  = (ATANH(0.998D0)*0.1D0)**2D0
             
             DATAPYCLAY(1)  = 0.000D0           *YFAC
             DATAPYCLAY(2)  = PU*DPIPE*0.23D0   *YFAC 
             DATAPYCLAY(3)  = PU*DPIPE*0.33D0   *YFAC
             DATAPYCLAY(4)  = PU*DPIPE*0.50D0   *YFAC
             DATAPYCLAY(5)  = PU*DPIPE*0.72D0   *YFAC
             DATAPYCLAY(6)  = PU*DPIPE          *YFAC
             DATAPYCLAY(7)  = PU*DPIPE          *YFAC
             DATAPYCLAY(8)  = Y1*DPIPE          *XFAC
             DATAPYCLAY(9)  = Y2*DPIPE          *XFAC
             DATAPYCLAY(10) = Y3*DPIPE          *XFAC
             DATAPYCLAY(11) = Y4*DPIPE          *XFAC
             DATAPYCLAY(12) = Y5*DPIPE          *XFAC
             DATAPYCLAY(13) = Y6*DPIPE          *XFAC
             DATAPYCLAY(14) = Y7*DPIPE          *XFAC
             
          ENDIF
          
          IF (SBC(4).EQ.1) THEN ! CLAY
             IF (NSUC.NE.1) THEN
C           --- API 2A  WSD ---
              
             XR    = 6D0*DPIPE/((GRMMA*DPIPE/C)+AJ)
             !IF (ALENGTHPY.GE.XR)THEN
             IF (ALENGTHPP.GT.XR)THEN     ! 2018
             PU    = (9D0*C)
             ELSEIF (ALENGTHPP.LE.XR)THEN ! 2018
             !ELSEIF (ALENGTHPY.LT.XR)THEN
             PU    = (3D0*C+GRMMA*ALENGTHPY+(AJ*C*ALENGTHPY/DPIPE))
             ENDIF
             YC    = 2.5D0*STRAIN*DPIPE
             
             DATAPYCLAY(1)  = 0.000D0           *YFAC
             DATAPYCLAY(2)  = PU*DPIPE*0.23D0   *YFAC 
             DATAPYCLAY(3)  = PU*DPIPE*0.33D0   *YFAC
             DATAPYCLAY(4)  = PU*DPIPE*0.50D0   *YFAC
             DATAPYCLAY(5)  = PU*DPIPE*0.72D0   *YFAC
             DATAPYCLAY(6)  = PU*DPIPE          *YFAC
             DATAPYCLAY(7)  = PU*DPIPE          *YFAC
             DATAPYCLAY(8)  = 0.00D0            *XFAC
             DATAPYCLAY(9)  = YC*0.1D0          *XFAC
             DATAPYCLAY(10) = YC*0.3D0          *XFAC
             DATAPYCLAY(11) = YC*1.0D0          *XFAC
             DATAPYCLAY(12) = YC*3.0D0          *XFAC
             DATAPYCLAY(13) = YC*8.0D0          *XFAC
             DATAPYCLAY(14) = YC*15D0           *XFAC
             ENDIF    
             
          CALL SOILDATA_CURVE (ISN,ISP,DATAPYCLAY,SBC,"PY_CLAY","WRIT") 
          ! FREE = FRICTION ANGLE
          ELSEIF (SBC(4).EQ.2) THEN ! CLEAN SAND
           FREE = 35D0 ! GIVE BY P-Y SAND DATA
          ELSEIF (SBC(4).EQ.3) THEN ! SILTY SAND
           FREE = 30D0 ! GIVE BY P-Y SAND DATA
          ELSEIF (SBC(4).EQ.4) THEN ! SANDY SILT
           FREE = 25D0 ! GIVE BY P-Y SAND DATA
          ELSEIF (SBC(4).EQ.5) THEN ! GRAVEL
           FREE = 20D0 ! GIVE BY P-Y SAND DATA
          ENDIF
 
          ENDIF
          ! FOR PRINT CLAY DATA GRAPH
                   IF (ICOUNT.EQ.1)THEN
                       WRITE (7102,300)
                       WRITE (7102,301)
                       WRITE (7102,302)
                       WRITE (7102,303)   
                       DO I =1,7
                       WRITE (7102,304) DATAPYCLAY(I+7),DATAPYCLAY(I)
                       ENDDO
                       WRITE (7102,200)
                       WRITE (7102,201)
                       WRITE (7102,201)
300                    FORMAT ('# Graph: "P-Y CLAY"')
301                    FORMAT ("#")      
302                    FORMAT ('# X: "Z" Y: "T"')
303                    FORMAT ("# Units: - -")  
304                    FORMAT (F12.5,2X,F12.5)
                       
                       !WRITE (7120,300)
305                    FORMAT ('P-Y CLAY')
                       !DO I =1,7
                       !WRITE (7120,304) DATAPYCLAY(I+7),DATAPYCLAY(I)
                       !ENDDO
                   ENDIF
      ENDIF
      
      IF (OPT.EQ."T-Z")THEN
c          AK     = SBC(11)
          NUNIT  = SBC(14)
          NSOILY = SBC(15)
          C      = SBC(16)
          GRAMMA = SBC(17)
          DELTA  = SBC(18)
          ANQ    = SBC(19)
          BETA   = SBC(20)
          TFAC   = SBC(21)
          ZFAC   = SBC(22)
          AOPT   = SBC(23)
       POTOL = 0.
       CALL PO_TZ (ISP,ISN,PO)
       IF (NUNIT.EQ.2)THEN ! ENGLISH UNIT
         IF (NSOILY.EQ.2) THEN
          IF (AOPT.EQ.1)THEN ! GRAVEL
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.6108652382) ! 0.8D0
             IF (A.GE.(2.4D0/(12D0*12D0))) THEN
             A = 2400D0/(12D0*12D0)
             ENDIF  
          ELSEIF (AOPT.EQ.2)THEN ! DENSE SAND
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.5235987756D0)
             IF (A.GE.(2000D0/(12D0*12D0))) THEN ! 0.8D0
             A = 2000D0/(12D0*12D0)
          ENDIF
          ELSEIF (AOPT.EQ.3)THEN ! DENSE SAND-SILT
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.436332313D0)
             IF (A.GE.(1700D0/(12D0*12D0))) THEN ! 0.8D0
             A = 1700D0/(12D0*12D0)
             ENDIF  
          ELSEIF (AOPT.EQ.4)THEN ! MEDIUM SAND-SILT
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.3490658504D0)
             IF (A.GE.(1400D0/(12D0*12D0))) THEN ! 0.8D0
             A = 1400D0/(12D0*12D0)
             ENDIF  
          ELSEIF (AOPT.EQ.5)THEN ! MEDIUM SILT
             !PO    = POTOL+GRAMMA*ALENGTHTZ
             A     = PO*TAN(0.2617993878D0)
             IF (A.GE.(1000D0/(12D0*12D0))) THEN ! 0.8D0
             A = 1000D0/(12D0*12D0)
             ENDIF
          ENDIF
       ENDIF
          IF (NSOILY.EQ.1)THEN ! CLAY
C             POTOL  = POTOL+GRAMMA*ALENGTH/2D0
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = C/PO
                IF (A.LE.1.0D0)THEN
                   DE = 0.5D0*(A**(-0.5D0))
                ELSEIF (A.GT.1.0D0)THEN
                   DE = 0.5D0*(A**(-0.25D0))
                ENDIF
          ENDIF   
          
       ELSEIF (NUNIT.EQ.1)THEN ! SI UNIT
      IF (NSOILY.EQ.2) THEN
       IF (AOPT.EQ.1)THEN ! GRAVEL
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.6108652382)
             IF (A.GE.115000D0) THEN
             A = 115000D0
             ENDIF  
          ELSEIF (AOPT.EQ.2)THEN ! DENSE SAND
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.5235987756D0)
             IF (A.GE.96000D0) THEN
             A =96000D0
          ENDIF
          ELSEIF (AOPT.EQ.3)THEN ! DENSE SAND-SILT
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.436332313D0)
             IF (A.GE.81000D0) THEN
             A =81000D0
             ENDIF  
          ELSEIF (AOPT.EQ.4)THEN ! MEDIUM SAND-SILT
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.3490658504D0)
             IF (A.GE.67000D0) THEN
             A =67000D0
             ENDIF  
          ELSEIF (AOPT.EQ.5)THEN ! MEDIUM SILT
             !PO     = POTOL+GRAMMA*ALENGTHTZ
             A      = PO*TAN(0.2617993878D0)
             IF (A.GE.47880D0) THEN
             A =47880D0
             ENDIF  
          ENDIF
      ENDIF
          IF (NSOILY.EQ.1)THEN ! CLAY
C             POTOL  = POTOL+GRAMMA*ALENGTH/2D0
             !PO = POTOL+GRAMMA*ALENGTHTZ
             A  = C/PO
                IF (A.LE.1.0D0)THEN
                   DE = 0.5D0*(A**(-0.5D0))
                ELSEIF (A.GT.1.0D0)THEN
                   DE = 0.5D0*(A**(-0.25D0))
                ENDIF
          ENDIF   
       ENDIF
       
       !POTOL  = POTOL+GRAMMA*ALENGTH
C       IF (PRERESULT(17).EQ.6)THEN ! SAND
           XDISP      = (PREDATATZSAND(5)-PREDATATZSAND(4))*0.002D0
           YDISP      = PREDATATZCLAY(10)
C       ELSEIF (PRERESULT(17).EQ.16)THEN ! CLAY
           
C       ENDIF
       
          IF (NSOILY.EQ.1) THEN
             IF (NTOPOPT.EQ.1)THEN
                 DATATZCLAY(1)  = 0.0D0   *DPIPE         *ZFAC
                 DATATZCLAY(2)  = 0.0016D0*DPIPE         *ZFAC
                 DATATZCLAY(3)  = 0.0031D0*DPIPE         *ZFAC
                 DATATZCLAY(4)  = 0.0057D0*DPIPE         *ZFAC
                 DATATZCLAY(5)  = 0.0080D0*DPIPE         *ZFAC
                 DATATZCLAY(6)  = 0.01D0  *DPIPE         *ZFAC
                 DATATZCLAY(7)  = 0.020D0 *DPIPE         *ZFAC
                 DATATZCLAY(8)  = 0.030D0 *DPIPE         *ZFAC
                 DATATZCLAY(9)  = 0.000D0      *TFAC
                 DATATZCLAY(10) = DE*C*0.303D0 *TFAC  
                 DATATZCLAY(11) = DE*C*0.5D0   *TFAC
                 DATATZCLAY(12) = DE*C*0.754D0 *TFAC
                 DATATZCLAY(13) = DE*C*0.902D0 *TFAC 
                 DATATZCLAY(14) = DE*C         *TFAC
                 DATATZCLAY(15) = 0.7D0*DE*C   *TFAC
                 DATATZCLAY(16) = 0.7D0*DE*C   *TFAC
 
                   IF (ICOUNT.EQ.1)THEN
                       WRITE (7102,100)
                       WRITE (7102,101)
                       WRITE (7102,102)
                       WRITE (7102,103)   
                       DO I =1,8
                       WRITE (7102,104) DATATZCLAY(I),DATATZCLAY(I+8)
                       ENDDO
                       WRITE (7102,200)
                       WRITE (7102,201)
                       WRITE (7102,201)
100                    FORMAT ('# Graph: "T-Z CLAY"')
101                    FORMAT ("#")      
102                    FORMAT ('# X: "Z" Y: "T"')
103                    FORMAT ("# Units: - -")  
104                    FORMAT (F12.5,2X,F12.5)
                   ENDIF
                   
                 CALL SOILDATA_CURVE (ISN,ISP,DATATZCLAY,SBC,"TZ_CLAY","WRIT") 
                 ELSE
             
                 DATATZCLAY(1)  = 0.0D0   *DPIPE       *ZFAC
                 DATATZCLAY(2)  = 0.0016D0*DPIPE       *ZFAC
                 DATATZCLAY(3)  = 0.0031D0*DPIPE       *ZFAC
                 DATATZCLAY(4)  = 0.0057D0*DPIPE       *ZFAC
                 DATATZCLAY(5)  = 0.0080D0*DPIPE       *ZFAC
                 DATATZCLAY(6)  = 0.01D0  *DPIPE       *ZFAC
                 DATATZCLAY(7)  = 0.020D0 *DPIPE       *ZFAC
                 DATATZCLAY(8)  = 0.030D0 *DPIPE       *ZFAC
                 DATATZCLAY(9)  = 0.000D0              *TFAC
                 DATATZCLAY(10) = (DE*C*0.303D0)       *TFAC 
                 DATATZCLAY(11) = (DE*C*0.5D0 )        *TFAC
                 DATATZCLAY(12) = (DE*C*0.754D0)       *TFAC
                 DATATZCLAY(13) = (DE*C*0.902D0)       *TFAC 
                 DATATZCLAY(14) = (DE*C )              *TFAC
                 DATATZCLAY(15) = (0.7D0*DE*C)         *TFAC
                 DATATZCLAY(16) = (0.7D0*DE*C)         *TFAC
    
                   IF (ICOUNT.EQ.1)THEN
                       WRITE (7102,105)
                       WRITE (7102,106)
                       WRITE (7102,107)
                       WRITE (7102,108)
                       DO I =1,8
                       WRITE (7102,109) DATATZCLAY(I),DATATZCLAY(I+8)
                       ENDDO
                       WRITE (7102,200)
                       WRITE (7102,201)
                       WRITE (7102,201)
105                    FORMAT ('# Graph: "T-Z CLAY"')
106                    FORMAT ("#")      
107                    FORMAT ('# X: "Z" Y: "T"')
108                    FORMAT ("# Units: - -")  
109                    FORMAT (F12.5,2X,F12.5)
                   ENDIF
                   
                 CALL SOILDATA_CURVE (ISN,ISP,DATATZCLAY,SBC,"TZ_CLAY","WRIT") 
                 ENDIF
          ENDIF
             
          IF (NSOILY.EQ.2)THEN

          IF (NTZZAND.EQ.1.EQ.NTOPOPT.EQ.0.0D0)THEN ! WITH NTOPOPT.EQ.0.0D0
             DATATZSAND(1) = 0.00D0                     *ZFAC
             DATATZSAND(2) = 0.10D0                     *ZFAC
             DATATZSAND(3) = 0.20D0                     *ZFAC
             DATATZSAND(4) = 0.00D0                     *TFAC
             DATATZSAND(5) = (A)                        *TFAC
             DATATZSAND(6) = (A)                        *TFAC
             
                  IF (ICOUNT.EQ.1)THEN
                       WRITE (7102,110)
                       WRITE (7102,111)
                       WRITE (7102,112)
                       WRITE (7102,113) 
                       DO I = 1,3
                       WRITE (7102,114) DATATZSAND(I),DATATZSAND(I+3)
                       ENDDO
                       WRITE (7102,200)
                       WRITE (7102,201)
                       WRITE (7102,201)
110                    FORMAT ('# Graph: "T-Z Sand"')
111                    FORMAT ("#")      
112                    FORMAT ('# X: "Z" Y: "T"')
113                    FORMAT ("# Units: - -")  
114                    FORMAT (F12.5,2X,F12.5)
                  ENDIF
                  CALL SOILDATA_CURVE (ISN,ISP,DATATZSAND,SBC,"TZ_SAND","WRIT") 

          ELSE
              IF (NUNIT.EQ.2)THEN ! ENGLISH UNIT (in)
              DATATZSAND(1) = 0.00D0                     *ZFAC
              DATATZSAND(2) = 0.10D0                     *ZFAC
              DATATZSAND(3) = 0.20D0                     *ZFAC
              DATATZSAND(4) = 0.00D0                     *TFAC
              DATATZSAND(5) = (A)                        *TFAC
              DATATZSAND(6) = (A)                        *TFAC
              ELSEIF (NUNIT.EQ.1)THEN ! SI UNIT (m)
              DATATZSAND(1) = 0.00D0                     *ZFAC
              DATATZSAND(2) = 0.00254D0                  *ZFAC
              DATATZSAND(3) = 0.00508D0                  *ZFAC
              DATATZSAND(4) = 0.00D0                     *TFAC
              DATATZSAND(5) = (A)                        *TFAC
              DATATZSAND(6) = (A)                        *TFAC
              ENDIF
             
                  IF (ICOUNT.EQ.1)THEN
                       WRITE (7102,120)
                       WRITE (7102,121)
                       WRITE (7102,122)
                       WRITE (7102,123)     
                       DO I = 1,3
                       WRITE (7102,124) DATATZSAND(I),DATATZSAND(I+3)
                       ENDDO
                       WRITE (7102,200)
                       WRITE (7102,201)
                       WRITE (7102,201)
120                    FORMAT ('# Graph: "T-Z Sand"')
121                    FORMAT ("#")      
122                    FORMAT ('# X: "Z" Y: "T"')
123                    FORMAT ("# Units: - -")  
124                    FORMAT (F12.5,2X,F12.5)
                  ENDIF
           CALL SOILDATA_CURVE (ISN,ISP,DATATZSAND,SBC,"TZ_SAND","WRIT") 
           ENDIF
      
          
          ENDIF
             
          ENDIF
      
      IF (OPT.EQ."Q-Z")THEN
          NUNIT  = SBC(14)
          NSOILY = SBC(15)
          C      = SBC(16)
          GRAMMA = SBC(17)
          DELTA  = SBC(18)
          ANQ    = SBC(19)
          BETA   = SBC(20)
          TFAC   = SBC(21)
          ZFAC   = SBC(22)
          AOPT   = SBC(23)
          IF (NSOILY.EQ.2)THEN
          IF (NUNIT.EQ.1)THEN ! ENGLISH UNIT
          IF (AOPT.EQ.1)THEN ! GRAVEL
             Q = POTOTAL*50D0
             IF (Q.GT.250D0*1000D0/(12D0*12D0))THEN ! UNIT TRANSFER KIPS/FT^3 > IB/IN^2
             Q = 250D0*144000D0
             ENDIF
          ELSEIF (AOPT.EQ.2)THEN ! DENSE SAND
             Q = POTOTAL*40D0
             IF (Q.GT.200D0*1000D0/(12D0*12D0))THEN
             Q = 200D0*144000D0
             ENDIF
          ELSEIF (AOPT.EQ.3)THEN ! DENSE SAND-SILT
             Q = POTOTAL*20D0
             IF (Q.GT.170D0*1000D0/(12D0*12D0))THEN
             Q = 170D0*144000D0
             ENDIF
          ELSEIF (AOPT.EQ.4)THEN ! MEDIUM SAND-SILT
             Q = POTOTAL*12D0
             IF (Q.GT.60D0*1000D0/(12D0*12D0))THEN
             Q = 60D0*144000D0
             ENDIF 
          ELSEIF (AOPT.EQ.5)THEN ! MEDIUM SILT
             Q = POTOTAL*8D0
             IF (Q.GT.40D0*1000D0/(12D0*12D0))THEN
             Q = 40D0*144000D0
             ENDIF
          ENDIF
          ENDIF
          ENDIF
          
          IF (NSOILY.EQ.1)THEN ! CLAY
             Q  = 9D0*C
             A  = C/POTOTAL
          ENDIF
          
          ELSEIF (NUNIT.EQ.2)THEN ! SI UNIT
          IF (NSOILY.EQ.2)THEN
          IF (AOPT.EQ.1)THEN ! GRAVEL
             Q = POTOTAL*50D0
             IF (Q.GT.12D0*(10**6D0))THEN
             Q = 12D0*(10**6D0)
             ENDIF
          ELSEIF (AOPT.EQ.2)THEN ! DENSE SAND
             Q = POTOTAL*40D0
             IF (Q.GT.10D0*(10**6D0))THEN
             Q = 10D0*(10**6D0)
             ENDIF
          ELSEIF (AOPT.EQ.3)THEN ! DENSE SAND-SILT
             Q = POTOTAL*20D0
             IF (Q.GT.10D0*(10**6D0))THEN
             Q = 10D0*(10**6D0)
             ENDIF
          ELSEIF (AOPT.EQ.4)THEN ! MEDIUM SAND-SILT
             Q = POTOTAL*12D0
             IF (Q.GT.5D0*(10**6D0))THEN
             Q = 5D0*(10**6D0)
             ENDIF 
          ELSEIF (AOPT.EQ.5)THEN ! MEDIUM SILT
             Q = POTOTAL*8D0
             IF (Q.GT.3D0*(10**6D0))THEN
             Q = 3D0*(10**6D0)
             ENDIF
          ENDIF
          IF (NSOILY.EQ.1)THEN ! CLAY
             Q  = 9D0*C
             A  = C/POTOTAL
          ENDIF
          ENDIF
      ENDIF
      END
C     ---------------------------------------- 
       SUBROUTINE APICURVE (YT,AMR,SBC,DATATZSAND,DATATZCLAY,DATAPYCLAY,STIFF,STIFF_RO,DEPTH,ALENGTHPY,NTSOIL,OPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPT
      DIMENSION SBC(37)
      DIMENSION X(16),Y(16)
      COMMON /MGRAV/ NGRAV 
      COMMON/LARGE_PILE/ NSUC
      DIMENSION YT(3),STIFF(3),AMR(3),STIFF_RO(3)
      DIMENSION DATATZSAND(6),DATATZCLAY(16)
      DIMENSION DATAPYCLAY(14)
      
      COMMON/COUNT/ICOUNT,NCOUNT
      
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
      COMMON / DISSOIL / DATADISP(10000)
      STIFF = 0.
      X     = 0.D0
      Y     = 0.D0
C      DO I = 1,16
C      X(I) = SBC(I+3)    
C      ENDDO
      
C      DO I = 1,16
C      Y(I) = SBC(I+19)     
C      ENDDO
      
      ! READ PROFILE DATA OF SOIL FOR P-Y CURVE
      IF (OPT.EQ."P-Y") THEN
      IF (SBC(1).EQ.1) THEN ! P-Y CURVE CLAY
         IF (SBC(4).EQ.1)THEN
          X(1:7) = DATAPYCLAY(8:14)
          Y(1:7) = DATAPYCLAY(1:7)
          NT     = 7
         ELSEIF (SBC(4).EQ.2)THEN ! P-Y CURVE SAND
          
         ENDIF
      ENDIF
      
      ENDIF
      
      ! READ PROFILE DATA OF SOIL FOR T-Z CURVE
      IF (OPT.EQ."T-Z") THEN
      IF (SBC(2).EQ.1)THEN
        IF (SBC(15).EQ.2)THEN ! SAND
          A      = SBC(23)
          IF (A.EQ.1)THEN
          X(1:3) = DATATZSAND(1:3)
          Y(1:3) = DATATZSAND(4:6)
          NT     = 3
          ELSEIF (A.EQ.2)THEN
          X(1:3) = DATATZSAND(1:3)
          Y(1:3) = DATATZSAND(4:6)
          NT     = 3 
          ELSEIF (A.EQ.3)THEN
          X(1:3) = DATATZSAND(1:3)
          Y(1:3) = DATATZSAND(4:6)
          NT     = 3
          ELSEIF (A.EQ.4)THEN
          X(1:3) = DATATZSAND(1:3)
          Y(1:3) = DATATZSAND(4:6)
          NT     = 3
          ELSEIF (A.EQ.5)THEN
          X(1:3) = DATATZSAND(1:3)
          Y(1:3) = DATATZSAND(4:6)
          NT     = 3
          ENDIF    
        ELSEIF (SBC(15).EQ.1)THEN ! CLAY
        X(1:8) = DATATZCLAY(1:8)
        Y(1:8) = DATATZCLAY(9:16)    
        NT     = 8    
      ENDIF
      ENDIF
      ENDIF
      
      ! CANCEL FOR GENERATE Q-Z CURVE 
      IF (OPT.EQ."Q-Z") THEN
          
      ENDIF
      
      ! ASSUME DISPLACENT FOR INTIAL CONDITION
      IF (ICOUNT.EQ.1) YT(1:3) = (MAXVAL(X)-MINVAL(X))/10000D0
      IF (NSUC.EQ.1)THEN ! SUCTION PILE (JEANJEAN 2017)
         IF (ICOUNT.EQ.1) AMR(1:3) = (MAXVAL(X)-MINVAL(X))/10000D0
      ENDIF
      
      DO I= 1,3
      IF (OPT.EQ."P-Y")THEN ! P-Y CURVE
      ! DISPLACEMENT    
      IF (NGRAV.EQ.1)THEN     ! X-DIRECTION
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = ABS(YT(3))
      ELSEIF (NGRAV.EQ.2)THEN ! Y-DIRECTION
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = ABS(YT(3))
      ELSEIF (NGRAV.EQ.3)THEN ! Z-DIRECTION
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = 0.0D0
      ENDIF
      
!      IF (NSUC.EQ.1)THEN ! SUCTION PILE (JEANJEAN 2017)
!         ! ROTATION 
!         IF (NGRAV.EQ.1)THEN 
!            IF (I.EQ.1) AOU_RO = 0.0D0
!            IF (I.EQ.2) AOU_RO = 0.0D0
!            IF (I.EQ.3) AOU_RO = 0.0D0
!         ELSEIF (NGRAV.EQ.2)THEN
!            IF (I.EQ.1) AOU_RO = 0.0D0
!            IF (I.EQ.2) AOU_RO = 0.0D0
!            IF (I.EQ.3) AOU_RO = 0.0D0
!         ELSEIF (NGRAV.EQ.3)THEN
!            IF (I.EQ.1) AOU_RO = 0.0D0
!            IF (I.EQ.2) AOU_RO = 0.0D0
!            IF (I.EQ.3) AOU_RO = 0.0D0
!         ENDIF
!      ENDIF
      
      ELSEIF (OPT.EQ."T-Z")THEN ! T-Z CURVE
          
      IF (NGRAV.EQ.1)THEN ! DISPLACEMENT
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = ABS(YT(3))
      ENDIF
      
      IF (NSUC.EQ.1)THEN
          IF (NGRAV.EQ.1)THEN ! ROTATION
            IF (I.EQ.1) AOU_RO = ABS(AMR(1))
            IF (I.EQ.2) AOU_RO = 0.0D0
            IF (I.EQ.3) AOU_RO = 0.0D0
         ELSEIF (NGRAV.EQ.2)THEN
            IF (I.EQ.1) AOU_RO = 0.0D0
            IF (I.EQ.2) AOU_RO = ABS(AMR(2))
            IF (I.EQ.3) AOU_RO = 0.0D0
         ELSEIF (NGRAV.EQ.3)THEN
            IF (I.EQ.1) AOU_RO = 0.0D0
            IF (I.EQ.2) AOU_RO = 0.0D0
            IF (I.EQ.3) AOU_RO = ABS(AMR(3))
         ENDIF
      ENDIF
      
      ELSEIF (SBC(3).EQ.1)THEN ! Q-Z CURVE
          
      IF (NGRAV.EQ.1)THEN
         IF (I.EQ.1) AOU = ABS(YT(1))
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.2)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = ABS(YT(2))
         IF (I.EQ.3) AOU = 0.0D0
      ELSEIF (NGRAV.EQ.3)THEN
         IF (I.EQ.1) AOU = 0.0D0
         IF (I.EQ.2) AOU = 0.0D0
         IF (I.EQ.3) AOU = ABS(YT(3))
      ENDIF
      
      ENDIF
  
      IF     (NT.GE.2.AND.AOU.GE.X(1).AND.AOU.LT.X(2))THEN ! DISPLACEMENT
             AIY = Y(2) - Y(1)
             AIX = X(2) - X(1)
             PP  = Y(1) + ABS((AIY/AIX)*(X(1)-AOU))
      ELSEIF (NT.GE.3.AND.AOU.GE.X(2).AND.AOU.LT.X(3))THEN
             AIY = Y(3) - Y(2)
             AIX = X(3) - X(2)
             PP  = Y(2) + ABS((AIY/AIX)*(X(2)-AOU))
      ELSEIF (NT.GE.4.AND.AOU.GE.X(3).AND.AOU.LT.X(4))THEN
             AIY = Y(4) - Y(3)
             AIX = X(4) - X(3)
             PP  = Y(3) + ABS((AIY/AIX)*(X(3)-AOU))
      ELSEIF (NT.GE.5.AND.AOU.GE.X(4).AND.AOU.LT.X(5))THEN
             AIY = Y(5) - Y(4)
             AIX = X(5) - X(4)
             PP  = Y(4) + ABS((AIY/AIX)*(X(4)-AOU))
      ELSEIF (NT.GE.6.AND.AOU.GE.X(5).AND.AOU.LT.X(6))THEN
             AIY = Y(6) - Y(5)
             AIX = X(6) - X(5)
             PP  = Y(5) + ABS((AIY/AIX)*(X(5)-AOU))
      ELSEIF (NT.GE.7.AND.AOU.GE.X(6).AND.AOU.LT.X(7))THEN
             AIY = Y(7) - Y(6)
             AIX = X(7) - X(6)
             PP  = Y(6) + ABS((AIY/AIX)*(X(6)-AOU))
      ELSEIF (NT.GE.8.AND.AOU.GE.X(7).AND.AOU.LT.X(8))THEN
             AIY = Y(8) - Y(7)
             AIX = X(8) - X(7)
             PP  = Y(7) + ABS((AIY/AIX)*(X(7)-AOU))
      ELSEIF (NT.GE.9.AND.AOU.GE.X(8).AND.AOU.LT.X(9))THEN
             AIY = Y(9) - Y(8)
             AIX = X(9) - X(8)
             PP  = Y(8) + ABS((AIY/AIX)*(X(8)-AOU))
      ELSEIF (NT.GE.10.AND.AOU.GE.X(9).AND.AOU.LT.X(10))THEN
             AIY = Y(10) - Y(9)
             AIX = X(10) - X(9)
             PP  = Y(9) + ABS((AIY/AIX)*(X(9)-AOU))
      ELSEIF (NT.GE.11.AND.AOU.GE.X(10).AND.AOU.LT.X(11))THEN
             AIY = Y(11) - Y(10)
             AIX = X(11) - X(10)
             PP  = Y(10) + ABS((AIY/AIX)*(X(10)-AOU))
      ELSEIF (NT.GE.12.AND.AOU.GE.X(11).AND.AOU.LT.X(12))THEN
             AIY = Y(12) - Y(11)
             AIX = X(12) - X(11)
             PP  = Y(11) + ABS((AIY/AIX)*(X(11)-AOU))
      ELSEIF (NT.GE.13.AND.AOU.GE.X(12).AND.AOU.LT.X(13))THEN
             AIY = Y(13) - Y(12)
             AIX = X(13) - X(12)
             PP  = Y(12) + ABS((AIX*(AOU-Y(12))/AIY))
      ELSEIF (NT.GE.14.AND.AOU.GE.X(13).AND.AOU.LT.X(14))THEN
             AIY = Y(14) - Y(13)
             AIX = X(14) - X(13)
             PP  = Y(13) + ABS((AIY/AIX)*(X(13)-AOU))
      ELSEIF (NT.GE.15.AND.AOU.GE.X(14).AND.AOU.LT.X(15))THEN
             AIY = Y(15) - Y(14)
             AIX = X(15) - X(14)
             PP  = Y(14) + ABS((AIY/AIX)*(X(14)-AOU))
      ELSEIF (NT.GE.15.AND.AOU.GE.X(15).AND.AOU.LT.X(16))THEN
             AIY = Y(16) - Y(15)
             AIX = X(16) - X(15)
             PP  = Y(15) + ABS((AIY/AIX)*(X(15)-AOU))
      ELSEIF (NT.GE.16)THEN
             PP  = Y(16)
      ELSE
             PP  =MAXVAL(Y)*1.0
      ENDIF
      
      IF (NSUC.EQ.1)THEN
      IF     (NT.GE.2.AND.AOU_RO.GE.X(1).AND.AOU_RO.LT.X(2))THEN ! ROTATION
             AIY = Y(2) - Y(1)
             AIX = X(2) - X(1)
             PP_RO  = Y(1) + ABS((AIY/AIX)*(X(1)-AOU_RO))
      ELSEIF (NT.GE.3.AND.AOU_RO.GE.X(2).AND.AOU_RO.LT.X(3))THEN
             AIY = Y(3) - Y(2)
             AIX = X(3) - X(2)
             PP_RO  = Y(2) + ABS((AIY/AIX)*(X(2)-AOU_RO))
      ELSEIF (NT.GE.4.AND.AOU_RO.GE.X(3).AND.AOU_RO.LT.X(4))THEN
             AIY = Y(4) - Y(3)
             AIX = X(4) - X(3)
             PP_RO  = Y(3) + ABS((AIY/AIX)*(X(3)-AOU_RO))
      ELSEIF (NT.GE.5.AND.AOU_RO.GE.X(4).AND.AOU_RO.LT.X(5))THEN
             AIY = Y(5) - Y(4)
             AIX = X(5) - X(4)
             PP_RO  = Y(4) + ABS((AIY/AIX)*(X(4)-AOU_RO))
      ELSEIF (NT.GE.6.AND.AOU_RO.GE.X(5).AND.AOU_RO.LT.X(6))THEN
             AIY = Y(6) - Y(5)
             AIX = X(6) - X(5)
             PP_RO  = Y(5) + ABS((AIY/AIX)*(X(5)-AOU_RO))
      ELSEIF (NT.GE.7.AND.AOU_RO.GE.X(6).AND.AOU_RO.LT.X(7))THEN
             AIY = Y(7) - Y(6)
             AIX = X(7) - X(6)
             PP_RO  = Y(6) + ABS((AIY/AIX)*(X(6)-AOU_RO))
      ELSEIF (NT.GE.8.AND.AOU_RO.GE.X(7).AND.AOU_RO.LT.X(8))THEN
             AIY = Y(8) - Y(7)
             AIX = X(8) - X(7)
             PP_RO  = Y(7) + ABS((AIY/AIX)*(X(7)-AOU_RO))
      ELSEIF (NT.GE.9.AND.AOU_RO.GE.X(8).AND.AOU_RO.LT.X(9))THEN
             AIY = Y(9) - Y(8)
             AIX = X(9) - X(8)
             PP_RO  = Y(8) + ABS((AIY/AIX)*(X(8)-AOU_RO))
      ELSEIF (NT.GE.10.AND.AOU_RO.GE.X(9).AND.AOU_RO.LT.X(10))THEN
             AIY = Y(10) - Y(9)
             AIX = X(10) - X(9)
             PP_RO  = Y(9) + ABS((AIY/AIX)*(X(9)-AOU_RO))
      ELSEIF (NT.GE.11.AND.AOU_RO.GE.X(10).AND.AOU_RO.LT.X(11))THEN
             AIY = Y(11) - Y(10)
             AIX = X(11) - X(10)
             PP_RO  = Y(10) + ABS((AIY/AIX)*(X(10)-AOU_RO))
      ELSEIF (NT.GE.12.AND.AOU_RO.GE.X(11).AND.AOU_RO.LT.X(12))THEN
             AIY = Y(12) - Y(11)
             AIX = X(12) - X(11)
             PP_RO  = Y(11) + ABS((AIY/AIX)*(X(11)-AOU_RO))
      ELSEIF (NT.GE.13.AND.AOU_RO.GE.X(12).AND.AOU_RO.LT.X(13))THEN
             AIY = Y(13) - Y(12)
             AIX = X(13) - X(12)
             PP_RO  = Y(12) + ABS((AIX*(AOU_RO-Y(12))/AIY))
      ELSEIF (NT.GE.14.AND.AOU_RO.GE.X(13).AND.AOU_RO.LT.X(14))THEN
             AIY = Y(14) - Y(13)
             AIX = X(14) - X(13)
             PP_RO  = Y(13) + ABS((AIY/AIX)*(X(13)-AOU_RO))
      ELSEIF (NT.GE.15.AND.AOU_RO.GE.X(14).AND.AOU_RO.LT.X(15))THEN
             AIY = Y(15) - Y(14)
             AIX = X(15) - X(14)
             PP_RO  = Y(14) + ABS((AIY/AIX)*(X(14)-AOU_RO))
      ELSEIF (NT.GE.15.AND.AOU_RO.GE.X(15).AND.AOU_RO.LT.X(16))THEN
             AIY = Y(16) - Y(15)
             AIX = X(16) - X(15)
             PP_RO  = Y(15) + ABS((AIY/AIX)*(X(15)-AOU_RO))
      ELSEIF (NT.GE.16)THEN
             PP_RO  = Y(16)
      ELSE
             PP_RO  =MAXVAL(Y)*1.0
      ENDIF
      ENDIF
      
      AROUND = 2D0*3.141592654D0*(DPIPE/2D0)
      IF (OPT.EQ."P-Y")THEN ! P-Y CURVE
      PP    = PP*AROUND
      !IF (NSUC.EQ.1) PP_RO = 0.0D0
      ELSEIF (OPT.EQ."T-Z")THEN ! T-Z CURVE
      PP    = PP*AROUND*DEPTH!*1.2D0
      IF (NSUC.EQ.1) PP_RO = PP_RO*AROUND*DEPTH
      ELSEIF (SBC(3).EQ.1)THEN ! Q-Z CURVE
          
      ENDIF
      ! STIFFNESS
      STIFF(I)    = ABS(PP)/(AOU)
      IF (AOU.EQ.0)    STIFF(I) = 0.0D0
      
      IF (NSUC.EQ.1)THEN ! LARGE PILE
         STIFF_RO(I) = ABS(PP_RO)/(AOU_RO)
         IF (AOU_RO.EQ.0) STIFF_RO(I) = 0.0D0
      ELSEIF (NSUC.NE.0)THEN
         STIFF_RO = 0.D0
      ENDIF
      
!      IF (I.EQ.2) WRITE(7101,'(E12.5)') PP
      ENDDO
      
      END
C     ----------------------------------------
      SUBROUTINE OVERBUND (POTOTAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION SBC(37)
      ! FOR CALCULATE OVERBUNDAL PRESSURE
      POTOTAL = 0.0D0
      CALL READSOILDATA ("CALT",NTSOIL,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
      
      DO 100 NISN = 2,NTSOIL,2
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN,ISN,PROPOUT)
      CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NISN-1,ISP,PROPOUT)

      IF (SBC(2).EQ.1)THEN
          CALL ELEMENTLENGTH (ISP,ISN,ALENGTH,ALENGTHPY,ALENGTHPP,NTSOIL,"QZ")
          GRAMMA  = SBC(17)
          POTOTAL = GRAMMA*ALENGTH + POTOTAL
      ENDIF
      
100   CONTINUE
      END
C     ----------------------------------------
      SUBROUTINE AXIALFORCE (IGM,AXIALFORCEOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      ALLOCATABLE AF1(:)

      AXIALFORCEOUT = 0.0D0
      IF(IGM.LE.0) RETURN
      
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)

C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
          
      ALLOCATE(AF1(N1))

      READ(NFL1,REC=IEL) AF1(1:N1) 
      AXIALFORCEOUT = AF1(1)

      DEALLOCATE(AF1)
        
      END
C     ----------------------------------------
      SUBROUTINE REARRAYNODE (NNODEAR,INDEX)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
CC      COMMON /STOREN/ IELEMENT(100000),N1(100000),N2(100000)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
C     SONGSAK OCT2019
      ALLOCATABLE NELEMENT1(:),NELEMENT2(:),NELEMENT3(:),NELEMENT4(:)
CC      DIMENSION NELEMENT1(500),NELEMENT2(500),NELEMENT3(500),NELEMENT4(500)
      
      DIMENSION NNODEAR(500)

      
C     PREVENT SWAPING NODE ON GID
C     NNODEAR    = TOTAL NUMBER OF NODE IN ONE GROUP
C     INDEX      = TOTAL NUMBER INDEX IN NNODEAR
          INDEX1 = 1
          INDEX2 = 1
          INDEX3 = 1
          INDEX4 = 1
           
CC          DO J = 2,INDEX,2 
CC             ! EACH 2 NODE WAS ONE LAYERS TO CLASSIFY THE ELEMENT
CC             DO 10 K = 1,NSN
CC                 ! FROM THE CONNECTIVITY N1 AND N2
CC                 ! IF NNODEAR IDENTICAL WITH CONNECTIVITY FIND ELEMENT
CC                 IF (NNODEAR(J-1).EQ.N1(K))THEN
CC                     NELEMENT1(INDEX1) = IELEMENT(K)
CC                     INDEX1 = INDEX1 + 1 
CC                 ENDIF
CC                 IF (NNODEAR(J-1).EQ.N2(K))THEN
CC                     NELEMENT2(INDEX2) = IELEMENT(K)
CC                     INDEX2 = INDEX2 + 1 
CC                 ENDIF
CC                 IF (NNODEAR(J).EQ.N1(K))THEN
CC                     NELEMENT3(INDEX3) = IELEMENT(K)
CC                     INDEX3 = INDEX3 + 1 
CC                 ENDIF
CC                 IF (NNODEAR(J).EQ.N2(K))THEN
CC                     NELEMENT4(INDEX4) = IELEMENT(K)
CC                     INDEX4 = INDEX4 + 1 
CC                 ENDIF
CC10            CONTINUE 
CC             ENDDO

C     ----------------------------------------------------------------  
C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT LOOP ... OCT2019  
      
C	CALL TOTAL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NGDAT,NGIDM,1) 
      
      ALLOCATE(NELEMENT1(NGIDM),NELEMENT2(NGIDM),NELEMENT3(NGIDM),NELEMENT4(NGIDM))
      
      DO J = 2,INDEX,2 
      ! EACH 2 NODE WAS ONE LAYERS TO CLASSIFY THE ELEMENT
          
      DO 10 IGM = 1,NGIDM    
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
              CALL CALLNUMNODE_F (IGM,N1,N2) !CALL FRAME CONNECTIVITY NODE
                 ! FROM THE CONNECTIVITY N1 AND N2
                 ! IF NNODEAR IDENTICAL WITH CONNECTIVITY FIND ELEMENT
                 IF (NNODEAR(J-1).EQ.N1)THEN
                     NELEMENT1(INDEX1) = IGM
                     INDEX1 = INDEX1 + 1 
                 ENDIF
                 IF (NNODEAR(J-1).EQ.N2)THEN
                     NELEMENT2(INDEX2) = IGM
                     INDEX2 = INDEX2 + 1 
                 ENDIF
                 IF (NNODEAR(J).EQ.N1)THEN
                     NELEMENT3(INDEX3) = IGM
                     INDEX3 = INDEX3 + 1 
                 ENDIF
                 IF (NNODEAR(J).EQ.N2)THEN
                     NELEMENT4(INDEX4) = IGM
                     INDEX4 = INDEX4 + 1 
                 ENDIF
          ENDIF
10    ENDDO
       
      ENDDO  
C     ----------------------------------------------------------------   
      DEALLOCATE(NELEMENT1,NELEMENT2,NELEMENT3,NELEMENT4)           
             
      END
C     ----------------------------------------
C     ----------------------------------------
C     ----------------------------------------
          SUBROUTINE INTIALSOIL 
          IMPLICIT REAL*8 (A-H,O-Z)
          IMPLICIT INTEGER*4 (I-N)
          COMMON / ITERATION / APSTIFF(3,50000),ATSTIFF(3,50000),AQSTIFF(3,50000),AUSTIFF(3,50000),ATSTIFF_RO(3,50000)
          COMMON / TZ/ PREDATATZSAND(6),PREDATATZCLAY(16)
          COMMON / SOILSTIFFNESS / SSTIFF(6,1000)
          APSTIFF = 0.
          ATSTIFF = 0.
          AQSTIFF = 0.
          AUSTIFF = 0.
          ATSTIFF_RO = 0.
          PREDATATZSAND = 0.
          PREDATATZCLAY = 0.
          SSTIFF = 0.
          END
C     ----------------------------------------
          SUBROUTINE SOILDATA_CURVE (ISN,ISP,DATA_CURVE,SBC,OPTDATA,OPT) 
          IMPLICIT REAL*8 (A-H,O-Z)
          IMPLICIT INTEGER*4 (I-N)
          CHARACTER*200 LHDSTD,LHDCOM,NAMEL
          CHARACTER*7 OPTDATA
          CHARACTER*4 OPT
          CHARACTER*10 SOIL_TYPE
          CHARACTER*16 SOIL_TYPE_TZ
          !COMMON /DESIGN/ AX(100000),AY(100000),AZ(100000),NN,NND(100000),NCOM
          COMMON /SOILOPT/ NSOIL,NNSO
          COMMON /NSOILT/ NNODEA(10000),NPROP(10000)
          COMMON /MGRAV/ NGRAV
          COMMON /SOIL_OUTPUT/ CURVE_DATA_PY(70,500),CURVE_DATA_TZ(70,500),CURVE_DATA_QZ(70,500),INDEX_STADA
          COMMON /STCAS/ ILC
          COMMON /LHEADER/ LHDSTD(500),LHDCOM(500)
          DIMENSION DATA_CURVE(30)
          DIMENSION NDATA_NODE(NNSO),NDATA_LOC(NNSO),NDATA_LOC_RE(NNSO),COOR_DATA(NNSO)
          DIMENSION SBC(37)
          COOR_DATA  = 0.
          NDATA_NODE = 0.
          NDATA_LOC  = 0.
          
          IF (OPT.EQ."WRIT") THEN ! STORAGE SOIL DATA
              
              ! DATA LOCATION
              DO JJ = 1,NNSO
              !IF (NGRAV.EQ.1)THEN
              !COOR_DATA(JJ) = AX(NNODEA(JJ))
              !ELSEIF (NGRAV.EQ.2)THEN
              !COOR_DATA(JJ) = AY(NNODEA(JJ))
              !ELSEIF (NGRAV.EQ.3)THEN
              !COOR_DATA(JJ) = AZ(NNODEA(JJ))
              !ENDIF
              IF (NGRAV.EQ.1)THEN
              CALL RELFILL('@XYZ',COOR_DATA(JJ),1,NNODEA(JJ),0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.2)THEN
	        CALL RELFILL('@XYZ',COOR_DATA(JJ),2,NNODEA(JJ),0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.3)THEN
              CALL RELFILL('@XYZ',COOR_DATA(JJ),3,NNODEA(JJ),0)  !GETTING NODAL COORDINATE	
              ENDIF
              ENDDO
                INDEX = 1
                NDATA_NODE(1:NNSO) = NNODEA(1:NNSO)
                DO 100 I =1,NNSO,2
                    IF (I.EQ.39)THEN
                        AAAA  = 1
                     ENDIF
                  IF (NNODEA(I).NE.0)THEN
                    LOCA_MAX   = MAXLOC(COOR_DATA,DIM=1)
                    VALUE_MAX  = MAXVAL(COOR_DATA)
                    VALUE_MIN  = MINVAL(COOR_DATA)
                    NDATA_MAX  = NDATA_NODE(LOCA_MAX)
                    NDATA_NODE(LOCA_MAX) = 0.0D0
                    COOR_DATA(LOCA_MAX)  = VALUE_MIN - 1D0
                    !COOR_DATA(II)  = 0.0D0
                    DO II = 1,NNSO
                        CALL RELFILL('@XYZ',AX_MAX,1,NDATA_MAX,0)  !GETTING NODAL COORDINATE
	                  CALL RELFILL('@XYZ',AY_MAX,2,NDATA_MAX,0)  !GETTING NODAL COORDINATE
                        CALL RELFILL('@XYZ',AZ_MAX,3,NDATA_MAX,0)  !GETTING NODAL COORDINATE	
                        CALL RELFILL('@XYZ',AOUT_X,1,NDATA_NODE(II),0)  !GETTING NODAL COORDINATE
	                  CALL RELFILL('@XYZ',AOUT_Y,2,NDATA_NODE(II),0)  !GETTING NODAL COORDINATE
                        CALL RELFILL('@XYZ',AOUT_Z,3,NDATA_NODE(II),0)  !GETTING NODAL COORDINATE	
                     IF (NDATA_MAX.EQ.NDATA_NODE(II)) THEN
                        NDATA_NODE(II) = 0.0D0
                        COOR_DATA(II)  = VALUE_MIN - 1D0
                        !COOR_DATA(II)  = 0.0D0
                     ELSEIF (NDATA_MAX.NE.NDATA_NODE(II))THEN
                         IF(NGRAV.EQ.1)THEN
                             !IF (AX(NDATA_MAX).EQ.AX(NDATA_NODE(II)))THEN
                             IF (AX_MAX.EQ.AOUT_X)THEN
                                NDATA_NODE(II) = 0.0D0
                                COOR_DATA(II)  = VALUE_MIN - 1D0
                             ENDIF
                         ELSEIF (NGRAV.EQ.2)THEN
                             IF (AY_MAX.EQ.AOUT_Y)THEN
                                NDATA_NODE(II) = 0.0D0
                                COOR_DATA(II)  = VALUE_MIN - 1D0
                             ENDIF
                          ELSEIF (NGRAV.EQ.3)THEN
                            IF (AZ_MAX.EQ.AOUT_Z)THEN
                                NDATA_NODE(II) = 0.0D0
                                COOR_DATA(II)  = VALUE_MIN - 1D0
                             ENDIF
                         ENDIF
                     ENDIF
                    ENDDO
                    ! STRATA LOCATION
                    NDATA_LOC(INDEX) = NDATA_MAX
                    INDEX = INDEX + 1
                    INDEX_STADA = INDEX - 1 ! INDEX STADA LOOP 
                    ENDIF
 100            CONTINUE
                
                !1INDEX_STADA = 0.0D0
                !1DO JJ = 1,NNSO
                !1IF (NDATA_LOC(JJ).NE.0) INDEX_STADA = INDEX_STADA + 1
                !1ENDDO
                
                NDATA_LOC_RE = 0.0D0
                INDEX = 1
                DO IJK = 1,NNSO
                IF (NDATA_LOC(IJK).NE.0)THEN
                    NDATA_LOC_RE(INDEX) = NDATA_LOC(IJK)
                    INDEX = INDEX + 1
                ENDIF
                ENDDO
                  
              IF (OPTDATA.EQ."PY_CLAY") THEN ! P-Y-CLAY
                  DO 200 I = 1,INDEX
                      IF (NGRAV.EQ.1)THEN
                      !A_LOC = AX(NDATA_LOC_RE(I))
                      !A_ISP = AX(ISP)
                      CALL RELFILL('@XYZ',A_LOC,1,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,1,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.2)THEN
                      !A_LOC = AY(NDATA_LOC_RE(I))
                      !A_ISP = AY(ISP)
                      CALL RELFILL('@XYZ',A_LOC,2,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,2,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.3)THEN
                      !A_LOC = AZ(NDATA_LOC_RE(I))
                      !A_ISP = AZ(ISP)
                      CALL RELFILL('@XYZ',A_LOC,3,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,3,ISP,0)              !GETTING NODAL COORDINATE
                      ENDIF
                      IF (A_LOC.EQ.A_ISP) THEN
                      NLOCATION = I
                      EXIT
                      ENDIF
200               CONTINUE
                  CURVE_DATA_PY(1,NLOCATION)    = ISP
                  CURVE_DATA_PY(2,NLOCATION)    = ISN
                  CURVE_DATA_PY(3,NLOCATION)    = 14D0
                  DO I = 1,14
                  CURVE_DATA_PY(I+3,NLOCATION) = DATA_CURVE(I)
                  ENDDO
              ELSEIF (OPTDATA.EQ."PY_SAND") THEN ! P-Y-SAND
                  DO 201 I = 1,INDEX
                      IF (NGRAV.EQ.1)THEN
                      !A_LOC = AX(NDATA_LOC_RE(I))
                      !A_ISP = AX(ISP)
                      CALL RELFILL('@XYZ',A_LOC,1,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,1,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.2)THEN
                      !A_LOC = AY(NDATA_LOC_RE(I))
                      !A_ISP = AY(ISP)
                      CALL RELFILL('@XYZ',A_LOC,2,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,2,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.3)THEN
                      !A_LOC = AZ(NDATA_LOC_RE(I))
                      !A_ISP = AZ(ISP)
                      CALL RELFILL('@XYZ',A_LOC,3,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,3,ISP,0)              !GETTING NODAL COORDINATE
                      ENDIF
                      IF (A_LOC.EQ.A_ISP) THEN
                      NLOCATION = I
                      EXIT
                      ENDIF
201               CONTINUE
                  CURVE_DATA_PY(1,NLOCATION)    = ISP
                  CURVE_DATA_PY(2,NLOCATION)    = ISN
                  CURVE_DATA_PY(3,NLOCATION)    = 30D0
                  DO I = 1,30
                  CURVE_DATA_PY(I+3,NLOCATION) = DATA_CURVE(I)
                  ENDDO 
              ELSEIF (OPTDATA.EQ."TZ_CLAY") THEN ! T-Z CLAY
                  DO 202 I = 1,INDEX
                      IF (NGRAV.EQ.1)THEN
                      !A_LOC = AX(NDATA_LOC_RE(I))
                      !A_ISP = AX(ISP)
                      CALL RELFILL('@XYZ',A_LOC,1,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,1,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.2)THEN
                      !A_LOC = AY(NDATA_LOC_RE(I))
                      !A_ISP = AY(ISP)
                      CALL RELFILL('@XYZ',A_LOC,2,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,2,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.3)THEN
                      !A_LOC = AZ(NDATA_LOC_RE(I))
                      !A_ISP = AZ(ISP)
                      CALL RELFILL('@XYZ',A_LOC,3,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,3,ISP,0)              !GETTING NODAL COORDINATE
                      ENDIF
                      IF (A_LOC.EQ.A_ISP) THEN
                      NLOCATION = I
                      EXIT
                      ENDIF
202               CONTINUE
                  CURVE_DATA_TZ(1,NLOCATION)    = ISP
                  CURVE_DATA_TZ(2,NLOCATION)    = ISN
                  CURVE_DATA_TZ(3,NLOCATION)    = 16D0
                  DO I = 1,16
                  CURVE_DATA_TZ(I+3,NLOCATION) = DATA_CURVE(I)
                  ENDDO  
              ELSEIF (OPTDATA.EQ."TZ_SAND") THEN ! T-Z FOR SAND
                  DO 203 I = 1,INDEX
                      IF (NGRAV.EQ.1)THEN
                      !A_LOC = AX(NDATA_LOC(I))
                      !A_ISP = AX(ISP)
                      CALL RELFILL('@XYZ',A_LOC,1,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,1,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.2)THEN
                      !A_LOC = AY(NDATA_LOC(I))
                      !A_ISP = AY(ISP)
                      CALL RELFILL('@XYZ',A_LOC,2,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,2,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.3)THEN
                      !A_LOC = AZ(NDATA_LOC(I))
                      !A_ISP = AZ(ISP)
                      CALL RELFILL('@XYZ',A_LOC,3,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,3,ISP,0)              !GETTING NODAL COORDINATE
                      ENDIF
                      IF (A_LOC.EQ.A_ISP) THEN
                      NLOCATION = I
                      EXIT
                      ENDIF
203               CONTINUE
                  CURVE_DATA_TZ(1,NLOCATION)    = ISP
                  CURVE_DATA_TZ(2,NLOCATION)    = ISN
                  CURVE_DATA_TZ(3,NLOCATION)    = 6D0
                  DO I = 1,6
                  CURVE_DATA_TZ(I+3,NLOCATION) = DATA_CURVE(I)
                  ENDDO
              ELSEIF (OPTDATA.EQ."USERDEF") THEN ! USER DEFIND FUNCTION
                  DO 204 I = 1,INDEX
                      IF (NGRAV.EQ.1)THEN
                      !A_LOC = AX(NDATA_LOC(I))
                      !A_ISP = AX(ISP)
                      CALL RELFILL('@XYZ',A_LOC,1,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,1,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.2)THEN
                      !A_LOC = AY(NDATA_LOC(I))
                      !A_ISP = AY(ISP)
                      CALL RELFILL('@XYZ',A_LOC,2,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,2,ISP,0)              !GETTING NODAL COORDINATE
                      ELSEIF (NGRAV.EQ.3)THEN
                      !A_LOC = AZ(NDATA_LOC(I))
                      !A_ISP = AZ(ISP)
                      CALL RELFILL('@XYZ',A_LOC,3,NDATA_LOC_RE(I),0)  !GETTING NODAL COORDINATE
                      CALL RELFILL('@XYZ',A_ISP,3,ISP,0)              !GETTING NODAL COORDINATE
                      ENDIF
                      IF (A_LOC.EQ.A_ISP) THEN
                      NLOCATION = I
                      EXIT
                      ENDIF
204               CONTINUE
                  IF (SBC(1).EQ.1)THEN
                    CURVE_DATA_PY(1,NLOCATION)    = ISP
                    CURVE_DATA_PY(2,NLOCATION)    = ISN
                    CURVE_DATA_PY(3,NLOCATION)    = 32D0   
                    DO I = 1,32
                    CURVE_DATA_PY(I,NLOCATION)    = SBC(I+3)    
                    ENDDO 
                  ELSEIF (SBC(2).EQ.2)THEN
                    CURVE_DATA_TZ(1,NLOCATION)    = ISP
                    CURVE_DATA_TZ(2,NLOCATION)    = ISN
                    CURVE_DATA_TZ(3,NLOCATION)    = 32D0   
                    DO I = 1,32
                    CURVE_DATA_TZ(I,NLOCATION)    = SBC(I+3)    
                    ENDDO   
                  ELSEIF (SBC(3).EQ.3)THEN
                    CURVE_DATA_QZ(1,NLOCATION)    = ISP
                    CURVE_DATA_QZ(2,NLOCATION)    = ISN
                    CURVE_DATA_QZ(3,NLOCATION)    = 32D0   
                    DO I = 1,32
                    CURVE_DATA_QZ(I,NLOCATION)    = SBC(I+3)    
                    ENDDO 
                  ENDIF
              ENDIF
              
          ELSEIF (OPT.EQ."PRIN") THEN  ! PRINT VALUE
          ! PRINT HEAD
          LENGTH = LEN_TRIM(LHDSTD(ILC))    
          NAMEL  = LHDSTD(ILC)
          WRITE (7120,1500) NAMEL(1:LENGTH)
1500      FORMAT ("NAME OF LOAD CASE:",2X,A)   
          ! SOIL PROPERTIES      
          WRITE (7120,2000) 
          WRITE (7120,3500)
          ! INTIAL VALUE FOR UNIT
          DO 3900 I = 1,100
               NODEINPUT = CURVE_DATA_PY(2,I)
               CALL READSOILDATA ("CALO",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               NUNIT  = SBC(14)
               IF (NUNIT.NE.0) EXIT
3900      CONTINUE 
          IF (NUNIT.EQ.1) WRITE (7120,3401) ! ENGLISH UNIT HEADER
          IF (NUNIT.EQ.2) WRITE (7120,3402) ! SI UNIT HEADER
           IF (NUNIT.NE.1.AND.NUNIT.NE.2) THEN
             WRITE (7120,3401)
             IF (CURVE_DATA_PY(1,I).EQ.0.AND.CURVE_DATA_PY(2,I).EQ.0) THEN
             WRITE (7120,1005)
             WRITE (7120,1003)
             WRITE (7120,1003)
             ENDIF
          ENDIF
            INDEX_WRITE = 1D0
           DO 3000 I = 1,INDEX_STADA
               ISP = CURVE_DATA_PY(1,I)
               ISN = CURVE_DATA_PY(2,I)
               !CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               CALL READSOILPROP (ISN,ISP,SBC)
               IF (NODEINPUT.EQ.0) GOTO 3000
                    ISAND  = SBC(5)   ! OPTION SELECTING OF THE SAND TYPE
                    GRMMA  = SBC(6)   ! UNIT WEIGHT OF SOIL
                    CA     = SBC(7)   ! UNDRAINED SHEAR STRENGTH
                    FREE   = SBC(8)   ! FRICTION ANGLE
                    AJ     = SBC(9)   ! EMPRIRICAL PARAMETER
                    AK     = SBC(10)  ! SUBGRADE REACTION
                    STRAIN = SBC(11)  ! SOIL STRAIN
                    XFAC   = SBC(12)  ! X AXIS FACTOR
                    YFAC   = SBC(13)  ! Y AXIS FACTOR
                    
                    IF (SBC(4).EQ.1)THEN
                    SOIL_TYPE = "   CLAY"
                    AK        = 0.0D0
                    FREE      = 0.0D0
                    ENDIF
                    IF (SBC(4).NE.1)THEN
                    IF (ISAND.EQ.1) SOIL_TYPE = "   SILT"
                    IF (ISAND.EQ.2) SOIL_TYPE = "SANDY SILT"
                    IF (ISAND.EQ.3) SOIL_TYPE = "SILTY SAND"
                    IF (ISAND.EQ.4) SOIL_TYPE = "CLEAN SAND"  
                    IF (ISAND.EQ.5) SOIL_TYPE = "   GRAVEL"
                    AJ     = 0.0D0
                    STRAIN = 0.0D0
                    CA     = 0.0D0
                    ENDIF
                        
                    !IF (NGRAV.EQ.1)THEN
                    !!AC_OUT_MAX = AX(CURVE_DATA_PY(1,I))   
                    !!AC_OUT_MIN = AX(CURVE_DATA_PY(2,I))  
                    !CALL RELFILL('@XYZ',AC_OUT_MAX,1,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                    !CALL RELFILL('@XYZ',AC_OUT_MIN,1,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                    !ELSEIF (NGRAV.EQ.2)THEN
                    !!AC_OUT_MAX = AY(CURVE_DATA_PY(1,I))   
                    !!AC_OUT_MIN = AY(CURVE_DATA_PY(2,I))   
                    !CALL RELFILL('@XYZ',AC_OUT_MAX,2,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                    !CALL RELFILL('@XYZ',AC_OUT_MIN,2,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                    !ELSEIF (NGRAV.EQ.3)THEN
                    !!AC_OUT_MAX = AZ(CURVE_DATA_PY(1,I))   
                    !!AC_OUT_MIN = AZ(CURVE_DATA_PY(2,I)) 
                    !CALL RELFILL('@XYZ',AC_OUT_MAX,3,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                    !CALL RELFILL('@XYZ',AC_OUT_MIN,3,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                    !ENDIF
                    IF (NGRAV.EQ.1)THEN
                    !AC_OUT_MAX = AX(CURVE_DATA_PY(1,I))   
                    !AC_OUT_MIN = AX(CURVE_DATA_PY(2,I)) 
                    NNODEMAX = CURVE_DATA_PY(1,I)
                    NNODEMIN = CURVE_DATA_PY(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.2)THEN
                    !AC_OUT_MAX = AY(CURVE_DATA_PY(1,I))   
                    !AC_OUT_MIN = AY(CURVE_DATA_PY(2,I))   
                    NNODEMAX = CURVE_DATA_PY(1,I)
                    NNODEMIN = CURVE_DATA_PY(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.3)THEN
                    !AC_OUT_MAX = AZ(CURVE_DATA_PY(1,I))   
                    !AC_OUT_MIN = AZ(CURVE_DATA_PY(2,I))  
                    NNODEMAX = CURVE_DATA_PY(1,I)
                    NNODEMIN = CURVE_DATA_PY(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ENDIF
                    AK = AK/1000D0
                    IF(AC_OUT_MIN.EQ.AC_OUT_MAX) GOTO 3000
                    IF(AC_OUT_MIN.EQ.0.0D0.AND.AC_OUT_MAX.EQ.0.0D0) GOTO 3000
               WRITE (7120,3501) INDEX_WRITE,AC_OUT_MIN,AC_OUT_MAX,SOIL_TYPE,CA,GRMMA,STRAIN,AJ,AK,FREE
               INDEX_WRITE = INDEX_WRITE + 1D0
3000       CONTINUE
          WRITE (7120,3502)
2000      FORMAT ("SOIL PROPERTIES FOR P-Y CURVE")    
3500      FORMAT ("STRATA    FROM DEPTH     TO DEPTH   SOIL-TYPE           C         UNIT-WEIGHT   SOIL-STRAIN        J
     1    MODULUS-SUBGRADE  FRICTION-ANGLE")       
3401  FORMAT ("             (m)          (m)                          N/m^2         N/m^3
     1                                    kN/m^3        DEGREE")      
3402  FORMAT ("             (in)          (in)                       lb/in^2       lb/in^3
     1                                 kips/in^3        DEGREE")  
3501  FORMAT (I5,2X,F12.3,2X,F12.3,4X,A10,2X,F12.5,2X,F12.5,2X,F12.5,2X,F12.5,2X,F12.5,2X,F12.5)   
3502  FORMAT ("") 
      
          ! P-Y CURVE     
          ! INTIAL VALUE FOR UNIT
          DO 4000 I = 1,INDEX_STADA
               ISP = CURVE_DATA_PY(1,I)
               ISN = CURVE_DATA_PY(2,I)
               !CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               CALL READSOILPROP (ISN,ISP,SBC)
               NUNIT  = SBC(14)
               IF (NUNIT.NE.0) EXIT
4000      CONTINUE
          
          WRITE (7120,1000)
          IF (NUNIT.EQ.1) WRITE (7120,1101) ! ENGLISH UNIT HEADER
          IF (NUNIT.EQ.2) WRITE (7120,1102) ! SI UNIT HEADER
          IF (NUNIT.NE.1.AND.NUNIT.NE.2) THEN
             WRITE (7120,1101)
             IF (CURVE_DATA_PY(1,I).EQ.0.AND.CURVE_DATA_PY(2,I).EQ.0) THEN
             WRITE (7120,1005)
             WRITE (7120,1003)
             WRITE (7120,1003)
             WRITE (7120,1003)
             ENDIF
          ENDIF
          INDEX_WRITE = 1D0
          DO 5000 I = 1,INDEX_STADA
               IF (CURVE_DATA_PY(1,I).EQ.0.AND.CURVE_DATA_PY(2,I).EQ.0) GOTO 5000
               ISP = CURVE_DATA_PY(1,I)
               ISN = CURVE_DATA_PY(2,I)
               !CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               CALL READSOILPROP (ISN,ISP,SBC)
               IF (NCODE.EQ.1)THEN  ! API AND JEAN JEAN CURVE
               XFAC   = SBC(12)  ! X AXIS FACTOR
               YFAC   = SBC(13)  ! Y AXIS FACTOR
               IF (CURVE_DATA_PY(3,I).EQ.14)THEN ! P-Y CLAY
                  !IF (NGRAV.EQ.1)THEN
                  !!AC_OUT_MAX = AX(CURVE_DATA_PY(1,I))   
                  !!AC_OUT_MIN = AX(CURVE_DATA_PY(2,I))  
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,1,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,1,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.2)THEN
                  !!AC_OUT_MAX = AY(CURVE_DATA_PY(1,I))   
                  !!AC_OUT_MIN = AY(CURVE_DATA_PY(2,I))  
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,2,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,2,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.3)THEN
                  !!AC_OUT_MAX = AZ(CURVE_DATA_PY(1,I))   
                  !!AC_OUT_MIN = AZ(CURVE_DATA_PY(2,I)) 
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,3,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,3,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                  !ENDIF
                  IF (NGRAV.EQ.1)THEN
                  !AC_OUT_MAX = AX(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AX(CURVE_DATA_PY(2,I)) 
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ELSEIF (NGRAV.EQ.2)THEN
                  !AC_OUT_MAX = AY(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AY(CURVE_DATA_PY(2,I))   
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ELSEIF (NGRAV.EQ.3)THEN
                  !AC_OUT_MAX = AZ(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AZ(CURVE_DATA_PY(2,I))  
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ENDIF
               WRITE (7120,1001) INDEX_WRITE,YFAC,AC_OUT_MIN,AC_OUT_MAX,CURVE_DATA_PY(4,I),CURVE_DATA_PY(11,I)
     1                           ,CURVE_DATA_PY(5,I),CURVE_DATA_PY(12,I),CURVE_DATA_PY(6,I),CURVE_DATA_PY(13,I)
     1                           ,CURVE_DATA_PY(7,I),CURVE_DATA_PY(14,I),CURVE_DATA_PY(8,I),CURVE_DATA_PY(15,I)
               WRITE (7120,1002) CURVE_DATA_PY(9,I),CURVE_DATA_PY(16,I),CURVE_DATA_PY(10,I),CURVE_DATA_PY(17,I)
               WRITE (7120,1003) ! SPACE
               INDEX_WRITE = INDEX_WRITE + 1D0
                  ELSEIF (CURVE_DATA_PY(3,I).EQ.30)THEN ! P-Y SAND
                  !IF (NGRAV.EQ.1)THEN
                  !!AC_OUT_MAX = AX(CURVE_DATA_PY(1,I))   
                  !!AC_OUT_MIN = AX(CURVE_DATA_PY(2,I))  
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,1,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,1,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.2)THEN
                  !!AC_OUT_MAX = AY(CURVE_DATA_PY(1,I))   
                  !!AC_OUT_MIN = AY(CURVE_DATA_PY(2,I))  
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,2,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,2,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.3)THEN
                  !!AC_OUT_MAX = AZ(CURVE_DATA_PY(1,I))   
                  !!AC_OUT_MIN = AZ(CURVE_DATA_PY(2,I))  
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,3,CURVE_DATA_PY(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,3,CURVE_DATA_PY(2,I),0)  !GETTING NODAL COORDINATE
                  !ENDIF
                  IF (NGRAV.EQ.1)THEN
                  !AC_OUT_MAX = AX(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AX(CURVE_DATA_PY(2,I)) 
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ELSEIF (NGRAV.EQ.2)THEN
                  !AC_OUT_MAX = AY(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AY(CURVE_DATA_PY(2,I))   
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ELSEIF (NGRAV.EQ.3)THEN
                  !AC_OUT_MAX = AZ(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AZ(CURVE_DATA_PY(2,I))  
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ENDIF
                  
               WRITE (7120,1001) INDEX_WRITE,YFAC,AC_OUT_MIN,AC_OUT_MAX,CURVE_DATA_PY(4,I),CURVE_DATA_PY(19,I)
     1                           ,CURVE_DATA_PY(5,I),CURVE_DATA_PY(20,I),CURVE_DATA_PY(6,I),CURVE_DATA_PY(21,I)
     1                           ,CURVE_DATA_PY(7,I),CURVE_DATA_PY(22,I),CURVE_DATA_PY(8,I),CURVE_DATA_PY(23,I)
               WRITE (7120,1004)  CURVE_DATA_PY(9,I),CURVE_DATA_PY(24,I)
     1                           ,CURVE_DATA_PY(10,I),CURVE_DATA_PY(25,I),CURVE_DATA_PY(11,I),CURVE_DATA_PY(26,I)
     1                           ,CURVE_DATA_PY(12,I),CURVE_DATA_PY(27,I),CURVE_DATA_PY(13,I),CURVE_DATA_PY(28,I)
               WRITE (7120,1004)  CURVE_DATA_PY(14,I),CURVE_DATA_PY(29,I)
     1                           ,CURVE_DATA_PY(15,I),CURVE_DATA_PY(30,I),CURVE_DATA_PY(16,I),CURVE_DATA_PY(31,I)
     1                           ,CURVE_DATA_PY(17,I),CURVE_DATA_PY(32,I),CURVE_DATA_PY(18,I),CURVE_DATA_PY(33,I)
               WRITE (7120,1003) ! SPACE
               INDEX_WRITE = INDEX_WRITE + 1D0
               ENDIF
               ELSEIF (NCODE.EQ.3) THEN ! USER DEFINED FUNCTION
               XFAC   = SBC(36)  ! X AXIS FACTOR
               YFAC   = SBC(37)  ! Y AXIS FACTOR
                  IF (NGRAV.EQ.1)THEN
                  !AC_OUT_MAX = AX(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AX(CURVE_DATA_PY(2,I)) 
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ELSEIF (NGRAV.EQ.2)THEN
                  !AC_OUT_MAX = AY(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AY(CURVE_DATA_PY(2,I))   
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ELSEIF (NGRAV.EQ.3)THEN
                  !AC_OUT_MAX = AZ(CURVE_DATA_PY(1,I))   
                  !AC_OUT_MIN = AZ(CURVE_DATA_PY(2,I))  
                  NNODEMAX = CURVE_DATA_PY(1,I)
                  NNODEMIN = CURVE_DATA_PY(2,I)
                  CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                  CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ENDIF
                  
               WRITE (7120,1001) INDEX_WRITE,YFAC,AC_OUT_MIN,AC_OUT_MAX,CURVE_DATA_PY(4,I),CURVE_DATA_PY(20,I)
     1                           ,CURVE_DATA_PY(5,I),CURVE_DATA_PY(21,I),CURVE_DATA_PY(6,I),CURVE_DATA_PY(22,I)
     1                           ,CURVE_DATA_PY(7,I),CURVE_DATA_PY(23,I),CURVE_DATA_PY(8,I),CURVE_DATA_PY(24,I)
               WRITE (7120,1004)  CURVE_DATA_PY(9,I),CURVE_DATA_PY(25,I)
     1                           ,CURVE_DATA_PY(10,I),CURVE_DATA_PY(26,I),CURVE_DATA_PY(11,I),CURVE_DATA_PY(27,I)
     1                           ,CURVE_DATA_PY(12,I),CURVE_DATA_PY(28,I),CURVE_DATA_PY(13,I),CURVE_DATA_PY(29,I)
               WRITE (7120,1004)  CURVE_DATA_PY(14,I),CURVE_DATA_PY(30,I)
     1                           ,CURVE_DATA_PY(15,I),CURVE_DATA_PY(31,I),CURVE_DATA_PY(16,I),CURVE_DATA_PY(32,I)
     1                           ,CURVE_DATA_PY(17,I),CURVE_DATA_PY(33,I),CURVE_DATA_PY(18,I),CURVE_DATA_PY(34,I)
               WRITE (7120,1006) CURVE_DATA_PY(19,I),CURVE_DATA_PY(35,I)
               WRITE (7120,1003) ! SPACE
               ENDIF
5000      CONTINUE

1005      FORMAT ("**** NO DATA EXIST ****")
1000      FORMAT ("STRATA FACTOR   FROM DEPTH    TO DEPTH         P             Y             P             Y             P
     1             Y             P             Y             P             Y")
1001      FORMAT (I5,2X,F5.3,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X
     1           ,F12.3,2X,F12.3)
1002      FORMAT (40X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X)
1003      FORMAT ("")  
1004      FORMAT (40X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X)
1101      FORMAT ("                  (m)            (m)           (N/m)        (m)           (N/m)         (m)         (N/m)
     1         (m)         (N/m)           (m)          (N/m)          (m)")
1102      FORMAT ("                  (in)           (in)        (lb/in)        (in)        (lb/in)         (in)        (lb/in)
     1         (in)        (lb/in)         (in)        (lb/in)         (in)")

1010      FORMAT ("STRATA FACTOR   FROM DEPTH    TO DEPTH         T             Z             T             Z             T
     1             Z             T             Z             T             Z")
1006      FORMAT (40X,F12.3,2X,F12.3)          
          ! SOIL PROPERTIES      
          WRITE (7120,2100) 
          WRITE (7120,3550)
          ! INTIAL VALUE FOR UNIT
          
          DO 3950 I = 1,INDEX_STADA
               NODEINPUT = CURVE_DATA_TZ(2,I)
               CALL READSOILDATA ("CALO",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               NUNIT  = SBC(14)
               IF (NUNIT.NE.0) EXIT
3950      CONTINUE 
          IF (NUNIT.EQ.1) WRITE (7120,3451) ! ENGLISH UNIT HEADER
          IF (NUNIT.EQ.2) WRITE (7120,3452) ! SI UNIT HEADER
           IF (NUNIT.NE.1.AND.NUNIT.NE.2) THEN
             WRITE (7120,3451)
             IF (CURVE_DATA_TZ(1,I).EQ.0.AND.CURVE_DATA_TZ(2,I).EQ.0) THEN
             WRITE (7120,1005) ! SPACE
             WRITE (7120,1003) ! SPEAC
             WRITE (7120,1003) ! SPACE
             ENDIF
          ENDIF
           INDEX_WRITE = 1D0  
           DO 3100 I = 1,INDEX_STADA
               ISP = CURVE_DATA_TZ(1,I)
               ISN = CURVE_DATA_TZ(2,I)
               !CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               CALL READSOILPROP (ISN,ISP,SBC)
               IF (ISN.EQ.0) GOTO 3100
                     NUNIT  = SBC(14) ! SOIL UNIT
                     NSOILY = SBC(15) ! SOIL TYPE
                     CA     = SBC(16) ! UNDRAIN SHEAR STRENGTH
                     GRAMMA = SBC(17) ! UNIT WEIGHT
                     DELTA  = SBC(18) ! FRICTION ANGLE
                     ANQ    = SBC(19) ! 
                     BETA   = SBC(20) !
                     TFAC   = SBC(21) ! T FACTOR
                     ZFAC   = SBC(22) ! Z FACTOR
                     AOPT   = SBC(23) ! SOIL TYPE
                    
                    IF (SBC(15).EQ.1)THEN
                    SOIL_TYPE_TZ = "      CLAY      "
                    DELTA = 0.D0
                    ENDIF
                    IF (SBC(15).NE.1)THEN
                    IF (AOPT.EQ.1) SOIL_TYPE_TZ = "     GRAVEL     "
                    IF (AOPT.EQ.2) SOIL_TYPE_TZ = "   DENSE SAND   "
                    IF (AOPT.EQ.3) SOIL_TYPE_TZ = "DENSE SAND-SILT "
                    IF (AOPT.EQ.4) SOIL_TYPE_TZ = "MEDIUM SAND-SILT"  
                    IF (AOPT.EQ.5) SOIL_TYPE_TZ = "   MEDIUM SILT  "
                    CA   = 0.0D0
                    ENDIF
                    
                        
                    IF (NGRAV.EQ.1)THEN
                    !AC_OUT_MAX = AX(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AX(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.2)THEN
                    !AC_OUT_MAX = AY(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AY(CURVE_DATA_TZ(2,I))
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.3)THEN
                    !AC_OUT_MAX = AZ(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AZ(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ENDIF
                  
               WRITE (7120,3551) INDEX_WRITE,AC_OUT_MIN,AC_OUT_MAX,SOIL_TYPE_TZ,CA,GRAMMA,DELTA
               INDEX_WRITE =  INDEX_WRITE  + 1
3100       CONTINUE
          WRITE (7120,3552)
2100      FORMAT ("SOIL PROPERTIES FOR T-Z CURVE")    
3550      FORMAT ("STRATA    FROM DEPTH     TO DEPTH        SOIL-TYPE            C         UNIT-WEIGHT  FRICTION-ANGLE" )       
3451  FORMAT ("             (m)           (m)                              N/m^2          N/m^3        DEGREE")      
3452  FORMAT ("             (in)          (in)                            lb/in^2        lb/in^3       DEGREE")  
3551  FORMAT (I5,2X,F12.3,2X,F12.3,4X,A16,2X,F12.5,2X,F12.5,2X,F12.5)   
3552  FORMAT ("") 
      
          ! T-Z CURVE
           DO 5500 I = 1,INDEX_STADA
               NODEINPUT = CURVE_DATA_TZ(2,I)
               CALL READSOILDATA ("CALO",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               NUNIT  = SBC(14)
               IF (NUNIT.NE.0.) EXIT
5500       CONTINUE
           
          WRITE (7120,1010) 
          IF (NUNIT.EQ.1) WRITE (7120,1101)
          IF (NUNIT.EQ.2) WRITE (7120,1102)
          IF (NUNIT.NE.1.AND.NUNIT.NE.2) THEN
             WRITE (7120,1101)
             IF (CURVE_DATA_TZ(1,I).EQ.0.AND.CURVE_DATA_TZ(2,I).EQ.0) THEN
             WRITE (7120,1005)
             WRITE (7120,1003) ! SPACE
             WRITE (7120,1003) ! SPACE
             WRITE (7120,1003) ! SPACE
             ENDIF
          ENDIF
          INDEX_WRITE = 1D0
          DO 6000 I = 1,INDEX_STADA
               IF (CURVE_DATA_TZ(1,I).EQ.0.AND.CURVE_DATA_TZ(2,I).EQ.0) GOTO 6000
               ISP = CURVE_DATA_TZ(1,I)
               ISN = CURVE_DATA_TZ(2,I)
               !CALL READSOILDATA ("CALL",NODEOUT,DFAC,NCODE,NT,SBC(1:37),NODEINPUT,ISN,PROPOUT)
               CALL READSOILPROP (ISN,ISP,SBC)
               IF (NCODE.EQ.1)THEN ! API CURVE
               XFAC   = SBC(12)  ! X AXIS FACTOR
               YFAC   = SBC(13)  ! Y AXIS FACTOR
               IF (CURVE_DATA_TZ(3,I).EQ.16)THEN ! T-Z CLAY
                  !IF (NGRAV.EQ.1)THEN
                  !!AC_OUT_MAX = AX(CURVE_DATA_TZ(1,I))   
                  !!AC_OUT_MIN = AX(CURVE_DATA_TZ(2,I))
                  !NNODEMAX = CURVE_DATA_TZ(1,I)
                  !NNODEMIN = CURVE_DATA_TZ(2,I)
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,1,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,1,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.2)THEN
                  !!AC_OUT_MAX = AY(CURVE_DATA_TZ(1,I))   
                  !!AC_OUT_MIN = AY(CURVE_DATA_TZ(2,I))   
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,2,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,2,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.3)THEN
                  !!AC_OUT_MAX = AZ(CURVE_DATA_TZ(1,I))   
                  !!AC_OUT_MIN = AZ(CURVE_DATA_TZ(2,I)) 
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,3,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  !  CALL RELFILL('@XYZ',AC_OUT_MIN,3,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  !ENDIF
                  
                  IF (NGRAV.EQ.1)THEN
                    !AC_OUT_MAX = AX(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AX(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.2)THEN
                    !AC_OUT_MAX = AY(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AY(CURVE_DATA_TZ(2,I))
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.3)THEN
                    !AC_OUT_MAX = AZ(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AZ(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ENDIF
               WRITE (7120,1001) INDEX_WRITE,YFAC,AC_OUT_MIN,AC_OUT_MAX,CURVE_DATA_TZ(4,I),CURVE_DATA_TZ(12,I)
     1                           ,CURVE_DATA_TZ(5,I),CURVE_DATA_TZ(13,I),CURVE_DATA_TZ(6,I),CURVE_DATA_TZ(14,I)
     1                           ,CURVE_DATA_TZ(7,I),CURVE_DATA_TZ(15,I),CURVE_DATA_TZ(8,I),CURVE_DATA_TZ(16,I)
               WRITE (7120,1051) CURVE_DATA_TZ(9,I),CURVE_DATA_TZ(17,I),CURVE_DATA_TZ(10,I),CURVE_DATA_TZ(18,I)
     1                           ,CURVE_DATA_TZ(11,I),CURVE_DATA_TZ(19,I)
               WRITE (7120,1003) ! SPACE
               INDEX_WRITE =  INDEX_WRITE  + 1
                    ELSEIF (CURVE_DATA_TZ(3,I).EQ.6)THEN ! T-Z SAND
                  ! IF (NGRAV.EQ.1)THEN
                  ! !AC_OUT_MAX = AX(CURVE_DATA_TZ(1,I))   
                  ! !AC_OUT_MIN = AX(CURVE_DATA_TZ(2,I))  
                  ! CALL RELFILL('@XYZ',AC_OUT_MAX,1,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  ! CALL RELFILL('@XYZ',AC_OUT_MIN,1,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  ! ELSEIF (NGRAV.EQ.2)THEN
                  ! !AC_OUT_MAX = AY(CURVE_DATA_TZ(1,I))   
                  ! !AC_OUT_MIN = AY(CURVE_DATA_TZ(2,I))   
                  ! CALL RELFILL('@XYZ',AC_OUT_MAX,2,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  ! CALL RELFILL('@XYZ',AC_OUT_MIN,2,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  ! ELSEIF (NGRAV.EQ.3)THEN
                  ! !AC_OUT_MAX = AZ(CURVE_DATA_TZ(1,I))   
                  ! !AC_OUT_MIN = AZ(CURVE_DATA_TZ(2,I))  
                  ! CALL RELFILL('@XYZ',AC_OUT_MAX,3,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  ! CALL RELFILL('@XYZ',AC_OUT_MIN,3,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  ! ENDIF
                        
                  IF (NGRAV.EQ.1)THEN
                    !AC_OUT_MAX = AX(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AX(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.2)THEN
                    !AC_OUT_MAX = AY(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AY(CURVE_DATA_TZ(2,I))
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.3)THEN
                    !AC_OUT_MAX = AZ(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AZ(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ENDIF
               WRITE (7120,1001) INDEX_WRITE,YFAC,AC_OUT_MIN,AC_OUT_MAX,CURVE_DATA_TZ(4,I),CURVE_DATA_TZ(7,I)
     1                           ,CURVE_DATA_TZ(5,I),CURVE_DATA_TZ(8,I),CURVE_DATA_TZ(6,I),CURVE_DATA_TZ(9,I)
               WRITE (7120,1003) ! SPACE
               INDEX_WRITE =  INDEX_WRITE  + 1
               ENDIF
               
               ELSEIF (NCODE.EQ.3)THEN ! USER DEFINED FUNCTION
               XFAC   = SBC(36)  ! X AXIS FACTOR
               YFAC   = SBC(37)  ! Y AXIS FACTOR
                  !IF (NGRAV.EQ.1)THEN
                  !!AC_OUT_MAX = AX(CURVE_DATA_TZ(1,I))   
                  !!AC_OUT_MIN = AX(CURVE_DATA_TZ(2,I))
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,1,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,1,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.2)THEN
                  !!AC_OUT_MAX = AY(CURVE_DATA_TZ(1,I))   
                  !!AC_OUT_MIN = AY(CURVE_DATA_TZ(2,I)) 
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,2,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,2,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  !ELSEIF (NGRAV.EQ.3)THEN
                  !!AC_OUT_MAX = AZ(CURVE_DATA_TZ(1,I))   
                  !!AC_OUT_MIN = AZ(CURVE_DATA_TZ(2,I))
                  !CALL RELFILL('@XYZ',AC_OUT_MAX,3,CURVE_DATA_TZ(1,I),0)  !GETTING NODAL COORDINATE
                  !CALL RELFILL('@XYZ',AC_OUT_MIN,3,CURVE_DATA_TZ(2,I),0)  !GETTING NODAL COORDINATE
                  !ENDIF
                 IF (NGRAV.EQ.1)THEN
                    !AC_OUT_MAX = AX(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AX(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,1,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,1,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.2)THEN
                    !AC_OUT_MAX = AY(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AY(CURVE_DATA_TZ(2,I))
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,2,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,2,NNODEMIN,0)  !GETTING NODAL COORDINATE
                    ELSEIF (NGRAV.EQ.3)THEN
                    !AC_OUT_MAX = AZ(CURVE_DATA_TZ(1,I))   
                    !AC_OUT_MIN = AZ(CURVE_DATA_TZ(2,I))  
                    NNODEMAX = CURVE_DATA_TZ(1,I)
                    NNODEMIN = CURVE_DATA_TZ(2,I)
                    CALL RELFILL('@XYZ',AC_OUT_MAX,3,NNODEMAX,0)  !GETTING NODAL COORDINATE
                    CALL RELFILL('@XYZ',AC_OUT_MIN,3,NNODEMIN,0)  !GETTING NODAL COORDINATE
                  ENDIF
                  
               WRITE (7120,1001) INDEX_WRITE,YFAC,AC_OUT_MIN,AC_OUT_MAX,CURVE_DATA_TZ(4,I),CURVE_DATA_TZ(20,I)
     1                           ,CURVE_DATA_TZ(5,I),CURVE_DATA_TZ(21,I),CURVE_DATA_TZ(6,I),CURVE_DATA_TZ(22,I)
     1                           ,CURVE_DATA_TZ(7,I),CURVE_DATA_TZ(23,I),CURVE_DATA_TZ(8,I),CURVE_DATA_TZ(24,I)
               WRITE (7120,1004)  CURVE_DATA_TZ(9,I),CURVE_DATA_TZ(25,I)
     1                           ,CURVE_DATA_TZ(10,I),CURVE_DATA_TZ(26,I),CURVE_DATA_TZ(11,I),CURVE_DATA_TZ(27,I)
     1                           ,CURVE_DATA_TZ(12,I),CURVE_DATA_TZ(28,I),CURVE_DATA_TZ(13,I),CURVE_DATA_TZ(29,I)
               WRITE (7120,1004)  CURVE_DATA_TZ(14,I),CURVE_DATA_TZ(30,I)
     1                           ,CURVE_DATA_TZ(15,I),CURVE_DATA_TZ(31,I),CURVE_DATA_TZ(16,I),CURVE_DATA_TZ(32,I)
     1                           ,CURVE_DATA_TZ(17,I),CURVE_DATA_TZ(33,I),CURVE_DATA_TZ(18,I),CURVE_DATA_TZ(34,I)
               WRITE (7120,1006) CURVE_DATA_TZ(19,I),CURVE_DATA_TZ(35,I)
               WRITE (7120,1003) ! SPACE
               INDEX_WRITE =  INDEX_WRITE  + 1
               ENDIF
6000      CONTINUE
              
1051      FORMAT (40X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X,F12.3,2X)
          
          WRITE (7120,"") 
          WRITE (7120,1600) 
1600      FORMAT ("**********************************")
          WRITE (7120,"") 
          ENDIF ! ENDIF PRINT

          END
          
          
C     ----------------------------------------------------------------
      SUBROUTINE READSOILPROP (ISN,ISP,SEBC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 NOPT
      DIMENSION SEBC(37),PROPOUT(3)  
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOILOPT/ NSOIL,NNSO
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
      COMMON /NSOILT/ NNODEA(10000),NPROP(10000)
       
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NA,NEF,NELTOT,NMV,MTYP,ISECT
      
      COMMON /SOILDATA/ NODEA(1000),SBC(37,1000),ITER(1000),NCODE(1000),PROP(2,1000),START_ELV(1000),END_ELV(1000)
     1                 ,NNODE_DIRECT(1000,100),NTOTAL_NODE(100)
      
      DO 100 I =1,NNSO,2
          NIS1 = NNODEA(I)
          NIS2 = NNODEA(I+1)
          IF (NIS1.EQ.ISN.AND.NIS2.EQ.ISP)THEN
             KPROP = NPROP(I) 
          ELSEIF (NIS1.EQ.ISP.AND.NIS2.EQ.ISN)THEN
             KPROP = NPROP(I)
          ENDIF
100   CONTINUE
      
      IF (NCODE(KPROP).EQ.1)THEN
      SEBC(1:37) = SBC(1:37,KPROP)
      ENDIF
      END

C     ----------------------------------------
      SUBROUTINE MODIFY_INPUT_SOIL_SPRING 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
      COMMON /MGRAV/ NGRAV  
      
      DIMENSION N_PILE_GROUP(1000,500)
      DIMENSION NINDEX(500),NSNODE(20)
      DIMENSION N_INDEX (500)
      DIMENSION LAYER(1000,500)
      DIMENSION NCODE_MODITY(1000)
      
      COMMON /SOILOPT/ NSOIL,NNSO
      
      COMMON /NSOILT/ NNODEA(10000),NPROP(10000)
      
      COMMON /CALENODE/ NODECALE1(1000),NODECALE2(1000)
      
      COMMON /SOILDATA/ NODEA(1000),SBC(37,1000),ITER(1000),NCODE(1000),PROP(2,1000),START_ELV(1000),END_ELV(1000)
     1                 ,NNODE_DIRECT(1000,100),NTOTAL_NODE(100)
      
      !READ (ITI,*) SOIL_ELV_BEGIN,SOIL_ELV_END
      N_PILE_GROUP = 0.
      NINDEX       = 0.
      
      CALL LOCATN('OGDM',KGDM,NGDAT,NGIDM,1)  !CALLING TOTAL NUMBER OF GID ELEMENT   ... SONGSAK OCT2019
      CALL INTFILL('%NUB',NGGEO,1,13,0)       !CALLING TOTAL NUMBER OF GEOMETRIC SET ... SONGSAK OCT2019

      INDEX = 1
      N_INDEX = 1D0
      LAYER = 0.
      NNODEA = 0.
      NPROP = 0.
      NCODE_MODITY = 0.
      INDEX_DIRECT = 0.
      DO 100 IJK = 1,NSOIL
      IF (NCODE(IJK).EQ.1) THEN 
      NCODE_MODITY(IJK) = 1D0 
      DO 200 IGM = 1,NGIDM   
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
               CALL INTFILL('OGDM',IGST,3,IGM,0) !CALL GEOMETRIC SET INDEX ... IGST
               CALL CALLNUMNODE_F (IGM,N1,N2) !CALL FRAME CONNECTIVITY NODE
               NSNODE = 0. 
               NSNODE(1) = N1
               NSNODE(2) = N2
               
              CALL NODEABOVE (NSNODE,NODEOUT)
              IF (NODEOUT.EQ.NSNODE(1)) NODEOUT1  = NODEOUT
              IF (NODEOUT.EQ.NSNODE(1)) NODEOUT2  = NSNODE(2)
              IF (NODEOUT.EQ.NSNODE(2)) NODEOUT2  = NSNODE(1)
              IF (NODEOUT.EQ.NSNODE(2)) NODEOUT1  = NODEOUT
              
              IF (NGRAV.EQ.1)THEN
              CALL RELFILL('@XYZ',AY_OUT1,1,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,1,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.2)THEN
              CALL RELFILL('@XYZ',AY_OUT1,2,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,2,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.3)THEN
              CALL RELFILL('@XYZ',AY_OUT1,3,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,3,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ENDIF
          
               IF (AY_OUT1.LE.START_ELV(IJK).AND.AY_OUT2.GE.END_ELV(IJK))THEN
                  N_PILE_GROUP(N_INDEX(IGST),IGST) = IGM ! NUMBER OF ELEMENT
                  LAYER(N_INDEX(IGST),IGST) = IJK
                  N_INDEX(IGST) = N_INDEX(IGST) + 1
               ENDIF  
          ENDIF
200   ENDDO
      ELSEIF (NCODE(IJK).EQ.2) THEN
      !NCODE_MODITY(IGST) = 2D0
      !N_PILE_GROUP(N_INDEX(IGST),1:IGST) = 0. ! NUMBER OF ELEMENT
      !LAYER(N_INDEX(IGST),1:IGST) = 0.
      !N_INDEX(1:IGST) = N_INDEX(1:IGST) + 1
      ELSEIF (NCODE(IJK).EQ.3) THEN
      IF (START_ELV(IJK).EQ.END_ELV(IJK)) EXIT
      NCODE_MODITY(IJK) = 3D0   
      
      DO II = 1,NGIDM
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,II,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
               CALL INTFILL('OGDM',IGST,3,II,0) !CALL GEOMETRIC SET INDEX ... IGST
               CALL CALLNUMNODE_F (II,N1,N2) !CALL FRAME CONNECTIVITY NODE
               NSNODE = 0. 
               NSNODE(1) = N1
               NSNODE(2) = N2
               
              CALL NODEABOVE (NSNODE,NODEOUT)
              IF (NODEOUT.EQ.NSNODE(1)) NODEOUT1  = NODEOUT
              IF (NODEOUT.EQ.NSNODE(1)) NODEOUT2  = NSNODE(2)
              IF (NODEOUT.EQ.NSNODE(2)) NODEOUT2  = NSNODE(1)
              IF (NODEOUT.EQ.NSNODE(2)) NODEOUT1  = NODEOUT
              
              IF (NGRAV.EQ.1)THEN
              CALL RELFILL('@XYZ',AY_OUT1,1,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,1,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.2)THEN
              CALL RELFILL('@XYZ',AY_OUT1,2,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,2,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.3)THEN
              CALL RELFILL('@XYZ',AY_OUT1,3,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,3,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ENDIF

          
               IF (AY_OUT1.LE.START_ELV(IJK).AND.AY_OUT2.GE.END_ELV(IJK))THEN
                  N_PILE_GROUP(N_INDEX(IGST),IGST) = II ! NUMBER OF ELEMENT
                  LAYER(N_INDEX(IGST),IGST) = IJK
                  N_INDEX(IGST) = N_INDEX(IGST) + 1
               ENDIF  
          ENDIF
      ENDDO
      
      ELSEIF (NCODE(IJK).EQ.4) THEN
      NCODE_MODITY(IJK) = 4D0   
      
      DO II = 1,NGIDM
C	CALL GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,II,0)  !IEG
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE  
          IF(ITYP.EQ.5) THEN !FOR FRAME ELEMENT
               CALL INTFILL('OGDM',IGST,3,II,0) !CALL GEOMETRIC SET INDEX ... IGST
               CALL CALLNUMNODE_F (II,N1,N2) !CALL FRAME CONNECTIVITY NODE
               NSNODE = 0. 
               NSNODE(1) = N1
               NSNODE(2) = N2
               
              CALL NODEABOVE (NSNODE,NODEOUT)
              IF (NODEOUT.EQ.NSNODE(1)) NODEOUT1  = NODEOUT
              IF (NODEOUT.EQ.NSNODE(1)) NODEOUT2  = NSNODE(2)
              IF (NODEOUT.EQ.NSNODE(2)) NODEOUT2  = NSNODE(1)
              IF (NODEOUT.EQ.NSNODE(2)) NODEOUT1  = NODEOUT
              
              IF (NGRAV.EQ.1)THEN
              CALL RELFILL('@XYZ',AY_OUT1,1,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,1,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.2)THEN
              CALL RELFILL('@XYZ',AY_OUT1,2,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,2,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ELSEIF (NGRAV.EQ.3)THEN
              CALL RELFILL('@XYZ',AY_OUT1,3,NODEOUT1,0)  !GETTING NODAL COORDINATE
	        CALL RELFILL('@XYZ',AY_OUT2,3,NODEOUT2,0)  !GETTING NODAL COORDINATE
              ENDIF

          
               IF (AY_OUT1.LE.START_ELV(IJK).AND.AY_OUT2.GE.END_ELV(IJK))THEN
                  N_PILE_GROUP(N_INDEX(IGST),IGST) = II ! NUMBER OF ELEMENT
                  LAYER(N_INDEX(IGST),IGST) = IJK
                  N_INDEX(IGST) = N_INDEX(IGST) + 1
               ENDIF  
          ENDIF
      ENDDO
      ENDIF
      
100   ENDDO
     
      INDEX_DIRECT = 1.
      INDEX_PO     = 1
      CALL INTFILL('%NUB',NGGEO,1,13,0)       !CALLING TOTAL NUMBER OF GEOMETRIC SET ... SONGSAK OCT2019
      DO JJ = 1,NGGEO
          DO II = 1, N_INDEX(JJ)-1
              IF (N_PILE_GROUP(II,JJ).NE.0D0)THEN
                  NSNODE = 0.
                  CALL CALLNUMNODE_F (N_PILE_GROUP(II,JJ),NSNODE(1),NSNODE(2))
                  CALL NODEABOVE (NSNODE,NODEOUT)
                  IF (NODEOUT.EQ.NSNODE(1)) NNODEA(INDEX)   =  NSNODE(1)
                  IF (NODEOUT.EQ.NSNODE(1)) NNODEA(INDEX+1) =  NSNODE(2)
                  IF (NODEOUT.EQ.NSNODE(2)) NNODEA(INDEX)   =  NSNODE(2)
                  IF (NODEOUT.EQ.NSNODE(2)) NNODEA(INDEX+1) =  NSNODE(1)
                  NPROP(INDEX) = LAYER(II,JJ)
                  NPROP(INDEX+1) = LAYER(II,JJ)
                  INDEX = INDEX + 2
              ELSEIF (N_PILE_GROUP(II,JJ).EQ.0D0)THEN
              NNODEA(INDEX) = -1.D0
              NPROP(INDEX)  = -1.D0
              INDEX = INDEX + 1
              ENDIF
          ENDDO
      ENDDO
      
      !DO JJ = 1,INDEX - 1
      !IF     (NGRAV.EQ.1) THEN  
      !ELSEIF (NGRAV.EQ.2) THEN    
      !ELSEIF (NGRAV.EQ.3) THEN 
      !ENDIF
      !ENDDO
      
      DO II = 1, NSOIL
          IF (NCODE(II).EQ.2)THEN
              !IJK = II*2D0 - INDEX_DIRECT
          NNODEA(INDEX) = NTOTAL_NODE(II) !NNODE_DIRECT(II)
          NPROP(INDEX)  = II
          INDEX = INDEX + 1D0
          ENDIF
      ENDDO

      
      
      NNSO = INDEX - 1
      
      TCOUNT = 1
      DO I = 2,NNSO,2
          NODECALE1(TCOUNT) = NNODEA(I)
          NODECALE2(TCOUNT) = NNODEA(I-1)
          TCOUNT       = TCOUNT + 1
      ENDDO
      
      
      ! ---- SKIP 09/2021 ----
      !DO I =1,NNSO
      !READ(ITI,*)  NNODEA(I),NPROP(I)
      !ENDDO
      
      ! ---- SKIP 09/2021 ----
       !TCOUNT = 1
       !DO I = 2,NNSO,2
       !    NODECALE1(TCOUNT) = NNODEA(I)
       !    NODECALE2(TCOUNT) = NNODEA(I-1)
       !    TCOUNT       = TCOUNT + 1
       !ENDDO
      
      END
!     --------------------------------------------------------------                
      SUBROUTINE DIRECT_STIFFNESS_SOIL_PILE (NISN,NNODE_DIREC,NNODE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N) 
      DIMENSION NNODE_DIREC(1000)
      COMMON /SOILDATA/ NODEA(1000),SBC(37,1000),ITER(1000),NCODE(1000),PROP(2,1000),START_ELV(1000),END_ELV(1000)
     1                 ,NNODE_DIRECT(1000,100),NTOTAL_NODE(100)
      COMMON /NSOILT/ NNODEA(10000),NPROP(10000)
      
      NNODE_DIREC = 0.
      NPROP_INDEX = NPROP(NISN)
      NNODE = NTOTAL_NODE(NPROP_INDEX)
      
      NNODE_DIREC(1:NNODE) = NNODE_DIRECT(1:NNODE,NPROP_INDEX)
      
      END
!     --------------------------------------------------------------                
      SUBROUTINE PO_TZ (ISP,ISN,PO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N) 
      COMMON /SOILOPT/ NSOIL,NNSO
      
      COMMON /NSOILT/ NNODEA(10000),NPROP(10000)
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
      COMMON /MGRAV/ NGRAV  
      
      COMMON /SOILDATA/ NODEA(1000),SBC(37,1000),ITER(1000),NCODE(1000),PROP(2,1000),START_ELV(1000),END_ELV(1000)
     1                 ,NNODE_DIRECT(1000,100),NTOTAL_NODE(100)
      
      DIMENSION NLAYER(100)
      DIMENSION NSNODE(NSN)
      
      NLAYER = 0.
      NSNODE = 0.
      NSNODE(1) = ISP
      NSNODE(2) = ISN
      IF (NGRAV.EQ.1)THEN
      CALL RELFILL('@XYZ',AVALE1,1,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AVALE2,1,ISN,0)  !GETTING NODAL COORDINATE
      ELSEIF (NGRAV.EQ.2)THEN
      CALL RELFILL('@XYZ',AVALE1,2,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AVALE2,2,ISN,0)  !GETTING NODAL COORDINATE
      ELSEIF (NGRAV.EQ.3)THEN
      CALL RELFILL('@XYZ',AVALE1,3,ISP,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',AVALE2,3,ISN,0)  !GETTING NODAL COORDINATE
      ENDIF
      
      INDEX = 1
      NLAYER = 0.
      DO I = 1,NSOIL
      AVG = (AVALE1+AVALE2)/2D0
          IF (AVG.LE.START_ELV(I))THEN
              IF (AVG.GE.END_ELV(I))THEN
              NLAYER(INDEX) = I
              INDEX = INDEX + 1
              ELSEIF (AVG.LE.END_ELV(I))THEN
              NLAYER(INDEX) = I
              INDEX = INDEX + 1    
              ENDIF
          ENDIF
      ENDDO
      
      PO = 0.0D0
      DO I = 1,INDEX-1
          CALL NODEABOVE (NSNODE,NODEOUT)
          IF (NODEOUT.EQ.NSNODE(1)) NMAX =  NSNODE(1)
          IF (NODEOUT.EQ.NSNODE(1)) NMIN =  NSNODE(2)
          IF (NODEOUT.EQ.NSNODE(2)) NMAX =  NSNODE(2)
          IF (NODEOUT.EQ.NSNODE(2)) NMIN =  NSNODE(1)
          
               IF (NGRAV.EQ.1)THEN
               CALL RELFILL('@XYZ',AVALE1,1,NMAX,0)  !GETTING NODAL COORDINATE
	         CALL RELFILL('@XYZ',AVALE2,1,NMIN,0)  !GETTING NODAL COORDINATE
               ELSEIF (NGRAV.EQ.2)THEN
               CALL RELFILL('@XYZ',AVALE1,2,NMAX,0)  !GETTING NODAL COORDINATE
	         CALL RELFILL('@XYZ',AVALE2,2,NMIN,0)  !GETTING NODAL COORDINATE
               ELSEIF (NGRAV.EQ.3)THEN
               CALL RELFILL('@XYZ',AVALE1,3,NMAX,0)  !GETTING NODAL COORDINATE
	         CALL RELFILL('@XYZ',AVALE2,3,NMIN,0)  !GETTING NODAL COORDINATE
               ENDIF
          IF (START_ELV(NLAYER(I)).GE.0.0D0.AND.END_ELV(NLAYER(I)).GE.0.0)THEN   !! + +
             IF (START_ELV(NLAYER(I)).GE.AVALE1.AND.END_ELV(NLAYER(I)).LE.AVALE1) THEN
             ALENGTH = START_ELV(NLAYER(I)) - AVALE1
                     GRAMMA = SBC(17,NLAYER(I))
                     IF (I.NE.INDEX-1) THEN
                     PO     = PO + GRAMMA*ALENGTH 
                     ELSEIF (I.EQ.INDEX-1) THEN ! LAST LAYER
                     PO     = PO + GRAMMA*ALENGTH + GRAMMA*(AVALE1-AVALE2)/2D0 
                     ENDIF
             ELSE
                  
             ALENGTH = START_ELV(NLAYER(I)) - END_ELV(NLAYER(I))
                     GRAMMA = SBC(17,NLAYER(I))
                     IF (I.NE.INDEX-1) THEN
                     PO     = PO + GRAMMA*ALENGTH 
                     ELSEIF (I.EQ.INDEX-1) THEN ! LAST LAYER
                     PO     = PO + GRAMMA*ALENGTH/2D0 
                     ENDIF
             ENDIF
          
          ELSEIF (START_ELV(NLAYER(I)).LE.0.0D0.AND.END_ELV(NLAYER(I)).LE.0.0)THEN !! - - 
              IF (ABS(START_ELV(NLAYER(I))).LE.ABS(AVALE2).AND.ABS(END_ELV(NLAYER(I))).GE.ABS(AVALE2)) THEN
              ALENGTH = ABS(ABS(START_ELV(NLAYER(I))) - ABS(AVALE1))
                     GRAMMA = SBC(17,NLAYER(I))
                     IF (I.NE.INDEX-1) THEN
                     PO     = PO + GRAMMA*ALENGTH 
                     ELSEIF (I.EQ.INDEX-1) THEN ! LAST LAYER
                     PO     = PO + GRAMMA*ALENGTH + GRAMMA*(ABS(AVALE2) - ABS(AVALE1))/2D0 
                     ENDIF
             ELSE
                  
             ALENGTH = START_ELV(NLAYER(I)) - END_ELV(NLAYER(I))
                     GRAMMA = SBC(17,NLAYER(I))
                     IF (I.NE.INDEX-1) THEN
                     PO     = PO + GRAMMA*ALENGTH 
                     ELSEIF (I.EQ.INDEX-1) THEN ! LAST LAYER
                     PO     = PO + GRAMMA*ALENGTH/2D0 
                     ENDIF
             ENDIF
        ELSEIF (START_ELV(NLAYER(I)).GE.0.0D0.AND.END_ELV(NLAYER(I)).LE.0.0)THEN  !! + -
               IF (START_ELV(NLAYER(I)).GE.AVALE1.AND.END_ELV(NLAYER(I)).LE.AVALE1) THEN
               ALENGTH = START_ELV(NLAYER(I)) - AVALE1
                     GRAMMA = SBC(17,NLAYER(I))
                     IF (I.NE.INDEX-1) THEN
                     PO     = PO + GRAMMA*ALENGTH 
                     ELSEIF (I.EQ.INDEX-1) THEN ! LAST LAYER
                     PO     = PO + GRAMMA*ALENGTH + GRAMMA*(AVALE1-AVALE2)/2D0 
                     ENDIF
             ELSE
                  
             ALENGTH = START_ELV(NLAYER(I)) - END_ELV(NLAYER(I))
                     GRAMMA = SBC(17,NLAYER(I))
                     IF (I.NE.INDEX-1) THEN
                     PO     = PO + GRAMMA*ALENGTH 
                     ELSEIF (I.EQ.INDEX-1) THEN ! LAST LAYER
                     PO     = PO + GRAMMA*ALENGTH/2D0 
                     ENDIF
             ENDIF
             
        ENDIF
          
      ENDDO
      
      END