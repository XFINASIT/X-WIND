      SUBROUTINE WRITEJOINTWAVEFORCE(MLE,CHXR,CHYR,CHZR,CHXR2,CHYR2,CHZR2,SEABED)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /WAVEWAITEDATAREAL/ CHXRI,CHYRI,CHZRI,CHXRI2,CHYRI2,CHZRI2,SEABEDX
      COMMON /WAVEWAITEDATAINTEGER/ MLEI
      
      MLEI = MLE 
      CHXRI = CHXR
      CHYRI = CHYR
      CHZRI = CHZR
      CHXRI2 = CHXR2
      CHYRI2 = CHYR2
      CHZRI2 = CHZR2
      SEABEDX = SEABED
      
      CALL JACOBIAN_ELIPFUNCTION

      
            
      END SUBROUTINE
!     ============================================================================================================ 
      SUBROUTINE WRITE_WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,VH1,VGV,VH2,YL,XYZ,MLE,IORRE,IRWAVE,
     1                ITIME,ILCAS)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      DIMENSION VH1(3),VGV(3),WFX(2),WFY(2),WFZ(2),WRTELEV(2),VH(2)
      
      COMMON /WAVEWAITEDATAREAL/ CHXRI,CHYRI,CHZRI,CHXRI2,CHYRI2,CHZRI2,SEABED
      COMMON /WAVEWAITEDATAINTEGER/ MLEI
      
      COMMON /MGRAV/ NGRAV 
      
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
      DIMENSION XYZ(NCO*NNM,NELE)
          
      WAVESURVACE = YL+SEABED
      
      IF (ISOLOP.EQ.1)THEN
         IF (ILCAS.EQ.1) THEN
         NEC = MLEI
         ELSEIF (ILCAS.GT.1D0)THEN
         NEC = ILCAS*NELE + MLEI
         ENDIF 
      ELSEIF (ISOLOP.EQ.5)THEN
         IF (ITIME.EQ.1) THEN
         NEC = MLEI
         ELSEIF (ITIME.GT.1D0)THEN
         NEC = NELE*ITIME + MLEI
         ENDIF 
      ELSEIF (ISOLOP.NE.0)THEN
         IF (ITIME.EQ.1) THEN
         NEC = MLEI
         ELSEIF (ITIME.GT.1D0)THEN
         NEC = NELE*ITIME + MLEI
         ENDIF 
      ENDIF
      
      if( CHZRI.GE. WAVESURVACE .AND. CHZRI2.GE. WAVESURVACE ) THEN
          
!      WRITE(5027,1) MLEI, CHZRI , 0 , 0 , 0 , 0 , 0 , 0  , CHZRI2 ! SKIP BY TOEY 03/23
      
      WRITE (812,REC=NEC) MLEI, CHZRI , 0 , 0 , 0 , 0 , 0 , 0  , CHZRI2 
      ELSE 

      CHANA = 3
      
      DO I = 1,2
          
      IF( NGRAV .EQ. 3)THEN
         IF(I.EQ.1)THEN 
         H1R = CHXRI*VH1(1)
         HGM = (CHZRI-SEABED)*VGV(3)
         HEM1 = HGM
         ELSEIF(I.EQ.2)THEN
         H1R = CHXRI2*VH1(1)
         HGM = (CHZRI2-SEABED)*VGV(3) 
         HGM2 = HGM
         ENDIF
      ELSEIF( NGRAV .EQ. 2)THEN
         IF(I.EQ.1)THEN 
         H1R = CHXRI
         HGM = (CHYRI-SEABED)*VGV(3)
         HEM1 = HGM
         ELSEIF(I.EQ.2)THEN
         H1R = CHXRI2
         HGM = (CHYRI2-SEABED)*VGV(2) 
         HGM2 = HGM
         ENDIF    
      ENDIF
                   
		 
         WRTELEV(I) = HGM
      CALL WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,IORRE,IRWAVE)
      
!      WFX(I) = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
!      WFY(I) = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
!      WFZ(I) = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2
      
       WFX(I) = WH1
       WFY(I) = WGV
       WFZ(I) = WH2
       
      ! CUTTING FORCE
      CALL CUTTING_FORCE(XYZ,'X',MLE,WFX(I))
      CALL CUTTING_FORCE(XYZ,'Y',MLE,WFY(I))
      CALL CUTTING_FORCE(XYZ,'Z',MLE,WFZ(I))
      
      IF (HGM.GT.YL) THEN
      WFX(I) = 0.0D0
      WFY(I) = 0.0D0
      WFZ(I) = 0.0D0
      ENDIF
      
      IF( IWAVE . EQ . 3 ) THEN
      
          APRECENT = 2.50D0
          
          SURFACEMAX = Yl
          
          AWAVE_ROTIO = HGM / Yl
          
          AFACTORWAVE = ( 1.0D0 + AWAVE_ROTIO * (APRECENT / 100.0D0 ) )
          
            WFX(I) = WFX(I) *AFACTORWAVE
            WFY(I) = WFY(I) *AFACTORWAVE
            WFZ(I) = WFZ(I) *AFACTORWAVE
      
      ELSEIF( IWAVE . EQ . 4 ) THEN
          
          STCvalue = 1.45
          Endvalue =0.3
          
          !=$S$1-(F3/13*$T$1)
          
          AFACTORWAVE = STCvalue - ( HGM / Yl ) * Endvalue
          
          
          IF( HGM.GT.0.4D0*Yl.AND. HGM.LE.0.5D0*Yl ) THEN
              
              AFACTORWAVE = AFACTORWAVE *  0.990D0
          
          ENDIF
          
          IF( HGM.GT.0.5D0*Yl.AND. HGM.LT.0.95D0*Yl ) THEN
              
              AFACTORWAVE = AFACTORWAVE *  0.970D0
          
          ENDIF
          
            WFX(I) = WFX(I) *AFACTORWAVE
            WFY(I) = WFY(I) *AFACTORWAVE
            WFZ(I) = WFZ(I) *AFACTORWAVE
          
      ENDIF
            
      
      ENDDO
      
      
      
!      WRITE(*,1) MLEI, HGM1 , WFX(1) , WFY(1) , WFZ(1) , WFX(2) , WFY(2) , WFZ(2) , HGM2
!      WRITE(5027,1) MLEI, HGM1 , WFX(1) , WFY(1) , WFZ(1) , WFX(2) , WFY(2) , WFZ(2) , HGM2 , YL
!      WRITE(5027,1) MLEI, WRTELEV(1) , WFX(1) , WFY(1) , WFZ(1) , WFX(2) , WFY(2) , WFZ(2) , WRTELEV(2) , YL ! SKIP BY TOEY 03/23
      
      WRITE (812,REC=NEC) MLEI, HGM1 , WFX(1) , WFY(1) , WFZ(1) , WFX(2) , WFY(2) , WFZ(2) , HGM2 , YL
      
      !WRITE (812,REC=NEC+NELE) MLEI, WRTELEV(1) , WFX(1) , WFY(1) , WFZ(1) , WFX(2) , WFY(2) , WFZ(2) , WRTELEV(2) , YL
1     FORMAT(I5,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,
     1 X,F12.4)

      ENDIF

            
      END SUBROUTINE
!     ============================================================================================================ 
      SUBROUTINE WRITEJOINTWAVEFORCEHEADER
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /WAVEWAITEDATAREAL/ CHXRI,CHYRI,CHZRI,CHXRI2,CHYRI2,CHZRI2
      COMMON /WAVEWAITEDATAINTEGER/ MLEI
      
!      WRITE(5027,1) ! SKIP BY TOEY 03/23
!1     FORMAT('ELEMENT   ELEVATION     WFX_JOINT1     WFY_JOINT1    WFZ_JOINT1   WFX_JOINT2    WFY_JOINT2    WFZ_JOINT2 ') ! SKIP BY TOEY 03/23
      
      END SUBROUTINE
!     ============================================================================================================ 
      SUBROUTINE PRNOUTWAVEFORCE(ISTEP)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*200 NAMEX
      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON /OSPEC/ ISTIME,IOUTSPEC,NFATIGUE
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !
C     ----------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,                      
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,                  
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC  
      
      COMMON/OFFSHORE_CASE/ LGEN,IOFFL,IORI_OFFSHORE
      COMMON/WAVE_DISTRIBUTED/ ICASE
      
      IF (ISOLOP.EQ.1) THEN
        IF (ISTEP.LE.LGEN) RETURN
      ENDIF
      IF (IOFFL.EQ.0) RETURN
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      JSTEP = ISTEP
      CALL PRNLCHD(NAMEX,NAML,JSTEP,INM,IND)
      WRITE (LFPRIN,101) NAMEX(1:NAML),JSTEP
      
101	FORMAT(/'Result "Wave Distributed Force"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints  "XWave"'/,
	2		'ComponentNames  "F-X" "F-Y" "F-Z"'/,
	3		'Values')

      CALL READNODALWAVEFORCE(ISTEP,ICASE)
      ICASE = ICASE + 1
      
      
      RETURN
      END SUBROUTINE
C     ----------------------------------------------------------------      
      SUBROUTINE READNODALWAVEFORCE(ISTEP,ICASE)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,                      
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,                  
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT   
            
      COMMON /MWRITETIONWAVEFORCE/ NWELEMENT
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON /OSPEC/ ISTIME,IOUTSPEC,NFATIGUE
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 ! 
C     ----------------------------------------------------------------
!      DIMENSION MELEMENT(NELE*5),WFOCX1(NELE*5),WFOCY1(NELE*5),WFOCZ1(NELE*5),WFOCX2(NELE*5),WFOCY2(NELE*5),WFOCZ2(NELE*5) 
!      DIMENSION NEWMELEMENT(NELE),ANEWWFOCX1(NELE),ANEWWFOCY1(NELE),ANEWWFOCZ1(NELE)
!      DIMENSION ANEWWFOCX2(NELE),ANEWWFOCY2(NELE),ANEWWFOCZ2(NELE)
      
!      ALLOCATABLE MELEMENT(:) , WFOCX1(:) , WFOCY1(:) , WFOCZ1(:) , WFOCX2(:) , WFOCY2(:) , WFOCZ2(:)
!      ALLOCATABLE NEWMELEMENT(:) , ANEWWFOCX1(:) , ANEWWFOCY1(:) , ANEWWFOCZ1(:) , ANEWWFOCX2(:) , ANEWWFOCY2(:) , ANEWWFOCZ2(:)
!      
!      ALLOCATE( MELEMENT(NELE*5) , WFOCX1(NELE*5) , WFOCY1(NELE*5) , WFOCZ1(NELE*5) , WFOCX2(NELE*5) , WFOCY2(NELE*5),
!     1  WFOCZ2(NELE*5), 
!     1  NEWMELEMENT(NELE)  , ANEWWFOCX1(NELE) , ANEWWFOCY1(NELE) , ANEWWFOCZ1(NELE) , ANEWWFOCX2(NELE) , ANEWWFOCY2(NELE) ,
!     1  ANEWWFOCZ2(NELE) )
      
      
!      REWIND(5027) ! SKIP BY TOEY 03/23
!      
!      READ(5027,*,IOSTAT = IOS) ! SKIP BY TOEY 03/23
!      
!      ICOUNT = 1
!      DO I = 1,NWELEMENT*5
!      READ(5027,*,END=100,IOSTAT = IOS) MLEI, HGM , WFX1 , WFY1 , WFZ1, WFX2 , WFY2 , WFZ2 ! SKIP BY TOEY 03/23
!      
!      MELEMENT(I) = MLEI
!      WFOCX1(I)   = WFX1
!      WFOCY1(I)   = WFY1
!      WFOCZ1(I)   = WFZ1
!      WFOCX2(I)   = WFX2
!      WFOCY2(I)   = WFY2
!      WFOCZ2(I)   = WFZ2
!      
!1     FORMAT(I5,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4,2X,F12.4)
!      
!      
!      !WRITE(LFPRIN,200) MLEI , WFX1 , WFY1 , WFZ1 !, WFX2 , WFY2 , WFZ2
!      !WRITE(LFPRIN,201) WFX2 , WFY2 , WFZ2
!      
!       ICOUNT = ICOUNT + 1
!       
!      ENDDO 
!       
!100   CONTINUE
200   FORMAT(2X,I5,10E15.6)
201   FORMAT(7X,10E15.6)
!      
!      ! REARANGE DATA
!      
!      NEWMELEMENT = -100
      
!      DO J = 1,NWELEMENT
!          DO 400 K=1,(ICOUNT - 1)
!           IF(J.EQ.MELEMENT(K))THEN
!               
!           !WRITE(LFPRIN,200) MELEMENT(K) , WFOCX1(K) , WFOCY1(K) , WFOCZ1(K)
!           !WRITE(LFPRIN,201)  WFOCX2(K) , WFOCY2(K) , WFOCZ2(K)
!               NEWMELEMENT(J) = MELEMENT(K) 
!               ANEWWFOCX1(J)   = WFOCX1(K)
!               ANEWWFOCY1(J)   = WFOCY1(K)
!               ANEWWFOCZ1(J)   = WFOCZ1(K)
!               ANEWWFOCX2(J)   = WFOCX2(K)
!               ANEWWFOCY2(J)   = WFOCY2(K)
!               ANEWWFOCZ2(J)   = WFOCZ2(K)
!               EXIT
!           ENDIF         
!400     CONTINUE
!      ENDDO
      
      ! WRITE DATA
      DO J = 1,NWELEMENT
          
      IF (ISOLOP.EQ.1)THEN
         IF (ICASE.EQ.1) THEN
         NEC = J
         ELSEIF (ICASE.GT.1D0)THEN
         NEC = ICASE*NELE + J
         ENDIF 
      ELSEIF (ISOLOP.EQ.5)THEN
         IF (ISTEP.EQ.1) THEN
         NEC = J
         ELSEIF (ISTEP.GT.1D0)THEN
         NEC = NELE*ISTEP + J
         ENDIF 
      ELSEIF (ISOLOP.NE.0)THEN
         IF (ISTEP.EQ.1) THEN
         NEC = J
         ELSEIF (ISTEP.GT.1D0)THEN
         NEC = NELE*ISTEP + J
         ENDIF 
      ENDIF
      MLEI = 0.
      HGM  = 0.
      WFX1 = 0.
      WFY1 = 0.
      WFZ1 = 0.
      WFX2 = 0.
      WFY2 = 0.
      WFZ2 = 0.
      READ (812,REC=NEC,IOSTAT = IOS) MLEI, HGM , WFX1 , WFY1 , WFZ1, WFX2 , WFY2 , WFZ2 ! SKIP BY TOEY 03/23
        !IF((J).EQ.NEWMELEMENT(J))THEN
            
            !WRITE(LFPRIN,200) NEWMELEMENT(J) ,ANEWWFOCX1(J),ANEWWFOCY1(J),ANEWWFOCZ1(J)
            !WRITE(LFPRIN,201) ANEWWFOCX2(J),ANEWWFOCY2(J),ANEWWFOCZ2(J) 
            
        !ELSE
            
            !WRITE(LFPRIN,200) J ,0,0,0
            !WRITE(LFPRIN,201) 0,0,0
            
        !ENDIF
      
            IF (MLEI.NE.0.D0)THEN
            WRITE(LFPRIN,200) J ,WFX1,WFY1,WFZ1
            WRITE(LFPRIN,201) WFX2,WFY2,WFZ2
            ELSEIF (MLEI.EQ.0.D0)THEN
            WRITE(LFPRIN,200) J ,0,0,0
            WRITE(LFPRIN,201) 0,0,0 
            ENDIF
          
      ENDDO
  
      
      !IF (NWELEMENT.EQ.0)THEN
      !    DO J = 1,NELE
      !      WRITE(LFPRIN,200) J ,0,0,0
      !      WRITE(LFPRIN,201) 0,0,0
      !    ENDDO
      !ENDIF

      WRITE(LFPRIN,300)
300   FORMAT('End Values')      
      WRITE(LFPRIN,*)
      
      
!      DEALLOCATE( MELEMENT , WFOCX1 , WFOCY1 , WFOCZ1 , WFOCX2 , WFOCY2 , WFOCZ2 ,
!     1 NEWMELEMENT , ANEWWFOCX1 , ANEWWFOCY1 , ANEWWFOCZ1  , ANEWWFOCX2 , ANEWWFOCY2 , ANEWWFOCZ2 )
      
      END SUBROUTINE
          
C     ----------------------------------------------------------------
