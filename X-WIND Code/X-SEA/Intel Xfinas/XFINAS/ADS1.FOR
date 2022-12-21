C --- SUBROUTINE ADS (INFO,ISTRAT,IOPT,IONED,IPRINT,IGRAD,NDV,NCON,     ADS00010
C    1 X,VLB,VUB,OBJ,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,NRWK,IWK,NRIWK)     ADS00020
C     DIMENSION X(NDV),VLB(NDV),VUB(NDV),G(NCON),IDG(NCON),IC(NCOLA),   ADS00030
C    1 DF(NDV),A(NRA,NCOLA),WK(NRWK),IWK(NRIWK)                         ADS00040
C     ------------------------------------------------------------------ADS00050
C --- MAIN CONTROL ROUTINE FOR THE ADS DESIGN OPTIMIZATION PROGRAM.     ADS00060
C     ------------------------------------------------------------------ADS00070
C     BY G. N. VANDERPLAATS                               JUNE 1982.    ADS00080
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS00090
C --- INPUT.                                                            ADS00100
C     INFO - CONTROL PARAMETER.                                         ADS00110
C        -2 - SET SCALAR PARAMETERS TO DEFAULTS, DEFINE WORK ARRAY      ADS00120
C             STORAGE AND EXIT.                                         ADS00130
C        -1 - FIRST CALL AFTER INFO = -2 CALL.  USER HAS OVER-RIDDEN    ADS00140
C             SOME DEFAULTS.  NOW PROCEED WITH OPTIMIZATION.            ADS00150
C         0 - FIRST CALL.  SET ALL SCALARS TO DEFAULTS, DEFINE STORAGE  ADS00160
C             AND PROCEED WITH OPTIMIZATION.                            ADS00170
C      GT.0 - OPTIMIZATION IS IN PROGRESS.  THIS IS A RE-CALL.          ADS00180
C     ISTRAT - OPTIMIZATION STRATEGY IDENTIFIER.                        ADS00190
C     IOPT - OPTIMIZER IDENTIFIER.                                      ADS00200
C     IONED - ONE-DIMENSIONAL SEARCH IDENTIFIER.                        ADS00210
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS00220
C     NCON - NUMBER OF CONSTRAINTS.  MAY BE ZERO.                       ADS00230
C     IGRAD - GRADIENT CALCULATION IDENTIFIER.                          ADS00240
C         0 - ALL GRADIENTS ARE CALCULATED USING FINITE DIFFERENCE.     ADS00250
C         1 - ALL GRADIENTS ARE SUPPLIED BY THE USER.                   ADS00260
C     IPRINT - PRINT CONTROL.                                           ADS00270
C     X - VECTOR OF DESIGN VARIABLES.  INITIAL DESIGN IF FIRST CALL.    ADS00280
C         PREVIOUSLY RETURNED DESIGN IS SUBSEQUENT CALL.  MINIMUM       ADS00290
C         DIMENSION IS NDV.                                             ADS00300
C     OBJ - VALUE OF THE OBJECTIVE FUNCTION ASSOCIATED WITH X.          ADS00310
C     VLB - LOWER BOUNDS ON THE DESIGN VARIABLES.  MINIMUM DIMENSION    ADS00320
C           IS NDV.                                                     ADS00330
C     VUB - UPPER BOUNDS ON THE DESIGN VARIABLES.  MINIMUM DIMENSION    ADS00340
C           IS NDV.                                                     ADS00350
C     G - CONSTRAINT VALUES ASSOCIATED WITH CURRENT DESIGN.  MINIMUM    ADS00360
C         DIMENSION IS NCON.                                            ADS00370
C     IDG - CONSTRAINT TYPE IDENTIFIER.  MINIMUM DIMENSION IS NCON.     ADS00380
C     IC - ARRAY IDENTIFYING CONSTRAINTS FOR WHICH GRADIENTS WERE       ADS00390
C          REQUIRED ON LAST RETURN FROM ADS. MINIMUM DIMENSION IS NCOLA.ADS00400
C     DF - GRADIENT OF OBJECTIVE.  MINIMUM DIMENSION IS NDV.            ADS00410
C     A - ARRAY CONTAINING GRADIENTS OF CONSTRAINTS IDENTIFIED BY ARRAY ADS00420
C         IC.  GRADIENTS ARE STORED IN COLUMNS OF A.  DIMENSIONED SIZE  ADS00430
C         IS NRA X NCOLA.                                               ADS00440
C     NRA - NUMBER OF ROWS IN ARRAY A.  NRA.GE.NDV.                     ADS00450
C     NCOLA - NUMBER OF COLUMNS OF ARRAY A.  NCOLA.GE.NGT.              ADS00460
C     WK - REAL WORK ARRAY.                                             ADS00470
C     NRWK - DIMENSIONED SIZE OF WK.                                    ADS00480
C     IWK INTEGER WORK ARRAY.                                           ADS00490
C     NRIWK - DIMENSIONED SIZE OF IWK.                                  ADS00500
C --- OUTPUT.                                                           ADS00510
C     INFO - REQUIRED INFORMATION TO BE PROVIDED ON NEXT CALL.          ADS00520
C       -2 - INITIALIZATION IS COMPLETE.  SET INFO=-1 AND CALL AGAIN    ADS00530
C            FOR OPTIMIZATION.                                          ADS00540
C        0 - OPTIMIZATION IS COMPLETE.                                  ADS00550
C        1 - EVALUATE OBJECTIVE AND CONSTRAINTS AND CALL AGAIN.         ADS00560
C        2 - EVALUATE GRADIENT OF OBJECTIVE AND THE NGT SUBSET OF       ADS00570
C            CONSTRAINTS AND CALL AGAIN.                                ADS00580
C     OBJ - CURRENT OR OPTIMUM OBJECTIVE FUNCTION VALUE.                ADS00590
C     X - VECTOR OF DESIGN VARIABLES.                                   ADS00600
C     G - VECTOR OF CONSTRAINTS.                                        ADS00610
C     NGT - NUMBER OF CURRENTLY ACTIVE/VIOLATED CONSTRAINTS.  IGNORE    ADS00620
C          IF INFO.LE.2.                                                ADS00630
C     IC - ARRAY IDENTIFYING THE NGT CONSTRAINTS FOR WHICH GRADIENTS    ADS00640
C          ARE REQUIRED.                                                ADS00650
C                                                                       ADS00660
      SUBROUTINE ADS (INFO,ISTRAT,IOPT,IONED,IPRINT,IGRAD,NDV,NCON,X,VLBADS00670
     1,VUB,OBJ,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,NRWK,IWK,NRIWK)           ADS00680
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

      DIMENSION X(1),VLB(1),VUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WK(1ADS00690
     1),IWK(1)                                                          ADS00700
C                                                                       ADS00710

C      print*,'> > > obj g1 g2=',obj,g(1),g(2)

      
C     PRINT CONTROLS.                                                   ADS00720
      IPGLOB=IPRINT/1000                                                ADS00730
      I=MOD(IPRINT,1000)                                                ADS00740
      IPRNT3=I/100                                                      ADS00750
      IPRNT2=(I-100*IPRNT3)/10                                          ADS00760
      IPRNT1=MOD(I,10)                                                  ADS00770
      NCOLA1=NCOLA                                                      ADS00780
      IF(ISTRAT.EQ.7) NCOLA1=NCOLA1-1                                   ADS00790
      NGT=MIN0(NCON,NCOLA1)                                             ADS00800
      IF (INFO.GT.0) GO TO 70                                           ADS00810
      IF (IPGLOB.EQ.0.OR.INFO.EQ.-1) GO TO 10                           ADS00820

C	SUPPRESS BY SONGSAK NOV2006
C      WRITE (6,360)                                                     ADS00830
C      WRITE (6,370)                                                     ADS00840
C      WRITE (6,380) ISTRAT,IOPT,IONED,IPRINT,IGRAD,NDV,NCON             ADS00850

10    CONTINUE                                                          ADS00860
C     CHECK FOR VALID COMBINATIONS OF OPTIONS.                          ADS00870
      IER=0                                                             ADS00880
C     ISTRAT/IOPT COMBINATION.                                          ADS00890
      IF(ISTRAT.LT.0.OR.ISTRAT.GT.9) WRITE(6,420)ISTRAT                 ADS00900
      IF(ISTRAT.LT.0.OR.ISTRAT.GT.9) IER=1                              ADS00910
      IF(IOPT.LT.0.OR.IOPT.GT.5) WRITE(6,430)IOPT                       ADS00920
      IF(IOPT.LT.0.OR.IOPT.GT.5) IER=1                                  ADS00930
      IF(IONED.LT.0.OR.IONED.GT.8) WRITE(6,440)IONED                    ADS00940
      IF(IONED.LT.0.OR.IONED.GT.8) IER=1                                ADS00950
      IF(ISTRAT.GT.0) GO TO 12                                          ADS00960
      IF(IOPT.LT.4.AND.NCON.EQ.0) GO TO 14                              ADS00970
      IF(IOPT.GE.4.AND.NCON.GT.0) GO TO 14                              ADS00980
      WRITE(6,390)ISTRAT,IOPT                                           ADS00990
      IER=1                                                             ADS01000
      GO TO 14                                                          ADS01010
12    IF(ISTRAT.GT.5) GO TO 13                                          ADS01020
      IF(IOPT.LT.4) GO TO 13                                            ADS01030
      WRITE(6,390)ISTRAT,IOPT                                           ADS01040
      IER=1                                                             ADS01050
      GO TO 14                                                          ADS01060
13    IF(ISTRAT.LT.6) GO TO 14                                          ADS01070
      IF(IOPT.GT.3) GO TO 14                                            ADS01080
      WRITE(6,390)ISTRAT,IOPT                                           ADS01090
      IER=1                                                             ADS01100
C     CHECK IOPT/IONED COMBINATION.                                     ADS01110
14    IF(IOPT.GT.3) GO TO 15                                            ADS01120
      IF(IONED.LE.4) GO TO 15                                           ADS01130
      WRITE(6,400)IOPT,IONED                                            ADS01140
      IER=1                                                             ADS01150
      GO TO 16                                                          ADS01160
15    IF(IOPT.LT.4) GO TO 16                                            ADS01170
      IF(IONED.GT.4) GO TO 16                                           ADS01180
      WRITE(6,400)IOPT,IONED                                            ADS01190
      IER=1                                                             ADS01200
16    IF(IER.GT.0) WRITE(6,410)                                         ADS01210
      IF(IER.GT.0) INFO=0                                               ADS01220
      IF(IER.GT.0) RETURN                                               ADS01230
      IFCALL=0                                                          ADS01240
      IGCALL=0                                                          ADS01250
      IF (INFO.EQ.-1) GO TO 50                                          ADS01260
C     INFO = -2 OR 0.  PROCEED WITH INITIALIZATION TO DEFAULTS.         ADS01270
C     ------------------------------------------------------------------ADS01280
C     FIRST CALL.  INITIALIZE ALL NEEDED INFORMATION.                   ADS01290
C     ------------------------------------------------------------------ADS01300
C                                                                       ADS01310
C     ZERO OUT WK AND IWK ARRAYS.                                       ADS01320
      DO 20 I=1,NRWK                                                    ADS01330
20    WK(I)=0.                                                          ADS01340
      DO 30 I=1,NRIWK                                                   ADS01350
30    IWK(I)=0.                                                         ADS01360
C     CHECK FOR SIDE CONSTRAINT VIOLATIONS.                             ADS01370
      DO 35 I=1,NDV                                                     ADS01380
      IF(X(I).GE.VLB(I)) GO TO 36                                       ADS01390
      WRITE(6,37)I,X(I),VLB(I),I                                        ADS01400
      X(I)=VLB(I)                                                       ADS01410
36    IF(X(I).LE.VUB(I)) GO TO 35                                       ADS01420
      WRITE(6,38)I,X(I),VUB(I),I                                        ADS01430
      X(I)=VUB(I)                                                       ADS01440
35    CONTINUE                                                          ADS01450
37    FORMAT (//5X,14H* * * VARIABLE,I5,                                ADS01460
     *29H IS LESS THAN ITS LOWER BOUND/5X,3HX =,1PE12.5,5X,5HVLB =,     ADS01470
     *1PE12.5/5X,35HX(I) IS SET EQUAL TO VLB(I) FOR I =,I5)             ADS01480
38    FORMAT(//5X,14H* * * VARIABLE,I5,                                 ADS01490
     *32H IS GREATER THAN ITS UPPER BOUND/5X,3HX =,1PE12.5,5X,5HVUB =,  ADS01500
     *1PE12.5/5X,35HX(I) IS SET EQUAL TO VUB(I) FOR I =,I5)             ADS01510
C --- IDENTIFY EXISTANCE OF EQUALITY CONSTRAINTS.                       ADS01520
      NEQCON=0                                                          ADS01530
      IF(NCON.EQ.0) GO TO 41                                            ADS01540
      DO 42 I=1,NCON                                                    ADS01550
      IF(IDG(I).LT.0) NEQCON=NEQCON+1                                   ADS01560
42    CONTINUE                                                          ADS01570
41    CONTINUE                                                          ADS01580
      IWK(71)=NEQCON                                                    ADS01590
C     SCALAR DEFAULTS.                                                  ADS01600
      CALL ADS005 (ISTRAT,IOPT,IONED,NDV,WK,IWK)                        ADS01610
C     WORK ARRAY STORAGE ALLOCATIONS.                                   ADS01620
      CALL ADS004 (ISTRAT,IOPT,IONED,IWK,NDV,NCON,NCOLA,NGT,NEQCON)     ADS01630
      NRB=MIN0(NCON,NCOLA)                                              ADS01640
      IF(NRB.LT.(2*NDV)) NRB=NRB+NDV                                    ADS01650
      IWK(67)=NRB                                                       ADS01660
C     IF DIMENSIONED SIZE IS TOO SMALL, PRINT WARNING COMENT.           ADS01670
      IF (IWK(150).LE.NRWK.AND.IWK(170).LE.NRIWK) GO TO 40              ADS01680
      WRITE (6,330)                                                     ADS01690
      WRITE (6,340) NRWK,IWK(150),NRIWK,IWK(170)                        ADS01700
40    CONTINUE                                                          ADS01710
      NSCFG=IWK(102)                                                    ADS01720
      NSCFX=IWK(105)                                                    ADS01730
C     IF INFO = -2.  INITIALIZATION IS COMPLETE.  SET INFO=-1 AND EXIT. ADS01740
      IF (INFO.EQ.-2) INFO=-1                                           ADS01750
      IF (INFO.LT.0) RETURN                                             ADS01760
C     INFO.NE.-2.  PROCEED.                                             ADS01770
50    CONTINUE                                                          ADS01780
C     IF DIMENSIONED SIZE IS STILL TOO SMALL, STOP THE OPTIMIZATION.    ADS01790
      IF (IWK(150).LE.NRWK.AND.IWK(170).LE.NRIWK) GO TO 60              ADS01800
      WRITE (6,350)                                                     ADS01810
      WRITE (6,340) NRWK,IWK(150),NRIWK,IWK(170)                        ADS01820
      INFO=0                                                            ADS01830
      IMAT=1                                                            ADS01840
      IWK(34)=IMAT                                                      ADS01850
      RETURN                                                            ADS01860
60    CONTINUE                                                          ADS01870
C     ----------------------------------------------------------------- ADS01880
C     INITIALIZATION FOR SCALE FACTORS.                                 ADS01890
C     ----------------------------------------------------------------- ADS01900
C     REPLACE THE SCALE FACTORS BY OVER-RIDEN VALUES.  IF NOT           ADS01910
C     OVER-RIDEN, THEN SET TO SCALE FACTORS TO 1.0                      ADS01920
      SCFO=WK(32)                                                       ADS01930
      ZRO=WK(37)                                                        ADS01940
      NSCFG=IWK(102)                                                    ADS01950
      NSVLB=IWK(103)                                                    ADS01960
      NSVUB=IWK(104)                                                    ADS01970
      NSCFX=IWK(105)                                                    ADS01980
      CALL ADS006 (NDV,NCON,VLB,VUB,ZRO,WK(NSCFX),WK(NSVLB),WK(NSVUB),SCADS01990
     1FO,WK(NSCFG))                                                     ADS02000
      WK(32)=SCFO                                                       ADS02010
      IGOTO=1                                                           ADS02020
      IWK(31)=IGOTO                                                     ADS02030
      IMAT=2                                                            ADS02040
      IWK(34)=IMAT                                                      ADS02050
C     INITIAL FUNCTION EVALUATION.                                      ADS02060
      INFO=1                                                            ADS02070
      RETURN                                                            ADS02080
C                                                                       ADS02090
70    CONTINUE                                                          ADS02100
C     UPDATE PRINT CONTROL VALUES IN WORK ARRAYS.                       ADS02110
      IWK(38)=IPRNT1                                                    ADS02120
      IWK(39)=IPRNT2                                                    ADS02130
      IWK(40)=IPRNT3                                                    ADS02140
C     ------------------------------------------------------------------ADS02150
C     GET SCALARS FROM WORK ARRAYS.                                     ADS02160
C     ------------------------------------------------------------------ADS02170
C     REAL.                                                             ADS02180
      CT=WK(3)                                                          ADS02190
      CTL=WK(4)                                                         ADS02200
      CTLMIN=WK(5)                                                      ADS02210
      CTMIN=WK(6)                                                       ADS02220
      DABOBJ=WK(8)                                                      ADS02230
      DABOBM=WK(9)                                                      ADS02240
      DABSTR=WK(10)                                                     ADS02250
      DELOBJ=WK(12)                                                     ADS02260
      DELOBM=WK(13)                                                     ADS02270
      EPSPEN=WK(19)                                                     ADS02280
      FDCH=WK(21)                                                       ADS02290
      FDCHM=WK(22)                                                      ADS02300
      OBJS=WK(78)                                                       ADS02310
      PENALT=WK(82)                                                     ADS02320
      PMLT=WK(38)                                                       ADS02330
      RP=WK(27)                                                         ADS02340
      SCFO=WK(32)                                                       ADS02350
      ZRO=WK(37)                                                        ADS02360
C     INTEGER.                                                          ADS02370
      IFCALL=IWK(28)                                                    ADS02380
      IGCALL=IWK(29)                                                    ADS02390
      IGOTO=IWK(31)                                                     ADS02400
      IPRNT1=IWK(38)                                                    ADS02410
      IPRNT2=IWK(39)                                                    ADS02420
      IPRNT3=IWK(40)                                                    ADS02430
      ISCAL=IWK(2)                                                      ADS02440
      ISFLAG=IWK(41)                                                    ADS02450
      JSCAL=IWK(50)                                                     ADS02460
      LGOTOZ=IWK(55)                                                    ADS02470
      MGOTOZ=IWK(57)                                                    ADS02480
      NAC=IWK(58)                                                       ADS02490
      NGOTOZ=IWK(64)                                                    ADS02500
      NVC=IWK(68)                                                       ADS02510
      NXFD=IWK(69)                                                      ADS02520
      NEQCON=IWK(71)                                                    ADS02530
C     ARRAY STARTING LOCATIONS.                                         ADS02540
      NSCFG=IWK(102)                                                    ADS02550
      NSVLB=IWK(103)                                                    ADS02560
      NSVUB=IWK(104)                                                    ADS02570
      NSCFX=IWK(105)                                                    ADS02580
      NXTMP=IWK(106)                                                    ADS02590
      NXTMP3=IWK(108)                                                   ADS02600
      NS=IWK(109)                                                       ADS02610
      NGL=IWK(114)                                                      ADS02620
      NG1=IWK(115)                                                      ADS02630
      NLAMDA=IWK(123)                                                   ADS02640
      NPARAM=IWK(124)                                                   ADS02650
      NPENLT=IWK(128)                                                   ADS02660
C     DEFAULT VALUE FOR DELOBJ, DABOBJ, DELOBM, DABOBM.                 ADS02670
      IF(IGOTO.EQ.1) CALL ADS010 (ISTRAT,NDV,OBJ,X,DELOBJ,DABOBJ,       ADS02680
     1DELOBM,DABOBM,DABSTR,ZRO)                                         ADS02690
      WK(8)=DABOBJ                                                      ADS02700
      WK(9)=DABOBM                                                      ADS02710
      WK(10)=DABSTR                                                     ADS02720
      WK(12)=DELOBJ                                                     ADS02730
      WK(13)=DELOBM                                                     ADS02740
C     OUTPUT VALUES OF PARAMETERS AND STORAGE REQUIREMENTS.             ADS02750
      IF (IPGLOB.GE.2.AND.IGOTO.EQ.1) CALL ADS403 (WK,NRWK,IWK,NRIWK)   ADS02760
      GO TO (80,100,120,200,240,270), IGOTO                             ADS02770
80    CONTINUE                                                          ADS02780
      ISFLAG=-1                                                         ADS02790
      INFO=2                                                            ADS02800
      IF (NCON.EQ.0) GO TO 90                                           ADS02810
C     SELECT THE NGT ACTIVE CONSTRAINTS.                                ADS02820
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS02830
90    CONTINUE                                                          ADS02840
C     ------------------------------------------------------------------ADS02850
C     CALCULATE SCALE FACTOR AND SCALE EACH VARIABLE.                   ADS02860
C     ------------------------------------------------------------------ADS02870
C     SAVE OBJECTIVE, DESIGN VARIABLE AND CONSTRAINT VALUES.            ADS02880
      OBJS=OBJ                                                          ADS02890
      CALL ADS503 (NDV,X,WK(NXTMP))                                     ADS02900
      IF (NCON.GT.0) CALL ADS503 (NCON,G,WK(NGL))                       ADS02910
C     GRADIENTS.  IF IGRAD.GT.0, GRADIENTS ARE USER-SUPPLIED.           ADS02920
      IGOTO=3                                                           ADS02930
      IF (IGRAD.GT.0.AND.INFO.EQ.2) IGCALL=IGCALL+1                     ADS02940
      IMAT=3                                                            ADS02950
      IF (IGRAD.GT.0) GO TO 320                                         ADS02960
C     IGRAD = 0.  FINITE DIFFERENCE GRADIENTS.                          ADS02970
      MGOTOZ=0                                                          ADS02980
100   CONTINUE                                                          ADS02990
      CALL ADS216 (MGOTOZ,OBJ,NDV,X,WK(NSVLB),WK(NSVUB),NCON,NGT,IC,G,WKADS03000
     1(NGL),FDCH,FDCHM,DF,A,NRA,WK(NXTMP),ZRO,WK(NXTMP3),NXFD,OBJS)     ADS03010
C     FINISHED?                                                         ADS03020
      INFO=1                                                            ADS03030
      IF (MGOTOZ.EQ.0) GO TO 110                                        ADS03040
C     NO.  EVALUATE OBJECTIVE AND CONSTRAINT FUNCTIONS AND CALL AGAIN.  ADS03050
      IGOTO=2                                                           ADS03060
      IFCALL=IFCALL+1                                                   ADS03070
      IMAT=4                                                            ADS03080
      GO TO 320                                                         ADS03090
110   CONTINUE                                                          ADS03100
      INFO=2                                                            ADS03110
120   CONTINUE                                                          ADS03120
C     GRADIENT CALCULATIONS COMPLETE. STORE NOMINAL CONSTRAINT VALUES   ADS03130
C     BACK IN OBJ, X AND G.                                             ADS03140
      OBJ=OBJS                                                          ADS03150
      CALL ADS503 (NDV,WK(NXTMP),X)                                     ADS03160
      IF (NCON.GT.0) CALL ADS503 (NCON,WK(NGL),G)                       ADS03170
C     (IGRAD ISCAL ISFLAG) = (JSCAL)   ,   *; 0 OR 1.                   ADS03180
C     (* 1 -1) = 1 ; SCALE FACTOR, SCALE X, G, OBJ, DF AND A.           ADS03190
C     (* 0 -1) = 7 ; SCALE X, G, OBJ, DF AND A.                         ADS03200
      IF (ISCAL.EQ.1) JSCAL=1                                           ADS03210
      IF (ISCAL.EQ.0) JSCAL=6                                           ADS03220
      CALL ADS007 (JSCAL,INFO,ISTRAT,NDV,NCON,NGT,IC,OBJ,X,G,DF,A,NRA,NCADS03230
     1OLA,WK(NSCFX),SCFO,WK(NSCFG),VLB,VUB,IPGLOB,WK(NXTMP),WK(NGL),ZRO,ADS03240
     2WK(NSVLB),WK(NSVUB))                                              ADS03250
      WK(32)=SCFO                                                       ADS03260
      WK(8)=WK(8)*SCFO                                                  ADS03270
      WK(9)=WK(9)*SCFO                                                  ADS03280
      WK(10)=WK(10)*SCFO                                                ADS03290
C     ------------------------------------------------------------------ADS03300
C     SCALE FACTORS ARE GAINED. SET ISFLAG 1.                           ADS03310
C     ------------------------------------------------------------------ADS03320
      ISFLAG=1                                                          ADS03330
C     REORDER IC AND A MATRIX, IF NECESSARY.                            ADS03340
      IF (ISTRAT.GE.1.AND.ISTRAT.LE.5) GO TO 130                        ADS03350
      IF (NCON.EQ.0) GO TO 130                                          ADS03360
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,CTLMIN,NAC,NVC)  ADS03370
130   CONTINUE                                                          ADS03380
      INFO=0                                                            ADS03390
      IF (ISTRAT.EQ.0) GO TO 150                                        ADS03400
140   CONTINUE                                                          ADS03410
C     ------------------------------------------------------------------ADS03420
C     STRATEGY LEVEL.                                                   ADS03430
C     ------------------------------------------------------------------ADS03440
C --- PENALTY PARAMETERS FOR EQUALITY CONSTRAINTS IF IOPT=4 OR 5.       ADS03450
      IF(NEQCON.EQ.0) GO TO 149                                         ADS03460
      IF(IOPT.NE.4.AND.IOPT.NE.5) GO TO 149                             ADS03470
      IF(INFO.GT.0) GO TO 143                                           ADS03480
      DO 144 ICALC=1,3                                                  ADS03490
144   CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS03500
     *ZRO,PENALT)                                                       ADS03510
      OBJ=OBJ+PENALT                                                    ADS03520
      GO TO 149                                                         ADS03530
143   ICALC=2                                                           ADS03540
      CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS03550
     *ZRO,PENALT)                                                       ADS03560
      OBJ=OBJ+PENALT                                                    ADS03570
      WK(82)=PENALT                                                     ADS03580
      IF(INFO.EQ.1) GO TO 149                                           ADS03590
      ICALC=3                                                           ADS03600
      CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS03610
     *ZRO,PENALT)                                                       ADS03620
149   CONTINUE                                                          ADS03630
      CALL ADS003 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,IGRAD,OBJ,X,WK(NSVLB)ADS03640
     1,WK(NSVUB),G,IDG,NGT,NAC,NVC,IC,DF,A,NRA,NCOLA,WK,IWK,ISFLAG)     ADS03650
      IF(NEQCON.EQ.0) GO TO 148                                         ADS03660
      IF(IOPT.NE.4.AND.IOPT.NE.5) GO TO 148                             ADS03670
      ICALC=4                                                           ADS03680
      CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS03690
     *ZRO,PENALT)                                                       ADS03700
      OBJ=OBJ-PENALT                                                    ADS03710
148   CONTINUE                                                          ADS03720
C     FINISHED?                                                         ADS03730
      IF (INFO.EQ.0) GO TO 310                                          ADS03740
C     NO. EVALUATE FUNCTIONS OR GRADIENTS AS REQUIRED.                  ADS03750
      LGOTOZ=1                                                          ADS03760
      GO TO 190                                                         ADS03770
150   IF (IOPT.EQ.0) GO TO 170                                          ADS03780
160   CONTINUE                                                          ADS03790
C     ------------------------------------------------------------------ADS03800
C     OPTIMIZATION LEVEL.                                               ADS03810
C     ------------------------------------------------------------------ADS03820
C --- ADD PENALTY FOR EQUALITY CONSTRAINTS IF IOPT=4 OR 5.              ADS03830
      IF(NEQCON.EQ.0) GO TO 169                                         ADS03840
      IF(IOPT.NE.4.AND.IOPT.NE.5) GO TO 169                             ADS03850
      IF(INFO.GT.0) GO TO 163                                           ADS03860
      DO 164 ICALC=1,3                                                  ADS03870
164   CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS03880
     *ZRO,PENALT)                                                       ADS03890
      OBJ=OBJ+PENALT                                                    ADS03900
      WK(82)=PENALT                                                     ADS03910
      GO TO 169                                                         ADS03920
163   ICALC=2                                                           ADS03930
      CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS03940
     *ZRO,PENALT)                                                       ADS03950
      OBJ=OBJ+PENALT                                                    ADS03960
      WK(82)=PENALT                                                     ADS03970
      IF(INFO.EQ.1) GO TO 169                                           ADS03980
      ICALC=3                                                           ADS03990
      CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS04000
     *ZRO,PENALT)                                                       ADS04010
169   CONTINUE                                                          ADS04020
      IWK(58)=NAC                                                       ADS04030
      IWK(68)=NVC                                                       ADS04040
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,WK(NSVLB),WK(NSADS04050
     1VUB),G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                          ADS04060
      IF(NEQCON.EQ.0) GO TO 168                                         ADS04070
      IF(IOPT.NE.4.AND.IOPT.NE.5) GO TO 168                             ADS04080
      ICALC=4                                                           ADS04090
      CALL ADS237(ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,WK(NPENLT),      ADS04100
     *ZRO,PENALT)                                                       ADS04110
      OBJ=OBJ-PENALT                                                    ADS04120
168   CONTINUE                                                          ADS04130
      NAC=IWK(58)                                                       ADS04140
      NVC=IWK(68)                                                       ADS04150
C     FINISHED?                                                         ADS04160
      IF (INFO.EQ.0) GO TO 310                                          ADS04170
C     NO. EVALUATE FUNCTIONS OR GRADIENTS AS REQUIRED.                  ADS04180
      LGOTOZ=2                                                          ADS04190
      GO TO 190                                                         ADS04200
170   IF (IONED.EQ.0) GO TO 310                                         ADS04210
180   CONTINUE                                                          ADS04220
C     ------------------------------------------------------------------ADS04230
C     ONE-DIMENSIONAL SEARCH LEVEL.                                     ADS04240
C     ------------------------------------------------------------------ADS04250
      ALPHA=WK(52)                                                      ADS04260
      ALPMAX=WK(57)                                                     ADS04270
      SLOPF=WK(85)                                                      ADS04280
      CALL ADS001 (INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,WK(NSVLB),          ADS04290
     1WK(NSVUB),G,IDG,NGT,IC,WK,IWK,SLOPF,ALPHA,ALPMAX)                 ADS04300
      WK(52)=ALPHA                                                      ADS04310
C     FINISHED?                                                         ADS04320
      IF (INFO.EQ.0) GO TO 310                                          ADS04330
C     NO.  EVALUATE FUNCTIONS.                                          ADS04340
      LGOTOZ=3                                                          ADS04350
190   CONTINUE                                                          ADS04360
C     ------------------------------------------------------------------ADS04370
C     FUNCTION OR GRADIENT EVALUATION.                                  ADS04380
C     ------------------------------------------------------------------ADS04390
      IF (INFO.GT.1) GO TO 220                                          ADS04400
C     IF INFO IS EQUAL TO 1, UN-SCALE X, G AND OBJ.                     ADS04410
      JSCAL=4                                                           ADS04420
      CALL ADS007 (JSCAL,INFO,ISTRAT,NDV,NCON,NGT,IC,OBJ,X,G,DF,A,NRA,NCADS04430
     1OLA,WK(NSCFX),SCFO,WK(NSCFG),VLB,VUB,IPGLOB,WK(NXTMP),WK(NGL),ZRO,ADS04440
     2WK(NSVLB),WK(NSVUB))                                              ADS04450
C     INFO = 1.  EVALUATE OBJECTIVE.                                    ADS04460
      IGOTO=4                                                           ADS04470
      IFCALL=IFCALL+1                                                   ADS04480
      IMAT=5                                                            ADS04490
      GO TO 320                                                         ADS04500
200   CONTINUE                                                          ADS04510
C     SCALE X, G AND OBJ.                                               ADS04520
      JSCAL=3                                                           ADS04530
      CALL ADS007 (JSCAL,INFO,ISTRAT,NDV,NCON,NGT,IC,OBJ,X,G,DF,A,NRA,NCADS04540
     1OLA,WK(NSCFX),SCFO,WK(NSCFG),VLB,VUB,IPGLOB,WK(NXTMP),WK(NGL),ZRO,ADS04550
     2WK(NSVLB),WK(NSVUB))                                              ADS04560
C     IF STRATEGY IS SUMT, EVALUATE PENALIZED FUNCTION.                 ADS04570
      IF (ISTRAT.LE.0.OR.ISTRAT.GE.6) GO TO 210                         ADS04580
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,WK(NLAMDA),WK(NPARADS04590
     1AM),NGT,IC)                                                       ADS04600
      OBJ=OBJ+PENALT                                                    ADS04610
      WK(82)=PENALT                                                     ADS04620

C       print*,'  > psudo penal=',obj,penalt

210   CONTINUE                                                          ADS04630
      GO TO (140,160,180), LGOTOZ                                       ADS04640
C                                                                       ADS04650
220   CONTINUE                                                          ADS04660
C     CALCULATE REQUIRED GRADIENTS.                                     ADS04670
C     SAVE OBJECTIVE, DESIGN VARIABLE AND CONSTRAINT VALUES.            ADS04680
      OBJS=OBJ                                                          ADS04690
      CALL ADS503 (NDV,X,WK(NXTMP))                                     ADS04700
      IF (NCON.GT.0) CALL ADS503 (NCON,G,WK(NGL))                       ADS04710
C     UN-SCALE X, G AND OBJ.                                            ADS04720
      JSCAL=4                                                           ADS04730
      CALL ADS007 (JSCAL,INFO,ISTRAT,NDV,NCON,NGT,IC,OBJ,X,G,DF,A,NRA,NCADS04740
     1OLA,WK(NSCFX),SCFO,WK(NSCFG),VLB,VUB,IPGLOB,WK(NXTMP),WK(NGL),ZRO,ADS04750
     2WK(NSVLB),WK(NSVUB))                                              ADS04760
C     SAVE X IN S.                                                      ADS04770
      CALL ADS503(NDV,X,WK(NS))                                         ADS04780
C     SAVE G IN G1.                                                     ADS04790
      IF (NCON.GT.0) CALL ADS503 (NCON,G,WK(NG1))                       ADS04800
C     GRADIENTS.  IF IGRAD.GT.0, GRADIENTS ARE USER-SUPPLIED.           ADS04810
      IGOTO=6                                                           ADS04820
      IF (IGRAD.GT.0.AND.INFO.EQ.2) IGCALL=IGCALL+1                     ADS04830
      IMAT=3                                                            ADS04840
      IF (IGRAD.GT.0) GO TO 320                                         ADS04850
C     IGRAD = 0.  FINITE DIFFERENCE GRADIENTS.                          ADS04860
      NGOTOZ=0                                                          ADS04870
240   CONTINUE                                                          ADS04880
      CALL ADS216 (NGOTOZ,OBJ,NDV,X,VLB,VUB,NCON,NGT,IC,G,WK(NG1),      ADS04890
     1FDCH,FDCHM,DF,A,NRA,WK(NS),ZRO,WK(NXTMP3),NXFD,OBJS)              ADS04900
C     FINISHED?                                                         ADS04910
      INFO=1                                                            ADS04920
      IF (NGOTOZ.EQ.0) GO TO 260                                        ADS04930
C     NO.  EVALUATE OBJECTIVE AND CONSTRAINT FUNCTIONS AND CALL AGAIN.  ADS04940
      IGOTO=5                                                           ADS04950
      IFCALL=IFCALL+1                                                   ADS04960
      IMAT=4                                                            ADS04970
      GO TO 320                                                         ADS04980
260   CONTINUE                                                          ADS04990
      INFO=2                                                            ADS05000
270   CONTINUE                                                          ADS05010
C     GRADIENT CALCULATIONS COMPLETE. STORE NOMINAL FUNCTION VALUES     ADS05020
C     BACK IN OBJ, X AND G.                                             ADS05030
      OBJ=OBJ*SCFO                                                      ADS05040
      CALL ADS503 (NDV,WK(NXTMP),X)                                     ADS05050
      IF (NCON.GT.0) CALL ADS503 (NCON,WK(NGL),G)                       ADS05060
C     ------------------------------------------------------------------ADS05070
C     CALCULATE SCALE FACTOR AND SCALE.                                 ADS05080
C     ------------------------------------------------------------------ADS05090
C     SCALE DF AND A.                                                   ADS05100
      JSCAL=2                                                           ADS05110
      CALL ADS007 (JSCAL,INFO,ISTRAT,NDV,NCON,NGT,IC,OBJ,X,G,DF,A,NRA,NCADS05120
     1OLA,WK(NSCFX),SCFO,WK(NSCFG),VLB,VUB,IPGLOB,WK(NXTMP),WK(NGL),ZRO,ADS05130
     2WK(NSVLB),WK(NSVUB))                                              ADS05140
      ISFLAG=1                                                          ADS05150
C     FOR F.D. OR M.F.D., RE-ORDER IC AND A MATRIX.                     ADS05160
      IF (ISTRAT.GE.1.AND.ISTRAT.LE.5) GO TO 305                        ADS05170
      IF(NCON.EQ.0) GO TO 305                                           ADS05180
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G                               ADS05190
     1,CT,CTMIN,CTL,CTLMIN,NAC,NVC)                                     ADS05200
305   CONTINUE                                                          ADS05210
C     GO BACK TO STRATEGY, OPTIMIZATION, OR ONE-DIMENSIONAL SEARCH.     ADS05220
      GO TO (140,160,180), LGOTOZ                                       ADS05230
310   CONTINUE                                                          ADS05240
C     ------------------------------------------------------------------ADS05250
C     OPTIMIZATION IS COMPLETE.                                         ADS05260
C     ------------------------------------------------------------------ADS05270
C     UN-SCALE X,OBJ AND G.                                             ADS05280
      JSCAL=4                                                           ADS05290
      CALL ADS007 (JSCAL,INFO,ISTRAT,NDV,NCON,NGT,IC,OBJ,X,G,DF,A,NRA,NCADS05300
     1OLA,WK(NSCFX),SCFO,WK(NSCFG),VLB,VUB,IPGLOB,WK(NXTMP),WK(NGL),ZRO,ADS05310
     2WK(NSVLB),WK(NSVUB))                                              ADS05320
C     PRINT FINAL RESULTS.                                              ADS05330
      IF (IPGLOB.GT.0) CALL ADS420 (NDV,NCON,X,VLB,VUB,G,OBJ,IFCALL,IGCAADS05340
     1LL,IGRAD)                                                         ADS05350
      IMAT=0                                                            ADS05360
320   CONTINUE                                                          ADS05370
C     ------------------------------------------------------------------ADS05380
C     PUT PARAMETERS IN WORK ARRAYS.                                    ADS05390
C     ------------------------------------------------------------------ADS05400
C     REAL.                                                             ADS05410
      WK(78)=OBJS                                                       ADS05420
      WK(82)=PENALT                                                     ADS05430
C     INTEGER.                                                          ADS05440
      IWK(28)=IFCALL                                                    ADS05450
      IWK(29)=IGCALL                                                    ADS05460
      IWK(31)=IGOTO                                                     ADS05470
      IWK(34)=IMAT                                                      ADS05480
      IWK(41)=ISFLAG                                                    ADS05490
      IWK(50)=JSCAL                                                     ADS05500
      IWK(55)=LGOTOZ                                                    ADS05510
      IWK(57)=MGOTOZ                                                    ADS05520
      IWK(58)=NAC                                                       ADS05530
      IWK(64)=NGOTOZ                                                    ADS05540
      IWK(68)=NVC                                                       ADS05550
      IWK(69)=NXFD                                                      ADS05560
      RETURN                                                            ADS05570
C     ------------------------------------------------------------------ADS05580
C     FORMATS                                                           ADS05590
C     ------------------------------------------------------------------ADS05600
330   FORMAT (/,5X,'<< WARNING >>')                                     ADS05610
340   FORMAT(7X,'DIMENSIONED SIZE OF WK OR IWK IS TOO SMALL !!',//,14X, ADS05620
     1'DIMENSIONED  REQUIRED',/,7X,'ARRAY     SIZE        SIZE',/,      ADS05630
     29X,'WK',4X,I6,6X,I6,/,8X,'IWK',4X,I6,6X,I6)                       ADS05640
350   FORMAT (/,5X,'<<DISASTER ERROR >>')                               ADS05650
360   FORMAT ('1',/////////,21X,                                        ADS05660
     1' AAAAA    DDDDDD     SSSSSS',/,21X,                              ADS05670
     2'A     A   D     D   S',/,21X,                                    ADS05680
     3'A     A   D     D   S',/,21X,                                    ADS05690
     4'AAAAAAA   D     D    SSSSS',/,21X,                               ADS05700
     5'A     A   D     D         S',/,21X,                              ADS05710
     6'A     A   D     D         S',/,21X,                              ADS05720
     7'A     A   DDDDDD    SSSSSS')                                     ADS05730
370   FORMAT (////,20X,'F O R T R A N   P R O G R A M',//,32X,'F O R',  ADS05740
     1//,9X,'A U T O M A T E D   D E S I G N   S Y N T H E S I S',//,   ADS05750
     225X,'V E R S I O N  2.01')                                        ADS05760
380   FORMAT (////////,5X,'CONTROL PARAMETERS',/,5X,'ISTRAT =',I5,4X,   ADS05770
     1'IOPT   =',I5,4X,'IONED  =',I5,4X,'IPRINT =',I5,/,5X,'IGRAD  =',  ADS05780
     2I5,4X,'NDV    =',I5,4X,'NCON   =',I5)                             ADS05790
390   FORMAT (//,5X,'* * * ISTRAT =',I3,' IOPT =',I3,                   ADS05800
     1' IS NOT A VALID COMBINATION')                                    ADS05810
400   FORMAT (//,5X,'* * * IOPT   =',I3,' IONED =',I3,                  ADS05820
     1' IS NOT A VALID COMBINATION')                                    ADS05830
410   FORMAT (/,5X,'* * * THE COMBINATION OF ISTRAT, IOPT AND IONED',   ADS05840
     1' IS NOT VALID',/,5X,'* * * OPTIMIZATION TERMINATED')             ADS05850
420   FORMAT (//,5X,'* * * ISTRAT =',I5,' IS NOT VALID')                ADS05860
430   FORMAT (//,5X,'* * * IOPT =',I5,' IS NOT VALID')                  ADS05870
440   FORMAT (//,5X,'* * * IONED =',I5,' IS NOT VALID')                 ADS05880
      END                                                               ADS05890
C --- SUBROUTINE ADS001 (INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,  ADS05900
C    1IDG,NAVC,IC,WK,IWK)                                               ADS05910
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),IDG(NCON),IC(NRIC),  ADS05920
C    * WK(NRWK),IWK(NRIWK)                                              ADS05930
C     ------------------------------------------------------------------ADS05940
C --- ROUTINE TO CONTROL ONE-DIMENSIONAL SEARCH.                        ADS05950
C     ------------------------------------------------------------------ADS05960
C     BY G. N. VANDERPLAATS                              JUNE 1982.     ADS05970
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS05980
C                                                                       ADS05990
C --- INPUT                                                             ADS06000
C                                                                       ADS06010
C --- OUTPUT                                                            ADS06020
      SUBROUTINE ADS001 (INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,  ADS06030
     1IDG,NAVC,IC,WK,IWK,SLOPF,ALPHA,ALPMAX)                            ADS06040
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),WK(1),IWK(1),    ADS06050
     1TEMP(10)                                                          ADS06060
C                                                                       ADS06070
C     ------------------------------------------------------------------ADS06080
C     SCALARS AND ARRAY LOCATIONS.                                      ADS06090
C     ------------------------------------------------------------------ADS06100
C                                                                       ADS06110
C     SCALARS.                                                          ADS06120
C     REAL.                                                             ADS06130
      SCFO=WK(32)                                                       ADS06140
      ALPBST=WK(51)                                                     ADS06150
      ALPL=WK(56)                                                       ADS06160
      ALPU=WK(60)                                                       ADS06170
      CTLMIN=WK(5)                                                      ADS06180
      CTMIN=WK(6)                                                       ADS06190
      DABALP=WK(7)                                                      ADS06200
      DELALP=WK(11)                                                     ADS06210
      EXTRAP=WK(20)                                                     ADS06220
      FBST=WK(64)                                                       ADS06230
      FL=WK(65)                                                         ADS06240
      FU=WK(66)                                                         ADS06250
      F1=WK(67)                                                         ADS06260
      F2=WK(68)                                                         ADS06270
      GBST=WK(69)                                                       ADS06280
      GMAX=WK(70)                                                       ADS06290
      GMCT=WK(72)                                                       ADS06300
      GMCTL=WK(73)                                                      ADS06310
      GMCTU=WK(74)                                                      ADS06320
      GMCT1=WK(75)                                                      ADS06330
      GMCT2=WK(76)                                                      ADS06340
      XL=WK(86)                                                         ADS06350
      XMULT=WK(36)                                                      ADS06360
      XU=WK(87)                                                         ADS06370
      X1=WK(88)                                                         ADS06380
      X2=WK(89)                                                         ADS06390
      ZRO=WK(37)                                                        ADS06400
C     INTEGER.                                                          ADS06410
      II=IWK(32)                                                        ADS06420
      II1=IWK(33)                                                       ADS06430
      IPRNT1=IWK(38)                                                    ADS06440
      ISLOPF=IWK(43)                                                    ADS06450
      ISLOPG=IWK(44)                                                    ADS06460
      IX=IWK(47)                                                        ADS06470
      JJ=IWK(49)                                                        ADS06480
      LGOTO=IWK(54)                                                     ADS06490
      MGOTO=IWK(56)                                                     ADS06500
      NITER=IWK(66)                                                     ADS06510
      ICAL=IWK(72)                                                      ADS06520
C                                                                       ADS06530
C     ARRAY STARTING LOCATIONS.                                         ADS06540
      NSCFG=IWK(102)                                                    ADS06550
      NSCAL=IWK(105)                                                    ADS06560
      NXTMP=IWK(106)                                                    ADS06570
      NS=IWK(109)                                                       ADS06580
      NXOLD=IWK(110)                                                    ADS06590
      NGU=IWK(113)                                                      ADS06600
      NGL=IWK(114)                                                      ADS06610
      NG1=IWK(115)                                                      ADS06620
      NG2=IWK(116)                                                      ADS06630
      NDG=IWK(117)                                                      ADS06640
      NXBST=IWK(125)                                                    ADS06650
C     ------------------------------------------------------------------ADS06660
C     WRITE OUTPUT INFORMATION AS REQUIRED.                             ADS06670
C     ------------------------------------------------------------------ADS06680
      IF (IPRNT1.EQ.0) GO TO 20                                         ADS06690
      IF (INFO.GT.0) GO TO 10                                           ADS06700
C     INITIAL WRITES.                                                   ADS06710
      WRITE (6,120)                                                     ADS06720
      IF (IONED.GE.1.AND.IONED.LE.4) WRITE (6,60)                       ADS06730
      IF (IONED.GE.5.AND.IONED.LE.8) WRITE (6,70)                       ADS06740
      IF (IONED.EQ.1.OR.IONED.EQ.5) WRITE (6,80) IONED                  ADS06750
      IF (IONED.EQ.2.OR.IONED.EQ.6) WRITE (6,90) IONED                  ADS06760
      IF (IONED.EQ.3.OR.IONED.EQ.7) WRITE (6,100) IONED                 ADS06770
      IF (IONED.EQ.4.OR.IONED.EQ.8) WRITE (6,110) IONED                 ADS06780
      Z=0.                                                              ADS06790
      WRITE (6,130) Z                                                   ADS06800
10    CONTINUE                                                          ADS06810
      OBJW=OBJ                                                          ADS06820
      IF(ISTRAT.EQ.0) OBJW=OBJ/SCFO                                     ADS06830
      WRITE (6,140) OBJW                                                ADS06840
      ICALC=1                                                           ADS06850
      IF (INFO.EQ.0) CALL ADS402 (ICALC,X,NDV,WK(NSCAL))                ADS06860
      ICALC=4                                                           ADS06870
      IF (IPRNT1.GT.1.AND.NCON.GT.0) CALL ADS402 (ICALC,G,NCON,WK(NSCFG)ADS06880
     1)                                                                 ADS06890
20    CONTINUE                                                          ADS06900
C                                                                       ADS06910
      IF (INFO.EQ.0) LGOTO=0                                            ADS06920
      IF (IONED.GT.4) GO TO 30                                          ADS06930
C     ------------------------------------------------------------------ADS06940
C     UNCONSTRAINED FUNCTION.                                           ADS06950
C     ------------------------------------------------------------------ADS06960
C                                                                       ADS06970
      CALL ADS116 (LGOTO,IONED,IPRNT1,NDV,NCON,OBJ,SLOPF,ISLOPF,ALPHA,ALADS06980
     1PMAX,XMULT,EXTRAP,ZRO,DELALP,DABALP,ALPL,ALPU,WK(NXOLD),SVLB,SVUB,ADS06990
     2WK(NS),X,G,WK(NGL),WK(NXBST),FBST,FL,FU,F1,F2,ALPBST,XL,XU,X1,X2,IADS07000
     3CAL,II,MGOTO,II1,IX,NITER,WK(NSCAL),WK(NSCFG),ISTRAT,SCFO)        ADS07010
      GO TO 40                                                          ADS07020
30    CONTINUE                                                          ADS07030
C     ------------------------------------------------------------------ADS07040
C     CONSTRAINED FUNCTION.                                             ADS07050
C     ------------------------------------------------------------------ADS07060
C                                                                       ADS07070
C     PUT TRANSFER INFORMATION IN TEMP.                                 ADS07080
      TEMP(1)=XL                                                        ADS07090
      TEMP(2)=XU                                                        ADS07100
      TEMP(3)=X1                                                        ADS07110
      TEMP(4)=X2                                                        ADS07120
      TEMP(5)=SCFO                                                      ADS07130
      CALL ADS115 (LGOTO,IONED,IPRNT1,NDV,NCON,OBJ,SLOPF,ISLOPF,ISLOPG,CADS07140
     1TMIN,CTLMIN,ALPHA,ALPMAX,XMULT,EXTRAP,ZRO,DELALP,DABALP,ALPL,ALPU,ADS07150
     2WK(NXOLD),SVLB,SVUB,WK(NS),X,G,WK(NDG),IDG,WK(NGL),WK(NG1),WK(NG2)ADS07160
     3,WK(NGU),WK(NXTMP),WK(NXBST),NAVC,IC,ALPBST,FBST,FL,FU,F1,F2,GBST,ADS07170
     4GMAX,GMCT,GMCTL,GMCTU,GMCT1,GMCT2,ICAL,II,MGOTO,JJ,II1,           ADS07180
     5IX,NITER,WK(NSCAL),WK(NSCFG),ISTRAT,TEMP)                         ADS07190
C     GET INFORMATION BACK FROM TEMP.                                   ADS07200
      XL=TEMP(1)                                                        ADS07210
      XU=TEMP(2)                                                        ADS07220
      X1=TEMP(3)                                                        ADS07230
      X2=TEMP(4)                                                        ADS07240
      SCFO=TEMP(5)                                                      ADS07250
40    CONTINUE                                                          ADS07260
C     SAVE SCALARS IN WORK ARRAYS.                                      ADS07270
C     REAL.                                                             ADS07280
      WK(51)=ALPBST                                                     ADS07290
      WK(56)=ALPL                                                       ADS07300
      WK(57)=ALPMAX                                                     ADS07310
      WK(60)=ALPU                                                       ADS07320
      WK(64)=FBST                                                       ADS07330
      WK(65)=FL                                                         ADS07340
      WK(66)=FU                                                         ADS07350
      WK(67)=F1                                                         ADS07360
      WK(68)=F2                                                         ADS07370
      WK(69)=GBST                                                       ADS07380
      WK(70)=GMAX                                                       ADS07390
      WK(72)=GMCT                                                       ADS07400
      WK(73)=GMCTL                                                      ADS07410
      WK(74)=GMCTU                                                      ADS07420
      WK(75)=GMCT1                                                      ADS07430
      WK(76)=GMCT2                                                      ADS07440
      WK(86)=XL                                                         ADS07450
      WK(87)=XU                                                         ADS07460
      WK(88)=X1                                                         ADS07470
      WK(89)=X2                                                         ADS07480
C     INTEGER.                                                          ADS07490
      IWK(32)=II                                                        ADS07500
      IWK(33)=II1                                                       ADS07510
      IWK(47)=IX                                                        ADS07520
      IWK(49)=JJ                                                        ADS07530
      IWK(54)=LGOTO                                                     ADS07540
      IWK(56)=MGOTO                                                     ADS07550
      IWK(66)=NITER                                                     ADS07560
      IWK(72)=ICAL                                                      ADS07570
      INFO=1                                                            ADS07580
      IF (LGOTO.EQ.0) INFO=0                                            ADS07590
      IF (IPRNT1.EQ.0) RETURN                                           ADS07600
      WRITE (6,130) ALPHA                                               ADS07610
      OBJW=OBJ                                                          ADS07620
      IF(ISTRAT.EQ.0) OBJW=OBJ/SCFO                                     ADS07630
      IF(INFO.EQ.0) WRITE (6,140) OBJW                                  ADS07640
      ICALC=1                                                           ADS07650
      CALL ADS402 (ICALC,X,NDV,WK(NSCAL))                               ADS07660
      IF (IPRNT1.LT.2) GO TO 50                                         ADS07670
      ICALC=4                                                           ADS07680
      IF (INFO.EQ.0.AND.NCON.GT.0) CALL ADS402 (ICALC,G,NCON,WK(NSCFG)) ADS07690
50    IF (INFO.EQ.0) WRITE (6,150)                                      ADS07700
      RETURN                                                            ADS07710
C     ------------------------------------------------------------------ADS07720
C     FORMATS.                                                          ADS07730
C     ------------------------------------------------------------------ADS07740
60    FORMAT (/,5X,'UNCONSTRAINED FUNCTION')                            ADS07750
70    FORMAT (/,5X,'CONSTRAINED FUNCTION')                              ADS07760
80    FORMAT (5X,'IONED =',I2,'  GOLDEN SECTION METHOD')                ADS07770
90    FORMAT (5X,'IONED =',I2,'  GOLDEN SECTION + POLYNOMIAL REFINEMENT'ADS07780
     1)                                                                 ADS07790
100   FORMAT (5X,'IONED =',I2,'  POLYNOMIAL FIT AFTER GETTING BOUNDS')  ADS07800
110   FORMAT (5X,'IONED =',I2,                                          ADS07810
     1'  POLYNOMIAL FIT WITHOUT FIRST GETTING BOUNDS')                  ADS07820
120   FORMAT(//,5X,'* * * * * BEGIN ONE-DIMENSIONAL SEARCH * * * * *')  ADS07830
130   FORMAT (/,5X,'ALPHA =',1PE13.5)                                   ADS07840
140   FORMAT (/,5X,'OBJ =',1PE13.5)                                     ADS07850
150   FORMAT(/,5X,'* * * * * END OF ONE-DIMENSIONAL SEARCH * * * * *')  ADS07860
      END                                                               ADS07870
C --- SUBROUTINE ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,         ADS07880
C    1SVLB,SVUB,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                     ADS07890
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),IDG(NCON),IC(NCOLA), ADS07900
C    1 DF(NDV),A(NRA,NCOLA),WK(NRWK),IWK(NRIWK)                         ADS07910
C     ------------------------------------------------------------------ADS07920
C --- ROUTINE TO CONTROL OPTIMIZATION.                                  ADS07930
C     ------------------------------------------------------------------ADS07940
C     BY G. N. VANDERPLAATS                               JUNE 1982.    ADS07950
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS07960
C --- INPUT                                                             ADS07970
C --- OUTPUT                                                            ADS07980
C                                                                       ADS07990
      SUBROUTINE ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUBADS08000
     1,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                              ADS08010
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WKADS08020
     1(1),IWK(1),TEMP(10),ITEMP(10)                                     ADS08030
C                                                                       ADS08040
C     ------------------------------------------------------------------ADS08050
C     SCALAR PARAMETERS.                                                ADS08060
C     ------------------------------------------------------------------ADS08070
      ALPHA=WK(52)                                                      ADS08080
      ALPMAX=WK(57)                                                     ADS08090
      ALPSAV=WK(58)                                                     ADS08100
      BETA=WK(61)                                                       ADS08110
      CT=WK(3)                                                          ADS08120
      CTL=WK(4)                                                         ADS08130
      CTLMIN=WK(5)                                                      ADS08140
      CTMIN=WK(6)                                                       ADS08150
      DABOBJ=WK(8)                                                      ADS08160
      DABOBM=WK(9)                                                      ADS08170
      DELOBJ=WK(12)                                                     ADS08180
      DELOBM=WK(13)                                                     ADS08190
      DF2=WK(63)                                                        ADS08200
      DLOBJ1=WK(15)                                                     ADS08210
      DLOBJ2=WK(16)                                                     ADS08220
      DX1=WK(17)                                                        ADS08230
      DX2=WK(18)                                                        ADS08240
      GMAX1=WK(71)                                                      ADS08250
      GMX1=WK(90)                                                       ADS08260
      OBJ1=WK(80)                                                       ADS08270
      PHI=WK(83)                                                        ADS08280
      SCFO=WK(32)                                                       ADS08290
      SLOPF=WK(85)                                                      ADS08300
      STOL=WK(34)                                                       ADS08310
      THETAZ=WK(35)                                                     ADS08320
      ZRO=WK(37)                                                        ADS08330
C     INTEGER.                                                          ADS08340
      ICNDIR=IWK(1)                                                     ADS08350
      ICOUNT=IWK(22)                                                    ADS08360
      IDAB=IWK(23)                                                      ADS08370
      IDEL=IWK(25)                                                      ADS08380
      IDIVCT=IWK(27)                                                    ADS08390
      IFCALL=IWK(28)                                                    ADS08400
      IGCALL=IWK(29)                                                    ADS08410
      IPRNT2=IWK(39)                                                    ADS08420
      ITER=IWK(45)                                                      ADS08430
      ITMAX=IWK(3)                                                      ADS08440
      ITRMOP=IWK(4)                                                     ADS08450
      IVAMET=IWK(46)                                                    ADS08460
      KGOTO=IWK(52)                                                     ADS08470
      KOUNT=IWK(53)                                                     ADS08480
      NAC=IWK(58)                                                       ADS08490
      NACS=IWK(59)                                                      ADS08500
      NAVC=IWK(61)                                                      ADS08510
      NGRAD=IWK(65)                                                     ADS08520
      N2GOTO=IWK(70)                                                    ADS08530
      NRB=IWK(67)                                                       ADS08540
      NVC=IWK(68)                                                       ADS08550
C     ARRAY STARTING LOCATIONS.                                         ADS08560
      NSCFG=IWK(102)                                                    ADS08570
      NSCAL=IWK(105)                                                    ADS08580
      NXTMP1=IWK(106)                                                   ADS08590
      NS=IWK(109)                                                       ADS08600
      NXOLD=IWK(110)                                                    ADS08610
      NB=IWK(111)                                                       ADS08620
      NGL=IWK(114)                                                      ADS08630
      NDG=IWK(117)                                                      ADS08640
      NC=IWK(118)                                                       ADS08650
      NTHETA=IWK(119)                                                   ADS08660
      NH=IWK(120)                                                       ADS08670
      NDX=IWK(121)                                                      ADS08680
      NDFG=IWK(122)                                                     ADS08690
      NGTMP=IWK(112)                                                    ADS08700
      NICQP=IWK(153)                                                    ADS08710
      NISIDC=IWK(152)                                                   ADS08720
      NITMP2=IWK(155)                                                   ADS08730
      IF (INFO.EQ.0) KGOTO=0                                            ADS08740
C     TURN ON ISLOPF.                                                   ADS08750
      IWK(43)=1                                                         ADS08760
      IF(ISTRAT.EQ.8) IWK(43)=0                                         ADS08770
C     TURN OFF ISLOPG.                                                  ADS08780
      IWK(44)=0                                                         ADS08790
      GO TO (10,20,20,30,40), IOPT                                      ADS08800
10    CONTINUE                                                          ADS08810
C     ------------------------------------------------------------------ADS08820
C     IOPT = 1.  FLETCHER-REEVES UNCONSTRAINED MINIMIZATION.            ADS08830
C     ------------------------------------------------------------------ADS08840
      IVAMET=3                                                          ADS08850
      CALL ADS204 (KGOTO,INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,ITMAX,ICNDIR,DADS08860
     1F,WK(NS),ZRO,DLOBJ1,DLOBJ2,DX1,DX2,WK(NXOLD),SVLB,SVUB,DELOBJ,DABOADS08870
     2BJ,ITRMOP,WK(NXTMP1),IPRNT2,WK,IWK,G,WK(NSCAL),WK(NH),WK(NDX),WK(NADS08880
     3DFG),IVAMET,IC,IDG,IWK(NICQP),SCFO,IFCALL,IGCALL,ALPHA,ALPSAV,DF2,ADS08890
     4OBJ1,ICOUNT,IDAB,IDEL,ITER,N2GOTO,ALPMAX,WK(NSCFG))               ADS08900
      GO TO 50                                                          ADS08910
20    CONTINUE                                                          ADS08920
C     ------------------------------------------------------------------ADS08930
C     IOPT = 2 OR 3.  DFP OR BFGS METHOD.                               ADS08940
C     ------------------------------------------------------------------ADS08950
      IVAMET=1                                                          ADS08960
      IF (IOPT.EQ.3) IVAMET=2                                           ADS08970
      CALL ADS235 (KGOTO,INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,ITMAX,ICNDIR,DADS08980
     1F,WK(NS),ZRO,DLOBJ1,DLOBJ2,DX1,DX2,WK(NXOLD),SVLB,SVUB,DELOBJ,DABOADS08990
     2BJ,ITRMOP,WK(NXTMP1),IPRNT2,WK,IWK,G,WK(NSCAL),WK(NH),WK(NDX),WK(NADS09000
     3DFG),IVAMET,IC,IDG,IWK(NISIDC),SCFO,IFCALL,IGCALL,ALPHA,ALPSAV,DF2ADS09010
     4,OBJ1,ICOUNT,IDAB,IDEL,IDIVCT,ITER,N2GOTO,ALPMAX,WK(NSCFG))       ADS09020
      GO TO 50                                                          ADS09030
30    CONTINUE                                                          ADS09040
C     ------------------------------------------------------------------ADS09050
C     IOPT = 4.  METHOD OF FEASIBLE DIRECTIONS.                         ADS09060
C     ------------------------------------------------------------------ADS09070
C     TURN ON ISLOPG.                                                   ADS09080
      IWK(44)=1                                                         ADS09090
C     PUT INFORMATION IN TEMP FOR TRANSFER (TO REDUCE NUMBER OF PARAM)  ADS09100
      TEMP(1)=ALPSAV                                                    ADS09110
      TEMP(2)=BETA                                                      ADS09120
      TEMP(3)=DF2                                                       ADS09130
      TEMP(4)=OBJ1                                                      ADS09140
      TEMP(5)=PHI                                                       ADS09150
      TEMP(6)=DLOBJ1                                                    ADS09160
      TEMP(7)=DLOBJ2                                                    ADS09170
      TEMP(8)=DX1                                                       ADS09180
      TEMP(9)=DX2                                                       ADS09190
      CALL ADS205 (KGOTO,INFO,IONED,NDV,NCON,ITMAX,ICNDIR,NGT,NAC,NVC,IPADS09200
     1RNT2,ITRMOP,OBJ,ZRO,CT,CTMIN,CTL,CTLMIN,STOL,THETAZ,DELOBJ,DABOBJ,ADS09210
     2SCFO,SLOPF,ALPHA,ALPMAX,X,SVLB,SVUB,DF,WK(NS),G,IC,IDG,WK(NDG),WK(ADS09220
     3NTHETA),WK(NC),IWK(NISIDC),WK(NXOLD),WK(NXTMP1),WK(NGL),WK(NSCAL),ADS09230
     4IWK(NITMP2),WK(NSCFG),A,NRA,NCOLA,WK(NB),NRB,WK,IWK,              ADS09240
     5ICOUNT,IDAB,IDEL,ITER,N2GOTO,NACS,NAVC,NGRAD,ISTRAT,TEMP)         ADS09250
C     GET INFORMATION BACK FROM TEMP.                                   ADS09260
      ALPSAV=TEMP(1)                                                    ADS09270
      BETA=TEMP(2)                                                      ADS09280
      DF2=TEMP(3)                                                       ADS09290
      OBJ1=TEMP(4)                                                      ADS09300
      PHI=TEMP(5)                                                       ADS09310
      DLOBJ1=TEMP(6)                                                    ADS09320
      DLOBJ2=TEMP(7)                                                    ADS09330
      DX1=TEMP(8)                                                       ADS09340
      DX2=TEMP(9)                                                       ADS09350
      GO TO 50                                                          ADS09360
40    CONTINUE                                                          ADS09370
C     ------------------------------------------------------------------ADS09380
C     IOPT = 5.  MODIFIED METHOD OF FEASIBLE DIRECTIONS.                ADS09390
C     ------------------------------------------------------------------ADS09400
C     TURN ON ISLOPG.                                                   ADS09410
      IWK(44)=1                                                         ADS09420
C     PUT INFO IN TEMP TO SHORTEN PARAMETER LIST.                       ADS09430
      TEMP(1)=ALPSAV                                                    ADS09440
      TEMP(2)=BETA                                                      ADS09450
      TEMP(3)=DF2                                                       ADS09460
      TEMP(4)=OBJ1                                                      ADS09470
      TEMP(5)=PHI                                                       ADS09480
      TEMP(6)=DLOBJ1                                                    ADS09490
      TEMP(7)=DLOBJ2                                                    ADS09500
      TEMP(8)=DX1                                                       ADS09510
      TEMP(9)=DX2                                                       ADS09520
      TEMP(10)=GMAX1                                                    ADS09530
      ITEMP(1)=ICOUNT                                                   ADS09540
      ITEMP(2)=IDAB                                                     ADS09550
      ITEMP(3)=IDEL                                                     ADS09560
      ITEMP(4)=ITER                                                     ADS09570
      ITEMP(5)=KOUNT                                                    ADS09580
      ITEMP(6)=LGOTO                                                    ADS09590
      ITEMP(7)=NACS                                                     ADS09600
      ITEMP(8)=NAVC                                                     ADS09610
      ITEMP(9)=NGRAD                                                    ADS09620
      ITEMP(10)=ISTRAT                                                  ADS09630
      CALL ADS231 (KGOTO,INFO,IONED,NDV,NCON,ITMAX,ICNDIR,NGT,NAC,NVC,IPADS09640
     1RNT2,ITRMOP,OBJ,ZRO,CT,CTMIN,CTL,CTLMIN,STOL,THETAZ,DELOBJ,DABOBJ,ADS09650
     2SCFO,SLOPF,ALPHA,ALPMAX,X,SVLB,SVUB,DF,WK(NS),G,IC,IDG,WK(NDG),WK(ADS09660
     3NTHETA),WK(NC),IWK(NISIDC),WK(NXOLD),WK(NXTMP1),WK(NGL),WK(NSCAL),ADS09670
     4IWK(NITMP2),WK(NSCFG),A,NRA,NCOLA,WK(NB),NRB,WK,IWK,WK(NGTMP),    ADS09680
     5ITEMP,TEMP,GMX1)                                                  ADS09690
C     GET INFO BACK FROM TEMP.                                          ADS09700
      ALPSAV=TEMP(1)                                                    ADS09710
      BETA=TEMP(2)                                                      ADS09720
      DF2=TEMP(3)                                                       ADS09730
      OBJ1=TEMP(4)                                                      ADS09740
      PHI=TEMP(5)                                                       ADS09750
      DLOBJ1=TEMP(6)                                                    ADS09760
      DLOBJ2=TEMP(7)                                                    ADS09770
      DX1=TEMP(8)                                                       ADS09780
      DX2=TEMP(9)                                                       ADS09790
      GMAX1=TEMP(10)                                                    ADS09800
      ICOUNT=ITEMP(1)                                                   ADS09810
      IDAB=ITEMP(2)                                                     ADS09820
      IDEL=ITEMP(3)                                                     ADS09830
      ITER=ITEMP(4)                                                     ADS09840
      KOUNT=ITEMP(5)                                                    ADS09850
      LGOTO=ITEMP(6)                                                    ADS09860
      NACS=ITEMP(7)                                                     ADS09870
      NAVC=ITEMP(8)                                                     ADS09880
      NGRAD=ITEMP(9)                                                    ADS09890
      ISTRAT=ITEMP(10)                                                  ADS09900
50    CONTINUE                                                          ADS09910
C     STORE CHANGED VALUES IN WORK ARRAYS.                              ADS09920
C     REAL.                                                             ADS09930
      WK(52)=ALPHA                                                      ADS09940
      WK(57)=ALPMAX                                                     ADS09950
      WK(58)=ALPSAV                                                     ADS09960
      WK(61)=BETA                                                       ADS09970
      WK(3)=CT                                                          ADS09980
      WK(4)=CTL                                                         ADS09990
      WK(63)=DF2                                                        ADS10000
      WK(71)=GMAX1                                                      ADS10010
      WK(80)=OBJ1                                                       ADS10020
      WK(83)=PHI                                                        ADS10030
      WK(85)=SLOPF                                                      ADS10040
      WK(90)=GMX1                                                       ADS10050
C     INTEGER.                                                          ADS10060
      IWK(22)=ICOUNT                                                    ADS10070
      IWK(23)=IDAB                                                      ADS10080
      IWK(25)=IDEL                                                      ADS10090
      IWK(27)=IDIVCT                                                    ADS10100
      IWK(28)=IFCALL                                                    ADS10110
      IWK(29)=IGCALL                                                    ADS10120
      IWK(45)=ITER                                                      ADS10130
      IWK(46)=IVAMET                                                    ADS10140
      IWK(52)=KGOTO                                                     ADS10150
      IWK(53)=KOUNT                                                     ADS10160
      IWK(58)=NAC                                                       ADS10170
      IWK(59)=NACS                                                      ADS10180
      IWK(61)=NAVC                                                      ADS10190
      IWK(65)=NGRAD                                                     ADS10200
      IWK(68)=NVC                                                       ADS10210
      IWK(70)=N2GOTO                                                    ADS10220
      IF (KGOTO.GT.0) GO TO 60                                          ADS10230
      IF (ISTRAT.LT.1.OR.ISTRAT.GT.5) GO TO 60                          ADS10240
      WK(12)=DMAX1(0.6*WK(12),DELOBM)                                   ADS10250
      WK(8)=DMAX1(0.6*WK(8),DABOBM)                                     ADS10260
60    CONTINUE                                                          ADS10270
      IF (KGOTO.EQ.0) INFO=0                                            ADS10280
      IF(INFO.EQ.0) RETURN                                              ADS10290
      IWK(43)=0                                                         ADS10300
      IWK(44)=0                                                         ADS10310
      RETURN                                                            ADS10320
      END                                                               ADS10330
C --- SUBROUTINE ADS003 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,IGRAD,         ADS10340
C    1 OBJ,X,SVLB,SVUB,G,IDG,NGT,NAC,NVC,IC,DF,A,NRA,NCOLA,WK,IWK,      ADS10350
C    2 ISFLAG)                                                          ADS10360
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),A(NRA,NCOLA),WK(NRWK)ADS10370
C    1,IWK(NRIWK)                                                       ADS10380
C     ------------------------------------------------------------------ADS10390
C --- ROUTINE TO CONTROL OPTIMIZATION STRATEGY.                         ADS10400
C     ------------------------------------------------------------------ADS10410
C     BY G. N. VANDERPLAATSJUNE 1982.                                   ADS10420
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS10430
C --- INPUT                                                             ADS10440
C                                                                       ADS10450
C --- OUTPUT                                                            ADS10460
      SUBROUTINE ADS003 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,IGRAD,OBJ,X,SVLADS10470
     1B,SVUB,G,IDG,NGT,NAC,NVC,IC,DF,A,NRA,NCOLA,WK,IWK,ISFLAG)         ADS10480
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WKADS10490
     1(1),IWK(1),TEMP2(10)                                              ADS10500
C                                                                       ADS10510
C     REAL.                                                             ADS10520
      ALAMDZ=WK(1)                                                      ADS10530
      ALPHA3=WK(53)                                                     ADS10540
      ALPHB3=WK(55)                                                     ADS10550
      ALPSV3=WK(59)                                                     ADS10560
      BETAMC=WK(2)                                                      ADS10570
      CT=WK(3)                                                          ADS10580
      CTL=WK(4)                                                         ADS10590
      CTLMIN=WK(5)                                                      ADS10600
      CTMIN=WK(6)                                                       ADS10610
      DABSTR=WK(10)                                                     ADS10620
      DELSTR=WK(14)                                                     ADS10630
      EPSPEN=WK(19)                                                     ADS10640
      GMAX=WK(70)                                                       ADS10650
      GMULTZ=WK(23)                                                     ADS10660
      OBJR=WK(77)                                                       ADS10670
      OBJZ=WK(79)                                                       ADS10680
      OBJ1B=WK(81)                                                      ADS10690
      PENALT=WK(82)                                                     ADS10700
      PSAIZ=WK(24)                                                      ADS10710
      RMULT=WK(25)                                                      ADS10720
      RMVLMT=WK(84)                                                     ADS10730
      RMVLMZ=WK(26)                                                     ADS10740
      RP=WK(27)                                                         ADS10750
      RPMAX=WK(28)                                                      ADS10760
      RPMULT=WK(29)                                                     ADS10770
      RPPMIN=WK(30)                                                     ADS10780
      RPPRIM=WK(31)                                                     ADS10790
      SCFO=WK(32)                                                       ADS10800
      SCLMIN=WK(33)                                                     ADS10810
      STOL=WK(34)                                                       ADS10820
      SLOPF=WK(85)                                                      ADS10830
      ZRO=WK(37)                                                        ADS10840
C     INTEGER.                                                          ADS10850
      IDAB3=IWK(25)                                                     ADS10860
      IDEL3=IWK(26)                                                     ADS10870
      IFCALL=IWK(28)                                                    ADS10880
      IGCALL=IWK(29)                                                    ADS10890
      IGMXC=IWK(30)                                                     ADS10900
      IMVCH=IWK(35)                                                     ADS10910
      IMVCM=IWK(36)                                                     ADS10920
      INFO2C=IWK(37)                                                    ADS10930
      IPRNT3=IWK(40)                                                    ADS10940
      ISLOP=IWK(42)                                                     ADS10950
      ITRMST=IWK(5)                                                     ADS10960
      JGOTO=IWK(48)                                                     ADS10970
      JONED=IWK(6)                                                      ADS10980
      JTER=IWK(51)                                                      ADS10990
      JTMAX=IWK(7)                                                      ADS11000
      NACT=IWK(60)                                                      ADS11010
      NAVC=IWK(61)                                                      ADS11020
      NGOTO=IWK(63)                                                     ADS11030
      NRB=IWK(67)                                                       ADS11040
      NRAA=NDV+1                                                        ADS11050
C     ARRAY STARTING LOCATIONS.                                         ADS11060
      NSCFG=IWK(102)                                                    ADS11070
      NSCAL=IWK(105)                                                    ADS11080
      NXTMP=IWK(106)                                                    ADS11090
      NXTMP2=IWK(107)                                                   ADS11100
      NXTMP3=IWK(108)                                                   ADS11110
      NS=IWK(109)                                                       ADS11120
      NXOLD=IWK(110)                                                    ADS11130
      NB=IWK(111)                                                       ADS11140
      NGL=IWK(114)                                                      ADS11150
      NC=IWK(118)                                                       ADS11160
      NTHETA=IWK(119)                                                   ADS11170
      NH=IWK(120)                                                       ADS11180
      NDX=IWK(121)                                                      ADS11190
      NDFG=IWK(122)                                                     ADS11200
      NLAMDA=IWK(123)                                                   ADS11210
      NPARAM=IWK(124)                                                   ADS11220
      NVLBTP=IWK(126)                                                   ADS11230
      NVUBTP=IWK(127)                                                   ADS11240
      NISIDC=IWK(152)                                                   ADS11250
      NICQP=IWK(153)                                                    ADS11260
      NITMP=IWK(154)                                                    ADS11270
      GO TO (10,20,30,40,50,60,70,80,85), ISTRAT                        ADS11280
C     ------------------------------------------------------------------ADS11290
C     ISTRAT = 1.  EXTERIOR PENALTY FUNCTION METHOD.                    ADS11300
C     ------------------------------------------------------------------ADS11310
10    CONTINUE                                                          ADS11320
      CALL ADS301 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,ADS11330
     1NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,WK(NXTMPADS11340
     2),CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,WK(NSCAL),WK(NLAMDA),WK(ADS11350
     3NPARAM),SCFO,ISFLAG,IFCALL,IGCALL,OBJZ,OBJ1B,PENALT,IDAB,IDEL,INFOADS11360
     42C,JTER,IWK(NITMP),IWK(NISIDC),WK(NSCFG))                         ADS11370
      GO TO 90                                                          ADS11380
20    CONTINUE                                                          ADS11390
C     ------------------------------------------------------------------ADS11400
C     ISTRAT = 2.  LINEAR EXTENDED INTERIOR PENALTY FUNCTION METHOD.    ADS11410
C     ------------------------------------------------------------------ADS11420
      CALL ADS303 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,ADS11430
     1NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,WK(NXTMPADS11440
     2),CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPMADS11450
     3ULT,WK(NSCAL),WK(NLAMDA),WK(NPARAM),SCFO,ISFLAG,IFCALL,IGCALL,OBJZADS11460
     4,OBJ1B,PENALT,IDAB,IDEL,INFO2C,JTER,IWK(NITMP),IWK(NISIDC),WK(NSCFADS11470
     5G))                                                               ADS11480
      GO TO 90                                                          ADS11490
30    CONTINUE                                                          ADS11500
C     ------------------------------------------------------------------ADS11510
C     ISTRAT = 3.  QUADRATIC EXTENDED INTERIOR PENALTY FUNCTION METHOD. ADS11520
C     ------------------------------------------------------------------ADS11530
      CALL ADS304 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,ADS11540
     1NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,WK(NXTMPADS11550
     2),CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPMADS11560
     3ULT,WK(NSCAL),WK(NLAMDA),WK(NPARAM),SCFO,ISFLAG,IFCALL,IGCALL,OBJZADS11570
     4,OBJ1B,PENALT,IDAB,IDEL,INFO2C,JTER,IWK(NITMP),IWK(NISIDC),WK(NSCFADS11580
     5G))                                                               ADS11590
      GO TO 90                                                          ADS11600
40    CONTINUE                                                          ADS11610
C     ------------------------------------------------------------------ADS11620
C     ISTRAT = 4.  CUBIC EXTENDED INTERIOR PENALTY FUNCTION METHOD.     ADS11630
C     ------------------------------------------------------------------ADS11640
      CALL ADS305 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,ADS11650
     1NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,WK(NXTMPADS11660
     2),CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPMADS11670
     3ULT,WK(NSCAL),WK(NLAMDA),WK(NPARAM),SCFO,ISFLAG,IFCALL,IGCALL,OBJZADS11680
     4,OBJ1B,PENALT,IDAB,IDEL,INFO2C,JTER,IWK(NITMP),IWK(NISIDC),WK(NSCFADS11690
     5G))                                                               ADS11700
      GO TO 90                                                          ADS11710
50    CONTINUE                                                          ADS11720
C     ------------------------------------------------------------------ADS11730
C     ISTRAT = 5.  AUGMENTED LAGRANGE MULTIPLIER METHOD.                ADS11740
C     ------------------------------------------------------------------ADS11750
      CALL ADS306 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,ADS11760
     1NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,WK(NXTMPADS11770
     2),CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPMADS11780
     3ULT,WK(NSCAL),ALAMDZ,WK(NLAMDA),WK(NPARAM),SCFO,ISFLAG,IFCALL,IGCAADS11790
     4LL,CK,OBJZ,OBJ1B,PENALT,IDAB,IDEL,INFO2C,JTER,IWK(NITMP),IWK(NISIDADS11800
     5C),WK(NSCFG))                                                     ADS11810
      GO TO 90                                                          ADS11820
60    CONTINUE                                                          ADS11830
C     ------------------------------------------------------------------ADS11840
C     ISTRAT = 6.  SEQUENTIAL LINEAR PROGRAMMING.                       ADS11850
C     ------------------------------------------------------------------ADS11860
      CALL ADS320 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,A,NCOLADS11870
     1A,NRA,X,IDG,WK,IWK,IWK(NICQP),SVLB,SVUB,WK(NS),WK(NLAMDA),WK(NDFG)ADS11880
     2,WK(NDX),ZRO,IPRNT3,RMVLMZ,WK(NSCAL),WK(NXTMP),SCLMIN,DELSTR,DABSTADS11890
     3R,ITRMST,JTMAX,SCFO,IWK(NITMP),CT,CTMIN,CTL,CTLMIN,NAC,NVC,WK(NXTMADS11900
     4P2),WK(NXTMP3),OBJR,OBJ1B,RMVLMT,IDAB,IDEL,IMVCM,JTER,IWK(NISIDC),ADS11910
     5WK(NSCFG))                                                        ADS11920
      GO TO 90                                                          ADS11930
70    CONTINUE                                                          ADS11940
C     ------------------------------------------------------------------ADS11950
C     ISTRAT = 7.  METHOD OF CENTERS.                                   ADS11960
C     ------------------------------------------------------------------ADS11970
      CALL ADS330 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,A,NCOLADS11980
     1A,NRA,X,IDG,WK,IWK,IWK(NICQP),SVLB,SVUB,WK(NS),WK(NLAMDA),WK(NPARAADS11990
     2M),WK(NVLBTP),WK(NVUBTP),ZRO,IPRNT3,RMVLMZ,WK(NSCAL),WK(NXTMP),SCLADS12000
     3MIN,DELSTR,DABSTR,ITRMST,JTMAX,BETAMC,SCFO,IWK(NITMP),CT,CTMIN,CTLADS12010
     4,CTLMIN,NAC,NVC,WK(NDX),WK(NXTMP2),WK(NXTMP3),GMAX,OBJR,OBJ1B,IDABADS12020
     5,IDEL,IGMXC,IMVCH,JTER,IWK(NISIDC),WK(NSCFG))                     ADS12030
      GO TO 90                                                          ADS12040
80    CONTINUE                                                          ADS12050
C     ------------------------------------------------------------------ADS12060
C     ISTRAT = 8.  VARIABLE METRIC METHOD FOR CONSTRAINED OPTIMIZATION. ADS12070
C     ------------------------------------------------------------------ADS12080
C     PUT INFO IN TEMP2 TO REDUCE NUMBER OF CALLING PARAMETERS.         ADS12090
      TEMP2(1)=CTMIN                                                    ADS12100
      TEMP2(2)=CTL                                                      ADS12110
      TEMP2(3)=CTLMIN                                                   ADS12120
      TEMP2(4)=GMULTZ                                                   ADS12130
      TEMP2(5)=ALPHA3                                                   ADS12140
      TEMP2(6)=ALPHB3                                                   ADS12150
      TEMP2(7)=ALPSV3                                                   ADS12160
      TEMP2(8)=OBJ1B                                                    ADS12170
      TEMP2(9)=SLOPF                                                    ADS12180
      TEMP2(10)=STOL                                                    ADS12190
      CALL ADS310 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,A,NCOLADS12200
     1A,NRA,WK(NH),X,IDG,WK,IWK,IWK(NICQP),SVLB,SVUB,WK(NS),WK(NLAMDA),WADS12210
     2K(NPARAM),WK(NXOLD),WK(NDFG),WK(NDX),ZRO,IPRNT3,PSAIZ,WK(NSCAL),WKADS12220
     3(NXTMP2),DELSTR,DABSTR,ITRMST,JTMAX,IWK(NISIDC),WK(NGL),WK(NXTMP),ADS12230
     4SCFO,IWK(NITMP),WK(NB),NRB,WK(NTHETA),CT,NAC,NVC,                 ADS12240
     5WK(NC),JONED,WK(NXTMP3),WK(NVLBTP),WK(NVUBTP),                    ADS12250
     6IDAB,IDEL,ISLOP,JTER,NGOTO,NACT,NAVC,TEMP2,WK(NSCFG))             ADS12260
C     GET INFO BACK FROM TEMP2.                                         ADS12270
      CTMIN=TEMP2(1)                                                    ADS12280
      CTL=TEMP2(2)                                                      ADS12290
      CTLMIN=TEMP2(3)                                                   ADS12300
      GMULTZ=TEMP2(4)                                                   ADS12310
      ALPHA3=TEMP2(5)                                                   ADS12320
      ALPHB3=TEMP2(6)                                                   ADS12330
      ALPSV3=TEMP2(7)                                                   ADS12340
      OBJ1B=TEMP2(8)                                                    ADS12350
      SLOPF=TEMP2(9)                                                    ADS12360
      STOL=TEMP2(10)                                                    ADS12370
      GO TO 90                                                          ADS12380
85    CONTINUE                                                          ADS12390
C     ------------------------------------------------------------------ADS12400
C     ISTRAT = 9.  SEQUENTIAL CONVEX PROGRAMMING.                       ADS12410
C     ------------------------------------------------------------------ADS12420
      CALL ADS325 (JGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,A,NCOLADS12430
     1A,NRA,X,IDG,WK,IWK,IWK(NICQP),SVLB,SVUB,WK(NS),WK(NLAMDA),WK(NDFG)ADS12440
     2,WK(NDX),ZRO,IPRNT3,WK(NSCAL),WK(NXTMP),SCLMIN,DELSTR,DABST       ADS12450
     3R,ITRMST,JTMAX,SCFO,IWK(NITMP),CT,CTMIN,CTL,CTLMIN,NAC,NVC,WK(NXTMADS12460
     4P2),WK(NXTMP3),OBJR,OBJ1B,IDAB,IDEL,JTER,IWK(NISIDC),             ADS12470
     5WK(NB),NRB,WK(NC),WK(NTHETA),WK(NGL),WK(NH),NRAA,RMVLMZ,RMVLMT,   ADS12480
     6WK(NSCFG))                                                        ADS12490
90    CONTINUE                                                          ADS12500
C     STORE PARAMETERS BACK IN WORK ARRAYS.                             ADS12510
C     REAL.                                                             ADS12520
      WK(53)=ALPHA3                                                     ADS12530
      WK(55)=ALPHB3                                                     ADS12540
      WK(59)=ALPSV3                                                     ADS12550
      WK(19)=EPSPEN                                                     ADS12560
      WK(70)=GMAX                                                       ADS12570
      WK(77)=OBJR                                                       ADS12580
      WK(79)=OBJZ                                                       ADS12590
      WK(81)=OBJ1B                                                      ADS12600
      WK(82)=PENALT                                                     ADS12610
      WK(84)=RMVLMT                                                     ADS12620
      WK(27)=RP                                                         ADS12630
      WK(31)=RPPRIM                                                     ADS12640
      WK(85)=SLOPF                                                      ADS12650
C     INTEGER.                                                          ADS12660
      IWK(25)=IDAB3                                                     ADS12670
      IWK(26)=IDEL3                                                     ADS12680
      IWK(28)=IFCALL                                                    ADS12690
      IWK(29)=IGCALL                                                    ADS12700
      IWK(30)=IGMXC                                                     ADS12710
      IWK(35)=IMVCH                                                     ADS12720
      IWK(37)=INFO2C                                                    ADS12730
      IWK(42)=ISLOP                                                     ADS12740
      IWK(48)=JGOTO                                                     ADS12750
      IWK(51)=JTER                                                      ADS12760
      IWK(60)=NACT                                                      ADS12770
      IWK(61)=NAVC                                                      ADS12780
      IWK(63)=NGOTO                                                     ADS12790
      IF (JGOTO.EQ.0) INFO=0                                            ADS12800
      RETURN                                                            ADS12810
      END                                                               ADS12820
C --- SUBROUTINE ADS004 (ISTRAT,IOPT,IONED,IWK,NDV,NCON,NCOLA,NGT,      ADS12830
C    *NEQCON)                                                           ADS12840
C     DIMENSION IWK(NRIWK)                                              ADS12850
C     ------------------------------------------------------------------ADS12860
C --- ROUTINE TO DEFINE WORK ARRAY STORAGE ALLOCATIONS.                 ADS12870
C     ------------------------------------------------------------------ADS12880
C     BY G. N. VANDERPLAATS                               JUNE 1982.    ADS12890
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS12900
C --- INPUT.                                                            ADS12910
C     ISTRAT - OPTIMIZATION STRATEGY IDENTIFIER.                        ADS12920
C     IOPT - OPTIMIZER NUMBER.                                          ADS12930
C     IONED - ONE-DIMENSIONAL SEARCH IDENTIFIER.                        ADS12940
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS12950
C     NCON - NUMBER OF CONSTRAINTS.                                     ADS12960
C     NCOLA - DIMENSIONED COLUMNS OF MATRIX A.                          ADS12970
C     NEQCON - NUMBER OF EQUALITY CONSTRAINTS.                          ADS12980
C                                                                       ADS12990
C --- OUTPUT.                                                           ADS13000
C     IWK(I), I=101,200 - STORAGE ALLOCATION POINTERS FOR WORK ARRAYS.  ADS13010
C                                                                       ADS13020
      SUBROUTINE ADS004 (ISTRAT,IOPT,IONED,IWK,NDV,NCON,NCOLA,NGT,      ADS13030
     $NEQCON)                                                           ADS13040
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION IWK(200)                                                ADS13050
C                                                                       ADS13060
C     ------------------------------------------------------------------ADS13070
C     REAL ARRAY STORAGE.                                               ADS13080
C     ------------------------------------------------------------------ADS13090
C                                                                       ADS13100
C     RSCALR, REAL SCALAR VARIABLES.                                    ADS13110
C     100 LOCATIONS PRESENTLY RESERVED.                                 ADS13120
C     IF METHOD OF CENTERS , NCON,NDV SHOULD BE ADDED 1.                ADS13130
      NRB=MIN0(NCON,NCOLA)                                              ADS13140
      IF(NRB.LT.(2*NDV)) NRB=NRB+NDV                                    ADS13150
      IWK(101)=1                                                        ADS13160
C     SCFG. SCALE FACTORS FOR CONSTRAINTS.                              ADS13170
      IWK(102)=101                                                      ADS13180
C     SVLB, SCALED VLB-VECTOR.                                          ADS13190
      IWK(103)=IWK(102)+NCON                                            ADS13200
C     SVUB, SCALED VUB-VECTOR.                                          ADS13210
      IWK(104)=IWK(103)+NDV+1                                           ADS13220
C     SCFX, SCALE FACTORS FOR DESIGN VARIABLES.                         ADS13230
      IWK(105)=IWK(104)+NDV+1                                           ADS13240
C     XTMP, TEMPORARY X-VECTOR.                                         ADS13250
      IWK(106)=IWK(105)+NDV+1                                           ADS13260
C     XTMP2, TEMPORARY X-VECTOR.                                        ADS13270
      IWK(107)=IWK(106)+NDV+1                                           ADS13280
C     XTMP3, TEMPORARY X-VECTOR.                                        ADS13290
      IWK(108)=IWK(107)+NDV+1                                           ADS13300
C     S, SEARCH DIRECTION.                                              ADS13310
      IWK(109)=IWK(108)+NDV+1                                           ADS13320
C     XOLD, X-VECTOR AT ALPHA=0.                                        ADS13330
      IWK(110)=IWK(109)+NDV+1                                           ADS13340
C     B, BTMP MATRIX USED TO FIND USABLE-FEASIBLE DIRECTION.            ADS13350
      IWK(111)=IWK(110)+NDV+1                                           ADS13360
C     GTMP, TEMPORARY STORAGE OF ACTIVE CONSTRAINTS.                    ADS13370
      IWK(112)=IWK(111)+NRB**2                                          ADS13380
C                                                                       ADS13390
      IF (IOPT.LT.4.AND.ISTRAT.LT.8) IWK(112)=IWK(111)                  ADS13400
C                                                                       ADS13410
C     GGU, TEMPORARY STORAGE FOR G-VECTOR.                              ADS13420
      IWK(113)=IWK(112)+NCON+NDV+1                                      ADS13430
C     GGL, TEMPORARY STORAGE FOR G-VECTOR.                              ADS13440
      IWK(114)=IWK(113)+NCON+1                                          ADS13450
C     GG1, TEMPORARY STORAGE FOR G-VECTOR.                              ADS13460
      IWK(115)=IWK(114)+NCON+NDV+1                                      ADS13470
C     GG2, TEMPORARY STORAGE FOR G-VECTOR.                              ADS13480
      IWK(116)=IWK(115)+NCON+1                                          ADS13490
C     DG, GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS W.R.T. ALPHA.         ADS13500
      IWK(117)=IWK(116)+NCON+1                                          ADS13510
C     C.  METHOD OF FEASIBLE DIRECTIONS SUB-PROBLEM.                    ADS13520
      IWK(118)=IWK(117)+NCON+1                                          ADS13530
C     THETA.  METHOD OF FEASIBLE DIRECTIONS.                            ADS13540
      IWK(119)=IWK(118)+NRB                                             ADS13550
C                                                                       ADS13560
      IF (IOPT.LT.4) IWK(119)=IWK(118)                                  ADS13570
C                                                                       ADS13580
C     H-MATRIX.                                                         ADS13590
      IWK(120)=IWK(119)+NRB                                             ADS13600
C                                                                       ADS13610
      IF (IOPT.LT.4) IWK(120)=IWK(119)                                  ADS13620
C                                                                       ADS13630
C     DX.                                                               ADS13640
      IWK(121)=IWK(120)+(NRB*(NRB+1))/2                                 ADS13650
      IF(ISTRAT.EQ.9) IWK(121)=IWK(120)+(NDV+1)*NGT                     ADS13660
C                                                                       ADS13670
      IF (IOPT.GT.3.AND.ISTRAT.LT.8) IWK(121)=IWK(120)                  ADS13680
C                                                                       ADS13690
C     DFG.                                                              ADS13700
      IWK(122)=IWK(121)+NRB                                             ADS13710
C     ALAMDA.                                                           ADS13720
      IWK(123)=IWK(122)+NDV+1                                           ADS13730
C     PPARAM.                                                           ADS13740
      NFORLG=MAX0(NCON,NDV+NRB)                                         ADS13750
      IWK(124)=IWK(123)+NFORLG+1                                        ADS13760
C     XBST.                                                             ADS13770
      IWK(125)=IWK(124)+NCON+1                                          ADS13780
C     VLBTMP.                                                           ADS13790
      IWK(126)=IWK(125)+NDV+1                                           ADS13800
C     VUBTMP.                                                           ADS13810
      IWK(127)=IWK(126)+NDV+1                                           ADS13820
C     PENLT.                                                            ADS13830
      IWK(128)=IWK(127)+NDV+1                                           ADS13840
C     NEXT AVAILABLE STORAGE.                                           ADS13850
      IWK(129)=IWK(128)                                                 ADS13860
      IF (NEQCON.GT.0) IWK(129)=IWK(128)+NCON                           ADS13870
C     STORAGE USED.                                                     ADS13880
      IWK(150)=IWK(129)                                                 ADS13890
C     ------------------------------------------------------------------ADS13900
C     INTEGER ARRAY STORAGE.                                            ADS13910
C     ------------------------------------------------------------------ADS13920
C                                                                       ADS13930
C     ISCALR, INTEGER SCALAR INFORMATION.                               ADS13940
C     RESERVE 100 STORAGE ALLOCATIONS FOR NOW.                          ADS13950
      IWK(151)=1                                                        ADS13960
C     ISIDC, SIDE CONSTRAINT IDENTIFIERS.                               ADS13970
      IWK(152)=IWK(151)+170                                             ADS13980
C     ICQP. IC FOR QP SUB=PROBLEM.                                      ADS13990
      IWK(153)=IWK(152)+NDV+1                                           ADS14000
C     ITMP, TEMPORARY INTEGER ARRAY USED IN STRATEGIES FOR IDG.         ADS14010
      IWK(154)=IWK(153)+NCON                                            ADS14020
C     ITMP2, TEMPORARY INTEGER ARRAY USED IN FEASIBLE DIRECTIONS.       ADS14030
      IWK(155)=IWK(154)+NCON                                            ADS14040
C     NEXT AVAILABLE STORAGE.                                           ADS14050
      IWK(156)=IWK(155)+NRB+NDV                                         ADS14060
C     STORAGE USED.                                                     ADS14070
      IWK(170)=IWK(156)                                                 ADS14080
      RETURN                                                            ADS14090
      END                                                               ADS14100
C --- SUBROUTINE ADS005 (ISTRAT,IOPT,IONED,NDV,WK,IWK)                  ADS14110
C     DIMENSION WK(NRWK),IWK(NRIWK)                                     ADS14120
C     ------------------------------------------------------------------ADS14130
C --- ROUTINE TO INITIALIZE SCALAR PARAMETERS TO THEIR DEFAULT VALUES.  ADS14140
C     ------------------------------------------------------------------ADS14150
C     BY G. N. VANDERPLAATS                             JUNE 1982.      ADS14160
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS14170
C                                                                       ADS14180
C --- INPUT                                                             ADS14190
C                                                                       ADS14200
C --- OUTPUT                                                            ADS14210
C     WK(I), I=1,50 - SCALAR REAL VARIABLES.                            ADS14220
C     IWK(I), I=101,150 - SCALAR INTEGER VARIABLES.                     ADS14230
C                                                                       ADS14240
      SUBROUTINE ADS005 (ISTRAT,IOPT,IONED,NDV,WK,IWK)                  ADS14250
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION WK(1),IWK(1)                                            ADS14260
C     DATA STATEMENTS FOR DEFAULTS.                                     ADS14270
C     ------------------------------------------------------------------ADS14280
C     REAL SCALARS.                                                     ADS14290
C     ------------------------------------------------------------------ADS14300
C     FIRST 50 LOCATIONS ARE RESERVED FOR USER-CONTROLED VALUES.        ADS14310
C     ( 1) ALAMDZ ( 2) BETAMC ( 3) CT     ( 4) CTL    ( 5) CTLMIN       ADS14320
C     ( 6) CTMIN  ( 7) DABALP ( 8) DABOBJ ( 9) DABOBM (10) DABSTR       ADS14330
C     (11) DELALP (12) DELOBJ (13) DELOBM (14) DELSTR (15) DLOBJ1       ADS14340
C     (16) DLOBJ2 (17) DX1    (18) DX2    (19) EPSPEN (20) EXTRAP       ADS14350
C     (21) FDCH   (22) FDCHM  (23) GMULTZ (24) PSAIZ  (25) RMULT        ADS14360
C     (26) RMVLMZ (27) RP     (28) RPMAX  (29) RPMULT (30) RPPMIN       ADS14370
C     (31) RPPRIM (32) SCFO   (33) SCLMIN (34) STOL   (35) THETAZ       ADS14380
C     (36) XMULT  (37) ZRO    (38) PMLT                                 ADS14390
      WK(1)=0.0                                                         ADS14400
      WK(2)=0.                                                          ADS14410
      WK(3)=-0.03                                                       ADS14420
      IF (IOPT.EQ.4) WK(3)=-.1                                          ADS14430
      WK(4)=-0.005                                                      ADS14440
      WK(5)=0.001                                                       ADS14450
      WK(6)=0.004                                                       ADS14460
      WK(7)=0.0001                                                      ADS14470
      IF (IONED.EQ.3.OR.IONED.EQ.8) WK(7)=.001                          ADS14480
      WK(8)=0.0                                                         ADS14490
      WK(9)=0.0                                                         ADS14500
      WK(10)=0.0                                                        ADS14510
      WK(11)=0.005                                                      ADS14520
      IF (IONED.EQ.3.OR.IONED.EQ.8) WK(11)=.05                          ADS14530
      WK(12)=0.001                                                      ADS14540
      WK(13)=0.01                                                       ADS14550
      WK(14)=0.001                                                      ADS14560
      WK(15)=0.1                                                        ADS14570
      WK(16)=1000.                                                      ADS14580
      WK(17)=0.01                                                       ADS14590
      WK(18)=0.2                                                        ADS14600
      WK(19)=-0.05                                                      ADS14610
      WK(20)=5.0                                                        ADS14620
      WK(21)=0.01                                                       ADS14630
      WK(22)=0.001                                                      ADS14640
      WK(23)=10.                                                        ADS14650
      WK(24)=0.95                                                       ADS14660
      WK(25)=5.                                                         ADS14670
      WK(26)=0.2                                                        ADS14680
      IF(ISTRAT.EQ.9) WK(26)=.4                                         ADS14690
      WK(27)=10.                                                        ADS14700
      WK(28)=1.0E+10                                                    ADS14710
      WK(29)=0.2                                                        ADS14720
      WK(30)=1.0E-10                                                    ADS14730
      WK(31)=100.                                                       ADS14740
      WK(32)=0.                                                         ADS14750
      WK(33)=0.001                                                      ADS14760
      WK(34)=0.001                                                      ADS14770
      WK(35)=0.1                                                        ADS14780
      WK(36)=2.618034                                                   ADS14790
      WK(37)=0.00001                                                    ADS14800
      WK(38)=10.                                                        ADS14810
C     ------------------------------------------------------------------ADS14820
C     INTEGER SCALARS.                                                  ADS14830
C     ------------------------------------------------------------------ADS14840
C     ( 1) ICNDIR ( 2) ISCAL  ( 3) ITMAX  ( 4) ITRMOP ( 5) ITRMST       ADS14850
C     ( 6) JONED  ( 7) JTMAX                                            ADS14860
      IWK(1)=NDV+1                                                      ADS14870
      IWK(2)=1                                                          ADS14880
      IWK(3)=40                                                         ADS14890
      IWK(4)=3                                                          ADS14900
      IWK(5)=2                                                          ADS14910
      IWK(6)=IONED                                                      ADS14920
      IWK(7)=20                                                         ADS14930
      RETURN                                                            ADS14940
      END                                                               ADS14950
C --- SUBROUTINE ADS006 (NDV,NCON,VLB,VUB,ZRO,SCFX,SVLB,SVUB,SCFO,SCFG) ADS14960
C     DIMENSION VLB(NDV),VUB(NDV),SVLB(NDV),SVUB(NDV),SCFX(NDV),        ADS14970
C    1SCFG(NCON)                                                        ADS14980
C     ------------------------------------------------------------------ADS14990
C --- ROUTINE TO INITIALIZE SCALE FACTORS.                              ADS15000
C     ------------------------------------------------------------------ADS15010
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS15020
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS15030
C --- INPUT                                                             ADS15040
C      NDV  - NUMBER OF DESIGN VARIABLES.                               ADS15050
C      NCON - NUMBER OF CONSTRAINTS.                                    ADS15060
C      VLB  - LOWER BOUND OF DESIGN VARIABLES.                          ADS15070
C      VUB  - UPPER BOUND OF DESIGN VARIABLES.                          ADS15080
C                                                                       ADS15090
      SUBROUTINE ADS006 (NDV,NCON,VLB,VUB,ZRO,SCFX,SVLB,SVUB,SCFO,SCFG) ADS15100
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION VLB(1),VUB(1),SVLB(1),SVUB(1),SCFX(1),SCFG(1)           ADS15110
C                                                                       ADS15120
      IF (SCFX(1).GT.ZRO) GO TO 20                                      ADS15130
      DO 10 I=1,NDV                                                     ADS15140
      SCFX(I)=1.                                                        ADS15150
10    CONTINUE                                                          ADS15160
20    CONTINUE                                                          ADS15170
      DO 30 I=1,NDV                                                     ADS15180
      SVLB(I)=VLB(I)*SCFX(I)                                            ADS15190
      SVUB(I)=VUB(I)*SCFX(I)                                            ADS15200
30    CONTINUE                                                          ADS15210
C     SCALE FACTOR FOR OBJECTIVE.                                       ADS15220
      IF (SCFO.LE.ZRO) SCFO=1.                                          ADS15230
      IF (NCON.EQ.0) RETURN                                             ADS15240
      IF (SCFG(1).GT.ZRO) RETURN                                        ADS15250
      DO 40 J=1,NCON                                                    ADS15260
      SCFG(J)=1.                                                        ADS15270
40    CONTINUE                                                          ADS15280
      RETURN                                                            ADS15290
      END                                                               ADS15300
      SUBROUTINE ADS007 (JSCAL,INFO,ISTRAT,NDV,NCON,NGT,IC,OBJ,X,G,DF,A,ADS15310
     1NRA,NCOLA,SCFX,SCFO,SCFG,VLB,VUB,IPRT,SXTMP,SGTMP,ZRO,SVLB,SVUB)  ADS15320
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION IC(1),X(1),G(1),DF(1),A(NRA,1),SCFX(1),SCFG(1),VLB(1),VUADS15330
     1B(1),SXTMP(1),SGTMP(1),SVLB(1),SVUB(1)                            ADS15340
C                                                                       ADS15350
        IF (JSCAL.LE.0.OR.JSCAL.GE.7) RETURN                            ADS15360
        GO TO (10,110,260,340,340,260), JSCAL                           ADS15370
10      CONTINUE                                                        ADS15380
C     ------------------------------------------------------------------ADS15390
C     CALCULATE SCALE FACTORS NEWLY.  JSCAL=1                           ADS15400
C     ------------------------------------------------------------------ADS15410
        DO 20 I=1,NDV                                                   ADS15420

	   TOLS = 0.1   !ADDED BY SONGSAK NOV06
         D=DMAX1(TOLS,ABS(X(I)))                                        ADS15430

         SCFX(I)=1./D                                                   ADS15440
         IF(SCFX(I).LT.0.001) SCFX(I)=0.001                             ADS15450
 20     CONTINUE                                                        ADS15460
        DO 30 I=1,NDV                                                   ADS15470
         SVLB(I)=SCFX(I)*VLB(I)                                         ADS15480
         SVUB(I)=SCFX(I)*VUB(I)                                         ADS15490
         IF(SVLB(I).LT.-1.0E+20) SVLB(I)=-1.0E+20                       ADS15500
         IF(SVUB(I).GT.+1.0E+20) SVUB(I)=+1.0E+20                       ADS15510
 30     CONTINUE                                                        ADS15520
C                                                                       ADS15530
        VF=0.                                                           ADS15540
        DO 40 I=1,NDV                                                   ADS15550
 40     VF=VF+DF(I)**2                                                  ADS15560
        VF=SQRT(VF)*FLOAT(NDV)                                          ADS15570

	  TOLS = 0.1   !ADDED BY SONGSAK NOV06
        D=DMAX1(TOLS,VF)                                                ADS15580

        SCFO=1./D                                                       ADS15590
        IF(SCFO.LT.0.001) SCFO=0.001                                    ADS15600
C                                                                       ADS15610
        IF(NGT.EQ.0.OR.NCON.EQ.0) GOTO 100                              ADS15620
        DO 60 K=1,NGT                                                   ADS15630
         J=IC(K)                                                        ADS15640
         VG=0.                                                          ADS15650
         DO 50 I=1,NDV                                                  ADS15660
 50      VG=VG+A(I,K)**2                                                ADS15670
         VG=SQRT(VG)*FLOAT(NDV)                                         ADS15680

	   TOLS = 0.001   !ADDED BY SONGSAK NOV06
         D=DMAX1(TOLS,VG)                                               ADS15690

         SCFG(J)=1./D                                                   ADS15700
C                                                                       ADS15710
         IF(SCFG(J).LT.1.0) SCFG(J)=1.001                               ADS15720
 60     CONTINUE                                                        ADS15730
C                                                                       ADS15740
100     IF(IPRT.LT.3) GOTO 260                                          ADS15750
        WRITE (6,410) SCFO                                              ADS15760
        II=(NDV+4)/5                                                    ADS15770
        I1=1                                                            ADS15780
        I2=5                                                            ADS15790
        DO 103 I=1,II                                                   ADS15800
         IF (I2.GT.NDV) I2=NDV                                          ADS15810
         WRITE (6,420) I1,(SCFX(K),K=I1,I2)                             ADS15820
         I1=I1+5                                                        ADS15830
         I2=I2+5                                                        ADS15840
103     CONTINUE                                                        ADS15850
        IF (NCON.EQ.0.OR.NGT.EQ.0) GOTO 260                             ADS15860
        WRITE (6,430)                                                   ADS15870
        WRITE (6,440) (IC(I),I=1,NGT)                                   ADS15880
        WRITE (6,450)                                                   ADS15890
        DO 106 I=1,NGT                                                  ADS15900
         J=IC(I)                                                        ADS15910
         SGTMP(I)=SCFG(J)                                               ADS15920
106     CONTINUE                                                        ADS15930
        WRITE (6,460) (SGTMP(I),I=1,NGT)                                ADS15940
        IF (NGT.LT.NCON) WRITE (6,470)                                  ADS15950
        GOTO 260                                                        ADS15960
C                                                                       ADS15970
110   CONTINUE                                                          ADS15980
C     ------------------------------------------------------------------ADS15990
C     CALCULATE ADDITIONAL SCALE FACTOR. JSCAL=2                        ADS16000
C     ------------------------------------------------------------------ADS16010
        IF (NGT.EQ.0.OR.NCON.EQ.0) GO TO 260                            ADS16020
        IPP=0                                                           ADS16030
        DO 150 K=1,NGT                                                  ADS16040
         J=IC(K)                                                        ADS16050
         R=ABS(1.-SCFG(J))                                              ADS16060
         IF (R.GT.0.0005) GO TO 150                                     ADS16070
         IF (IPP.EQ.0.AND.IPRT.GE.3) WRITE (6,480)                      ADS16080
         IPP=1                                                          ADS16090
         VG=0.                                                          ADS16100
         DO 145 I=1,NDV                                                 ADS16110
145      VG=VG+A(I,K)**2                                                ADS16120
         VG=SQRT(VG)*FLOAT(NDV)                                         ADS16130

	   TOLS = 0.001   !ADDED BY SONGSAK NOV06
         D=DMAX1(TOLS,VG)                                               ADS16140

         SCFG(J)=1./D                                                   ADS16150
         IF(SCFG(J).LT.1.0) SCFG(J)=1.001                               ADS16160
         G(J)=SCFG(J)*G(J)                                              ADS16170
         IF(IPRT.GE.3) WRITE(6,490) J,SCFG(J)                           ADS16180
150     CONTINUE                                                        ADS16190
C                                                                       ADS16200
260     CONTINUE                                                        ADS16210
C     ------------------------------------------------------------------ADS16220
C     SCALE.  JSCAL=1,2,3,6                                             ADS16230
C     ------------------------------------------------------------------ADS16240
      IF (JSCAL.EQ.2) GO TO 290                                         ADS16250
      OBJ=SCFO*OBJ                                                      ADS16260
      DO 270 I=1,NDV                                                    ADS16270
      X(I)=SCFX(I)*X(I)                                                 ADS16280
270   CONTINUE                                                          ADS16290
      IF (NCON.EQ.0) GO TO 290                                          ADS16300
      DO 280 J=1,NCON                                                   ADS16310
      G(J)=SCFG(J)*G(J)                                                 ADS16320
280   CONTINUE                                                          ADS16330
290   CONTINUE                                                          ADS16340
      IF (JSCAL.EQ.3) RETURN                                            ADS16350
      DO 300 I=1,NDV                                                    ADS16360
      DF(I)=SCFO*DF(I)/SCFX(I)                                          ADS16370
300   CONTINUE                                                          ADS16380
      IF (NGT.EQ.0.OR.NCON.EQ.0) GO TO 320                              ADS16390
      DO 310 I=1,NDV                                                    ADS16400
      DO 310 K=1,NGT                                                    ADS16410
      J=IC(K)                                                           ADS16420
      A(I,K)=SCFG(J)*A(I,K)/SCFX(I)                                     ADS16430
310   CONTINUE                                                          ADS16440
320   CONTINUE                                                          ADS16450
      RETURN                                                            ADS16460
C                                                                       ADS16470
340   CONTINUE                                                          ADS16480
C     ------------------------------------------------------------------ADS16490
C     UN-SCALE.  JSCAL=4                                                ADS16500
C     ------------------------------------------------------------------ADS16510
      DO 350 I=1,NDV                                                    ADS16520
      X(I)=X(I)/SCFX(I)                                                 ADS16530
350   CONTINUE                                                          ADS16540
      IF (JSCAL.EQ.5) RETURN                                            ADS16550
      OBJ=OBJ/SCFO                                                      ADS16560
      IF (NCON.EQ.0) RETURN                                             ADS16570
      DO 360 J=1,NCON                                                   ADS16580
      G(J)=G(J)/SCFG(J)                                                 ADS16590
360   CONTINUE                                                          ADS16600
      RETURN                                                            ADS16610
C                                                                       ADS16620
C     ------------------------------------------------------------------ADS16630
C                            FORMATS                                    ADS16640
C     ------------------------------------------------------------------ADS16650
410   FORMAT (//,5X,37HSCALING INFORMATION CALCULATED BY ADS,//,5X,     ADS16660
     139HOBJECTIVE FUNCTION SCALE FACTOR, SCFO =,1PE12.5,//,5X,         ADS16670
     229HDESIGN VARIABLE SCALE FACTORS)                                 ADS16680
420   FORMAT (I10,5(1PE12.5))                                           ADS16690
430   FORMAT (/,5X,24HCONSTRAINT SCALE FACTORS,//,5X,                   ADS16700
     118HCONSTRAINT NUMBERS)                                            ADS16710
440   FORMAT (5I12)                                                     ADS16720
450   FORMAT (/,5X,'SCALE FACTORS')                                     ADS16730
460   FORMAT (5X,1PE12.5,1PE12.5,1PE12.5,1PE12.5,1PE12.5)               ADS16740
470   FORMAT (/,5X,'ALL OTHER SCALE FACTORS = 1.0')                     ADS16750
480   FORMAT (/,5X,                                                     ADS16760
     1'ADDITIONAL CONSTRAINT SCALE FACTORS CALCULATED BY ADS',/,5X,     ADS16770
     1'CONSTRAINT   SCALE FACTOR')                                      ADS16780
490   FORMAT (5X,I6,6X,1PE12.5)                                         ADS16790
      END                                                               ADS16800
C --- SUBROUTINE ADS008 (NDV,ISTRAT,NCON,NGT,A,NRA,IC,G,IDG,DF,ALAMDA,  ADS16810
C    *PPARAM,RP,EPSPEN)                                                 ADS16820
C     DIMENSION G(NCON),IDG(NCON),A(NRA,NCOLA),IC(NCOLA),DF(NDV),       ADS16830
C    *ALAMDA(NCON),PPARAM(NCON)                                         ADS16840
C     ------------------------------------------------------------------ADS16850
C --- ROUTINE TO CALCULATE THE GRADIENTS OF PSEUDO OBJECTIVE FUNCTION.  ADS16860
C     ------------------------------------------------------------------ADS16870
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS16880
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS16890
C --- INPUT                                                             ADS16900
C      NDV    ; NUMBER OF DESIGN VARIABLES.                             ADS16910
C      ISTRAT ; OPTIMIZATION STRATEGY IDENTIFIER.                       ADS16920
C      NGT   ; NUMBER OF ACTIVE CONSTRAINTS.                            ADS16930
C      A      ; ARRAY CONTAINING GRADIENTS OF CONSTRAINTS IDENTIFIED BY ADS16940
C               ARRAY IC.                                               ADS16950
C      NRA    ; NUMBER OF ROWS IN ARRAY A.                              ADS16960
C      IC     ; ARRAY IDENTIFYING NGT ACTIVE CONSTRAINTS.               ADS16970
C      G      ; CONSTRAINT VALUES ASSOCIATED WITH CURRENT DESIGN.       ADS16980
C      IDG    ; CONSTRAINT TYPE IDENTIFIER.                             ADS16990
C      ALAMDA ; ARRAY FOR LAGRANGE MULTIPLIER.                          ADS17000
C      PPARAM ; ARRAY FOR PENALTY PARAMETER.                            ADS17010
C      RP     ; PENALTY PARAMETER FOR EXTERIOR PENALTY FUNCTION METHOD. ADS17020
C      EPSPEN ; TRANSITION POINT.                                       ADS17030
C --- OUTPUT                                                            ADS17040
C      DF     ; GRADIENTS OF PSEUDO PBJECTIVE FUNCTION.                 ADS17050
      SUBROUTINE ADS008 (NDV,ISTRAT,NCON,NGT,A,NRA,IC,G,IDG,DF,ALAMDA,PPADS17060
     1ARAM,RP,EPSPEN)                                                   ADS17070
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IDG(1),A(NRA,1),IC(1),DF(1),ALAMDA(1),PPARAM(1)    ADS17080
C                                                                       ADS17090
      GO TO (10,30,30,30,110), ISTRAT                                   ADS17100
10    CONTINUE                                                          ADS17110
C     ------------------------------------------------------------------ADS17120
C     ISTRAT=1 . EXTERIOR PENALTY.                                      ADS17130
C     ------------------------------------------------------------------ADS17140
      RPEX=2.*RP                                                        ADS17150
      DO 20 K=1,NGT                                                     ADS17160
      J=IC(K)                                                           ADS17170
      IF (IDG(J).GE.0.AND.G(J).LE.0.) GO TO 20                          ADS17180
      GMULT=RPEX*G(J)                                                   ADS17190
      CALL ADS508 (NDV,K,GMULT,A,NRA,DF)                                ADS17200
20    CONTINUE                                                          ADS17210
      RETURN                                                            ADS17220
C                                                                       ADS17230
30    CONTINUE                                                          ADS17240
C     ------------------------------------------------------------------ADS17250
C     ISTRAT=2,3,4 . INTERIOR PENALTY.                                  ADS17260
C     ------------------------------------------------------------------ADS17270
      IF (ISTRAT.LE.3) GO TO 50                                         ADS17280
      AA=-100.                                                          ADS17290
      DO 40 I=1,NCON                                                    ADS17300
      IF (G(I).GT.AA) AA=G(I)                                           ADS17310
40    CONTINUE                                                          ADS17320
      AA=AA/EPSPEN                                                      ADS17330
      SHAPE=2./3.                                                       ADS17340
      IF (AA.GT.SHAPE) AA=SHAPE                                         ADS17350
      SHAPE=1./(3.*(1.-AA))                                             ADS17360
50    CONTINUE                                                          ADS17370
      C=-1./EPSPEN**2                                                   ADS17380
      DO 100 K=1,NGT                                                    ADS17390
      J=IC(K)                                                           ADS17400
      IF (IDG(J).GE.0) GO TO 60                                         ADS17410
      GMULT=2.*PPARAM(J)*G(J)                                           ADS17420
      GO TO 90                                                          ADS17430
60    CONTINUE                                                          ADS17440
      IF (G(J).GT.EPSPEN) GO TO 70                                      ADS17450
      GMULT=1./G(J)**2                                                  ADS17460
      GMULT=GMULT*PPARAM(J)                                             ADS17470
      GO TO 90                                                          ADS17480
70    CONTINUE                                                          ADS17490
      GMULT=-1.                                                         ADS17500
      IF (ISTRAT.EQ.2) GO TO 80                                         ADS17510
      AA=G(J)/EPSPEN-1.                                                 ADS17520
      GMULT=GMULT+2.*AA                                                 ADS17530
      IF (ISTRAT.EQ.3) GO TO 80                                         ADS17540
      GMULT=GMULT+3.*SHAPE*AA*AA                                        ADS17550
80    GMULT=C*GMULT*PPARAM(J)                                           ADS17560
90    CONTINUE                                                          ADS17570
      CALL ADS508 (NDV,K,GMULT,A,NRA,DF)                                ADS17580
100   CONTINUE                                                          ADS17590
      RETURN                                                            ADS17600
C                                                                       ADS17610
110   CONTINUE                                                          ADS17620
C     ------------------------------------------------------------------ADS17630
C     ISTRAT=5 . AUGMENTED LAGRANGIAN MULTIPLIER.                       ADS17640
C     ------------------------------------------------------------------ADS17650
      DO 130 K=1,NGT                                                    ADS17660
      J=IC(K)                                                           ADS17670
      GMULT=ALAMDA(J)+PPARAM(J)*G(J)                                    ADS17680
      IF (IDG(J).LT.0) GO TO 120                                        ADS17690
      IF (GMULT.LE.0.) GO TO 130                                        ADS17700
120   CONTINUE                                                          ADS17710
      CALL ADS508 (NDV,K,GMULT,A,NRA,DF)                                ADS17720
130   CONTINUE                                                          ADS17730
      RETURN                                                            ADS17740
      END                                                               ADS17750
C --- SUBROUTINE ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,CTLMIN,NAC,ADS17760
C    *NVC)                                                              ADS17770
C     DIMENSION A(NRA,NCOLA),IC(NCOLA),IDG(NCON),G(NCON)                ADS17780
C     ------------------------------------------------------------------ADS17790
C --- ROUTINE TO REORDER IC AND A MATRIX                                ADS17800
C     ------------------------------------------------------------------ADS17810
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS17820
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS17830
C --- INPUT                                                             ADS17840
C      NDV  - NUMBER OF DESIGN VARIABLES.                               ADS17850
C      A    - ARRAY CONTAINING GRADIENTS OF CONSTRAINTS IDENTIFIED BY   ADS17860
C             ARRAY IC.                                                 ADS17870
C      NRA  - NUMBER OF ROWS IN ARRAY A.                                ADS17880
C      NGT  - NUMBER OF CONSTRAINT GRADIENTS STORED IN A.               ADS17890
C      IC   - ARRAY IDENTIFYING CONSTRAINTS THAT ARE ACTIVE.            ADS17900
C      IDG  - CONSTRAINTS TYPE IDENTIFIER.                              ADS17910
C      G    - CONSTRAINT VALUES.                                        ADS17920
C      CT   - CONSTRAINT TOLORANCE DEFINING ACTIVE NONLINEAR            ADS17930
C             CONSTRAINTS.                                              ADS17940
C      CTMIN- CONSTRAINT TOLORANCE DEFINING VIOLATION OF NONLINEAR      ADS17950
C             CONSTRAINTS.                                              ADS17960
C      CTL  - CONSTRAINT TOLORANCE DEFINING ACTIVE LINEAR CONSTRAINTS.  ADS17970
C      CTLMIN- CONSTRAINT TOLORANCE DEFINING VIOLATION OF LINEAR        ADS17980
C             CONSTRAINTS.                                              ADS17990
C --- OUTPUT                                                            ADS18000
C      A    - REORDERED ARRAY A.                                        ADS18010
C      IC   - REORDERED ARRAY IC.                                       ADS18020
C      NAC  - NUMBER OF ACTIVE CONSTRAINTS.                             ADS18030
C      NVC  - NUMBER OF VIOLATED CONSTRAINTS.                           ADS18040
C                                                                       ADS18050
C --- NOTE                                                              ADS18060
C     ORDERING OF A IS EQUALITY CONSTRAINTS, ACTIVE AND VIOLATED        ADS18070
C     INEQUALITY CONSTRAINTS, REMAINING CONSTRAINTS IN OREDER OF        ADS18080
C     DECREASING CRITICALITY.                                           ADS18090
C                                                                       ADS18100
      SUBROUTINE ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,CTLMIN,NAC,ADS18110
     1NVC)                                                              ADS18120
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(NRA,1),IC(1),IDG(1),G(1)                              ADS18130
C                                                                       ADS18140
      NAC=0                                                             ADS18150
      NVC=0                                                             ADS18160
      DO 90 K=1,NGT                                                     ADS18170
      ICK=IC(K)                                                         ADS18180
      C1=CT                                                             ADS18190
      C2=CTMIN                                                          ADS18200
      IF (IABS(IDG(ICK)).LT.2) GO TO 30                                 ADS18210
C     LINEAR CONSTRAINT.                                                ADS18220
      C1=CTL                                                            ADS18230
      C2=CTLMIN                                                         ADS18240
30    CONTINUE                                                          ADS18250
      GMAX=G(ICK)                                                       ADS18260
      IF (GMAX.LT.C1) GO TO 40                                          ADS18270
C     THIS CONSTRAINT IS ACTIVE OR VIOLATED.  NO NEED TO SHIFT.         ADS18280
      IF (GMAX.GE.C1.AND.GMAX.LE.C2) NAC=NAC+1                          ADS18290
      IF (GMAX.GT.C2) NVC=NVC+1                                         ADS18300
      GO TO 90                                                          ADS18310
40    CONTINUE                                                          ADS18320
      IF (K.EQ.NGT) GO TO 90                                            ADS18330
C     SEARCH FOR NEXT ACTIVE OR VIOLATED CONSTRAINT.                    ADS18340
      K1=K+1                                                            ADS18350
      IK=0                                                              ADS18360
      DO 60 J=K1,NGT                                                    ADS18370
      ICJ=IC(J)                                                         ADS18380
      IK=J                                                              ADS18390
      C1=CT                                                             ADS18400
      C2=CTMIN                                                          ADS18410
      IF (IABS(IDG(ICJ)).LT.2) GO TO 50                                 ADS18420
C     LINEAR CONSTRAINT.                                                ADS18430
      C1=CTL                                                            ADS18440
      C2=CTLMIN                                                         ADS18450
50    CONTINUE                                                          ADS18460
      GG=G(ICJ)                                                         ADS18470
      IF (GG.GE.C1.AND.GG.LE.C2) NAC=NAC+1                              ADS18480
      IF (GG.GT.C2) NVC=NVC+1                                           ADS18490
      IF (GG.GE.C1) GO TO 70                                            ADS18500
60    CONTINUE                                                          ADS18510
      GO TO 90                                                          ADS18520
C     SHIFT.                                                            ADS18530
70    ICS=IC(K)                                                         ADS18540
      IC(K)=IC(IK)                                                      ADS18550
      IC(IK)=ICS                                                        ADS18560
      DO 80 I=1,NDV                                                     ADS18570
      TMP=A(I,K)                                                        ADS18580
      A(I,K)=A(I,IK)                                                    ADS18590
      A(I,IK)=TMP                                                       ADS18600
80    CONTINUE                                                          ADS18610
90    CONTINUE                                                          ADS18620
      NN=NAC+NVC+1                                                      ADS18630
      IF (NN.GE.NGT) RETURN                                             ADS18640
      NGT1=NGT-1                                                        ADS18650
      DO 120 K=NN,NGT1                                                  ADS18660
      ICK=IC(K)                                                         ADS18670
      GMAX=G(ICK)                                                       ADS18680
      K1=K+1                                                            ADS18690
      IK=0                                                              ADS18700
      DO 100 J=K1,NGT                                                   ADS18710
      ICJ=IC(J)                                                         ADS18720
      IF (GMAX.GT.G(ICJ)) GO TO 100                                     ADS18730
      IK=J                                                              ADS18740
      GMAX=G(ICK)                                                       ADS18750
100   CONTINUE                                                          ADS18760
      IF (IK.EQ.0) GO TO 120                                            ADS18770
      ICS=IC(K)                                                         ADS18780
      IC(K)=IC(IK)                                                      ADS18790
      IC(IK)=ICS                                                        ADS18800
      DO 110 I=1,NDV                                                    ADS18810
      TMP=A(I,K)                                                        ADS18820
      A(I,K)=A(I,IK)                                                    ADS18830
110   A(I,IK)=TMP                                                       ADS18840
120   CONTINUE                                                          ADS18850
      RETURN                                                            ADS18860
      END                                                               ADS18870
C --- SUBROUTINE ADS010 (ISTRAT,NDV,OBJ,X,DELOBJ,DABOBJ,DELOBM,DABOBM,DAADS18880
C    1BSTR,ZRO)                                                         ADS18890
C     DIMENSION X(NDV)                                                  ADS18900
C     ------------------------------------------------------------------ADS18910
C --- ROUTINE TO CALCULATE PARAMETERS FOR CONVERGENCE CRITERIA.         ADS18920
C     ------------------------------------------------------------------ADS18930
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY  1983   ADS18940
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS18950
C --- INPUT                                                             ADS18960
C      ISTRAT - OPTIMIZATION STRATEGY IDENTIFIER.                       ADS18970
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS18980
C      OBJ    - OBJECTIVE.                                              ADS18990
C      X      - DESIGN VARIABLES.                                       ADS19000
C      DELOBJ - RELATIVE CHANGE CRITERION FOR OPTIMIZATION LEVEL.       ADS19010
C      DABOBJ - ABSOLUTE CHANGE CRITERION FOR OPTIMIZATION LEVEL.       ADS19020
C      DELOBM - MINIMUM RELATIVE CHANGE CRITERION FOR OPTIMIZATION LEVELADS19030
C               IN SUMT METHOD.                                         ADS19040
C      DABOBM - MINIMUM ABSOLUTE CHANGE CRITERION FOR OPTIMIZATION LEVELADS19050
C               IN SUMT METHOD.                                         ADS19060
C      DABSTR - ABSOLUTE CHANGE CRITERION FOR STRATEGY LEVEL.           ADS19070
C --- OUTPUT                                                            ADS19080
C                                                                       ADS19090
      SUBROUTINE ADS010 (ISTRAT,NDV,OBJ,X,DELOBJ,DABOBJ,DELOBM,DABOBM,DAADS19100
     1BSTR,ZRO)                                                         ADS19110
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1)                                                    ADS19120
C                                                                       ADS19130
C     RELATIVE CHANGE.                                                  ADS19140
      IF (DELOBM.LE.ZRO) DELOBM=.01                                     ADS19150
      IF (DELOBJ.LE.ZRO) DELOBJ=.001                                    ADS19160
C     ABSOLUTE CHANGE.                                                  ADS19170
      ABSOB=ABS(OBJ)                                                    ADS19180
      DABA=0.001*ABSOB                                                  ADS19190
      IF (DABA.LE.ZRO) DABA=.001                                        ADS19200
      IF (ISTRAT.NE.7) GO TO 20                                         ADS19210
C     IF STRATEGY IS METHOD OF CENTER, THEN CALCULATE THE AVERAGE VALUE ADS19220
C     OF INITIAL DESIGN VARIABLES.                                      ADS19230
      XT=0.                                                             ADS19240
      DO 10 I=1,NDV                                                     ADS19250
      ABSX=ABS(X(I))                                                    ADS19260
      XT=XT+ABSX                                                        ADS19270
10    CONTINUE                                                          ADS19280
      XT=XT/NDV                                                         ADS19290
      IF (XT.LT.1.) XT=1.                                               ADS19300
      DABA=0.001*XT                                                     ADS19310
20    CONTINUE                                                          ADS19320
      IF(DABOBJ.LE.ZRO) DABOBJ=DABA                                     ADS19330
      IF (DABOBM.LE.ZRO) DABOBM=2.*DABA                                 ADS19340
      IF (DABSTR.LE.ZRO) DABSTR=DABA                                    ADS19350
      RETURN                                                            ADS19360
      END                                                               ADS19370



C --- SUBROUTINE ADS101(ICALC,X1,F1,FP1,X2,F2,AZ,A1)                    ADS19380
C     ------------------------------------------------------------------ADS19390
C --- LINEAR POLYNOMIAL APPROXIMATION.                                  ADS19400
C     ------------------------------------------------------------------ADS19410
C     BY G. N. VANDERPLAATS                            AUG. 1980        ADS19420
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS19430
C                                                                       ADS19440
C --- INPUT                                                             ADS19450
C     ICALC - CALCULATION CONTROL.  IF ICALC=1 DO ONE-POINT             ADS19460
C             INTERPOLATION.  IF ICALC=2, DO TWO-POINT INTERPOLATION.   ADS19470
C     X1,F1,FP1 - FIRST X-F PAIR AND SLOPE.                             ADS19480
C                 F1 AND YP1 ARE FUNCTION VALUE AND SLOPE AT X=X1.      ADS19490
C                 FP1 IS NOT USED IF ICALC=2.                           ADS19500
C     X2,F2 - SECOND POINT AT WHICH FUNCTION IS GIVEN.  F2=FUNCTION     ADS19510
C             VALUE AT X=X2.  X2 AND F2 ARE NOT USED IF ICALC=1.        ADS19520
C                                                                       ADS19530
C --- OUTPUT                                                            ADS19540
C     AZ,A1 - POLYNOMIAL COEFFICIENTS OF LINEAR INTERPOLATING           ADS19550
C             POLYNOMIAL.                                               ADS19560
C                                                                       ADS19570
C     APPROXIMATING FUNCTION IS F = AZ + A1*X.                          ADS19580
C                                                                       ADS19590
      SUBROUTINE ADS101 (ICALC,X1,F1,FP1,X2,F2,AZ,A1)                   ADS19600
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS19610
C     EVALUATE COEFFICIENTS.                                            ADS19620
      IF (ICALC.EQ.1) A1=FP1                                            ADS19630
      IF (ICALC.GT.1) A1=(F2-F1)/(X2-X1)                                ADS19640
      AZ=F1-A1*X1                                                       ADS19650
      RETURN                                                            ADS19660
      END                                                               ADS19670
C --- SUBROUTINE ADS102(ICALC,X1,F1,FP1,X2,F2,X3,F3,AZ,A1,A2)           ADS19680
C     ------------------------------------------------------------------ADS19690
C --- QUADRATIC POLYNOMIAL APPROXIMATION.                               ADS19700
C     ------------------------------------------------------------------ADS19710
C     BY G. N. VANDERPLAATS                           AUG. 1980         ADS19720
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS19730
C                                                                       ADS19740
C --- INPUT                                                             ADS19750
C     ICALC - CALCULATION CONTROL.  IF ICALC=1, DO TWO-POINT            ADS19760
C             APPROXIMATION.  IF ICALC=2, DO THREE-POINT APPROXIMATION. ADS19770
C     X1,F1,FP1 - FIRST X-F PAIR AND SLOPE.                             ADS19780
C                 F1 AND FP1 ARE FUNCTION VALUE AND SLOPE AT X=X1.      ADS19790
C                 FP1 IS NOT USED IF ICALC=2.                           ADS19800
C     X2,F2 - SECOND X-F PAIR.                                          ADS19810
C     X3,F3 - THIRD X-F PAIR.  NOT USED IF ICALC=1.                     ADS19820
C                                                                       ADS19830
C --- OUTPUT                                                            ADS19840
C     AZ,A1,A2 - COEFFICIENTS OF INTERPOLATING POLYNOMIAL.              ADS19850
C                                                                       ADS19860
C     APPROXIMATING FUNCTION IS                                         ADS19870
C         F = AZ + A1*X + A2*X**2                                       ADS19880
C                                                                       ADS19890
      SUBROUTINE ADS102 (ICALC,X1,F1,FP1,X2,F2,X3,F3,AZ,A1,A2)          ADS19900
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS19910
C     EVALUATE COEFFICIENTS.                                            ADS19920
      IF (ICALC.GT.1) GO TO 10                                          ADS19930
C     TWO-POINT APPROXIMATION, GIVEN (X1,F1,FP1), (X2,F2).              ADS19940
      A2=((F2-F1)/(X2-X1)-FP1)/(X2-X1)                                  ADS19950
      A1=FP1-2.*A2*X1                                                   ADS19960
      AZ=F1-A1*X1-A2*X1*X1                                              ADS19970
      RETURN                                                            ADS19980
C     THREE-POINT APPROXIMATION, GIVEN (X1,F1), (X2,F2), (X3,F3).       ADS19990
10    CONTINUE                                                          ADS20000
      A1=(F2-F1)/(X2-X1)                                                ADS20010
      A2=((F3-F1)/(X3-X1)-A1)/(X3-X2)                                   ADS20020
      A1=A1-(X2+X1)*A2                                                  ADS20030
      AZ=F1-A1*X1-A2*X1*X1                                              ADS20040
      RETURN                                                            ADS20050
      END                                                               ADS20060
C --- SUBROUTINE ADS103 (ICALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,AZ,A1,A2,ADS20070
C    * A3)                                                              ADS20080
C     ------------------------------------------------------------------ADS20090
C --- CUBIC POLYNOMIAL APPROXIMATION.                                   ADS20100
C     ------------------------------------------------------------------ADS20110
C     BY G. N. VANDERPLAATS                         AUG. 1980           ADS20120
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS20130
C                                                                       ADS20140
C --- INPUT                                                             ADS20150
C     ICALC - CALCULATION CONTROL.  IF ICALC=1, DO THREE-POINT          ADS20160
C             CURVE FIT.  IF ICALC=2, DO FOUR-POINT CURVE FIT.          ADS20170
C     X1,F1,FP1 - FIRST X-F PAIR AND SLOPE.                             ADS20180
C                 F1 AND FP1 ARE FUNCTION VALUE AND SLOPE.              ADS20190
C                 IF ICALC=2, FP1 IS NOT USED.                          ADS20200
C     X2,F2 - SECOND X-F PAIR.                                          ADS20210
C     X3,F3 - THIRD X-F PAIR.                                           ADS20220
C     X4,F4 - FOURTH X-F PAIR.  NOT USED IF ICALC=1.                    ADS20230
C     ZRO - SMALL POSITIVE NUMBER TO INDICATE NUMERICAL ZERO.           ADS20240
C         INDICATOR OF MACHINE ACCURACY.                                ADS20250
C                                                                       ADS20260
C --- OUTPUT                                                            ADS20270
C     AZ,A1,A2,A3 - COEFFICIENTS OF INTERPOLATING POLYNOMIAL.           ADS20280
C                                                                       ADS20290
C     APPROXIMATING FUNCTION IS                                         ADS20300
C         F = AZ + A1*X + A2*X**2 + A3*X**3                             ADS20310
C                                                                       ADS20320
      SUBROUTINE ADS103 (ICALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,AZ,A1,A2,ADS20330
     1A3)                                                               ADS20340
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS20350
C     EVALUTE COEFFICIENTS.                                             ADS20360
C     CONSTANTS.                                                        ADS20370
      X21=X2-X1                                                         ADS20380
      X31=X3-X1                                                         ADS20390
      X32=X3-X2                                                         ADS20400
      X11=X1**2                                                         ADS20410
      X22=X2**2                                                         ADS20420
      X33=X3**2                                                         ADS20430
      IF (ICALC.GT.1) GO TO 10                                          ADS20440
C     THREE-POINT APPROXIMATION, GIVEN (X1,F1,FP1), (X2,F2), (X3,F3).   ADS20450
      A3=(X21*(F3-F1)/X31-X31*(F2-F1)/X21+X32*FP1)/(X21*X31*X32)        ADS20460
      A2=((F2-F1)/X21-FP1)/X21-(2.*X1+X2)*A3                            ADS20470
      A1=FP1-2.*A2*X1-3.*A3*X1*X1                                       ADS20480
      AZ=F1-X1*(A1+X1*(A2+A3*X1))                                       ADS20490
      RETURN                                                            ADS20500
C     FOUR-POINT APPROXIMATION, GIVEN (X1,F1), (X2,F2), (X3,F3),        ADS20510
C     (X4,F4).                                                          ADS20520
10    CONTINUE                                                          ADS20530
      X42=X4-X2                                                         ADS20540
      X41=X4-X1                                                         ADS20550
      Q1=X33*X3*X21-X22*X2*X31+X11*X1*X32                               ADS20560
      Q2=(X4**3)*X21-X22*X2*X41+X11*X1*X42                              ADS20570
      Q3=X32*X21*X31                                                    ADS20580
      Q4=X42*X21*X41                                                    ADS20590
      Q5=F3*X21-F2*X31+F1*X32                                           ADS20600
      Q6=F4*X21-F2*X41+F1*X42                                           ADS20610
      DENOM=Q2*Q3-Q1*Q4                                                 ADS20620
      A3=0.                                                             ADS20630
      ANUM=Q3*Q6-Q4*Q5                                                  ADS20640
      CHK=ZRO*ANUM                                                      ADS20650
      IF (ABS(CHK).LT.ABS(DENOM).AND.ABS(DENOM).GT.1.0E-15) A3=ANUM/DENOADS20660
     1M                                                                 ADS20670
      A2=(Q5-Q1*A3)/Q3                                                  ADS20680
      A1=(F2-F1)/X21-(X22*X2-X11*X1)*A3/X21-(X1+X2)*A2                  ADS20690
      AZ=F1-X1*(A1+X1*(A2+X1*A3))                                       ADS20700
      RETURN                                                            ADS20710
      END                                                               ADS20720
C --- SUBROUTINE ADS104(ICALC,AZ,A1,A2,A3,ZRO,XZ1,XZ2,XZ3,IER)          ADS20730
C     ------------------------------------------------------------------ADS20740
C --- FIND ZERO(S) OF POLYNOMIAL TO 3RD ORDER.                          ADS20750
C     ------------------------------------------------------------------ADS20760
C     BY G. N. VANDERPLAATS                             AUG. 1980       ADS20770
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS20780
C                                                                       ADS20790
C --- INPUT                                                             ADS20800
C     ICALC - CALCULATION CONTROL.                                      ADS20810
C             ICALC=1, FIRST ORDER POLYNOMIAL.  FIND ONLY FIRST ZERO.   ADS20820
C             ICALC=2, SECOND ORDER POLYNOMIAL.  FIND TWO ZEROES.       ADS20830
C             ICALC=3, THIRD ORDER POLYNOMIAL.  FIND THREE ZEROES IF    ADS20840
C             ALL ARE REAL.                                             ADS20850
C     AZ, A1, A2, A3 - COEFFICIENTS OF THE POLYNOMIAL,                  ADS20860
C         Y = AZ + A1*X + A2*X**2 + A3*X**3.                            ADS20870
C     ZRO - SMALL POSITIVE NUMBER TO INDICATE NUMERICAL ZERO.           ADS20880
C         INDICATOR OF MACHINE ACCURACY.  ALSO THE PRECISION WITH WHICH ADS20890
C         FIRST ZERO IS FOUND IN CUBIC INTERPOLATION.                   ADS20900
C                                                                       ADS20910
C --- OUTPUT                                                            ADS20920
C     XZ1, XZ2, XZ3 - THE VALUES OF X AT WHICH Y = 0.                   ADS20930
C         ONLY REAL VALUES ARE FOUND.                                   ADS20940
C         IF ICALC = 1, ONLY XZ1 IS FOUND.                              ADS20950
C         IF ICALC = 2, XZ1 AND XZ2 ARE FOUND.                          ADS20960
C         IF ICALC = 3, XZ1, XZ2 AND XZ3 ARE FOUND.                     ADS20970
C     IER - ERROR CODE.  IF IER=1, NEWTON'S METHOD DID NOT CONVERGE     ADS20980
C         IN 50 ITERATIONS.  THE RESULTS MAY NOT BE VALID.              ADS20990
C --- NOTE.  IF XZ1, XZ2 OR XZ3 ARE NOT FOUND, OR ARE IMAGINARY,        ADS21000
C     A VALUE OF -1.0E+20 IS RETURNED.                                  ADS21010
      SUBROUTINE ADS104 (ICALC,AZ,A1,A2,A3,ZRO,XZ1,XZ2,XZ3,IER)         ADS21020
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS21030
      ZRO10=10.*ZRO                                                     ADS21040
C     INITIALIZE OUTPUT PARAMETERS.                                     ADS21050
      XZ1=-1.0E+20                                                      ADS21060
      XZ2=-1.0E+20                                                      ADS21070
      XZ3=-1.0E+20                                                      ADS21080
      IF (ICALC.GT.1) GO TO 20                                          ADS21090
10    CONTINUE                                                          ADS21100
C     ICALC = 1.                                                        ADS21110
C     LINEAR APPROXIMATION.  FIND XZ1 ONLY.                             ADS21120
      IF (ABS(A1).GT.ZRO) XZ1=-AZ/A1                                    ADS21130
      RETURN                                                            ADS21140
20    IF (ICALC.GT.2) GO TO 40                                          ADS21150
C     ICALC = 2.                                                        ADS21160
30    CONTINUE                                                          ADS21170
C     QUADRATIC APPROXIMATION.  FIND XZ1 AND XZ2.                       ADS21180
C     IF QUADRATIC TERM IS SMALL, GO USE LINEAR INTERPOLATION.          ADS21190
      IF (ABS(A2).LT.ZRO10) GO TO 10                                    ADS21200
C     Y=0.                                                              ADS21210
      B=A1**2-4.*A2*AZ                                                  ADS21220
      IF (B.LT.-1.0E-6) RETURN                                          ADS21230
      IF (B.LT.0.) B=0.                                                 ADS21240
      IF (ABS(A2).GT.ZRO) XX=-.5*A1/A2                                  ADS21250
      IF (ABS(A2).GT.ZRO) YY=.5*SQRT(B)/A2                              ADS21260
      X1=XX+YY                                                          ADS21270
      X2=XX-YY                                                          ADS21280
      XZ1=DMIN1(X1,X2)                                                  ADS21290
      XZ2=DMAX1(X1,X2)                                                  ADS21300
      RETURN                                                            ADS21310
40    CONTINUE                                                          ADS21320
C     ICALC = 3.                                                        ADS21330
C     CUBIC APPROXIMATION.                                              ADS21340
C     FIND XZ1, XZ2, XZ3 IF ALL ARE REAL.                               ADS21350
C     IF CUBIC TERM IS SMALL, GO DO QUADRATIC INTERPOLATION.            ADS21360
      IF (ABS(A3).LT.ZRO10) GO TO 30                                    ADS21370
C     FIND XZ1 BY NEWTON'S METHOD.                                      ADS21380
C     INITIAL ESTIMATE.                                                 ADS21390
      X=0.                                                              ADS21400
      IF (ABS(A1).GT.ZRO) X=-AZ/A1                                      ADS21410
C     ITERATE.                                                          ADS21420
      KOUNT=0                                                           ADS21430
50    Y=AZ+X*(A1+X*(A2+X*A3))                                           ADS21440
      YP=A1+X*(2.*A2+3.*X*A3)                                           ADS21450
      IF (ABS(YP).LT.1.) YP=0.9*YP/(ABS(YP)+ZRO)                        ADS21460
      IF (ABS(YP).LT.ZRO) YP=ZRO                                        ADS21470
      DX=-Y/YP                                                          ADS21480

	TOLS = 0.01   !ADDED BY SONGSAK NOV06
      DEL=DMAX1(ABS(X),TOLS)                                            ADS21490

      DEL=ABS(DX)/DEL                                                   ADS21500
      KOUNT=KOUNT+1                                                     ADS21510
      IF (DEL.LT.ZRO.OR.KOUNT.GT.50) GO TO 60                           ADS21520
      X=X+DX                                                            ADS21530
      GO TO 50                                                          ADS21540
C     IF NON-CONVERGENCE, SET ERROR CODE TO 1.                          ADS21550
60    IER=1                                                             ADS21560
      IF (DEL.LT.ZRO) IER=0                                             ADS21570
      XZ1=X                                                             ADS21580
C     REDUCE POLYNOMIAL TO SECOND ORDER AND SOLVE FOR REMAINING ROOTS.  ADS21590
      B3=A3                                                             ADS21600
      B2=A2+XZ1*A3                                                      ADS21610
      B1=A1+XZ1*B2                                                      ADS21620
C     LINEAR ESTIMATE.                                                  ADS21630
      IF (ABS(B2).GT.ZRO) XZ2=-B1/B2                                    ADS21640
C     IF QUADRATIC TERM IS SMALL, END.                                  ADS21650
      IF (ABS(B3).LT.ZRO10) GO TO 80                                    ADS21660
      AC=B2*B2-4.*B3*B1                                                 ADS21670
C     IF ROOTS ARE COMPLEX, END SEARCH.                                 ADS21680
      IF (AC.GE.0.) GO TO 70                                            ADS21690
      XZ2=-1.0E+20                                                      ADS21700
      RETURN                                                            ADS21710
70    AC=SQRT(AC)                                                       ADS21720
      IF (ABS(B3).GT.ZRO) XZ2=(-B2+AC)/(2.*B3)                          ADS21730
      IF (ABS(B3).GT.ZRO) XZ3=(-B2-AC)/(2.*B3)                          ADS21740
C     PUT XZ1, XZ2 AND XZ3 IN ASCENDING ORDER.                          ADS21750
80    CONTINUE                                                          ADS21760
      IF (XZ2.GT.XZ1) GO TO 90                                          ADS21770
      TMP=XZ1                                                           ADS21780
      XZ1=XZ2                                                           ADS21790
      XZ2=TMP                                                           ADS21800
90    IF (XZ3.GT.XZ1) GO TO 100                                         ADS21810
      TMP=XZ1                                                           ADS21820
      XZ1=XZ3                                                           ADS21830
      XZ3=TMP                                                           ADS21840
100   IF (XZ2.LT.XZ3) RETURN                                            ADS21850
      TMP=XZ3                                                           ADS21860
      XZ3=XZ2                                                           ADS21870
      XZ2=TMP                                                           ADS21880
      RETURN                                                            ADS21890
      END                                                               ADS21900
C --- SUBROUTINE ADS105(ICALC,A1,A2,A3,ZRO,XMN1,XMX1)                   ADS21910
C     ------------------------------------------------------------------ADS21920
C --- FIND MINIMUM(S) OF POLYNOMIAL TO 3RD ORDER.                       ADS21930
C     ------------------------------------------------------------------ADS21940
C     BY G. N. VANDERPLAATS                             AUG. 1980       ADS21950
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS21960
C                                                                       ADS21970
C --- INPUT                                                             ADS21980
C     ICALC - CALCULATION CONTROL.                                      ADS21990
C             ICALC=1, FIRST ORDER POLYNOMIAL.  NO MINIMUM OR           ADS22000
C             MAXIMUM EXISTS.  RETURN.                                  ADS22010
C             ICALC=2, SECOND ORDER POLYNOMIAL.  FIND MINIMUM OR        ADS22020
C                      MAXIMUM.                                         ADS22030
C             ICALC=3, THIRD ORDER POLYNOMIAL.  FIND MINIMUM AND MAXIMUMADS22040
C                 IF THEY ARE REAL.                                     ADS22050
C     A1, A2, A3 - COEFFICIENTS OF THE POLYNOMIAL,                      ADS22060
C         Y = AZ + A1*X + A2*X**2 + A3*X**3.                            ADS22070
C     ZRO - SMALL NUMBER TO INDICATE NUMERICAL ZERO.                    ADS22080
C         INDICATOR OF MACHINE ACCURACY.                                ADS22090
C     NOTE - AZ IS NOT NEED IN CALCULATIONS, SO IS NOT BROUGHT          ADS22100
C         THROUGH THE PARAMETER LIST.                                   ADS22110
C                                                                       ADS22120
C --- OUTPUT                                                            ADS22130
C     XMN1, XMX1 - VALUE OF X AT WHICH Y IS MINIMUM AND/OR MAXIMUM.     ADS22140
C         IF ICALC = 1, XMN1 AND XMX1 ARE NOT FOUND.                    ADS22150
C         IF ICALC = 2, XMN1 OR XMX1 IS FOUND, BUT NOT BOTH.            ADS22160
C         IF ICALC = 3, XMN1 AND XMX1 ARE FOUND.                        ADS22170
C                                                                       ADS22180
C --- NOTE.  IF XMN1 OR XMX1 IS NOT FOUND, OR IS IMAGINARY, A           ADS22190
C     VALUE OF -1.0E+20 IS RETURNED.                                    ADS22200
C                                                                       ADS22210
      SUBROUTINE ADS105 (ICALC,A1,A2,A3,ZRO,XMN1,XMX1)                  ADS22220
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS22230
      ZRO10=10.*ZRO                                                     ADS22240
C     INITIALIZE OUTPUT PARAMETERS.                                     ADS22250
      XMN1=-1.0E+20                                                     ADS22260
      XMX1=-1.0E+20                                                     ADS22270
C     ICALC = 1.                                                        ADS22280
C     LINEAR APPROXIMATION.  NO MINIMUM OR MAXIMUM EXISTS.  RETURN.     ADS22290
      IF (ICALC.LT.2) RETURN                                            ADS22300
      IF (ICALC.GT.2) GO TO 20                                          ADS22310
C     ICALC = 2.                                                        ADS22320
C     QUADRATIC APPROXIMATION.  FIND EITHER XMN1 OR XMX1.               ADS22330
10    CONTINUE                                                          ADS22340
C     IF A2 IS NEAR ZERO, NO SOLUTION.                                  ADS22350
      IF (ABS(A2).LT.ZRO) RETURN                                        ADS22360
      XX=-.5*A1/A2                                                      ADS22370
      IF (A2.GT.0.) XMN1=XX                                             ADS22380
      IF (A2.LT.0.) XMX1=XX                                             ADS22390
      RETURN                                                            ADS22400
20    CONTINUE                                                          ADS22410
C     ICALC = 3.                                                        ADS22420
C     IF A3 IS NEAR ZERO, USE QUADRATIC INTERPOLATION.                  ADS22430
      IF (ABS(A3).LT.ZRO10) GO TO 10                                    ADS22440
C     CUBIC APPROXIMATION.  FIND XMN1 AND XMX1.                         ADS22450
      B=A2**2-3.*A1*A3                                                  ADS22460
      IF (B.LT.-1.0E-6) GO TO 30                                        ADS22470
      IF (B.LT.0.) B=0.                                                 ADS22480
      B=SQRT(B)                                                         ADS22490
      XMN1=(-A2+B)/(3.*A3)                                              ADS22500
      XMX1=(-A2-B)/(3.*A3)                                              ADS22510
30    CONTINUE                                                          ADS22520
      RETURN                                                            ADS22530
      END                                                               ADS22540
C --- SUBROUTINE ADS106(N,A,X,Y)                                        ADS22550
C     ------------------------------------------------------------------ADS22560
C --- ROUTINE TO EVALUATE N-TH ORDER POLYNOMIAL, Y=F(X).                ADS22570
C     ------------------------------------------------------------------ADS22580
C     BY G. N. VANDERPLAATS                              AUG. 1980      ADS22590
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS22600
C                                                                       ADS22610
C --- INPUT                                                             ADS22620
C     N - ORDER OF POLYNOMIAL.                                          ADS22630
C     A - ARRAY CONTAINING COEFFICIENTS OF POLYNOMIAL, STARTING         ADS22640
C         WITH A(1)=A0, A(2)=A1, ETC.                                   ADS22650
C         MINIMUM DIMENSION OF A IS N.                                  ADS22660
C     X - POINT AT WHICH THE POLYNOMIAL IS TO BE EVALUATED.             ADS22670
C                                                                       ADS22680
C --- OUTPUT                                                            ADS22690
C     Y - VALUE OF THE POLYNOMIAL AT X.                                 ADS22700
C         Y = A(1) + A(2)*X + A(3)*(X**2) + .... + A(N+1)*(X**N).       ADS22710
C                                                                       ADS22720
      SUBROUTINE ADS106 (N,A,X,Y)                                       ADS22730
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(1)                                                    ADS22740
C     EVALUATE Y.                                                       ADS22750
      Y=A(1)                                                            ADS22760
      IF (N.LT.1) RETURN                                                ADS22770
      XX=X                                                              ADS22780
      DO 10 I=1,N                                                       ADS22790
      IP1=I+1                                                           ADS22800
      Y=Y+XX*A(IP1)                                                     ADS22810
10    XX=XX*X                                                           ADS22820
      RETURN                                                            ADS22830
      END                                                               ADS22840
C --- SUBROUTINE ADS108(ICALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,XMN1,     ADS22850
C    * FMN1)                                                            ADS22860
C     ------------------------------------------------------------------ADS22870
C --- ROUTINE TO FIND THE MINIMUM OF A FUNCTION F(X) BY                 ADS22880
C     POLYNOMIAL INTERPOLATION.                                         ADS22890
C     ------------------------------------------------------------------ADS22900
C     BY G. N. VANDERPLAATS                            DEC. 1980        ADS22910
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS22920
C                                                                       ADS22930
C --- INPUT                                                             ADS22940
C     ICALC - CALCULATION CONTROL                                       ADS22950
C         1 - 2-POINT QUADRATIC INTERPOLATION.                          ADS22960
C         2 - 3-POINT QUADRATIC INTERPOLATION.                          ADS22970
C         3 - 3-POINT CUBIC INTERPOLATION.                              ADS22980
C         4 - 4-POINT CUBIC INTERPOLATION.                              ADS22990
C     X1,F1...X4,F4 - FOUR X-COORDINATES WITH CORRESPONDING F-VALUES.   ADS23000
C     FP1 - SLOPE AT X1.                                                ADS23010
C     ZRO - SMALL POSITIVE NUMBER TO INDICATE NUMERICAL ZERO.           ADS23020
C         INDICATOR OF MACHINE ACCURACY.                                ADS23030
C                                                                       ADS23040
C --- OUTPUT                                                            ADS23050
C     XMN1 - X AT WHICH F IS ESTIMATED TO BE MINIMUM.                   ADS23060
C     FMN1 - POLYNOMIAL APPROXIMATION TO F AT XMN1.                     ADS23070
C                                                                       ADS23080
C --- NOTE.  IF XMN1 CANNOT BE FOUND OR IS IMAGINARY, XMN1 AND FMN1     ADS23090
C     ARE SET TO -1.0E+20.                                              ADS23100
C                                                                       ADS23110
C --- EXTERNALS                                                         ADS23120
C     ADS102, ADS103, ADS105, ADS106.                                   ADS23130
C                                                                       ADS23140
      SUBROUTINE ADS108 (ICALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,XMN1,FMN1ADS23150
     1)                                                                 ADS23160
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(4)                                                    ADS23170
      IF (ICALC.GT.2) GO TO 10                                          ADS23180
C     QUADRATIC INTERPOLATION.                                          ADS23190
C     COEFFICIENTS.                                                     ADS23200
      CALL ADS102 (ICALC,X1,F1,FP1,X2,F2,X3,F3,AZ,A1,A2)                ADS23210
C     GO FIND MINIMUM.                                                  ADS23220
      A3=0.                                                             ADS23230
      JCALC=2                                                           ADS23240
      GO TO 20                                                          ADS23250
10    CONTINUE                                                          ADS23260
C     CUBIC INTERPOLATION.                                              ADS23270
C     COEFFICIENTS.                                                     ADS23280
      JCALC=ICALC-2                                                     ADS23290
      CALL ADS103 (JCALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,AZ,A1,A2,A3)   ADS23300
C     GO FIND MINIMUM.                                                  ADS23310
      JCALC=3                                                           ADS23320
20    CALL ADS105 (JCALC,A1,A2,A3,ZRO,XMN1,XMX1)                        ADS23330
C     CHECK FOR SOLUTION.  IF NO, RETURN.                               ADS23340
      FMN1=-1.0E+20                                                     ADS23350
      IF (XMN1.LT.-0.9E+20) RETURN                                      ADS23360
C     EVALUATE POLYNOMIAL APPROXIMATION TO FMN1.                        ADS23370
      N=3                                                               ADS23380
      IF (ICALC.LT.3) N=2                                               ADS23390
      A(1)=AZ                                                           ADS23400
      A(2)=A1                                                           ADS23410
      A(3)=A2                                                           ADS23420
      A(4)=A3                                                           ADS23430
      CALL ADS106 (N,A,XMN1,FMN1)                                       ADS23440
      RETURN                                                            ADS23450
      END                                                               ADS23460
C --- SUBROUTINE ADS109(ICALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,          ADS23470
C    * XZ1,Z1,ZP1,XZ2,Z2,ZP2,XZ3,Z3,ZP3,IER)                            ADS23480
C     DIMENSION AA(4)                                                   ADS23490
C     ------------------------------------------------------------------ADS23500
C --- ROUTINE TO PERFORM POLYNOMIAL INTERPOLATION UP TO 3RD-ORDER       ADS23510
C     FOR THE ZEROES OF A FUNCTION.                                     ADS23520
C     ------------------------------------------------------------------ADS23530
C     BY G. N. VANDERPLAATS                              JULY 1981      ADS23540
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS23550
C                                                                       ADS23560
C --- INPUT                                                             ADS23570
C     ICALC - CALCULATION CONTROL.                                      ADS23580
C         1 = 1-POINT LINEAR INTERPOLATION.                             ADS23590
C         2 = 2-POINT LINEAR INTERPOLATION.                             ADS23600
C         3 = 2-POINT QUADRATIC INTERPOLATION.                          ADS23610
C         4 = 3-POINT QUADRATIC INTERPOLATION.                          ADS23620
C         5 = 3-POINT CUBIC INTERPOLATION.                              ADS23630
C         6 = 4-POINT CUBIC INTERPOLATION.                              ADS23640
C     X1,F1,X2,F2,X3,F3,X4,F4 - X-F PAIRS OF DATA.                      ADS23650
C     FP1 - SLOPE AT X1.                                                ADS23660
C     ZRO - SMALL POSITIVE NUMBER TO INDICATE NUMERICAL ZERO.           ADS23670
C         INDICATOR OF MACHINE ACCURACY.                                ADS23680
C --- NOTE                                                              ADS23690
C     ONLY SUFFICIENT DATA NEED BE GIVEN TO PERFORM THE CALCULATIONS    ADS23700
C         REQUIRED BY ICALC.  REMAINING DATA MAY BE INPUT AS 0.0.       ADS23710
C                                                                       ADS23720
C --- OUTPUT                                                            ADS23730
C     XZ1,Z1,ZP1....XZ3,Z3,ZP3 - PROPOSED X-VALUES AT WHICH F=0,        ADS23740
C         AND THEIR CORRESPONDING FUNCTION VALUES AND SLOPES (BASED     ADS23750
C         ON THE POLYNOMIAL INTERPOLATION).                             ADS23760
C     IER - ERROR CODE.  IF IER=1, NEWTON'S METHOD DID NOT CONVERGE     ADS23770
C         IN 50 ITERATIONS.  THE RESULTS MAY NOT BE VALID.              ADS23780
C --- NOTE                                                              ADS23790
C     IF ICALC=1 OR 2, ONLY Z1 AND ZP1 ARE OUTPUT.                      ADS23800
C     IF ICALC=3 OR 4, Z1, ZP1, Z2 AND ZP2 ARE OUTPUT.                  ADS23810
C     IF ICALC=5 OR 6, Z1 THROUGH ZP3 ARE OUTPUT.                       ADS23820
C                                                                       ADS23830
C --- EXTERNALS                                                         ADS23840
C     ADS101, ADS102, ADS103, ADS104.                                   ADS23850
C                                                                       ADS23860
      SUBROUTINE ADS109 (ICALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,XZ1,Z1,ZPADS23870
     11,XZ2,Z2,ZP2,XZ3,Z3,ZP3,IER)                                      ADS23880
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION AA(4)                                                   ADS23890
      A2=0.                                                             ADS23900
      A3=0.                                                             ADS23910
      GO TO (10,10,20,20,30,30), ICALC                                  ADS23920
10    CONTINUE                                                          ADS23930
C     LINEAR INTERPOLATION.                                             ADS23940
      JCALC=ICALC                                                       ADS23950
      CALL ADS101 (JCALC,X1,F1,FP1,X2,F2,AZ,A1)                         ADS23960
      GO TO 40                                                          ADS23970
C     QUADRATIC INTERPOLATION.                                          ADS23980
20    CONTINUE                                                          ADS23990
      JCALC=ICALC-2                                                     ADS24000
      CALL ADS102 (JCALC,X1,F1,FP1,X2,F2,X3,F3,AZ,A1,A2)                ADS24010
      GO TO 40                                                          ADS24020
30    CONTINUE                                                          ADS24030
C     CUBIC INTERPOLATION.                                              ADS24040
      JCALC=ICALC-4                                                     ADS24050
      CALL ADS103 (JCALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,AZ,A1,A2,A3)   ADS24060
40    CONTINUE                                                          ADS24070
C     FIND ZEROES.                                                      ADS24080
      JCALC=(ICALC+1)/2                                                 ADS24090
      CALL ADS104 (JCALC,AZ,A1,A2,A3,ZRO,XZ1,XZ2,XZ3,IER)               ADS24100
C     EVALUATE SLOPE OF POLYNOMIAL AT ZEROES.                           ADS24110
      Z1=-1.0E+20                                                       ADS24120
      Z2=-1.0E+20                                                       ADS24130
      Z3=-1.0E+20                                                       ADS24140
      ZP1=-1.0E+20                                                      ADS24150
      ZP2=-1.0E+20                                                      ADS24160
      ZP3=-1.0E+20                                                      ADS24170
      AA(1)=AZ                                                          ADS24180
      AA(2)=A1                                                          ADS24190
      AA(3)=A2                                                          ADS24200
      AA(4)=A3                                                          ADS24210
      N=(ICALC+1)/2                                                     ADS24220
      IF (XZ1.GT.-1.0E+20) CALL ADS106 (N,AA,XZ1,Z1)                    ADS24230
      IF (XZ2.GT.-1.0E+20) CALL ADS106 (N,AA,XZ2,Z2)                    ADS24240
      IF (XZ3.GT.-1.0E+20) CALL ADS106 (N,AA,XZ3,Z3)                    ADS24250
      IF (XZ1.GT.-1.0E+20) CALL ADS110 (N,AA,XZ1,ZP1)                   ADS24260
      IF (XZ2.GT.-1.0E+20) CALL ADS110 (N,AA,XZ2,ZP2)                   ADS24270
      IF (XZ3.GT.-1.0E+20) CALL ADS110 (N,AA,XZ3,ZP3)                   ADS24280
      RETURN                                                            ADS24290
      END                                                               ADS24300
C --- SUBROUTINE ADS110(N,A,X,YP)                                       ADS24310
C     DIMENSION A(N+1)                                                  ADS24320
C     ------------------------------------------------------------------ADS24330
C --- ROUTINE TO EVALUATE SLOPE OF N-TH ORDER POLYNOMIAL, Y=F(X).       ADS24340
C     ------------------------------------------------------------------ADS24350
C     BY G. N. VANDERPLAATS                              AUG. 1980      ADS24360
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS24370
C                                                                       ADS24380
C --- INPUT                                                             ADS24390
C     N - ORDER OF POLYNOMIAL.                                          ADS24400
C     A - ARRAY CONTAINING COEFFICIENTS OF POLYNOMIAL, STARTING         ADS24410
C         WITH A(1)=A0, A(2)=A1, ETC.                                   ADS24420
C         MINIMUM DIMENSION OF A IS N.                                  ADS24430
C     X - POINT AT WHICH THE POLYNOMIAL IS TO BE EVALUATED.             ADS24440
C                                                                       ADS24450
C --- OUTPUT                                                            ADS24460
C     YP - VALUE OF THE SLOPE OF THE POLYNOMIAL AT X.                   ADS24470
C      YP = A(2) + 2*A(3)*X + 3*A(4)*(X**2) + ... + N*A(N+1)*(X**(N-1)).ADS24480
C                                                                       ADS24490
      SUBROUTINE ADS110 (N,A,X,YP)                                      ADS24500
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(1)                                                    ADS24510
C     EVALUATE YP.                                                      ADS24520
      YP=0.                                                             ADS24530
      IF (N.LT.1) RETURN                                                ADS24540
      YP=A(2)                                                           ADS24550
      IF (N.LT.2) RETURN                                                ADS24560
      XX=X                                                              ADS24570
      DO 10 I=2,N                                                       ADS24580
      IP1=I+1                                                           ADS24590
      YP=YP+FLOAT(I)*A(IP1)*XX                                          ADS24600
10    XX=XX*X                                                           ADS24610
      RETURN                                                            ADS24620
      END                                                               ADS24630
C --- SUBROUTINE ADS111 (NCON,NX,ISLOPG,CTMIN,CTLMIN,ZRO,X1,X2,X3,X4,   ADS24640
C    1DG,GGL,GG1,GG2,GGU,IDG,ALPL,ALPU)                                 ADS24650
C     ------------------------------------------------------------------ADS24660
C --- POLYNOMIAL INTERPOLATION FOR CONSTRAINT BOUNDARIES.               ADS24670
C     ------------------------------------------------------------------ADS24680
C     BY G. N. VANDERPLAATS                                 ***** 1983  ADS24690
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS24700
      SUBROUTINE ADS111 (NCON,NX,ISLOPG,CTMIN,CTLMIN,ZRO,X1,X2,X3,X4,DG,ADS24710
     1GGL,GG1,GG2,GGU,IDG,ALPL,ALPU)                                    ADS24720
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION DG(1),GGL(1),GG1(1),GG2(1),GGU(1),IDG(1)                ADS24730
      ALPL1=ALPL                                                        ADS24740
      ALPU1=ALPU                                                        ADS24750
C     BEGIN INTERPOLATION.                                              ADS24760
      DO 40 I=1,NCON                                                    ADS24770
C     IS THIS CONSTRAINT INITIALLY ACTIVE OR VIOLATED?                  ADS24780
      FP1=0.                                                            ADS24790
      IAV=0                                                             ADS24800
      IF (ISLOPG.EQ.0) GO TO 10                                         ADS24810
      FP1=DG(I)                                                         ADS24820
      IF (ABS(FP1).GT.0.01) IAV=I                                       ADS24830
10    CONTINUE                                                          ADS24840
C     CONSTRAINT VIOLATION THRESHOLD.                                   ADS24850
      CC=CTMIN                                                          ADS24860
      IF (IABS(IDG(I)).GT.1) CC=CTLMIN                                  ADS24870
C     FUNCTION VALUES.                                                  ADS24880
      F1=GGL(I)                                                         ADS24890
      IF (IABS(IDG(I)).GT.1.AND.ABS(F1).LE.CC) GO TO 40                 ADS24900
      F2=GG1(I)                                                         ADS24910
      F3=0.                                                             ADS24920
      F4=0.                                                             ADS24930
      IF (NX.GE.3) F3=GG2(I)                                            ADS24940
      IF (NX.GE.4) F4=GGU(I)                                            ADS24950
C     IGNORE THIS CONSTRAINT IF ALL ARE WITHIN A TOLORANCE CC OF ZERO.  ADS24960
      AF1=ABS(F1)                                                       ADS24970
      AF2=ABS(F2)                                                       ADS24980
      AF3=ABS(F3)                                                       ADS24990
      AF4=ABS(F4)                                                       ADS25000
      AF=DMAX1(AF1,AF2,AF3,AF4)                                         ADS25010
      IF (AF.LE.CC) GO TO 40                                            ADS25020
      IF(IAV.EQ.0.AND.NX.LT.2) GO TO 40                                 ADS25030
      ICALC=2*(NX-1)                                                    ADS25040
      IF (ICALC.GT.6) ICALC=6                                           ADS25050
      IF (IAV.GT.0.AND.NX.LT.4) ICALC=ICALC+1                           ADS25060
      CALL ADS109 (ICALC,X1,F1,FP1,X2,F2,X3,F3,X4,F4,ZRO,XZ1,Z1,ZP1,XZ2,ADS25070
     1Z2,ZP2,XZ3,Z3,ZP3,IER)                                            ADS25080
      DO 30 II=1,3                                                      ADS25090
      XZ=XZ1                                                            ADS25100
      ZP=ZP1                                                            ADS25110
      IF (II.EQ.2) XZ=XZ2                                               ADS25120
      IF (II.EQ.2) ZP=ZP2                                               ADS25130
      IF (II.EQ.3) XZ=XZ3                                               ADS25140
      IF (II.EQ.3) ZP=ZP3                                               ADS25150
      IF (ABS(XZ).GT.0.9E+20) GO TO 30                                  ADS25160
      IF (F1.LE.CC) GO TO 20                                            ADS25170
C     INITIALLY INFEASIBLE CASE, BUT DECREASING.                        ADS25180
      IF (XZ.GT.ALPL.AND.ZP.LT.0.) ALPL=XZ                              ADS25190
      GO TO 30                                                          ADS25200
20    CONTINUE                                                          ADS25210
      IF (XZ.LT.ALPL1) GO TO 30                                         ADS25220
C     INITIALLY FEASIBLE, AND F1.LT.0.                                  ADS25230
      IF (XZ.LT.ALPU.AND.ZP.GT.0.) ALPU=XZ                              ADS25240
      IF (FP1.GT.ZRO) GO TO 30                                          ADS25250
C     INITIALLY FEASIBLE.  SECOND INTERCEPT.                            ADS25260
      IF (XZ.LT.ALPU.AND.(F1.LT.0.0.AND.ZP.GT.0.)) ALPU=XZ              ADS25270
30    CONTINUE                                                          ADS25280
40    CONTINUE                                                          ADS25290
      IF (ALPL.GT.ALPU1) ALPL=ALPU1                                     ADS25300
      RETURN                                                            ADS25310
      END                                                               ADS25320
C --- SUBROUTINE ADS112 (NDV,X,S,SVLB,SVUB,ZRO,ALPMAX)                  ADS25330
C     DIMENSION X(NDV),S(NDV),SVLB(NDV),SVUB(NDV)                       ADS25340
C     ------------------------------------------------------------------ADS25350
C --- ROUTINE TO FIND ALPMAX SO NDV SIDE CONSTRAINTS ARE ACTIVE.        ADS25360
C     ------------------------------------------------------------------ADS25370
C     BY G. N. VANDERPLAATS                             AUG. 1982       ADS25380
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS25390
C                                                                       ADS25400
C --- INPUT                                                             ADS25410
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS25420
C     X - CURRENT VECTOR OF DESIGN VARIABLES.  MINIMUM DIMENSION IS NDV.ADS25430
C     S - SEARCH DIRECTION.  MINIMUM DIMENSION IS NDV.                  ADS25440
C     SVLB - LOWER BOUNDS ON X.  MINIMUM DIMENSION IS NDV.              ADS25450
C     SVUB - UPPER BOUNDS ON X.  MINIMUM DIMENSION IS NDV.              ADS25460
C     ZRO - SMALL POSITIVE NUMBER INDICATING NUMERICAL ZERO.            ADS25470
C                                                                       ADS25480
C --- OUTPUT                                                            ADS25490
C     ALPMAX - MOVE, ALPHA SUCH THAT NDV SIDE CONSTRAINTS ARE ACTIVE.   ADS25500
C                                                                       ADS25510
      SUBROUTINE ADS112 (NDV,X,S,SVLB,SVUB,ZRO,ALPMAX)                  ADS25520
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),S(1),SVLB(1),SVUB(1)                               ADS25530
C                                                                       ADS25540
      ALPMAX=0.                                                         ADS25550
      DO 10 I=1,NDV                                                     ADS25560
      IF (ABS(S(I)).LE.ZRO) GO TO 10                                    ADS25570
      ALP=(SVUB(I)-X(I))/S(I)                                           ADS25580
      IF (ALP.GT.ALPMAX) ALPMAX=ALP                                     ADS25590
      ALP=(SVLB(I)-X(I))/S(I)                                           ADS25600
      IF (ALP.GT.ALPMAX) ALPMAX=ALP                                     ADS25610
10    CONTINUE                                                          ADS25620
      RETURN                                                            ADS25630
      END                                                               ADS25640
C --- SUBROUTINE ADS115(IGOTO,IONED,IPRNT1,NDV,NCON,OBJ,SLOPF,ISLOPF,   ADS25650
C    *ISLOPG,CTMIN,CTLMIN,ALPMAX,XMULT,EXTRAP,ZRO,DELALP,DABALP,ALPL,   ADS25660
C    * ALPU,XOLD,SVLB,SVUB,S,X,G,DG,IDG,GGL,GG1,GG2,XTMP,XBST,NAVC,IC)  ADS25670
C     DIMENSION X(NDV),S(NDV),XOLD(NDV),SVLB(NDV),SVUB(NDV),G(NCON),    ADS25680
C    * DG(NACG),IDG(NCON),GGL(NCON),GG1(NCON),GG2(NCON),XTMP(NDV),      ADS25690
C    * XBST(NDV),IC(NAVC)                                               ADS25700
C     ------------------------------------------------------------------ADS25710
C --- ROUTINE TO CONTROL ONE-DIMENSIONAL SEARCH FOR CONSTRAINED         ADS25720
C     FUNCTIONS.                                                        ADS25730
C     ------------------------------------------------------------------ADS25740
C     BY G. N. VANDERPLAATS                          MAY, 1982.         ADS25750
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS25760
C                                                                       ADS25770
C --- INPUT.                                                            ADS25780
C     IGOTO - POINTER.  INITIALLY INPUT AS ZERO.  ON SUBSEQUENT CALLS,  ADS25790
C         IGOTO = VALUE RETURNED FROM LAST CALL.                        ADS25800
C     IONED - METHOD                                                    ADS25810
C         5 - USE GOLDEN SECTION METHOD.                                ADS25820
C         6 - GOLDEN SECTION FOLLOWED BY POLYNOMIAL INTERPOLATION.      ADS25830
C         7 - POLYNOMIAL INTERPOLATION WITH BOUNDS.                     ADS25840
C         8 - POLYNOMIAL INTERPOLATION/EXTRAPOLATION WITHOUT FIRST      ADS25850
C             GETTING BOUNDS.                                           ADS25860
C     IPRNT1 - PRINT CONTROL.                                           ADS25870
C         0 - NO PRINT.                                                 ADS25880
C         2 - PRINT TITLE, METHOD, EACH PROPOSED ALPHA, OBJ AND GMAX.   ADS25890
C         3 - PRINT ALL OF ABOVE PLUS CONTROL PARAMETERS, S-VECTOR,     ADS25900
C             X-VECTOR AND G-VECTOR.                                    ADS25910
C         NOTE; IPRNT1 IS USUALLY USED FOR DEBUG PRINTING ONLY.         ADS25920
C     NDV - NUMBER OF DESIGN VARIABLES CONTAINED IN X.                  ADS25930
C     NCON - NUMBER OF CONSTRAINT VALUES STORED IN G-VECTOR.            ADS25940
C     OBJ - ON FIRST CALL, OBJ = OBJECTIVE FUNCTION CORRESPONDING TO    ADS25950
C         ALPHA = 0.  ON SUBSEQUENT CALLS, OBJ = OBJECTIVE FUNCTION     ADS25960
C         CORRESPONDING TO X-VECTOR RETURNED FROM LAST CALL.            ADS25970
C     SLOPF - SLOPF OF OBJ.VS.ALPHA CURVE.  NOT USED IF ISLOPF=0.       ADS25980
C     ISLOPF _ TELLS IF SLOPE OF OBJ.VS.ALPHA CURVE IS AVAILABLE.  IF   ADS25990
C         ISLOPF = 1, SLOPF IS AVAILABLE FOR USE IN POLYNOMIAL          ADS26000
C         INTERPOLATION.                                                ADS26010
C     ISLOPG - TELLS IF SLOPE OF G(I).VS.ALPHA IS AVAILABLE.  IF        ADS26020
C         ISLOPG = 1, DG CONTAINS THE DESIRED GRADIENTS.                ADS26030
C     CTMIN - NONLINEAR CONSTRAINT TOLORANCE.  IF G(I).GT.CTMIN,        ADS26040
C         CONSTRAINT I IS VIOLATED.                                     ADS26050
C     CTLMIN - LINEAR CONSTRAINT TOLORANCE.  IF G(I).GT.CTLMIN,         ADS26060
C         CONSTRAINT IS VIOLATED.                                       ADS26070
C     ALPHA - ON FIRST CALL, ALPHA = PROPOSED MOVE PARAMETER.  ON       ADS26080
C         SUBSEQUENT CALLS, ALPHA = ALPHA RETURNED FROM PREVIOUS CALL.  ADS26090
C     ALPMAX - MAXIMUM ALPHA THAT WILL BE CONSIDERED.                   ADS26100
C     XMULT - MULTIPLIER FOR FINDING BRACKETS ON THE MINIMUM OF OBJ.    ADS26110
C         XMULT = 1.5 + SQRT(5)/2 IS RECOMMENDED, SINCE THIS WILL GIVE  ADS26120
C         ONE INTERIOR POINT FOR THE GOLDEN SECTION METHOD.             ADS26130
C     EXTRAP - THE FACTOR BY WHICH ALPHA MAY BE EXTRAPOLATED IN         ADS26140
C         POLYNOMIAL FIT WITHOUT FIRST GETTING BOUNDS.                  ADS26150
C     ZRO - IDENTIFIES NUMERICAL ZERO, SAY 1.0E-5.  A MEASURE OF        ADS26160
C         MACHINE ACCURACY.                                             ADS26170
C     DELALP - RELATIVE CONVERGENCE CRITERION FOR GOLDEN SECTION.       ADS26180
C     DABALP - ABSOLUTE CONVERGENCE CRITERION FOR GOLDEN SECTION.       ADS26190
C     XOLD - VECTOR FOR STORAGE OF DESIGN VARIABLES AT ALPHA = 0.       ADS26200
C         MINIMUM DIMENSION IS NDV.                                     ADS26210
C     SVLB - LOWER BOUNDS ON THE VARIABLES.  IF ANY COMPONENT OF        ADS26220
C         X VIOLATES SVLB, IT IS SET EQUAL TO ITS LOWER BOUND.  MINIMUM ADS26230
C         DIMENSION IS NDV.                                             ADS26240
C     SVUB - UPPER BOUNDS ON THE VARIABLES.  IF ANY COMPONENT OF        ADS26250
C         X VIOLATES SVUB, IT IS SET EQUAL TO ITS UPPER BOUND.  MINIMUM ADS26260
C         DIMENSION IS NDV.                                             ADS26270
C     S - SEARCH DIRECTION.  MINIMUM DIMENSION IS NDV.                  ADS26280
C     G - CONSTRAINT VECTOR.  MINIMUM DIMENSION IS NCON.                ADS26290
C     DG - GRADIENTS OF CURRENTLY ACTIVE/VIOLATED CONSTRAINTS.          ADS26300
C         MINIMUM DIMENSION IS NACG.  IF THIS INFORMATION IS NOT TO BE  ADS26310
C         USED, NACG SHOULD BE DEFINED AS ZERO, EVEN THOUGH CONSTRAINTS ADS26320
C         MAY BE ACTIVE/VIOLATED.                                       ADS26330
C     IDG - CONSTRAINT TYPE IDENTIFIER.  MINIMUM DIMENSION IS NCON.     ADS26340
C     GGL, GG1, GG2 - TEMPORARY STORAGE OF G-VECTOR.  MINIMUM DIMENSION ADS26350
C         FOR EACH IS NCON.                                             ADS26360
C     XTMP - TEMPORARY STORAGE OF THE X-VECTOR.  MINIMUM DIMENSION IS   ADS26370
C         NDV.                                                          ADS26380
C     NAVC - NUMBER OA ACTIVE AND VIOLATED CONSTRAINTS AT ALPHA=0.      ADS26390
C     IC - VECTOR IDENTIFYING ACTIVE AND VIOLATED CONSTRAINTS.          ADS26400
C     XBST - TEMPORARY STORAGE OF BEST X-VECTOR.                        ADS26410
C                                                                       ADS26420
C --- OUTPUT.                                                           ADS26430
C     IGOTO - POINTER.  IF IGOTO.GT.0, EVALUATE OBJ AND CALL ADS115     ADS26440
C         AGAIN.  IF IGOTO = 0, THE ONE-DIMENSIONAL SEARCH IS COMPLETE. ADS26450
C     OBJ - MINIMUM OBJECTIVE FUNCTION VALUE WHEN IGOTO = 0.            ADS26460
C     ALPHA - CURRENT MOVE PARAMETER IF IGOTO.GT.0.  OPTIMUM MOVE       ADS26470
C         PARAMETER IF IGOTO = 0.                                       ADS26480
C     ALPL - BEST KNOWN LOWER BOUND ON ALPHA.                           ADS26490
C     ALPU - BEST KNOWN UPPER BOUND ON ALPHA.                           ADS26500
C     X - VECTOR OF DESIGN VARIABLES FOR THE CURRENT VALUE OF ALPHA     ADS26510
C         IF IGOTO.GT.0.  IN THIS CASE EVALUATE OBJ AND CALL AGAIN.     ADS26520
C         VECTOR OF DESIGN VARIABLES CORRESPONDING TO THE OPTIMUM       ADS26530
C         ALPHA IF IGOGO = 0.  MINIMUM DIMENSION IS NDV.                ADS26540
C     G - VECTOR OF CONSTRAINT VALUES ASSOCIATED WITH OPTIMUM ALPHA.    ADS26550
C                                                                       ADS26560
C --- EXTERNALS ADS102, ADS103, ADS105, ADS503, ADS206, ADS121,         ADS26570
C         ADS126.                                                       ADS26580
      SUBROUTINE ADS115 (IGOTO,IONED,IPRNT1,NDV,NCON,OBJ,SLOPF,ISLOPF,ISADS26590
     1LOPG,CTMIN,CTLMIN,ALPHA,ALPMAX,XMULT,EXTRAP,ZRO,DELALP,DABALP,ALPLADS26600
     2,ALPU,XOLD,SVLB,SVUB,S,X,G,DG,IDG,GGL,GG1,GG2,GGU,XTMP,XBST,NAVC,IADS26610
     3C,ALPBST,FBST,FL,FU,F1,F2,GBST,GMAX,GMCT,GMCTL,GMCTU,GMCT1,GMCT2, ADS26620
     4ICAL,II,JGOTO,JJ,II1,IX,NITER,SCFX,SCFG,ISTRAT,TEMP)              ADS26630
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),S(1),XOLD(1),SVLB(1),SVUB(1),G(1),DG(1),IDG(1),GGL(ADS26640
     11),GG1(1),GG2(1),GGU(1),XTMP(1),XBST(1),IC(1),SCFX(1),SCFG(1),    ADS26650
     2TEMP(10)                                                          ADS26660
      XL=TEMP(1)                                                        ADS26670
      XU=TEMP(2)                                                        ADS26680
      X1=TEMP(3)                                                        ADS26690
      X2=TEMP(4)                                                        ADS26700
      SCFO=TEMP(5)                                                      ADS26710
C                                                                       ADS26720
C     FIND MAXIMUM CONSTRAINT TOLORANCE VALUE.                          ADS26730
      CALL ADS502 (NCON,G,IDG,CTMIN,CTLMIN,GMCT,IGMAX)                  ADS26740
      IEVAL=1                                                           ADS26750
      IF (IGOTO.GT.0) GO TO 40                                          ADS26760
C     IGOTO = 0, FIRST CALL.                                            ADS26770
      X1=-10.*ZRO                                                       ADS26780
      X2=-10.*ZRO                                                       ADS26790
      XU=-10.*ZRO                                                       ADS26800
      FBST=OBJ                                                          ADS26810
      ALPBST=0.                                                         ADS26820
      GBST=GMCT                                                         ADS26830
      ICAL=0                                                            ADS26840
C     SAVE X-VECTOR IN XOLD AND XBST.                                   ADS26850
      CALL ADS503 (NDV,X,XOLD)                                          ADS26860
      CALL ADS503 (NDV,X,XBST)                                          ADS26870
      IF (IPRNT1.GT.1) WRITE (6,360) SLOPF,ALPHA,ALPMAX,XMULT,EXTRAP,ZROADS26880
     1,DELALP,DABALP,ISLOPF,ISLOPG                                      ADS26890
      ICALC=13                                                          ADS26900
      IF (IPRNT1.GT.1) CALL ADS402 (ICALC,DG,NCON,TMP1)                 ADS26910
40    IGOTO1=IGOTO+1                                                    ADS26920
      ICAL=ICAL+1                                                       ADS26930
      IF (ICAL.EQ.1) GO TO 80                                           ADS26940
C     ------------------------------------------------------------------ADS26950
C     IDENTIFY BEST DESIGN.                                             ADS26960
C     ------------------------------------------------------------------ADS26970
C     IS THIS DESIGN LESS INFEASIBLE?                                   ADS26980
      IF (GBST.GT.0.0.AND.GMCT.LT.GBST) GO TO 50                        ADS26990
C     NO.  IS THIS DESIGN FEASIBLE WITH LOWER OBJECTIVE?                ADS27000
      IF (OBJ.LT.FBST.AND.GMCT.LE.0.) GO TO 50                          ADS27010
C     NO. THIS IS NOT THE BEST DESIGN.                                  ADS27020
      GO TO 80                                                          ADS27030
C     SAVE THIS DESIGN.                                                 ADS27040
50    FBST=OBJ                                                          ADS27050
      ALPBST=ALPHA                                                      ADS27060
      GBST=GMCT                                                         ADS27070
      CALL ADS503 (NDV,X,XBST)                                          ADS27080
      IF (GBST.GT.0.) GO TO 80                                          ADS27090
C     IF THIS DESIGN ENCOUNTERED A NEW CONSTRAINT AND IT IS ACTIVE      ADS27100
C     WITHIN A TOLORANCE OF CTMIN (CTLMIN) THEN END ONE-D SEARCH.       ADS27110
      DO 70 I=1,NCON                                                    ADS27120
      I1=IDG(I)                                                         ADS27130
      CC=CTMIN                                                          ADS27140
      IF (IABS(I1).EQ.1) CC=CTLMIN                                      ADS27150
      CC=ABS(G(I))-CC                                                   ADS27160
      IF (CC.GT.0.) GO TO 70                                            ADS27170
C     ACTIVE CONSTRAINT.  IS IT NEW?                                    ADS27180
      IF (NAVC.EQ.0) GO TO 290                                          ADS27190
      DO 60 J=1,NAVC                                                    ADS27200
      IF (IC(J).EQ.I) GO TO 70                                          ADS27210
60    CONTINUE                                                          ADS27220
C     NEW CONSTRAINT.  GO PICK BEST DESIGN.                             ADS27230
      GO TO 290                                                         ADS27240
70    CONTINUE                                                          ADS27250
80    CONTINUE                                                          ADS27260
      GO TO (90,100,120,290,200,270), IGOTO1                            ADS27270
C     IGOTO = 0.  FIRST CALL.                                           ADS27280
90    CONTINUE                                                          ADS27290
C     DO WE NEED TO BRACKET THE MINIMUM?                                ADS27300
      IF (IONED.GT.7) GO TO 260                                         ADS27310
C     YES.                                                              ADS27320
C     ------------------------------------------------------------------ADS27330
C     FIND BRACKETS ON MINIMUM.                                         ADS27340
C     ------------------------------------------------------------------ADS27350
      JGOTO=0                                                           ADS27360
      II=1                                                              ADS27370
      IF (IPRNT1.GT.0) WRITE (6,370)                                    ADS27380
100   CONTINUE                                                          ADS27390
      CALL ADS125 (JGOTO,NCON,OBJ,GMCT,ALPHA,XMULT,ALPMAX,XL,X1,X2,XU,FLADS27400
     1,F1,F2,FU,GMCTL,GMCT1,GMCT2,GMCTU,G,GGL,GG1,GG2,GGU)              ADS27410
      II=II+1                                                           ADS27420
      IF (JGOTO.EQ.0) GO TO 110                                         ADS27430
C     UPDATE X.                                                         ADS27440
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS27450
C     EVALUATE OBJECTIVE AND CONSTRAINT FUNCTIONS AND CALL AGAIN.       ADS27460
      IGOTO=1                                                           ADS27470
      TEMP(1)=XL                                                        ADS27480
      TEMP(2)=XU                                                        ADS27490
      TEMP(3)=X1                                                        ADS27500
      TEMP(4)=X2                                                        ADS27510
      TEMP(5)=SCFO                                                      ADS27520
      RETURN                                                            ADS27530
110   CONTINUE                                                          ADS27540
C     IF ONLY TWO DESIGNS ARE AVAILABLE, PUT G IN GGU.                  ADS27550
      IF(II.EQ.2) CALL ADS503(NCON,G,GGU)                               ADS27560
      IF (IPRNT1.GT.1) WRITE (6,380) XL,XU                              ADS27570
C     ------------------------------------------------------------------ADS27580
C     MINIMUM IS BRACKETED.  REFINE SOLUTION AS REQUIRED.               ADS27590
C     ------------------------------------------------------------------ADS27600
      IF (IONED.GT.6) GO TO 170                                         ADS27610
C     ------------------------------------------------------------------ADS27620
C     IONED = 5 OR 6.  REFINE SOLUTION USING THE GOLDEN SECTION METHOD. ADS27630
C     ------------------------------------------------------------------ADS27640
      JGOTO=0                                                           ADS27650
      IF (IPRNT1.GT.0) WRITE (6,390)                                    ADS27660
      NPTS=ICAL                                                         ADS27670
      IF (NPTS.GT.4) NPTS=4                                             ADS27680
      XX=XMULT-2.618034                                                 ADS27690
      IF (ABS(XX).GT.0.0001) NPTS=2                                     ADS27700
      XX=XU-ALPMAX                                                      ADS27710
      IF (ABS(XX).LT.0.0001) NPTS=2                                     ADS27720
      IF (NPTS.NE.3) GO TO 120                                          ADS27730
      X1=X2                                                             ADS27740
      F1=F2                                                             ADS27750
      GMCT1=GMCT2                                                       ADS27760
      CALL ADS503 (NCON,GG2,GG1)                                        ADS27770
120   CONTINUE                                                          ADS27780
      CALL ADS123 (JGOTO,NCON,ICAL,NPTS,OBJ,GMCT,ALPHA,ALPMAX,DELALP,DABADS27790
     1ALP,ZRO,XL,X1,X2,XU,FL,F1,F2,FU,GMCTL,GMCT1,GMCT2,GMCTU,IDL,ID1,IDADS27800
     22,IDU,G,GGL,GG1,GG2,GGU,II1,IX,NITER)                             ADS27810
C     FINISHED?                                                         ADS27820
      IF (JGOTO.EQ.0) GO TO 130                                         ADS27830
C     NO. UPDATE X.                                                     ADS27840
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS27850
C     EVALUATE OBJECTIVE FUNCTION AND CALL AGAIN.                       ADS27860
      IGOTO=2                                                           ADS27870
      TEMP(1)=XL                                                        ADS27880
      TEMP(2)=XU                                                        ADS27890
      TEMP(3)=X1                                                        ADS27900
      TEMP(4)=X2                                                        ADS27910
      TEMP(5)=SCFO                                                      ADS27920
      RETURN                                                            ADS27930
130   CONTINUE                                                          ADS27940
      IF (IPRNT1.GT.0) WRITE (6,400)                                    ADS27950
C     IF IONED = 5, SOLUTION IS COMPLETE.  GO FIND BEST DESIGN.         ADS27960
C     IF IONED = 6, DO POLYNOMIAL INTERPOLATION.                        ADS27970
      IGOTO=0                                                           ADS27980
      IF (IONED.EQ.5) GO TO 290                                         ADS27990
C     ------------------------------------------------------------------ADS28000
C     IONED = 6.  POLYNOMIAL INTERPOLATION USING LAST FOUR GOLDEN       ADS28010
C     SECTION POINTS.                                                   ADS28020
C     ------------------------------------------------------------------ADS28030
      IF (IPRNT1.GT.0) WRITE (6,410)                                    ADS28040
      GMIN=DMIN1(GMCTL,GMCT1,GMCT2,GMCTU)                               ADS28050
      GMAX=DMAX1(GMCTL,GMCT1,GMCT2,GMCTU)                               ADS28060
      ALPL=XL                                                           ADS28070
      ALPU=XU                                                           ADS28080
      ALPHA=ALPU                                                        ADS28090
      IF (GMIN.GT.0.0.OR.GMAX.LT.0.) GO TO 140                          ADS28100
C     SOME DESIGNS ARE FEASIBLE AND SOME INFEASIBLE.  UPDATE BOUNDS.    ADS28110
      II=4                                                              ADS28120
      CALL ADS111 (NCON,II,ISLOPG,CTMIN,CTLMIN,ZRO,XL,X1,X2,XU,DG,GGL,GGADS28130
     11,GG2,GGU,IDG,ALPL,ALPU)                                          ADS28140
      ALPHA=ALPL                                                        ADS28150
      IF (GMCTL.LT.0.) ALPHA=ALPU                                       ADS28160
      IF (ALPL.GT.ALPU) ALPHA=.5*(ALPL+ALPU)                            ADS28170
C     GO FIND ALPHA TO MINIMIZE THE OBJECTIVE.                          ADS28180
      GO TO 150                                                         ADS28190
140   CONTINUE                                                          ADS28200
C     NO CONSTRAINTS ARE VIOLATED OR ALL ARE VIOLATED.                  ADS28210
      IF (GMIN.LT.0.) GO TO 150                                         ADS28220
C     ALL ARE VIOLATED.  INTERPOLATE FOR MINIMUM VIOLATION.             ADS28230
      ICALC=2                                                           ADS28240
      CALL ADS103 (ICALC,XL,GMCTL,YP1,X1,GMCT1,X2,GMCT2,XU,GMCTU,ZRO,AZ,ADS28250
     1A1,A2,A3)                                                         ADS28260
      GO TO 160                                                         ADS28270
C     MINIMUM OF OBJECTIVE FUNCTION.                                    ADS28280
150   ICALC=2                                                           ADS28290
      CALL ADS103 (ICALC,XL,FL,YP1,X1,F1,X2,F2,XU,FU,ZRO,AZ,A1,A2,A3)   ADS28300
160   ICALC=3                                                           ADS28310
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,ALP,XMX1)                         ADS28320
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS28330
      IF (IPRNT1.GT.1) WRITE (6,420) XL,FL,X1,F1,X2,F2,XU,FU,AZ,A1,A2,A3ADS28340
     1,ALPHA                                                            ADS28350
C     REJECT IF ALPHA IS NOT UNIQUE.  GO PICK BEST AVAILABLE.           ADS28360
      DXL=ABS(XL-ALPHA)                                                 ADS28370
      DX1=ABS(ALPHA-X1)                                                 ADS28380
      DX2=ABS(ALPHA-X2)                                                 ADS28390
      DXU=ABS(ALPHA-XU)                                                 ADS28400
      DX=DMIN1(DXL,DX1,DX2,DXU)                                         ADS28410
      IEVAL=0                                                           ADS28420
      IF (DX.LE.ZRO) GO TO 290                                          ADS28430
C     TRY THIS DESIGN.                                                  ADS28440
C     UPDATE X.                                                         ADS28450
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS28460
C     EVALUATE OBJECTIVE FUNCTION AND CALL AGAIN.                       ADS28470
      IGOTO=3                                                           ADS28480
      TEMP(1)=XL                                                        ADS28490
      TEMP(2)=XU                                                        ADS28500
      TEMP(3)=X1                                                        ADS28510
      TEMP(4)=X2                                                        ADS28520
      TEMP(5)=SCFO                                                      ADS28530
      RETURN                                                            ADS28540
C     GO PICK BEST DESIGN UPON RETURN.                                  ADS28550
170   CONTINUE                                                          ADS28560
      IF (IONED.GT.7) GO TO 260                                         ADS28570
C     ------------------------------------------------------------------ADS28580
C     IONED = 7.  POLYNOMIAL INTERPOLATION WITH BOUNDS.                 ADS28590
C     ------------------------------------------------------------------ADS28600
      ALPU1=XU                                                          ADS28610
      JJ=0                                                              ADS28620
      II=ICAL                                                           ADS28630
      IF (II.GT.4) II=4                                                 ADS28640
      IF (II.GT.2) GO TO 180                                            ADS28650
C     II=2.  SET GG1 = GGU.                                             ADS28660
      CALL ADS503 (NCON,GGU,GG1)                                        ADS28670
      X1=XU                                                             ADS28680
      GMCT1=GMCTU                                                       ADS28690
      F1=FU                                                             ADS28700
      X2=XU                                                             ADS28710
      GO TO 220                                                         ADS28720
180   CONTINUE                                                          ADS28730
      IF (II.GT.3) GO TO 190                                            ADS28740
C     II=3. SET GG1=GG2 AND GG2=GGU.                                    ADS28750
      CALL ADS503 (NCON,GG2,GG1)                                        ADS28760
      CALL ADS503 (NCON,GGU,GG2)                                        ADS28770
      X1=X2                                                             ADS28780
      GMCT1=GMCT2                                                       ADS28790
      F1=F2                                                             ADS28800
      X2=XU                                                             ADS28810
      GMCT2=GMCTU                                                       ADS28820
      F2=FU                                                             ADS28830
      GO TO 220                                                         ADS28840
190   CONTINUE                                                          ADS28850
C     II=4.  SET G=GGU.                                                 ADS28860
      CALL ADS503 (NCON,GGU,G)                                          ADS28870
      GO TO 220                                                         ADS28880
200   CONTINUE                                                          ADS28890
C     IF FINISHED, GO FIND BEST DESIGN.                                 ADS28900
      IF (JJ.GE.3) GO TO 290                                            ADS28910
      IF (II.GT.3) GO TO 210                                            ADS28920
C     II=3, ALPHA=X2.                                                   ADS28930
      X2=ALPHA                                                          ADS28940
      F2=OBJ                                                            ADS28950
      GMCT2=GMCT                                                        ADS28960
      CALL ADS503 (NCON,G,GG2)                                          ADS28970
      GO TO 220                                                         ADS28980
210   CONTINUE                                                          ADS28990
C     II=4, ALPHA=XU.                                                   ADS29000
      XU=ALPHA                                                          ADS29010
      FU=OBJ                                                            ADS29020
      GMCTU=GMCT                                                        ADS29030
      CALL ADS503 (NCON,G,GGU)                                          ADS29040
220   CONTINUE                                                          ADS29050
      ALPL=XL                                                           ADS29060
      ALPU=ALPU1                                                        ADS29070
      GMIN=DMIN1(GMCTL,GMCT1,GMCT2,GMCTU)                               ADS29080
      ALPHA=.25*(XL+X1+X2+XU)                                           ADS29090
      IF (GMIN.LE.0.) GO TO 230                                         ADS29100
C     ALL DESIGNS ARE INFEASIBLE.  FIND MINIMUM CONSTRAINT VIOLATION.   ADS29110
      IF (II.LT.3) GO TO 250                                            ADS29120
      ICALC=2                                                           ADS29130
      IF (II.GT.3) ICALC=4                                              ADS29140
      CALL ADS108 (ICALC,XL,GMCTL,YP1,X1,GMCT1,X2,GMCT2,XU,GMCTU,ZRO,ALPADS29150
     1,YMN1)                                                            ADS29160
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS29170
      GO TO 250                                                         ADS29180
230   CONTINUE                                                          ADS29190
C     A FEASIBLE DESIGN EXISTS.  FIND NEW BOUNDS ON ALPHA.              ADS29200
      CALL ADS111 (NCON,II,ISLOPG,CTMIN,CTLMIN,ZRO,XL,X1,X2,XU,DG,GGL,GGADS29210
     11,GG2,G,IDG,ALPL,ALPU)                                            ADS29220
      IF (ABS(ALPU-ALPU1).GT.ZRO) ALPHA=ALPU                            ADS29230
C     FIND ALPHA TO MINIMIZE OBJECTIVE.                                 ADS29240
      ALP=-1.                                                           ADS29250
      IF (II.EQ.2.AND.ISLOPF.EQ.0) GO TO 240                            ADS29260
C     2-POINT QUADRATIC FIT.                                            ADS29270
      ICALC=1                                                           ADS29280
C     3-POINT QUADRATIC FIT.                                            ADS29290
      IF (II.EQ.3.AND.ISLOPF.EQ.0) ICALC=2                              ADS29300
C     3-POINT CUBIC FIT.                                                ADS29310
      IF (II.EQ.3.AND.ISLOPF.GT.0) ICALC=3                              ADS29320
C     4-POINT CUBIC FIT.                                                ADS29330
      IF (II.GT.3) ICALC=4                                              ADS29340
      CALL ADS108 (ICALC,XL,FL,SLOPF,X1,F1,X2,F2,XU,FU,ZRO,ALP,YMN1)    ADS29350
240   CONTINUE                                                          ADS29360
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS29370
      IF (FU.GT.FL.AND.GMCTL.GT.0.) ALPHA=ALPL                          ADS29380
      IF (ALPL.GT.ALPU) ALPHA=.5*(ALPL+ALPU)                            ADS29390
250   CONTINUE                                                          ADS29400
C     INSURE THAT ALPHA IS UNIQUE.  IF NOT, GO PICK BEST DESIGN.        ADS29410
      DXL=ABS(ALPHA-XL)                                                 ADS29420
      DX1=ABS(ALPHA-X1)                                                 ADS29430
      DX2=ABS(ALPHA-X2)                                                 ADS29440
      DXU=ABS(ALPHA-XU)                                                 ADS29450
      DX=DMIN1(DXL,DX1,DX2,DXU)                                         ADS29460
      IEVAL=0                                                           ADS29470
      IF (DX.LE.ZRO) GO TO 290                                          ADS29480
      JJ=II-1                                                           ADS29490
      IF (ISLOPF.GT.0) JJ=JJ+1                                          ADS29500
      II=II+1                                                           ADS29510
C     UPDATE X.                                                         ADS29520
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS29530
C     EVALUATE OBJECTIVE AND CONSTRAINT FUNCTIONS AND CALL AGAIN.       ADS29540
      IGOTO=4                                                           ADS29550
      TEMP(1)=XL                                                        ADS29560
      TEMP(2)=XU                                                        ADS29570
      TEMP(3)=X1                                                        ADS29580
      TEMP(4)=X2                                                        ADS29590
      TEMP(5)=SCFO                                                      ADS29600
      RETURN                                                            ADS29610
260   CONTINUE                                                          ADS29620
C     ------------------------------------------------------------------ADS29630
C     IONED = 9.  POLYNOMIAL SEARCH WITHOUT FIRST GETTING BOUNDS.       ADS29640
C     ------------------------------------------------------------------ADS29650
      IF (IPRNT1.GT.0) WRITE (6,340)                                    ADS29660
      JGOTO=0                                                           ADS29670
270   CONTINUE                                                          ADS29680
      CALL ADS119 (JGOTO,NCON,XMULT,EXTRAP,ALPMAX,SLOPF,ISLOPF,ISLOPG,ZRADS29690
     1O,CTMIN,CTLMIN,GMCT,G,DG,GGL,GG1,GG2,GGU,IDG,ALPHA,OBJ,XL,X1,X2,XUADS29700
     2,GMCTL,GMCT1,GMCT2,GMCTU,FL,F1,F2,FU)                             ADS29710
C     INSURE THAT ALPHA IS UNIQUE.  IF NOT, GO PICK BEST DESIGN.        ADS29720
      DXL=ABS(ALPHA-XL)                                                 ADS29730
      DX1=ABS(ALPHA-X1)                                                 ADS29740
      DX2=ABS(ALPHA-X2)                                                 ADS29750
      DXU=ABS(ALPHA-XU)                                                 ADS29760
      DX=DMIN1(DXL,DX1,DX2,DXU)                                         ADS29770
      IF (DX.LE.ZRO) GO TO 290                                          ADS29780
      IF (JGOTO.EQ.0) GO TO 280                                         ADS29790
C     UPDATE X.                                                         ADS29800
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS29810
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS29820
      IGOTO=5                                                           ADS29830
      TEMP(1)=XL                                                        ADS29840
      TEMP(2)=XU                                                        ADS29850
      TEMP(3)=X1                                                        ADS29860
      TEMP(4)=X2                                                        ADS29870
      TEMP(5)=SCFO                                                      ADS29880
      RETURN                                                            ADS29890
280   CONTINUE                                                          ADS29900
290   CONTINUE                                                          ADS29910
C     ------------------------------------------------------------------ADS29920
C     PICK BEST OF ALL DESIGNS CONSIDERED.  FIND ALPL AND ALPU.         ADS29930
C     ------------------------------------------------------------------ADS29940
C     RETRIEVE XBST.                                                    ADS29950
      CALL ADS503 (NDV,XBST,X)                                          ADS29960
      OBJ=FBST                                                          ADS29970
      GMCT=GBST                                                         ADS29980
      DXZ=ABS(ALPHA-ALPBST)                                             ADS29990
      IF (DXZ.LE.ZRO.AND.IEVAL.EQ.1) GO TO 330                          ADS30000
      DXL=ABS(ALPBST-XL)                                                ADS30010
      DX1=ABS(ALPBST-X1)                                                ADS30020
      DX2=ABS(ALPBST-X2)                                                ADS30030
      DXU=ABS(ALPBST-XU)                                                ADS30040
C     COMPARE WITH XL.                                                  ADS30050
      IF (DXL.GT.DXZ) GO TO 300                                         ADS30060
C     XL IS BETTER.                                                     ADS30070
      DXZ=DXL                                                           ADS30080
      CALL ADS503 (NCON,GGL,G)                                          ADS30090
C     COMPARE WITH X1.                                                  ADS30100
300   IF (DX1.GT.DXZ) GO TO 310                                         ADS30110
C     X1 IS BETTER.                                                     ADS30120
      DXZ=DX1                                                           ADS30130
      CALL ADS503 (NCON,GG1,G)                                          ADS30140
C     COMPARE WITH X2.                                                  ADS30150
310   IF (DX2.GT.DXZ) GO TO 320                                         ADS30160
C     X2 IS BETTER.                                                     ADS30170
      DXZ=DX2                                                           ADS30180
      CALL ADS503 (NCON,GG2,G)                                          ADS30190
C     COMPARE WITH XU.                                                  ADS30200
320   IF (DXU.GT.DXZ) GO TO 330                                         ADS30210
C     XU IS BETTER.                                                     ADS30220
      CALL ADS503 (NCON,GGU,G)                                          ADS30230
330   CONTINUE                                                          ADS30240
      ALPHA=ALPBST                                                      ADS30250
C     FIND ALPL AND ALPU.                                               ADS30260
      ALPL=DMIN1(XL,X1,X2,XU)                                           ADS30270
      ALPU=DMAX1(XL,X1,X2,XU)                                           ADS30280
      IF (XL.GT.ALPL.AND.XL.LT.ALPHA) ALPL=XL                           ADS30290
      IF (X1.GT.ALPL.AND.X1.LT.ALPHA) ALPL=X1                           ADS30300
      IF (X2.GT.ALPL.AND.X2.LT.ALPHA) ALPL=X2                           ADS30310
      IF (XU.GT.ALPL.AND.XU.LT.ALPHA) ALPL=XU                           ADS30320
      IF (XL.LT.ALPU.AND.XL.GT.ALPHA) ALPU=XL                           ADS30330
      IF (X1.LT.ALPU.AND.X1.GT.ALPHA) ALPU=X1                           ADS30340
      IF (X2.LT.ALPU.AND.X2.GT.ALPHA) ALPU=X2                           ADS30350
      IF (XU.LT.ALPU.AND.XU.GT.ALPHA) ALPU=XU                           ADS30360
      IGOTO=0                                                           ADS30370
      TEMP(1)=XL                                                        ADS30380
      TEMP(2)=XU                                                        ADS30390
      TEMP(3)=X1                                                        ADS30400
      TEMP(4)=X2                                                        ADS30410
      TEMP(5)=SCFO                                                      ADS30420
      IF (IPRNT1.LE.0) RETURN                                           ADS30430
      OBJW=OBJ                                                          ADS30440
      IF(ISTRAT.EQ.0) OBJW=OBJ/SCFO                                     ADS30450
      WRITE (6,350) ALPL,ALPU,ALPHA,OBJW                                ADS30460
      RETURN                                                            ADS30470
C     ------------------------------------------------------------------ADS30480
C     FORMATS.                                                          ADS30490
C     ------------------------------------------------------------------ADS30500
340   FORMAT (/,5X,'BEGIN POLYNOMIAL INTERPOLATION WITHOUT FIRST ',     ADS30510
     1'GETTING BOUNDS')                                                 ADS30520
350   FORMAT (/,5X,'RESULTS OF ONE-DIMENSIONAL SEARCH',/,6X,'ALPL =',   ADS30530
     11PE13.5,2X,'ALPU =',1PE13.5,/,5X,'ALPHA =',1PE13.5,3X,'OBJ =',    ADS30540
     21PE13.5)                                                          ADS30550
360   FORMAT (/,5X,'CONTROL PARAMETERS',/,5X,'SLOPF, ALPHA, ALPMAX,',   ADS30560
     1' XMULT,  EXTRAP',/,5X,5(1PE13.5),/,5X,'ZRO, DELALP, DABALP',/,5X,ADS30570
     23(1PE13.5),/,5X,'ISLOPF, ISLOPG',2I5)                             ADS30580
370   FORMAT (/,5X,'FIND BRACKETS ON MINIMUM')                          ADS30590
380   FORMAT (/,5X,'BRACKETS, XL =',1PE13.5,2X,'XU =',1PE13.5)          ADS30600
390   FORMAT (/,5X,'BEGIN GOLDEN SECTION SEARCH')                       ADS30610
400   FORMAT (/,5X,'GOLDEN SECTION SEARCH IS COMPLETE')                 ADS30620
410   FORMAT (/,5X,'BEGIN CUBIC INTERPOLATION ON RESULTS OF GOLDEN',    ADS30630
     1' SECTION SEARCH')                                                ADS30640
420   FORMAT (/,8X,'ALPHA',11X,'F',4(/,5X,2(1PE13.5)),/,12X,'A0',11X,   ADS30650
     1'A1',11X,'A2',11X,'A3',/,5X,4(1PE13.5),/,5X,'ESTIMATED ALPHA =',  ADS30660
     21PE13.5,2X,'(NOT USED IF LT.XL OR GT.XU)')                        ADS30670
      END                                                               ADS30680
C --- SUBROUTINE ADS116(IGOTO,IONED,IPRNT1,NDV,NCON,OBJ,SLOPE,ISLOPE,   ADS30690
C    * ALPHA,ALPMAX,XMULT,EXTRAP,ZRO,DELALP,DABALP,ALPL,ALPU,           ADS30700
C    * XOLD,SVLB,SVUB,S,X,G,GGL,XSAV)                                   ADS30710
C     DIMENSION X(NDV),S(NDV),XOLD(NDV),SVLB(NDV),SVUB(NDV),G(NCON),    ADS30720
C    * GGL(NCON),XSAV(NDV)                                              ADS30730
C     ------------------------------------------------------------------ADS30740
C --- ROUTINE TO CONTROL ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED       ADS30750
C     FUNCTIONS.                                                        ADS30760
C     ------------------------------------------------------------------ADS30770
C     BY G. N. VANDERPLAATS                                 ***** 1983  ADS30780
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS30790
      SUBROUTINE ADS116 (IGOTO,IONED,IPRNT1,NDV,NCON,OBJ,SLOPE,ISLOPE,ALADS30800
     1PHA,ALPMAX,XMULT,EXTRAP,ZRO,DELALP,DABALP,ALPL,ALPU,XOLD,SVLB,SVUBADS30810
     2,S,X,G,GGL,XSAV,FBST,FL,FU,F1,F2,ALPBST,XL,XU,X1,X2,ICAL,II,JGOTO,ADS30820
     3II1,IX,NITER,SCFX,SCFG,ISTRAT,SCFO)                               ADS30830
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),S(1),XOLD(1),SVLB(1),SVUB(1),G(1),GGL(1),XSAV(1),  ADS30840
     1SCFX(1),SCFG(1)                                                   ADS30850
      ISLOPE=0                                                          ADS30860
      IF (IGOTO.GT.0) GO TO 10                                          ADS30870
C     IGOTO = 0, FIRST CALL.                                            ADS30880
      ICAL=1                                                            ADS30890
      ALPL=0.                                                           ADS30900
      ALPU=ALPMAX                                                       ADS30910
      ALPBST=0.                                                         ADS30920
      FBST=OBJ                                                          ADS30930
      CALL ADS503 (NDV,X,XSAV)                                          ADS30940
      IF (NCON.GT.0) CALL ADS503 (NCON,G,GGL)                           ADS30950
C     SAVE X-VECTOR IN XOLD.                                            ADS30960
      CALL ADS503 (NDV,X,XOLD)                                          ADS30970
      IF (IPRNT1.GT.1) WRITE (6,140) SLOPE,ALPHA,ALPMAX,XMULT,EXTRAP,ZROADS30980
     1,DELALP,DABALP                                                    ADS30990
      GO TO 20                                                          ADS31000
10    CONTINUE                                                          ADS31010
C     ------------------------------------------------------------------ADS31020
C     IDENTIFY BEST DESIGN.                                             ADS31030
C     ------------------------------------------------------------------ADS31040
      ICAL=ICAL+1                                                       ADS31050
      IF (OBJ.GT.FBST) GO TO 20                                         ADS31060
      IF (ALPBST.GT.ALPL.AND.ALPBST.LT.ALPHA) ALPL=ALPBST               ADS31070
      IF (ALPBST.LT.ALPU.AND.ALPBST.GT.ALPHA) ALPU=ALPBST               ADS31080
      ALPBST=ALPHA                                                      ADS31090
      FBST=OBJ                                                          ADS31100
      CALL ADS503 (NDV,X,XSAV)                                          ADS31110
      IF (NCON.GT.0) CALL ADS503 (NCON,G,GGL)                           ADS31120
20    IGOTO1=IGOTO+1                                                    ADS31130
      GO TO (30,40,60,120,90,110), IGOTO1                               ADS31140
C     IGOTO = 0.  FIRST CALL.                                           ADS31150
30    CONTINUE                                                          ADS31160
C     DO WE NEED TO BRACKET THE MINIMUM?                                ADS31170
      IF (IONED.GT.3) GO TO 100                                         ADS31180
C     YES.                                                              ADS31190
C     ------------------------------------------------------------------ADS31200
C     FIND BRACKETS ON MINIMUM.                                         ADS31210
C     ------------------------------------------------------------------ADS31220
      JGOTO=0                                                           ADS31230
      IF (IPRNT1.GT.1) WRITE (6,160)                                    ADS31240
40    CONTINUE                                                          ADS31250
      CALL ADS121 (JGOTO,OBJ,ALPHA,XMULT,ALPMAX,XL,X1,X2,XU,FL,F1,F2,FU)ADS31260
      IF (JGOTO.EQ.0) GO TO 50                                          ADS31270
C     UPDATE X.                                                         ADS31280
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS31290
C     EVALUATE OBJECTIVE FUNCTION AND CALL AGAIN.                       ADS31300
      IGOTO=1                                                           ADS31310
      RETURN                                                            ADS31320
50    CONTINUE                                                          ADS31330
      IF (IPRNT1.GT.1) WRITE (6,170) XL,XU                              ADS31340
C     ------------------------------------------------------------------ADS31350
C     MINIMUM IS BRACKETED.  REFINE SOLUTION AS REQUIRED.               ADS31360
C     ------------------------------------------------------------------ADS31370
      IF (IONED.GT.2) GO TO 80                                          ADS31380
C     ------------------------------------------------------------------ADS31390
C     IONED = 1 OR 2.  REFINE SOLUTION USING THE GOLDEN SECTION METHOD. ADS31400
C     ------------------------------------------------------------------ADS31410
      JGOTO=0                                                           ADS31420
      IF (IPRNT1.GT.1) WRITE (6,180)                                    ADS31430
      NPTS=ICAL                                                         ADS31440
      XX=XMULT-2.618034                                                 ADS31450
      IF (ABS(XX).GT.0.0001) NPTS=2                                     ADS31460
      IF (NPTS.NE.3) GO TO 60                                           ADS31470
      X1=X2                                                             ADS31480
      F1=F2                                                             ADS31490
60    CONTINUE                                                          ADS31500
      CALL ADS126 (JGOTO,NPTS,OBJ,ALPHA,ALPMAX,DELALP,DABALP,ZRO,XL,X1,XADS31510
     12,XU,FL,F1,F2,FU,II1,IX,NITER)                                    ADS31520
C     FINISHED?                                                         ADS31530
      IF (JGOTO.EQ.0) GO TO 70                                          ADS31540
C     NO. UPDATE X.                                                     ADS31550
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS31560
C     EVALUATE OBJECTIVE FUNCTION AND CALL AGAIN.                       ADS31570
      IGOTO=2                                                           ADS31580
      RETURN                                                            ADS31590
70    CONTINUE                                                          ADS31600
      IF (IPRNT1.GT.1) WRITE (6,190)                                    ADS31610
C     IF IONED = 1, SOLUTION IS COMPLETE.                               ADS31620
C     IF IONED = 2, DO POLYNOMIAL INTERPOLATION.                        ADS31630
      IGOTO=0                                                           ADS31640
      IF (IONED.EQ.1) GO TO 120                                         ADS31650
C     ------------------------------------------------------------------ADS31660
C     IONED = 2.  POLYNOMIAL INTERPOLATION USING LAST FOUR GOLDEN       ADS31670
C     SECTION POINTS.                                                   ADS31680
C     ------------------------------------------------------------------ADS31690
      IF (IPRNT1.GT.1) WRITE (6,200)                                    ADS31700
      ICALC=2                                                           ADS31710
      CALL ADS103 (ICALC,XL,FL,YP1,X1,F1,X2,F2,XU,FU,ZRO,AZ,A1,A2,A3)   ADS31720
      ICALC=3                                                           ADS31730
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,ALPHA,XMX1)                       ADS31740
      IF (IPRNT1.GT.1) WRITE (6,210) XL,FL,X1,F1,X2,F2,XU,FU,AZ,A1,A2,A3ADS31750
     1,ALPHA                                                            ADS31760
C     REJECT IF OUT OF BOUNDS.  GO PICK BEST AVAILABLE.                 ADS31770
      IGOTO=0                                                           ADS31780
      IF (ALPHA.LT.XL.OR.ALPHA.GT.XU) GO TO 120                         ADS31790
C     TRY THIS DESIGN.                                                  ADS31800
C     UPDATE X.                                                         ADS31810
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS31820
C     EVALUATE OBJECTIVE FUNCTION AND CALL AGAIN.                       ADS31830
      IGOTO=3                                                           ADS31840
      RETURN                                                            ADS31850
C     GO PICK BEST DESIGN UPON RETURN.                                  ADS31860
80    CONTINUE                                                          ADS31870
      IF (IONED.GT.3) GO TO 100                                         ADS31880
C     ------------------------------------------------------------------ADS31890
C     IONED = 3.  POLYNOMIAL INTERPOLATION WITH BOUNDS.                 ADS31900
C     ------------------------------------------------------------------ADS31910
      II=ICAL                                                           ADS31920
      IF (II.GT.4) II=4                                                 ADS31930
C     II = NUMBER OF POINTS FOR INTERPOLATION.                          ADS31940
C         2 - XL,X1                                                     ADS31950
C         3 - XL,X1,X2                                                  ADS31960
C         4 - XL,X1,X2,XU                                               ADS31970
      JGOTO=0                                                           ADS31980
      IF (II.NE.3) GO TO 90                                             ADS31990
C     II=3                                                              ADS32000
      X1=X2                                                             ADS32010
      F1=F2                                                             ADS32020
90    CONTINUE                                                          ADS32030
      CALL ADS117 (JGOTO,II,XL,X1,X2,XU,FL,F1,F2,FU,SLOPE,ISLOPE,ALPHA,OADS32040
     1BJ,ZRO)                                                           ADS32050
      IF (JGOTO.EQ.0) GO TO 120                                         ADS32060
C     UPDATE X.                                                         ADS32070
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS32080
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS32090
      IGOTO=4                                                           ADS32100
      RETURN                                                            ADS32110
100   CONTINUE                                                          ADS32120
C     ------------------------------------------------------------------ADS32130
C     IONED = 4.  POLYNOMIAL SEARCH WITHOUT FIRST GETTING BOUNDS.       ADS32140
C     ------------------------------------------------------------------ADS32150
      IF (IPRNT1.GT.1) WRITE (6,130)                                    ADS32160
      JGOTO=0                                                           ADS32170
110   CONTINUE                                                          ADS32180
      CALL ADS118 (JGOTO,ZRO,XMULT,EXTRAP,ALPMAX,SLOPE,ISLOPE,ALPHA,OBJ,ADS32190
     1FL,FU,F1,F2,XL,XU,X1,X2)                                          ADS32200
      IF (JGOTO.EQ.0) GO TO 120                                         ADS32210
C     UPDATE X.                                                         ADS32220
      CALL ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                        ADS32230
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS32240
      IGOTO=5                                                           ADS32250
      RETURN                                                            ADS32260
120   CONTINUE                                                          ADS32270
C     ------------------------------------------------------------------ADS32280
C     SEARCH COMPLETE.                                                  ADS32290
C     ------------------------------------------------------------------ADS32300
C     WRITE LAST OBJ IF PRINT IS TURNED ON.                             ADS32310
      OBJW=OBJ                                                          ADS32320
      IF(ISTRAT.EQ.0) OBJW=OBJ/SCFO                                     ADS32330
      IF (IGOTO.EQ.3.AND.IPRNT1.GT.1) WRITE (6,150) OBJW                ADS32340
      IF (IGOTO.EQ.6.AND.IPRNT1.GT.1) WRITE (6,150) OBJW                ADS32350
      IF (IGOTO.EQ.10.AND.IPRNT1.GT.1) WRITE (6,150) OBJW               ADS32360
      IF (IPRNT1.GT.0) WRITE (6,220) ICAL                               ADS32370
      ALPHA=ALPBST                                                      ADS32380
      OBJ=FBST                                                          ADS32390
C     RETRIEVE ALPBST.                                                  ADS32400
      CALL ADS503 (NDV,XSAV,X)                                          ADS32410
C     RETRIEVE CONSTRAINT VALUES.                                       ADS32420
      IF (NCON.GT.0) CALL ADS503 (NCON,GGL,G)                           ADS32430
      IGOTO=0                                                           ADS32440
      RETURN                                                            ADS32450
C     ------------------------------------------------------------------ADS32460
C     FORMATS.                                                          ADS32470
C     ------------------------------------------------------------------ADS32480
130   FORMAT (/,5X,'BEGIN POLYNOMIAL INTERPOLATION WITHOUT FIRST',      ADS32490
     1' GETTING BOUNDS')                                                ADS32500
140   FORMAT (/,5X,'CONTROL PARAMETERS',/,5X,'SLOPE, ALPHA, ALPMAX,',   ADS32510
     1' XMULT,  EXTRAP',/,5X,5(1PE13.5),/,5X,'ZRO, DELALP, DABALP',//5X,ADS32520
     23(1PE13.5))                                                       ADS32530
150   FORMAT (/,5X,'OBJ =',1PE13.5)                                     ADS32540
160   FORMAT (/,5X,'FIND BRACKETS ON MINIMUM')                          ADS32550
170   FORMAT (/,5X,'BRACKETS, XL =',1PE13.5,2X,'XU =',1PE13.5)          ADS32560
180   FORMAT (/,5X,'BEGIN GOLDEN SECTION SEARCH')                       ADS32570
190   FORMAT (/,5X,'GOLDEN SECTION SEARCH IS COMPLETE')                 ADS32580
200   FORMAT (/,5X,'BEGIN CUBIC INTERPOLATION ON RESULTS OF GOLDEN',    ADS32590
     1' SECTION SEARCH')                                                ADS32600
210   FORMAT (/,5X,'ALPHA',11X,'F',4(/,5X,2(1PE13.5)),/,12X,'A0',11X,   ADS32610
     1'A1',11X,'A2',11X,'A3',/,5X,4(1PE13.5),/,5X,'ESTIMATED ALPHA =',  ADS32620
     21PE13.5,2X,'(NOT USED IF LT.XL OR GT.XU)')                        ADS32630
220   FORMAT (/,5X,'FUNCTION CALLS =',I5)                               ADS32640
      END                                                               ADS32650
C --- SUBROUTINE ADS117 (IGOTO,II,XL,X1,X2,XU,FL,F1,F2,FU,SLOPE,ISLOPE, ADS32660
C    * ALPHA,OBJ,ZRO)                                                   ADS32670
C     ------------------------------------------------------------------ADS32680
C --- POLYNOMIAL INTERPOLATION FOR MINIMUM OF UNCONSTRAINED FUNCTION,   ADS32690
C     GIVEN INITIAL BOUNDS.                                             ADS32700
C     ------------------------------------------------------------------ADS32710
C     BY G. N. VANDERPLAATS                            JUNE 1982.       ADS32720
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS32730
C                                                                       ADS32740
C --- INPUT                                                             ADS32750
C                                                                       ADS32760
C --- OUTPUT                                                            ADS32770
      SUBROUTINE ADS117 (IGOTO,II,XL,X1,X2,XU,FL,F1,F2,FU,SLOPE,ISLOPE,AADS32780
     1LPHA,OBJ,ZRO)                                                     ADS32790
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS32800
      IGOTO1=IGOTO+1                                                    ADS32810
      GO TO (10,30,80,100), IGOTO1                                      ADS32820
10    CONTINUE                                                          ADS32830
      IF (II.GT.2) GO TO 40                                             ADS32840
C --- II=2.  DO QUADRATIC INTERPOLATION IF ISLOPE.GT.0.  ELSE           ADS32850
C     PICK DEFAULT ESTIMATE FOR X1.                                     ADS32860
      FA=FL+2.*ABS(FL)                                                  ADS32870
      X1=.5*XU                                                          ADS32880
      IF (FU.GT.FA) X1=.1*XU                                            ADS32890
      IF (ISLOPE.EQ.0) GO TO 20                                         ADS32900
C     QUADRATIC INTERPOLATION.                                          ADS32910
      ICALC=1                                                           ADS32920
      CALL ADS102 (ICALC,XL,FL,SLOPE,XU,FU,XX,YY,AZ,A1,A2)              ADS32930
      ICALC=2                                                           ADS32940
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,XMN1,XMX1)                        ADS32950
      IF (XMN1.GT.XL.AND.XMN1.LT.XU) X1=XMN1                            ADS32960
20    CONTINUE                                                          ADS32970
      ALPHA=X1                                                          ADS32980
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS32990
      IGOTO=1                                                           ADS33000
      RETURN                                                            ADS33010
30    CONTINUE                                                          ADS33020
      F1=OBJ                                                            ADS33030
40    CONTINUE                                                          ADS33040
      IF (II.GT.3) GO TO 90                                             ADS33050
C --- II=3.  DO 3-POINT INTERPOLATION.                                  ADS33060
      IF (ISLOPE.GT.0) GO TO 50                                         ADS33070
C     3-POINT QUADRATIC COEFFICIENTS.                                   ADS33080
      ICALC=2                                                           ADS33090
      CALL ADS102 (ICALC,XL,FL,YP1,X1,F1,XU,FU,AZ,A1,A2)                ADS33100
      GO TO 60                                                          ADS33110
50    CONTINUE                                                          ADS33120
C     3-POINT CUBIC INTERPOLATION.                                      ADS33130
      ICALC=1                                                           ADS33140
      CALL ADS103 (ICALC,XL,FL,SLOPE,X1,F1,XU,FU,XX,YY,ZRO,AZ,A1,A2,A3) ADS33150
      ICALC=3                                                           ADS33160
60    CONTINUE                                                          ADS33170
C     FIND MINIMUM.                                                     ADS33180
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,XMN1,XMX1)                        ADS33190
      X2=.5*(XU+X1)                                                     ADS33200
      IF (F1.LT.FL) GO TO 70                                            ADS33210
C     X1 IS UPPER BOUND.  PICK BETTER X2.                               ADS33220
      FA=FL+2.*ABS(FL)                                                  ADS33230
      X2=.5*(XL+X1)                                                     ADS33240
      IF (F1.GT.FA) X2=XL+.1*(X1-XL)                                    ADS33250
70    CONTINUE                                                          ADS33260
      IF (XMN1.GT.XL.AND.XMN1.LT.XU) X2=XMN1                            ADS33270
      ALPHA=X2                                                          ADS33280
C     CHECK TO BE SURE THIS ALPHA IS UNIQUE.                            ADS33290
      DXL=ABS(X2-XL)                                                    ADS33300
      DX1=ABS(X2-X1)                                                    ADS33310
      DXU=ABS(X2-XU)                                                    ADS33320
      DX=DMIN1(DXL,DX1,DXU)                                             ADS33330
      IF (DX.LE.ZRO) GO TO 100                                          ADS33340
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS33350
      IGOTO=2                                                           ADS33360
      RETURN                                                            ADS33370
80    CONTINUE                                                          ADS33380
      F2=OBJ                                                            ADS33390
      IF (ISLOPE.GT.0) GO TO 100                                        ADS33400
90    CONTINUE                                                          ADS33410
C --- II=4.  DO 4-POINT CUBIC INTERPOLATION.                            ADS33420
      ICALC=2                                                           ADS33430
      CALL ADS103 (ICALC,XL,FL,SLOPE,X1,F1,X2,F2,XU,FU,ZRO,AZ,A1,A2,A3) ADS33440
      ICALC=3                                                           ADS33450
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,ALPHA,XMX1)                       ADS33460
C     IF ALPHA IS OUT OF BOUNDS, EXIT.                                  ADS33470
      IF (ALPHA.LE.XL.OR.ALPHA.GE.XU) GO TO 100                         ADS33480
C     IF THIS ALPHA IS NOT UNIQUE, EXIT.                                ADS33490
      DXL=ABS(ALPHA-XL)                                                 ADS33500
      DX1=ABS(ALPHA-X1)                                                 ADS33510
      DX2=ABS(ALPHA-X2)                                                 ADS33520
      DXU=ABS(ALPHA-XU)                                                 ADS33530
      DX=DMIN1(DXL,DX1,DX2,DXU)                                         ADS33540
      IF (DX.LE.ZRO) GO TO 100                                          ADS33550
      IGOTO=3                                                           ADS33560
      RETURN                                                            ADS33570
100   CONTINUE                                                          ADS33580
C     EXIT.                                                             ADS33590
      IGOTO=0                                                           ADS33600
      RETURN                                                            ADS33610
      END                                                               ADS33620
C --- SUBROUTINE ADS118 (IGOTO,ZRO,XMULT,EXTRAP,ALPMAX,SLOPE,           ADS33630
C    * ISLOPE,ALPHA,OBJ)                                                ADS33640
C     ------------------------------------------------------------------ADS33650
C --- POLYNOMIAL INTERPOLATION/EXTRAPOLATION FOR UNCONSTRAINED MINIMUM  ADS33660
C     WITHOUT FIRST GETTING BOUNDS.                                     ADS33670
C     ------------------------------------------------------------------ADS33680
C     BY G. N. VANDERPLAATS                             JUNE 1982.      ADS33690
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS33700
C                                                                       ADS33710
C --- INPUT                                                             ADS33720
C                                                                       ADS33730
C --- OUTPUT                                                            ADS33740
      SUBROUTINE ADS118 (IGOTO,ZRO,XMULT,EXTRAP,ALPMAX,SLOPE,ISLOPE,ALPHADS33750
     1A,OBJ,FL,FU,F1,F2,XL,XU,X1,X2)                                    ADS33760
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS33770
      IGOTO1=IGOTO+1                                                    ADS33780
      GO TO (10,20,40,70,80), IGOTO1                                    ADS33790
10    CONTINUE                                                          ADS33800
C     FIRST CALL.  EVALUATE FUNCTION AT X1=ALPHA AND CALL AGAIN.        ADS33810
      XL=0.                                                             ADS33820
      FL=OBJ                                                            ADS33830
      X1=ALPHA                                                          ADS33840
      IGOTO=1                                                           ADS33850
      RETURN                                                            ADS33860
20    CONTINUE                                                          ADS33870
      F1=OBJ                                                            ADS33880
C     ------------------------------------------------------------------ADS33890
C     TWO-POINT INTERPOLATION.                                          ADS33900
C     ------------------------------------------------------------------ADS33910
C                                                                       ADS33920
C     DO 2-POINT QUADRATIC INTERPOLATION IF ISLOPE.GT.0.  ELSE USE      ADS33930
C     XMULT*X1 OR 0.5*X1 TO GET NEW ALPHA.                              ADS33940
      X2=X1*XMULT                                                       ADS33950
      IF (F1.GE.FL) X2=0.5*X1                                           ADS33960
      IF (ISLOPE.EQ.0) GO TO 30                                         ADS33970
C     2-POINT QUADRATIC COEFFICIENTS.                                   ADS33980
      ICALC=1                                                           ADS33990
      CALL ADS102 (ICALC,XL,FL,SLOPE,X1,F1,XX,YY,AZ,A1,A2)              ADS34000
      ICALC=2                                                           ADS34010
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,XMN1,XMX1)                        ADS34020
      X2=XMN1                                                           ADS34030
      IF (XMN1.LE.0.0) X2=.5*X1                                         ADS34040
      XX=EXTRAP*X1                                                      ADS34050
      IF (XMN1.GT.XX) X2=XX                                             ADS34060
30    CONTINUE                                                          ADS34070
      IF (X2.GT.ALPMAX) X2=ALPMAX                                       ADS34080
C     CHECK TO BE SURE THIS ALPHA IS UNIQUE.  IF NOT, EXIT.             ADS34090
      DX=ABS(X2-X1)                                                     ADS34100
      IF (DX.LE.ZRO) GO TO 80                                           ADS34110
      ALPHA=X2                                                          ADS34120
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS34130
      IGOTO=2                                                           ADS34140
      RETURN                                                            ADS34150
40    CONTINUE                                                          ADS34160
      F2=OBJ                                                            ADS34170
      IF (X2.GE.(0.99999*ALPMAX)) GO TO 80                              ADS34180
C     ------------------------------------------------------------------ADS34190
C     3-POINT INTERPOLATION.                                            ADS34200
C     ------------------------------------------------------------------ADS34210
      IF (ISLOPE.GT.0) GO TO 50                                         ADS34220
C     3-POINT QUADRATIC FIT.                                            ADS34230
      ICALC=2                                                           ADS34240
      CALL ADS102 (ICALC,XL,FL,SLOPE,X1,F1,X2,F2,AZ,A1,A2)              ADS34250
      GO TO 60                                                          ADS34260
50    CONTINUE                                                          ADS34270
C     3-POINT CUBIC COEFFICIENTS.                                       ADS34280
      ICALC=1                                                           ADS34290
      CALL ADS103 (ICALC,XL,FL,SLOPE,X1,F1,X2,F2,XX,YY,ZRO,AZ,A1,A2,A3) ADS34300
      ICALC=3                                                           ADS34310
60    CONTINUE                                                          ADS34320
C     FIND MINIMUM.                                                     ADS34330
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,XU,XMX1)                          ADS34340
C     IF ISLOPE.GT.0 AND XU.LT.XL, EXIT.  ELSE PICK NEW POINT.          ADS34350
      IF (ISLOPE.GT.0.AND.XU.LE.XL) GO TO 80                            ADS34360
      XX=X2*EXTRAP                                                      ADS34370
      IF (XU.GT.XX) XU=XX                                               ADS34380
      IF (XU.GT.ALPMAX) XU=ALPMAX                                       ADS34390
      IF (XU.LE.0.) XU=.5*DMIN1(X1,X2)                                  ADS34400
C     CHECK TO BE SURE THIS ALPHA IS UNIQUE.  IF NOT, EXIT.             ADS34410
      DX1=ABS(XU-X1)                                                    ADS34420
      DX2=ABS(XU-X2)                                                    ADS34430
      DX=DMIN1(DX1,DX2)                                                 ADS34440
      IF (DX.LE.ZRO) GO TO 80                                           ADS34450
      ALPHA=XU                                                          ADS34460
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS34470
      IGOTO=3                                                           ADS34480
      RETURN                                                            ADS34490
70    CONTINUE                                                          ADS34500
      FU=OBJ                                                            ADS34510
      IF (XU.GE.(0.99999*ALPMAX)) GO TO 80                              ADS34520
      IF (ISLOPE.GT.0) GO TO 80                                         ADS34530
C     ------------------------------------------------------------------ADS34540
C     4-POINT CUBIC FIT.                                                ADS34550
C     ------------------------------------------------------------------ADS34560
      ICALC=2                                                           ADS34570
      CALL ADS103 (ICALC,XL,FL,SLOPE,X1,F1,X2,F2,XU,FU,ZRO,AZ,A1,A2,A3) ADS34580
      ICALC=3                                                           ADS34590
      CALL ADS105 (ICALC,A1,A2,A3,ZRO,ALPHA,XMX1)                       ADS34600
C     IF ALPHA.LT.0, EXIT.                                              ADS34610
      IF (ALPHA.LT.0.) GO TO 80                                         ADS34620
C     LIMIT ALPHA TO XU*EXTRAP                                          ADS34630
      XX=XU*EXTRAP                                                      ADS34640
      IF (ALPHA.GT.XX) ALPHA=XX                                         ADS34650
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS34660
C     CHECK TO BE SURE THIS ALPHA IS UNIQUE.  IF NOT, EXIT.             ADS34670
      DX1=ABS(ALPHA-X1)                                                 ADS34680
      DX2=ABS(ALPHA-X2)                                                 ADS34690
      DXU=ABS(ALPHA-XU)                                                 ADS34700
      DX=DMIN1(DX1,DX2,DXU)                                             ADS34710
      IF (DX.LE.ZRO) GO TO 80                                           ADS34720
C     EVALUATE FUNCTION AND CALL AGAIN.                                 ADS34730
      IGOTO=4                                                           ADS34740
      RETURN                                                            ADS34750
80    CONTINUE                                                          ADS34760
C     EXIT.                                                             ADS34770
      IGOTO=0                                                           ADS34780
      RETURN                                                            ADS34790
      END                                                               ADS34800
C --- SUBROUTINE ADS119 (IGOTO,NCON,XMULT,EXTRAP,ALPMAX,SLOPE,ISLOPE,   ADS34810
C    1ISLOPG,ZRO,CTMIN,CTLMIN,GMCT,G,DG,GGL,GG1,GG2,GGU,                ADS34820
C    2IDG,ALPHA,OBJ,XL,X1,X2,XU,GMCTL,GMCT1,GMCT2,GMCTU,FL,F1,F2,FU)    ADS34830
C     ------------------------------------------------------------------ADS34840
C --- ROUTINE TO INTERPOLATE/EXTRAPOLATE FOR CONSTRAINED MINIMUM        ADS34850
C     WITHOUT FIRST GETTING BOUNDS.                                     ADS34860
C     ------------------------------------------------------------------ADS34870
C     BY G. N. VANDERPLAATS                               JUNE 1982.    ADS34880
C     NAVAL POSTGRADUATE SCHOOL, MONETERY, CALIF.                       ADS34890
C                                                                       ADS34900
C --- INPUT                                                             ADS34910
C                                                                       ADS34920
C --- OUTPUT                                                            ADS34930
      SUBROUTINE ADS119 (IGOTO,NCON,XMULT,EXTRAP,ALPMAX,SLOPE,ISLOPE,ISLADS34940
     1OPG,ZRO,CTMIN,CTLMIN,GMCT,G,DG,GGL,GG1,GG2,GGU,IDG,ALPHA,OBJ,XL,X1ADS34950
     2,X2,XU,GMCTL,GMCT1,GMCT2,GMCTU,FL,F1,F2,FU)                       ADS34960
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),DG(1),GGL(1),GG1(1),GG2(1),IDG(1),GGU(1)           ADS34970
      IGOTO1=IGOTO+1                                                    ADS34980
      GO TO (10,20,30,60,90), IGOTO1                                    ADS34990
10    CONTINUE                                                          ADS35000
C     FIRST CALL.  EVALUATE FUNCTIONS AT X1=ALPHA AND CALL AGAIN.       ADS35010
      XL=0.                                                             ADS35020
      FL=OBJ                                                            ADS35030
      GMCTL=GMCT                                                        ADS35040
      CALL ADS503 (NCON,G,GGL)                                          ADS35050
      IGOTO=1                                                           ADS35060
      RETURN                                                            ADS35070
20    CONTINUE                                                          ADS35080
C     ------------------------------------------------------------------ADS35090
C     2-POINT INTERPOLATION/EXTRAPOLATION.                              ADS35100
C     ------------------------------------------------------------------ADS35110
      X1=ALPHA                                                          ADS35120
      F1=OBJ                                                            ADS35130
      GMCT1=GMCT                                                        ADS35140
      CALL ADS503 (NCON,G,GG1)                                          ADS35150
      ALPL=0.                                                           ADS35160
      ALPU=EXTRAP*X1                                                    ADS35170
      IF (GMCT1.GT.GMCTL.AND.GMCT1.GE.0.) ALPU=X1                       ADS35180
C     UPDATE BOUNDS.                                                    ADS35190
      II=2                                                              ADS35200
      CALL ADS111 (NCON,II,ISLOPG,CTMIN,CTLMIN,ZRO,XL,X1,X2,XU,DG,GGL,GGADS35210
     11,GG2,G,IDG,ALPL,ALPU)                                            ADS35220
C     MINIMUM OBJECTIVE.                                                ADS35230
      ALP=-1.                                                           ADS35240
      ICALC=1                                                           ADS35250
      IF (ISLOPE.GT.0) CALL ADS108 (ICALC,XL,FL,SLOPE,X1,F1,X2,F2,XU,FU,ADS35260
     1ZRO,ALP,YMN1)                                                     ADS35270
C     PICK BEST MOVE.                                                   ADS35280
      ALPHA=ALPU                                                        ADS35290
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS35300
      IF (F1.GT.FL.AND.GMCTL.GT.0.) ALPHA=ALPL                          ADS35310
      IF (ALPL.GT.ALPU) ALPHA=.5*(ALPL+ALPU)                            ADS35320
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS35330
C     CHECK TO BE SURE THIS ALPHA IS UNIQUE.                            ADS35340
      DXL=ABS(ALPHA)                                                    ADS35350
      DX1=ABS(ALPHA-X1)                                                 ADS35360
C     IF NOT UNIQUE AND ISLOPE.GT.0, EXIT.                              ADS35370
      DX=DMIN1(DXL,DX1)                                                 ADS35380
      IF (DX.LE.ZRO.AND.ISLOPE.GT.0) GO TO 90                           ADS35390
C     IF ALPHA IS NOT UNIQUE, PICK AVERAGE.                             ADS35400
      IF (DX.LE.ZRO) ALPHA=.5*X1                                        ADS35410
C     EVALUATE FUNCTIONS AND CALL AGAIN.                                ADS35420
      IGOTO=2                                                           ADS35430
      RETURN                                                            ADS35440
30    CONTINUE                                                          ADS35450
C     ------------------------------------------------------------------ADS35460
C     3-POINT INTERPOLATION/EXTRAPOLATION.                              ADS35470
C     ------------------------------------------------------------------ADS35480
      X2=ALPHA                                                          ADS35490
      F2=OBJ                                                            ADS35500
      GMCT2=GMCT                                                        ADS35510
      CALL ADS503 (NCON,G,GG2)                                          ADS35520
C     BOUNDS.                                                           ADS35530
      ALPL=0.                                                           ADS35540
      ALPU=EXTRAP*DMAX1(XL,X1,X2)                                       ADS35550
      IF (GMCT2.GT.GMCT1.AND.GMCT1.GT.0.) ALPU=X2                       ADS35560
      GMIN=DMIN1(GMCTL,GMCTL,GMCT2)                                     ADS35570
      IF (ALPU.GT.X2.OR.GMIN.LT.0.) GO TO 40                            ADS35580
C     ALL DESIGNS ARE INFEASIBLE AND X2 IS UPPER BOUND.                 ADS35590
C     INTERPOLATE FOR MINIMUM CONSTRAINT VIOLATION.                     ADS35600
      ICALC=2                                                           ADS35610
      CALL ADS108 (ICALC,XL,GMCTL,YP1,X1,GMCT1,X2,GMCT2,XU,GMCTU,ZRO,ALPADS35620
     1,YMN1)                                                            ADS35630
      ALPHA=.5*X1                                                       ADS35640
      IF (GMCT2.LT.GMCTL) ALPHA=.5*(X1+X2)                              ADS35650
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS35660
      GO TO 50                                                          ADS35670
40    CONTINUE                                                          ADS35680
C     UPDATE BOUNDS.                                                    ADS35690
      II=3                                                              ADS35700
      CALL ADS111 (NCON,II,ISLOPG,CTMIN,CTLMIN,ZRO,XL,X1,X2,XU,DG,GGL,GGADS35710
     11,GG2,G,IDG,ALPL,ALPU)                                            ADS35720
C     MINIMIZE OBJECTIVE.                                               ADS35730
      ALP=-1.                                                           ADS35740
      ICALC=2                                                           ADS35750
      IF (ISLOPE.GT.0) ICALC=3                                          ADS35760
      CALL ADS108 (ICALC,XL,FL,SLOPE,X1,F1,X2,F2,XU,FU,ZRO,ALP,YMN1)    ADS35770
C     PICK BEST MOVE.                                                   ADS35780
      ALPHA=ALPU                                                        ADS35790
      IF (F2.GT.FL.AND.GMCTL.GT.0.) ALPHA=ALPL                          ADS35800
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS35810
      IF (ALPL.GT.ALPU) ALPHA=.5*(ALPL+ALPU)                            ADS35820
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS35830
C     CHECK TO BE SURE ALPHA IS UNIQUE.                                 ADS35840
50    DXL=ABS(ALPHA)                                                    ADS35850
      DX1=ABS(ALPHA-X1)                                                 ADS35860
      DX2=ABS(ALPHA-X2)                                                 ADS35870
      DX=DMIN1(DXL,DX1,DX2)                                             ADS35880
C     IF NOT UNIQUE, EXIT.                                              ADS35890
      IF (DX.LE.ZRO) GO TO 90                                           ADS35900
C     EVALUATE FUNCTIONS AND CALL AGAIN.                                ADS35910
      IGOTO=3                                                           ADS35920
      RETURN                                                            ADS35930
60    CONTINUE                                                          ADS35940
C     ------------------------------------------------------------------ADS35950
C     4-POINT INTERPOLATION/EXTRAPOLATION.                              ADS35960
C     ------------------------------------------------------------------ADS35970
      IF (ISLOPE.GT.0.AND.ISLOPG.GT.0) GO TO 90                         ADS35980
      XU=ALPHA                                                          ADS35990
      FU=OBJ                                                            ADS36000
      GMCTU=GMCT                                                        ADS36010
      CALL ADS503 (NCON,G,GGU)                                          ADS36020
C     BOUNDS.                                                           ADS36030
      ALPL=0.                                                           ADS36040
      ALPU=EXTRAP*DMAX1(XL,X1,X2,XU)                                    ADS36050
      IF (GMCTU.GT.GMCT2.AND.GMCT2.GT.0.) ALPU=XU                       ADS36060
      GMIN=DMIN1(GMCTL,GMCT1,GMCT2,GMCTU)                               ADS36070
      IF (ALPU.GT.XU.OR.GMIN.LT.0.) GO TO 70                            ADS36080
C     ALL DESIGNS ARE INFEASIBLE AND XU IS UPPER BOUND.                 ADS36090
C     INTERPOLATE FOR MINIMUM CONSTRAINT VIOLATION.                     ADS36100
      ICALC=4                                                           ADS36110
      CALL ADS108 (ICALC,XL,GMCTL,YP1,X1,GMCT1,X2,GMCT2,XU,GMCTU,ZRO,ALPADS36120
     1,YMN1)                                                            ADS36130
      ALPHA=.5*(X2+XU)                                                  ADS36140
      IF (GMCT1.LT.GMCTU) ALPHA=.5*(X1+X2)                              ADS36150
      IF (GMCTL.LT.GMCT2) ALPHA=.5*X1                                   ADS36160
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS36170
      GO TO 80                                                          ADS36180
70    CONTINUE                                                          ADS36190
C     UPDATE BOUNDS.                                                    ADS36200
      II=4                                                              ADS36210
      CALL ADS111 (NCON,II,ISLOPG,CTMIN,CTLMIN,ZRO,XL,X1,X2,XU,DG,GGL,GGADS36220
     11,GG2,G,IDG,ALPL,ALPU)                                            ADS36230
C     MINIMIZE OBJECTIVE.                                               ADS36240
      ALP=-1.                                                           ADS36250
      ICALC=4                                                           ADS36260
      CALL ADS108 (ICALC,XL,FL,SLOPE,X1,F1,X2,F2,XU,FU,ZRO,ALP,YMN1)    ADS36270
C     PICK BEST MOVE.                                                   ADS36280
      ALPHA=ALPU                                                        ADS36290
      IF (FU.GT.FL.AND.GMCTL.GT.0.) ALPHA=ALPL                          ADS36300
      IF (ALP.GT.ALPL.AND.ALP.LT.ALPU) ALPHA=ALP                        ADS36310
      IF (ALPL.GT.ALPU) ALPHA=.5*(ALPL+ALPU)                            ADS36320
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS36330
C     CHECK TO SEE IF THIS ALPHA IS UNIQUE.                             ADS36340
80    DXL=ABS(ALPHA)                                                    ADS36350
      DX1=ABS(ALPHA-X1)                                                 ADS36360
      DX2=ABS(ALPHA-X2)                                                 ADS36370
      DXU=ABS(ALPHA-XU)                                                 ADS36380
      DX=DMIN1(DXL,DX1,DX2,DXU)                                         ADS36390
C     IF NOT UNIQUE, EXIT.                                              ADS36400
      IF (DX.LE.ZRO) GO TO 90                                           ADS36410
C     EVALUATE FUNCTIONS AND CALL AGAIN.                                ADS36420
      IGOTO=4                                                           ADS36430
      RETURN                                                            ADS36440
90    CONTINUE                                                          ADS36450
C     ------------------------------------------------------------------ADS36460
C     SEARCH IS COMPLETE.  EXIT.                                        ADS36470
C     ------------------------------------------------------------------ADS36480
      IGOTO=0                                                           ADS36490
      RETURN                                                            ADS36500
      END                                                               ADS36510
C --- SUBROUTINE ADS121 (JGOTO,OBJ,ALPHA,XMULT,ALPMAX,XL,X1,XU,FL,F1,FU)ADS36520
C     ------------------------------------------------------------------ADS36530
C --- ROUTINE TO FIND BRACKETS ON MINIMUM IN UNCONSTRAINED ONE-         ADS36540
C     DIMENSIONAL SEARCH                                                ADS36550
C     ------------------------------------------------------------------ADS36560
C     BY G. N. VANDERPLAATS                             NOV. 1980       ADS36570
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS36580
C                                                                       ADS36590
C --- INPUT                                                             ADS36600
C     JGOTO - INITIALLY ZERO TO INDICATE FIRST CALL.  SUBSEQUENTLY      ADS36610
C          JGOTO IS LEFT AT THE VALUE RETURNED ON PREVIOUS CALL.        ADS36620
C     ICAL - NUMBER OF FUNCTION VALUE THIS CALL CORRESPONDS TO.         ADS36630
C     OBJ - VALUE OF THE OBJECTIVE FUNCTION CORRESPONDING TO X.         ADS36640
C     ALPHA - INITIAL ESTIMATE FOR THE MOVE PARAMETER OR THE VALUE      ADS36650
C         RETURNED FROM THE LAST CALL TO ASDS013.                       ADS36660
C     XMULT - MULTIPLIER TO UPDATE THE BOUNDS ON ALPHA.                 ADS36670
C     ALPMAX - MAXIMUM ALPHA ALLOWED.  IF ALPHA EXCEEDS THIS,           ADS36680
C         THE ROUTINE TEMINNATES AND THE SOLUTION IS ASSUMED            ADS36690
C         TO BE UNBOUNDED.  SEE NOTE BELOW.                             ADS36700
C                                                                       ADS36710
C --- OUTPUT                                                            ADS36720
C     JGOTO - IF ZERO, BRACKETS HAVE BEEN FOUND.  ELSE, EVALUATE        ADS36730
C         OBJ AND CALL ASD013 AGAIN.                                    ADS36740
C     ALPHA - MOVE PARAMETER FOR NEW X.                                 ADS36750
C     XL, X1, XU - ALPHA VALUES FOR LATER USE IN THE GOLDEN             ADS36760
C         SECTION METHOD.  IGNORE X1 IF ICAL.LT.3.                      ADS36770
C     FL, F1, FU - FUNCTION VALUES CORRESPONDING TO XL, X1 AND XU.      ADS36780
C                                                                       ADS36790
C --- NOTE                                                              ADS36800
C     IF AT THE END, ALPHA EXCEEDS ALPMAX, THIS IS AN ERROR CODE.       ADS36810
C     THE OBJECTIVE CORRESPONDING TO THIS ALPHA IS NOT EVALUATED, BUT   ADS36820
C         THE XI AND FI ARRAYS CONTAIN USEFUL INFORMATION               ADS36830
C         (ASSUMING THE INPUT INITIAL ALPHA WAS LESS THAN ALPMAX).      ADS36840
C                                                                       ADS36850
C --- NOTE                                                              ADS36860
C     IF GOLDEN SECTION IS TO BE USED, IT IS SUGGESTED THAT             ADS36870
C     XMULT=(3+SQRT(5))/2.                                              ADS36880
C     THIS WILL RESULT IN EVALUATION OF F AT XL, X1 AND XU.             ADS36890
C                                                                       ADS36900
      SUBROUTINE ADS121 (JGOTO,OBJ,ALPHA,XMULT,ALPMAX,XL,X1,X2,XU,FL,F1,ADS36910
     1F2,FU)                                                            ADS36920
C                                                                       ADS36930
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      IF (JGOTO.GT.0) GO TO 10                                          ADS36940
C     FIRST CALL.                                                       ADS36950
      XL=0.                                                             ADS36960
      FL=OBJ                                                            ADS36970
      X1=0.                                                             ADS36980
      F1=OBJ                                                            ADS36990
      X2=0.                                                             ADS37000
      F2=OBJ                                                            ADS37010
C     CHECK FOR INPUT ALPHA GREATER THAN ALPMAX.                        ADS37020
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS37030
      XU=ALPHA                                                          ADS37040
      JGOTO=1                                                           ADS37050
C     EVALUATE OBJ.                                                     ADS37060
      RETURN                                                            ADS37070
10    CONTINUE                                                          ADS37080
C     WE HAVE AT LEAST TWO POINTS.  IS XU THE UPPER BOUND?              ADS37090
      FU=OBJ                                                            ADS37100
      IF (ALPHA.GE.(0.99999*ALPMAX)) GO TO 20                           ADS37110
      IF (OBJ.GE.F1) GO TO 20                                           ADS37120
C     NO.  PICK NEW ALPHA.                                              ADS37130
      ALPHA=X2+XMULT*(XU-X2)                                            ADS37140
      XL=X1                                                             ADS37150
      X1=X2                                                             ADS37160
      X2=XU                                                             ADS37170
      FL=F1                                                             ADS37180
      F1=F2                                                             ADS37190
      F2=FU                                                             ADS37200
C     CHECK FOR UNBOUNDED SOLUTION.                                     ADS37210
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS37220
      XU=ALPHA                                                          ADS37230
C     EVALUATE OBJ.                                                     ADS37240
      JGOTO=1                                                           ADS37250
      RETURN                                                            ADS37260
20    CONTINUE                                                          ADS37270
C     XU IS UPPER BOUND UNLESS ALPHA EXCEEDS ALPMAX.  SOLUTION IS       ADS37280
C     COMPLETE.                                                         ADS37290
      JGOTO=0                                                           ADS37300
      RETURN                                                            ADS37310
      END                                                               ADS37320
C --- SUBROUTINE ADS122(XL,X1,X2,XU)                                    ADS37330
C     ------------------------------------------------------------------ADS37340
C --- ROUTINE TO CALCULATE INITIAL INTERIOR POINTS, X1 AND X2 FOR       ADS37350
C     FINDING THE MINIMUM OF A FUNCTION USING THE GOLDEN SECTION METHOD.ADS37360
C     ------------------------------------------------------------------ADS37370
C     BY G. N. VANDERPLAATS                           AUG. 1980         ADS37380
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS37390
C                                                                       ADS37400
C --- INPUT                                                             ADS37410
C     XL - LOWER BOUND ON SOLUTION.                                     ADS37420
C     XU - UPPER BOUND ON SOLUTION.                                     ADS37430
C                                                                       ADS37440
C --- OUTPUT                                                            ADS37450
C     X1 - FIRST INTERMEDIATE POINT.  X1=GOLD*(XU-XL).                  ADS37460
C     X2 - SECOND INTERMEDIATE POINT.  X2=XU-GOLD*(XU-XL).              ADS37470
C          GOLD IS THE GOLDEN SECTION NUMBER (3-SQRT(5))/2.             ADS37480
C                                                                       ADS37490
      SUBROUTINE ADS122 (XL,X1,X2,XU)                                   ADS37500
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      GOLDX=0.38196601*(XU-XL)                                          ADS37510
      X1=XL+GOLDX                                                       ADS37520
      X2=XU-GOLDX                                                       ADS37530
      RETURN                                                            ADS37540
      END                                                               ADS37550
                                                                        ADS37560
C --- SUBROUTINE ADS123(IGOTO,NCON,ICAL,NPTS,OBJ,GMAX,ALPHA,ALPMAX,     ADS37570
C    *DELALP,DABALP,ZRO,XL,X1,X2,XU,FL,F1,F2,FU,GL,G1,G2,GU,            ADS37580
C    *IDL,ID1,ID2,IDU,G,GGL,GG1,GG2,GGU)                                ADS37590
C     DIMENSION G(NCON),GGL(NCON),GG1(NCON),GG2(NCON),GGU(NCON)         ADS37600
C     ------------------------------------------------------------------ADS37610
C --- ROUTINE TO PERFORM CONSTRAINED ONE-DIMENSIONAL SEARCH USING       ADS37620
C     THE GOLDEN SECTION METHOD AFTER INITIAL BOUNDS HAVE BEEN FOUND.   ADS37630
C     ------------------------------------------------------------------ADS37640
C     BY G. N. VANDERPLAATS                               NOV. 1980.    ADS37650
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS37660
C                                                                       ADS37670
C --- INPUT                                                             ADS37680
C     IGOTO - CONTROL PARAMETER.  ON FIRST CALL IGOTO=0, ON             ADS37690
C         SUBSEQUENT CALLS, IGOTO = VALUE RETURNED FROM PREVIOUS CALL.  ADS37700
C     ICAL - ANALYSIS NUMBER ASSOCIATED WITH CURRENT DESIGN.            ADS37710
C     NPTS - NUMBER OF POINTS INITIALLY AVAILABLE. MUST BE GE.2.        ADS37720
C     OBJ - OBJECTIVE FUNCTION VALUE CORRESPONDING TO X.                ADS37730
C     GMAX - MAXIMUM CONSTRAINT VALUE FOR CURRENT DESIGN.               ADS37740
C     ALPHA - INITIAL SEARCH PARAMETER.                                 ADS37750
C     ALPMAX - MAX ALPHA THAT WILL BE ALLOWED.                          ADS37760
C     DELALP - RELATIVE CONVERGENCE TOLERANCE ON ALPHA.                 ADS37770
C     DABALP - ABSOLUTE CONVERGENCE TOLERANCE ON ALPHA.  THE            ADS37780
C         CONVERGENCE CRITERION GIVING THE LARGEST TOLERANCE WILL       ADS37790
C         BE USED.                                                      ADS37800
C     ZRO - SMALL POSITIVE NUMBER TO INDICATE NUMERICAL ZERO.           ADS37810
C         INDICATOR OF MACHINE ACCURACY.                                ADS37820
C     XL-XU - UP TO FOUR ALPHA VALUES, THE FIRST AND LAST OF WHICH      ADS37830
C         BRACKET THE MINIMUM.                                          ADS37840
C     FL-FU - FUNCTIONS ASSOCIATED WITH XL-XU.                          ADS37850
C     GL-GU - GMAX VALUES ASSOCIATED WITH XL-XU.                        ADS37860
C     IDL-IDU - ANALYSIS NUMBERS ASSOCIATED WITH XL-XU.                 ADS37870
C                                                                       ADS37880
C --- OUTPUT                                                            ADS37890
C     IGOTO - CALCULATION PARAMETER.  IF IGOTO = 0, SOLUTION IS         ADS37900
C         COMPLETE.  ELSE CALCULATE OBJ AND CALL ADS123 AGAIN.          ADS37910
C     OBJ - OPTIMUM OBJECTIVE IF IGOTO = 0.                             ADS37920
C     ALPHA - OPTIMUM MOVE PARAMETER IF IGOTO = 0.                      ADS37930
C     XL-XU - UPDATED VALUES OF XL-XU.                                  ADS37940
C     FL-FU - UPDATED VALUES OF FL-FU.                                  ADS37950
C     GL-GU - UPDATED VALUES OF GL-GU.                                  ADS37960
C     IDL-IDU - UPDATED VALUES OF IDL-IDU.                              ADS37970
C                                                                       ADS37980
C --- EXTERNALS                                                         ADS37990
C     ADS122, ADS124                                                    ADS38000
C                                                                       ADS38010
      SUBROUTINE ADS123 (IGOTO,NCON,ICAL,NPTS,OBJ,GMAX,ALPHA,ALPMAX,DELAADS38020
     1LP,DABALP,ZRO,XL,X1,X2,XU,FL,F1,F2,FU,GL,G1,G2,GU,IDL,ID1,ID2,IDU,ADS38030
     2G,GGL,GG1,GG2,GGU,II,IX,NITER)                                    ADS38040
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),GGL(1),GG1(1),GG2(1),GGU(1)                        ADS38050
      IGOTO1=IGOTO+1                                                    ADS38060
      GO TO (10,20,40,50), IGOTO1                                       ADS38070
10    CONTINUE                                                          ADS38080
C     EVALUATE INTERIOR POINT(S) AS REQUIRED.                           ADS38090
C     INTERIOR POINTS.                                                  ADS38100
      CALL ADS122 (XL,X1,X2,XU)                                         ADS38110
      IF (NPTS.GT.2) GO TO 30                                           ADS38120
C     ICAL=2.  EVALUATE F1 THEN F2.                                     ADS38130
      ALPHA=X1                                                          ADS38140
      IGOTO=1                                                           ADS38150
      RETURN                                                            ADS38160
20    CONTINUE                                                          ADS38170
      ID1=ICAL                                                          ADS38180
      F1=OBJ                                                            ADS38190
      G1=GMAX                                                           ADS38200
      CALL ADS503 (NCON,G,GG1)                                          ADS38210
30    CONTINUE                                                          ADS38220
C     EVALUATE F2.                                                      ADS38230
      ALPHA=X2                                                          ADS38240
      IGOTO=2                                                           ADS38250
      RETURN                                                            ADS38260
40    CONTINUE                                                          ADS38270
      ID2=ICAL                                                          ADS38280
      F2=OBJ                                                            ADS38290
      G2=GMAX                                                           ADS38300
      CALL ADS503 (NCON,G,GG2)                                          ADS38310
C     BRACKETS HAVE BEEN FOUND AND INTERIOR POINTS EVALUATED.           ADS38320
C     UPDATE TO CONVERGENCE.                                            ADS38330
      DAB=DELALP*(XU-XL)                                                ADS38340
      IF (DAB.LT.DABALP) DAB=DABALP                                     ADS38350
      IF (DAB.LT.ZRO) DAB=ZRO                                           ADS38360
      EPS=DAB/(XU-XL)                                                   ADS38370
      AITER=-2.0780869*DLOG(EPS)                                        ADS38380
      NITER=AITER                                                       ADS38390
      II=0                                                              ADS38400
      IX=0                                                              ADS38410
50    II=II+1                                                           ADS38420
      IF (IX.NE.1) GO TO 60                                             ADS38430
C     IX=1.                                                             ADS38440
      ID1=ICAL                                                          ADS38450
      F1=OBJ                                                            ADS38460
      G1=GMAX                                                           ADS38470
      CALL ADS503 (NCON,G,GG1)                                          ADS38480
      GO TO 70                                                          ADS38490
60    IF (IX.NE.2) GO TO 70                                             ADS38500
C     IX=2.                                                             ADS38510
      ID2=ICAL                                                          ADS38520
      F2=OBJ                                                            ADS38530
      G2=GMAX                                                           ADS38540
      CALL ADS503 (NCON,G,GG2)                                          ADS38550
70    CONTINUE                                                          ADS38560
      IF (II.GT.NITER) GO TO 80                                         ADS38570
C     UPDATE BOUNDS.                                                    ADS38580
      CALL ADS124 (XL,X1,X2,XU,FL,F1,F2,FU,GL,G1,G2,GU,IDL,ID1,ID2,IDU,ZADS38590
     1RO,IX,NCON,G,GGL,GG1,GG2,GGU)                                     ADS38600
      ALPHA=X1                                                          ADS38610
      IF (IX.GT.1) ALPHA=X2                                             ADS38620
C     EVALUATE FUNCTION.                                                ADS38630
      IGOTO=3                                                           ADS38640
      RETURN                                                            ADS38650
80    CONTINUE                                                          ADS38660
C     SOLUTION IS COMPLETE.                                             ADS38670
      IGOTO=0                                                           ADS38680
      RETURN                                                            ADS38690
      END                                                               ADS38700
C --- SUBROUTINE ADS124(XL,X1,X2,XU,FL,F1,F2,FU,GL,G1,G2,GU,IDL,ID1,    ADS38710
C     ID2,IDU,GTOL,IX,NCON,G,GGL,GG1,GG2,GGU)                           ADS38720
C     DIMENSION G(NCON),GGL(NCON),GG1(NCON),GG2(NCON),GGU(NCON)         ADS38730
C     ------------------------------------------------------------------ADS38740
C --- ROUTINE TO UPDATE THE BOUNDS ON THE MINIMUM AND ESTIMATE A NEW    ADS38750
C     INTERIOR POINT FOR THE GOLDEN SECTION METHOD FOR CONSTRAINED      ADS38760
C     FUNCTION.                                                         ADS38770
C     ------------------------------------------------------------------ADS38780
C     BY G. N. VANDERPLAATS                             AUG. 1980       ADS38790
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS38800
C                                                                       ADS38810
C --- INPUT                                                             ADS38820
C     XL-XU - FOUR ALPHA VALUES FOR GOLDEN SECTION METHOD.              ADS38830
C     FL-FU - FUNCTION VALUES ASSOCIATED WITH XL-XU.                    ADS38840
C     GL-GU - GMAX VALUES ASSOCIATED WITH XL-XU.                        ADS38850
C     IDL-IDU - ANALYSIS NUMBERS ASSOCIATED WITH XL-XU.                 ADS38860
C     GTOL - TOLORANCE DETERMINING ACCEPTABILITY OF CONSTRAINT VALUES.  ADS38870
C                                                                       ADS38880
C --- OUTPUT                                                            ADS38890
C     XL-XU - UPDATED VALUES OF XL-XU.                                  ADS38900
C     FL-FU - UPDATED VALUES OF FL-FU.                                  ADS38910
C         IF IX=1, F1 IS NOT UPDATED.  IF IX=2, F2 IS NOT UPDATED.      ADS38920
C     GL-GU - UPDATED VALUES OF GL-GU.                                  ADS38930
C       IF IX=1, G1 IS NOT UPDATED.  IF IX=2, G2 IS NOT UPDATED.        ADS38940
C     IDL-IDU - UPDATED VALUES OF IDL-IDU.                              ADS38950
C         IF IX=1, ID1 IS NOT UPDATED.  IF IX=2, ID2 IS NOT UPDATED.    ADS38960
C     IX - NUMBER IDENTIFYING NEW INTERIOR POINT WHICH WAS CREATED.     ADS38970
C          IX=1 OR 2, CORRESPONDING TO X1 OR X2.  THE CORRESPONDING     ADS38980
C         FUNCTION VALUE, F1 OR F2 MUST BE CALCULATED LATER.            ADS38990
C                                                                       ADS39000
      SUBROUTINE ADS124 (XL,X1,X2,XU,FL,F1,F2,FU,GL,G1,G2,GU,IDL,ID1,ID2ADS39010
     1,IDU,GTOL,IX,NCON,G,GGL,GG1,GG2,GGU)                              ADS39020
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),GGL(1),GG1(1),GG2(1),GGU(1)                        ADS39030
C     GOLDEN SECTION FRACTION.                                          ADS39040
      GOLD=0.38196601                                                   ADS39050
      IF (IX.EQ.1) CALL ADS503 (NCON,G,GG1)                             ADS39060
      IF (IX.EQ.2) CALL ADS503 (NCON,G,GG2)                             ADS39070
C     PICK BEST POINT AND REDUCE BOUNDS.                                ADS39080
C     ARE G1 AND G2 FEASIBLE.                                           ADS39090
      IF (G1.LE.GTOL.AND.G2.LE.GTOL) GO TO 10                           ADS39100
C     NO. G1 AND/OR G2 IS INFEASIBLE.                                   ADS39110
C     IF G1 IS FEASIBLE, G2 IS NOT.  THUS X2 IS NEW UPPER BOUND.        ADS39120
      IF (G1.LT.GTOL) GO TO 30                                          ADS39130
C     G1 IS INFEASIBLE.                                                 ADS39140
C     FIND MINIMUM CONSTRAINT VIOLATION.                                ADS39150
      IF (G2.GT.G1) GO TO 30                                            ADS39160
C     G2 IS LEAST VIOLATED CONSTRAINT.  X1 IS NEW LOWER BOUND.          ADS39170
      GO TO 20                                                          ADS39180
10    IF (F2.GT.F1) GO TO 30                                            ADS39190
C     F2.LT.F1 .... X1 = NEW LOWER BOUND.                               ADS39200
20    XL=X1                                                             ADS39210
      FL=F1                                                             ADS39220
      GL=G1                                                             ADS39230
      CALL ADS503 (NCON,GG1,GGL)                                        ADS39240
      IDL=ID1                                                           ADS39250
      X1=X2                                                             ADS39260
      F1=F2                                                             ADS39270
      G1=G2                                                             ADS39280
      CALL ADS503 (NCON,GG2,GG1)                                        ADS39290
      ID1=ID2                                                           ADS39300
      X2=XU-GOLD*(XU-XL)                                                ADS39310
      IX=2                                                              ADS39320
      RETURN                                                            ADS39330
30    CONTINUE                                                          ADS39340
C     F1.LT.F2 .... X2 = NEW UPPER BOUND.                               ADS39350
      XU=X2                                                             ADS39360
      FU=F2                                                             ADS39370
      GU=G2                                                             ADS39380
      CALL ADS503 (NCON,GG2,GGU)                                        ADS39390
      IDU=ID2                                                           ADS39400
      X2=X1                                                             ADS39410
      F2=F1                                                             ADS39420
      G2=G1                                                             ADS39430
      CALL ADS503 (NCON,GG1,GG2)                                        ADS39440
      ID2=ID1                                                           ADS39450
      X1=XL+GOLD*(XU-XL)                                                ADS39460
      IX=1                                                              ADS39470
      RETURN                                                            ADS39480
      END                                                               ADS39490
C --- SUBROUTINE ADS125(JGOTO,NCON,ICAL,OBJ,GMCT,ALPHA,XMULT,ALPMAX,    ADS39500
C     *XL,X1,XU,FL,F1,FU,GL,G1,GU,IDL,ID1,IDU,G,GGL,GG1,GGU)            ADS39510
C     DIMENSION G(NCON),GGL(NCON),GG1(NCON),GGU(NCON)                   ADS39520
C     ------------------------------------------------------------------ADS39530
C --- ROUTINE TO FIND BRACKETS ON MINIMUM IN CONSTRAINED ONE-           ADS39540
C     DIMENSIONAL SEARCH                                                ADS39550
C     ------------------------------------------------------------------ADS39560
C     BY G. N. VANDERPLAATS                             NOV. 1980       ADS39570
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS39580
C                                                                       ADS39590
C --- INPUT                                                             ADS39600
C     JGOTO - INITIALLY ZERO TO INDICATE FIRST CALL.  SUBSEQUENTLY      ADS39610
C          JGOTO IS LEFT AT THE VALUE RETURNED ON PREVIOUS CALL.        ADS39620
C     ICAL - ANALYSIS NUMBER ASSOCIATED WITH CURRENT DESIGN.            ADS39630
C     OBJ - VALUE OF THE OBJECTIVE FUNCTION CORRESPONDING TO X.         ADS39640
C     GMCT - MAXIMUM OF (CONSTRAINT VALUE - CONSTRAINT VIOLATION        ADS39650
C         TOLORANCE) FOR THE CURRENT DESIGN.                            ADS39660
C     ALPHA - INITIAL ESTIMATE FOR THE MOVE PARAMETER OR THE VALUE      ADS39670
C         RETURNED FROM THE LAST CALL TO ASDS030.                       ADS39680
C     XMULT - MULTIPLIER TO UPDATE THE BOUNDS ON ALPHA.                 ADS39690
C     ALPMAX - MAXIMUM ALPHA ALLOWED.  IF ALPHA EXCEEDS THIS,           ADS39700
C         THE ROUTINE TEMINNATES AND THE SOLUTION IS ASSUMED            ADS39710
C         TO BE UNBOUNDED.  SEE NOTE BELOW.                             ADS39720
C     GGL,GG1,GGU - CONSTRAINT VECTORS.                                 ADS39730
C                                                                       ADS39740
C --- OUTPUT                                                            ADS39750
C     JGOTO - IF ZERO, BRACKETS HAVE BEEN FOUND.  ELSE, EVALUATE        ADS39760
C         OBJ AND CALL ASD013 AGAIN.                                    ADS39770
C     ALPHA - MOVE PARAMETER CORRESPONDING TO NEW X.                    ADS39780
C     XL-XU - THREE ALPHA VALUES FOR LATER USE IN GOLDEN SECTION METHOD.ADS39790
C     FL-FU - FUNCTION VALUES ASSOCIATED WITH XL-XU.                    ADS39800
C     GL-GU - GMCT VALUES ASSOCIATED WITH XL-XU.                        ADS39810
C                                                                       ADS39820
C --- NOTE                                                              ADS39830
C     IF AT THE END, ALPHA EXCEEDS ALPMAX, THIS IS AN ERROR CODE.       ADS39840
C         THE OBJECTIVE CORRESPONDING TO THIS ALPHA IS NOT EVALUATED,   ADS39850
C         BUT THE XI AND FI ARRAYS CONTAIN USEFUL INFORMATION           ADS39860
C         (ASSUMING THE INPUT INITIAL ALPHA WAS LESS THAN ALPMAX).      ADS39870
C                                                                       ADS39880
C --- NOTE                                                              ADS39890
C     IF GOLDEN SECTION IS TO BE USED, IT IS SUGGESTED THAT             ADS39900
C     XMULT=(3+SQRT(5))/2.                                              ADS39910
C     THIS WILL RESULT IN EVALUATION OF F AT XL, X1 AND XU.             ADS39920
C                                                                       ADS39930
      SUBROUTINE ADS125 (JGOTO,NCON,OBJ,GMCT,ALPHA,XMULT,ALPMAX,XL,X1,X2ADS39940
     1,XU,FL,F1,F2,FU,GL,G1,G2,GU,G,GGL,GG1,GG2,GGU)                    ADS39950
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),GGL(1),GG1(1),GG2(1),GGU(1)                        ADS39960
C                                                                       ADS39970
      IF (JGOTO.GT.0) GO TO 10                                          ADS39980
C     FIRST CALL.  INITIALIZE XL, FL, GL.                               ADS39990
      XL=0.                                                             ADS40000
      FL=OBJ                                                            ADS40010
      GL=GMCT                                                           ADS40020
      X1=0.                                                             ADS40030
      F1=FL                                                             ADS40040
      G1=GL                                                             ADS40050
      X2=0.                                                             ADS40060
      F2=FL                                                             ADS40070
      G2=GL                                                             ADS40080
      CALL ADS503 (NCON,G,GGL)                                          ADS40090
      CALL ADS503 (NCON,G,GG1)                                          ADS40100
      CALL ADS503 (NCON,G,GG2)                                          ADS40110
C     CHECK FOR INPUT ALPHA GREATER THAN ALPMAX.                        ADS40120
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS40130
C     SET XU=ALPHA.                                                     ADS40140
      XU=ALPHA                                                          ADS40150
      JGOTO=1                                                           ADS40160
C     EVALUATE OBJ.                                                     ADS40170
      RETURN                                                            ADS40180
10    CONTINUE                                                          ADS40190
      FU=OBJ                                                            ADS40200
      GU=GMCT                                                           ADS40210
      CALL ADS503 (NCON,G,GGU)                                          ADS40220
      IF (ALPHA.GE.(0.999999*ALPMAX)) GO TO 40                          ADS40230
      IF (GU.LE.0.0.AND.G1.LE.0.0) GO TO 20                             ADS40240
C     SOME CONSTRAINT IS VIOLATED.                                      ADS40250
C     IS A NEW CONSTRAINT CRITICAL OR THE MAXIMUM CONSTRAINT VALUE      ADS40260
C     INCREASING.                                                       ADS40270
      IF (GU.GT.G1) GO TO 40                                            ADS40280
C     NO. IF CONSTRAINT GU IS VIOLATED, GO INCREASE ALPHA.              ADS40290
      IF (GU.GT.0.) GO TO 30                                            ADS40300
20    CONTINUE                                                          ADS40310
C     CHECK FOR INCREASING FUNCTION VALUE.                              ADS40320
      IF (OBJ.GT.F1) GO TO 40                                           ADS40330
30    CONTINUE                                                          ADS40340
C     NO.                                                               ADS40350
C     PICK NEW ALPHA.                                                   ADS40360
      ALPHA=X2+XMULT*(XU-X2)                                            ADS40370
      XL=X1                                                             ADS40380
      X1=X2                                                             ADS40390
      X2=XU                                                             ADS40400
      FL=F1                                                             ADS40410
      F1=F2                                                             ADS40420
      F2=FU                                                             ADS40430
      GL=G1                                                             ADS40440
      G1=G2                                                             ADS40450
      G2=GU                                                             ADS40460
      CALL ADS503 (NCON,GG1,GGL)                                        ADS40470
      CALL ADS503 (NCON,GG2,GG1)                                        ADS40480
      CALL ADS503 (NCON,GGU,GG2)                                        ADS40490
C     CHECK FOR UNBOUNDED SOLUTION.                                     ADS40500
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS40510
      XU=ALPHA                                                          ADS40520
C     EVALUATE OBJ.                                                     ADS40530
      JGOTO=1                                                           ADS40540
      RETURN                                                            ADS40550
40    CONTINUE                                                          ADS40560
C     XU IS UPPER BOUND UNLESS ALPHA EXCEEDS ALPMAX.                    ADS40570
C     SOLUTION IS COMPLETE.                                             ADS40580
      JGOTO=0                                                           ADS40590
      RETURN                                                            ADS40600
      END                                                               ADS40610
C --- SUBROUTINE ADS126(IGOTO,NPTS,OBJ,ALPHA,ALPMAX,DELALP,DABALP,ZRO,X,ADS40620
C    *XL,X1,X2,XU,FL,F1,F2,FU)                                          ADS40630
C     ------------------------------------------------------------------ADS40640
C --- ROUTINE TO PERFORM UNCONSTRAINED ONE-DIMENSIONAL SEARCH USING     ADS40650
C     THE GOLDEN SECTION METHOD.                                        ADS40660
C     ------------------------------------------------------------------ADS40670
C     BY G. N. VANDERPLAATS                               NOV. 1980.    ADS40680
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS40690
C                                                                       ADS40700
C --- INPUT                                                             ADS40710
C     IGOTO - CONTROL PARAMETER.  ON FIRST CALL IGOTO=0, ON             ADS40720
C         SUBSEQUENT CALLS, IGOTO = VALUE RETURNED FROM PREVIOUS CALL.  ADS40730
C     ICAL - FUNCTION NUMBER ASSOCIATED WITH OBJ.                       ADS40740
C     OBJ - OBJECTIVE FUNCTION VALUE CORRESPONDING TO X.                ADS40750
C     ALPHA - INITIAL SEARCH PARAMETER.                                 ADS40760
C     ALPMAX - MAX ALPHA THAT WILL BE ALLOWED.                          ADS40770
C     DELALP - RELATIVE CONVERGENCE TOLERANCE ON ALPHA.                 ADS40780
C     ZRO - SMALL POSITIVE NUMBER TO INDICATE NUMERICAL ZERO.           ADS40790
C         INDICATOR OF MACHINE ACCURACY.                                ADS40800
C                                                                       ADS40810
C --- OUTPUT                                                            ADS40820
C     IGOTO - CALCULATION PARAMETER.  IF IGOTO = 0, SOLUTION IS         ADS40830
C         COMPLETE.  ELSE CALCULATE OBJ AND CALL ADS126 AGAIN.          ADS40840
C     OBJ - OPTIMUM OBJECTIVE IF IGOTO = 0.                             ADS40850
C     ALPHA - OPTIMUM MOVE PARAMETER IF IGOTO = 0.                      ADS40860
C                                                                       ADS40870
C --- NOTE                                                              ADS40880
C     IF ALPHA EXCEEDS ALPMAX ON RETURN, THE SOLUTION                   ADS40890
C     IS ASSUMED UNBOUNDED AND ALPHA=ALPMAX IS USED.                    ADS40900
C                                                                       ADS40910
C --- EXTERNALS                                                         ADS40920
C     ADS122, ADS127                                                    ADS40930
C                                                                       ADS40940
      SUBROUTINE ADS126 (IGOTO,NPTS,OBJ,ALPHA,ALPMAX,DELALP,DABALP,ZRO,XADS40950
     1L,X1,X2,XU,FL,F1,F2,FU,II,IX,NITER)                               ADS40960
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      IGOTO1=IGOTO+1                                                    ADS40970
      GO TO (10,20,40,50), IGOTO1                                       ADS40980
C     INTERIOR POINTS.                                                  ADS40990
10    CALL ADS122 (XL,X1,X2,XU)                                         ADS41000
      IF (NPTS.GT.2) GO TO 30                                           ADS41010
C     NPTS=2.  EVALUATE F1 THEN F2.                                     ADS41020
      ALPHA=X1                                                          ADS41030
      IGOTO=1                                                           ADS41040
      RETURN                                                            ADS41050
20    CONTINUE                                                          ADS41060
      F1=OBJ                                                            ADS41070
30    CONTINUE                                                          ADS41080
C     EVALUATE F2.                                                      ADS41090
      ALPHA=X2                                                          ADS41100
      IGOTO=2                                                           ADS41110
      RETURN                                                            ADS41120
40    CONTINUE                                                          ADS41130
      F2=OBJ                                                            ADS41140
C     BRACKETS HAVE BEEN FOUND AND INTERIOR POINTS EVALUATED.           ADS41150
C     UPDATE TO CONVERGENCE.                                            ADS41160
      DAB=DELALP*(XU-XL)                                                ADS41170
      IF (DAB.LT.DABALP) DAB=DABALP                                     ADS41180
      IF (DAB.LT.ZRO) DAB=ZRO                                           ADS41190
      EPS=DAB/(XU-XL)                                                   ADS41200
      AITER=-2.0780869*DLOG(EPS)                                        ADS41210
      NITER=AITER                                                       ADS41220
      II=0                                                              ADS41230
      IX=0                                                              ADS41240
50    II=II+1                                                           ADS41250
      IF (IX.EQ.1) F1=OBJ                                               ADS41260
      IF (IX.EQ.2) F2=OBJ                                               ADS41270
      IF (II.GT.NITER) GO TO 60                                         ADS41280
C     UPDATE BOUNDS.                                                    ADS41290
      CALL ADS127 (XL,X1,X2,XU,FL,F1,F2,FU,IX)                          ADS41300
      ALPHA=X1                                                          ADS41310
      IF (IX.GT.1) ALPHA=X2                                             ADS41320
C     EVALUATE FUNCTION.                                                ADS41330
      IGOTO=3                                                           ADS41340
      RETURN                                                            ADS41350
60    CONTINUE                                                          ADS41360
C     SOLUTION IS COMPLETE.                                             ADS41370
C     PICK BEST DESIGN.                                                 ADS41380
C     START WITH XL.                                                    ADS41390
      ALPHA=XL                                                          ADS41400
      OBJ=FL                                                            ADS41410
      IF (F1.GT.OBJ) GO TO 70                                           ADS41420
      ALPHA=X1                                                          ADS41430
      OBJ=F1                                                            ADS41440
70    IF (F2.GT.OBJ) GO TO 80                                           ADS41450
      ALPHA=X2                                                          ADS41460
      OBJ=F2                                                            ADS41470
80    IF (FU.GT.OBJ) GO TO 90                                           ADS41480
      ALPHA=XU                                                          ADS41490
      OBJ=FU                                                            ADS41500
90    CONTINUE                                                          ADS41510
      IGOTO=0                                                           ADS41520
      RETURN                                                            ADS41530
      END                                                               ADS41540
C --- SUBROUTINE ADS127(XL,X1,X2,XU,FL,F1,F2,FU,IX)                     ADS41550
C     ------------------------------------------------------------------ADS41560
C --- ROUTINE TO UPDATE THE BOUNDS ON THE MINIMUM AND ESTIMATE A NEW    ADS41570
C     INTERIOR POINT FOR THE GOLDEN SECTION METHOD.  UNCONSTRAINED      ADS41580
C     FUNCTION.                                                         ADS41590
C     ------------------------------------------------------------------ADS41600
C     BY G. N. VANDERPLAATS                             AUG. 1980       ADS41610
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS41620
C                                                                       ADS41630
C --- INPUT                                                             ADS41640
C     XL,FL - LOWER BOUND AND FUNCTION VALUE.                           ADS41650
C     X1,F1 - FIRST INTERMEDIATE POINT AND FUNCTION VALUE.              ADS41660
C     X2,F2 - SECOND INTERMEDIATE POINT AND FUNCTION VALUE.             ADS41670
C     XU,FU - UPPER BOUND AND FUNCTION VALUE.                           ADS41680
C                                                                       ADS41690
C --- OUTPUT                                                            ADS41700
C     XL,FL - UPDATED XL AND FL.                                        ADS41710
C     X1,F1 - UPDATED X1 AND F1 IF IX=2.  X1 ONLY IF IX=1.              ADS41720
C     X2,F2 - UPDATED X2 AND F2 IF IX=1.  X2 ONLY IF IX=2.              ADS41730
C     XU,FU - UPDATED XU AND FU.                                        ADS41740
C     IX - NUMBER IDENTIFYING NEW INTERIOR POINT WHICH WAS CREATED.     ADS41750
C          IX=1 OR 2, CORRESPONDING TO X1 OR X2.  THE CORRESPONDING     ADS41760
C         FUNCTION VALUE, F1 OR F2 MUST BE CALCULATED LATER.            ADS41770
C                                                                       ADS41780
      SUBROUTINE ADS127 (XL,X1,X2,XU,FL,F1,F2,FU,IX)                    ADS41790
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     GOLDEN SECTION FRACTION.                                          ADS41800
      GOLD=0.38196601                                                   ADS41810
C     PICK BEST POINT AND REDUCE BOUNDS.                                ADS41820
      IF (F2.GT.F1) GO TO 10                                            ADS41830
C     F2.LT.F1 .... X1 = NEW LOWER BOUND.                               ADS41840
      XL=X1                                                             ADS41850
      FL=F1                                                             ADS41860
      X1=X2                                                             ADS41870
      F1=F2                                                             ADS41880
      X2=XU-GOLD*(XU-XL)                                                ADS41890
      IX=2                                                              ADS41900
      RETURN                                                            ADS41910
10    CONTINUE                                                          ADS41920
C     F1.LT.F2 .... X2 = NEW UPPER BOUND.                               ADS41930
      XU=X2                                                             ADS41940
      FU=F2                                                             ADS41950
      X2=X1                                                             ADS41960
      F2=F1                                                             ADS41970
      X1=XL+GOLD*(XU-XL)                                                ADS41980
      IX=1                                                              ADS41990
      RETURN                                                            ADS42000
      END                                                               ADS42010
      
