C --- SUBROUTINE ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC) ADS42020
C     DIMENSION G(NCON),IDG(NCON),IC(NCOLA)                             ADS42030
C     ------------------------------------------------------------------ADS42040
C --- ROUTINE TO IDENTIFY THE NGT MOST CRITICAL CONSTRAINTS             ADS42050
C     ------------------------------------------------------------------ADS42060
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY  1983   ADS42070
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS42080
C                                                                       ADS42090
C --- INPUT                                                             ADS42100
C     NCON - NUMBER OF CONSTRAINTS.                                     ADS42110
C     NGT - NUMBER OF MOST CRITICAL CONSTRAINTS TO BE IDENTIFIED.       ADS42120
C     G - VECTOR OF CONSTRAINT VALUES.  MINIMUM DIMENSION IS NCON.      ADS42130
C     IDG - VECTOR IDENTIFYING CONSTRAINT TYPE.                         ADS42140
C         MINIMUM DIMENSION IS NCON.                                    ADS42150
C         -2 - LINEAR EQUALITY CONSTRAINT.                              ADS42160
C         -1 - NONLINEAR EQUALITY CONSTRAINT.                           ADS42170
C         0,1 - NONLINEAR INEQUALITY CONSTRAINT.                        ADS42180
C         2 - LINEAR INEQUALITY CONSTRAINT.                             ADS42190
C                                                                       ADS42200
C --- OUTPUT.                                                           ADS42210
C     IC - VECTOR CONTAINING NUMBER OF NGT ACTIVE CONSTRAINTS.          ADS42220
C                                                                       ADS42230
      SUBROUTINE ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC) ADS42240
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IDG(1),IC(1)                                       ADS42250
C                                                                       ADS42260
      NAC=0                                                             ADS42270
      NVC=0                                                             ADS42280
C     IDENTIFY ALL EQUALITY CONSTRAINTS.                                ADS42290
      N=0                                                               ADS42300
      DO 10 I=1,NGT                                                     ADS42310
10    IC(I)=0                                                           ADS42320
      DO 20 I=1,NCON                                                    ADS42330
      IF (IDG(I).GE.0) GO TO 20                                         ADS42340
      C1=CTMIN                                                          ADS42350
      IF (IDG(I).EQ.-2) C1=CTLMIN                                       ADS42360
      IF (ABS(G(I)).LE.C1) NAC=NAC+1                                    ADS42370
      IF (G(I).GT.C1) NVC=NVC+1                                         ADS42380
      N=N+1                                                             ADS42390
      IC(N)=I                                                           ADS42400
20    CONTINUE                                                          ADS42410
      IF (N.GE.NGT) RETURN                                              ADS42420
C     INEQUALITY CONSTRAINTS.                                           ADS42430
C     FIND ACTIVE AND VIOLATED G(I)                                     ADS42440
      DO 60 I=1,NCON                                                    ADS42450
      IF (N.EQ.0) GO TO 40                                              ADS42460
C     IS THIS CONSTRAINT ALREADY IN IC?                                 ADS42470
      DO 30 K=1,N                                                       ADS42480
      IF (I.EQ.IC(K)) GO TO 60                                          ADS42490
30    CONTINUE                                                          ADS42500
40    CONTINUE                                                          ADS42510
      C1=CT                                                             ADS42520
      C2=CTMIN                                                          ADS42530
      I1=IDG(I)                                                         ADS42540
      IF (IABS(I1).LT.2) GO TO 50                                       ADS42550
C     LINEAR CONSTRAINT.                                                ADS42560
      C1=CTL                                                            ADS42570
      C2=CTLMIN                                                         ADS42580
50    G1=G(I)                                                           ADS42590
      IF (G1.LT.C1) GO TO 60                                            ADS42600
      IF (G1.GE.C1.AND.G1.LE.C2) NAC=NAC+1                              ADS42610
      IF (G1.GT.C2) NVC=NVC+1                                           ADS42620
      N=N+1                                                             ADS42630
      IC(N)=I                                                           ADS42640
60    CONTINUE                                                          ADS42650
C     FIND REMAINING NGT-N MOST CRITICAL CONSTRAINTS BY MAXIMUM VALUE.  ADS42660
      M=NGT-N                                                           ADS42670
      IF (M.LE.0) RETURN                                                ADS42680
      DO 100 II=1,M                                                     ADS42690
      GMAX=-1.0E+20                                                     ADS42700
      IMX=0                                                             ADS42710
C     FIND MAXIMUM G(I).                                                ADS42720
      DO 90 I=1,NCON                                                    ADS42730
      IF (N.EQ.0) GO TO 80                                              ADS42740
      DO 70 K=1,N                                                       ADS42750
      IF (IC(K).EQ.I) GO TO 90                                          ADS42760
70    CONTINUE                                                          ADS42770
80    IF (G(I).LT.GMAX) GO TO 90                                        ADS42780
      IMX=I                                                             ADS42790
      GMAX=G(I)                                                         ADS42800
90    CONTINUE                                                          ADS42810
      IF (IMX.EQ.0) GO TO 110                                           ADS42820
      N=N+1                                                             ADS42830
      IC(N)=IMX                                                         ADS42840
100   CONTINUE                                                          ADS42850
110   CONTINUE                                                          ADS42860
      RETURN                                                            ADS42870
      END                                                               ADS42880
C --- SUBROUTINE ADS202 (NDB,ZRO,B,NRB)                                 ADS42890
C     DIMENSION B(NRB,NRB)                                              ADS42900
C     ------------------------------------------------------------------ADS42910
C --- ROUTINE TO INVERT MATRIX B AND STORE BACK IN B.                   ADS42920
C     ------------------------------------------------------------------ADS42930
C     BY G. N. VANDERPLAATS                               AUG. 1983     ADS42940
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS42950
C                                                                       ADS42960
C --- INPUT                                                             ADS42970
C     NDB - DIMENSIONED SIZE OF THE MATRIX CONTAINED IN B.              ADS42980
C     ZRO - NUMERICAL ZERO FOR COMPUTATIONAL PURPOSES.                  ADS42990
C     B - MATRIX TO BE INVERTED.                                        ADS43000
C     NRB - DIMENSION OF B.                                             ADS43010
C                                                                       ADS43020
C --- OUTPUT                                                            ADS43030
C     B - INVERSE OF B, CALCULATED BY ELIMINATION WITHOUT PIVOT SEARCH. ADS43040
C                                                                       ADS43050
      SUBROUTINE ADS202 (NDB,ZRO,B,NRB)                                 ADS43060
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION B(NRB,1)                                                ADS43070
C                                                                       ADS43080
      ZRO10=10.*ZRO                                                     ADS43090
C     BEGIN INVERSION.                                                  ADS43100
      DO 60 I=1,NDB                                                     ADS43110
C     IF ZERO ROW, IGNORE.                                              ADS43120
      BMX=B(I,I)                                                        ADS43130
      IF (ABS(BMX).GT.ZRO10) GO TO 20                                   ADS43140
C     ZERO OUT THIS ROW AND COLUMN.                                     ADS43150
      DO 10 J=1,NDB                                                     ADS43160
      B(I,J)=0.                                                         ADS43170
10    B(J,I)=0.                                                         ADS43180
      GO TO 60                                                          ADS43190
C     PIVOT.                                                            ADS43200
20    BMX=1./BMX                                                        ADS43210
      DO 30 J=1,NDB                                                     ADS43220
30    B(I,J)=BMX*B(I,J)                                                 ADS43230
      DO 50 II=1,NDB                                                    ADS43240
      IF (II.EQ.I) GO TO 50                                             ADS43250
      B1=B(II,I)                                                        ADS43260
      DO 40 J=1,NDB                                                     ADS43270
40    B(II,J)=B(II,J)-B1*B(I,J)                                         ADS43280
      B(II,I)=-BMX*B1                                                   ADS43290
50    CONTINUE                                                          ADS43300
      B(I,I)=BMX                                                        ADS43310
60    CONTINUE                                                          ADS43320
      RETURN                                                            ADS43330
      END                                                               ADS43340
C --- SUBROUTINE ADS203 (NDV,NAVC,NACS,A,NRA,B,NRB,XTMP,GTMP,ISIDC)     ADS43350
C     DIMENSION A(NRA,NAVC,B(NRB,NRB),XTMP(NDV),GTMP(NACT)              ADS43360
C     DIMENSION ISIDC(NDV)                                              ADS43370
C     ------------------------------------------------------------------ADS43380
C --- ROUTINE TO CALCULATE DELTA-X TO MOVE BACK TO CONSTRAINT           ADS43390
C     BOUNDARIES IN ROBUST METHOD OF FEASIBLE DIRECTIONS.               ADS43400
C     ------------------------------------------------------------------ADS43410
C     BY G. N. VANDERPLAATS                                AUG. 1983    ADS43420
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS43430
C                                                                       ADS43440
C --- INPUT                                                             ADS43450
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS43460
C     NAVC - NUMBER OF ACTIVE CONSTRAINTS.                              ADS43470
C     NACS - NUMBER OF ACTIVE SIDE CONSTRAINTS.                         ADS43480
C     A - MATRIX, THE COLUMNS OF WHICH CONTAIN GRADIENTS OF THE NAVC    ADS43490
C         ACTIVE CONSTRAINTS.                                           ADS43500
C     NRA - DIMENSIONED ROWS OF A.                                      ADS43510
C     B - INVERSE OF A-TRANS TIMES A.                                   ADS43520
C     NRB - DIMENSIONED ROWS OF B.                                      ADS43530
C     GTMP - CURRENT VALUES OF ACTIVE CONSTRAINTS.                      ADS43540
C     ISIDC - INTEGER ARRAY IDENTIFYING ACTIVE SIDE CONSTRAINTS.        ADS43550
C                                                                       ADS43560
C --- OUTPUT.                                                           ADS43570
C     XTMP - DELTA-X VECTOR TO DRIVE CONSTRAINTS TO ZERO.               ADS43580
C                                                                       ADS43590
      SUBROUTINE ADS203 (NDV,NAVC,NACS,A,NRA,B,NRB,XTMP,GTMP,ISIDC)     ADS43600
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(NRA,1),B(NRB,1),XTMP(1),GTMP(1),ISIDC(1)              ADS43610
C                                                                       ADS43620
      NACT=NAVC+NACS                                                    ADS43630
C     STORE G IN XTMP FOR NOW.                                          ADS43640
      DO 10 I=1,NACT                                                    ADS43650
10    XTMP(I)=GTMP(I)                                                   ADS43660
C     PUT B*XTMP IN GTMP.                                               ADS43670
      DO 30 I=1,NACT                                                    ADS43680
      GG=0.                                                             ADS43690
      DO 20 J=1,NACT                                                    ADS43700
20    GG=GG+B(I,J)*XTMP(J)                                              ADS43710
30    GTMP(I)=GG                                                        ADS43720
C     PUT -A*GTMP IN XTMP.                                              ADS43730
      DO 40 I=1,NDV                                                     ADS43740
40    XTMP(I)=0.                                                        ADS43750
      DO 50 J=1,NAVC                                                    ADS43760
      GG=GTMP(J)                                                        ADS43770
      DO 50 I=1,NDV                                                     ADS43780
50    XTMP(I)=XTMP(I)-A(I,J)*GG                                         ADS43790
C     LINEAR CONSTRAINT TERMS.                                          ADS43800
      IF (NACS.EQ.0) GO TO 80                                           ADS43810
      NN=NAVC                                                           ADS43820
      DO 70 I=1,NDV                                                     ADS43830
      II=ISIDC(I)                                                       ADS43840
      IF (II.EQ.0) GO TO 70                                             ADS43850
      IF (II.EQ.1) GO TO 60                                             ADS43860
C     LOWER BOUND IS ACTIVE.                                            ADS43870
      NN=NN+1                                                           ADS43880
      XTMP(I)=XTMP(I)+GTMP(NN)                                          ADS43890
      IF (II.LT.0) GO TO 70                                             ADS43900
C     UPPER BOUND IS ACTIVE.                                            ADS43910
60    NN=NN+1                                                           ADS43920
      XTMP(I)=XTMP(I)-GTMP(NN)                                          ADS43930
70    CONTINUE                                                          ADS43940
80    CONTINUE                                                          ADS43950
      RETURN                                                            ADS43960
      END                                                               ADS43970
C --- SUBROUTINE ADS204 (KGOTO,INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,ITMAX,  ADS43980
C    1 ICNDIR,DF,S,ZRO,DLOBJ1,DLOBJ2,DX1,DX2,XOLD,SVLB,SVUB,DELOBJ,     ADS43990
C    2 DABOBJ,ITRM,XTMP,IPRNT2,WK,IWK,G,SCFX,H,DX,DFG,IVAMET,IC,        ADS44000
C    3 IDG,ISIDC,SCFO,IFCALL,IGCALL,ALPMAX)                             ADS44010
C     DIMENSION X(NDV),S(NDV),DF(NDV),XOLD(NDV),SVLB(NDV),SVUB(NDV),    ADS44020
C    1XTMP(NDV),WK(NRWK),IWK(NRIWK),G(NCON),SCFX(NDV),DX(NDV),DFG(NDV), ADS44030
C    2H(NDV*(NDV+1)/2),IC(NCOLA),IDG(NCON),ISIDC(NDV)                   ADS44040
C     ------------------------------------------------------------------ADS44050
C --- ROUTINE TO PERFORM FLETCHER-REEVES UNCONSTRAINED MINIMIZATION.    ADS44060
C     ------------------------------------------------------------------ADS44070
C     BY G. N. VANDERPLAATS                        FEB. 1982            ADS44080
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS44090
C --- INPUT                                                             ADS44100
C      KGOTO  - CONTROL PQRAMETER.                                      ADS44110
C      INFO   - CONTROL PARAMETER.                                      ADS44120
C      ISTRAT - STRATEGY IDENTIFIER.                                    ADS44130
C      IONED  - ONE D. SEARCH IDENTIFIER.                               ADS44140
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS44150
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS44160
C      OBJ    - OBJECTIVE.                                              ADS44170
C      X      - DESIGN VARIABLES.                                       ADS44180
C      ITMAX  - MAXIMUM ITERATION NUMBER.                               ADS44190
C      ICNDIR - ITERATION NUMBER RETURN TO STEETEST DESCENT METHOD.     ADS44200
C      DF     - GRADIENT OF OBJECTIVE.                                  ADS44210
C      S      - ARRAY FOR SEARCH DIRECTION.                             ADS44220
C      ZRO    - ZERO.                                                   ADS44230
C      DLOBJ1,DLOBJ2,DX1,DX2 - PARAMETERS TO ESTIMATE INITIAL ALPHA.    ADS44240
C      XOLD   - FORMER DESIGN VARIABLE.                                 ADS44250
C      SVLB,SVUB - LOWER AND UPPER BOUND OF DESIGN VARIABLES.           ADS44260
C      DELOBJ,DABOBJ,ITRM - PARAMETERS FOR CONVERGENCE CROTERIA.        ADS44270
C      XTMP   - TEMPORARY ARRAY FOR DESIGN VARIABLES.                   ADS44280
C      IPRNT2 - PRINT CONTROL IDENTIFIER FOR OPTIMIZATION LEVEL.        ADS44290
C      WK,IWK - REAL AND INTEGER WORK ARRAY.                            ADS44300
C      G      - CONSTRAINT VALUE.                                       ADS44310
C      SCFX   - SCALE FACTOR FOR DESIGN VARIABLES.                      ADS44320
C      H,DX,DFG - ARRAY FOR VARIABLE METRIC METHOD.                     ADS44330
C      IVAMET - VARIABLE METRIC METHOD CONTROL IDENTIFIER.              ADS44340
C      IC     - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                   ADS44350
C      IDG    - CONSTRAINT TYPE IDENTIFIER.                             ADS44360
C      ISIDC  - ARRAY FOR IDENTIFY ACTIVE SIDE CONSTRAINTS.             ADS44370
C      SCFO   - SCALE FACTOR OF OBJECTIVE.                              ADS44380
C      IFCALL,IGCALL - NUMBER OF FUNCTION AND GRADIENTS CALL.           ADS44390
C --- OUTPUT                                                            ADS44400
C      OBJ    - OPTIMUM OBJECTIVE.                                      ADS44410
C      X      - OPTIMUM DESIGN VARIABLES.                               ADS44420
C      G      - CONSTRAIN VALUE CORESPOND TO  OPTIMUM DESIGN VARIABLES. ADS44430
C                                                                       ADS44440
      SUBROUTINE ADS204 (KGOTO,INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,ITMAX,ICADS44450
     1NDIR,DF,S,ZRO,DLOBJ1,DLOBJ2,DX1,DX2,XOLD,SVLB,SVUB,DELOBJ,DABOBJ,IADS44460
     2TRM,XTMP,IPRNT2,WK,IWK,G,SCFX,H,DX,DFG,IVAMET,IC,IDG,ISIDC,SCFO,IFADS44470
     3CALL,IGCALL,ALPHA,ALPSAV,DF2,OBJ1,ICOUNT,IDAB,IDEL,ITER,LGOTO,ALPMADS44480
     4AX,SCFG)                                                          ADS44490
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),S(1),DF(1),XOLD(1),SVLB(1),SVUB(1),XTMP(1),WK(1),IWADS44500
     1K(1),G(1),SCFX(1),DX(1),DFG(1),H(1),IC(1),IDG(1),ISIDC(1),SCFG(1) ADS44510
C                                                                       ADS44520
      KGOTO1=KGOTO+1                                                    ADS44530
      GO TO (10,30,70), KGOTO1                                          ADS44540
10    CONTINUE                                                          ADS44550
C     INITIAL CALL.                                                     ADS44560
C     PRINT AS REQUIRED.                                                ADS44570
      IF (IPRNT2.GT.0.AND.ISTRAT.EQ.0) WRITE (6,190)                    ADS44580
      IF (IPRNT2.GT.0.AND.ICNDIR.NE.1) WRITE (6,200)                    ADS44590
      IF (IPRNT2.GT.0.AND.ICNDIR.EQ.1) WRITE (6,210)                    ADS44600
C     SET UP PARAMETERS.                                                ADS44610
      ITER=0                                                            ADS44620
      ICOUNT=0                                                          ADS44630
      IDEL=0                                                            ADS44640
      IDAB=0                                                            ADS44650
      ALPSAV=0.                                                         ADS44660
      IF (IPRNT2.EQ.0) GO TO 20                                         ADS44670
      OBJW=OBJ                                                          ADS44680
      IF (ISTRAT.EQ.0) OBJW=OBJ/SCFO                                    ADS44690
      WRITE (6,220) OBJW                                                ADS44700
      ICALC=1                                                           ADS44710
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS44720
      ICALC=4                                                           ADS44730
      IF (NCON.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                    ADS44740
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS44750
C     ------------------------------------------------------------------ADS44760
C     OPTIMIZATION.  BEGIN NEW ITERATION.                               ADS44770
C     ------------------------------------------------------------------ADS44780
20    CONTINUE                                                          ADS44790
      ITER=ITER+1                                                       ADS44800
      IF (IPRNT2.GT.2) WRITE (6,230) ITER                               ADS44810
      OBJ1=OBJ                                                          ADS44820
C     CALCULATE GRADIENT OF OBJECTIVE.                                  ADS44830
      KGOTO=1                                                           ADS44840
      INFO=2                                                            ADS44850
      IF (ISTRAT.EQ.0.AND.ITER.EQ.1) GO TO 30                           ADS44860
      IF (ITER.GT.1.AND.ALPHA.LE.ZRO) GO TO 30                          ADS44870
      RETURN                                                            ADS44880
C                                                                       ADS44890
30    CONTINUE                                                          ADS44900
      OBJ=OBJ1                                                          ADS44910
      ICALC=3                                                           ADS44920
      IF (IPRNT2.GT.4) CALL ADS402 (ICALC,DF,NDV,SCFX)                  ADS44930
C     IDENTIFY ACTIVE SIDE CONSTRAINTS.                                 ADS44940
      CALL ADS509 (NDV,X,SVLB,SVUB,ZRO,ISIDC,NACS)                      ADS44950
C     SEARCH DIRECTION.                                                 ADS44960
      ICOUNT=ICOUNT+1                                                   ADS44970
      IF (ICOUNT.GT.ICNDIR) ICOUNT=1                                    ADS44980
C     STEEPEST DESCENT.                                                 ADS44990
40    IF (ICOUNT.EQ.1) JVAMET=0                                         ADS45000
      IF (ICOUNT.EQ.1) CALL ADS236 (NDV,JVAMET,DF,DX,DFG,H,DF2,S,S2,SCALADS45010
     1S,ZRO,ISIDC)                                                      ADS45020
C     CONJUGATE DIRECTION.                                              ADS45030
      IF (ICOUNT.GT.1) CALL ADS236 (NDV,IVAMET,DF,DX,DFG,H,DF2,S,S2,SCALADS45040
     1S,ZRO,ISIDC)                                                      ADS45050
C     IF ALL VARIABLES ARE EQUAL TO THEIR LIMIT VALUES , THEN TERMINATE.ADS45060
      IF (DF2.LE.ZRO) GO TO 100                                         ADS45070
      ICALC=2                                                           ADS45080
      IF (IPRNT2.GT.3) CALL ADS402 (ICALC,S,NDV,SCFX)                   ADS45090
C     ONE-DIMENSIONAL SEARCH.                                           ADS45100
C     SLOPE OF F(ALPHA) FUNCTION.                                       ADS45110
      CALL ADS501 (NDV,DF,S,SLOPF)                                      ADS45120
      ZROM=-ZRO                                                         ADS45130
      IF (SLOPF.LT.ZROM.OR.ICOUNT.EQ.1) GO TO 50                        ADS45140
C     SLOPE IS BAD (POSITIVE).  GO USE STEEPEST DESCENT.                ADS45150
      ICOUNT=1                                                          ADS45160
      GO TO 40                                                          ADS45170
50    CONTINUE                                                          ADS45180
C     IF FIRST ITERATION, ESTIMATE ALPHA.                               ADS45190
      IF (ITER.GT.1) GO TO 60                                           ADS45200
      CALL ADS213 (OBJ,SLOPF,DLOBJ1,DLOBJ2,ALPHA)                       ADS45210
      CALL ADS214 (NDV,X,S,DX1,DX2,ALP1)                                ADS45220
      IF (ALP1.LT.ALPHA) ALPHA=ALP1                                     ADS45230
60    CONTINUE                                                          ADS45240
      ALP1=ALPSAV                                                       ADS45250
      ALPSAV=ALPHA                                                      ADS45260
      ALPHA=.5*(ALPHA+ALP1)                                             ADS45270
      IF (ALPHA.LT.ZRO) ALPHA=ZRO                                       ADS45280
C     ALPMAX.                                                           ADS45290
      CALL ADS112 (NDV,X,S,SVLB,SVUB,ZRO,ALPMAX)                        ADS45300
      IF (ALPMAX.LE.ZRO) GO TO 100                                      ADS45310
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS45320
C     ------------------------------------------------------------------ADS45330
C     DO ONE-DIMENSIONAL SEARCH.                                        ADS45340
C     ------------------------------------------------------------------ADS45350
      LGOTO=0                                                           ADS45360
      IF (IPRNT2.GT.3) WRITE (6,120) ALPHA                              ADS45370
70    CONTINUE                                                          ADS45380
      NACG=0                                                            ADS45390
      CALL ADS001 (LGOTO,ISTRAT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,   ADS45400
     1NACG,IC,WK,IWK,SLOPF,ALPHA,ALPMAX)                                ADS45410
      IF (LGOTO.EQ.0) GO TO 80                                          ADS45420
C     EVALUATE OBJECTIVE FUNCTION.                                      ADS45430
      KGOTO=2                                                           ADS45440
      INFO=1                                                            ADS45450
      RETURN                                                            ADS45460
C                                                                       ADS45470
80    CONTINUE                                                          ADS45480
C     ONE-DIMENSIONAL SEARCH IS COMPLETE.                               ADS45490
      OBJW=OBJ                                                          ADS45500
      IF (ISTRAT.EQ.0) OBJW=OBJ/SCFO                                    ADS45510
      IF (IPRNT2.EQ.2) WRITE (6,240) ITER,OBJW                          ADS45520
      IF (IPRNT2.EQ.3) WRITE (6,130) OBJW                               ADS45530
      IF (IPRNT2.GE.4) WRITE (6,250) ALPHA,OBJW                         ADS45540
      IF (IPRNT2.LT.2) GO TO 90                                         ADS45550
      ICALC=1                                                           ADS45560
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS45570
      IF (IPRNT2.LT.3) GO TO 90                                         ADS45580
      ICALC=4                                                           ADS45590
      IF (NCON.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                    ADS45600
90    CONTINUE                                                          ADS45610
C     ------------------------------------------------------------------ADS45620
C     CHECK CONVERGENCE.                                                ADS45630
C     ------------------------------------------------------------------ADS45640
      CALL ADS211 (OBJ,OBJ1,DELOBJ,DABOBJ,DEL1,DAB1,IDEL,IDAB)          ADS45650
      IF (ITER.GE.ITMAX) GO TO 100                                      ADS45660
      IF (IDEL.LT.ITRM.AND.IDAB.LT.ITRM) GO TO 20                       ADS45670
100   CONTINUE                                                          ADS45680
C     ------------------------------------------------------------------ADS45690
C     OPTIMIZATION IS COMPLETE.                                         ADS45700
C     ------------------------------------------------------------------ADS45710
      IF (IPRNT2.EQ.0) GO TO 110                                        ADS45720
      WRITE (6,260)                                                     ADS45730
      OBJW=OBJ                                                          ADS45740
      IF (ISTRAT.EQ.0) OBJW=OBJ/SCFO                                    ADS45750
      WRITE (6,130) OBJW                                                ADS45760
      ICALC=1                                                           ADS45770
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS45780
      ICALC=4                                                           ADS45790
      IF (NCON.NE.0) CALL ADS402 (ICALC,G,NCON,SCFG)                    ADS45800
      WRITE (6,140)                                                     ADS45810
      IF (ITER.GE.ITMAX) WRITE (6,150)                                  ADS45820
      IF (IDEL.GE.ITRM) WRITE (6,160) ITRM                              ADS45830
      IF (IDAB.GE.ITRM) WRITE (6,170) ITRM                              ADS45840
      IF (DF2.LE.ZRO) WRITE (6,180)                                     ADS45850
110   CONTINUE                                                          ADS45860
      KGOTO=0                                                           ADS45870
      RETURN                                                            ADS45880
C     ------------------------------------------------------------------ADS45890
C     FORMATS                                                           ADS45900
C     ------------------------------------------------------------------ADS45910
120   FORMAT (/,5X,'PROPOSED ALPHA =',1PE12.5)                          ADS45920
130   FORMAT (/,5X,'OBJECTIVE =',1PE12.5)                               ADS45930
140   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS45940
150   FORMAT (/,5X,'NUMBER OF ITERATIONS = ITMAX')                      ADS45950
160   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS45960
     1' CONSECUTIVE ITERATIONS')                                        ADS45970
170   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS45980
     1' CONSECUTIVE ITERATIONS')                                        ADS45990
180   FORMAT (/,5X,'KUHN-TUCKER CONDITIONS ARE SATISFIED')              ADS46000
190   FORMAT ('1')                                                      ADS46010
200   FORMAT (//,5X,75('-'),/,15X,                                      ADS46020
     1'IOPT = 1; FLETCHER-REEVES CONJUGATE DIIRECTIONS METHOD',/,       ADS46030
     25X,75('-'))                                                       ADS46040
210   FORMAT (//,5X,75('-'),/,21X,                                      ADS46050
     1'IOPT = 1; ICNDIR = 1; STEEPEST DESCENT METHOD',/,5X,75('-'))     ADS46060
220   FORMAT (/,5X,'INITIAL DESIGN',//,5X,'OBJ =',1PE12.5)              ADS46070
230   FORMAT (//,5X,'-- BEGIN ITERATION',I5)                            ADS46080
240   FORMAT (//,5X,'-- ITERATION',I5,3X,'OBJ =',1PE12.5)               ADS46090
250   FORMAT (/,5X,'CALCULATED ALPHA =',1PE12.5,//,5X,'OBJECTIVE =',    ADS46100
     11PE12.5)                                                          ADS46110
260   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS')                     ADS46120
      END                                                               ADS46130
C --- SUBROUTINE ADS205 (IGOTO,INFO,IONED,NDV,NCON,ITMAX,ICNDIR,NGT,    ADS46140
C    1 NAC,NVC,IPRNT2,ITRMOP,OBJ,ZRO,CT,CTMIN,CTL,CTLMIN,STOL,          ADS46150
C    2 THETAZ,DELOBJ,DABOBJ,SCFO,SLOPF,ALPHA,ALPMAX,X,SVLB,SVUB,DF,S,G  ADS46160
C    3 ,IC,IDG,DG,THETA,C,ISIDC,XOLD,XTMP,GGL,SCFX,ITMP,SCFG,A,NRA,NCOLAADS46170
C    4 ,B,NRB,WK,IWK)                                                   ADS46180
C     ------------------------------------------------------------------ADS46190
C --- ROUTINE TO PERFORM CONSTRAINED MINIMIZATION BY THE METHOD OF      ADS46200
C     FEASIBLE DIRECTIONS.                                              ADS46210
C     ------------------------------------------------------------------ADS46220
C     BY G. N. VANDERPLAATS                                  FEB. 1982  ADS46230
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS46240
C --- INPUT                                                             ADS46250
C     IDG(I) - CONSTRAINT TYPE INDICATOR.                               ADS46260
C       0,1 - NONLINEAR INEQUALITY CONSTRAINT.                          ADS46270
C         2 - LINEAR INEQUALITY CONSTRAINT.                             ADS46280
C        -1 - NONLINEAR EQUALITY CONSTRAINT.                            ADS46290
C        -2 - LINEAR EQUALITY CONSTRAINT.                               ADS46300
C                                                                       ADS46310
C --- OUTPUT                                                            ADS46320
C                                                                       ADS46330
      SUBROUTINE ADS205 (IGOTO,INFO,IONED,NDV,NCON,ITMAX,ICNDIR,NGT,NAC,ADS46340
     1NVC,IPRNT2,ITRMOP,OBJ,ZRO,CT,CTMIN,CTL,CTLMIN,STOL,THETAZ,DELOBJ,DADS46350
     2ABOBJ,SCFO,SLOPF,ALPHA,ALPMAX,X,SVLB,SVUB,DF,S,G,IC,IDG,DG,THETA,CADS46360
     3,ISIDC,XOLD,XTMP,GGL,SCFX,ITMP,SCFG,A,NRA,NCOLA,B,NRB,WK,IWK,     ADS46370
     4ICOUNT,IDAB,IDEL,ITER,LGOTO,NACS,NAVC,NGRAD,ISTRAT,TEMP)          ADS46380
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),DF(1),S(1),G(1),IC(1),IDG(1),DG(1),ADS46390
     1THETA(1),C(1),ISIDC(1),XOLD(1),XTMP(1),GGL(1),SCFX(1),ITMP(1),SCFGADS46400
     2(1),A(NRA,1),B(NRB,1),WK(1),IWK(1),TEMP(10)                       ADS46410
C                                                                       ADS46420
C     GET INFO FROM TEMP                                                ADS46430
      ALPSAV=TEMP(1)                                                    ADS46440
      BETA=TEMP(2)                                                      ADS46450
      DF2=TEMP(3)                                                       ADS46460
      OBJ1=TEMP(4)                                                      ADS46470
      PHI=TEMP(5)                                                       ADS46480
      DLOBJ1=TEMP(6)                                                    ADS46490
      DLOBJ2=TEMP(7)                                                    ADS46500
      DX1=TEMP(8)                                                       ADS46510
      DX2=TEMP(9)                                                       ADS46520
	GMAX1=TEMP(10)                             !PUT THIS LINE BY SONGSAK NOV06
      ZRO10=10.*ZRO                                                     ADS46530
      IGOTO1=IGOTO+1                                                    ADS46540
      GO TO (10,60,240), IGOTO1                                         ADS46550
10    CONTINUE                                                          ADS46560
C     ------------------------------------------------------------------ADS46570
C     INITIAL CALL.  SET UP PARAMETERS.                                 ADS46580
C     ------------------------------------------------------------------ADS46590
C     PRINT AS REQUIRED.                                                ADS46600
      IF (ISTRAT.EQ.0.AND.IPRNT2.GT.0) WRITE (6,350)                    ADS46610
      IF (IPRNT2.GT.0) WRITE (6,460)                                    ADS46620
C     INITIALIZE PARAMETERS.                                            ADS46630
      ITER=0                                                            ADS46640
      ICOUNT=0                                                          ADS46650
      ALPSAV=0.                                                         ADS46660
      PHI=5.                                                            ADS46670
      IDEL=0                                                            ADS46680
      IDAB=0                                                            ADS46690
      NGRAD=1                                                           ADS46700
C     INITIALIZE S-VECTOR AND ISIDC-VECTOR TO ZERO.                     ADS46710
      DO 20 I=1,NDV                                                     ADS46720
      ISIDC(I)=0                                                        ADS46730
20    S(I)=0.                                                           ADS46740
      IF (IPRNT2.EQ.0) GO TO 30                                         ADS46750
      OBJA=OBJ/SCFO                                                     ADS46760
      WRITE (6,560) OBJA                                                ADS46770
      ICALC=1                                                           ADS46780
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS46790
      ICALC=5                                                           ADS46800
      CALL ADS402 (ICALC,SVLB,NDV,SCFX)                                 ADS46810
      ICALC=6                                                           ADS46820
      CALL ADS402 (ICALC,SVUB,NDV,SCFX)                                 ADS46830
      ICALC=4                                                           ADS46840
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS46850
30    CONTINUE                                                          ADS46860
C     ------------------------------------------------------------------ADS46870
C     OPTIMIZATION.  BEGIN NEW ITERATION.                               ADS46880
C     ------------------------------------------------------------------ADS46890
40    CONTINUE                                                          ADS46900
      ITER=ITER+1                                                       ADS46910
      OBJ1=OBJ                                                          ADS46920
      IF (IPRNT2.GT.2) WRITE (6,570) ITER                               ADS46930
      IF (NGRAD.EQ.0) GO TO 60                                          ADS46940
C     ------------------------------------------------------------------ADS46950
C     CALCULATE GRADIENT OF OBJECTIVE AND ACTIVE AND/OR VIOLATED        ADS46960
C     CONSTRAINTS.                                                      ADS46970
C     ------------------------------------------------------------------ADS46980
      IF (ITER.GT.1) GO TO 50                                           ADS46990
C     GRADIENTS ARE ALREADY AVAILABLE.  SHIFT A AND IC AS REQUIRED.     ADS47000
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,CTLMIN,NAC,NVC)  ADS47010
      GO TO 60                                                          ADS47020
50    CONTINUE                                                          ADS47030
C     DETERMINE THE ACTIVE CONSTRAINT SET.                              ADS47040
C     GENERAL CONSTRAINTS, G.                                           ADS47050
      IF (NGT.LT.NCON) CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIADS47060
     1N,NAC,NVC)                                                        ADS47070
      IF (NGT.EQ.NCON) CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,ADS47080
     1CTLMIN,NAC,NVC)                                                   ADS47090
C     EVALUATE GRADIENTS OF OBJECTIVE AND ALL ACTIVE/VIOLATED           ADS47100
C     CONSTRAINTS.                                                      ADS47110
      IGOTO=1                                                           ADS47120
      INFO=2                                                            ADS47130
C     PUT INFO IN TEMP FOR TRANSFER                                     ADS47140
      TEMP(1)=ALPSAV                                                    ADS47150
      TEMP(2)=BETA                                                      ADS47160
      TEMP(3)=DF2                                                       ADS47170
      TEMP(4)=OBJ1                                                      ADS47180
      TEMP(5)=PHI                                                       ADS47190
      TEMP(6)=DLOBJ1                                                    ADS47200
      TEMP(7)=DLOBJ2                                                    ADS47210
      TEMP(8)=DX1                                                       ADS47220
      TEMP(9)=DX2                                                       ADS47230
      RETURN                                                            ADS47240
60    CONTINUE                                                          ADS47250
      NAVC=NAC+NVC                                                      ADS47260
C     IDENTIFY ACTIVE SIDE CONSTRAINTS.                                 ADS47270
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS47280
      IF (IPRNT2.LT.4) GO TO 90                                         ADS47290
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS47300
      WRITE (6,470) NAC,NVC                                             ADS47310
      IF (NAVC.GT.0) WRITE (6,480)                                      ADS47320
      IF (NAVC.GT.0) WRITE (6,490) (IC(I),I=1,NAVC)                     ADS47330
      WRITE (6,500) NACS                                                ADS47340
      IF (NACS.EQ.0) GO TO 90                                           ADS47350
      N=0                                                               ADS47360
      DO 80 I=1,NDV                                                     ADS47370
      IF (ISIDC(I).EQ.0) GO TO 80                                       ADS47380
      IF (ISIDC(I).EQ.1) GO TO 70                                       ADS47390
C     LOWER BOUND IS ACTIVE.                                            ADS47400
      N=N+1                                                             ADS47410
      ITMP(N)=-I                                                        ADS47420
70    IF (ISIDC(I).LT.0) GO TO 80                                       ADS47430
C     UPPER BOUND IS ACTIVE.                                            ADS47440
      N=N+1                                                             ADS47450
      ITMP(N)=I                                                         ADS47460
80    CONTINUE                                                          ADS47470
      WRITE (6,510)                                                     ADS47480
      WRITE (6,490) (ITMP(I),I=1,NACS)                                  ADS47490
90    CONTINUE                                                          ADS47500
      IF (IPRNT2.LT.5) GO TO 110                                        ADS47510
C     WRITE GRADIENT OF OBJECTIVE.                                      ADS47520
      ICALC=3                                                           ADS47530
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS47540
      IF (NAVC.EQ.0) GO TO 110                                          ADS47550
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS47560
      ISCAL=0                                                           ADS47570
      DO 100 J=1,NAVC                                                   ADS47580
      WRITE (6,520) IC(J)                                               ADS47590
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS47600
100   CONTINUE                                                          ADS47610
110   CONTINUE                                                          ADS47620
C     ------------------------------------------------------------------ADS47630
C     SEARCH DIRECTION.                                                 ADS47640
C     ------------------------------------------------------------------ADS47650
      NACT=NAVC+NACS                                                    ADS47660
      IF (NACT.GT.0) GO TO 130                                          ADS47670
C     NO ACTIVE OR VIOLATED CONSTRAINTS.  FIND STEEPEST DESCENT OR      ADS47680
C     FLETCHER-REEVES CONJUGATE DIRECTION.                              ADS47690
      ICOUNT=ICOUNT+1                                                   ADS47700
      IF (IDEL.EQ.1.OR.IDAB.EQ.1) ICOUNT=1                              ADS47710
      IF (ICOUNT.GT.ICNDIR) ICOUNT=1                                    ADS47720
C     STEEPEST DESCENT.                                                 ADS47730
120   IF (ICOUNT.EQ.1) SCALS=0.                                         ADS47740
      IF (ICOUNT.EQ.1) DF2=1.                                           ADS47750
C     CONJUGATE DIRECTION.                                              ADS47760
      ICALC=3                                                           ADS47770
      CALL ADS236 (NDV,ICALC,DF,DM1,DM2,DM3,DF2,S,S2,SCALS,ZRO,ISIDC)   ADS47780
C     SLOPE OF F(ALPHA) FUNCTION.                                       ADS47790
      CALL ADS501 (NDV,DF,S,SLOPF)                                      ADS47800
      ZROM=-ZRO                                                         ADS47810
      IF (SLOPF.LT.ZROM.OR.ICOUNT.EQ.1) GO TO 160                       ADS47820
C     SLOPE IS BAD (POSITIVE).  GO USE STEEPEST DESCENT.                ADS47830
      ICOUNT=1                                                          ADS47840
      GO TO 120                                                         ADS47850
130   CONTINUE                                                          ADS47860
C     CONSTRAINTS ARE ACTIVE/VIOLATED.  FIND USABLE-FEASIBLE            ADS47870
C     SEARCH DIRECTION.                                                 ADS47880
      ICOUNT=0                                                          ADS47890
      IF (NVC.EQ.0) PHI=5.                                              ADS47900
      PHI=10.*PHI                                                       ADS47910
      IF (PHI.GT.1000.) PHI=1000.                                       ADS47920
C     CALCULATE SCALE FACTOR FOR GRADIENT OF OBJECTIVE.                 ADS47930
      N2=1                                                              ADS47940
      CALL ADS207 (NDV,N2,DF,NDV,ZRO,SCALF)                             ADS47950
C     CALCULATE NORMALIZATION FACTORS FOR CONSTRAINTS.                  ADS47960
C     STORE FACTORS IN GGL.                                             ADS47970
      IF (NAVC.GT.0) CALL ADS207 (NDV,NAVC,A,NRA,ZRO,GGL)               ADS47980
C     CALCULATE USABLE-FEASIBLE SEARCH DIRECTION.                       ADS47990
      CALL ADS217 (NDV,NAC,NVC,NACS,CT,CTMIN,CTLMIN,THETAZ,ZRO,PHI,GGL,GADS48000
     1,DF,A,NRA,IDG,IC,ISIDC,ITMP,B,NRB,SCALF,S,THETA,C,BETA,SCALS)     ADS48010
      IF (IPRNT2.GT.3) WRITE (6,530) BETA,SCALS                         ADS48020
C     ------------------------------------------------------------------ADS48030
C     CHECK KUHN-TUCKER CONDITIONS.                                     ADS48040
C     ------------------------------------------------------------------ADS48050
C     IF NDV SIDE CONSTRAINTS ARE ACTIVE AND BETA=0, TERMINATE.         ADS48060
      IF (NACS.EQ.NDV.AND.BETA.LE.ZRO10) GO TO 300                      ADS48070
      IF (SCALS.GT.STOL.AND.BETA.GT.ZRO10) GO TO 150                    ADS48080
C     EITHER SCALS OR BETA INDICATES K-T CONDITIONS ARE SATISFIED.      ADS48090
C     CHECK CONSTRAINT VALUES.                                          ADS48100
      ALPHA=0.                                                          ADS48110
      DO 140 I=1,NAVC                                                   ADS48120
      II=IC(I)                                                          ADS48130
      JJ=IABS(IDG(II))                                                  ADS48140
      CC=CTMIN                                                          ADS48150
      IF (JJ.GT.1) CC=CTLMIN                                            ADS48160
C     IF G.GT.CC, GO CHECK CONVERGENCE.                                 ADS48170
      IF (ABS(G(II)).GT.CC) GO TO 250                                   ADS48180
140   CONTINUE                                                          ADS48190
C     K-T CONDITIONS ARE SATISFIED.  GO EXIT.                           ADS48200
      GO TO 300                                                         ADS48210
150   CONTINUE                                                          ADS48220
C     SLOPE OF F(ALPHA) FUNCTION.                                       ADS48230
      CALL ADS501 (NDV,DF,S,SLOPF)                                      ADS48240
160   CONTINUE                                                          ADS48250
      IF (IPRNT2.LT.4) GO TO 170                                        ADS48260
      ICALC=2                                                           ADS48270
      CALL ADS402 (ICALC,S,NDV,SCFX)                                    ADS48280
170   CONTINUE                                                          ADS48290
C     ------------------------------------------------------------------ADS48300
C     ESTIMATE ALPHA.                                                   ADS48310
C     ------------------------------------------------------------------ADS48320
      IF (ITER.GT.1) GO TO 180                                          ADS48330
C     FIRST ITERATION ESTIMATE.                                         ADS48340
      CALL ADS213(OBJ,SLOPF,DLOBJ1,DLOBJ2,ALPHA)                        ADS48350
      CALL ADS214(NDV,X,S,DX1,DX2,ALP1)                                 ADS48360
      IF(ALP1.LT.ALPHA) ALPHA=ALP1                                      ADS48370
      ALPSAV=ALPHA                                                      ADS48380
180   CONTINUE                                                          ADS48390
      OBJ1=OBJ                                                          ADS48400
      ALP1=ALPHA                                                        ADS48410
      ALPHA=.5*(ALPHA+ALPSAV)                                           ADS48420
      ALPSAV=ALP1                                                       ADS48430
      IF (ALPHA.LE.0.001) ALPHA=0.001                                   ADS48440
C     ESTIMATE ALPMAX TO ENCOUNTER NDV SIDE CONSTRAINTS.                ADS48450
      CALL ADS112 (NDV,X,S,SVLB,SVUB,ZRO,ALPMAX)                        ADS48460
      IF (ALPMAX.LE.ZRO) GO TO 300                                      ADS48470
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS48480
C     ------------------------------------------------------------------ADS48490
C     CALCULATE DG/DALPHA VECTOR AND MOVE TO ENCOUNTER NEW CONSTRAINT,  ADS48500
C     OR OVERCOME CONSTRAINT VIOLATION.                                 ADS48510
C     ------------------------------------------------------------------ADS48520
      ALP1=5.*ALPHA                                                     ADS48530
      ALP2=0.                                                           ADS48540
      II=0                                                              ADS48550
      DO 220 I=1,NCON                                                   ADS48560
      DG(I)=0.                                                          ADS48570
      DO 190 J=1,NGT                                                    ADS48580
      J1=J                                                              ADS48590
      IF (IC(J).EQ.I) GO TO 200                                         ADS48600
190   CONTINUE                                                          ADS48610
C     THIS GRADIENT IS NOT AVAILABLE.                                   ADS48620
      GO TO 220                                                         ADS48630
200   CONTINUE                                                          ADS48640
C     GRADIENT IS AVAILABLE.  CALCULATE DG(I)                           ADS48650
      DO 210 J=1,NDV                                                    ADS48660
210   DG(I)=DG(I)+A(J,J1)*S(J)                                          ADS48670
      IF (ABS(DG(I)).LE.0.001) GO TO 220                                ADS48680
      ALP=-G(I)/DG(I)                                                   ADS48690
      IF (ALP.LT.ALP1.AND.ALP.GT.ZRO.AND.J1.GT.NAVC) ALP1=ALP           ADS48700
      IF (ALP.GT.ALP2.AND.J1.LE.NAVC) ALP2=ALP                          ADS48710
      II=1                                                              ADS48720
220   CONTINUE                                                          ADS48730
      IF(NVC.GT.0.AND.ALP2.GT.0.) ALPHA=ALP2                            ADS48740
      IF (II.EQ.1.AND.ALP1.LT.ALP2) ALPHA=ALP1                          ADS48750
C     LIMIT FIRST ESTIMATE TO MOVE TO SIDE CONSTRAINT.                  ADS48760
      ALP=1.0E+20                                                       ADS48770
      DO 230 I=1,NDV                                                    ADS48780
      IF (ABS(S(I)).LE.ZRO) GO TO 230                                   ADS48790
      D1=(SVLB(I)-X(I))/S(I)                                            ADS48800
      IF (D1.GT.ZRO.AND.D1.LT.ALP) ALP=D1                               ADS48810
      D1=(SVUB(I)-X(I))/S(I)                                            ADS48820
      IF (D1.GT.ZRO.AND.D1.LT.ALP) ALP=D1                               ADS48830
230   CONTINUE                                                          ADS48840
      IF (ALP.LT.ALPHA) ALPHA=ALP                                       ADS48850
C     ------------------------------------------------------------------ADS48860
C     DO ONE-DIMENSIONAL SEARCH.                                        ADS48870
C     ------------------------------------------------------------------ADS48880
      LGOTO=0                                                           ADS48890
      IF (IPRNT2.GE.4) WRITE (6,540) ALPHA                              ADS48900
240   CONTINUE                                                          ADS48910
      CALL ADS001 (LGOTO,ISTRAT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,   ADS48920
     1NAVC,IC,WK,IWK,SLOPF,ALPHA,ALPMAX)                                ADS48930
      IF (LGOTO.EQ.0) GO TO 250                                         ADS48940
C     CALCULATE FUNCTION VALUES.                                        ADS48950
      IGOTO=2                                                           ADS48960
      INFO=1                                                            ADS48970
C     PUT INFO IN TEMP FOR TRANSFER                                     ADS48980
      TEMP(1)=ALPSAV                                                    ADS48990
      TEMP(2)=BETA                                                      ADS49000
      TEMP(3)=DF2                                                       ADS49010
      TEMP(4)=OBJ1                                                      ADS49020
      TEMP(5)=PHI                                                       ADS49030
      TEMP(6)=DLOBJ1                                                    ADS49040
      TEMP(7)=DLOBJ2                                                    ADS49050
      TEMP(8)=DX1                                                       ADS49060
      TEMP(9)=DX2                                                       ADS49070
      RETURN                                                            ADS49080
250   CONTINUE                                                          ADS49090
C     ONE-DIMENSIONAL SEARCH IS COMPLETE.                               ADS49100
      OBJA=OBJ/SCFO                                                     ADS49110
      IF (IPRNT2.EQ.2) WRITE (6,580) ITER,OBJA                          ADS49120
      IF (IPRNT2.GE.4) WRITE (6,590) ALPHA                              ADS49130
      IF (IPRNT2.GT.2) WRITE (6,550) OBJA                               ADS49140
      ICALC=1                                                           ADS49150
      IF (IPRNT2.GE.2) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS49160
      ICALC=4                                                           ADS49170
      IF (IPRNT2.GT.2) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS49180
C     ------------------------------------------------------------------ADS49190
C     CHECK CONVERGENCE.                                                ADS49200
C     ------------------------------------------------------------------ADS49210
      NGRAD=1                                                           ADS49220
      CALL ADS211 (OBJ,OBJ1,DELOBJ,DABOBJ,DEL1,DAB1,IDEL,IDAB)          ADS49230
      IF (ITER.GE.ITMAX) GO TO 300                                      ADS49240
      IF (NVC.GT.0.AND.ITER.LT.10) IDEL=0                               ADS49250

C	IF (NVD.GT.0.AND.ITER.LT.10) IDAB=0                               ADS49260 
	IF (NDV.GT.0.AND.ITER.LT.10) IDAB=0   !CHANGE FROM NVD TO NDV BY SONGSAK NOV06

      IF (IDEL.GE.ITRMOP.OR.IDAB.GE.ITRMOP) GO TO 300                   ADS49270
      IF (IDEL.GE.1.OR.IDAB.GE.1) GO TO 260                             ADS49280
      IF (ALPHA.LE.ZRO10) GO TO 260                                     ADS49290
      II=MAX0(IDEL,IDAB)+1                                              ADS49300
      IF (II.LT.ITRMOP) GO TO 40                                        ADS49310
C     REDUCE CONSTRAINT TOLERANCE.                                      ADS49320
260   CONTINUE                                                          ADS49330
      DCT=.5*(CT+CTMIN)                                                 ADS49340
      DCTL=.5*(CTL+CTLMIN)                                              ADS49350
      CT=CT-DCT                                                         ADS49360
      IF (CT.GT.-CTMIN) CT=-CTMIN                                       ADS49370
      CTL=CTL-DCTL                                                      ADS49380
      IF (CTL.GT.-CTLMIN) CTL=-CTLMIN                                   ADS49390
      IF (IPRNT2.GT.3) WRITE (6,360) CT,CTL                             ADS49400
C     IF LAST MOVE WAS GT.ZRO10, GO START NEW ITERATION.                ADS49410
      IF (ALPHA.GT.ZRO10.OR.NAVC.EQ.0) GO TO 40                         ADS49420
C     ------------------------------------------------------------------ADS49430
C     LAST ALPHA.LE.ZRO10.  DELETE CONSTRAINTS THAT ARE NO LONGER       ADS49440
C     CRITICAL AND CONTINUE OPTIMIZATION WITHOUT GETTING NEW GRADIENTS. ADS49450
C     ------------------------------------------------------------------ADS49460
      IF (NGT.GE.NCON) GO TO 290                                        ADS49470
C     CHECK FOR NEW CRITICAL CONSTRAINTS.                               ADS49480
      DO 280 I=1,NCON                                                   ADS49490
      I1=IDG(I)                                                         ADS49500
      IF (I1.LT.0) GO TO 280                                            ADS49510
      CC=CT                                                             ADS49520
      IF (I1.GT.1) CC=CTL                                               ADS49530
      IF (G(I).LT.CC) GO TO 280                                         ADS49540
C     THIS CONSTRAINT IS CRITICAL.  IS IT ALREADY IN IC?                ADS49550
      DO 270 J=1,NGT                                                    ADS49560
      IF (IC(J).EQ.I) GO TO 280                                         ADS49570
270   CONTINUE                                                          ADS49580
C     NEW CRITICAL CONSTRAINT.  MUST EVALUATE GRADIENTS.                ADS49590
      NGRAD=1                                                           ADS49600
      GO TO 40                                                          ADS49610
280   CONTINUE                                                          ADS49620
290   CONTINUE                                                          ADS49630
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,CTLMIN,NAC,NVC)  ADS49640
      NAVC=NAC+NVC                                                      ADS49650
C     WE HAVE GRADIENTS OF ALL CRITICAL CONSTRAINTS.                    ADS49660
C     CONTINUE OPTIMIZATION WITHOUT GETTING NEW GRADIENTS.              ADS49670
      NGRAD=0                                                           ADS49680
      GO TO 40                                                          ADS49690
300   CONTINUE                                                          ADS49700
C     ------------------------------------------------------------------ADS49710
C     OPTIMIZATION IS COMPLETE.                                         ADS49720
C     ------------------------------------------------------------------ADS49730
      IF (IPRNT2.LE.0) GO TO 340                                        ADS49740
C     PRINT FINAL RESULTS.                                              ADS49750
      WRITE (6,380)                                                     ADS49760
      WRITE (6,390) ITER                                                ADS49770
      OBJA=OBJ/SCFO                                                     ADS49780
      WRITE (6,550) OBJA                                                ADS49790
C     X-VECTOR.                                                         ADS49800
      ICALC=1                                                           ADS49810
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS49820
C     G-VECTOR.                                                         ADS49830
      ICALC=4                                                           ADS49840
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS49850
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS49860
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS49870
C     SIDE CONSTRAINTS.                                                 ADS49880
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS49890
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS49900
      WRITE (6,370) CT,CTL                                              ADS49910
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS49920
      NAVC=NAC+NVC                                                      ADS49930
      WRITE (6,470) NAC,NVC                                             ADS49940
      IF (NAVC.GT.0) WRITE (6,480)                                      ADS49950
      IF (NAVC.GT.0) WRITE (6,490) (ITMP(I),I=1,NAVC)                   ADS49960
      WRITE (6,500) NACS                                                ADS49970
      IF (NACS.EQ.0) GO TO 330                                          ADS49980
      N=0                                                               ADS49990
      DO 320 I=1,NDV                                                    ADS50000
      IF (ISIDC(I).EQ.0) GO TO 320                                      ADS50010
      IF (ISIDC(I).EQ.1) GO TO 310                                      ADS50020
C     LOWER BOUND IS ACTIVE.                                            ADS50030
      N=N+1                                                             ADS50040
      ITMP(N)=-I                                                        ADS50050
310   IF (ISIDC(I).LT.0) GO TO 320                                      ADS50060
C     UPPER BOUND IS ACTIVE.                                            ADS50070
      N=N+1                                                             ADS50080
      ITMP(N)=I                                                         ADS50090
320   CONTINUE                                                          ADS50100
      WRITE (6,510)                                                     ADS50110
      WRITE (6,490) (ITMP(I),I=1,NACS)                                  ADS50120
330   WRITE (6,400)                                                     ADS50130
      IF (ITER.GE.ITMAX) WRITE (6,410)                                  ADS50140
      IF (IDEL.GE.ITRMOP) WRITE (6,420) ITRMOP                          ADS50150
      IF (IDAB.GE.ITRMOP) WRITE (6,430) ITRMOP                          ADS50160
      IF (SCALS.LE.STOL) WRITE (6,440) SCALS,STOL                       ADS50170
      IF (BETA.LE.ZRO10.AND.NACT.GT.0) WRITE (6,450) BETA,ZRO10         ADS50180
340   CONTINUE                                                          ADS50190
C     PUT FINAL INFO IN TEMP FOR TRANSFER                               ADS50200
      TEMP(1)=ALPSAV                                                    ADS50210
      TEMP(2)=BETA                                                      ADS50220
      TEMP(3)=DF2                                                       ADS50230
      TEMP(4)=OBJ1                                                      ADS50240
      TEMP(5)=PHI                                                       ADS50250
      TEMP(6)=DLOBJ1                                                    ADS50260
      TEMP(7)=DLOBJ2                                                    ADS50270
      TEMP(8)=DX1                                                       ADS50280
      TEMP(9)=DX2                                                       ADS50290
      TEMP(10)=GMAX1                                                    ADS50300
      IGOTO=0                                                           ADS50310
      RETURN                                                            ADS50320
C     ------------------------------------------------------------------ADS50330
C     FORMATS.                                                          ADS50340
C     ------------------------------------------------------------------ADS50350
350   FORMAT ('1')                                                      ADS50360
360   FORMAT (/,5X,'NEW CONSTRAINT TOLERANCE',2X,'CT =',1PE13.5,2X,     ADS50370
     1'CTL =',1PE13.5)                                                  ADS50380
370   FORMAT (/,5X,'CONSTRAINT TOLERANCE, CT =',1PE13.5,2X,'CTL =',     ADS50390
     11PE13.5)                                                          ADS50400
380   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS')                     ADS50410
390   FORMAT (/,5X,'NUMBER OF ITERATIONS =',I5)                         ADS50420
400   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS50430
410   FORMAT (/,5X,'NUMBER OF ITERATIONS = ITMAX')                      ADS50440
420   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS50450
     1' CONSECUTIVE ITERATIONS')                                        ADS50460
430   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS50470
     1' CONSECUTIVE ITERATIONS')                                        ADS50480
440   FORMAT (/,5X,'MAXIMUM K-T RESIDUAL =',1PE13.5,' IS LESS THAN',    ADS50490
     11PE13.5)                                                          ADS50500
450   FORMAT (/,5X,'KUHN-TUCKER PARAMETER, BETA =',1PE13.5,             ADS50510
     1' IS LESS THAN',1PE13.5)                                          ADS50520
460   FORMAT (//,5X,75('-'),/,23X,                                      ADS50530
     1'IOPT = 4; METHOD OF FEASIBLE DIRECTIONS',/,5X,75('-'))           ADS50540
470   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS50550
     1' VIOLATED CONSTRAINTS')                                          ADS50560
480   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS50570
490   FORMAT (5X,10I5)                                                  ADS50580
500   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS50590
510   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS50600
520   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS50610
530   FORMAT (/,5X,'K-T PARAMETERS, BETA =',1PE13.5,'  MAX. RESIDUAL =',ADS50620
     11PE13.5)                                                          ADS50630
540   FORMAT (/,5X,'PROPOSED ALPHA =',1PE13.5)                          ADS50640
550   FORMAT (/,5X,11HOBJECTIVE =,1PE12.5)                              ADS50650
560   FORMAT (/,5X,'-- INITIAL DESIGN',//,5X,'OBJ =',1PE12.5)           ADS50660
570   FORMAT (//,5X,'-- BEGIN ITERATION NUMBER',I5)                     ADS50670
580   FORMAT (//,5X,'-- ITERATION',I5,5X,'OBJ =',1PE12.5)               ADS50680
590   FORMAT (/,5X,'CALCULATED ALPHA =',1PE13.5)                        ADS50690
      END                                                               ADS50700
C --- SUBROUTINE ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                  ADS50710
C     DIMENSION XOLD(NDV),S(NDV),SVLB(NDV),SVUB(NDV),X(NDV)             ADS50720
C     ------------------------------------------------------------------ADS50730
C --- ROUTINE TO UPDATE THE X-VECTOR SUBJECT TO BOUNDS SVLB AND SVUB.   ADS50740
C     ------------------------------------------------------------------ADS50750
C     BY G. N. VANDERPLAATS                          SEPT. 1981         ADS50760
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS50770
C                                                                       ADS50780
C --- INPUT                                                             ADS50790
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS50800
C     ALPHA - MOVE PARAMETER IN ONE-DIMENSIONAL SEARCH.                 ADS50810
C     XOLD - OLD X-VECTOR.  MINIMUM DIMENSION IS NDV.                   ADS50820
C     S - SEARCH DIRECTION.  MINIMUM DIMENSION IS NDV.                  ADS50830
C     SVLB - LOWER BOUNDS ON X.  MINIMUM DIMENSION IS NDV.              ADS50840
C     SVUB - UPPER BOUNDS ON X.  MINIMUM DIMENSION IS NDV.              ADS50850
C                                                                       ADS50860
C --- OUTPUT                                                            ADS50870
C     X - NEW VECTOR OF DESIGN VARIABLES.  X = XOLD + ALPHA*S.          ADS50880
C         IF X(I).LT.SVLB(I), X(I)=SVLB(I).                             ADS50890
C         IF X(I).GT.SVUB(I), X(I)=SVUB(I).                             ADS50900
C         MINIMUM DIMENSION OF X IS NDV.                                ADS50910
C                                                                       ADS50920
      SUBROUTINE ADS206 (NDV,ALPHA,XOLD,S,SVLB,SVUB,X)                  ADS50930
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION XOLD(1),S(1),SVLB(1),SVUB(1),X(1)                       ADS50940
C                                                                       ADS50950
      DO 10 I=1,NDV                                                     ADS50960
      X(I)=XOLD(I)+ALPHA*S(I)                                           ADS50970
      IF (X(I).GT.SVUB(I)) X(I)=SVUB(I)                                 ADS50980
      IF (X(I).LT.SVLB(I)) X(I)=SVLB(I)                                 ADS50990
10    CONTINUE                                                          ADS51000
      RETURN                                                            ADS51010
      END                                                               ADS51020
C --- SUBROUTINE ADS207(N1,N2,A,NRA,ZRO,GGL)                            ADS51030
C     DIMENSION A(N1,N2),GGL(N2)                                        ADS51040
C     ------------------------------------------------------------------ADS51050
C --- ROUTINE TO FIND THE MAXIMUM COMPONENT MAGNITUDE OF EACH VECTOR    ADS51060
C     CONTAINED IN ARRAY A.                                             ADS51070
C     ------------------------------------------------------------------ADS51080
C     BY G. N. VANDERPLAATS                              AUG 1983       ADS51090
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS51100
C                                                                       ADS51110
C --- INPUT                                                             ADS51120
C     N1 - NUMBER OR VARIABLES (ROWS).                                  ADS51130
C     N2 - NUMBER OF VECTORS (COLUMNS).                                 ADS51140
C     A - ARRAY CONTAINING VECTORS TO BE NORMALIZED.                    ADS51150
C         MINIMUM DIMENSION IS (N1,N2).                                 ADS51160
C     NRA - DIMENSIONED ROWS OF ARRAY A IN CALLING PROGRAM.  NRA MUST   ADS51170
C         BE AT LEAST EQUAL TO N1.                                      ADS51180
C     ZRO - TOLERANCE IDENTIFYING ZERO ON THE COMPUTER.  A SMALL        ADS51190
C         POSITIVE NUMBER IDENTIFYING MACHINE PRECISION.                ADS51200
C                                                                       ADS51210
C --- OUTPUT                                                            ADS51220
C     GGL - ARRAY CONTAINING NORMALIZATION FACTORS.                     ADS51230
C         MINIMUM DIMENSION IS N2.                                      ADS51240
C                                                                       ADS51250
C --- NOTE                                                              ADS51260
C     IF A COLUMN OF A CONTAINS ALL ZEROS, THE NORMALIZATION FACTOR     ADS51270
C         IS SET TO ZRO.                                                ADS51280
C                                                                       ADS51290
      SUBROUTINE ADS207 (N1,N2,A,NRA,ZRO,GGL)                           ADS51300
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(NRA,1),GGL(1)                                         ADS51310
C                                                                       ADS51320
C     BEGIN CALCULATION OF NORMALIZATION FACTORS.                       ADS51330
      DO 20 J=1,N2                                                      ADS51340
C     MAGNITUDE OF MAXIMUM COMPONENT.                                   ADS51350
      AMXI=ZRO                                                          ADS51360
      DO 10 I=1,N1                                                      ADS51370
      AI=ABS(A(I,J))                                                    ADS51380
      IF (AI.GT.AMXI) AMXI=AI                                           ADS51390
10    CONTINUE                                                          ADS51400
      GGL(J)=1./AMXI                                                    ADS51410
20    CONTINUE                                                          ADS51420
      RETURN                                                            ADS51430
      END                                                               ADS51440
C --- SUBROUTINE ADS209(NDV,NAVC,NACS,A,NRA,B,NRB,ISIDC)                ADS51450
C     DIMENSION A(NRA,NCOLA),B(NRB,NRB),ISIDC(NDV)                      ADS51460
C     ------------------------------------------------------------------ADS51470
C --- ROUTINE TO CALCULATE B = A-TRANS TIMES A FOR USE IN THE ROBUST    ADS51480
C     FEASIBLE DIRECTIONS METHOD.                                       ADS51490
C     ------------------------------------------------------------------ADS51500
C     BY G. N. VANDERPLAATS                              AUG. 1983      ADS51510
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS51520
C                                                                       ADS51530
C --- INPUT                                                             ADS51540
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS51550
C     NAVC - NUMBER OF ACTIVE/VIOLATED CONSTRAINTS.                     ADS51560
C     NACS - NUMBER OF ACTIVE SIDE CONSTRAINTS.                         ADS51570
C     A - ARRAY CONTAINING GRADIENTS OF CONSTRAINTS STORED IN COLUMNS.  ADS51580
C     NRA - DIMENSIONED ROWS OF A.                                      ADS51590
C     ISIDC - VECTOR IDENTIFYING ACTIVE SIDE CONSTRAINTS.               ADS51600
C                                                                       ADS51610
C --- OUTPUT                                                            ADS51620
C     B = A-TRANSPOSE TIMES A, AUGMENTED TO INCLUDE GRADIENTS OF ACTIVE ADS51630
C         SIDE CONSTRAINTS.                                             ADS51640
C                                                                       ADS51650
      SUBROUTINE ADS209 (NDV,NAVC,NACS,A,NRA,B,NRB,ISIDC)               ADS51660
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(NRA,1),B(NRB,1),ISIDC(1)                              ADS51670
C                                                                       ADS51680
      NACT=NAVC+NACS                                                    ADS51690
      IF (NAVC.EQ.0) GO TO 70                                           ADS51700
C     ------------------------------------------------------------------ADS51710
C     BUILD A-TRANS TIMES A, LOWER TRIANGLE.                            ADS51720
C     ------------------------------------------------------------------ADS51730
      DO 20 I=1,NAVC                                                    ADS51740
      DO 20 J=I,NAVC                                                    ADS51750
      B(J,I)=0.                                                         ADS51760
C     DEL-G TRANS TIMES DEL-G.                                          ADS51770
      DO 10 K=1,NDV                                                     ADS51780
10    B(J,I)=B(J,I)+A(K,I)*A(K,J)                                       ADS51790
20    CONTINUE                                                          ADS51800
C     ------------------------------------------------------------------ADS51810
C     SIDE CONSTRAINT TERMS.                                            ADS51820
C     ------------------------------------------------------------------ADS51830
      IF (NACS.EQ.0) GO TO 110                                          ADS51840
      J=NAVC                                                            ADS51850
      DO 60 II=1,NDV                                                    ADS51860
      ID=ISIDC(II)                                                      ADS51870
      IF (ID.EQ.0) GO TO 60                                             ADS51880
      IF (ID.EQ.1) GO TO 40                                             ADS51890
C     LOWER BOUND IS ACTIVE.                                            ADS51900
      J=J+1                                                             ADS51910
      DO 30 I=1,NAVC                                                    ADS51920
30    B(J,I)=-A(II,I)                                                   ADS51930
      IF (ID.LT.0) GO TO 60                                             ADS51940
40    CONTINUE                                                          ADS51950
C     UPPER BOUND IS ACTIVE.                                            ADS51960
      J=J+1                                                             ADS51970
      DO 50 I=1,NAVC                                                    ADS51980
50    B(J,I)=A(II,I)                                                    ADS51990
60    CONTINUE                                                          ADS52000
70    CONTINUE                                                          ADS52010
C     ------------------------------------------------------------------ADS52020
C     TERMS IN B UNIQUE TO SIDE CONSTRAINTS.                            ADS52030
C     ------------------------------------------------------------------ADS52040
      J=NAVC                                                            ADS52050
      DO 100 II=1,NDV                                                   ADS52060
      ID=ISIDC(II)                                                      ADS52070
      IF (ID.EQ.0) GO TO 100                                            ADS52080
C     ACTIVE SIDE CONSTRAINT(S).                                        ADS52090
      JD=IABS(ID)                                                       ADS52100
      DO 90 K=1,JD                                                      ADS52110
      J=J+1                                                             ADS52120
      DO 80 I=J,NACT                                                    ADS52130
      B(I,J)=0.                                                         ADS52140
      IF (ID.EQ.2) B(I,J)=-1.                                           ADS52150
80    CONTINUE                                                          ADS52160
90    B(J,J)=1.                                                         ADS52170
100   CONTINUE                                                          ADS52180
110   CONTINUE                                                          ADS52190
C     ------------------------------------------------------------------ADS52200
C     IMPOSE SYMMETRY.                                                  ADS52210
C     ------------------------------------------------------------------ADS52220
      DO 120 I=1,NACT                                                   ADS52230
      DO 120 J=I,NACT                                                   ADS52240
120   B(I,J)=B(J,I)                                                     ADS52250
      RETURN                                                            ADS52260
      END                                                               ADS52270
C --- SUBROUTINE ADS211(OBJ,OBJ1,DEL,DAB,DEL1,DAB1,IDEL,IDAB)           ADS52280
C     ------------------------------------------------------------------ADS52290
C --- ROUTINE TO UPDATE CONVERGENCE PARAMETERS IDEL AND IDAB TO         ADS52300
C     IDENTIFY CONVERGENCE OF THE OBJECTIVE FUNCTION, OBJ.              ADS52310
C     ------------------------------------------------------------------ADS52320
C     BY G. N. VANDERPLAATS                  SEPT. 1980                 ADS52330
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS52340
C                                                                       ADS52350
C --- INPUT                                                             ADS52360
C     OBJ - PRESENT VALUE OF THE OBJECTIVE FUNCTION.                    ADS52370
C     OBJ1 - VALUE OF THE OBJECTIVE FUNCTION ON THE PREVIOUS ITERATION. ADS52380
C     DEL - RELATIVE CONVERGENCE TOLERANCE.                             ADS52390
C         CONVERGENCE IF 1.-ABS(OBJ)/ABS(OBJ1).LE.DEL.                  ADS52400
C     DAB - ABSOLUTE CONVERGENCE CRITERION.                             ADS52410
C         CONVERGED IF ABS(OBJ)-ABS(OBJ1).LE.DAB.                       ADS52420
C     IDEL - CONVERGENCE IDENTIFIER.  IDEL=0 ON FIRST CALL.  ON         ADS52430
C         SUBSEQUENT CALLS, IDEL=THE VALUE RETURNED FROM PREVIOUS CALL. ADS52440
C     IDAB - CONVERGENCE IDENTIFIER.  IDAB=0 ON FIRST CALL.  ON         ADS52450
C         SUBSEQUENT CALLS, IDAB=THE VALUE RETURNED FROM PREVIOUS CALL. ADS52460
C                                                                       ADS52470
C --- OUTPUT                                                            ADS52480
C     DEL1 - RELATIVE CHANGE IN THE OBJECTIVE FUNCTION.                 ADS52490
C     DAB1 - ABSOLUTE CHANGE IN THE OBJECTIVE FUNCTION.                 ADS52500
C     IDEL - CONVERGENCE IDENTIFIER.                                    ADS52510
C         IDEL = 0 IF DEL CONVERGENCE CRITERION IS NOT SATISFIED.       ADS52520
C         IDEL = IDEL + 1 IF DEL CONVERGENCE CRITERION IS SATISFIED.    ADS52530
C     IDAB - CONVERGENCE IDENTIFIER.                                    ADS52540
C         IDAB = 0 IS DAB CONVERGENCE CRITERION IS NOT SATISFIED.       ADS52550
C         IDAB = IDAB + 1 IF DAB CONVERGENCE CRITERION IS SATISFIED.    ADS52560
C                                                                       ADS52570
      SUBROUTINE ADS211 (OBJ,OBJ1,DEL,DAB,DEL1,DAB1,IDEL,IDAB)          ADS52580
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	OBKK = ABS(OBJ1)        !CHANGED BY SONGSAK NOV06
	TNORM = 1.0E-10         !CHANGED BY SONGSAK NOV06
      DENOM=DMAX1(OBKK,TNORM)                                           ADS52590


      DEL1=1.-ABS(OBJ)/DENOM                                            ADS52600
      IDEL=IDEL+1                                                       ADS52610
      IF (ABS(DEL1).GT.DEL) IDEL=0                                      ADS52620
      DAB1=ABS(OBJ1)-ABS(OBJ)                                           ADS52630
      IDAB=IDAB+1                                                       ADS52640
      IF (ABS(DAB1).GT.DAB) IDAB=0                                      ADS52650
      RETURN                                                            ADS52660
      END                                                               ADS52670
C --- SUBROUTINE ADS213(OBJ,DFDALP,DLOBJ1,DLOBJ2,ALPHA)                 ADS52680
C     ------------------------------------------------------------------ADS52690
C --- ROUTINE TO ESTIMATE AN INITIAL STEP, ALPHA, TO CHANGE OBJ BY THE  ADS52700
C     FRACTION DLOBJ1 OR THE MAGNITUDE DLOBJ2, WHICHEVER IS LESS.       ADS52710
C     ------------------------------------------------------------------ADS52720
C     BY G. N. VANDERPLAATS                         NOV. 1980           ADS52730
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS52740
C                                                                       ADS52750
C --- INPUT                                                             ADS52760
C     OBJ - PRESENT VALUE OF THE OBJECTIVE FUNCTION.                    ADS52770
C     DFDALP - THE SLOPE OF THE OBJ VS. ALPHA CURVE AT ALPHA=0.         ADS52780
C         DFDALP = THE VECTOR PRODUCT OF DF AND S.                      ADS52790
C     DLOBJ1 - RELATIVE CHANGE DESIRED IN THE OBJECTIVE.                ADS52800
C     DLOBJ2 - MINIMUM ABSOLUTE CHANGE DESIRED IN THE OBJECTIVE.        ADS52810
C                                                                       ADS52820
C --- OUTPUT                                                            ADS52830
C     ALPHA - INITIAL STEP FOR ONE-DIMENSIONAL SEARCH TO REDUCE THE     ADS52840
C         OBJECTIVE BY THE FRACTION DLOBJ1 OR THE ABSOLUTE AMOUNT       ADS52850
C         DLOBJ2, WHICHEVER IS LESS.                                    ADS52860
C                                                                       ADS52870
      SUBROUTINE ADS213 (OBJ,DFDALP,DLOBJ1,DLOBJ2,ALPHA)                ADS52880
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DENOM=DFDALP                                                      ADS52890
      IF (ABS(DFDALP).LT.1.0E-20) DENOM=1.0E-20                         ADS52900
      ANUM=(DLOBJ1*OBJ)                                                 ADS52910
      IF (ANUM.LT.0.01) ANUM=0.01                                       ADS52920
      ALP1=ANUM/DENOM                                                   ADS52930
      ALPHA=ABS(ALP1)                                                   ADS52940
      ALP1=DLOBJ2/DENOM                                                 ADS52950
      ALP1=ABS(ALP1)                                                    ADS52960
      IF (ALP1.LT.ALPHA) ALPHA=ALP1                                     ADS52970
      RETURN                                                            ADS52980
      END                                                               ADS52990
C --- SUBROUTINE ADS214(NDV,X,S,DX1,DX2,ALPHA)                          ADS53000
C     DIMENSION X(NDV),S(NDV)                                           ADS53010
C     ------------------------------------------------------------------ADS53020
C --- ROUTINE TO ESTIMATE AN INITIAL STEP, ALPHA, TO CHANGE SOME        ADS53030
C     COMPONENT OF X BY THE GREATEST OF DX1*X(I) AND DX2.               ADS53040
C     ------------------------------------------------------------------ADS53050
C     BY G. N. VANDERPLAATS                            NOV. 1980        ADS53060
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS53070
C                                                                       ADS53080
C --- INPUT                                                             ADS53090
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS53100
C     X - PRESENT DESIGN VARIABLES.  MINIMUM DIMENSION OF X IS NDV.     ADS53110
C     S - SEARCH DIRECTION.  MINIMUM DIMENSION OF S IS NDV.             ADS53120
C     DX1 - RELATIVE CHANGE DESIRED IN SOME COMPONENT OF X.             ADS53130
C     DX2 - MINIMUM ABSOLUTE CHANGE DESIRED IN SOME COMPONENT OF X.     ADS53140
C                                                                       ADS53150
C --- OUTPUT                                                            ADS53160
C     ALPHA - INITIAL PROPOSED STEP IN THE ONE-DIMENSIONAL SEARCH.      ADS53170
C                                                                       ADS53180
C     NOTE - FOR EACH X(I), THE MAXIMUM OF X(I)*DX1/S(I) AND DX2/S(I)   ADS53190
C         IN MAGNITUDE IS FOUND.  FROM AMONG THESE NDV VALUES FOR ALPHA,ADS53200
C         THE SMALLEST IS CHOSEN AS THE PROPOSED MOVE.                  ADS53210
C                                                                       ADS53220
      SUBROUTINE ADS214 (NDV,X,S,DX1,DX2,ALPHA)                         ADS53230
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),S(1)                                               ADS53240
      ALPHA=1.0E+20                                                     ADS53250
      DO 10 I=1,NDV                                                     ADS53260
      SI=ABS(S(I))                                                      ADS53270
C     PREVENT DIVIDE BY ZERO.                                           ADS53280
      IF (SI.LT.1.0E-20) SI=1.0E-20                                     ADS53290
      DX=X(I)                                                           ADS53300
C     MOVE TO CHANGE X(I) BY THE FRACTION DX1.                          ADS53310
      ALP=DX*DX1/SI                                                     ADS53320
      ALP=ABS(ALP)                                                      ADS53330
C     MOVE TO CHANGE X(I) BY THE AMOUNT DX2.                            ADS53340
      ALP2=DX2/SI                                                       ADS53350
      ALP2=ABS(ALP2)                                                    ADS53360
C     PICK MAXIMUM OF ALP AND ALP2.                                     ADS53370
      IF (ALP2.GT.ALP) ALP=ALP2                                         ADS53380
C     PICK MINIMUM OF ALPHA AND ALP.                                    ADS53390
      IF (ALP.LT.ALPHA) ALPHA=ALP                                       ADS53400
10    CONTINUE                                                          ADS53410
      RETURN                                                            ADS53420
      END                                                               ADS53430
C --- SUBROUTINE ADS216 (JGOTO,OBJ,NDV,X,SVLB,SVUB,NCON,NCG,IC,G,GOLD,  ADS53440
C    1 FDCH,FDCHM,DF,A,NRA,XSAV,ZRO,S,NXFD,OBJS)                        ADS53450
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),IC(NCG),G(NCON),GOLD(NCON),  ADS53460
C    1 DF(NDV),A(NRA,NCOLA),XSAV(NDV),S(NDV)                            ADS53470
C     ------------------------------------------------------------------ADS53480
C --- ROUTINE TO CALCULATE GRADIENT OF FUNCTION  SUBJECT TO SIDE        ADS53490
C     CONSTRAINTS AND A SPECIFIED SET OF GENERAL CONSTRAINTS BY FIRST   ADS53500
C     FORWARD FINITE DIFFERENCE                                         ADS53510
C     ------------------------------------------------------------------ADS53520
C     BY G. N. VANDERPLAATS                      NOV. 1980              ADS53530
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS53540
C --- INPUT                                                             ADS53550
C     JGOTO - MUST BE INPUT AS ZERO ON FIRST CALL TO ADS216.  ON        ADS53560
C         SUBSEQUENT CALLS, JGOTO SHOULD BE WHAT WAS RETURNED FROM      ADS53570
C         PREVIOUS CALL.                                                ADS53580
C     OBJ - VALUE OF THE OBJECTIVE FUNCTION ASSOCIATED WITH X.          ADS53590
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS53600
C     X - VECTOR OF DESIGN VARIABLES.  MINIMUM DIMENSION IS NDV.        ADS53610
C     SVLB - LOWER BOUNDS ON X.  MINIMUM DIMENSION IS NDV.              ADS53620
C     SVUB - UPPER BOUNDS ON X.  MINIMUM DIMENSION IS NDV.              ADS53630
C     NCON - NUMBER OF CONSTRAINT FUNCTIONS STORED IN G AND GOLD.       ADS53640
C     NCG - NUMBER OF CONSTRAINTS FOR WHICH GRADIENTS ARE TO BE         ADS53650
C         EVALUATED.  NCG MAY BE ZERO.                                  ADS53660
C     IC - LOCATIONS OF CONSTRAINTS IN G FOR WHICH GRADIENTS ARE TO     ADS53670
C     BE EVALUATED.  MINIMUM DIMENSION IS NCG.                          ADS53680
C     G - ARRAY OF CONSTRAINT VALUES ASSOCIATED WITH THE CURRENT        ADS53690
C         X-VECTOR.  MINIMUM DIMENSION IS NCON.                         ADS53700
C     GOLD - ARRAY OF CONSTRAINT VALUES ASSOCIATED WITH THE NOMINAL     ADS53710
C         X-VECTOR.  MINIMUM DIMENSION IS NCON.                         ADS53720
C     NRA - NUMBER OF ROWS IN ARRAY A.                                  ADS53730
C     FDCH - FINITE DIFFERENCE STEP AS A FRACTION OF X(I).              ADS53740
C     FDCHM - MINIMUM ABSOLUTE STEP FOR FINITE DIFFERENCE.              ADS53750
C     XSAV - OLD VALUES OF DESIGN VARIABLES.                            ADS53760
C     S    - WORK ARRAY TO SAVE DX.                                     ADS53770
C --- OUTPUT                                                            ADS53780
C     JGOTO - CALCULATE CONTROL.  IF JGOTO IS GREATER THAN ZERO,        ADS53790
C     EVALUATE OBJ FOR THE CURRENT VALUE OF X AND CALL ADS216 AGAIN.    ADS53800
C         IF JGOTO = 0, CALCULATION OF FINITE DIFFERENCE GRADIENTS IS   ADS53810
C         COMPLETE.                                                     ADS53820
C     X - NEW PERTURBED X VECTOR.  IF JGOTO = 0, X IS THE ORIGINAL      ADS53830
C         X-VECTOR.                                                     ADS53840
C     DF - VECTOR CONTAINING FIRST PARTIAL DERIVATIVES.  DF IS COMPLETE ADS53850
C         WHEN JGOTO IS RETURNED AS ZERO.                               ADS53860
C     A - ARRAY CONTAINING GRADIENTS OF THE SPECIFIED CONSTRAINTS.      ADS53870
C         THE REQUIRED DIMENSIONS ARE (NRA,NCG).                        ADS53880
C     NXFD - DESIGN VARIABLE ASSOCIATED WITH CURRENT F.D. STEP.         ADS53890
C      OBJS - NOMINAL OBJECTIVE.                                        ADS53900
C                                                                       ADS53910
      SUBROUTINE ADS216 (JGOTO,OBJ,NDV,X,SVLB,SVUB,NCON,NCG,IC,G,GOLD,FDADS53920
     1CH,FDCHM,DF,A,NRA,XSAV,ZRO,S,NXFD,OBJS)                           ADS53930
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),DF(1),IC(1),G(1),GOLD(1),A(NRA,1),SVLB(1),SVUB(1),XADS53940
     1SAV(1),S(1)                                                       ADS53950
      IF (JGOTO.GT.0) GO TO 10                                          ADS53960
C     FIRST CALL.                                                       ADS53970
      OBJS=OBJ                                                          ADS53980
      NXFD=1                                                            ADS53990
C     FINITE DIFFERENCE STEP.                                           ADS54000
      DX=FDCH*ABS(XSAV(1))                                              ADS54010
      IF (DX.LT.FDCHM) DX=FDCHM                                         ADS54020
      X1=XSAV(1)+DX                                                     ADS54030
      IF (X1.GT.SVUB(1)) DX=-DX                                         ADS54040
      S(NXFD)=DX                                                        ADS54050
      X(1)=XSAV(1)+DX                                                   ADS54060
C     EVALUATE OBJECTIVE AND CONSTRAINTS.                               ADS54070
      JGOTO=1                                                           ADS54080
      RETURN                                                            ADS54090
10    CONTINUE                                                          ADS54100
C     GRADIENT OF OBJECTIVE AND CONSTRAINTS.                            ADS54110
      DF(NXFD)=(OBJ-OBJS)/S(NXFD)                                       ADS54120
      IF (ABS(DF(NXFD)).LT.ZRO) DF(NXFD)=0.                             ADS54130
      IF (NCG.EQ.0) GO TO 30                                            ADS54140
      DO 20 J=1,NCG                                                     ADS54150
      JX=IC(J)                                                          ADS54160
      A(NXFD,J)=G(JX)-GOLD(JX)                                          ADS54170
20    CONTINUE                                                          ADS54180
30    X(NXFD)=XSAV(NXFD)                                                ADS54190
      NXFD=NXFD+1                                                       ADS54200
      JGOTO=0                                                           ADS54210
C     TERMINATE IF SOLUTION IS COMPLETE.                                ADS54220
      OBJ=OBJS                                                          ADS54230
      IF (NXFD.GT.NDV) GO TO 40                                         ADS54240
C     CONTINUE WITH FINITE DIFFERENCE GRADIENT COMPUTATIONS.            ADS54250
      DX=FDCH*ABS(XSAV(NXFD))                                           ADS54260
      IF (DX.LT.FDCHM) DX=FDCHM                                         ADS54270
      X1=XSAV(NXFD)+DX                                                  ADS54280
      IF (X1.GT.SVUB(NXFD)) DX=-DX                                      ADS54290
      S(NXFD)=DX                                                        ADS54300
      X(NXFD)=XSAV(NXFD)+DX                                             ADS54310
C     EVALUATE OBJECTIVE FUNCTION.                                      ADS54320
      OBJ=OBJS                                                          ADS54330
      JGOTO=1                                                           ADS54340
      RETURN                                                            ADS54350
40    CONTINUE                                                          ADS54360
      NXFD=0                                                            ADS54370
      IF (NCG.EQ.0) RETURN                                              ADS54380
      DO 90 J=1,NCG                                                     ADS54390
      NZRO=0                                                            ADS54400
      DO 50 I=1,NDV                                                     ADS54410
      DG=ABS(A(I,J))                                                    ADS54420
      IF (DG.LE.ZRO) NZRO=NZRO+1                                        ADS54430
50    CONTINUE                                                          ADS54440
      IF (NZRO.EQ.NDV) GO TO 70                                         ADS54450
      DO 60 I=1,NDV                                                     ADS54460
      DG=ABS(A(I,J))                                                    ADS54470
      IF (DG.LE.ZRO) A(I,J)=0.                                          ADS54480
      A(I,J)=A(I,J)/S(I)                                                ADS54490
60    CONTINUE                                                          ADS54500
      GO TO 90                                                          ADS54510
70    CONTINUE                                                          ADS54520
      DO 80 I=1,NDV                                                     ADS54530
      A(I,J)=A(I,J)/S(I)                                                ADS54540
80    CONTINUE                                                          ADS54550
90    CONTINUE                                                          ADS54560
      RETURN                                                            ADS54570
      END                                                               ADS54580
C --- SUBROUTINE ADS217 (NDV,NAC,NVC,NACS,CT,CTMIN,CTLMIN,THETAZ,ZRO,   ADS54590
C    *PHI,GGL,G,DF,A,NRA,IDG,IC,ISIDC,ITMP,B,NRB,SCALF,S,THETA,BETA,    ADS54600
C    *RSDMX)                                                            ADS54610
C     DIMENSION G(NCON),DF(1),A(NDV+1,NAC+1),IDG(NCON),IC(NAC),         ADS54620
C    * S(1),THETA(NAC),C(NACT+1),ISIDC(NDV),ITMP(NACT),GGL(NACT)        ADS54630
C     NACT=NAC+NVC+NACS                                                 ADS54640
C     ------------------------------------------------------------------ADS54650
C --- ROUTINE TO SOLVE DIRECTION FINDING PROBLEM IN MODIFIED METHOD OF  ADS54660
C     FEASIBLE DIRECTIONS.                                              ADS54670
C     ------------------------------------------------------------------ADS54680
C     BY G. N. VANDERPLAATS                            JAN 1981         ADS54690
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS54700
C                                                                       ADS54710
C --- INPUT                                                             ADS54720
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS54730
C     NAC - NUMBER OF ACTIVE CONSTRAINTS.                               ADS54740
C     NVC - NUMBER OF VIOLATED CONSTRAINTS.                             ADS54750
C     NACS - NUMBER OF ACTIVE SIDE CONSTRAINTS.                         ADS54760
C     CT - CONSTRAINT TOLERANCE DEFINING WHEN A NONLINEAR INEQUALITY    ADS54770
C         CONSTRAINT IS ACTIVE.  CT IS A SMALL NEGATIVE NUMBER.         ADS54780
C     CTMIN - CONSTRAINT TOLERANCE DEFINING WHEN A NONLINEAR INEQUALITY ADS54790
C         CONSTRAINT IS VIOLATED.  CTMIN IS A SMALL POSITIVE NUMBER.    ADS54800
C     CTLMIN - CONSTRAINT TOLERANCE DEFINING WHEN A LINEAR INEQUALITY   ADS54810
C         CONSTRAINT IS VIOLATED.  CTLMIN IS A SMALL POSITIVE NUMBER.   ADS54820
C     THETAZ - NOMINAL PUSH-OFF FACTOR.  THE PUSH-OFF FACTOR WHEN       ADS54830
C         CONSTRAINT G(I) IS PRECISELY ZERO.  THE ACTUAL PUSH-OFF FACTORADS54840
C         IS A QUADRATIC FUNCTION OF G(I).                              ADS54850
C     ZRO - TOLERANCE TO IDENTIFY ZERO VALUE.  ESTIMATE OF REASONALBE   ADS54860
C         PRECISION FOR A GIVEN COMPUTER.                               ADS54870
C     PHI - PARTICIPATION COEFFICIENT FOR MODIFIED METHOD OF FEASIBLE   ADS54880
C         DIRECTIONS.                                                   ADS54890
C     GGL - NORMALIZATION FACTORS FOR CONSTRAINT GRADIENTS.             ADS54900
C     G - VECTOR OF CONSTRAINT VALUES.  MINIMUM DIMENSION IS NCON.      ADS54910
C     DF - GRADIENT OF THE OBJECTIVE FUNCTION.  MINIMUM DIMENSION       ADS54920
C         IS NDV.                                                       ADS54930
C     A - ARRAY CONTAINING GRADIENTS OF ACTIVE AND VIOLATED             ADS54940
C         CONSTRAINTS.  THE CONSTRAINT GRADIENTS ARE CONTAINED IN THE   ADS54950
C         FIRST NAC COLUMNS OF A.  MINIMUM DIMENSONS ARE (NDV+1,NACG)   ADS54960
C     NRA - DIMENSIONED NUMBER OF ROWS IN ARRAY A IN CALLING PROGRAM.   ADS54970
C     IDG - CONSTRAINT TYPE IDENTIFIER.                                 ADS54980
C         IDG(I) = 0, NONLINEAR INEQUALITY CONSTRAINT.                  ADS54990
C         IDG(I) = 1, LINEAR INEQUALITY CONSTRAINT.                     ADS55000
C     IC - ACTIVE AND VIOLATED CONSTRAINT IDENTIFIER.  IC(I) IS THE     ADS55010
C         LOCATION IN VECTOR G OF THE ITH ACTIVE OR VIOLATED CONSTRAINT.ADS55020
C     ISIDC - IDENTIFIES NACS ACTIVE SIDE CONSTRAINTS.                  ADS55030
C         MINIMUM DIMENSION IS NDV.                                     ADS55040
C     ITMP - WORK ARRAY.  MINIMUM DIMENSION IS NACT.                    ADS55050
C     B - WORK ARRAY.  MINIMUM DIMENSION IS (NRB,NRB).                  ADS55060
C     NRB - DIMENSIONED ROWS OF B.  MINIMUM VALUE = NACT+1.             ADS55070
C     SCALF - SCALING FACTOR FOR GRADIENT OF OBJECTIVE.                 ADS55080
C                                                                       ADS55090
C --- OUTPUT                                                            ADS55100
C     S - SEARCH DIRECTION.  MINIMUM DIMENSION IS NDV.                  ADS55110
C     THETA - VECTOR CONTAINING VALUES OF THE PUSH-OFF FACTORS.         ADS55120
C         MINIMUM DIMENSION IS NAC.                                     ADS55130
C     C - ESTIMATE OF THE LAGRANGE MULTIPLIERS.                         ADS55140
C     BETA - KUHN-TUCKER PARAMETER.                                     ADS55150
C     RSDMX - MAXIMUM RESIDUAL FROM THE KUHN-TUCKER CONDITIONS.         ADS55160
C                                                                       ADS55170
C --- EXTERNALS                                                         ADS55180
C     ADS201,ADS501,ADS507,ADS218,ADS220,ADS221,ADS222                  ADS55190
C                                                                       ADS55200
C --- NOTE                                                              ADS55210
C     NORM OF S VECTOR USED HERE IS S-TRANSPOSE TIMES S.LE.1.           ADS55220
C     IF NVC = 0 FIND DIRECTION BY ZOUTENDIJK'S METHOD.  OTHERWISE      ADS55230
C     FIND MODIFIED DIRECTION.                                          ADS55240
C                                                                       ADS55250
C --- NOTE                                                              ADS55260
C     IF THIS ROUTINE IS USED TO CALCULATE THE LAGRANGE MULTIPLIERS     ADS55270
C     ONLY, THEN THETAZ=0 SHOULD BE INPUT.                              ADS55280
C                                                                       ADS55290
      SUBROUTINE ADS217 (NDV,NAC,NVC,NACS,CT,CTMIN,CTLMIN,THETAZ,ZRO,PHIADS55300
     1,GGL,G,DF,A,NRA,IDG,IC,ISIDC,ITMP,B,NRB,SCALF,S,THETA,C,BETA,RSDMXADS55310
     2)                                                                 ADS55320
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),DF(1),A(NRA,1),IDG(1),IC(1),S(1),C(1),ISIDC(1),THETADS55330
     1A(1),ITMP(1),B(NRB,1),GGL(1)                                      ADS55340
C                                                                       ADS55350
C     NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS.                        ADS55360
      NAVC=NAC+NVC                                                      ADS55370
      NAVC1=NAVC+1                                                      ADS55380
C     TOTAL NUMBER OF ACTIVE/VIOLATED CONSTRAINTS.                      ADS55390
      NACT=NAVC+NACS                                                    ADS55400
C     ------------------------------------------------------------------ADS55410
C     CALCULATE USABLE-FEASIBLE DIRECTION.                              ADS55420
C     ------------------------------------------------------------------ADS55430
C     CALCULATE PUSH-OFF FACTORS, THETA.                                ADS55440
      IF (NAVC.GT.0) CALL ADS221 (NAVC,CT,CTMIN,CTLMIN,THETAZ,G,IC,IDG,TADS55450
     1HETA)                                                             ADS55460
C     ------------------------------------------------------------------ADS55470
C     BUILD COEFICIENT MATRIX FOR SPECIAL SUB-PROBLEM.                  ADS55480
C     ------------------------------------------------------------------ADS55490
C     BUILD A-TRANS TIMES DF TERMS OF C AND DF*DF TERM OF C.            ADS55500
      NACT1=NACT+1                                                      ADS55510
      DO 10 I=1,NACT1                                                   ADS55520
10    C(I)=0.                                                           ADS55530
      DO 60 J=1,NDV                                                     ADS55540
      D=DF(J)*SCALF                                                     ADS55550
      C(NACT1)=C(NACT1)+D*D                                             ADS55560
      IF (NAVC.EQ.0) GO TO 30                                           ADS55570
C     GENERAL CONSTRAINTS.                                              ADS55580
      DO 20 I=1,NAVC                                                    ADS55590
20    C(I)=C(I)+D*A(J,I)*GGL(I)                                         ADS55600
30    IF (NAVC.GE.NACT) GO TO 60                                        ADS55610
C     SIDE CONSTRAINTS.                                                 ADS55620
      I=NAVC                                                            ADS55630
      DO 50 II=1,NDV                                                    ADS55640
      ID=ISIDC(II)                                                      ADS55650
      IF (ID.EQ.0) GO TO 50                                             ADS55660
      IF (ID.EQ.1) GO TO 40                                             ADS55670
C     LOWER BOUND.                                                      ADS55680
      I=I+1                                                             ADS55690
      IF (II.EQ.J) C(I)=C(I)-D                                          ADS55700
      IF (ID.LT.0) GO TO 50                                             ADS55710
40    CONTINUE                                                          ADS55720
C     UPPER BOUND.                                                      ADS55730
      I=I+1                                                             ADS55740
      IF (II.EQ.J) C(I)=C(I)+D                                          ADS55750
50    CONTINUE                                                          ADS55760
60    CONTINUE                                                          ADS55770
      IF (NAVC.EQ.0) GO TO 130                                          ADS55780
C     BUILD A-TRANS TIMES A, LOWER TRIANGLE.                            ADS55790
      DO 80 I=1,NAVC                                                    ADS55800
      GA=GGL(I)                                                         ADS55810
      DO 80 J=I,NAVC                                                    ADS55820
      B(J,I)=0.                                                         ADS55830
C     DEL-G TRANS TIMES DEL-G.                                          ADS55840
      GG=GA*GGL(J)                                                      ADS55850
      DO 70 K=1,NDV                                                     ADS55860
70    B(J,I)=B(J,I)-A(K,I)*A(K,J)*GG                                    ADS55870
80    CONTINUE                                                          ADS55880
C     SIDE CONSTRAINT TERMS.                                            ADS55890
      IF (NACS.EQ.0) GO TO 130                                          ADS55900
      J=NAVC                                                            ADS55910
      DO 120 II=1,NDV                                                   ADS55920
      ID=ISIDC(II)                                                      ADS55930
      IF (ID.EQ.0) GO TO 120                                            ADS55940
      IF (ID.EQ.1) GO TO 100                                            ADS55950
C     LOWER BOUND IS ACTIVE.                                            ADS55960
      J=J+1                                                             ADS55970
      DO 90 I=1,NAVC                                                    ADS55980
90    B(J,I)=A(II,I)*GGL(I)                                             ADS55990
      IF (ID.LT.0) GO TO 120                                            ADS56000
100   CONTINUE                                                          ADS56010
C     UPPER BOUND IS ACTIVE.                                            ADS56020
      J=J+1                                                             ADS56030
      DO 110 I=1,NAVC                                                   ADS56040
110   B(J,I)=-A(II,I)*GGL(I)                                            ADS56050
120   CONTINUE                                                          ADS56060
130   CONTINUE                                                          ADS56070
      IF (NACS.EQ.0) GO TO 170                                          ADS56080
C     TERMS IN B FOR SIDE CONSTRAINTS.                                  ADS56090
      J=NAVC                                                            ADS56100
      DO 160 II=1,NDV                                                   ADS56110
      ID=ISIDC(II)                                                      ADS56120
      IF (ID.EQ.0) GO TO 160                                            ADS56130
C     ACTIVE SIDE CONSTRAINT(S).                                        ADS56140
      JD=IABS(ID)                                                       ADS56150
      DO 150 K=1,JD                                                     ADS56160
      J=J+1                                                             ADS56170
      DO 140 I=J,NACT                                                   ADS56180
      B(I,J)=0.                                                         ADS56190
      IF (ID.EQ.2) B(I,J)=1.                                            ADS56200
140   CONTINUE                                                          ADS56210
150   B(J,J)=-1.                                                        ADS56220
160   CONTINUE                                                          ADS56230
170   CONTINUE                                                          ADS56240
      IF (NAVC.EQ.0) GO TO 190                                          ADS56250
C     SUBTRACT THETA TIMES THETA-TRANS FROM A.                          ADS56260
      DO 180 I=1,NAVC                                                   ADS56270
      DO 180 J=1,I                                                      ADS56280
180   B(I,J)=B(I,J)-THETA(I)*THETA(J)                                   ADS56290
190   CONTINUE                                                          ADS56300
C     ARE CONSTRAINTS VIOLATED.                                         ADS56310
      IF (NVC.GT.0) GO TO 210                                           ADS56320
C     ------------------------------------------------------------------ADS56330
C     NO CONSTRAINTS ARE VIOLATED. BUILD NACT+1 ROW OF B MATRIX AND     ADS56340
C     FINAL C-VECTOR.                                                   ADS56350
C     ------------------------------------------------------------------ADS56360
      NDB=NACT+1                                                        ADS56370
      SAVC=C(NDB)                                                       ADS56380
      DO 200 J=1,NDB                                                    ADS56390
      B(NDB,J)=-C(J)                                                    ADS56400
      C(J)=0.                                                           ADS56410
      IF (J.GT.NAVC) GO TO 200                                          ADS56420
      B(NDB,J)=B(NDB,J)-THETA(J)                                        ADS56430
      C(J)=-THETA(J)                                                    ADS56440
200   CONTINUE                                                          ADS56450
      B(NDB,NDB)=-SAVC-1.                                               ADS56460
      C(NDB)=-1.                                                        ADS56470
      GO TO 240                                                         ADS56480
210   CONTINUE                                                          ADS56490
C     ------------------------------------------------------------------ADS56500
C     MODIFIED METHOD OF FEASIBLE DIRECTIONS.                           ADS56510
C     ------------------------------------------------------------------ADS56520
C     B IS COMPLETE.  BUILD FINAL C-VECTOR.                             ADS56530
      NDB=NACT                                                          ADS56540
      IF (NAVC.EQ.0) GO TO 230                                          ADS56550
      DO 220 I=1,NAVC                                                   ADS56560
      C(I)=C(I)-PHI*THETA(I)                                            ADS56570
220   CONTINUE                                                          ADS56580
230   CONTINUE                                                          ADS56590
240   CONTINUE                                                          ADS56600
C     IMPOSE SYMMETRY.                                                  ADS56610
      DO 250 I=1,NDB                                                    ADS56620
      DO 250 J=I,NDB                                                    ADS56630
250   B(I,J)=B(J,I)                                                     ADS56640
C     ------------------------------------------------------------------ADS56650
C     SOLVE SPECIAL L. P. PROBLEM                                       ADS56660
C     ------------------------------------------------------------------ADS56670
      CALL ADS218 (NDB,B,NRB,C,ITMP,ZRO)                                ADS56680
C     ------------------------------------------------------------------ADS56690
C     CALCULATE RESULTING DIRECTION VECTOR, S.                          ADS56700
C     ------------------------------------------------------------------ADS56710
C     USABLE-FEASIBLE DIRECTION                                         ADS56720
C     MULTIPLIER ON DF.                                                 ADS56730
      AMULT=C(NDB)                                                      ADS56740
      IF (NVC.GT.0) AMULT=1.                                            ADS56750
      AMULT=AMULT*SCALF                                                 ADS56760
      DO 270 I=1,NDV                                                    ADS56770
      S(I)=-AMULT*DF(I)                                                 ADS56780
      IF (NAVC.EQ.0) GO TO 270                                          ADS56790
C     GENERAL CONSTRAINT CASE.                                          ADS56800
      DO 260 J=1,NAVC                                                   ADS56810
260   S(I)=S(I)-A(I,J)*C(J)*GGL(J)                                      ADS56820
270   CONTINUE                                                          ADS56830
      IF (NACS.EQ.0) GO TO 300                                          ADS56840
C     SIDE CONSTRAINT.                                                  ADS56850
      J=NAVC                                                            ADS56860
      DO 290 II=1,NDV                                                   ADS56870
      ID=ISIDC(II)                                                      ADS56880
      IF (ID.EQ.0) GO TO 290                                            ADS56890
C     ACTIVE SIDE CONSTRAINT(S).                                        ADS56900
      IF (ID.EQ.1) GO TO 280                                            ADS56910
C     LOWER BOUND IS ACTIVE.                                            ADS56920
      J=J+1                                                             ADS56930
      S(II)=S(II)+C(J)                                                  ADS56940
      IF (ID.LT.0) GO TO 290                                            ADS56950
280   CONTINUE                                                          ADS56960
C     UPPER BOUND IS ACTIVE.                                            ADS56970
      J=J+1                                                             ADS56980
      S(II)=S(II)-C(J)                                                  ADS56990
290   CONTINUE                                                          ADS57000
300   CONTINUE                                                          ADS57010
C     ------------------------------------------------------------------ADS57020
C     KUHN-TUCKER PARAMETER, BETA.                                      ADS57030
C     ------------------------------------------------------------------ADS57040
C     IF CONSTRAINTS ARE VIOLATED, SET BETA = 1.                        ADS57050
      BETA=1.                                                           ADS57060
      IF (NVC.GT.0) GO TO 320                                           ADS57070
      BETA=1.-C(NDB)                                                    ADS57080
      IF (NAVC.EQ.0) GO TO 320                                          ADS57090
C     SUBTRACT THETA-TRANSPOSE TIMES U FROM BETA.                       ADS57100
      DO 310 I=1,NAVC                                                   ADS57110
310   BETA=BETA-THETA(I)*C(I)                                           ADS57120
320   CONTINUE                                                          ADS57130
C     ------------------------------------------------------------------ADS57140
C     NORMALIZE S.                                                      ADS57150
C     ------------------------------------------------------------------ADS57160
      CALL ADS507 (NDV,S,ZRO,RSDMX,IER)                                 ADS57170
C     ------------------------------------------------------------------ADS57180
C     SET NEAR ZERO COMPONENTS OF S TO PRECISE ZERO.                    ADS57190
C     ------------------------------------------------------------------ADS57200
      DO 330 I=1,NDV                                                    ADS57210
      IF (ABS(S(I)).LE.ZRO) S(I)=0.                                     ADS57220
330   CONTINUE                                                          ADS57230
      IF (NACS.EQ.0) GO TO 350                                          ADS57240
C     ------------------------------------------------------------------ADS57250
C     CHECK TO AVOID SIDE CONSTRAINT VIOLATIONS.                        ADS57260
C     ------------------------------------------------------------------ADS57270
      DO 340 I=1,NDV                                                    ADS57280
      I1=ISIDC(I)                                                       ADS57290
      IF (I1.GT.0.AND.S(I).GT.0.) S(I)=0.                               ADS57300
      IF (I1.LT.0.AND.S(I).LT.0.) S(I)=0.                               ADS57310
340   CONTINUE                                                          ADS57320
350   CONTINUE                                                          ADS57330
C     ------------------------------------------------------------------ADS57340
C     IF S IS NULL VECTOR, SET BETA TO ZERO.                            ADS57350
C     ------------------------------------------------------------------ADS57360
      IF (RSDMX.LE.ZRO) BETA=0.                                         ADS57370
C     ------------------------------------------------------------------ADS57380
C     UNSCALE LAGRANGE MULTIPLIERS.                                     ADS57390
C     ------------------------------------------------------------------ADS57400
      IF (NAVC.EQ.0) GO TO 370                                          ADS57410
      DO 360 I=1,NAVC                                                   ADS57420
360   C(I)=C(I)*GGL(I)/SCALF                                            ADS57430
370   IF (NACS.EQ.0) GO TO 390                                          ADS57440
      NAVC1=NAVC+1                                                      ADS57450
      DO 380 I=NAVC1,NACT                                               ADS57460
380   C(I)=C(I)/SCALF                                                   ADS57470
390   CONTINUE                                                          ADS57480
      RETURN                                                            ADS57490
      END                                                               ADS57500
C --- SUBROUTINE ADS218 (NDB,B,NRB,C,ITMP,ZRO)                          ADS57510
C     DIMENSION B(NRB,NRB),C(NRB),ITMP(NRB)                             ADS57520
C     ------------------------------------------------------------------ADS57530
C --- ROUTINE TO SOLVE SPECIAL LINEAR PROBLEM FOR IMPOSING S-TRANSPOSE  ADS57540
C     TIMES S.LE.1 BOUNDS IN THE MODIFIED METHOD OF FEASIBLE DIRECTIONS.ADS57550
C     ------------------------------------------------------------------ADS57560
C     BY G. N. VANDERPLAATS                             JAN 1981        ADS57570
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS57580
C                                                                       ADS57590
C --- INPUT                                                             ADS57600
C     NDB - SIZE OF ARRAY STORED IN B.  NDB = NAC+1 WHERE NAC IS THE    ADS57610
C         NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS.                    ADS57620
C     B - ARRAY OF COEFFICIENTS FOR THE SPECIAL L.P. PROBLEM IN         ADS57630
C         THE METHOD OF FEASIBLE DIRECTIONS.  MINIMUM DIMENSONS         ADS57640
C         ARE (NRB,NRB).                                                ADS57650
C     NRB - DIMENSIONED SIZE OF ARRAY B IN THE CALLING PROGRAM.         ADS57660
C     C - RIGHT HAND SIDE OF SPECIAL L.P. PROBLEM.  MINIMUM DIMENSION   ADS57670
C         IS NDB.                                                       ADS57680
C     ITMP - WORK ARRAY.  MINIMUM DIMENSION IS NDB.                     ADS57690
C     ZRO - NUMERICAL ZERO FOR COMPUTATIONAL PURPOSES.                  ADS57700
C                                                                       ADS57710
C --- OUTPUT                                                            ADS57720
C     C - MODIFIED RIGHT HAND SIDE OF THE SPECIAL L.P. EQUATION.        ADS57730
C         MINIMUM DIMENSION IS NDB.                                     ADS57740
C                                                                       ADS57750
C     NOTE - IF THE ITERATION PROCESS DOES NOT CONVERGE IN 5*NDB        ADS57760
C         ITERATIONS, THAT RESULT IS USED.  NO ERROR MESSAGE IS PRINTED.ADS57770
C                                                                       ADS57780
C     REF.  'STRUCTURAL OPTIMIZATION BY METHODS OF FEASIBLE DIRECTIONS',ADS57790
C     G. N. VANDERPLAATS AND F. MOSES, JOURNAL OF COMPUTERS             ADS57800
C     AND STRUCTURES, VOL 3, PP 739-755, 1973.                          ADS57810
C     FORM OF L. P. IS BX=C WHERE 1ST NDB COMPONENTS OF X CONTAIN VECTORADS57820
C     U AND LAST NDB COMPONENTS CONTAIN VECTOR V.  CONSTRAINTS ARE      ADS57830
C     U.GE.0, V.GE.0, AND U-TRANSPOSE TIMES V = 0.                      ADS57840
C                                                                       ADS57850
      SUBROUTINE ADS218 (NDB,B,NRB,C,ITMP,ZRO)                          ADS57860
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION B(NRB,1),C(1),ITMP(1)                                   ADS57870
C     VECTOR ITMP IDENTIFIES THE SET OF BASIC VARIABLES.                ADS57880
C     CHOOSE INITIAL BASIC VARIABLES AS V, AND INITIALIZE VECTOR ITMP   ADS57890
C     INITIALIZE ITMP.                                                  ADS57900
      DO 10 I=1,NDB                                                     ADS57910
10    ITMP(I)=0                                                         ADS57920
      ITER=0                                                            ADS57930
      EPS=-ZRO                                                          ADS57940
      NMAX=5*NDB                                                        ADS57950
C     BEGIN NEW ITERATION                                               ADS57960
20    ITER=ITER+1                                                       ADS57970
C     FIND MAX. C(I)/B(I,I) FOR I=1,NDB.                                ADS57980
      CBMAX=0.                                                          ADS57990
      ICHK=0                                                            ADS58000
      DO 30 I=1,NDB                                                     ADS58010
      C1=C(I)                                                           ADS58020
      BI=B(I,I)                                                         ADS58030
      IF (BI.GT.EPS.OR.C1.GT.0.) GO TO 30                               ADS58040
      CB=C1/BI                                                          ADS58050
      IF (CB.LE.CBMAX) GO TO 30                                         ADS58060
      ICHK=I                                                            ADS58070
      CBMAX=CB                                                          ADS58080
30    CONTINUE                                                          ADS58090
      IF (CBMAX.LT.ZRO.OR.ITER.GT.NMAX) GO TO 70                        ADS58100
      IF (ICHK.EQ.0) GO TO 70                                           ADS58110
C     UPDATE VECTOR ITMP.                                               ADS58120
      JJ=ITMP(ICHK)                                                     ADS58130
      IF (JJ.EQ.0) ITMP(ICHK)=ICHK                                      ADS58140
      IF (JJ.GT.0) ITMP(ICHK)=0                                         ADS58150
C     PIVOT OF B(ICHK,ICHK)                                             ADS58160
      BB=1./B(ICHK,ICHK)                                                ADS58170
      DO 40 J=1,NDB                                                     ADS58180
40    B(ICHK,J)=BB*B(ICHK,J)                                            ADS58190
      C(ICHK)=CBMAX                                                     ADS58200
      B(ICHK,ICHK)=BB                                                   ADS58210
C     ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE        ADS58220
C     COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE.             ADS58230
      DO 60 I=1,NDB                                                     ADS58240
      IF (I.EQ.ICHK) GO TO 60                                           ADS58250
      BB1=B(I,ICHK)                                                     ADS58260
      B(I,ICHK)=0.                                                      ADS58270
      DO 50 J=1,NDB                                                     ADS58280
50    B(I,J)=B(I,J)-BB1*B(ICHK,J)                                       ADS58290
      C(I)=C(I)-BB1*CBMAX                                               ADS58300
60    CONTINUE                                                          ADS58310
      GO TO 20                                                          ADS58320
70    CONTINUE                                                          ADS58330
      NER=0                                                             ADS58340
C     STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR         ADS58350
C     TEMPORARY STORAGE                                                 ADS58360
      DO 80 I=1,NDB                                                     ADS58370
      B(I,1)=C(I)                                                       ADS58380
80    CONTINUE                                                          ADS58390
      DO 90 I=1,NDB                                                     ADS58400
      C(I)=0.                                                           ADS58410
      J=ITMP(I)                                                         ADS58420
      IF (J.GT.0) C(I)=B(J,1)                                           ADS58430
      IF (C(I).LT.0.) C(I)=0.                                           ADS58440
90    CONTINUE                                                          ADS58450
      RETURN                                                            ADS58460
      END                                                               ADS58470
C --- SUBROUTINE ADS221(NAVC,CT,CTMIN,CTLMIN,THETAZ,G,IC,IDG,THETA)     ADS58480
C     DIMENSION G(NCON),IC(NAVC),IDG(NCON),THETA(NAVC)                  ADS58490
C     ------------------------------------------------------------------ADS58500
C --- ROUTINE TO CALCULATE THE PUSH-OFF FACTORS THETA(I) FOR THE        ADS58510
C     METHOD OF FEASIBLE DIRECTIONS.  THETA(I) IS QUADRATIC IN G.       ADS58520
C     ------------------------------------------------------------------ADS58530
C     BY G. N. VANDERPLAATS                            JAN 1981         ADS58540
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS58550
C                                                                       ADS58560
C --- INPUT                                                             ADS58570
C     NAVC - NUMBER OF ACTIVE/VIOLATED CONSTRAINTS.                     ADS58580
C     CT - CONSTRAINT TOLERANCE DEFINING WHEN A NONLINEAR INEQUALITY    ADS58590
C         CONSTRAINT IS ACTIVE.  CT IS A SMALL NEGATIVE NUMBER.         ADS58600
C     CTMIN - CONSTRAINT TOLERANCE DEFINING WHEN A NONLINEAR INEQUALITY ADS58610
C         CONSTRAINT IS VIOLATED.  CTMIN IS A SMALL POSITIVE NUMBER.    ADS58620
C     CTLMIN - CONSTRAINT TOLERANCE DEFINIGNG  WHEN A LINEAR INEQUALITY ADS58630
C         CONSTRAINT IS VIOLATED.  CTLMIN IS A SMALL POSITIVE NUMBER.   ADS58640
C     THETAZ - NOMINAL PUSH-OFF FACTOR.  THE PUSH-OFF FACTOR WHEN       ADS58650
C         CONSTRAINT G(I) IS PRECISELY ZERO.  THE ACTUAL PUSH-OFF FACTORADS58660
C         IS A QUADRATIC FUNCTION OF G(I).                              ADS58670
C     G - VECTOR OF CONSTRAINT VALUES.  MINIMUM DIMENSION IS NCON.      ADS58680
C     IC - ACTIVE AND VIOLATED CONSTRAINT IDENTIFIER.  IC(I) IS THE     ADS58690
C         LOCATION IN VECTOR G OF THE ITH ACTIVE OR VIOLATED CONSTRAINT.ADS58700
C     IDG - CONSTRAINT TYPE IDENTIFIER.                                 ADS58710
C     -2 - LINEAR EQUALITY CONSTRAINT.                                  ADS58720
C     -1 - NONLINEAR EQUALITY CONSTRAINT.                               ADS58730
C     0,1 - NONLINEAR INEQUALITY CONSTRAINT.                            ADS58740
C     2 - LINEAR INEQUALITY CONSTRAINT.                                 ADS58750
C         MINIMUM DIMENSION IS NCON.                                    ADS58760
C                                                                       ADS58770
C --- OUTPUT                                                            ADS58780
C     THETA - VECTOR OF PUSH-OFF FACTORS FOR ALL CONSTRAINTS.           ADS58790
C         MINIMUM DIMENSION IS NAVC.                                    ADS58800
C                                                                       ADS58810
      SUBROUTINE ADS221 (NAVC,CT,CTMIN,CTLMIN,THETAZ,G,IC,IDG,THETA)    ADS58820
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IC(1),IDG(1),THETA(1)                              ADS58830
C                                                                       ADS58840
C     CONSTANTS.                                                        ADS58850
      THT10=10.*THETAZ                                                  ADS58860
C     CALCULATE THETAS.                                                 ADS58870
      DO 20 I=1,NAVC                                                    ADS58880
      NCI=IC(I)                                                         ADS58890
      NCJ=IDG(NCI)                                                      ADS58900
      NCJ=IABS(NCJ)                                                     ADS58910
C     ASSUME NONLINEAR CONSTRAINT FOR NOW.                              ADS58920
      GI=G(NCI)                                                         ADS58930
      CTC=CT                                                            ADS58940
      CTV=CTMIN                                                         ADS58950
C     IS CONSTRAINT ACTUALLY LINEAR.                                    ADS58960
      IF (NCJ.LT.2) GO TO 10                                            ADS58970
C     YES.                                                              ADS58980
      CTC=-CTLMIN                                                       ADS58990
      CTV=CTLMIN                                                        ADS59000
10    CONTINUE                                                          ADS59010
      GG=GI/CTC-1.                                                      ADS59020
C     THETA VALUE.                                                      ADS59030
      THT=0.                                                            ADS59040
C     NON-LINEAR PUSH-OFF.                                              ADS59050
      THT1=THETAZ*GG*GG                                                 ADS59060
      IF (NCJ.EQ.0.OR.NCJ.EQ.1) THT=THT1                                ADS59070
      IF (GI.GT.CTV) THT=THT1                                           ADS59080
C     LIMIT THETA TO 10*THETAZ.                                         ADS59090
      IF (THT.GT.THT10) THT=THT10                                       ADS59100
      THETA(I)=THT                                                      ADS59110
20    CONTINUE                                                          ADS59120
      RETURN                                                            ADS59130
      END                                                               ADS59140
C --- SUBROUTINE ADS223(NDV,NACS,CTMIN,X,SVLB,SVUB,ISIDC)               ADS59150
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),ISIDC(NDV)                   ADS59160
C     ------------------------------------------------------------------ADS59170
C --- ROUTINE TO IDENTIFY ACTIVE SIDE CONSTRAINTS                       ADS59180
C     ------------------------------------------------------------------ADS59190
C     BY G. N. VANDERPLAATS                               JAN 1981      ADS59200
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS59210
C                                                                       ADS59220
C --- INPUT                                                             ADS59230
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS59240
C     CTMIN - TOLERANCE TO IDENTIFY ACTIVE SIDE CONSTRAINT.             ADS59250
C     X - VECTOR OF DESIGN VARIABLES.  MINIMUM DIMENSION IS NDV.        ADS59260
C     SVLB - VECTOR CONTAINING LOWER BOUNDS ON DESIGN VARIABLES, X.     ADS59270
C         MINIMUM DIMENSION IS NDV.                                     ADS59280
C     SVUB - VECTOR CONTAINING UPPER BOUNDS ON DESIGN VARIABLES, X.     ADS59290
C         MINIMUM DIMENSION IS NDV.                                     ADS59300
C                                                                       ADS59310
C --- OUTPUT                                                            ADS59320
C     NACS - THE NUMBER OF ACTIVE/VIOLATED SIDE CONSTRAINTS.            ADS59330
C     ISIDC - VECTOR IDENTIFYING ACTIVE/VIOLATED CONSTRAINTS.           ADS59340
C         MINIMUM DIMENSION IS NACS.                                    ADS59350
C                                                                       ADS59360
C --- NOTE                                                              ADS59370
C     IF SVLB(I) IS ACTIVE, ISIDC(I)=-1.                                ADS59380
C     IF SVUB(I) IS ACTIVE, ISIDC(I)=1.                                 ADS59390
C     IF BOTH SVLB(I) AND SVUB(I) ARE ACTIVE, ISIDC(I)=2.               ADS59400
C                                                                       ADS59410
      SUBROUTINE ADS223 (NDV,NACS,CTMIN,X,SVLB,SVUB,ISIDC)              ADS59420
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),ISIDC(1)                           ADS59430
C                                                                       ADS59440
C     BEGIN SEARCH FOR ACTIVE CONSTRAINTS.                              ADS59450
      NACS=0                                                            ADS59460
      DO 30 I=1,NDV                                                     ADS59470
      II=0                                                              ADS59480
      JJ=0                                                              ADS59490
C     LOWER BOUND.                                                      ADS59500
      DX=X(I)-SVLB(I)                                                   ADS59510
      IF (DX.GT.CTMIN) GO TO 10                                         ADS59520
C     ACTIVE LOWER BOUND.                                               ADS59530
      NACS=NACS+1                                                       ADS59540
      II=1                                                              ADS59550
      JJ=-1                                                             ADS59560
10    CONTINUE                                                          ADS59570
C     UPPER BOUND.                                                      ADS59580
      DX=SVUB(I)-X(I)                                                   ADS59590
      IF (DX.GT.CTMIN) GO TO 20                                         ADS59600
C     ACTIVE UPPER BOUND.                                               ADS59610
      NACS=NACS+1                                                       ADS59620
      JJ=1                                                              ADS59630
      IF (II.GT.0) JJ=2                                                 ADS59640
20    CONTINUE                                                          ADS59650
      ISIDC(I)=JJ                                                       ADS59660
30    CONTINUE                                                          ADS59670
      RETURN                                                            ADS59680
      END                                                               ADS59690
C --- SUBROUTINE ADS231 (IGOTO,INFO,IONED,NDV,NCON,ITMAX,ICNDIR,NGT,NAC,ADS59700
C    1 NVC,IPRNT2,ITRMOP,OBJ,ZRO,CT,CTMIN,CTL,CTLMIN,STOL,THETAZ,       ADS59710
C    2 DELOBJ,DABOBJ,SCFO,SLOPF,ALPHA,ALPMAX,X,SVLB,SVUB,DF,S,G,IC,IDG, ADS59720
C    3 DG,THETA,C,ISIDC,XOLD,XTMP,GGL,SCFX,ITMP,SCFG,A,NRA,NCOLA,B,NRB, ADS59730
C    4 WK,IWK)                                                          ADS59740
C     ------------------------------------------------------------------ADS59750
C --- ROUTINE TO PERFORM CONSTRAINED MINIMIZATION BY THE MODIFIED       ADS59760
C     METHOD OF FEASIBLE DIRECTIONS.                                    ADS59770
C     ------------------------------------------------------------------ADS59780
C     BY G. N. VANDERPLAATS                                  AUG. 1983  ADS59790
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS59800
C --- INPUT                                                             ADS59810
C     IDG(I) - CONSTRAINT TYPE INDICATOR.                               ADS59820
C       0,1 - NONLINEAR INEQUALITY CONSTRAINT.                          ADS59830
C         2 - LINEAR INEQUALITY CONSTRAINT.                             ADS59840
C        -1 - NONLINEAR EQUALITY CONSTRAINT.                            ADS59850
C        -2 - LINEAR EQUALITY CONSTRAINT.                               ADS59860
C                                                                       ADS59870
C --- OUTPUT                                                            ADS59880
C                                                                       ADS59890
      SUBROUTINE ADS231 (IGOTO,INFO,IONED,NDV,NCON,ITMAX,ICNDIR,NGT,NAC,ADS59900
     1NVC,IPRNT2,ITRMOP,OBJ,ZRO,CT,CTMIN,CTL,CTLMIN,STOL,THETAZ,DELOBJ,DADS59910
     2ABOBJ,SCFO,SLOPF,ALPHA,ALPMAX,X,SVLB,SVUB,DF,S,G,IC,IDG,DG,THETA,CADS59920
     3,ISIDC,XOLD,XTMP,GGL,SCFX,ITMP,SCFG,A,NRA,NCOLA,B,NRB,WK,IWK,GTMP,ADS59930
     4ITEMP,TEMP,GMX1)                                                  ADS59940
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),DF(1),S(1),G(1),IC(1),IDG(1),DG(1),ADS59950
     1THETA(1),C(1),ISIDC(1),XOLD(1),XTMP(1),GGL(1),SCFX(1),ITMP(1),SCFGADS59960
     2(1),A(NRA,1),B(NRB,1),WK(1),IWK(1),GTMP(1),ITEMP(10),TEMP(10)     ADS59970
C                                                                       ADS59980
C     GET INFO FROM TEMP                                                ADS59990
      ALPSAV=TEMP(1)                                                    ADS60000
      BETA=TEMP(2)                                                      ADS60010
      DF2=TEMP(3)                                                       ADS60020
      OBJ1=TEMP(4)                                                      ADS60030
      PHI=TEMP(5)                                                       ADS60040
      DLOBJ1=TEMP(6)                                                    ADS60050
      DLOBJ2=TEMP(7)                                                    ADS60060
      DX1=TEMP(8)                                                       ADS60070
      DX2=TEMP(9)                                                       ADS60080
      GMAX1=TEMP(10)                                                    ADS60090
      ICOUNT=ITEMP(1)                                                   ADS60100
      IDAB=ITEMP(2)                                                     ADS60110
      IDEL=ITEMP(3)                                                     ADS60120
      ITER=ITEMP(4)                                                     ADS60130
      KOUNT=ITEMP(5)                                                    ADS60140
      LGOTO=ITEMP(6)                                                    ADS60150
      NACS=ITEMP(7)                                                     ADS60160
      NAVC=ITEMP(8)                                                     ADS60170
      NGRAD=ITEMP(9)                                                    ADS60180
      ISTRAT=ITEMP(10)                                                  ADS60190
      ZRO10=10.*ZRO                                                     ADS60200
      IGOTO1=IGOTO+1                                                    ADS60210
      GO TO (10,60,330), IGOTO1                                         ADS60220
10    CONTINUE                                                          ADS60230
C     ------------------------------------------------------------------ADS60240
C     INITIAL CALL.  SET UP PARAMETERS.                                 ADS60250
C     ------------------------------------------------------------------ADS60260
C     PRINT AS REQUIRED.                                                ADS60270
      IF (ISTRAT.EQ.0.AND.IPRNT2.GT.0) WRITE (6,510)                    ADS60280
      IF (IPRNT2.GT.0) WRITE (6,610)                                    ADS60290
C     INITIALIZE PARAMETERS.                                            ADS60300
      ITER=0                                                            ADS60310
      ICOUNT=0                                                          ADS60320
      PHI=5.                                                            ADS60330
      IDEL=0                                                            ADS60340
      IDAB=0                                                            ADS60350
      NGRAD=1                                                           ADS60360
C     INITIALIZE S-VECTOR AND ISIDC-VECTOR TO ZERO.                     ADS60370
      DO 20 I=1,NDV                                                     ADS60380
      ISIDC(I)=0                                                        ADS60390
20    S(I)=0.                                                           ADS60400
      IF (IPRNT2.EQ.0) GO TO 30                                         ADS60410
      OBJA=OBJ/SCFO                                                     ADS60420
      WRITE (6,710) OBJA                                                ADS60430
      ICALC=1                                                           ADS60440
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS60450
      ICALC=5                                                           ADS60460
      CALL ADS402 (ICALC,SVLB,NDV,SCFX)                                 ADS60470
      ICALC=6                                                           ADS60480
      CALL ADS402 (ICALC,SVUB,NDV,SCFX)                                 ADS60490
      ICALC=4                                                           ADS60500
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS60510
30    CONTINUE                                                          ADS60520
C     ------------------------------------------------------------------ADS60530
C     OPTIMIZATION.  BEGIN NEW ITERATION.                               ADS60540
C     ------------------------------------------------------------------ADS60550
40    CONTINUE                                                          ADS60560
      ITER=ITER+1                                                       ADS60570
      OBJ1=OBJ                                                          ADS60580
      IF (IPRNT2.GT.2) WRITE (6,720) ITER                               ADS60590
      IF (NGRAD.EQ.0) GO TO 60                                          ADS60600
C     ------------------------------------------------------------------ADS60610
C     CALCULATE GRADIENT OF OBJECTIVE AND ACTIVE AND/OR VIOLATED        ADS60620
C     CONSTRAINTS.                                                      ADS60630
C     ------------------------------------------------------------------ADS60640
      IF (ITER.GT.1) GO TO 50                                           ADS60650
C     GRADIENTS ARE ALREADY AVAILABLE.  SHIFT A AND IC AS REQUIRED.     ADS60660
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,CTLMIN,NAC,NVC)  ADS60670
      GO TO 60                                                          ADS60680
50    CONTINUE                                                          ADS60690
C     DETERMINE THE ACTIVE CONSTRAINT SET.                              ADS60700
C     GENERAL CONSTRAINTS, G.                                           ADS60710
      IF (NGT.LT.NCON) CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIADS60720
     1N,NAC,NVC)                                                        ADS60730
      IF (NGT.EQ.NCON) CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,ADS60740
     1CTLMIN,NAC,NVC)                                                   ADS60750
C     EVALUATE GRADIENTS OF OBJECTIVE AND ALL ACTIVE/VIOLATED           ADS60760
C     CONSTRAINTS.                                                      ADS60770
      IGOTO=1                                                           ADS60780
      INFO=2                                                            ADS60790
C     PUT INFO IN TEMP FOR TRANSFER                                     ADS60800
      TEMP(1)=ALPSAV                                                    ADS60810
      TEMP(2)=BETA                                                      ADS60820
      TEMP(3)=DF2                                                       ADS60830
      TEMP(4)=OBJ1                                                      ADS60840
      TEMP(5)=PHI                                                       ADS60850
      TEMP(6)=DLOBJ1                                                    ADS60860
      TEMP(7)=DLOBJ2                                                    ADS60870
      TEMP(8)=DX1                                                       ADS60880
      TEMP(9)=DX2                                                       ADS60890
      TEMP(10)=GMAX1                                                    ADS60900
      ITEMP(1)=ICOUNT                                                   ADS60910
      ITEMP(2)=IDAB                                                     ADS60920
      ITEMP(3)=IDEL                                                     ADS60930
      ITEMP(4)=ITER                                                     ADS60940
      ITEMP(5)=KOUNT                                                    ADS60950
      ITEMP(6)=LGOTO                                                    ADS60960
      ITEMP(7)=NACS                                                     ADS60970
      ITEMP(8)=NAVC                                                     ADS60980
      ITEMP(9)=NGRAD                                                    ADS60990
      ITEMP(10)=ISTRAT                                                  ADS61000
      RETURN                                                            ADS61010
60    CONTINUE                                                          ADS61020
      NAVC=NAC+NVC                                                      ADS61030
C     IDENTIFY ACTIVE SIDE CONSTRAINTS.                                 ADS61040
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS61050
      IF (IPRNT2.LT.4) GO TO 90                                         ADS61060
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS61070
      WRITE (6,620) NAC,NVC                                             ADS61080
      IF (NAVC.GT.0) WRITE (6,630)                                      ADS61090
      IF (NAVC.GT.0) WRITE (6,640) (IC(I),I=1,NAVC)                     ADS61100
      WRITE (6,650) NACS                                                ADS61110
      IF (NACS.EQ.0) GO TO 90                                           ADS61120
      N=0                                                               ADS61130
      DO 80 I=1,NDV                                                     ADS61140
      IF (ISIDC(I).EQ.0) GO TO 80                                       ADS61150
      IF (ISIDC(I).EQ.1) GO TO 70                                       ADS61160
C     LOWER BOUND IS ACTIVE.                                            ADS61170
      N=N+1                                                             ADS61180
      ITMP(N)=-I                                                        ADS61190
70    IF (ISIDC(I).LT.0) GO TO 80                                       ADS61200
C     UPPER BOUND IS ACTIVE.                                            ADS61210
      N=N+1                                                             ADS61220
      ITMP(N)=I                                                         ADS61230
80    CONTINUE                                                          ADS61240
      WRITE (6,660)                                                     ADS61250
      WRITE (6,640) (ITMP(I),I=1,NACS)                                  ADS61260
90    CONTINUE                                                          ADS61270
      IF (IPRNT2.LT.5) GO TO 110                                        ADS61280
C     WRITE GRADIENT OF OBJECTIVE.                                      ADS61290
      ICALC=3                                                           ADS61300
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS61310
      IF (NAVC.EQ.0) GO TO 110                                          ADS61320
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS61330
      ISCAL=0                                                           ADS61340
      DO 100 J=1,NAVC                                                   ADS61350
      WRITE (6,670) IC(J)                                               ADS61360
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS61370
100   CONTINUE                                                          ADS61380
110   CONTINUE                                                          ADS61390
C     ------------------------------------------------------------------ADS61400
C     SEARCH DIRECTION.                                                 ADS61410
C     ------------------------------------------------------------------ADS61420
      NACT=NAVC+NACS                                                    ADS61430
      IF (NACT.GT.0) GO TO 130                                          ADS61440
C     NO ACTIVE OR VIOLATED CONSTRAINTS.  FIND STEEPEST DESCENT OR      ADS61450
C     FLETCHER-REEVES CONJUGATE DIRECTION.                              ADS61460
      ICOUNT=ICOUNT+1                                                   ADS61470
      IF (IDEL.EQ.1.OR.IDAB.EQ.1) ICOUNT=1                              ADS61480
      IF (ICOUNT.GT.ICNDIR) ICOUNT=1                                    ADS61490
C     STEEPEST DESCENT.                                                 ADS61500
120   IF (ICOUNT.EQ.1) SCALS=0.                                         ADS61510
      IF (ICOUNT.EQ.1) DF2=1.                                           ADS61520
C     CONJUGATE DIRECTION.                                              ADS61530
      ICALC=3                                                           ADS61540
      CALL ADS236 (NDV,ICALC,DF,DM1,DM2,DM3,DF2,S,S2,SCALS,ZRO,ISIDC)   ADS61550
C     SLOPE OF F(ALPHA) FUNCTION.                                       ADS61560
      CALL ADS501 (NDV,DF,S,SLOPF)                                      ADS61570
      ZROM=-ZRO                                                         ADS61580
      IF (SLOPF.LT.ZROM.OR.ICOUNT.EQ.1) GO TO 160                       ADS61590
C     SLOPE IS BAD (POSITIVE).  GO USE STEEPEST DESCENT.                ADS61600
      ICOUNT=1                                                          ADS61610
      GO TO 120                                                         ADS61620
130   CONTINUE                                                          ADS61630
C     CONSTRAINTS ARE ACTIVE/VIOLATED.  FIND USABLE-FEASIBLE            ADS61640
C     SEARCH DIRECTION.                                                 ADS61650
      ICOUNT=0                                                          ADS61660
      IF (NVC.EQ.0) PHI=5.                                              ADS61670
      PHI=10.*PHI                                                       ADS61680
      IF (PHI.GT.1000.) PHI=1000.                                       ADS61690
C     CALCULATE SCALE FACTOR FOR GRADIENT OF OBJECTIVE.                 ADS61700
      NN=1                                                              ADS61710
      CALL ADS207 (NDV,NN,DF,NDV,ZRO,SCALF)                             ADS61720
C     CALCULATE NORMALIZATION FACTORS FOR CONSTRAINTS.                  ADS61730
C     STORE FACTORS IN GGL.                                             ADS61740
      IF (NAVC.GT.0) CALL ADS207 (NDV,NAVC,A,NRA,ZRO,GGL)               ADS61750
C     CALCULATE USABLE-FEASIBLE SEARCH DIRECTION.                       ADS61760
      THETA1=THETAZ                                                     ADS61770
      THETAZ=0.                                                         ADS61780
      IF (NVC.GT.0) THETAZ=THETA1                                       ADS61790
      CALL ADS217 (NDV,NAC,NVC,NACS,CT,CTMIN,CTLMIN,THETAZ,ZRO,PHI,GGL,GADS61800
     1,DF,A,NRA,IDG,IC,ISIDC,ITMP,B,NRB,SCALF,S,THETA,C,BETA,RSDMX)     ADS61810
      THETAZ=THETA1                                                     ADS61820
      IF (IPRNT2.GT.3) WRITE (6,680) BETA,RSDMX                         ADS61830
C     ------------------------------------------------------------------ADS61840
C     CHECK KUHN-TUCKER CONDITIONS.                                     ADS61850
C     ------------------------------------------------------------------ADS61860
      IF (NVC.GT.0.AND.BETA.GT.0.01) GO TO 150                          ADS61870
C     IF NDV SIDE CONSTRAINTS ARE ACTIVE AND BETA=0, TERMINATE.         ADS61880
      IF (NACS.EQ.NDV.AND.BETA.LE.0.01) GO TO 460                       ADS61890
      IF (RSDMX.GT.STOL.AND.BETA.GT.ZRO10) GO TO 150                    ADS61900
C     EITHER RSDMX OR BETA INDICATES K-T CONDITIONS ARE SATISFIED.      ADS61910
C     CHECK CONSTRAINT VALUES.                                          ADS61920
      ALPHA=0.                                                          ADS61930
      DO 140 I=1,NAVC                                                   ADS61940
      II=IC(I)                                                          ADS61950
      JJ=IABS(IDG(II))                                                  ADS61960
      CC=CTMIN                                                          ADS61970
      IF (JJ.GT.1) CC=CTLMIN                                            ADS61980
C     IF G.GT.CC, GO CHECK CONVERGENCE.                                 ADS61990
      IF (ABS(G(II)).GT.CC) GO TO 420                                   ADS62000
140   CONTINUE                                                          ADS62010
C     K-T CONDITIONS ARE SATISFIED.  GO EXIT.                           ADS62020
      GO TO 460                                                         ADS62030
150   CONTINUE                                                          ADS62040
C     SLOPE OF F(ALPHA) FUNCTION.                                       ADS62050
      CALL ADS501 (NDV,DF,S,SLOPF)                                      ADS62060
160   CONTINUE                                                          ADS62070
      IF (IPRNT2.LT.4) GO TO 170                                        ADS62080
      ICALC=2                                                           ADS62090
      CALL ADS402 (ICALC,S,NDV,SCFX)                                    ADS62100
170   CONTINUE                                                          ADS62110
C     ------------------------------------------------------------------ADS62120
C     ESTIMATE ALPHA.                                                   ADS62130
C     ------------------------------------------------------------------ADS62140
      IF (ITER.GT.1) GO TO 180                                          ADS62150
C     FIRST ITERATION ESTIMATE.                                         ADS62160
      CALL ADS213 (OBJ,SLOPF,DLOBJ1,DLOBJ2,ALPHA)                       ADS62170
      CALL ADS214 (NDV,X,S,DX1,DX2,ALP1)                                ADS62180
      IF (ALP1.LT.ALPHA) ALPHA=ALP1                                     ADS62190
      ALPSAV=ALPHA                                                      ADS62200
180   CONTINUE                                                          ADS62210
      OBJ1=OBJ                                                          ADS62220
      ALP1=ALPHA                                                        ADS62230
      ALPHA=.5*(ALPHA+ALPSAV)                                           ADS62240
      ALPSAV=ALP1                                                       ADS62250
      IF (ALPHA.LE.0.001) ALPHA=0.001                                   ADS62260
C     ESTIMATE ALPMAX TO ENCOUNTER NDV SIDE CONSTRAINTS.                ADS62270
      CALL ADS112 (NDV,X,S,SVLB,SVUB,ZRO,ALPMAX)                        ADS62280
      IF (ALPMAX.LE.ZRO) GO TO 460                                      ADS62290
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS62300
C     ------------------------------------------------------------------ADS62310
C     CALCULATE DG/DALPHA VECTOR AND MOVE TO ENCOUNTER NEW CONSTRAINT,. ADS62320
C     OR OVERCOME CONSTRAINT VIOLATION.                                 ADS62330
C     ------------------------------------------------------------------ADS62340
      ALP1=5.*ALPHA                                                     ADS62350
      ALP2=-1.                                                          ADS62360
      II=0                                                              ADS62370
      DO 220 I=1,NCON                                                   ADS62380
      DG(I)=0.                                                          ADS62390
      DO 190 J=1,NGT                                                    ADS62400
      J1=J                                                              ADS62410
      IF (IC(J).EQ.I) GO TO 200                                         ADS62420
190   CONTINUE                                                          ADS62430
C     THIS GRADIENT IS NOT AVAILABLE.                                   ADS62440
      GO TO 220                                                         ADS62450
200   CONTINUE                                                          ADS62460
C     GRADIENT IS AVAILABLE.  CALCULATE DG(I)                           ADS62470
      DO 210 J=1,NDV                                                    ADS62480
210   DG(I)=DG(I)+A(J,J1)*S(J)                                          ADS62490
      IF (ABS(DG(I)).LE.0.001) GO TO 220                                ADS62500
      ALP=-G(I)/DG(I)                                                   ADS62510
      IF (ALP.LT.ALP1.AND.ALP.GT.ZRO.AND.J1.GT.NAVC) ALP1=ALP           ADS62520
      IF(ALP.GT.ALP2.AND.J1.LE.NAVC) ALP2=ALP                           ADS62530
      II=1                                                              ADS62540
220   CONTINUE                                                          ADS62550
      IF(NVC.GT.0.AND.ALP2.GT.0.) ALPHA=ALP2                            ADS62560
      IF (ALP1.LT.ALP2.AND.II.EQ.1) ALPHA=ALP1                          ADS62570
C     LIMIT FIRST ESTIMATE TO MOVE TO SIDE CONSTRAINT.                  ADS62580
      ALP=1.0E+20                                                       ADS62590
      DO 230 I=1,NDV                                                    ADS62600
      IF (ABS(S(I)).LE.ZRO) GO TO 230                                   ADS62610
      CC=(SVLB(I)-X(I))/S(I)                                            ADS62620
      IF (CC.GT.ZRO.AND.CC.LT.ALP) ALP=CC                               ADS62630
      CC=(SVUB(I)-X(I))/S(I)                                            ADS62640
      IF (CC.GT.ZRO.AND.CC.LT.ALP) ALP=CC                               ADS62650
230   CONTINUE                                                          ADS62660
      IF (ALP.LT.ALPHA) ALPHA=ALP                                       ADS62670
C     ------------------------------------------------------------------ADS62680
C     CREATE B-MATRIX AND INVERT.                                       ADS62690
C     ------------------------------------------------------------------ADS62700
      IF (NVC.GT.0.OR.NACT.EQ.0) GO TO 320                              ADS62710
C     B-MATRIX.                                                         ADS62720
      CALL ADS209 (NDV,NAVC,NACS,A,NRA,B,NRB,ISIDC)                     ADS62730
C     ZERO OUT ROWS AND COLUMNS OF B FOR DESIGNS LEAVING A CONSTRAINT.  ADS62740
      IF (NAVC.EQ.0) GO TO 260                                          ADS62750
      DO 250 I=1,NAVC                                                   ADS62760
      I1=IC(I)                                                          ADS62770
      IF (DG(I1).GT.-0.01) GO TO 250                                    ADS62780
C     ZERO THIS ROW AND COLUMN.                                         ADS62790
      DO 240 J=1,NACT                                                   ADS62800
      B(I,J)=0.                                                         ADS62810
240   B(J,I)=0.                                                         ADS62820
250   CONTINUE                                                          ADS62830
260   CONTINUE                                                          ADS62840
      IF (NACS.EQ.0) GO TO 310                                          ADS62850
C     SIDE CONSTRAINTS.                                                 ADS62860
      NN=NAVC                                                           ADS62870
      DO 300 I=1,NDV                                                    ADS62880
      ID=ISIDC(I)                                                       ADS62890
      IF (ID.EQ.0) GO TO 300                                            ADS62900
      IF (ID.EQ.1) GO TO 280                                            ADS62910
C     LOWER BOUND.                                                      ADS62920
      NN=NN+1                                                           ADS62930
      IF (S(I).LE.ZRO) GO TO 280                                        ADS62940
C     LEAVING THIS SIDE.  ZERO ROW AND COLUMN OF B.                     ADS62950
      DO 270 J=1,NACT                                                   ADS62960
      B(NN,J)=0.                                                        ADS62970
270   B(J,NN)=0.                                                        ADS62980
280   NACT=NAVC+NACS                                                    ADS62990
      IF (ID.LT.0) GO TO 300                                            ADS63000
C     UPPER BOUND.                                                      ADS63010
      NN=NN+1                                                           ADS63020
      IF (S(I).GE.-ZRO) GO TO 300                                       ADS63030
C     LEAVING THIS SIDE.  ZERO ROW AND COLUMN OF B.                     ADS63040
      DO 290 J=1,NACT                                                   ADS63050
      B(NN,J)=0.                                                        ADS63060
290   B(J,NN)=0.                                                        ADS63070
300   CONTINUE                                                          ADS63080
310   CONTINUE                                                          ADS63090
C     INVERT B.                                                         ADS63100
      CALL ADS202 (NACT,ZRO,B,NRB)                                      ADS63110
320   CONTINUE                                                          ADS63120
C     ------------------------------------------------------------------ADS63130
C     DO ONE-DIMENSIONAL SEARCH.                                        ADS63140
C     ------------------------------------------------------------------ADS63150
      LGOTO=0                                                           ADS63160
      KCOUNT=0                                                          ADS63170
      IF (IPRNT2.GE.4) WRITE (6,690) ALPHA                              ADS63180
      GMX1=1.0E+20                                                      ADS63190
      GO TO 410                                                         ADS63200
330   CONTINUE                                                          ADS63210
      NACT=NAVC+NACS                                                    ADS63220
      IF (NVC.GT.0.OR.NACT.EQ.0) GO TO 410                              ADS63230
      KCOUNT=KCOUNT+1                                                   ADS63240
      IF (KCOUNT.GT.10) GO TO 410                                       ADS63250
      IQ=0                                                              ADS63260
      GMX=0.                                                            ADS63270
C     BUILD GTMP.                                                       ADS63280
      IF (NAVC.EQ.0) GO TO 350                                          ADS63290
      DO 340 I=1,NAVC                                                   ADS63300
      I1=IC(I)                                                          ADS63310
      I2=IDG(I1)                                                        ADS63320
      CC=CTMIN                                                          ADS63330
      IF (IABS(I2).GT.1) CC=CTLMIN                                      ADS63340
      GI=G(I1)                                                          ADS63350
      IF (DG(I1).LT.-0.002) GI=0.                                       ADS63360
      GIA=ABS(GI)                                                       ADS63370
      IF (GIA.GT.CC) IQ=1                                               ADS63380
      IF (GIA.GT.GMX) GMX=GIA                                           ADS63390
340   GTMP(I)=GI                                                        ADS63400
350   CONTINUE                                                          ADS63410
C     SIDE CONSTRAINTS.                                                 ADS63420
      IF (NACS.EQ.0) GO TO 380                                          ADS63430
      NN=NAVC                                                           ADS63440
      DO 370 I=1,NDV                                                    ADS63450
      I1=ISIDC(I)                                                       ADS63460
      IF (I1.EQ.0) GO TO 370                                            ADS63470
      IF (I1.EQ.1) GO TO 360                                            ADS63480
C     LOWER BOUND IS ACTIVE.                                            ADS63490
      NN=NN+1                                                           ADS63500
      GTMP(NN)=SVLB(I)-X(I)                                             ADS63510
      GIA=ABS(GTMP(NN))                                                 ADS63520
      IF (GIA.GT.GMX) GMX=GIA                                           ADS63530
      IF (GIA.GT.CTLMIN) IQ=1                                           ADS63540
      IF (I1.LT.0) GO TO 370                                            ADS63550
360   CONTINUE                                                          ADS63560
C     UPPER BOUND IS ACTIVE.                                            ADS63570
      NN=NN+1                                                           ADS63580
      GTMP(NN)=X(I)-SVUB(I)                                             ADS63590
      GIA=ABS(GTMP(NN))                                                 ADS63600
      IF (GIA.GT.GMX) GMX=GIA                                           ADS63610
      IF (GIA.GT.CTLMIN) IQ=1                                           ADS63620
370   CONTINUE                                                          ADS63630
380   CONTINUE                                                          ADS63640
      GMX2=.9*GMX1                                                      ADS63650
      IF (GMX.GT.GMX2.OR.IQ.EQ.0) GO TO 410                             ADS63660
      GMX1=GMX                                                          ADS63670
C     CALCULATE DELTA-X.                                                ADS63680
      CALL ADS203 (NDV,NAVC,NACS,A,NRA,B,NRB,XTMP,GTMP,ISIDC)           ADS63690
C     ITERATE.                                                          ADS63700
C     UPDATE X.                                                         ADS63710
      XMX=0.                                                            ADS63720
      DO 390 I=1,NDV                                                    ADS63730
      IF (ABS(XTMP(I)).GT.XMX) XMX=ABS(XTMP(I))                         ADS63740
390   CONTINUE                                                          ADS63750
      IF (XMX.LE.ZRO) GO TO 410                                         ADS63760
      DO 400 I=1,NDV                                                    ADS63770
      X(I)=X(I)+XTMP(I)                                                 ADS63780
      IF (X(I).LT.SVLB(I)) X(I)=SVLB(I)                                 ADS63790
      IF (X(I).GT.SVUB(I)) X(I)=SVUB(I)                                 ADS63800
400   CONTINUE                                                          ADS63810
      ITEMP(1)=ICOUNT                                                   ADS63820
      ITEMP(2)=IDAB                                                     ADS63830
      ITEMP(3)=IDEL                                                     ADS63840
      ITEMP(4)=ITER                                                     ADS63850
      ITEMP(5)=KOUNT                                                    ADS63860
      ITEMP(6)=LGOTO                                                    ADS63870
      ITEMP(7)=NACS                                                     ADS63880
      ITEMP(8)=NAVC                                                     ADS63890
      ITEMP(9)=NGRAD                                                    ADS63900
      ITEMP(10)=ISTRAT                                                  ADS63910
      RETURN                                                            ADS63920
410   CONTINUE                                                          ADS63930
      GMX1=1.0E+20                                                      ADS63940
      CALL ADS001 (LGOTO,ISTRAT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,   ADS63950
     1NAVC,IC,WK,IWK,SLOPF,ALPHA,ALPMAX)                                ADS63960
      IF (LGOTO.EQ.0) GO TO 420                                         ADS63970
C     CALCULATE FUNCTION VALUES.                                        ADS63980
      KCOUNT=0                                                          ADS63990
      IGOTO=2                                                           ADS64000
      INFO=1                                                            ADS64010
C     PUT INFO BACK IN TEMP FOR TRANSFER.                               ADS64020
      TEMP(1)=ALPSAV                                                    ADS64030
      TEMP(2)=BETA                                                      ADS64040
      TEMP(3)=DF2                                                       ADS64050
      TEMP(4)=OBJ1                                                      ADS64060
      TEMP(5)=PHI                                                       ADS64070
      TEMP(6)=DLOBJ1                                                    ADS64080
      TEMP(7)=DLOBJ2                                                    ADS64090
      TEMP(8)=DX1                                                       ADS64100
      TEMP(9)=DX2                                                       ADS64110
      TEMP(10)=GMAX1                                                    ADS64120
      ITEMP(1)=ICOUNT                                                   ADS64130
      ITEMP(2)=IDAB                                                     ADS64140
      ITEMP(3)=IDEL                                                     ADS64150
      ITEMP(4)=ITER                                                     ADS64160
      ITEMP(5)=KOUNT                                                    ADS64170
      ITEMP(6)=LGOTO                                                    ADS64180
      ITEMP(7)=NACS                                                     ADS64190
      ITEMP(8)=NAVC                                                     ADS64200
      ITEMP(9)=NGRAD                                                    ADS64210
      ITEMP(10)=ISTRAT                                                  ADS64220
      RETURN                                                            ADS64230
420   CONTINUE                                                          ADS64240
C     ONE-DIMENSIONAL SEARCH IS COMPLETE.                               ADS64250
      OBJA=OBJ/SCFO                                                     ADS64260
      IF (IPRNT2.EQ.2) WRITE (6,730) ITER,OBJA                          ADS64270
      IF (IPRNT2.GE.4) WRITE (6,740) ALPHA                              ADS64280
      IF (IPRNT2.GT.2) WRITE (6,700) OBJA                               ADS64290
      ICALC=1                                                           ADS64300
      IF (IPRNT2.GE.2) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS64310
      ICALC=4                                                           ADS64320
      IF (IPRNT2.GT.2) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS64330
C     ------------------------------------------------------------------ADS64340
C     CHECK CONVERGENCE.                                                ADS64350
C     ------------------------------------------------------------------ADS64360
      NGRAD=1                                                           ADS64370
      CALL ADS211 (OBJ,OBJ1,DELOBJ,DABOBJ,DEL1,DAB1,IDEL,IDAB)          ADS64380
      IF (ITER.GE.ITMAX) GO TO 460                                      ADS64390
      IF (NVC.GT.0.AND.ITER.LT.10) IDEL=0                               ADS64400
      IF (NVC.GT.0.AND.ITER.LT.10) IDAB=0                               ADS64410
      IF (IDEL.GE.ITRMOP.OR.IDAB.GE.ITRMOP) GO TO 460                   ADS64420
C     REDUCE CT AND CTL IF NO PROGRESS IS BEING MADE.                   ADS64430
      IF (IDEL.EQ.0.AND.IDAB.EQ.0) GO TO 425                            ADS64440
      CT=.5*CT                                                          ADS64450
      IF(ABS(CT).LT.CTMIN) CT=-CTMIN                                    ADS64460
      CTL=.5*CTL                                                        ADS64470
      IF(ABS(CTL).LT.CTLMIN) CTL=-CTLMIN                                ADS64480
425   CONTINUE                                                          ADS64490
C     IF LAST MOVE WAS GT.ZRO10, GO START NEW ITERATION.                ADS64500
      IF (ALPHA.GT.ZRO10.OR.NAVC.EQ.0) GO TO 40                         ADS64510
C     ------------------------------------------------------------------ADS64520
C     LAST ALPHA.LE.ZRO10.  DELETE CONSTRAINTS THAT ARE NO LONGER       ADS64530
C     CRITICAL AND CONTINUE OPTIMIZATION WITHOUT GETTING NEW GRADIENTS. ADS64540
C     ------------------------------------------------------------------ADS64550
      IF (NGT.GE.NCON) GO TO 450                                        ADS64560
C     CHECK FOR NEW CRITICAL CONSTRAINTS.                               ADS64570
      DO 440 I=1,NCON                                                   ADS64580
      I1=IDG(I)                                                         ADS64590
      IF (I1.LT.0) GO TO 440                                            ADS64600
      CC=CT                                                             ADS64610
      IF (I1.GT.1) CC=CTL                                               ADS64620
      IF (G(I).LT.CC) GO TO 440                                         ADS64630
C     THIS CONSTRAINT IS CRITICAL.  IS IT ALREADY IN IC?                ADS64640
      DO 430 J=1,NGT                                                    ADS64650
      IF (IC(J).EQ.I) GO TO 440                                         ADS64660
430   CONTINUE                                                          ADS64670
C     NEW CRITICAL CONSTRAINT.  MUST EVALUATE GRADIENTS.                ADS64680
      NGRAD=1                                                           ADS64690
      GO TO 40                                                          ADS64700
440   CONTINUE                                                          ADS64710
450   CONTINUE                                                          ADS64720
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CT,CTMIN,CTL,CTLMIN,NAC,NVC)  ADS64730
      NAVC=NAC+NVC                                                      ADS64740
C     WE HAVE GRADIENTS OF ALL CRITICAL CONSTRAINTS.                    ADS64750
C     CONTINUE OPTIMIZATION WITHOUT GETTING NEW GRADIENTS.              ADS64760
      NGRAD=0                                                           ADS64770
      GO TO 40                                                          ADS64780
460   CONTINUE                                                          ADS64790
C     ------------------------------------------------------------------ADS64800
C     OPTIMIZATION IS COMPLETE.                                         ADS64810
C     ------------------------------------------------------------------ADS64820
      IF (IPRNT2.LE.0) GO TO 500                                        ADS64830
C     PRINT FINAL RESULTS.                                              ADS64840
      WRITE (6,530)                                                     ADS64850
      WRITE (6,540) ITER                                                ADS64860
      OBJA=OBJ/SCFO                                                     ADS64870
      WRITE (6,700) OBJA                                                ADS64880
C   X-VECTOR.                                                           ADS64890
      ICALC=1                                                           ADS64900
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS64910
C     G-VECTOR.                                                         ADS64920
      ICALC=4                                                           ADS64930
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS64940
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS64950
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS64960
C     SIDE CONSTRAINTS.                                                 ADS64970
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS64980
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS64990
      WRITE (6,520) CT,CTL                                              ADS65000
C     WRITE ACTIVE/VIOLATED COSNTRAINT NUMBERS.                         ADS65010
      NAVC=NAC+NVC                                                      ADS65020
      WRITE (6,620) NAC,NVC                                             ADS65030
      IF (NAVC.GT.0) WRITE (6,630)                                      ADS65040
      IF (NAVC.GT.0) WRITE (6,640) (ITMP(I),I=1,NAVC)                   ADS65050
      WRITE (6,650) NACS                                                ADS65060
      IF (NACS.EQ.0) GO TO 490                                          ADS65070
      N=0                                                               ADS65080
      DO 480 I=1,NDV                                                    ADS65090
      IF (ISIDC(I).EQ.0) GO TO 480                                      ADS65100
      IF (ISIDC(I).EQ.1) GO TO 470                                      ADS65110
C     LOWER BOUND IS ACTIVE.                                            ADS65120
      N=N+1                                                             ADS65130
      ITMP(N)=-I                                                        ADS65140
470   IF (ISIDC(I).LT.0) GO TO 480                                      ADS65150
C     UPPER BOUND IS ACTIVE.                                            ADS65160
      N=N+1                                                             ADS65170
      ITMP(N)=I                                                         ADS65180
480   CONTINUE                                                          ADS65190
      WRITE (6,660)                                                     ADS65200
      WRITE (6,640) (ITMP(I),I=1,NACS)                                  ADS65210
490   WRITE (6,550)                                                     ADS65220
      IF (ITER.GE.ITMAX) WRITE (6,560)                                  ADS65230
      IF (IDEL.GE.ITRMOP) WRITE (6,570) ITRMOP                          ADS65240
      IF (IDAB.GE.ITRMOP) WRITE (6,580) ITRMOP                          ADS65250
      CC=.001                                                           ADS65260
      IF (RSDMX.LE.0.001) WRITE (6,590) RSDMX,CC                        ADS65270
      NACT=NAVC+NACS                                                    ADS65280
      IF (BETA.LE.0.001.AND.NACT.GT.0) WRITE (6,600) BETA,CC            ADS65290
500   CONTINUE                                                          ADS65300
C     PUT INFO BACK IN TEMP FOR TRANSFER.                               ADS65310
      TEMP(1)=ALPSAV                                                    ADS65320
      TEMP(2)=BETA                                                      ADS65330
      TEMP(3)=DF2                                                       ADS65340
      TEMP(4)=OBJ1                                                      ADS65350
      TEMP(5)=PHI                                                       ADS65360
      TEMP(6)=DLOBJ1                                                    ADS65370
      TEMP(7)=DLOBJ2                                                    ADS65380
      TEMP(8)=DX1                                                       ADS65390
      TEMP(9)=DX2                                                       ADS65400
      TEMP(10)=GMAX1                                                    ADS65410
      ITEMP(1)=ICOUNT                                                   ADS65420
      ITEMP(2)=IDAB                                                     ADS65430
      ITEMP(3)=IDEL                                                     ADS65440
      ITEMP(4)=ITER                                                     ADS65450
      ITEMP(5)=KOUNT                                                    ADS65460
      ITEMP(6)=LGOTO                                                    ADS65470
      ITEMP(7)=NACS                                                     ADS65480
      ITEMP(8)=NAVC                                                     ADS65490
      ITEMP(9)=NGRAD                                                    ADS65500
      ITEMP(10)=ISTRAT                                                  ADS65510
      IGOTO=0                                                           ADS65520
      RETURN                                                            ADS65530
C     ------------------------------------------------------------------ADS65540
C     FORMATS.                                                          ADS65550
C     ------------------------------------------------------------------ADS65560
510   FORMAT ('1')                                                      ADS65570
520   FORMAT (/,5X,'CONSTRAINT TOLERANCE, CT =',1PE13.5,'  CTL =',      ADS65580
     11PE13.5)                                                          ADS65590
530   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS')                     ADS65600
540   FORMAT (/,5X,'NUMBER OF ITERATIONS =',I5)                         ADS65610
550   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS65620
560   FORMAT (/,5X,'NUMBER OF ITERATIONS = ITMAX')                      ADS65630
570   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS65640
     1' CONSECUTIVE ITERATIONS')                                        ADS65650
580   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS65660
     1' CONSECUTIVE ITERATIONS')                                        ADS65670
590   FORMAT (/,5X,'MAXIMUM K-T RESIDUAL =',1PE13.5,' IS LESS THAN',    ADS65680
     11PE13.5)                                                          ADS65690
600   FORMAT (/,5X,'KUHN-TUCKER PARAMETER, BETA =',1PE13.5,             ADS65700
     1' IS LESS THAN',1PE13.5)                                          ADS65710
610   FORMAT (//,5X,75('-'),/,19X,                                      ADS65720
     1'IOPT = 5; MODIFIED METHOD OF FEASIBLE DIRECTIONS',/,5X,75('-'))  ADS65730
620   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS65740
     1' VIOLATED CONSTRAINTS')                                          ADS65750
630   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS65760
640   FORMAT (5X,10I5)                                                  ADS65770
650   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS65780
660   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS65790
670   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS65800
680   FORMAT (/,5X,'K-T PARAMETERS, BETA =',1PE13.5,'  MAX. RESIDUAL =',ADS65810
     11PE13.5)                                                          ADS65820
690   FORMAT (/,5X,'PROPOSED ALPHA =',1PE13.5)                          ADS65830
700   FORMAT (/,5X,'OBJECTIVE =',1PE12.5)                               ADS65840
710   FORMAT (/,5X,'-- INITIAL DESIGN',//,5X,'OBJ =',1PE12.5)           ADS65850
720   FORMAT (//,5X,'-- BEGIN ITERATION NUMBER',I5)                     ADS65860
730   FORMAT (//,5X,'-- ITERATION',I5,5X,'OBJ =',1PE12.5)               ADS65870
740   FORMAT (/,5X,'CALCULATED ALPHA =',1PE13.5)                        ADS65880
      END                                                               ADS65890
C --- SUBROUTINE ADS235 (IGOTO,INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,ITMAX   ADS65900
C    1 ,ICNDIR,DF,S,ZRO,DLOBJ1,DLOBJ2,DX1,DX2,XOLD,SVLB,SVUB,DELOBJ,    ADS65910
C    2 DABOBJ,ITRM,XTMP,IPRNT2,WK,IWK,G,SCFX,H,DX,DFG,IVAMET,IC,        ADS65920
C    3 IDG,ISIDC,SCFO,IFCALL,IGCALL,ALPMAX)                             ADS65930
C     DIMENSION X(NDV),S(NDV),DF(NDV),XOLD(NDV),SVLB(NDV),SVUB(NDV),    ADS65940
C    1 XTMP(NDV),WK(NRWK),IWK(NRIWK),G(NCON),SCFX(NDV),DX(NDV),DFG(NDV),ADS65950
C    2 H(NDV*(NDV+1)/2),IC(NCOLA),IDG(NCON),ISIDC(NDV)                  ADS65960
C     ------------------------------------------------------------------ADS65970
C --- ROUTINE TO PERFORM UNCONSTRAINED MINIMIZATION BY VARIABLE METRIC  ADS65980
C     METHODS.                                                          ADS65990
C     ------------------------------------------------------------------ADS66000
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 NOV. 1982   ADS66010
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS66020
C --- INPUT                                                             ADS66030
C      IGOTO  - CONTROL PQRAMETER.                                      ADS66040
C      INFO   - CONTROL PARAMETER.                                      ADS66050
C      ISTRAT - STRATEGY IDENTIFIER.                                    ADS66060
C      IONED  - ONE D. SEARCH IDENTIFIER.                               ADS66070
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS66080
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS66090
C      OBJ    - OBJECTIVE.                                              ADS66100
C      X      - DESIGN VARIABLES.                                       ADS66110
C      ITMAX  - MAXIMUM ITERATION NUMBER.                               ADS66120
C      ICNDIR - ITERATION NUMBER RETURN TO STEEPEST DESCENT METHOD.     ADS66130
C      DF     - GRADIENT OF OBJECTIVE.                                  ADS66140
C      S      - ARRAY FOR SEARCH DIRECTION.                             ADS66150
C      ZRO    - ZERO.                                                   ADS66160
C      DLOBJ1,DLOBJ2,DX1,DX2 - PARAMETERS TO ESTIMATE INITIAL ALPHA.    ADS66170
C      XOLD   - FORMER DESIGN VARIABLE.                                 ADS66180
C      SVLB,SVUB - LOWER AND UPPER BOUND OF DESIGN VARIABLES.           ADS66190
C      DELOBJ,DABOBJ,ITRM - PARAMETERS FOR CONVERGENCE CROTERIA.        ADS66200
C      XTMP   - TEMPORARY ARRAY FOR DESIGN VARIABLES.                   ADS66210
C      IPRNT2 - PRINT CONTROL IDENTIFIER FOR OPTIMIZATION LEVEL.        ADS66220
C      WK,IWK - REAL AND INTEGER WORK ARRAY.                            ADS66230
C      G      - CONSTRAINT VALUE.                                       ADS66240
C      SCFX   - SCALE FACTOR FOR DESIGN VARIABLES.                      ADS66250
C      H      - ARRAY FOR HESSIAN MATRIX.                               ADS66260
C      DX     - ARRAY FOR THE DIFFERENCE BETWEEN X AND XOLD.            ADS66270
C      DFG    - ARRAY FOR THE DIFFERENCE BETWEEN DF AND DFOLD.          ADS66280
C      IVAMET - VARIABLE METRIC METHOD CONTROL IDENTIFIER.              ADS66290
C      IC     - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                   ADS66300
C      IDG    - CONSTRAINT TYPE IDENTIFIER.                             ADS66310
C      ISIDC  - ARRAY FOR IDENTIFY ACTIVE SIDE CONSTRAINTS.             ADS66320
C      SCFO   - SCALE FACTOR OF OBJECTIVE.                              ADS66330
C      IFCALL,IGCALL - NUMBER OF FUNCTION AND GRADIENT CALL.            ADS66340
C --- OUTPUT                                                            ADS66350
C      OBJ    - OPTIMUM OBJECTIVE.                                      ADS66360
C      X      - OPTIMUM DESIGN VARIABLES.                               ADS66370
C      G      - CONSTRAIN VALUE CORESPOND TO  OPTIMUM DESIGN VARIABLES. ADS66380
C                                                                       ADS66390
      SUBROUTINE ADS235 (IGOTO,INFO,ISTRAT,IONED,NDV,NCON,OBJ,X,ITMAX,ICADS66400
     1NDIR,DF,S,ZRO,DLOBJ1,DLOBJ2,DX1,DX2,XOLD,SVLB,SVUB,DELOBJ,DABOBJ,IADS66410
     2TRM,XTMP,IPRNT2,WK,IWK,G,SCFX,H,DX,DFG,IVAMET,IC,IDG,ISIDC,SCFO,IFADS66420
     3CALL,IGCALL,ALPHA,ALPSAV,DF2,OBJ1,ICOUNT,IDAB,IDEL,IDIVCT,ITER,LGOADS66430
     4TO,ALPMAX,SCFG)                                                   ADS66440
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),S(1),DF(1),XOLD(1),SVLB(1),SVUB(1),XTMP(1),WK(1),IWADS66450
     1K(1),G(1),SCFX(1),DX(1),DFG(1),H(1),IC(1),IDG(1),ISIDC(1),SCFG(1) ADS66460
C                                                                       ADS66470
      IGOTO1=IGOTO+1                                                    ADS66480
      GO TO (10,40,170), IGOTO1                                         ADS66490
10    CONTINUE                                                          ADS66500
C     INITIAL CALL.                                                     ADS66510
C     PRINT AS REQUIRED.                                                ADS66520
      IF (IPRNT2.GT.0.AND.ISTRAT.EQ.0) WRITE (6,280)                    ADS66530
      IF (ICNDIR.EQ.1) GO TO 20                                         ADS66540
      IF (IPRNT2.GT.0.AND.IVAMET.LE.1) WRITE (6,290)                    ADS66550
      IF (IPRNT2.GT.0.AND.IVAMET.EQ.2) WRITE (6,300)                    ADS66560
20    IF (IPRNT2.GT.0.AND.ICNDIR.EQ.1) WRITE (6,310)                    ADS66570
C     SET UP PARAMETERS.                                                ADS66580
      ITER=0                                                            ADS66590
      ICOUNT=0                                                          ADS66600
      IDEL=0                                                            ADS66610
      IDAB=0                                                            ADS66620
      ALPSAV=0.                                                         ADS66630
C     IDENTIFIER FOR CONTROL OF DIVIDE CHECK ERROR DUE TO DFG , DX.     ADS66640
      IDIVCT=0                                                          ADS66650
      IF (IPRNT2.EQ.0) GO TO 30                                         ADS66660
      OBJW=OBJ                                                          ADS66670
      IF (ISTRAT.EQ.0) OBJW=OBJ/SCFO                                    ADS66680
      WRITE (6,320) OBJW                                                ADS66690
      ICALC=1                                                           ADS66700
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS66710
      ICALC=4                                                           ADS66720
      IF (NCON.NE.0) CALL ADS402 (ICALC,G,NCON,SCFG)                    ADS66730
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS66740
C     ------------------------------------------------------------------ADS66750
C     OPTIMIZATION.  BEGIN NEW ITERATION.                               ADS66760
C     ------------------------------------------------------------------ADS66770
30    CONTINUE                                                          ADS66780
      ITER=ITER+1                                                       ADS66790
      IF (IPRNT2.GT.2) WRITE (6,340) ITER                               ADS66800
      OBJ1=OBJ                                                          ADS66810
C     CALCULATE GRADIENT OF OBJECTIVE.                                  ADS66820
      IGOTO=1                                                           ADS66830
      INFO=2                                                            ADS66840
      IF (ISTRAT.EQ.0.AND.ITER.EQ.1) GO TO 40                           ADS66850
      IF (ITER.GT.1.AND.ALPHA.LE.ZRO) GO TO 40                          ADS66860
      RETURN                                                            ADS66870
40    CONTINUE                                                          ADS66880
      OBJ=OBJ1                                                          ADS66890
      ICALC=3                                                           ADS66900
      IF (IPRNT2.GT.4) CALL ADS402 (ICALC,DF,NDV,SCFX)                  ADS66910
C     IDENTIFY ACTIVE SIDE CONSTRAINTS.                                 ADS66920
      CALL ADS509 (NDV,X,SVLB,SVUB,ZRO,ISIDC,NACS)                      ADS66930
C     SEARCH DIRECTION.                                                 ADS66940
      ICOUNT=ICOUNT+1                                                   ADS66950
      IF (ICOUNT.GT.ICNDIR) ICOUNT=1                                    ADS66960
C     STEEPEST DESCENT.                                                 ADS66970
50    IF (ICOUNT.EQ.1) JVAMET=0                                         ADS66980
      IF (ICOUNT.EQ.1) CALL ADS236 (NDV,JVAMET,DF,DX,DFG,H,DF2,S,S2,SCALADS66990
     1S,ZRO,ISIDC)                                                      ADS67000
C     DIRECTION BY VARIABLE METRIC METHODS.                             ADS67010
C     CALCULATE DFG.                                                    ADS67020
      IF (ICOUNT.LE.1) GO TO 130                                        ADS67030
      SQY=0.                                                            ADS67040
      SQDX=0.                                                           ADS67050
      DO 110 I=1,NDV                                                    ADS67060
      DFG(I)=DF(I)-DFG(I)                                               ADS67070
      SQY=SQY+DFG(I)*DFG(I)                                             ADS67080
      IF (ISIDC(I).EQ.0) GO TO 100                                      ADS67090
      DXOL=ABS(DX(I)-SVLB(I))                                           ADS67100
      DXOU=ABS(SVUB(I)-DX(I))                                           ADS67110
      VLIM=DMIN1(DXOL,DXOU)                                             ADS67120
      IF (VLIM.GT.ZRO) GO TO 70                                         ADS67130
      IF (DXOL.LE.ZRO.AND.DF(I).LE.-ZRO) GO TO 60                       ADS67140
      IF (DXOU.LE.ZRO.AND.DF(I).GE.ZRO) GO TO 60                        ADS67150
      GO TO 70                                                          ADS67160
60    ICOUNT=1                                                          ADS67170
      GO TO 50                                                          ADS67180
70    CONTINUE                                                          ADS67190
      NA=I*(I-1)/2                                                      ADS67200
      DX(I)=0.                                                          ADS67210
      J1=NA+1                                                           ADS67220
      J2=NA+I                                                           ADS67230
      DO 80 J=J1,J2                                                     ADS67240
80    H(J)=0.                                                           ADS67250
      JZ=NA+I                                                           ADS67260
      J1=I+1                                                            ADS67270
      DO 90 J=J1,NDV                                                    ADS67280
      JZ=JZ+J-1                                                         ADS67290
90    H(JZ)=0.                                                          ADS67300
      GO TO 110                                                         ADS67310
100   CONTINUE                                                          ADS67320
      DX(I)=X(I)-DX(I)                                                  ADS67330
      SQDX=SQDX+DX(I)*DX(I)                                             ADS67340
110   CONTINUE                                                          ADS67350
      SQY=SQRT(SQY)                                                     ADS67360
      SQDX=SQRT(SQDX)                                                   ADS67370
      IF (SQY.GT.ZRO.AND.SQDX.GT.ZRO) GO TO 120                         ADS67380
      IDIVCT=IDIVCT+1                                                   ADS67390
      IF (IDIVCT.GE.2) GO TO 200                                        ADS67400
      ICOUNT=1                                                          ADS67410
      GO TO 50                                                          ADS67420
120   CONTINUE                                                          ADS67430
      IF (ICOUNT.EQ.2) IDIVCT=0                                         ADS67440
      ICALC=8                                                           ADS67450
      IF (IPRNT2.GT.4) CALL ADS402 (ICALC,DX,NDV,SCFX)                  ADS67460
      ICALC=9                                                           ADS67470
      IF (IPRNT2.GT.4) CALL ADS402 (ICALC,DFG,NDV,SCFX)                 ADS67480
      IF (ICOUNT.GT.1) CALL ADS236 (NDV,IVAMET,DF,DX,DFG,H,DF2,S,S2,SCALADS67490
     1S,ZRO,ISIDC)                                                      ADS67500
130   CONTINUE                                                          ADS67510
      IF (IPRNT2.GE.5) CALL ADS404 (NDV,H)                              ADS67520
      ICALC=2                                                           ADS67530
      IF (IPRNT2.GT.3) CALL ADS402 (ICALC,S,NDV,SCFX)                   ADS67540
C     ONE-DIMENSIONAL SEARCH.                                           ADS67550
C     SLOPE OF F(ALPHA) FUNCTION.                                       ADS67560
      CALL ADS501 (NDV,DF,S,SLOPF)                                      ADS67570
      ZROM=-ZRO                                                         ADS67580
      IF (SLOPF.LT.ZROM.OR.ICOUNT.EQ.1) GO TO 150                       ADS67590
C     SLOPE IS BAD (POSITIVE).  RETURN TO STEEPEST DESCENT.             ADS67600
C     IF ALL VARIABLES ARE EQUAL TO THEIR LIMIT VALUES , THEN TERMINATE.ADS67610
      ICOUNT=1                                                          ADS67620
      IISID=0                                                           ADS67630
      DO 140 I=1,NDV                                                    ADS67640
      IF (ISIDC(I).EQ.0) GO TO 50                                       ADS67650
140   CONTINUE                                                          ADS67660
      IISID=1                                                           ADS67670
      GO TO 200                                                         ADS67680
150   CONTINUE                                                          ADS67690
C     IF FIRST ITERATION, ESTIMATE ALPHA.                               ADS67700
      IF (ITER.GT.1) GO TO 160                                          ADS67710
      CALL ADS213 (OBJ,SLOPF,DLOBJ1,DLOBJ2,ALPHA)                       ADS67720
      CALL ADS214 (NDV,X,S,DX1,DX2,ALP1)                                ADS67730
      IF (ALP1.LT.ALPHA) ALPHA=ALP1                                     ADS67740
160   CONTINUE                                                          ADS67750
      ALP1=ALPSAV                                                       ADS67760
      ALPSAV=ALPHA                                                      ADS67770
      ALPHA=.5*(ALPHA+ALP1)                                             ADS67780
      IF (ALPHA.LT.ZRO) ALPHA=ZRO                                       ADS67790
C     ALPMAX.                                                           ADS67800
      CALL ADS112 (NDV,X,S,SVLB,SVUB,ZRO,ALPMAX)                        ADS67810
      IF (ALPMAX.LE.ZRO) GO TO 200                                      ADS67820
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS67830
C     ------------------------------------------------------------------ADS67840
C     DO ONE-DIMENSIONAL SEARCH.                                        ADS67850
C     ------------------------------------------------------------------ADS67860
C     SAVE XOLD IN DX.                                                  ADS67870
      CALL ADS503 (NDV,X,DX)                                            ADS67880
C     SAVE DF IN DFG.                                                   ADS67890
      CALL ADS503 (NDV,DF,DFG)                                          ADS67900
      LGOTO=0                                                           ADS67910
      IF (IPRNT2.GT.3) WRITE (6,220) ALPHA                              ADS67920
170   CONTINUE                                                          ADS67930
      NACG=0                                                            ADS67940
      CALL ADS001 (LGOTO,ISTRAT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDG,   ADS67950
     1NACG,IC,WK,IWK,SLOPF,ALPHA,ALPMAX)                                ADS67960
      IF (LGOTO.EQ.0) GO TO 180                                         ADS67970
C     EVALUATE OBJECTIVE FUNCTION.                                      ADS67980
      IGOTO=2                                                           ADS67990
      INFO=1                                                            ADS68000
      RETURN                                                            ADS68010
C                                                                       ADS68020
180   CONTINUE                                                          ADS68030
C     ONE-DIMENSIONAL SEARCH IS COMPLETE.                               ADS68040
C     RE-SET COUNTER IF ALPHA=0.                                        ADS68050
      IF (ALPHA.LE.ZRO) ICOUNT=0                                        ADS68060
      OBJW=OBJ                                                          ADS68070
      IF (ISTRAT.EQ.0) OBJW=OBJ/SCFO                                    ADS68080
      IF (IPRNT2.EQ.2) WRITE (6,350) ITER,OBJW                          ADS68090
      IF (IPRNT2.EQ.3) WRITE (6,330) OBJW                               ADS68100
      IF (IPRNT2.GE.4) WRITE (6,360) ALPHA,OBJW                         ADS68110
      IF (IPRNT2.LT.2) GO TO 190                                        ADS68120
      ICALC=1                                                           ADS68130
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS68140
      IF (IPRNT2.LT.3) GO TO 190                                        ADS68150
      ICALC=4                                                           ADS68160
      IF (NCON.NE.0) CALL ADS402 (ICALC,G,NCON,SCFG)                    ADS68170
190   CONTINUE                                                          ADS68180
C     ------------------------------------------------------------------ADS68190
C     CHECK CONVERGENCE.                                                ADS68200
C     ------------------------------------------------------------------ADS68210
      CALL ADS211 (OBJ,OBJ1,DELOBJ,DABOBJ,DEL1,DAB1,IDEL,IDAB)          ADS68220
      IF (ITER.GE.ITMAX) GO TO 200                                      ADS68230
      IF (IDEL.LT.ITRM.AND.IDAB.LT.ITRM) GO TO 30                       ADS68240
200   CONTINUE                                                          ADS68250
C     ------------------------------------------------------------------ADS68260
C     OPTIMIZATION IS COMPLETE.                                         ADS68270
C     ------------------------------------------------------------------ADS68280
      IF (IPRNT2.EQ.0) GO TO 210                                        ADS68290
      WRITE (6,370)                                                     ADS68300
      OBJW=OBJ                                                          ADS68310
      IF (ISTRAT.EQ.0) OBJW=OBJ/SCFO                                    ADS68320
      WRITE (6,330) OBJW                                                ADS68330
      ICALC=1                                                           ADS68340
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS68350
      ICALC=4                                                           ADS68360
      IF (NCON.NE.0) CALL ADS402 (ICALC,G,NCON,SCFG)                    ADS68370
      WRITE (6,230)                                                     ADS68380
      IF (ITER.GE.ITMAX) WRITE (6,240)                                  ADS68390
      IF (IDEL.GE.ITRM) WRITE (6,250) ITRM                              ADS68400
      IF (IDAB.GE.ITRM) WRITE (6,260) ITRM                              ADS68410
      IF (IDIVCT.GE.2.OR.IISID.EQ.1) WRITE (6,270)                      ADS68420
      IF (ALPMAX.LE.ZRO) WRITE (6,270)                                  ADS68430
210   CONTINUE                                                          ADS68440
      IGOTO=0                                                           ADS68450
      RETURN                                                            ADS68460
C     ------------------------------------------------------------------ADS68470
C     FORMATS                                                           ADS68480
C     ------------------------------------------------------------------ADS68490
220   FORMAT (/,5X,'PROPOSED ALPHA =',1PE12.5)                          ADS68500
230   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS68510
240   FORMAT (/,5X,'NUMBER OF ITERATIONS = ITMAX')                      ADS68520
250   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS68530
     1' CONSECUTIVE ITERATIONS')                                        ADS68540
260   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS68550
     1' CONSECUTIVE ITERATIONS')                                        ADS68560
270   FORMAT (/,5X,'KUHN-TUCKER CONDITIONS ARE SATISFIED')              ADS68570
280   FORMAT ('1')                                                      ADS68580
290   FORMAT (//,5X,75('-'),/,20X,                                      ADS68590
     1'IOPT = 2; VARIABLE METRIC METHOD -- DFP ---',/,5X,75('-'))       ADS68600
300   FORMAT (//,5X,75('-'),/,20X,                                      ADS68610
     1'IOPT = 3; VARIABLE METRIC METHOD -- BFGS ---',/,5X,75('-'))      ADS68620
310   FORMAT (//,5X,75('-'),/,25X,' ICNDIR = 1; STEEPEST DESCENT METHOD'ADS68630
     1,/,5X,75('-'))                                                    ADS68640
320   FORMAT (//,5X,'INITIAL DESIGN',//,5X,'OBJ =',1PE12.5)             ADS68650
330   FORMAT (/,5X,'OBJECTIVE =',1PE12.5)                               ADS68660
340   FORMAT (//,5X,'-- BEGIN ITERATION',I3)                            ADS68670
350   FORMAT (//,5X,'-- ITERATION',I5,3X,'OBJ =',1PE12.5)               ADS68680
360   FORMAT (/,5X,'CALCULATED ALPHA =',1PE12.5,//,5X,'OBJECTIVE =',    ADS68690
     11PE12.5)                                                          ADS68700
370   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS')                     ADS68710
      END                                                               ADS68720
C --- SUBROUTINE ADS236 (NDV,ICONT,DF,P,Y,H,DF2,S,S2,SCAL,ZRO,ISIDC)    ADS68730
C     DIMENSION DF(NDV),P(NDV),Y(NDV),H(NDV*(NDV+1)/2),S(NDV),ISIDC(NDV)ADS68740
C     ------------------------------------------------------------------ADS68750
C --- ROUTINE FOR VARIABLE METRIC METHODS.                              ADS68760
C     ------------------------------------------------------------------ADS68770
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 OCT. 1982   ADS68780
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS68790
C --- INPUT                                                             ADS68800
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS68810
C     ICONT - METHOD CONTROL. ICONT=1;DFP , ICONT=2;BFGS , ICONT=3;FR , ADS68820
C         ICONT=0;STEEPEST DESCENT.                                     ADS68830
C     DF - GRADIENT OF THE OBJECTIVE FUNCTION. MINIMUM DIMENSION IS NDV.ADS68840
C     P - DIFFERENCE OF X AND XOLD. MINIMUM DIMENSION IS NDV.           ADS68850
C     Y - DIFFERENCE OF DF AND DFOLD. MINIMUM DIMENSION IS NDV.         ADS68860
C     H - THE FORMER INVERSE HESSIAN MATRIX APPROXIMATED. NDV * NDV     ADS68870
C         MATRIX, UPPER TRIANGLE ONLY, STORED COLUMN BY COLUMN.         ADS68880
C     DF2 - DF**2 FROM THE PREVIOUS OPTIMIZATION ITERATION.             ADS68890
C     S - SEARCH DIRECTION FROM THE PREVIOUS OPTIMIZATION ITERATION.    ADS68900
C     SCAL - SCALING FACTOR FROM THE PREVIOUS OPTIMIZATION ITERATION.   ADS68910
C     ISIDC - =0 ; DESIGN VARIABLE IS NOT SIDE LIMIT.                   ADS68920
C             =1 ; DESIGN VARIABLE IS UPPER SIDE LIMIT.                 ADS68930
C            =-1 ; DESIGN VARIABLE IS LOWER SIDE LIMIT.                 ADS68940
C --- OUTPUT                                                            ADS68950
C     S - NEW SEARCH DIRECTION.                                         ADS68960
C     H - NEW APPROXIMATION TO INVERSE HESSIAN MATRIX.                  ADS68970
C     S2 - S**2, SQUARE OF MAGNITUDE OF S.                              ADS68980
C     SCAL - SCALING FACTOR FOR S.                                      ADS68990
C     DF2 - NEW DF2.                                                    ADS69000
C                                                                       ADS69010
      SUBROUTINE ADS236 (NDV,ICONT,DF,P,Y,H,DF2,S,S2,SCAL,ZRO,ISIDC)    ADS69020
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION DF(1),P(1),Y(1),H(1),S(1),ISIDC(1)                      ADS69030
      IF (ICONT.EQ.0) GO TO 10                                          ADS69040
      IF (ICONT.EQ.3) GO TO 60                                          ADS69050
      CALL ADS504 (NDV,H,Y,S)                                           ADS69060
      CALL ADS501 (NDV,P,Y,SIGMA)                                       ADS69070
      CALL ADS501 (NDV,Y,S,TAU)                                         ADS69080
      SIGMB=1./SIGMA                                                    ADS69090
      TBU=1./TAU                                                        ADS69100
      IF (ICONT.EQ.1) GO TO 20                                          ADS69110
      IF (ICONT.EQ.2) GO TO 40                                          ADS69120
10    CONTINUE                                                          ADS69130
C     ------------------------------------------------------------------ADS69140
C     STEEPEST DESCENT FORMULATION.                                     ADS69150
C     ------------------------------------------------------------------ADS69160
      CALL ADS506 (NDV,H)                                               ADS69170
      CALL ADS504 (NDV,H,DF,S)                                          ADS69180
      GO TO 90                                                          ADS69190
20    CONTINUE                                                          ADS69200
C     ------------------------------------------------------------------ADS69210
C     D.F.P. FORMULATION.                                               ADS69220
C     ------------------------------------------------------------------ADS69230
      NA=0                                                              ADS69240
      DO 30 J=1,NDV                                                     ADS69250
      DO 30 I=1,J                                                       ADS69260
      NA=NA+1                                                           ADS69270
      H(NA)=H(NA)+SIGMB*P(I)*P(J)-TBU*S(I)*S(J)                         ADS69280
30    CONTINUE                                                          ADS69290
      CALL ADS504 (NDV,H,DF,S)                                          ADS69300
      GO TO 90                                                          ADS69310
40    CONTINUE                                                          ADS69320
C     ------------------------------------------------------------------ADS69330
C     B.F.G.S. FORMULATION.                                             ADS69340
C     ------------------------------------------------------------------ADS69350
      BB=(1.+TAU*SIGMB)*SIGMB                                           ADS69360
      NA=0                                                              ADS69370
      DO 50 J=1,NDV                                                     ADS69380
      DO 50 I=1,J                                                       ADS69390
      NA=NA+1                                                           ADS69400
50    H(NA)=H(NA)+BB*P(I)*P(J)-SIGMB*(S(I)*P(J)+S(J)*P(I))              ADS69410
      CALL ADS504 (NDV,H,DF,S)                                          ADS69420
      GO TO 90                                                          ADS69430
60    CONTINUE                                                          ADS69440
C     ------------------------------------------------------------------ADS69450
C     FLETCHER-REEVES FORMULATION.                                      ADS69460
C     ------------------------------------------------------------------ADS69470
      DFDF=0.                                                           ADS69480
      DO 70 I=1,NDV                                                     ADS69490
      DFD=DF(I)*DF(I)                                                   ADS69500
      IF (ISIDC(I).GE.1.AND.DF(I).GE.0.) DFD=0.                         ADS69510
      IF (ISIDC(I).LE.-1.AND.DF(I).LE.0.) DFD=0.                        ADS69520
70    DFDF=DFDF+DFD                                                     ADS69530
C     FLETCHER-REEVES MULTIPLIER.                                       ADS69540
      BETA=SCAL*DFDF/DF2                                                ADS69550
C     UPDATE DF2.                                                       ADS69560
      DF2=DFDF                                                          ADS69570
      IF (DF2.LE.ZRO) RETURN                                            ADS69580
C     UPDATE S.                                                         ADS69590
      DO 80 I=1,NDV                                                     ADS69600
      S(I)=DF(I)-BETA*S(I)                                              ADS69610
80    CONTINUE                                                          ADS69620
90    CONTINUE                                                          ADS69630
C     ------------------------------------------------------------------ADS69640
C     CALCULATION OF DIRECTION ; S.                                     ADS69650
C     ------------------------------------------------------------------ADS69660
C     CALCULATE S**2 AND MAX ABS(S(I))                                  ADS69670
      S2=0.                                                             ADS69680
      SCAL=ZRO                                                          ADS69690
      DO 100 I=1,NDV                                                    ADS69700
      S(I)=-S(I)                                                        ADS69710
      IF (ISIDC(I).EQ.-1.AND.S(I).LT.0.) S(I)=0.                        ADS69720
      IF (ISIDC(I).EQ.1.AND.S(I).GT.0.) S(I)=0.                         ADS69730
      IF (ISIDC(I).EQ.2) S(I)=0.                                        ADS69740
      SI=ABS(S(I))                                                      ADS69750
      S2=S2+SI*SI                                                       ADS69760
      IF (SI.GT.SCAL) SCAL=SI                                           ADS69770
100   CONTINUE                                                          ADS69780
      IF (ICONT.NE.0) GO TO 120                                         ADS69790
      DF2=0.                                                            ADS69800
      DO 110 I=1,NDV                                                    ADS69810
110   DF2=DF2+S(I)*S(I)                                                 ADS69820
120   CONTINUE                                                          ADS69830
C     ------------------------------------------------------------------ADS69840
C     NORMALIZATION.                                                    ADS69850
C     ------------------------------------------------------------------ADS69860
      SD=1./SCAL                                                        ADS69870
      DO 130 I=1,NDV                                                    ADS69880
130   S(I)=SD*S(I)                                                      ADS69890
      S2=S2*SD*SD                                                       ADS69900
      RETURN                                                            ADS69910
      END                                                               ADS69920
C     SUBROUTINE ADS237 (ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,            ADS69930
C    *G,PENLT,ZRO,PENALT)                                               ADS69940
C --------------------------------------------------------------------- ADS69950
C     ROUTINE TO CREATE PENALTY FOR EQUALITY CONSTRAINTS OR ADD PENALTY ADS69960
C     TO OBJECTIVE.  USED IN FEASIBLE DIRECTIONS METHODS.               ADS69970
C --------------------------------------------------------------------- ADS69980
C     BY G. N. VANDERPLAATS                              MAY 1985       ADS69990
C     UCSB                                                              ADS70000
C                                                                       ADS70010
C --- INPUT                                                             ADS70020
C     ICALC - CALCULATION CONTROL                                       ADS70030
C         1 - CREATE PENALTY MULTIPLIER FOR EACH EQUALITY CONSTRAINT.   ADS70040
C         2 - ADD PENALTY TO OBJ AND MULTIPLY G BY -1 IF NEEDED.        ADS70050
C         3 - MULTIPLY GRADIENT OF CONSTRIANT BY -1 IF NEEDED.          ADS70060
C             MODIFY DF TO INCLUDE PENALTY.                             ADS70070
C                                                                       ADS70080
      SUBROUTINE ADS237 (ICALC,NDV,PMLT,DF,A,NRA,IC,NGT,IDG,G,          ADS70090
     *PENLT,ZRO,PENALT)                                                 ADS70100
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION DF(1),A(NRA,1),IC(1),IDG(1),G(1),PENLT(1)               ADS70110
      IF(ICALC.GT.1) GO TO 30                                           ADS70120
C --------------------------------------------------------------------  ADS70130
C     1: CREATE PENALTY MULTIPLIERS.                                    ADS70140
C --------------------------------------------------------------------  ADS70150
      DO 10 J=1,NGT                                                     ADS70160
      K=IC(J)                                                           ADS70170
      IF(IDG(K).GE.0) GO TO 10                                          ADS70180
C     THIS IS AN EQUALITY CONSTRIANT.                                   ADS70190
C     NORM OF DH.                                                       ADS70200
      N1=1                                                              ADS70210
      CALL ADS501 (NDV,A(N1,J),A(N1,J),H2)                              ADS70220
      HNORM=SQRT(H2)                                                    ADS70230
C     NORM OF DF.                                                       ADS70240
      CALL ADS501 (NDV,DF,DF,DF2)                                       ADS70250
      DFNORM=SQRT(DF2)                                                  ADS70260
C     DH DOT DF.                                                        ADS70270
      CALL ADS501 (NDV,A(N1,J),DF,HDF)                                  ADS70280
      ZK=1.                                                             ADS70290
      IF(ABS(HDF).GT.ZRO) ZK=-HDF/ABS(HDF)                              ADS70300
      PENLT(K)=PMLT*DFNORM*ZK/HNORM                                     ADS70310
      IF(ABS(PENLT(K)).LT.ZRO) PENLT(K)=ZRO                             ADS70320
10    CONTINUE                                                          ADS70330
      RETURN                                                            ADS70340
30    IF(ICALC.EQ.3) GO TO 60                                           ADS70350
C --------------------------------------------------------------------  ADS70360
C     2 OR 4: CREATE PENALTY AND MULT. G BY -1 IF NEEDED.               ADS70370
C --------------------------------------------------------------------  ADS70380
      PENALT=0.                                                         ADS70390
      DO 40 J=1,NGT                                                     ADS70400
      K=IC(J)                                                           ADS70410
      IF(IDG(K).GE.0) GO TO 40                                          ADS70420
C     THIS IS AN EQUALITY CONSTRIANT.                                   ADS70430
      AMULT=PENLT(K)                                                    ADS70440
      GK=G(K)                                                           ADS70450
      IF(AMULT.LT.0.) G(K)=-GK                                          ADS70460
      IF(ICALC.EQ.4.AND.AMULT.LT.0.) GK=-GK                             ADS70470
      PENALT=PENALT-AMULT*GK                                            ADS70480
40    CONTINUE                                                          ADS70490
      RETURN                                                            ADS70500
C --------------------------------------------------------------------  ADS70510
C     3: MULTIPLY DH BY -1 IF NEEDED.  MODIFY DF TO INCLUDE PENALTY.    ADS70520
C --------------------------------------------------------------------  ADS70530
60    CONTINUE                                                          ADS70540
      DO 70 J=1,NGT                                                     ADS70550
      K=IC(J)                                                           ADS70560
      IF(IDG(K).GE.0) GO TO 70                                          ADS70570
C     EQUALITY CONSTRAINT.                                              ADS70580
C     MODIFY DF.                                                        ADS70590
      FCT=PENLT(K)                                                      ADS70600
      DO 62 I=1,NDV                                                     ADS70610
62    DF(I)=DF(I)-FCT*A(I,J)                                            ADS70620
      IF(FCT.GT.0.) GO TO 70                                            ADS70630
C     MULTIPLY DH BY -1.                                                ADS70640
      DO 65 I=1,NDV                                                     ADS70650
65    A(I,J)=-A(I,J)                                                    ADS70660
70    CONTINUE                                                          ADS70670
      RETURN                                                            ADS70680
      END                                                               ADS70690
C --- SUBROUTINE ADS301 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,          ADS70700
C    1SVLB,SVUB,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX, ADS70710
C    2RP,RMULT,XTMP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,            ADS70720
C    3SCFX,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL)                     ADS70730
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),IDG(NCON),IC(NCOLA), ADS70740
C    1DF(NDV),A(NRA,NCOLA),WK(NRWK),IWK(NRIWK),XTMP(NDV),SCFX(NDV),     ADS70750
C    2ALAMDA(NCON),PPARAM(NCON)                                         ADS70760
C     ------------------------------------------------------------------ADS70770
C --- ISTRAT = 1.  PERFORM SEQUENTIAL UNCONSTRAINED MINIMIZATION USING  ADS70780
C     THE EXTERIOR PENALTY FUNCTION METHOD.                             ADS70790
C     ------------------------------------------------------------------ADS70800
C     BY G. N. VANDERPLAATS                           JULY 1982.        ADS70810
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS70820
C --- INPUT.                                                            ADS70830
C      IGOTO  - CONTROL PARAMETER.                                      ADS70840
C      INFO   - CONTROL PARAMETER.                                      ADS70850
C      IOPT   - OPTIMIZATION IDENTIFIER.                                ADS70860
C      IONED  - ONE D. SEARCH IDENTIFIER.                               ADS70870
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS70880
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS70890
C      OBJ    - OBJECTIVE.                                              ADS70900
C      X      - CURRENT VALUE OF DESIGN VARIABLES.                      ADS70910
C      SVLB,SVUB - LOWER AND UPPER BOUND ON DESIGN VARIABLES.           ADS70920
C      G      - CONSTRAINT VALUES.                                      ADS70930
C      IDG    - CONSTRAINT TYPE IDENYIFIER.                             ADS70940
C      NGT   - NUMBER OF ACTIVE CONSTRAINTS.                            ADS70950
C      IC     - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                   ADS70960
C      DF     - GRADIENT OF OBKECTIVE.                                  ADS70970
C      A      - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS.       ADS70980
C      NRA    - NUMBER OF ROWS IN ARRAY A.                              ADS70990
C      NCOLA  - NUMBER OF COLUMNS IN ARRAY A.                           ADS71000
C      WK,IWK - REAL,INTEGER WORK ARRAY.                                ADS71010
C      IPRNT3 - PRINT CONTROL IDENTIFIER FOR STRATEGY LEVEL.            ADS71020
C      CT,CTL - PARAMETERS TO IDENTIFY ACTIVE CONSTRAINTS.              ADS71030
C      ITRMST,CTMIN,CTLMIN,DELSTR,DABSTR - PARAMETERS FOR CONVERGENCE   ADS71040
C               CRITERIA.                                               ADS71050
C      JTMAX  - MAXIMUM NUMBER OF ITERATIONS.                           ADS71060
C      RP     - PENALTY PARAMETER FOR EXTERIOR PENALTY FUNCTION.        ADS71070
C      RMULT  - FACTOR FOR RP.                                          ADS71080
C      XTMP   - TEMPORARY ARRAY FOR X.                                  ADS71090
C      RPMAX  - MAXIMUM VALUE OF RP.                                    ADS71100
C      RPPMIN - MINIMUM VALUE OF RPPRIM.                                ADS71110
C      EPSPEN - TRANSITION POINT.                                       ADS71120
C      RPPRIM - PENALTY PARAMETER FOR INTERIOR PENALTY FUNCTION.        ADS71130
C      RPMULT - FACTOR FOR RPPRIM.                                      ADS71140
C      SCFX   - SCALE FACTOR FOR X.                                     ADS71150
C      ALAMDA - LAGRANGE MULTIPLIER.                                    ADS71160
C      PPARAM - PENALTY PARAMETER OF EACH CONSTRAINT.                   ADS71170
C      SCFO   - SCALE FACTOR OF OBJECTIVE.                              ADS71180
C      IFCALL - NUMBER OF FUNCTION EVALUATION.                          ADS71190
C      IGCALL - NUMBER OF GRADIENT EVALUATION.                          ADS71200
C --- OUTPUT.                                                           ADS71210
C      OBJ    - OPTIMUM OBJECTIVE.                                      ADS71220
C      X      - OPTIMUM VALUE OF DESIGN VARIABLES.                      ADS71230
C      G      - CONSTRAINT VALUE CORRESPOND TO OPTIMUM DESIGN VARIABLES.ADS71240
C      ISFLAG - SCALE CONTROL PARAMETER. =1; ADDITIONAL SCALE FACTOR IS ADS71250
C               CALCULATED. =0; ONLY GRADIENTS ARE REQUIRED.            ADS71260
C                                                                       ADS71270
