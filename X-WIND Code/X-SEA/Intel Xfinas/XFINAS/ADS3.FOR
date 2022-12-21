      SUBROUTINE ADS301 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,ADS71280
     1G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,XTADS71290
     2MP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,SCFX,ALAMDA,PPARAM,SCFOADS71300
     3,ISFLAG,IFCALL,IGCALL,OBJZ,OBJZ1,PENALT,IDAB,IDEL,INFO2C,JTER,ITMPADS71310
     4,ISIDC,SCFG)                                                      ADS71320
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WKADS71330
     1(1),IWK(1),XTMP(1),SCFX(1),ALAMDA(1),PPARAM(1),ITMP(1),ISIDC(1),  ADS71340
     2SCFG(1)                                                           ADS71350
C                                                                       ADS71360
      ISTRAT=1                                                          ADS71370
      IGOTO1=IGOTO+1                                                    ADS71380
      GO TO (10,50,90), IGOTO1                                          ADS71390
10    CONTINUE                                                          ADS71400
      INFO=0                                                            ADS71410
      IF (IPRNT3.GT.0) WRITE (6,290)                                    ADS71420
C     ------------------------------------------------------------------ADS71430
C     FIRST CALL.  INITIALIZE PARAMETERS.                               ADS71440
C     ------------------------------------------------------------------ADS71450
      JTER=0                                                            ADS71460
      IDEL=0                                                            ADS71470
      IDAB=0                                                            ADS71480
      INFO2C=0                                                          ADS71490
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS71500
      IF (IPRNT3.GT.0) WRITE (6,300)                                    ADS71510
C     CALCULATE TRUE OBJECTIVE.                                         ADS71520
      OBJZ=OBJ                                                          ADS71530
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS71540
     1IC)                                                               ADS71550
      OBJ=OBJZ+PENALT                                                   ADS71560
      OBJW=OBJZ/SCFO                                                    ADS71570
      IF (IPRNT3.GT.0) WRITE (6,310) OBJ,PENALT,OBJW                    ADS71580
      ICALC=1                                                           ADS71590
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS71600
      ICALC=4                                                           ADS71610
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS71620
C     ------------------------------------------------------------------ADS71630
C     BEGIN NEW ITERATION.                                              ADS71640
C     ------------------------------------------------------------------ADS71650
20    JTER=JTER+1                                                       ADS71660
C     SAVE OBJZ.                                                        ADS71670
      OBJZ1=OBJZ                                                        ADS71680
      IF (IPRNT3.GT.0) WRITE (6,320) JTER,RP                            ADS71690
      ICALC=3                                                           ADS71700
      IF (IPRNT3.LT.5) GO TO 40                                         ADS71710
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS71720
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS71730
      ISCAL=0                                                           ADS71740
      DO 30 J=1,NGT                                                     ADS71750
      WRITE (6,150) IC(J)                                               ADS71760
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS71770
30    CONTINUE                                                          ADS71780
40    CONTINUE                                                          ADS71790
C     ------------------------------------------------------------------ADS71800
C     OPTIMIZE UNCONSTRAINED PSEUDO-OBJECTIVE FUNCTION.                 ADS71810
C     ------------------------------------------------------------------ADS71820
50    CONTINUE                                                          ADS71830
      IF (INFO.NE.2) GO TO 60                                           ADS71840
      CALL ADS008 (NDV,ISTRAT,NCON,NGT,A,NRA,IC,G,IDG,DF,ALAMDA,PPARAM,RADS71850
     1P,EPSPEN)                                                         ADS71860
60    CONTINUE                                                          ADS71870
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDGADS71880
     1,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                                    ADS71890
C     FINISHED?                                                         ADS71900
      IF (INFO.EQ.2) INFO2C=INFO2C+1                                    ADS71910
      IF (INFO2C.EQ.1.AND.INFO.EQ.2) GO TO 50                           ADS71920
C     NO. EVALUATE FUNCTION OR GRADIENT AND CALL AGAIN.                 ADS71930
      IGOTO=1                                                           ADS71940
      IF (INFO.EQ.1) RETURN                                             ADS71950
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS71960
     1IC)                                                               ADS71970
C     TRUE OBJECTIVE.                                                   ADS71980
      OBJZ=OBJ-PENALT                                                   ADS71990
      IF (INFO.EQ.0) GO TO 70                                           ADS72000
      OBJ=OBJZ                                                          ADS72010
      ISFLAG=0                                                          ADS72020
      RETURN                                                            ADS72030
70    CONTINUE                                                          ADS72040
C     ------------------------------------------------------------------ADS72050
C     UNCONSTRAINED MINIMIZATION IS COMPLETE. UPDATE OBJ AND RP.        ADS72060
C     ------------------------------------------------------------------ADS72070
      IF (IPRNT3.GT.3) WRITE (6,160)                                    ADS72080
      OBJW=OBJZ/SCFO                                                    ADS72090
      IF (IPRNT3.GT.0) WRITE (6,310) OBJ,PENALT,OBJW                    ADS72100
      ICALC=1                                                           ADS72110
      IF (IPRNT3.GE.2) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS72120
      ICALC=4                                                           ADS72130
      IF (IPRNT3.GE.3) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS72140
C     ------------------------------------------------------------------ADS72150
C     CHECK CONVERGENCE.                                                ADS72160
C     ------------------------------------------------------------------ADS72170
      CALL ADS211 (OBJZ,OBJZ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)        ADS72180
C     MAXIMUM CONSTRAINT VIOLATION.                                     ADS72190
      GMCT=-1.                                                          ADS72200
      DO 80 I=1,NCON                                                    ADS72210
      CC=CTMIN                                                          ADS72220
      IF (IDG(I).GT.1) CC=CTLMIN                                        ADS72230
      GC=G(I)-CC                                                        ADS72240
      IF (GC.GT.GMCT) GMCT=GC                                           ADS72250
80    CONTINUE                                                          ADS72260
      IF (GMCT.GT.0.) IDEL=0                                            ADS72270
      IF (IDEL.GE.ITRMST.OR.IDAB.GE.ITRMST) GO TO 100                   ADS72280
      RP=RP*RMULT                                                       ADS72290
      IF (RP.GT.RPMAX.AND.IDEL.GT.0) GO TO 100                          ADS72300
      IF (RP.GT.RPMAX) RP=RPMAX                                         ADS72310
      IF (JTER.GE.JTMAX) GO TO 100                                      ADS72320
C     EVALUATE GRADIENT AND CALL AGAIN.                                 ADS72330
      OBJ=OBJZ                                                          ADS72340
      INFO=2                                                            ADS72350
      ISFLAG=1                                                          ADS72360
      INFO2C=0                                                          ADS72370
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS72380
      IGOTO=2                                                           ADS72390
      RETURN                                                            ADS72400
90    CONTINUE                                                          ADS72410
C     NEW PENALIZED OBJECTIVE.                                          ADS72420
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS72430
     1IC)                                                               ADS72440
      OBJ=OBJZ+PENALT                                                   ADS72450
      GO TO 20                                                          ADS72460
100   CONTINUE                                                          ADS72470
C     ------------------------------------------------------------------ADS72480
C     OPTIMIZATION IS COMPLETE.                                         ADS72490
C     ------------------------------------------------------------------ADS72500
      INFO=0                                                            ADS72510
      IGOTO=0                                                           ADS72520
      OBJ=OBJZ                                                          ADS72530
      IF (IPRNT3.EQ.0) GO TO 140                                        ADS72540
      OBJW=OBJ/SCFO                                                     ADS72550
      WRITE (6,170) JTER,OBJW                                           ADS72560
C     X-VECTOR.                                                         ADS72570
      ICALC=1                                                           ADS72580
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS72590
C     G-VECTOR.                                                         ADS72600
      ICALC=4                                                           ADS72610
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS72620
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS72630
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS72640
C     SIDE CONSTRAINTS.                                                 ADS72650
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS72660
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS72670
      WRITE (6,180) CT,CTL                                              ADS72680
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS72690
      NAVC=NAC+NVC                                                      ADS72700
      WRITE (6,240) NAC,NVC                                             ADS72710
      IF (NAVC.GT.0) WRITE (6,250)                                      ADS72720
      IF (NAVC.GT.0) WRITE (6,260) (ITMP(I),I=1,NAVC)                   ADS72730
      WRITE (6,270) NACS                                                ADS72740
      IF (NACS.EQ.0) GO TO 130                                          ADS72750
      N=0                                                               ADS72760
      DO 120 I=1,NDV                                                    ADS72770
      IF (ISIDC(I).EQ.0) GO TO 120                                      ADS72780
      IF (ISIDC(I).EQ.1) GO TO 110                                      ADS72790
C     LOWER BOUND IS ACTIVE.                                            ADS72800
      N=N+1                                                             ADS72810
      ITMP(N)=-I                                                        ADS72820
110   IF (ISIDC(I).LT.0) GO TO 120                                      ADS72830
C     UPPER BOUND IS ACTIVE.                                            ADS72840
      N=N+1                                                             ADS72850
      ITMP(N)=I                                                         ADS72860
120   CONTINUE                                                          ADS72870
      WRITE (6,280)                                                     ADS72880
      WRITE (6,260) (ITMP(I),I=1,NACS)                                  ADS72890
130   WRITE (6,190)                                                     ADS72900
      IF (JTER.GE.JTMAX) WRITE (6,200)                                  ADS72910
      IF (IDEL.GE.ITRMST) WRITE (6,210) ITRMST                          ADS72920
      IF (IDAB.GE.ITRMST) WRITE (6,220) ITRMST                          ADS72930
      IF (RP.GE.RPMAX.AND.IDEL.GT.0) WRITE (6,230) RPMAX                ADS72940
140   CONTINUE                                                          ADS72950
      RETURN                                                            ADS72960
C     ------------------------------------------------------------------ADS72970
C     FORMATS.                                                          ADS72980
C     ------------------------------------------------------------------ADS72990
150   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS73000
160   FORMAT (/,5X,'UNCONSTRAINED MINIMIZATIN IS COMPLETE')             ADS73010
170   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS73020
     1'NUMBER OF UNCONSTRAINED MINIMIZATIONS =',I5,//,                  ADS73030
     25X,11HOBJECTIVE =,1PE12.5)                                        ADS73040
180   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,'  CTL =',      ADS73050
     11PE12.5)                                                          ADS73060
190   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS73070
200   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS73080
210   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS73090
     1' CONSECUTIVE ITERATIONS')                                        ADS73100
220   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS73110
     1' CONSECUTIVE ITERATIONS')                                        ADS73120
230   FORMAT (/,5X,'PENALTY PERAMETER EXCEEDS',1PE12.5)                 ADS73130
240   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS73140
     1' VIOLATED CONSTRAINTS')                                          ADS73150
250   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS73160
260   FORMAT (5X,10I5)                                                  ADS73170
270   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS73180
280   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS73190
290   FORMAT ('1',4X,75('-'),/,20X,                                     ADS73200
     1'ISTRAT = 1; EXTERIOR PENALTY FUNCTION METHOD',/,5X,75('-'))      ADS73210
300   FORMAT (/,5X,'-- INITIAL DESIGN')                                 ADS73220
310   FORMAT (/,5X,'SCALED OBJ =',1PE12.5,5X,'PENALT =',1PE12.5,/,5X,   ADS73230
     1'TRUE OBJ   =',1PE12.5)                                           ADS73240
320   FORMAT (///,5X,'-- BEGIN SUMT ITERATION',I5,5X,'RP =',1PE12.5)    ADS73250
      END                                                               ADS73260
C --- SUBROUTINE ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAADS73270
C    1M,NGT,IC)                                                         ADS73280
C     DIMENSION G(NCON),IDG(NCON),ALAMDA(NCON),PPARAM(NCON),IC(NCOLA)   ADS73290
C     ------------------------------------------------------------------ADS73300
C --- ROUTINE TO CALCULATE THE PENALTY IN PENALTY FUNCTION METHODS.     ADS73310
C     ------------------------------------- ----------------------------ADS73320
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 FEB. 1982   ADS73330
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS73340
C --- INPUT                                                             ADS73350
C     ISTRAT - PENALTY FUNCTION METHOD.                                 ADS73360
C         0 - NONE. SET PENALT=0 AND RETURN.                            ADS73370
C         1 - EXTERIOR PENALTY.                                         ADS73380
C         2 - LINEAR EXTENDED INTERIOR PENALTY.                         ADS73390
C         3 - QUADRATIC EXTENDED INTERIOR PENALTY.                      ADS73400
C         4 - QUBIC EXTENDED INTERIOR PENALTY.                          ADS73410
C         5 - AUGMENTED LAGRANGEAN MULTIPLIER.                          ADS73420
C     .GT.5 - NONE.                                                     ADS73430
C     NCON   - NUMBER OF CONSTRAINTS.                                   ADS73440
C     G      - CONSTRAINTS.                                             ADS73450
C     IDG    - CONSTRAINT TYPE IDENTIFIER.                              ADS73460
C     EPSPEN - TRANSITION POINT.                                        ADS73470
C     RP     - PENALTY PARAMETER FOR EXTERIOR PENALTY FUNCTION METHOD.  ADS73480
C     ALAMDA - LAGRANGE MULTIPLIER.                                     ADS73490
C     PPARAM - PENALTY PARAMETER.                                       ADS73500
C     NGT   - NUMBER OF ACTIVE CONSTRAINTS.                             ADS73510
C     IC     - ARRAY IDENTIFYING NGT ACTIVE CONSTRAINTS.                ADS73520
C --- OUTPUT                                                            ADS73530
C     PENALT - PENALTY.                                                 ADS73540
C                                                                       ADS73550
      SUBROUTINE ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAADS73560
     1M,NGT,IC)                                                         ADS73570
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IDG(1),ALAMDA(1),PPARAM(1),IC(1)                   ADS73580
C                                                                       ADS73590
      PENALT=0.                                                         ADS73600
      IF (ISTRAT.LE.0.OR.ISTRAT.GT.5) RETURN                            ADS73610
      IF (NCON.LE.0) RETURN                                             ADS73620
      GO TO (10,30,80,130,190), ISTRAT                                  ADS73630
C                                                                       ADS73640
10    CONTINUE                                                          ADS73650
C     ------------------------------------------------------------------ADS73660
C     ISTRAT = 1.  EXTERIOR PENALTY FUNCTION.                           ADS73670
C     ------------------------------------------------------------------ADS73680
      DO 20 K=1,NGT                                                     ADS73690
      I=IC(K)                                                           ADS73700
      IF (IDG(I).GE.0.AND.G(I).LE.0.) GO TO 20                          ADS73710
C     EQUALITY CONSTRAINT OR VIOLATED INEQUALITY CONSTRAINT.            ADS73720
      PENALT=PENALT+G(I)*G(I)                                           ADS73730
20    CONTINUE                                                          ADS73740
      PENALT=RP*PENALT                                                  ADS73750
      GO TO 230                                                         ADS73760
C                                                                       ADS73770
30    CONTINUE                                                          ADS73780
C     ------------------------------------------------------------------ADS73790
C     ISTRAT = 2.  LINEAR EXTENDED PENALTY.                             ADS73800
C     ------------------------------------------------------------------ADS73810
      DO 70 K=1,NGT                                                     ADS73820
      I=IC(K)                                                           ADS73830
      IF (IDG(I).GE.0) GO TO 40                                         ADS73840
      A=G(I)*G(I)                                                       ADS73850
      GO TO 60                                                          ADS73860
40    CONTINUE                                                          ADS73870
      IF (G(I).GT.EPSPEN) GO TO 50                                      ADS73880
C     G(I).LE.EPSPEN. USE INTERIOR PENALTY.                             ADS73890
      A=-G(I)                                                           ADS73900
      A=1./A                                                            ADS73910
      GO TO 60                                                          ADS73920
50    CONTINUE                                                          ADS73930
C     G(I).GT.EPSPEN. USE EXTENDED PENALTY.                             ADS73940
      A=G(I)/EPSPEN-1.                                                  ADS73950
      A=(-A+1.)/(-EPSPEN)                                               ADS73960
60    PENALT=PENALT+A*PPARAM(I)                                         ADS73970
70    CONTINUE                                                          ADS73980
      GO TO 230                                                         ADS73990
C                                                                       ADS74000
80    CONTINUE                                                          ADS74010
C     ------------------------------------------------------------------ADS74020
C     ISTRAT = 3.  QUADRATIC EXTENDED INTERIOR PENALTY.                 ADS74030
C     ------------------------------------------------------------------ADS74040
      DO 120 K=1,NGT                                                    ADS74050
      I=IC(K)                                                           ADS74060
      IF (IDG(I).GE.0) GO TO 90                                         ADS74070
      A=G(I)*G(I)                                                       ADS74080
      GO TO 110                                                         ADS74090
90    CONTINUE                                                          ADS74100
      IF (G(I).GT.EPSPEN) GO TO 100                                     ADS74110
C     G(I).LE.EPSPEN. USE INTERIOR PENALTY.                             ADS74120
      A=-G(I)                                                           ADS74130
      A=1./A                                                            ADS74140
      GO TO 110                                                         ADS74150
100   CONTINUE                                                          ADS74160
C     G(I).GT.EPSPEN. USE EXTENDED PENALTY.                             ADS74170
      A=G(I)/EPSPEN-1.                                                  ADS74180
      A=(A**2-A+1.)/(-EPSPEN)                                           ADS74190
110   PENALT=PENALT+A*PPARAM(I)                                         ADS74200
120   CONTINUE                                                          ADS74210
      GO TO 230                                                         ADS74220
C                                                                       ADS74230
130   CONTINUE                                                          ADS74240
C     ------------------------------------------------------------------ADS74250
C     ISTRAT = 4.  VARIABLE PENALTY.                                    ADS74260
C     ------------------------------------------------------------------ADS74270
      A=-100.                                                           ADS74280
      DO 140 I=1,NCON                                                   ADS74290
      IF (G(I).GT.A) A=G(I)                                             ADS74300
140   CONTINUE                                                          ADS74310
      A=A/EPSPEN                                                        ADS74320
      SHAPE=2./3.                                                       ADS74330
      IF (A.GT.SHAPE) A=SHAPE                                           ADS74340
      SHAPE=1./(3.*(1.-A))                                              ADS74350
      DO 180 K=1,NGT                                                    ADS74360
      I=IC(K)                                                           ADS74370
      IF (IDG(I).GE.0) GO TO 150                                        ADS74380
      A=G(I)*G(I)                                                       ADS74390
      GO TO 170                                                         ADS74400
150   CONTINUE                                                          ADS74410
      IF (G(I).GT.EPSPEN) GO TO 160                                     ADS74420
C     G(I).LE.EPSPEN. USE INTERIOR PENALTY.                             ADS74430
      A=-G(I)                                                           ADS74440
      A=1./A                                                            ADS74450
      GO TO 170                                                         ADS74460
160   CONTINUE                                                          ADS74470
C     G(I).GT.EPSPEN. USE EXTENDED PENALTY.                             ADS74480
      A=G(I)/EPSPEN-1.                                                  ADS74490
      A=(SHAPE*A**3+A**2-A+1.)/(-EPSPEN)                                ADS74500
170   PENALT=PENALT+A*PPARAM(I)                                         ADS74510
180   CONTINUE                                                          ADS74520
      GO TO 230                                                         ADS74530
C                                                                       ADS74540
190   CONTINUE                                                          ADS74550
C     ------------------------------------------------------------------ADS74560
C     ISTRAT = 5.  AUGMENTED LAGRANGEAN MULTIPLIER.                     ADS74570
C     ------------------------------------------------------------------ADS74580
      DO 220 K=1,NGT                                                    ADS74590
      I=IC(K)                                                           ADS74600
      IF (IDG(I).LT.0) GO TO 200                                        ADS74610
      V=ALAMDA(I)+PPARAM(I)*G(I)                                        ADS74620
      IF (V.LE.0.) V=0.                                                 ADS74630
      PENAD=(V*V-ALAMDA(I)**2)/(2.*PPARAM(I))                           ADS74640
      IF (G(I).LE.-0.1) PENAD=0.                                        ADS74650
      GO TO 210                                                         ADS74660
200   CONTINUE                                                          ADS74670
      PENAD=ALAMDA(I)*G(I)+0.5*PPARAM(I)*G(I)*G(I)                      ADS74680
210   CONTINUE                                                          ADS74690
      PENALT=PENALT+PENAD                                               ADS74700
220   CONTINUE                                                          ADS74710
C                                                                       ADS74720
230   CONTINUE                                                          ADS74730
C     ------------------------------------------------------------------ADS74740
C     ADD PENALTY, IF INITIALLY NON-ACTIVE CONSTRAINTS ARE VIOLATED.    ADS74750
C     ------------------------------------------------------------------ADS74760
      DO 250 I=1,NCON                                                   ADS74770
      DO 240 K=1,NGT                                                    ADS74780
      IF (I.EQ.IC(K)) GO TO 250                                         ADS74790
240   CONTINUE                                                          ADS74800
      IF (G(I).LE.0.) GO TO 250                                         ADS74810
      PENALT=PENALT+RP*G(I)*G(I)                                        ADS74820
250   CONTINUE                                                          ADS74830
      RETURN                                                            ADS74840
      END                                                               ADS74850
C --- SUBROUTINE ADS303 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,          ADS74860
C    1SVLB,SVUB,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX, ADS74870
C    2RP,RMULT,XTMP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,     ADS74880
C    3EPSPEN,RPPRIM,RPMULT,SCFX,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL)ADS74890
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),IDG(NCON),IC(NCOLA), ADS74900
C    1DF(NDV),A(NRA,NCOLA),WK(NRWK),IWK(NRIWK),XTMP(NDV),SCFX(NDV),     ADS74910
C    2ALAMDA(NCON),PPARAM(NCON)                                         ADS74920
C     ------------------------------------------------------------------ADS74930
C --- ISTRAT = 2.  PERFORM SEQUENTIAL UNCONSTRAINED MINIMIZATION USING  ADS74940
C     THE LINEAR EXTENDED INTERIOR PENALTY FUNCTION METHOD.             ADS74950
C     ------------------------------------------------------------------ADS74960
C     BY G. N. VANDERPLAATS                           JULY 1982.        ADS74970
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS74980
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS74990
C --- INPUT.                                                            ADS75000
C      IGOTO  - CONTROL PARAMETER.                                      ADS75010
C      INFO   - CONTROL PARAMETER.                                      ADS75020
C      IOPT   - OPTIMIZATION IDENTIFIER.                                ADS75030
C      IONED  - ONE D. SEARCH IDENTIFIER.                               ADS75040
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS75050
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS75060
C      OBJ    - OBJECTIVE.                                              ADS75070
C      X      - CURRENT VALUE OF DESIGN VARIABLES.                      ADS75080
C      SVLB,SVUB - LOWER AND UPPER BOUND ON DESIGN VARIABLES.           ADS75090
C      G      - CONSTRAINT VALUES.                                      ADS75100
C      IDG    - CONSTRAINT TYPE IDENYIFIER.                             ADS75110
C      NGT   - NUMBER OF ACTIVE CONSTRAINTS.                            ADS75120
C      IC     - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                   ADS75130
C      DF     - GRADIENT OF OBKECTIVE.                                  ADS75140
C      A      - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS.       ADS75150
C      NRA    - NUMBER OF ROWS IN ARRAY A.                              ADS75160
C      NCOLA  - NUMBER OF COLUMNS IN ARRAY A.                           ADS75170
C      WK,IWK - REAL,INTEGER WORK ARRAY.                                ADS75180
C      IPRNT3 - PRINT CONTROL IDENTIFIER FOR STRATEGY LEVEL.            ADS75190
C      CT,CTL - PARAMETERS TO IDENTIFY ACTIVE CONSTRAINTS.              ADS75200
C      ITRMST,CTMIN,CTLMIN,DELSTR,DABSTR - PARAMETERS FOR CONVERGENCE   ADS75210
C               CRITERIA.                                               ADS75220
C      JTMAX  - MAXIMUM NUMBER OF ITERATIONS.                           ADS75230
C      RP     - PENALTY PARAMETER FOR EXTERIOR PENALTY FUNCTION.        ADS75240
C      RMULT  - FACTOR FOR RP.                                          ADS75250
C      XTMP   - TEMPORARY ARRAY FOR X.                                  ADS75260
C      RPMAX  - MAXIMUM VALUE OF RP.                                    ADS75270
C      RPPMIN - MINIMUM VALUE OF RPPRIM.                                ADS75280
C      EPSPEN - TRANSITION POINT.                                       ADS75290
C      RPPRIM - PENALTY PARAMETER FOR INTERIOR PENALTY FUNCTION.        ADS75300
C      RPMULT - FACTOR FOR RPPRIM.                                      ADS75310
C      SCFX   - SCALE FACTOR FOR X.                                     ADS75320
C      ALAMDA - LAGRANGE MULTIPLIER.                                    ADS75330
C      PPARAM - PENALTY PARAMETER OF EACH CONSTRAINT.                   ADS75340
C      SCFO   - SCALE FACTOR OF OBJECTIVE.                              ADS75350
C      IFCALL - NUMBER OF FUNCTION EVALUATION.                          ADS75360
C      IGCALL - NUMBER OF GRADIENT EVALUATION.                          ADS75370
C --- OUTPUT.                                                           ADS75380
C      OBJ    - OPTIMUM OBJECTIVE.                                      ADS75390
C      X      - OPTIMUM VALUE OF DESIGN VARIABLES.                      ADS75400
C      G      - CONSTRAINT VALUE CORRESPOND TO OPTIMUM DESIGN VARIABLES.ADS75410
C      ISFLAG - SCALE CONTROL PARAMETER. =1; ADDITIONAL SCALE FACTOR IS ADS75420
C               CALCULATED. =0; ONLY GRADIENTS ARE REQUIRED.            ADS75430
C                                                                       ADS75440
      SUBROUTINE ADS303 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,ADS75450
     1G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,XTADS75460
     2MP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPADS75470
     3MULT,SCFX,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL,OBJZ,OBJZ1,PENALADS75480
     4T,IDAB,IDEL,INFO2C,JTER,ITMP,ISIDC,SCFG)                          ADS75490
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WKADS75500
     1(1),IWK(1),XTMP(1),SCFX(1),ALAMDA(1),PPARAM(1),ITMP(1),ISIDC(1),  ADS75510
     2SCFG(1)                                                           ADS75520
C                                                                       ADS75530
      ISTRAT=2                                                          ADS75540
      IGOTO1=IGOTO+1                                                    ADS75550
      GO TO (10,50,90), IGOTO1                                          ADS75560
10    CONTINUE                                                          ADS75570
      INFO=0                                                            ADS75580
      IF (IPRNT3.GT.0) WRITE (6,290)                                    ADS75590
C     ------------------------------------------------------------------ADS75600
C     FIRST CALL.  INITIALIZE PARAMETERS.                               ADS75610
C     ------------------------------------------------------------------ADS75620
      JTER=0                                                            ADS75630
      IDEL=0                                                            ADS75640
      IDAB=0                                                            ADS75650
      INFO2C=0                                                          ADS75660
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS75670
      IF (IPRNT3.GT.0) WRITE (6,300)                                    ADS75680
C     CALCULATE TRUE OBJECTIVE.                                         ADS75690
      CALL ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                           ADS75700
      OBJZ=OBJ                                                          ADS75710
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS75720
     1IC)                                                               ADS75730
      OBJ=OBJZ+PENALT                                                   ADS75740
      OBJW=OBJZ/SCFO                                                    ADS75750
      IF (IPRNT3.GT.0) WRITE (6,310) OBJ,PENALT,OBJW                    ADS75760
      ICALC=1                                                           ADS75770
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS75780
      ICALC=4                                                           ADS75790
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS75800
C     ------------------------------------------------------------------ADS75810
C     BEGIN NEW ITERATION.                                              ADS75820
C     ------------------------------------------------------------------ADS75830
20    JTER=JTER+1                                                       ADS75840
C     SAVE OBJZ.                                                        ADS75850
      OBJZ1=OBJZ                                                        ADS75860
      IF (IPRNT3.GT.0) WRITE (6,320) JTER,RPPRIM,EPSPEN                 ADS75870
      ICALC=3                                                           ADS75880
      IF (IPRNT3.LT.5) GO TO 40                                         ADS75890
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS75900
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS75910
      ISCAL=0                                                           ADS75920
      DO 30 J=1,NGT                                                     ADS75930
      WRITE (6,150) IC(J)                                               ADS75940
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS75950
30    CONTINUE                                                          ADS75960
40    CONTINUE                                                          ADS75970
C     ------------------------------------------------------------------ADS75980
C     OPTIMIZE UNCONSTRAINED PSEUDO-OBJECTIVE FUNCTION.                 ADS75990
C     ------------------------------------------------------------------ADS76000
50    CONTINUE                                                          ADS76010
      IF (INFO.NE.2) GO TO 60                                           ADS76020
      CALL ADS008 (NDV,ISTRAT,NCON,NGT,A,NRA,IC,G,IDG,DF,ALAMDA,PPARAM,RADS76030
     1P,EPSPEN)                                                         ADS76040
60    CONTINUE                                                          ADS76050
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDGADS76060
     1,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                                    ADS76070
C     FINISHED?                                                         ADS76080
      IF (INFO.EQ.2) INFO2C=INFO2C+1                                    ADS76090
      IF (INFO2C.EQ.1.AND.INFO.EQ.2) GO TO 50                           ADS76100
C     NO. EVALUATE FUNCTION OR GRADIENT AND CALL AGAIN.                 ADS76110
      IGOTO=1                                                           ADS76120
      IF (INFO.EQ.1) RETURN                                             ADS76130
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS76140
     1IC)                                                               ADS76150
C     TRUE OBJECTIVE.                                                   ADS76160
      OBJZ=OBJ-PENALT                                                   ADS76170
      IF (INFO.EQ.0) GO TO 70                                           ADS76180
      OBJ=OBJZ                                                          ADS76190
      ISFLAG=0                                                          ADS76200
      RETURN                                                            ADS76210
70    CONTINUE                                                          ADS76220
C     ------------------------------------------------------------------ADS76230
C     UNCONSTRAINED MINIMIZATION IS COMPLETE. UPDATE OBJ AND RP.        ADS76240
C     ------------------------------------------------------------------ADS76250
      IF (IPRNT3.GT.3) WRITE (6,160)                                    ADS76260
      OBJW=OBJZ/SCFO                                                    ADS76270
      IF (IPRNT3.GT.0) WRITE (6,310) OBJ,PENALT,OBJW                    ADS76280
      ICALC=1                                                           ADS76290
      IF (IPRNT3.GE.2) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS76300
      ICALC=4                                                           ADS76310
      IF (IPRNT3.GE.3) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS76320
C     ------------------------------------------------------------------ADS76330
C     CHECK CONVERGENCE.                                                ADS76340
C     ------------------------------------------------------------------ADS76350
      CALL ADS211 (OBJZ,OBJZ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)        ADS76360
C     MAXIMUM CONSTRAINT VIOLATION.                                     ADS76370
      GMCT=-1.                                                          ADS76380
      DO 80 I=1,NCON                                                    ADS76390
      CC=CTMIN                                                          ADS76400
      IF (IDG(I).GT.1) CC=CTLMIN                                        ADS76410
      GC=G(I)-CC                                                        ADS76420
      IF (GC.GT.GMCT) GMCT=GC                                           ADS76430
80    CONTINUE                                                          ADS76440
      IF (GMCT.GT.0.) IDEL=0                                            ADS76450
      IF (IDEL.GE.ITRMST.OR.IDAB.GE.ITRMST) GO TO 100                   ADS76460
      RPPRIM=RPMULT*RPPRIM                                              ADS76470
      IF (RPPRIM.LT.RPPMIN.AND.IDEL.GT.0) GO TO 100                     ADS76480
      IF (RPPRIM.LT.RPPMIN) RPPRIM=RPPMIN                               ADS76490
      RP=RP*RMULT                                                       ADS76500
      IF (RP.GT.RPMAX) RP=RPMAX                                         ADS76510
      CALL ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                           ADS76520
      EPSPEN=EPSPEN*(RPMULT)**0.5                                       ADS76530
      IF (JTER.GE.JTMAX) GO TO 100                                      ADS76540
C     EVALUATE GRADIENT AND CALL AGAIN.                                 ADS76550
      OBJ=OBJZ                                                          ADS76560
      INFO=2                                                            ADS76570
      ISFLAG=1                                                          ADS76580
      INFO2C=0                                                          ADS76590
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS76600
      IGOTO=2                                                           ADS76610
      RETURN                                                            ADS76620
90    CONTINUE                                                          ADS76630
C     NEW PENALIZED OBJECTIVE.                                          ADS76640
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS76650
     1IC)                                                               ADS76660
      OBJ=OBJZ+PENALT                                                   ADS76670
      GO TO 20                                                          ADS76680
100   CONTINUE                                                          ADS76690
C     ------------------------------------------------------------------ADS76700
C     OPTIMIZATION IS COMPLETE.                                         ADS76710
C     ------------------------------------------------------------------ADS76720
      INFO=0                                                            ADS76730
      IGOTO=0                                                           ADS76740
      OBJ=OBJZ                                                          ADS76750
      IF (IPRNT3.EQ.0) GO TO 140                                        ADS76760
      OBJW=OBJ/SCFO                                                     ADS76770
      WRITE (6,170) JTER,OBJW                                           ADS76780
C     X-VECTOR.                                                         ADS76790
      ICALC=1                                                           ADS76800
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS76810
C     G-VECTOR.                                                         ADS76820
      ICALC=4                                                           ADS76830
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS76840
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS76850
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS76860
C     SIDE CONSTRAINTS.                                                 ADS76870
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS76880
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS76890
      WRITE (6,180) CT,CTL                                              ADS76900
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS76910
      NAVC=NAC+NVC                                                      ADS76920
      WRITE (6,240) NAC,NVC                                             ADS76930
      IF (NAVC.GT.0) WRITE (6,250)                                      ADS76940
      IF (NAVC.GT.0) WRITE (6,260) (ITMP(I),I=1,NAVC)                   ADS76950
      WRITE (6,270) NACS                                                ADS76960
      IF (NACS.EQ.0) GO TO 130                                          ADS76970
      N=0                                                               ADS76980
      DO 120 I=1,NDV                                                    ADS76990
      IF (ISIDC(I).EQ.0) GO TO 120                                      ADS77000
      IF (ISIDC(I).EQ.1) GO TO 110                                      ADS77010
C     LOWER BOUND IS ACTIVE.                                            ADS77020
      N=N+1                                                             ADS77030
      ITMP(N)=-I                                                        ADS77040
110   IF (ISIDC(I).LT.0) GO TO 120                                      ADS77050
C     UPPER BOUND IS ACTIVE.                                            ADS77060
      N=N+1                                                             ADS77070
      ITMP(N)=I                                                         ADS77080
120   CONTINUE                                                          ADS77090
      WRITE (6,280)                                                     ADS77100
      WRITE (6,260) (ITMP(I),I=1,NACS)                                  ADS77110
130   WRITE (6,190)                                                     ADS77120
      IF (JTER.GE.JTMAX) WRITE (6,200)                                  ADS77130
      IF (IDEL.GE.ITRMST) WRITE (6,210) ITRMST                          ADS77140
      IF (IDAB.GE.ITRMST) WRITE (6,220) ITRMST                          ADS77150
      IF (RPPRIM.LE.RPPMIN) WRITE (6,230) RPPMIN                        ADS77160
140   CONTINUE                                                          ADS77170
      RETURN                                                            ADS77180
C     ------------------------------------------------------------------ADS77190
C     FORMATS.                                                          ADS77200
C     ------------------------------------------------------------------ADS77210
150   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS77220
160   FORMAT (/,5X,'UNCONSTRAINED MINIMIZATION IS COMPLETE')            ADS77230
170   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS77240
     1'NUMBER OF UNCONSTRAINED MINIMIZATIONS =',I5,//,5X,               ADS77250
     2'OBJECTIVE =',1PE12.5)                                            ADS77260
180   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,'  CTL =',      ADS77270
     11PE12.5)                                                          ADS77280
190   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS77290
200   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS77300
210   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS77310
     1' CONSECUTIVE ITERATIONS')                                        ADS77320
220   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS77330
     1' CONSECUTIVE ITERATIONS')                                        ADS77340
230   FORMAT (/,5X,'PENALTY PARAMETER IS LESS THAN',1PE12.5)            ADS77350
240   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS77360
     1' VIOLATED CONSTRAINTS')                                          ADS77370
250   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS77380
260   FORMAT (5X,10I5)                                                  ADS77390
270   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS77400
280   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS77410
290   FORMAT ('1',4X,75('-'),/,16X,                                     ADS77420
     1'ISTRAT = 2: LINEAR EXTENDED PENALTY FUNCTION METHOD',/,5X,       ADS77430
     275('-'))                                                          ADS77440
300   FORMAT (//,5X,'-- INITIAL DESIGN')                                ADS77450
310   FORMAT (/,5X,'SCALED OBJ =',1PE12.5,5X,'PENALT =',1PE12.5,/,5X,   ADS77460
     1'TRUE OBJ =',1PE12.5)                                             ADS77470
320   FORMAT (//,5X,'-- BEGIN SUMT ITERATION',I5,5X,'RPPRIM =',1PE12.5, ADS77480
     15X,'EPSPEN =',1PE12.5)                                            ADS77490
      END                                                               ADS77500
C --- SUBROUTINE ADS304 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,          ADS77510
C    1SVLB,SVUB,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX, ADS77520
C    2RP,RMULT,XTMP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,     ADS77530
C    3EPSPEN,RPPRIM,RPMULT,SCFX,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL)ADS77540
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),IDG(NCON),IC(NCOLA), ADS77550
C    1DF(NDV),A(NRA,NCOLA),WK(NRWK),IWK(NRIWK),XTMP(NDV),SCFX(NDV),     ADS77560
C    2ALAMDA(NCON),PPARAM(NCON)                                         ADS77570
C     ------------------------------------------------------------------ADS77580
C --- ISTRAT = 3.  PERFORM SEQUENTIAL UNCONSTRAINED MINIMIZATION USING  ADS77590
C     THE QUADRATIC EXTENDED INTERIOR PENALTY FUNCTION METHOD.          ADS77600
C     ------------------------------------------------------------------ADS77610
C     BY G. N. VANDERPLAATS                           AUG. 1982.        ADS77620
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS77630
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS77640
C --- INPUT.                                                            ADS77650
C      IGOTO  - CONTROL PARAMETER.                                      ADS77660
C      INFO   - CONTROL PARAMETER.                                      ADS77670
C      IOPT   - OPTIMIZATION IDENTIFIER.                                ADS77680
C      IONED  - ONE D. SEARCH IDENTIFIER.                               ADS77690
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS77700
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS77710
C      OBJ    - OBJECTIVE.                                              ADS77720
C      X      - CURRENT VALUE OF DESIGN VARIABLES.                      ADS77730
C      SVLB,SVUB - LOWER AND UPPER BOUND ON DESIGN VARIABLES.           ADS77740
C      G      - CONSTRAINT VALUES.                                      ADS77750
C      IDG    - CONSTRAINT TYPE IDENYIFIER.                             ADS77760
C      NGT   - NUMBER OF ACTIVE CONSTRAINTS.                            ADS77770
C      IC     - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                   ADS77780
C      DF     - GRADIENT OF OBKECTIVE.                                  ADS77790
C      A      - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS.       ADS77800
C      NRA    - NUMBER OF ROWS IN ARRAY A.                              ADS77810
C      NCOLA  - NUMBER OF COLUMNS IN ARRAY A.                           ADS77820
C      WK,IWK - REAL,INTEGER WORK ARRAY.                                ADS77830
C      IPRNT3 - PRINT CONTROL IDENTIFIER FOR STRATEGY LEVEL.            ADS77840
C      CT,CTL - PARAMETER TO IDENTIFY ACTIVE CONSTRAINTS.               ADS77850
C      ITRMST,CTMIN,CTLMIN,DELSTR,DABSTR - PARAMETERS FOR CONVERGENCE   ADS77860
C               CRITERIA.                                               ADS77870
C      JTMAX  - MAXIMUM NUMBER OF ITERATIONS.                           ADS77880
C      RP     - PENALTY PARAMETER FOR EXTERIOR PENALTY FUNCTION.        ADS77890
C      RMULT  - FACTOR FOR RP.                                          ADS77900
C      XTMP   - TEMPORARY ARRAY FOR X.                                  ADS77910
C      RPMAX  - MAXIMUM VALUE OF RP.                                    ADS77920
C      RPPMIN - MINIMUM VALUE OF RPPRIM.                                ADS77930
C      EPSPEN - TRANSITION POINT.                                       ADS77940
C      RPPRIM - PENALTY PARAMETER FOR INTERIOR PENALTY FUNCTION.        ADS77950
C      RPMULT - FACTOR FOR RPPRIM.                                      ADS77960
C      SCFX   - SCALE FACTOR FOR X.                                     ADS77970
C      ALAMDA - LAGRANGE MULTIPLIER.                                    ADS77980
C      PPARAM - PENALTY PARAMETER OF EACH CONSTRAINT.                   ADS77990
C      SCFO   - SCALE FACTOR OF OBJECTIVE.                              ADS78000
C      IFCALL - NUMBER OF FUNCTION EVALUATION.                          ADS78010
C      IGCALL - NUMBER OF GRADIENT EVALUATION.                          ADS78020
C --- OUTPUT.                                                           ADS78030
C      OBJ    - OPTIMUM OBJECTIVE.                                      ADS78040
C      X      - OPTIMUM VALUE OF DESIGN VARIABLES.                      ADS78050
C      G      - CONSTRAINT VALUE CORRESPOND TO OPTIMUM DESIGN VARIABLES.ADS78060
C      ISFLAG - SCALE CONTROL PARAMETER. =1; ADDITIONAL SCALE FACTOR IS ADS78070
C               CALCULATED. =0; ONLY GRADIENTS ARE REQUIRED.            ADS78080
C                                                                       ADS78090
      SUBROUTINE ADS304 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,ADS78100
     1G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,XTADS78110
     2MP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPADS78120
     3MULT,SCFX,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL,OBJZ,OBJZ1,PENALADS78130
     4T,IDAB,IDEL,INFO2C,JTER,ITMP,ISIDC,SCFG)                          ADS78140
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WKADS78150
     1(1),IWK(1),XTMP(1),SCFX(1),ALAMDA(1),PPARAM(1),ITMP(1),ISIDC(1),  ADS78160
     2SCFG(1)                                                           ADS78170
C                                                                       ADS78180
      ISTRAT=3                                                          ADS78190
      IGOTO1=IGOTO+1                                                    ADS78200
      GO TO (10,50,90), IGOTO1                                          ADS78210
10    CONTINUE                                                          ADS78220
      INFO=0                                                            ADS78230
      IF (IPRNT3.GT.0) WRITE (6,280)                                    ADS78240
C     ------------------------------------------------------------------ADS78250
C     FIRST CALL.  INITIALIZE PARAMETERS.                               ADS78260
C     ------------------------------------------------------------------ADS78270
      JTER=0                                                            ADS78280
      IDEL=0                                                            ADS78290
      IDAB=0                                                            ADS78300
      INFO2C=0                                                          ADS78310
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS78320
      IF (IPRNT3.GT.0) WRITE (6,290)                                    ADS78330
C     CALCULATE TRUE OBJECTIVE.                                         ADS78340
      CALL ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                           ADS78350
      OBJZ=OBJ                                                          ADS78360
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS78370
     1IC)                                                               ADS78380
      OBJ=OBJZ+PENALT                                                   ADS78390
      OBJW=OBJZ/SCFO                                                    ADS78400
      IF (IPRNT3.GT.0) WRITE (6,300) OBJ,PENALT,OBJW                    ADS78410
      ICALC=1                                                           ADS78420
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS78430
      ICALC=4                                                           ADS78440
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS78450
C     ------------------------------------------------------------------ADS78460
C     BEGIN NEW ITERATION.                                              ADS78470
C     ------------------------------------------------------------------ADS78480
20    JTER=JTER+1                                                       ADS78490
C     SAVE OBJZ.                                                        ADS78500
      OBJZ1=OBJZ                                                        ADS78510
      IF (IPRNT3.GT.0) WRITE (6,310) JTER,RPPRIM,EPSPEN                 ADS78520
      ICALC=3                                                           ADS78530
      IF (IPRNT3.LT.5) GO TO 40                                         ADS78540
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS78550
C     WRITE GRADIENT OF RETAINED CONSTRAINTS.                           ADS78560
      ISCAL=0                                                           ADS78570
      DO 30 J=1,NGT                                                     ADS78580
      WRITE (6,150) IC(J)                                               ADS78590
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS78600
30    CONTINUE                                                          ADS78610
40    CONTINUE                                                          ADS78620
C     ------------------------------------------------------------------ADS78630
C     OPTIMIZE UNCONSTRAINED PSEUDO-OBJECTIVE FUNCTION.                 ADS78640
C     ------------------------------------------------------------------ADS78650
50    CONTINUE                                                          ADS78660
      IF (INFO.NE.2) GO TO 60                                           ADS78670
      CALL ADS008 (NDV,ISTRAT,NCON,NGT,A,NRA,IC,G,IDG,DF,ALAMDA,PPARAM,RADS78680
     1P,EPSPEN)                                                         ADS78690
60    CONTINUE                                                          ADS78700
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDGADS78710
     1,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                                    ADS78720
C     FINISHED?                                                         ADS78730
      IF (INFO.EQ.2) INFO2C=INFO2C+1                                    ADS78740
      IF (INFO2C.EQ.1.AND.INFO.EQ.2) GO TO 50                           ADS78750
C     NO. EVALUATE FUNCTION OR GRADIENT AND CALL AGAIN.                 ADS78760
      IGOTO=1                                                           ADS78770
      IF (INFO.EQ.1) RETURN                                             ADS78780
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS78790
     1IC)                                                               ADS78800
C     TRUE OBJECTIVE.                                                   ADS78810
      OBJZ=OBJ-PENALT                                                   ADS78820
      IF (INFO.EQ.0) GO TO 70                                           ADS78830
      OBJ=OBJZ                                                          ADS78840
      ISFLAG=0                                                          ADS78850
      RETURN                                                            ADS78860
70    CONTINUE                                                          ADS78870
C     ------------------------------------------------------------------ADS78880
C     UNCONSTRAINED MINIMIZATION IS COMPLETE. UPDATE OBJ AND RP.        ADS78890
C     ------------------------------------------------------------------ADS78900
      OBJW=OBJZ/SCFO                                                    ADS78910
      IF (IPRNT3.GT.0) WRITE (6,300) OBJ,PENALT,OBJW                    ADS78920
      ICALC=1                                                           ADS78930
      IF (IPRNT3.GE.2) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS78940
      ICALC=4                                                           ADS78950
      IF (IPRNT3.GE.3) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS78960
C     ------------------------------------------------------------------ADS78970
C     CHECK CONVERGENCE.                                                ADS78980
C     ------------------------------------------------------------------ADS78990
      CALL ADS211 (OBJZ,OBJZ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)        ADS79000
C     MAXIMUM CONSTRAINT VIOLATION.                                     ADS79010
      GMCT=-1.                                                          ADS79020
      DO 80 I=1,NCON                                                    ADS79030
      CC=CTMIN                                                          ADS79040
      IF (IDG(I).GT.1) CC=CTLMIN                                        ADS79050
      GC=G(I)-CC                                                        ADS79060
      IF (GC.GT.GMCT) GMCT=GC                                           ADS79070
80    CONTINUE                                                          ADS79080
      IF (GMCT.GT.0.) IDEL=0                                            ADS79090
      IF (IDEL.GE.ITRMST.OR.IDAB.GE.ITRMST) GO TO 100                   ADS79100
      RPPRIM=RPMULT*RPPRIM                                              ADS79110
      IF (RPPRIM.LT.RPPMIN.AND.IDEL.GT.0) GO TO 100                     ADS79120
      IF (RPPRIM.LT.RPPMIN) RPPRIM=RPPMIN                               ADS79130
      RP=RP*RMULT                                                       ADS79140
      IF (RP.GT.RPMAX) RP=RPMAX                                         ADS79150
      CALL ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                           ADS79160
      EPSPEN=EPSPEN*(RPMULT)**0.45                                      ADS79170
      IF (JTER.GE.JTMAX) GO TO 100                                      ADS79180
C     EVALUATE GRADIENT AND CALL AGAIN.                                 ADS79190
      OBJ=OBJZ                                                          ADS79200
      INFO=2                                                            ADS79210
      ISFLAG=1                                                          ADS79220
      INFO2C=0                                                          ADS79230
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS79240
      IGOTO=2                                                           ADS79250
      RETURN                                                            ADS79260
90    CONTINUE                                                          ADS79270
C     NEW PENALIZED OBJECTIVE.                                          ADS79280
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS79290
     1IC)                                                               ADS79300
      OBJ=OBJZ+PENALT                                                   ADS79310
      GO TO 20                                                          ADS79320
100   CONTINUE                                                          ADS79330
C     ------------------------------------------------------------------ADS79340
C     OPTIMIZATION IS COMPLETE.                                         ADS79350
C     ------------------------------------------------------------------ADS79360
      INFO=0                                                            ADS79370
      IGOTO=0                                                           ADS79380
      OBJ=OBJZ                                                          ADS79390
      IF (IPRNT3.EQ.0) GO TO 140                                        ADS79400
      OBJW=OBJ/SCFO                                                     ADS79410
      WRITE (6,160) JTER,OBJW                                           ADS79420
C     X-VECTOR.                                                         ADS79430
      ICALC=1                                                           ADS79440
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS79450
C     G-VECTOR.                                                         ADS79460
      ICALC=4                                                           ADS79470
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS79480
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS79490
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS79500
C     SIDE CONSTRAINTS.                                                 ADS79510
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS79520
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS79530
      WRITE (6,170) CT,CTL                                              ADS79540
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS79550
      NAVC=NAC+NVC                                                      ADS79560
      WRITE (6,230) NAC,NVC                                             ADS79570
      IF (NAVC.GT.0) WRITE (6,240)                                      ADS79580
      IF (NAVC.GT.0) WRITE (6,250) (ITMP(I),I=1,NAVC)                   ADS79590
      WRITE (6,260) NACS                                                ADS79600
      IF (NACS.EQ.0) GO TO 130                                          ADS79610
      N=0                                                               ADS79620
      DO 120 I=1,NDV                                                    ADS79630
      IF (ISIDC(I).EQ.0) GO TO 120                                      ADS79640
      IF (ISIDC(I).EQ.1) GO TO 110                                      ADS79650
C     LOWER BOUND IS ACTIVE.                                            ADS79660
      N=N+1                                                             ADS79670
      ITMP(N)=-I                                                        ADS79680
110   IF (ISIDC(I).LT.0) GO TO 120                                      ADS79690
C     UPPER BOUND IS ACTIVE.                                            ADS79700
      N=N+1                                                             ADS79710
      ITMP(N)=I                                                         ADS79720
120   CONTINUE                                                          ADS79730
      WRITE (6,270)                                                     ADS79740
      WRITE (6,250) (ITMP(I),I=1,NACS)                                  ADS79750
130   WRITE (6,180)                                                     ADS79760
      IF (JTER.GE.JTMAX) WRITE (6,190)                                  ADS79770
      IF (IDEL.GE.ITRMST) WRITE (6,200) ITRMST                          ADS79780
      IF (IDAB.GE.ITRMST) WRITE (6,210) ITRMST                          ADS79790
      IF (RPPRIM.LE.RPPMIN) WRITE (6,220) RPPMIN                        ADS79800
140   CONTINUE                                                          ADS79810
      RETURN                                                            ADS79820
C     ------------------------------------------------------------------ADS79830
C     FORMATS.                                                          ADS79840
C     ------------------------------------------------------------------ADS79850
150   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS79860
160   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS79870
     1'NUMBER OF UNCONSTRAINED MINIMIZATIONS =',I5,//,5X,'OBJECTIVE =', ADS79880
     21PE12.5)                                                          ADS79890
170   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,2X,'CTL =',     ADS79900
     11PE12.5)                                                          ADS79910
180   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS79920
190   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS79930
200   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS79940
     1' CONSECUTIVE ITERATIONS')                                        ADS79950
210   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS79960
     1' CONSECUTIVE ITERATIONS')                                        ADS79970
220   FORMAT (/,5X,'PENALTY PARAMETER IS LESS THAN',1PE12.5)            ADS79980
230   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS79990
     1' VIOLATED CONSTRAINTS')                                          ADS80000
240   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS80010
250   FORMAT (5X,10I5)                                                  ADS80020
260   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS80030
270   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS80040
280   FORMAT ('1',4X,75('-'),/,15X,                                     ADS80050
     1'ISTRAT = 3: QUADRATIC EXTENDED PENALTY FUNCTION METHOD',/,       ADS80060
     15X,75('-'))                                                       ADS80070
290   FORMAT (/,5X,'-- INITIAL DESIGN')                                 ADS80080
300   FORMAT (/,5X,'SCALED OBJ =',1PE12.5,5X,'PENALT =',1PE12.5,/,5X,   ADS80090
     1'TRUE OBJ =',1PE12.5)                                             ADS80100
310   FORMAT (/,5X,'-- BEGIN SUMT ITERATION',I5,5X,'RPPRIM =',1PE12.5,  ADS80110
     15X,'EPSPEN =',1PE12.5)                                            ADS80120
      END                                                               ADS80130
C --- SUBROUTINE ADS305 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,          ADS80140
C    1SVLB,SVUB,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX, ADS80150
C    2RP,RMULT,XTMP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,     ADS80160
C    3EPSPEN,RPPRIM,RPMULT,SCFX,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL)ADS80170
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),IDG(NCON),IC(NCOLA), ADS80180
C    1DF(NDV),A(NRA,NCOLA),WK(NRWK),IWK(NRIWK),XTMP(NDV),SCFX(NDV),     ADS80190
C    2ALAMDA(NCON),PPARAM(NCON)                                         ADS80200
C     ------------------------------------------------------------------ADS80210
C --- ISTRAT = 4.  PERFORM SEQUENTIAL UNCONSTRAINED MINIMIZATION USING  ADS80220
C     THE CUBIC EXTENDED INTERIOR PENALTY FUNCTION METHOD.              ADS80230
C     ------------------------------------------------------------------ADS80240
C     BY G. N. VANDERPLAATS                           AUG. 1982.        ADS80250
C     MODIFIED BY H.SUGIMOTO                          NOV. 1982.        ADS80260
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS80270
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS80280
C --- INPUT.                                                            ADS80290
C      IGOTO  - CONTROL PARAMETER.                                      ADS80300
C      INFO   - CONTROL PARAMETER.                                      ADS80310
C      IOPT   - OPTIMIZATION IDENTIFIER.                                ADS80320
C      IONED  - ONE D. SEARCH IDENTIFIER.                               ADS80330
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS80340
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS80350
C      OBJ    - OBJECTIVE.                                              ADS80360
C      X      - CURRENT VALUE OF DESIGN VARIABLES.                      ADS80370
C      SVLB,SVUB - LOWER AND UPPER BOUND ON DESIGN VARIABLES.           ADS80380
C      G      - CONSTRAINT VALUES.                                      ADS80390
C      IDG    - CONSTRAINT TYPE IDENYIFIER.                             ADS80400
C      NGT   - NUMBER OF ACTIVE CONSTRAINTS.                            ADS80410
C      IC     - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                   ADS80420
C      DF     - GRADIENT OF OBKECTIVE.                                  ADS80430
C      A      - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS.       ADS80440
C      NRA    - NUMBER OF ROWS IN ARRAY A.                              ADS80450
C      NCOLA  - NUMBER OF COLUMNS IN ARRAY A.                           ADS80460
C      WK,IWK - REAL,INTEGER WORK ARRAY.                                ADS80470
C      IPRNT3 - PRINT CONTROL IDENTIFIER FOR STRATEGY LEVEL.            ADS80480
C      CT,CTL - PARAMETERS TO IDENTIFY ACTIVE CONSTRAINTS.              ADS80490
C      ITRMST,CTMIN,CTLMIN,DELSTR,DABSTR - PARAMETERS FOR CONVERGENCE   ADS80500
C               CRITERIA.                                               ADS80510
C      JTMAX  - MAXIMUM NUMBER OF ITERATIONS.                           ADS80520
C      RP     - PENALTY PARAMETER FOR EXTERIOR PENALTY FUNCTION.        ADS80530
C      RMULT  - FACTOR FOR RP.                                          ADS80540
C      XTMP   - TEMPORARY ARRAY FOR X.                                  ADS80550
C      RPMAX  - MAXIMUM VALUE OF RP.                                    ADS80560
C      RPPMIN - MINIMUM VALUE OF RPPRIM.                                ADS80570
C      EPSPEN - TRANSITION POINT.                                       ADS80580
C      RPPRIM - PENALTY PARAMETER FOR INTERIOR PENALTY FUNCTION.        ADS80590
C      RPMULT - FACTOR FOR RPPRIM.                                      ADS80600
C      SCFX   - SCALE FACTOR FOR X.                                     ADS80610
C      ALAMDA - LAGRANGE MULTIPLIER.                                    ADS80620
C      PPARAM - PENALTY PARAMETER OF EACH CONSTRAINT.                   ADS80630
C      SCFO   - SCALE FACTOR OF OBJECTIVE.                              ADS80640
C      IFCALL - NUMBER OF FUNCTION EVALUATION.                          ADS80650
C      IGCALL - NUMBER OF GRADIENT EVALUATION.                          ADS80660
C --- OUTPUT.                                                           ADS80670
C      OBJ    - OPTIMUM OBJECTIVE.                                      ADS80680
C      X      - OPTIMUM VALUE OF DESIGN VARIABLES.                      ADS80690
C      G      - CONSTRAINT VALUE CORRESPOND TO OPTIMUM DESIGN VARIABLES.ADS80700
C      ISFLAG - SCALE CONTROL PARAMETER. =1; ADDITIONAL SCALE FACTOR IS ADS80710
C               CALCULATED. =0; ONLY GRADIENTS ARE REQUIRED.            ADS80720
C                                                                       ADS80730
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS80740
      SUBROUTINE ADS305 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,ADS80750
     1G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,XTADS80760
     2MP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPADS80770
     3MULT,SCFX,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL,OBJZ,OBJZ1,PENALADS80780
     4T,IDAB,IDEL,INFO2C,JTER,ITMP,ISIDC,SCFG)                          ADS80790
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WKADS80800
     1(1),IWK(1),XTMP(1),SCFX(1),ALAMDA(1),PPARAM(1),ITMP(1),ISIDC(1),  ADS80810
     2SCFG(1)                                                           ADS80820
C                                                                       ADS80830
      ISTRAT=4                                                          ADS80840
      IGOTO1=IGOTO+1                                                    ADS80850
      GO TO (10,50,90), IGOTO1                                          ADS80860
10    CONTINUE                                                          ADS80870
      INFO=0                                                            ADS80880
      IF (IPRNT3.GT.0) WRITE (6,290)                                    ADS80890
C     ------------------------------------------------------------------ADS80900
C     FIRST CALL.  INITIALIZE PARAMETERS.                               ADS80910
C     ------------------------------------------------------------------ADS80920
      JTER=0                                                            ADS80930
      IDEL=0                                                            ADS80940
      IDAB=0                                                            ADS80950
      INFO2C=0                                                          ADS80960
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS80970
      IF (IPRNT3.GT.0) WRITE (6,300)                                    ADS80980
C     CALCULATE TRUE OBJECTIVE.                                         ADS80990
      CALL ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                           ADS81000
      OBJZ=OBJ                                                          ADS81010
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS81020
     1IC)                                                               ADS81030
      OBJ=OBJZ+PENALT                                                   ADS81040
      OBJW=OBJZ/SCFO                                                    ADS81050
      IF (IPRNT3.GT.0) WRITE (6,310) OBJ,PENALT,OBJW                    ADS81060
      ICALC=1                                                           ADS81070
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS81080
      ICALC=4                                                           ADS81090
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS81100
C     ------------------------------------------------------------------ADS81110
C     BEGIN NEW ITERATION.                                              ADS81120
C     ------------------------------------------------------------------ADS81130
20    JTER=JTER+1                                                       ADS81140
C     SAVE OBJZ.                                                        ADS81150
      OBJZ1=OBJZ                                                        ADS81160
      IF (IPRNT3.GT.0) WRITE (6,320) JTER,RPPRIM,EPSPEN                 ADS81170
      ICALC=3                                                           ADS81180
      IF (IPRNT3.LT.5) GO TO 40                                         ADS81190
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS81200
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS81210
      ISCAL=0                                                           ADS81220
      DO 30 J=1,NGT                                                     ADS81230
      WRITE (6,150) IC(J)                                               ADS81240
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS81250
30    CONTINUE                                                          ADS81260
40    CONTINUE                                                          ADS81270
C     ------------------------------------------------------------------ADS81280
C     OPTIMIZE UNCONSTRAINED PSEUDO-OBJECTIVE FUNCTION.                 ADS81290
C     ------------------------------------------------------------------ADS81300
50    CONTINUE                                                          ADS81310
      IF (INFO.NE.2) GO TO 60                                           ADS81320
      CALL ADS008 (NDV,ISTRAT,NCON,NGT,A,NRA,IC,G,IDG,DF,ALAMDA,PPARAM,RADS81330
     1P,EPSPEN)                                                         ADS81340
60    CONTINUE                                                          ADS81350
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDGADS81360
     1,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                                    ADS81370
C     FINISHED?                                                         ADS81380
      IF (INFO.EQ.2) INFO2C=INFO2C+1                                    ADS81390
      IF (INFO2C.EQ.1.AND.INFO.EQ.2) GO TO 50                           ADS81400
C     NO. EVALUATE FUNCTION OR GRADIENT AND CALL AGAIN.                 ADS81410
      IGOTO=1                                                           ADS81420
      IF (INFO.EQ.1) RETURN                                             ADS81430
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS81440
     1IC)                                                               ADS81450
C     TRUE OBJECTIVE.                                                   ADS81460
      OBJZ=OBJ-PENALT                                                   ADS81470
      IF (INFO.EQ.0) GO TO 70                                           ADS81480
      OBJ=OBJZ                                                          ADS81490
      ISFLAG=0                                                          ADS81500
      RETURN                                                            ADS81510
70    CONTINUE                                                          ADS81520
C     ------------------------------------------------------------------ADS81530
C     UNCONSTRAINED MINIMIZATION IS COMPLETE. UPDATE OBJ AND RP.        ADS81540
C     ------------------------------------------------------------------ADS81550
      IF (IPRNT3.GT.3) WRITE (6,160)                                    ADS81560
      OBJW=OBJZ/SCFO                                                    ADS81570
      IF (IPRNT3.GT.0) WRITE (6,310) OBJ,PENALT,OBJW                    ADS81580
      ICALC=1                                                           ADS81590
      IF (IPRNT3.GE.2) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS81600
      ICALC=4                                                           ADS81610
      IF (IPRNT3.GE.3) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS81620
C     ------------------------------------------------------------------ADS81630
C     CHECK CONVERGENCE.                                                ADS81640
C     ------------------------------------------------------------------ADS81650
      CALL ADS211 (OBJZ,OBJZ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)        ADS81660
C     MAXIMUM CONSTRAINT VIOLATION.                                     ADS81670
      GMCT=-1.                                                          ADS81680
      DO 80 I=1,NCON                                                    ADS81690
      CC=CTMIN                                                          ADS81700
      IF (IDG(I).GT.1) CC=CTLMIN                                        ADS81710
      GC=G(I)-CC                                                        ADS81720
      IF (GC.GT.GMCT) GMCT=GC                                           ADS81730
80    CONTINUE                                                          ADS81740
      IF (GMCT.GT.0.) IDEL=0                                            ADS81750
      IF (IDEL.GE.ITRMST.OR.IDAB.GE.ITRMST) GO TO 100                   ADS81760
      RPPRIM=RPMULT*RPPRIM                                              ADS81770
      IF (RPPRIM.LT.RPPMIN.AND.IDEL.GT.0) GO TO 100                     ADS81780
      IF (RPPRIM.LT.RPPMIN) RPPRIM=RPPMIN                               ADS81790
      RP=RP*RMULT                                                       ADS81800
      IF (RP.GT.RPMAX) RP=RPMAX                                         ADS81810
      CALL ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                           ADS81820
      EPSPEN=EPSPEN*(RPMULT)**0.40                                      ADS81830
      IF (JTER.GE.JTMAX) GO TO 100                                      ADS81840
C     EVALUATE GRADIENT AND CALL AGAIN.                                 ADS81850
      OBJ=OBJZ                                                          ADS81860
      INFO=2                                                            ADS81870
      ISFLAG=1                                                          ADS81880
      INFO2C=0                                                          ADS81890
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS81900
      IGOTO=2                                                           ADS81910
      RETURN                                                            ADS81920
90    CONTINUE                                                          ADS81930
C     NEW PENALIZED OBJECTIVE.                                          ADS81940
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS81950
     1IC)                                                               ADS81960
      OBJ=OBJZ+PENALT                                                   ADS81970
      GO TO 20                                                          ADS81980
100   CONTINUE                                                          ADS81990
C     ------------------------------------------------------------------ADS82000
C     OPTIMIZATION IS COMPLETE.                                         ADS82010
C     ------------------------------------------------------------------ADS82020
      INFO=0                                                            ADS82030
      IGOTO=0                                                           ADS82040
      OBJ=OBJZ                                                          ADS82050
      IF (IPRNT3.EQ.0) GO TO 140                                        ADS82060
      OBJW=OBJ/SCFO                                                     ADS82070
      WRITE (6,170) JTER,OBJW                                           ADS82080
C     X-VECTOR.                                                         ADS82090
      ICALC=1                                                           ADS82100
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS82110
C     G-VECTOR.                                                         ADS82120
      ICALC=4                                                           ADS82130
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS82140
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS82150
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS82160
C     SIDE CONSTRAINTS.                                                 ADS82170
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS82180
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS82190
      WRITE (6,180) CT,CTL                                              ADS82200
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS82210
      NAVC=NAC+NVC                                                      ADS82220
      WRITE (6,240) NAC,NVC                                             ADS82230
      IF (NAVC.GT.0) WRITE (6,250)                                      ADS82240
      IF (NAVC.GT.0) WRITE (6,260) (ITMP(I),I=1,NAVC)                   ADS82250
      WRITE (6,270) NACS                                                ADS82260
      IF (NACS.EQ.0) GO TO 130                                          ADS82270
      N=0                                                               ADS82280
      DO 120 I=1,NDV                                                    ADS82290
      IF (ISIDC(I).EQ.0) GO TO 120                                      ADS82300
      IF (ISIDC(I).EQ.1) GO TO 110                                      ADS82310
C     LOWER BOUND IS ACTIVE.                                            ADS82320
      N=N+1                                                             ADS82330
      ITMP(N)=-I                                                        ADS82340
110   IF (ISIDC(I).LT.0) GO TO 120                                      ADS82350
C     UPPER BOUND IS ACTIVE.                                            ADS82360
      N=N+1                                                             ADS82370
      ITMP(N)=I                                                         ADS82380
120   CONTINUE                                                          ADS82390
      WRITE (6,280)                                                     ADS82400
      WRITE (6,260) (ITMP(I),I=1,NACS)                                  ADS82410
130   WRITE (6,190)                                                     ADS82420
      IF (JTER.GE.JTMAX) WRITE (6,200)                                  ADS82430
      IF (IDEL.GE.ITRMST) WRITE (6,210) ITRMST                          ADS82440
      IF (IDAB.GE.ITRMST) WRITE (6,220) ITRMST                          ADS82450
      IF (RPPRIM.LE.RPPMIN) WRITE (6,230) RPPMIN                        ADS82460
140   CONTINUE                                                          ADS82470
      RETURN                                                            ADS82480
C     ------------------------------------------------------------------ADS82490
C     FORMATS.                                                          ADS82500
C     ------------------------------------------------------------------ADS82510
150   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS82520
160   FORMAT (/,5X,'UNCONSTRAINED MINIMIZATION IS COMPLETE')            ADS82530
170   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS82540
     1'NUMBER OF UNCONSTRAINED MINIMIZATIONS =',I5,//,5X,'OBJECTIVE =', ADS82550
     21PE12.5)                                                          ADS82560
180   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,'  CTL =',      ADS82570
     11PE12.5)                                                          ADS82580
190   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS82590
200   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS82600
210   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS82610
     1' CONSECUTIVE ITERATIONS')                                        ADS82620
220   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS82630
     1' CONSECUTIVE ITERATIONS')                                        ADS82640
230   FORMAT (/,5X,'PENALTY PARAMETER IS LESS THAN',1PE12.5)            ADS82650
240   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS82660
     1' VIOLATED CONSTRAINTS')                                          ADS82670
250   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS82680
260   FORMAT (5X,10I5)                                                  ADS82690
270   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS82700
280   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS82710
290   FORMAT ('1',4X,75('-'),/,16X,                                     ADS82720
     1'ISTRAT = 4: CUBIC EXTENDED PENALTY FUNCTION METHOD',/,5X,75('-'))ADS82730
300   FORMAT (/,5X,'-- INITIAL DESIGN')                                 ADS82740
310   FORMAT (/,5X,'SCALED OBJ =',1PE12.5,5X,'PENALT =',1PE12.5,/,5X,   ADS82750
     1'TRUE OBJ =',1PE12.5)                                             ADS82760
320   FORMAT (/,5X,'-- BEGIN SUMT ITERATION',I5,5X,'RPPRIM =',1PE12.5,  ADS82770
     15X,'EPSPEN =',1PE12.5)                                            ADS82780
      END                                                               ADS82790
C --- SUBROUTINE ADS306 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,          ADS82800
C    1SVLB,SVUB,G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX, ADS82810
C    2RP,RMULT,XTMP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,     ADS82820
C    3EPSPEN,RPPRIM,RPMULT,SCFX,ALAMDZ,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,ADS82830
C    4IGCALL)                                                           ADS82840
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),G(NCON),IDG(NCON),IC(NCOLA), ADS82850
C    1DF(NDV),A(NRA,NCOLA),WK(NRWK),IWK(NRIWK),XTMP(NDV),SCFX(NDV),     ADS82860
C    2ALAMDA(NCON),PPARAM(NCON)                                         ADS82870
C     ------------------------------------------------------------------ADS82880
C --- ISTRAT = 5.  PERFORM SEQUENTIAL UNCONSTRAINED MINIMIZATION USING  ADS82890
C     THE AUGMENTED LAGRANGEAN MULTIPLIER METHOD.                       ADS82900
C     ------------------------------------------------------------------ADS82910
C     BY G. N. VANDERPLAATS                           AUG. 1982.        ADS82920
C     MODIFIED BY H.SUGIMOTO                     NOV. DEC. 1982.        ADS82930
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS82940
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS82950
C --- INPUT.                                                            ADS82960
C      IGOTO  - CONTROL PARAMETER.                                      ADS82970
C      INFO   - CONTROL PARAMETER.                                      ADS82980
C      IOPT   - OPTIMIZATION IDENTIFIER.                                ADS82990
C      IONED  - ONE D. SEARCH IDENTIFIER.                               ADS83000
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS83010
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS83020
C      OBJ    - OBJECTIVE.                                              ADS83030
C      X      - CURRENT VALUE OF DESIGN VARIABLES.                      ADS83040
C      SVLB,SVUB - LOWER AND UPPER BOUND ON DESIGN VARIABLES.           ADS83050
C      G      - CONSTRAINT VALUES.                                      ADS83060
C      IDG    - CONSTRAINT TYPE IDENYIFIER.                             ADS83070
C      NGT   - NUMBER OF ACTIVE CONSTRAINTS.                            ADS83080
C      IC     - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                   ADS83090
C      DF     - GRADIENT OF OBKECTIVE.                                  ADS83100
C      A      - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS.       ADS83110
C      NRA    - NUMBER OF ROWS IN ARRAY A.                              ADS83120
C      NCOLA  - NUMBER OF COLUMNS IN ARRAY A.                           ADS83130
C      WK,IWK - REAL,INTEGER WORK ARRAY.                                ADS83140
C      IPRNT3 - PRINT CONTROL IDENTIFIER FOR STRATEGY LEVEL.            ADS83150
C      CT,CTL - PARAMETERS TO IDENTIFY ACTIVE CONSTRAINTS.              ADS83160
C      ITRMST,CTMIN,CTLMIN,DELSTR,DABSTR - PARAMETERS FOR CONVERGENCE   ADS83170
C               CRITERIA.                                               ADS83180
C      JTMAX  - MAXIMUM NUMBER OF ITERATIONS.                           ADS83190
C      RP     - PENALTY PARAMETER FOR EXTERIOR PENALTY FUNCTION.        ADS83200
C      RMULT  - FACTOR FOR RP.                                          ADS83210
C      XTMP   - TEMPORARY ARRAY FOR X.                                  ADS83220
C      RPMAX  - MAXIMUM VALUE OF RP.                                    ADS83230
C      RPPMIN - MINIMUM VALUE OF RPPRIM.                                ADS83240
C      EPSPEN - TRANSITION POINT.                                       ADS83250
C      RPPRIM - PENALTY PARAMETER FOR INTERIOR PENALTY FUNCTION.        ADS83260
C      RPMULT - FACTOR FOR RPPRIM.                                      ADS83270
C      SCFX   - SCALE FACTOR FOR X.                                     ADS83280
C      ALAMDZ - INITIAL VALUE OF LGRANGE MULTIPLIER.                    ADS83290
C      ALAMDA - LAGRANGE MULTIPLIER.                                    ADS83300
C      PPARAM - PENALTY PARAMETER OF EACH CONSTRAINT.                   ADS83310
C      SCFO   - SCALE FACTOR OF OBJECTIVE.                              ADS83320
C      IFCALL - NUMBER OF FUNCTION EVALUATION.                          ADS83330
C      IGCALL - NUMBER OF GRADIENT EVALUATION.                          ADS83340
C --- OUTPUT.                                                           ADS83350
C      OBJ    - OPTIMUM OBJECTIVE.                                      ADS83360
C      X      - OPTIMUM VALUE OF DESIGN VARIABLES.                      ADS83370
C      G      - CONSTRAINT VALUE CORRESPOND TO OPTIMUM DESIGN VARIABLES.ADS83380
C      ISFLAG - SCALE CONTROL PARAMETER. =1; ADDITIONAL SCALE FACTOR IS ADS83390
C               CALCULATED. =0; ONLY GRADIENTS ARE REQUIRED.            ADS83400
C                                                                       ADS83410
      SUBROUTINE ADS306 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,ADS83420
     1G,IDG,NGT,IC,DF,A,NRA,NCOLA,WK,IWK,IPRNT3,ITRMST,JTMAX,RP,RMULT,XTADS83430
     2MP,CT,CTL,CTMIN,CTLMIN,DELSTR,DABSTR,RPMAX,RPPMIN,EPSPEN,RPPRIM,RPADS83440
     3MULT,SCFX,ALAMDZ,ALAMDA,PPARAM,SCFO,ISFLAG,IFCALL,IGCALL,CK,OBJZ,OADS83450
     4BJZ1,PENALT,IDAB,IDEL,INFO2C,JTER,ITMP,ISIDC,SCFG)                ADS83460
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),G(1),IDG(1),IC(1),DF(1),A(NRA,1),WKADS83470
     1(1),IWK(1),XTMP(1),SCFX(1),ALAMDA(1),PPARAM(1),ITMP(1),ISIDC(1),  ADS83480
     2SCFG(1)                                                           ADS83490
C                                                                       ADS83500
      ISTRAT=5                                                          ADS83510
      BETA=.1                                                           ADS83520
      IGOTO1=IGOTO+1                                                    ADS83530
      GO TO (10,60,100), IGOTO1                                         ADS83540
10    CONTINUE                                                          ADS83550
      INFO=0                                                            ADS83560
      IF (IPRNT3.GT.0) WRITE (6,300)                                    ADS83570
C     ------------------------------------------------------------------ADS83580
C     FIRST CALL.  INITIALIZE PARAMETERS.                               ADS83590
C     ------------------------------------------------------------------ADS83600
      JTER=0                                                            ADS83610
      IDEL=0                                                            ADS83620
      IDAB=0                                                            ADS83630
      INFO2C=0                                                          ADS83640
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS83650
C     DEFINE THE INITIAL VALUE OF LAGRANGEAN MULTIPLIER.                ADS83660
      CK=1.0E10                                                         ADS83670
      DO 20 I=1,NCON                                                    ADS83680
      PPARAM(I)=RP                                                      ADS83690
      IF (G(I).GE.0.1) PPARAM(I)=RMULT*PPARAM(I)                        ADS83700
      IF (G(I).GT.0.0) PPARAM(I)=RMULT*PPARAM(I)                        ADS83710
      ALAMDA(I)=ALAMDZ                                                  ADS83720
20    CONTINUE                                                          ADS83730
      IF (IPRNT3.GT.0) WRITE (6,310)                                    ADS83740
C     CALCULATE TRUE OBJECTIVE.                                         ADS83750
      OBJZ=OBJ                                                          ADS83760
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS83770
     1IC)                                                               ADS83780
      OBJ=OBJZ+PENALT                                                   ADS83790
      OBJW=OBJZ/SCFO                                                    ADS83800
      IF (IPRNT3.GT.0) WRITE (6,320) OBJ,PENALT,OBJW                    ADS83810
      ICALC=1                                                           ADS83820
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS83830
      ICALC=4                                                           ADS83840
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS83850
C     ------------------------------------------------------------------ADS83860
C     BEGIN NEW ITERATION.                                              ADS83870
C     ------------------------------------------------------------------ADS83880
30    JTER=JTER+1                                                       ADS83890
C     SAVE OBJZ.                                                        ADS83900
      OBJZ1=OBJZ                                                        ADS83910
      IF (IPRNT3.GT.0) WRITE (6,330) JTER,RP                            ADS83920
      ICALC=3                                                           ADS83930
      IF (IPRNT3.LT.5) GO TO 50                                         ADS83940
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS83950
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS83960
      ISCAL=0                                                           ADS83970
      DO 40 J=1,NGT                                                     ADS83980
      WRITE (6,160) IC(J)                                               ADS83990
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS84000
40    CONTINUE                                                          ADS84010
50    CONTINUE                                                          ADS84020
C     ------------------------------------------------------------------ADS84030
C     OPTIMIZE UNCONSTRAINED PSEUDO-OBJECTIVE FUNCTION.                 ADS84040
C     ------------------------------------------------------------------ADS84050
      INFO=0                                                            ADS84060
      ICALC=10                                                          ADS84070
      IF (IPRNT3.GE.4) CALL ADS402 (ICALC,ALAMDA,NCON,SCFX)             ADS84080
      ICALC=11                                                          ADS84090
      IF (IPRNT3.GE.4) CALL ADS402 (ICALC,PPARAM,NCON,SCFX)             ADS84100
60    CONTINUE                                                          ADS84110
      IF (INFO.NE.2) GO TO 70                                           ADS84120
      CALL ADS008 (NDV,ISTRAT,NCON,NGT,A,NRA,IC,G,IDG,DF,ALAMDA,PPARAM,RADS84130
     1P,EPSPEN)                                                         ADS84140
70    CONTINUE                                                          ADS84150
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NCON,OBJ,X,SVLB,SVUB,G,IDGADS84160
     1,NGT,IC,DF,A,NRA,NCOLA,WK,IWK)                                    ADS84170
C     FINISHED?                                                         ADS84180
      IF (INFO.EQ.2) INFO2C=INFO2C+1                                    ADS84190
      IF (INFO2C.EQ.1.AND.INFO.EQ.2) GO TO 60                           ADS84200
C     NO. EVALUATE FUNCTION OR GRADIENT AND CALL AGAIN.                 ADS84210
      IGOTO=1                                                           ADS84220
      IF (INFO.EQ.1) RETURN                                             ADS84230
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS84240
     1IC)                                                               ADS84250
C     TRUE OBJECTIVE.                                                   ADS84260
      OBJZ=OBJ-PENALT                                                   ADS84270
      IF (INFO.EQ.0) GO TO 80                                           ADS84280
      OBJ=OBJZ                                                          ADS84290
      ISFLAG=0                                                          ADS84300
      RETURN                                                            ADS84310
80    CONTINUE                                                          ADS84320
C     ------------------------------------------------------------------ADS84330
C     UNCONSTRAINED MINIMIZATION IS COMPLETE. UPDATE OBJ AND RP.        ADS84340
C     ------------------------------------------------------------------ADS84350
      IF (IPRNT3.GT.3) WRITE (6,170)                                    ADS84360
      OBJW=OBJZ/SCFO                                                    ADS84370
      IF (IPRNT3.GT.0) WRITE (6,320) OBJ,PENALT,OBJW                    ADS84380
      ICALC=1                                                           ADS84390
      IF (IPRNT3.GE.2) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS84400
      ICALC=4                                                           ADS84410
      IF (IPRNT3.GE.3) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS84420
C     ------------------------------------------------------------------ADS84430
C     CHECK CONVERGENCE.                                                ADS84440
C     ------------------------------------------------------------------ADS84450
      CALL ADS211 (OBJZ,OBJZ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)        ADS84460
C     MAXIMUM CONSTRAINT VIOLATION.                                     ADS84470
      GMCT=-1.                                                          ADS84480
      DO 90 I=1,NCON                                                    ADS84490
      CC=CTMIN                                                          ADS84500
      IF (IDG(I).GT.1) CC=CTLMIN                                        ADS84510
      GC=G(I)-CC                                                        ADS84520
      IF (GC.GT.GMCT) GMCT=GC                                           ADS84530
90    CONTINUE                                                          ADS84540
      IF (GMCT.GT.0.) IDEL=0                                            ADS84550
      IF (IDEL.GE.ITRMST.OR.IDAB.GE.ITRMST) GO TO 110                   ADS84560
C     UPDATE THE LAGRANGEAN MULTIPLIER.                                 ADS84570
      CALL ADS307 (NCON,BETA,CK,G,PPARAM,ALAMDA,IDG,RMULT,RPMAX)        ADS84580
      IF (JTER.GT.3.AND.CK.LT.0.004) GO TO 110                          ADS84590
      IF (JTER.GE.JTMAX) GO TO 110                                      ADS84600
C     EVALUATE GRADIENT AND CALL AGAIN.                                 ADS84610
      OBJ=OBJZ                                                          ADS84620
      INFO=2                                                            ADS84630
      ISFLAG=1                                                          ADS84640
      INFO2C=0                                                          ADS84650
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS84660
      IGOTO=2                                                           ADS84670
      RETURN                                                            ADS84680
100   CONTINUE                                                          ADS84690
C     NEW PENALIZED OBJECTIVE.                                          ADS84700
      CALL ADS302 (ISTRAT,NCON,G,IDG,EPSPEN,RP,PENALT,ALAMDA,PPARAM,NGT,ADS84710
     1IC)                                                               ADS84720
      OBJ=OBJZ+PENALT                                                   ADS84730
      GO TO 30                                                          ADS84740
110   CONTINUE                                                          ADS84750
C     ------------------------------------------------------------------ADS84760
C     OPTIMIZATION IS COMPLETE.                                         ADS84770
C     ------------------------------------------------------------------ADS84780
      INFO=0                                                            ADS84790
      IGOTO=0                                                           ADS84800
      OBJ=OBJZ                                                          ADS84810
      IF (IPRNT3.EQ.0) GO TO 150                                        ADS84820
      OBJW=OBJ/SCFO                                                     ADS84830
      WRITE (6,180) JTER,OBJW                                           ADS84840
C     X-VECTOR.                                                         ADS84850
      ICALC=1                                                           ADS84860
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS84870
C     G-VECTOR.                                                         ADS84880
      ICALC=4                                                           ADS84890
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS84900
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS84910
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS84920
C     SIDE CONSTRAINTS.                                                 ADS84930
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS84940
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS84950
      WRITE (6,200) CT,CTL                                              ADS84960
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS84970
      NAVC=NAC+NVC                                                      ADS84980
      WRITE (6,250) NAC,NVC                                             ADS84990
      IF (NAVC.GT.0) WRITE (6,260)                                      ADS85000
      IF (NAVC.GT.0) WRITE (6,270) (ITMP(I),I=1,NAVC)                   ADS85010
      WRITE (6,280) NACS                                                ADS85020
      IF (NACS.EQ.0) GO TO 140                                          ADS85030
      N=0                                                               ADS85040
      DO 130 I=1,NDV                                                    ADS85050
      IF (ISIDC(I).EQ.0) GO TO 130                                      ADS85060
      IF (ISIDC(I).EQ.1) GO TO 120                                      ADS85070
C     LOWER BOUND IS ACTIVE.                                            ADS85080
      N=N+1                                                             ADS85090
      ITMP(N)=-I                                                        ADS85100
120   IF (ISIDC(I).LT.0) GO TO 130                                      ADS85110
C     UPPER BOUND IS ACTIVE.                                            ADS85120
      N=N+1                                                             ADS85130
      ITMP(N)=I                                                         ADS85140
130   CONTINUE                                                          ADS85150
      WRITE (6,290)                                                     ADS85160
      WRITE (6,270) (ITMP(I),I=1,NACS)                                  ADS85170
140   WRITE (6,210)                                                     ADS85180
      IF (JTER.GE.JTMAX) WRITE (6,220)                                  ADS85190
      IF (JTER.GT.3.AND.CK.LT.0.004) WRITE (6,190)                      ADS85200
      IF (IDEL.GE.ITRMST) WRITE (6,230) ITRMST                          ADS85210
      IF (IDAB.GE.ITRMST) WRITE (6,240) ITRMST                          ADS85220
150   CONTINUE                                                          ADS85230
      RETURN                                                            ADS85240
C     ------------------------------------------------------------------ADS85250
C                      FORMAT                                           ADS85260
C     ------------------------------------------------------------------ADS85270
160   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS85280
170   FORMAT (/,5X,'UNCONSTRAINED MINIMIZATION IS COMPLETE')            ADS85290
180   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS85300
     1'NUMBER OF UNCONSTRAINED MINIMIZATIONS =',I5,//,5X,'OBJECTIVE =', ADS85310
     21PE12.5)                                                          ADS85320
190   FORMAT (/,5X,'KUHN-TUCKER PARAMETER, CK, IS LESS THAN 0.004')     ADS85330
200   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,'  CTL =',      ADS85340
     11PE12.5)                                                          ADS85350
210   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS85360
220   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS85370
230   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS85380
     1' CONSECUTIVE ITERATIONS')                                        ADS85390
240   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS85400
     1' CONSECUTIVE ITERATIONS')                                        ADS85410
250   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS85420
     1' VIOLATED CONSTRAINTS')                                          ADS85430
260   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS85440
270   FORMAT (5X,10I5)                                                  ADS85450
280   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS85460
290   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS85470
300   FORMAT ('1',4X,75('-'),/,18X,                                     ADS85480
     1'ISTRAT = 5: AUGMENTED LAGRANGE MULTIPLIER METHOD',/,5X,75('-'))  ADS85490
310   FORMAT (/,5X,'-- INITIAL DESIGN')                                 ADS85500
320   FORMAT (/,5X,'SCALED OBJ =',1PE12.5,5X,'PENALT =',1PE12.5,/,5X,   ADS85510
     1'TRUE OBJ   =',1PE12.5)                                           ADS85520
330   FORMAT (//,5X,'-- BEGIN ALM ITERATION',I5,5X,'RP =',1PE12.5)      ADS85530
      END                                                               ADS85540
C --- SUBROUTINE ADS307 (NCON,BETA,CK,G,PPARAM,ALAMDA,IDG,RMULT,RPMAX)  ADS85550
C     DIMENSION G(NCON),ALAMDA(NCON),PPARAM(NCON),IDG(NCON)             ADS85560
C     ------------------------------------------------------------------ADS85570
C --- ROUTINE TO UPDATE THE LAGRANGE MULTIPLIER IN AUGMENTED LAGRANGE   ADS85580
C     METHOD.                                                           ADS85590
C     ------------------------------------------------------------------ADS85600
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 NOV. 1983   ADS85610
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS85620
C --- INPUT                                                             ADS85630
C      NCON  - NUMBER OF CONSTRAINTS.                                   ADS85640
C      BETA,CK - PARAMETERS TO UPDATE LAGRANGE MULTIPLIER.              ADS85650
C      G     - CONSTRAINT VALUES.                                       ADS85660
C      PPARAM - PENALTY PARAMETER OF EACH CONSTRAINT.                   ADS85670
C      ALAMDA - LAGRANGE MULTIPLIER.                                    ADS85680
C      IDG    - CONSTRAINT TYPE IDENTIFIER.                             ADS85690
C      RMULT  - FACTOR FOR PPARAM.                                      ADS85700
C      RPMAX  - MAXIMUM VALUE OF PPARAM.                                ADS85710
C --- OUTPUT                                                            ADS85720
C      PPARAM - UPDATED PENALTY PARAMETER.                              ADS85730
C      ALAMDA - UPDATED LAGRANGE PARAMETER.                             ADS85740
C                                                                       ADS85750
      SUBROUTINE ADS307 (NCON,BETA,CK,G,PPARAM,ALAMDA,IDG,RMULT,RPMAX)  ADS85760
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),ALAMDA(1),PPARAM(1),IDG(1)                         ADS85770
      CRIT=0.001                                                        ADS85780
      BCK=BETA*CK                                                       ADS85790
      CKM=0.                                                            ADS85800
      CKN=0.                                                            ADS85810
      DO 20 I=1,NCON                                                    ADS85820
      IF (IDG(I).LT.0) GO TO 10                                         ADS85830
      GMAX=DMAX1(G(I),-ALAMDA(I)/PPARAM(I))                             ADS85840
      GMAX=ABS(GMAX)                                                    ADS85850
      IF (GMAX.GT.CKM) CKM=GMAX                                         ADS85860
      GO TO 20                                                          ADS85870
10    CONTINUE                                                          ADS85880
      HMAX=ABS(G(I))                                                    ADS85890
      IF (HMAX.GT.CKN) CKN=HMAX                                         ADS85900
20    CONTINUE                                                          ADS85910
      IF (CKM.GT.CK.OR.CKN.GT.CK) GO TO 50                              ADS85920
      CK=DMAX1(CKM,CKN)                                                 ADS85930
      DO 40 I=1,NCON                                                    ADS85940
      AGI=ABS(G(I))                                                     ADS85950
      GI=G(I)                                                           ADS85960
      IF (AGI.GE.CRIT) GO TO 30                                         ADS85970
      GIC=GI/CRIT                                                       ADS85980
      GI=CRIT*GIC*GIC*GIC*(2.-GIC*GIC)                                  ADS85990
30    CONTINUE                                                          ADS86000
      ALAM=ALAMDA(I)+PPARAM(I)*GI                                       ADS86010
      IF (IDG(I).GE.0.AND.ALAM.LT.0.0) ALAM=0.                          ADS86020
      ALAMDA(I)=ALAM                                                    ADS86030
40    CONTINUE                                                          ADS86040
C                                                                       ADS86050
50    CONTINUE                                                          ADS86060
      DO 80 I=1,NCON                                                    ADS86070
      IF (IDG(I).LT.0) GO TO 60                                         ADS86080
      GMAX=DMAX1(G(I),-ALAMDA(I)/PPARAM(I))                             ADS86090
      GMAX=ABS(GMAX)                                                    ADS86100
      GO TO 70                                                          ADS86110
60    CONTINUE                                                          ADS86120
      GMAX=ABS(G(I))                                                    ADS86130
70    CONTINUE                                                          ADS86140
      IF (GMAX.GT.BCK.OR.G(I).GE.0.01) PPARAM(I)=RMULT*PPARAM(I)        ADS86150
      IF (PPARAM(I).GT.RPMAX) PPARAM(I)=RPMAX                           ADS86160
80    CONTINUE                                                          ADS86170
      RETURN                                                            ADS86180
      END                                                               ADS86190
C --- SUBROUTINE ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                     ADS86200
C     DIMENSION IDG(NCON),PPARAM(NCON)                                  ADS86210
C     ------------------------------------------------------------------ADS86220
C --- ROUTINE TO CALCULATE THE PENALTY PARAMETER(PPARAM).               ADS86230
C     ------------------------------------------------------------------ADS86240
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 NOV. 1982   ADS86250
C                                                MODIFIED 1983/JULY.    ADS86260
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS86270
C --- INPUT                                                             ADS86280
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS86290
C      IDG    - CONSTRAINT TYPE IDENTIFIER.                             ADS86300
C      RPPRIM - PENALTY PARAMETER FOR INEQUALITY CONSTRAINTS.           ADS86310
C      RP     - PENALTY PARAMETER FOR EQUALITY CONSTRAINTS.             ADS86320
C --- OUTPUT                                                            ADS86330
C      PPARAM - PENALTY PARAMETER OF EACH CONSTRAINT.                   ADS86340
C                                                                       ADS86350
      SUBROUTINE ADS308 (NCON,IDG,RPPRIM,RP,PPARAM)                     ADS86360
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION IDG(1),PPARAM(1)                                        ADS86370
      DO 10 I=1,NCON                                                    ADS86380
      PQR=RPPRIM                                                        ADS86390
      IF (IDG(I).LT.0) PQR=RP                                           ADS86400
      PPARAM(I)=PQR                                                     ADS86410
10    CONTINUE                                                          ADS86420
      RETURN                                                            ADS86430
      END                                                               ADS86440
C --- SUBROUTINE ADS310 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,ADS86450
C    1A,NCOLA,NRA,B,X,IDG,WK,IWK,ICQP,SVLB,SVUB,S,RAG,GMULT,XOLD,P,Y,   ADS86460
C    2ZRO,IPRNT3,PSAIZ,SCFX,XSUB,DELSTR,DABSTR,ITRMST,JTMAX,            ADS86470
C    3ISIDC,GGL,XTMP,SCFO,ITMP,BTMP,NRB,THETA,CT,CTMIN,CTL,CTLMIN,NAC,  ADS86480
C    4NVC,C,GMULTZ,JONED,DFSUB,STOL)                                    ADS86490
C     DIMENSION G(NCON),IC(NCOLA),DF(NDV),A(NRA,NCOLA),B(NDV*(NDV+1)/2),ADS86500
C    1X(NDV),IDG(NCON),WK(NRWK),IWK(NRIWK),ICQP(NCOLA),SVLB(NDV),       ADS86510
C    2SVUB(NDV),S(NDV+1),RAG(NCOLA),GMULT(NCON),XOLD(NDV),P(NDV),Y(NDV),ADS86520
C    3SCFX(NDV),XSUB(NCOLA),ISIDC(NDV),GGL(NCOLA),XTMP(NDV),ITMP(NCOLA),ADS86530
C    4BTMP(NRB,NRB),THETA(NCOLA),C(NCOLA+NDV+1),DFSUB(NDV)              ADS86540
C     ------------------------------------------------------------------ADS86550
C --- ROUTINE FOR PERFORM THE VARIABLE METRIC METHOD WITH CONSTRAINTS.  ADS86560
C     ------------------------------------------------------------------ADS86570
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                MARCH 1983   ADS86580
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS86590
C --- INPUT                                                             ADS86600
C     INFO  - CONTROL PARAMETER.                                        ADS86610
C     IOPT  - OPTIMIZER IDENTIFIER.                                     ADS86620
C     IONED - ONE DIMENSIONAL SEARCH IDENTIFIER.                        ADS86630
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS86640
C     NCON  - NUMBER OF CONSTRAINTS.                                    ADS86650
C     OBJ   - VALUE OF PRIMAL OBJECTIVE FUNCTION.                       ADS86660
C     G     - CURRENT VALUES OF CONSTRAINTS.                            ADS86670
C     NGT   - TOTAL NUMBER OF GRADIENTS CALCULATED.                     ADS86680
C     IC    - ARRAY IDENTIFYING CONSTRAINTS THAT ARE ACTIVE.            ADS86690
C     DF    - GRADIENT OF PRIMAL OBJECTIVE FUNCTION.                    ADS86700
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS86710
C             PROBLEM.                                                  ADS86720
C     NCOLA - MAXIMUM NUMBER OF ACTIVE CONSTRAINTS.                     ADS86730
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS86740
C     B     - HESSIAN MATRIX.                                           ADS86750
C     X     - INITIAL VALLUE FOR DESIGN VARIABLES.                      ADS86760
C     IDG   - CONSTRAINT TYPE IDENTIFIER.                               ADS86770
C     WK    - REAL WORK ARRAY.                                          ADS86780
C     IWK   - INTEGER WORK ARRAY.                                       ADS86790
C     ICQP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR Q.P. SUB-PROBLEM.ADS86800
C     SVLB,SVUB - UPPER AND LOWER LIMITS FOR PRIMAL DESIGN VARIABLES.   ADS86810
C     S     - WORK ARRAY (NDV).                                         ADS86820
C     RAG   - LAGRANGE MULTIPLIER.                                      ADS86830
C     GMULT - PENALTY PARAMETER.                                        ADS86840
C     XOLD  - ARRAY FOR SAVING THE OLD 'X'.                             ADS86850
C     P     - TRUE DIFFERENCE BETWEEN FORMER AND CURRENT GRADIENTS OF   ADS86860
C             OBJECTIVE FUNCTION.                                       ADS86870
C     Y     - DIFFERENCE BETWEEN FORMER AND CURRENT DESIGN VARIABLES.   ADS86880
C     ZRO   - VALUE OF ZERO.                                            ADS86890
C     IPRNT3 - PRINT CONTROL FOR STRATEGY LEVEL.                        ADS86900
C     PSAIZ - THE FACTOR TO AVOID INCONSISTENCE BETWEEN CONSTRAINTS.    ADS86910
C     SCFX  - SCALE FACTOR FOR DESIGN VARIABLES.                        ADS86920
C     XSUB  - DESIGN VARIABLES OF SUB-PROBLEM.                          ADS86930
C     DELSTR - CONVERGENCE CRITERIOR OF RELATIVE CHANGE.                ADS86940
C              0.0010 IS RECOMMENDED.                                   ADS86950
C     DABSTR - CONVERGENCE CRITERIOR OF ABSOLUTE CHANGE.                ADS86960
C     ITRMST - MINIMUM NUMBER OF CONVERGENCE CRITERIOR SATISFACTION.    ADS86970
C     JONED  - ONE D.-SEARCH IDENTIFIER FOR PSEUDO-OBJECTIVE.           ADS86980
C     JTMAX  - MAXIMUM NUMBER OF ITERATION.                             ADS86990
C     SCFO   - SCALE FACTOR OF OBJECTIVE.                               ADS87000
C     ITMP,BTMP,THETA - WORK ARRAY FOR ADS217.                          ADS87010
C     NRB - NUMBER OF ROWS IN ARRAY BTMP.                               ADS87020
C     GMULTZ - INITIAL VALUE OF PENALTY PARAMETER.                      ADS87030
C     DFSUB  - GRADIENTS OF OBJECTIVE OF SUB-PROBLEM.                   ADS87040
C --- OUTPUT                                                            ADS87050
C     X     - OPTIMUM VALUE FOR DESIGN VARIABLES.                       ADS87060
C                                                                       ADS87070
      SUBROUTINE ADS310 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,ADS87080
     1A,NCOLA,NRA,B,X,IDG,WK,IWK,ICQP,SVLB,SVUB,S,RAG,GMULT,XOLD,P,Y,ZROADS87090
     2,IPRNT3,PSAIZ,SCFX,XSUB,DELSTR,DABSTR,ITRMST,JTMAX,ISIDC,GGL,XTMP,ADS87100
     3SCFO,ITMP,BTMP,NRB,THETA,CT,NAC,NVC,C,JONED,                      ADS87110
     4DFSUB,VLBST,VUBST,IDAB,IDEL,ISLOP,JTER,LGOTO,NACT,NAVC,TEMP2,     ADS87120
     5SCFG)                                                             ADS87130
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IC(1),DF(1),A(NRA,1),B(1),X(1),IDG(1),WK(1),IWK(1),ADS87140
     1ICQP(1),SVLB(1),SVUB(1),S(1),RAG(1),GMULT(1),XOLD(1),P(1),Y(1),SCFADS87150
     2X(1),XSUB(1),ISIDC(1),GGL(1),XTMP(1),ITMP(1),BTMP(NRB,1),THETA(1),ADS87160
     3C(1),DFSUB(1),VLBST(1),VUBST(1),TEMP2(10),SCFG(1)                 ADS87170
C                                                                       ADS87180
C     GET INFO OUT OF TEMP2.                                            ADS87190
      CTMIN=TEMP2(1)                                                    ADS87200
      CTL=TEMP2(2)                                                      ADS87210
      CTLMIN=TEMP2(3)                                                   ADS87220
      GMULTZ=TEMP2(4)                                                   ADS87230
      ALPHA=TEMP2(5)                                                    ADS87240
      ALPHB=TEMP2(6)                                                    ADS87250
      ALPSAV=TEMP2(7)                                                   ADS87260
      OBJ1=TEMP2(8)                                                     ADS87270
      SLOPF=TEMP2(9)                                                    ADS87280
      STOL=TEMP2(10)                                                    ADS87290
      IGOTO1=IGOTO+1                                                    ADS87300
      GO TO (10,220,290,370), IGOTO1                                    ADS87310
10    CONTINUE                                                          ADS87320
      IF (IPRNT3.GT.0) WRITE (6,600)                                    ADS87330
      JTER=0                                                            ADS87340
      IDEL=0                                                            ADS87350
      IDAB=0                                                            ADS87360
      NAVC=NAC+NVC                                                      ADS87370
C     INITIALLY HESSIAN IS UNIT MATRIX.                                 ADS87380
      CALL ADS506 (NDV,B)                                               ADS87390
      DO 20 J=1,NCON                                                    ADS87400
      GMULT(J)=GMULTZ                                                   ADS87410
20    CONTINUE                                                          ADS87420
C     INITIAL FUNCTION EVALUATION IS AVAILABLE.                         ADS87430
      IF (IPRNT3.EQ.0) GO TO 30                                         ADS87440
      WRITE (6,610)                                                     ADS87450
      OBJW=OBJ/SCFO                                                     ADS87460
      WRITE (6,620) OBJW                                                ADS87470
      ICALC=1                                                           ADS87480
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS87490
      ICALC=4                                                           ADS87500
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS87510
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS87520
30    CONTINUE                                                          ADS87530
C     ------------------------------------------------------------------ADS87540
C     ITERATION START.                                                  ADS87550
C     ------------------------------------------------------------------ADS87560
40    CONTINUE                                                          ADS87570
      JTER=JTER+1                                                       ADS87580
      IF (IPRNT3.GT.0) WRITE (6,630) JTER                               ADS87590
      IF (IPRNT3.GT.3) WRITE (6,670) NGT                                ADS87600
      IF (IPRNT3.GT.3) WRITE (6,680) (IC(I),I=1,NGT)                    ADS87610
      OBJ1=OBJ                                                          ADS87620
      ICALC=3                                                           ADS87630
      IF (IPRNT3.LT.5) GO TO 60                                         ADS87640
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS87650
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS87660
      ISCAL=1                                                           ADS87670
      DO 50 J=1,NGT                                                     ADS87680
      WRITE (6,420) IC(J)                                               ADS87690
      CALL ADS401 (A(1,J),NDV,SCFG,ISCAL)                               ADS87700
50    CONTINUE                                                          ADS87710
60    CONTINUE                                                          ADS87720
      NBOUT=(NDV+4)/5                                                   ADS87730
      IF (IPRNT3.GT.4) CALL ADS404 (NDV,B)                              ADS87740
C     ----------------------------------------------------------------- ADS87750
C     CHECK CONVERGENCE BY RESIDUAL OF KUHN TUCKER CONDITION.           ADS87760
C     ----------------------------------------------------------------- ADS87770
      CTTMP=-CTMIN                                                      ADS87780
      CTLTMP=-CTLMIN                                                    ADS87790
C     DETERMINE ACTIVE/VIOLATED CONSTRAINTS SET.                        ADS87800
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CTTMP,CTMIN,CTLTMP,CTLMIN,NACTADS87810
     1,NVCT)                                                            ADS87820
      NAVCT=NACT+NVCT                                                   ADS87830
      IF (NVCT.GT.0) GO TO 70                                           ADS87840
C     IDENTIFY ACTIVE SIDE CONSTRAINTS.                                 ADS87850
      CALL ADS509 (NDV,X,SVLB,SVUB,ZRO,ISIDC,NACS)                      ADS87860
      IF(NAVCT.EQ.0) GO TO 61                                           ADS87870
      CALL ADS207 (NDV,NAVCT,A,NRA,ZRO,GGL)                             ADS87880
61    N2=1                                                              ADS87890
      CALL ADS207 (NDV,N2,DF,NDV,ZRO,SCALF)                             ADS87900
      IF(NACS.EQ.0.AND.NAVCT.EQ.0) GO TO 70                             ADS87910
C     CALCULATE THE RESIDUAL OF KUHN-TUCKER CONDITION.                  ADS87920
      THZ=0.                                                            ADS87930
      DUM=1.                                                            ADS87940
      NDUM=0                                                            ADS87950
      CALL ADS217 (NDV,NAVCT,NDUM,NACS,DUM,DUM,DUM,THZ,ZRO,DUM,GGL,G,DF,ADS87960
     1A,NRA,IDG,IC,ISIDC,ITMP,BTMP,NRB,SCALF,S,THETA,RAG,BETA,RSDMX)    ADS87970
      IF (IPRNT3.GE.4) WRITE (6,560) RSDMX                              ADS87980
      IF (RSDMX.LT.0.001) GO TO 300                                     ADS87990
70    CONTINUE                                                          ADS88000
C     ------------------------------------------------------------------ADS88010
C     SOLVE Q.P. SUB-PROBLEM.                                           ADS88020
C     ------------------------------------------------------------------ADS88030
      DO 80 K=1,NGT                                                     ADS88040
      J=IC(K)                                                           ADS88050
      Y(K)=G(J)                                                         ADS88060
      IF (IDG(J).LT.0.OR.G(J).GE.0.) Y(K)=PSAIZ*G(J)                    ADS88070
      ITMP(K)=IDG(J)                                                    ADS88080
80    CONTINUE                                                          ADS88090
      CALL ADS503 (NDV,SVLB,VLBST)                                      ADS88100
      CALL ADS503 (NDV,SVUB,VUBST)                                      ADS88110
      CALL ADS371 (IOPT,IONED,NDV,NGT,X,XSUB,A,NRA,B,DF,Y,ICQP,RAG,DFSUBADS88120
     1,SVLB,SVUB,ITMP,NCOLA,WK,IWK,P)                                   ADS88130
      CALL ADS503 (NDV,VLBST,SVLB)                                      ADS88140
      CALL ADS503 (NDV,VUBST,SVUB)                                      ADS88150
      IF (IPRNT3.GT.3) WRITE (6,650)                                    ADS88160
      ICALC=2                                                           ADS88170
      IF (IPRNT3.GT.3) CALL ADS402 (ICALC,XSUB,NDV,SCFX)                ADS88180
C     IF ALL S(I) ARE EQUAL TO ZERO , THEN REGARD AS CONVERGENCE.       ADS88190
      SMX=0.                                                            ADS88200
      DO 90 I=1,NDV                                                     ADS88210
      SS=ABS(XSUB(I))                                                   ADS88220
      IF (SS.GE.SMX) SMX=SS                                             ADS88230
90    CONTINUE                                                          ADS88240
      IF (SMX.LE.STOL) GO TO 300                                        ADS88250
C     ----------------------------------------------------------------- ADS88260
C     DETERMINE LAGRANGE MULTIPLIER.                                    ADS88270
C     ----------------------------------------------------------------- ADS88280
C     ONLY ACTIVE CONSTRAINTS OF Q.P. SUB-PROBLEM ARE USED TO CALCULATE ADS88290
C     THE LAGRANGE MULTIPLIER.                                          ADS88300
      CALL ADS009 (NDV,A,NRA,NGT,ICQP,ITMP,RAG,CT,CTMIN,CTL,CTLMIN,NACQ,ADS88310
     1NVCQ)                                                             ADS88320
      NAVCQ=NACQ+NVCQ                                                   ADS88330
      IF (IPRNT3.GT.3.AND.NAVCQ.EQ.0) WRITE (6,690)                     ADS88340
      IF (IPRNT3.GT.3.AND.NAVCQ.GT.0) WRITE (6,700) NAVCQ               ADS88350
      IF (IPRNT3.GT.3.AND.NAVCQ.GT.0) WRITE (6,710) (IC(ICQP(I)),I=1,NAVADS88360
     1CQ)                                                               ADS88370
      DO 100 I=1,NGT                                                    ADS88380
      RAG(I)=0.                                                         ADS88390
100   CONTINUE                                                          ADS88400
      IF (NAVCQ.EQ.0) GO TO 130                                         ADS88410
      DO 110 I=1,NDV                                                    ADS88420
      XTMP(I)=X(I)+XSUB(I)                                              ADS88430
110   CONTINUE                                                          ADS88440
C     IDENTIFY ACTIVE SIDE CONSTRAINTS.                                 ADS88450
      CALL ADS509 (NDV,XTMP,SVLB,SVUB,ZRO,ISIDC,NACS)                   ADS88460
      CALL ADS207 (NDV,NAVCQ,A,NRA,ZRO,GGL)                             ADS88470
      N2=1                                                              ADS88480
      CALL ADS207 (NDV,N2,DFSUB,NDV,ZRO,SCALF)                          ADS88490
C     CALCULATE LAGRANGE MULTIPLIER.                                    ADS88500
      THZ=0.                                                            ADS88510
      DUM=1.                                                            ADS88520
      NDUM=0                                                            ADS88530
      CALL ADS217 (NDV,NAVCQ,NDUM,NACS,DUM,DUM,DUM,THZ,ZRO,DUM,GGL,G,DFSADS88540
     1UB,A,NRA,IDG,ICQP,ISIDC,ITMP,BTMP,NRB,SCALF,S,THETA,C,BETA,RSDMX) ADS88550
      DO 120 K=1,NAVCQ                                                  ADS88560
      J=ICQP(K)                                                         ADS88570
      RAG(J)=C(K)                                                       ADS88580
120   CONTINUE                                                          ADS88590
130   CONTINUE                                                          ADS88600
      ICALC=10                                                          ADS88610
      IF (IPRNT3.GT.3) CALL ADS402 (ICALC,RAG,NGT,SCFX)                 ADS88620
C     CALCULAT                                                          ADS88630
C     REVISION OF HESSIAN MATRIX.                                       ADS88640
      DO 150 I=1,NDV                                                    ADS88650
      Y(I)=DF(I)                                                        ADS88660
      IF (NAVCQ.EQ.0) GO TO 150                                         ADS88670
      DO 140 KK=1,NAVCQ                                                 ADS88680
      K=ICQP(KK)                                                        ADS88690
      Y(I)=Y(I)+A(I,KK)*RAG(K)                                          ADS88700
140   CONTINUE                                                          ADS88710
150   CONTINUE                                                          ADS88720
C     ------------------------------------------------------------------ADS88730
C     CALCULATION OF 'GMULT' FOR ONE-D. SEARCH.                         ADS88740
C     ------------------------------------------------------------------ADS88750
      DO 180 J=1,NCON                                                   ADS88760
      DO 160 I=1,NAVCQ                                                  ADS88770
      II=ICQP(I)                                                        ADS88780
      IF (J.EQ.IC(II)) GO TO 170                                        ADS88790
160   CONTINUE                                                          ADS88800
      GMULT(J)=0.5*GMULT(J)                                             ADS88810
      GO TO 180                                                         ADS88820
170   CONTINUE                                                          ADS88830
      GMULT(J)=DMAX1(GMULT(J),RAG(II))                                  ADS88840
180   CONTINUE                                                          ADS88850
      ICALC=11                                                          ADS88860
      IF (IPRNT3.GT.3) CALL ADS402 (ICALC,GMULT,NCON,SCFX)              ADS88870
C     ------------------------------------------------------------------ADS88880
C     SAVE OLD DESIGN VARIABLES AND CALCULATE THE GRADIENT OF TEMPORARY ADS88890
C     OBJECTIVE FUNCTION.                                               ADS88900
C     ------------------------------------------------------------------ADS88910
      ISLOP=0                                                           ADS88920
      CALL ADS503 (NDV,X,XOLD)                                          ADS88930
C     CALCULATE GRADIENT OF PSEUDO-OBJECTIVE.                           ADS88940
      CALL ADS376 (NDV,NGT,G,IDG,IC,ICQP,DF,GMULT,A,NRA,S)              ADS88950
      ICALC=12                                                          ADS88960
      IF (IPRNT3.GT.3) CALL ADS402 (ICALC,S,NDV,SCFX)                   ADS88970
C     ------------------------------------------------------------------ADS88980
C     SLOPE OF TEMPORARY OBJECTIVE FUNCTION.                            ADS88990
C     ------------------------------------------------------------------ADS89000
      CALL ADS501 (NDV,S,XSUB,SLOPF)                                    ADS89010
      ZROM=-ZRO                                                         ADS89020
      IF (SLOPF.LT.ZROM) GO TO 190                                      ADS89030
      IF (GMAX.GE.0.0) GO TO 190                                        ADS89040
C     STEEPEST DESCENT METHOD.                                          ADS89050
      IF (IPRNT3.GT.3) WRITE (6,590) SLOPF                              ADS89060
      ICONT=3                                                           ADS89070
      SCALT=0.                                                          ADS89080
      CALL ADS509 (NDV,X,SVLB,SVUB,ZRO,ISIDC,NACS)                      ADS89090
      DUM=1.                                                            ADS89100
      CALL ADS236 (NDV,ICONT,S,XTMP,XTMP,B,DUM,XSUB,DM2,SCALT,ZRO,ISIDC)ADS89110
      CALL ADS501 (NDV,S,XSUB,SLOPF)                                    ADS89120
190   CONTINUE                                                          ADS89130
      WK(40)=1.0                                                        ADS89140
      IF (JTER.GT.1) GO TO 200                                          ADS89150
      ALPHA=0.10                                                        ADS89160
      ALPSAV=ALPHA                                                      ADS89170
200   CONTINUE                                                          ADS89180
      ALP1=ALPSAV                                                       ADS89190
      ALPSAV=ALPHA                                                      ADS89200
      ALPHA=0.5*(ALPHA+ALP1)                                            ADS89210
      ALPHB=ALPHA                                                       ADS89220
      IF (ALPHA.LT.ZRO) ALPHA=ZRO                                       ADS89230
C     ALPMAX.                                                           ADS89240
      CALL ADS112 (NDV,X,S,SVLB,SVUB,ZRO,ALPMAX)                        ADS89250
      IF (ALPMAX.LE.ZRO) GO TO 250                                      ADS89260
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX                                 ADS89270
      DO 210 I=1,NDV                                                    ADS89280
      S(I)=XSUB(I)                                                      ADS89290
210   CONTINUE                                                          ADS89300
      LGOTO=0                                                           ADS89310
      IF (IPRNT3.GT.3) WRITE (6,570) ALPHA                              ADS89320
C     ------------------------------------------------------------------ADS89330
C     DO ONE-DIMENSIONAL SEARCH                                         ADS89340
C     ------------------------------------------------------------------ADS89350
220   CONTINUE                                                          ADS89360
      CALL ADS375 (NCON,IDG,G,GMULT,PENALT)                             ADS89370
      FAI=OBJ+PENALT                                                    ADS89380
      DO 230 J=1,NCON                                                   ADS89390
230   G(J)=G(J)-0.20                                                    ADS89400
      CALL ADS001 (LGOTO,ISTRAT,JONED,NDV,NCON,FAI,X,SVLB,SVUB,G,IDG,   ADS89410
     1NAVC,IC,WK,IWK,SLOPF,ALPHA,ALPMAX)                                ADS89420
      DO 240 J=1,NCON                                                   ADS89430
240   G(J)=G(J)+0.20                                                    ADS89440
      IF (LGOTO.EQ.0) GO TO 250                                         ADS89450
C     EVALUATE OBJECTIVE.                                               ADS89460
      IGOTO=1                                                           ADS89470
      INFO=1                                                            ADS89480
C     PUT INFO IN TEMP2 FOR TRANSFER.                                   ADS89490
      TEMP2(1)=CTMIN                                                    ADS89500
      TEMP2(2)=CTL                                                      ADS89510
      TEMP2(3)=CTLMIN                                                   ADS89520
      TEMP2(4)=GMULTZ                                                   ADS89530
      TEMP2(5)=ALPHA                                                    ADS89540
      TEMP2(6)=ALPHB                                                    ADS89550
      TEMP2(7)=ALPSAV                                                   ADS89560
      TEMP2(8)=OBJ1                                                     ADS89570
      TEMP2(9)=SLOPF                                                    ADS89580
      TEMP2(10)=STOL                                                    ADS89590
      RETURN                                                            ADS89600
250   CONTINUE                                                          ADS89610
      IF (IPRNT3.GT.3) WRITE (6,640) ALPHA                              ADS89620
      IF (ALPHA.LE.ZRO) GO TO 260                                       ADS89630
      CALL ADS375 (NCON,IDG,G,GMULT,PENALT)                             ADS89640
      OBJ=FAI-PENALT                                                    ADS89650
      GO TO 290                                                         ADS89660
260   CONTINUE                                                          ADS89670
      IF (ISLOP.EQ.1) GO TO 270                                         ADS89680
      ISLOP=1                                                           ADS89690
      ALPHA=0.05*ALPHB                                                  ADS89700
      IF (ALPHA.LT.ZRO) ALPHA=ZRO                                       ADS89710
      WK(39)=ALPHA                                                      ADS89720
      LGOTO=0                                                           ADS89730
      GO TO 220                                                         ADS89740
270   CONTINUE                                                          ADS89750
      ALPHA=0.01                                                        ADS89760
      DO 280 I=1,NDV                                                    ADS89770
      X(I)=XOLD(I)-ALPHA*S(I)                                           ADS89780
280   CONTINUE                                                          ADS89790
C     ------------------------------------------------------------------ADS89800
C     CALCULATE THE VALUES OF OBJECTIVE FUNCTION AND CONSTRAINTS FOR    ADS89810
C     NEW DESIGN VARIABLES.                                             ADS89820
C     ------------------------------------------------------------------ADS89830
      IGOTO=2                                                           ADS89840
      INFO=1                                                            ADS89850
C     PUT INFO IN TEMP2 FOR TRANSFER.                                   ADS89860
      TEMP2(1)=CTMIN                                                    ADS89870
      TEMP2(2)=CTL                                                      ADS89880
      TEMP2(3)=CTLMIN                                                   ADS89890
      TEMP2(4)=GMULTZ                                                   ADS89900
      TEMP2(5)=ALPHA                                                    ADS89910
      TEMP2(6)=ALPHB                                                    ADS89920
      TEMP2(7)=ALPSAV                                                   ADS89930
      TEMP2(8)=OBJ1                                                     ADS89940
      TEMP2(9)=SLOPF                                                    ADS89950
      TEMP2(10)=STOL                                                    ADS89960
      RETURN                                                            ADS89970
290   CONTINUE                                                          ADS89980
      OBJW=OBJ/SCFO                                                     ADS89990
      IF (IPRNT3.GT.1) WRITE (6,620) OBJW                               ADS90000
      ICALC=1                                                           ADS90010
      IF (IPRNT3.GT.1) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS90020
      ICALC=4                                                           ADS90030
      IF (IPRNT3.GT.2) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS90040
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS90050
      IF (IPRNT3.GT.3) WRITE (6,660) GMAX                               ADS90060
C     ----------------------------------------------------------------- ADS90070
C     CHECK FOR CONVERGENCE.                                            ADS90080
C     ----------------------------------------------------------------- ADS90090
      CALL ADS211 (OBJ,OBJ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)          ADS90100
      IF (JTER.GE.JTMAX) GO TO 300                                      ADS90110
      IF (IDEL.LT.ITRMST.AND.IDAB.LT.ITRMST) GO TO 350                  ADS90120
C     OPTIMIZATION IS COMPLETE.                                         ADS90130
300   CONTINUE                                                          ADS90140
      IF (IPRNT3.EQ.0) GO TO 340                                        ADS90150
      OBJW=OBJ/SCFO                                                     ADS90160
      WRITE (6,430) JTER,OBJW                                           ADS90170
C     X-VECTOR.                                                         ADS90180
      ICALC=1                                                           ADS90190
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS90200
C     G-VECTOR.                                                         ADS90210
      ICALC=4                                                           ADS90220
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS90230
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS90240
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS90250
C     SIDE CONSTRAINTS.                                                 ADS90260
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS90270
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS90280
      WRITE (6,440) CT,CTL                                              ADS90290
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS90300
      NAVC=NAC+NVC                                                      ADS90310
      WRITE (6,510) NAC,NVC                                             ADS90320
      IF (NAVC.GT.0) WRITE (6,520)                                      ADS90330
      IF (NAVC.GT.0) WRITE (6,530) (ITMP(I),I=1,NAVC)                   ADS90340
      WRITE (6,540) NACS                                                ADS90350
      IF (NACS.EQ.0) GO TO 330                                          ADS90360
      N=0                                                               ADS90370
      DO 320 I=1,NDV                                                    ADS90380
      IF (ISIDC(I).EQ.0) GO TO 320                                      ADS90390
      IF (ISIDC(I).EQ.1) GO TO 310                                      ADS90400
C     LOWER BOUND IS ACTIVE.                                            ADS90410
      N=N+1                                                             ADS90420
      ITMP(N)=-I                                                        ADS90430
310   IF (ISIDC(I).LT.0) GO TO 320                                      ADS90440
C     UPPER BOUND IS ACTIVE.                                            ADS90450
      N=N+1                                                             ADS90460
      ITMP(N)=I                                                         ADS90470
320   CONTINUE                                                          ADS90480
      WRITE (6,550)                                                     ADS90490
      WRITE (6,530) (ITMP(I),I=1,NACS)                                  ADS90500
330   WRITE (6,450)                                                     ADS90510
      IF (JTER.GE.JTMAX) WRITE (6,460)                                  ADS90520
      IF (IDEL.GE.ITRMST) WRITE (6,470) ITRMST                          ADS90530
      IF (IDAB.GE.ITRMST) WRITE (6,480) ITRMST                          ADS90540
      IF (RSDMX.LE.0.001) WRITE (6,490) RSDMX                           ADS90550
      IF (SMX.LE.STOL) WRITE (6,500) SMX,STOL                           ADS90560
340   CONTINUE                                                          ADS90570
      INFO=0                                                            ADS90580
      IGOTO=0                                                           ADS90590
      RETURN                                                            ADS90600
350   CONTINUE                                                          ADS90610
C     SAVE THE INFORMATION ABOUT ACTIVE CONSTRAINTS SET.                ADS90620
      DO 360 J=1,NGT                                                    ADS90630
      ICQP(J)=IC(J)                                                     ADS90640
360   CONTINUE                                                          ADS90650
C     ------------------------------------------------------------------ADS90660
C     DETERMINE THE ACTIVE CONSTRAINTS SET.                             ADS90670
C     ------------------------------------------------------------------ADS90680
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS90690
      NAVC=NAC+NVC                                                      ADS90700
C     CALCULATE GRADIENTS OF OBJECTIVE FUNCTION AND CONSTRAINTS.        ADS90710
      IGOTO=3                                                           ADS90720
      INFO=2                                                            ADS90730
C     PUT INFO IN TEMP2 FOR TRANSFER.                                   ADS90740
      TEMP2(1)=CTMIN                                                    ADS90750
      TEMP2(2)=CTL                                                      ADS90760
      TEMP2(3)=CTLMIN                                                   ADS90770
      TEMP2(4)=GMULTZ                                                   ADS90780
      TEMP2(5)=ALPHA                                                    ADS90790
      TEMP2(6)=ALPHB                                                    ADS90800
      TEMP2(7)=ALPSAV                                                   ADS90810
      TEMP2(8)=OBJ1                                                     ADS90820
      TEMP2(9)=SLOPF                                                    ADS90830
      TEMP2(10)=STOL                                                    ADS90840
      RETURN                                                            ADS90850
370   CONTINUE                                                          ADS90860
C     ------------------------------------------------------------------ADS90870
C     CALCULATE CURRENT GRADIENTS OF TEMPORARY OBJECTIVE FUNCTION AND   ADS90880
C     THE DIFFERENCE OF THE VALUES.                                     ADS90890
C     ------------------------------------------------------------------ADS90900
      DO 410 I=1,NDV                                                    ADS90910
      P(I)=DF(I)                                                        ADS90920
      DO 400 J=1,NGT                                                    ADS90930
      JJ=J                                                              ADS90940
      DO 380 K=1,NGT                                                    ADS90950
      KK=K                                                              ADS90960
      IF (IC(J).EQ.ICQP(K)) GO TO 390                                   ADS90970
380   CONTINUE                                                          ADS90980
      GO TO 400                                                         ADS90990
390   P(I)=P(I)+A(I,JJ)*RAG(KK)                                         ADS91000
400   CONTINUE                                                          ADS91010
      P(I)=P(I)-Y(I)                                                    ADS91020
      Y(I)=X(I)-XOLD(I)                                                 ADS91030
410   CONTINUE                                                          ADS91040
C     CALCULATE THE TRANSFORMED DIFFERENCE.                             ADS91050
      CALL ADS377 (NDV,P,Y,B,S)                                         ADS91060
C     PERFORM THE REVISION OF HESSIAN MATRIX.                           ADS91070
      CALL ADS378 (NDV,P,Y,B,S,ZRO,ILL)                                 ADS91080
      IF (ILL.NE.0.AND.IPRNT3.GT.3) WRITE (6,580)                       ADS91090
      GO TO 40                                                          ADS91100
C     ------------------------------------------------------------------ADS91110
C     FORMATS.                                                          ADS91120
C     ------------------------------------------------------------------ADS91130
420   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS91140
430   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS91150
     1'NUMBER OF CONSTRAINED MINIMIZATIONS =',I5,//,5X,'OBJECTIVE =',   ADS91160
     21PE12.5)                                                          ADS91170
440   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',E12.5,'  CTL =',E12.5)  ADS91180
450   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS91190
460   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS91200
470   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS91210
     1' CONSECUTIVE ITERATIONS')                                        ADS91220
480   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS91230
     1' CONSECUTIVE ITERATIONS')                                        ADS91240
490   FORMAT (/,5X,'MAXIMUM K-T RESIDUAL =',E12.5,' IS LESS THAN 0.001')ADS91250
500   FORMAT (/,5X,'MAXIMUM S-VECTOR COMPONENT =',1PE12.5,              ADS91260
     1' IS LESS THAN',1PE12.5)                                          ADS91270
510   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS91280
     1' VIOLATED CONSTRAINTS')                                          ADS91290
520   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS91300
530   FORMAT (5X,10I5)                                                  ADS91310
540   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS91320
550   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS91330
560   FORMAT (/,5X,'MAXIMUM RESIDUAL OF KUHN TUCKER CONDITION =',       ADS91340
     11PE12.5)                                                          ADS91350
570   FORMAT (/,5X,'ESTIMATED ALPHA  =',F9.4)                           ADS91360
580   FORMAT (/,5X,'REVISION OF HESSIAN FAILED')                        ADS91370
590   FORMAT (/,5X,'STEEPEST DESCENT METHOD',5X,1PE12.4)                ADS91380
600   FORMAT ('1',4X,75('-'),/,20X,                                     ADS91390
     1'ISTRAT = 8: SEQUENTIAL QUADRATIC PROGRAMMING',/,5X,75('-'))      ADS91400
610   FORMAT (/,5X,'-- INITIAL DESIGN')                                 ADS91410
620   FORMAT (/,5X,'OBJ =',1PE12.5)                                     ADS91420
630   FORMAT (//,5X,'-- BEGIN SQP INTERATION',I5)                       ADS91430
640   FORMAT (/,5X,'CALCULATED ALPHA =',F9.4)                           ADS91440
650   FORMAT (/,5X,'Q.P. IS SOLVED')                                    ADS91450
660   FORMAT (/,5X,'GMAX =',1PE12.4)                                    ADS91460
670   FORMAT (/,5X,'RETAINED CONSTRAINTS (',I4,' )')                    ADS91470
680   FORMAT (5X,15I5)                                                  ADS91480
690   FORMAT (/,5X,'NO ACTIVE/VIOLATED CONSTRAINTS IN SUB-PROBLEM')     ADS91490
700   FORMAT (/,5X,'ACTIVE CONSTRAINTS OF SUB-PROBLEM (',I4,' )')       ADS91500
710   FORMAT (5X,15I5)                                                  ADS91510
      END                                                               ADS91520
C --- SUBROUTINE ADS320 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF ADS91530
C    1,A,NCOLA,NRA,X,IDG,WK,IWK,ICLP,SVLB,SVUB,S,GS,VLBST,VUBST,ZRO,    ADS91540
C    2IPRNT3,RMVLMZ,SCFX,XTMP,SCLMIN,DELSTR,DABSTR,ITRMST,JTMAX,        ADS91550
C    3SCFO,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC,XSUB,DFSUB)                 ADS91560
C     DIMENSION G(NCON),IC(NCOLA),DF(NDV),A(NRA,NCOLA),X(NDV),IDG(NCON),ADS91570
C    1 WK(NRWK),IWK(NRIWK),ICLP(NCOLA),SVLB(NDV),SVUB(NDV),S(NDV),      ADS91580
C    2 GS(NCOLA),VLBST(NDV),VUBST(NDV),SCFX(NDV),XTMP(NDV),ITMP(NCOLA), ADS91590
C    3 XSUB(NDV),DFSUB(NDV)                                             ADS91600
C     ------------------------------------------------------------------ADS91610
C --- ROUTINE FOR PERFORM SEQUENTIAL LINEAR PROGRAMMING.                ADS91620
C     ------------------------------------------------------------------ADS91630
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                APRIL 1983   ADS91640
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS91650
C --- INPUT                                                             ADS91660
C     INFO  - CONTROL PARAMETER.                                        ADS91670
C     IOPT  - OPTIMIZER IDENTIFIER.                                     ADS91680
C     IONED - ONE DIMENSIONAL SEARCH IDENTIFIER.                        ADS91690
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS91700
C     NCON  - NUMBER OF CONSTRAINTS.                                    ADS91710
C     OBJ   - VALUE OF PRIMAL OBJECTIVE FUNCTION.                       ADS91720
C     G     - CURRENT VALUES OF CONSTRAINTS.                            ADS91730
C     NGT   - TOTAL NUMBER OF GRADIENTS CALCULATED.                     ADS91740
C     IC    - ARRAY IDENTIFYING CONSTRAINTS THAT ARE ACTIVE.            ADS91750
C     DF    - GRADIENT OF PRIMAL OBJECTIVE FUNCTION.                    ADS91760
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS91770
C             PROBLEM.                                                  ADS91780
C     NCOLA - MAXIMUM NUMBER OF ACTIVE CONSTRAINTS.                     ADS91790
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS91800
C     X     - INITIAL VALLUE FOR DESIGN VARIABLES.                      ADS91810
C     IDG   - CONSTRAINT TYPE IDENTIFIER.                               ADS91820
C     WK    - REAL WORK ARRAY.                                          ADS91830
C     IWK   - INTEGER WORK ARRAY.                                       ADS91840
C     ICLP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR L.P. SUB-PROBLEM.ADS91850
C     SVLB,SVUB - UPPER AND LOWER LIMITS FOR PRIMAL DESIGN VARIABLES.   ADS91860
C     S     - WORK ARRAY (NDV).                                         ADS91870
C     VUBST,VLBST - ARRAY TO STORE VLB,VUB.                             ADS91880
C     GS    - ARRAY TO STORE THE CURRENT VALUE OF CONSTRAINTS.          ADS91890
C     ZRO   - VALUE OF ZERO.                                            ADS91900
C     IPRNT3 - PRINT CONTROL FOR STRATEGY LEVEL.                        ADS91910
C     RMVLMZ - INITIAL VALUE OF MOVE LIMIT.                             ADS91920
C     SCFX  - SCALING FACTOR. ARRAY FOR QP SUB-PROBLEM.                 ADS91930
C     XTMP  - ARRAY FOR WORK.                                           ADS91940
C     SCLMIN - MINIMUM VALUE OF SCALING FACTOR.                         ADS91950
C     DELSTR - CONVERGENCE CRITERIOR OF RELATIVE CHANGE.                ADS91960
C     DABSTR - CONVERGENCE CRITERIOR OF ABSOLUTE CHANGE.                ADS91970
C     ITRMST - MINIMUM NUMBER OF CONVERGENCE CRITERIOR SATISFACTION.    ADS91980
C     JTMAX  - MAXIMUM NUMBER OF ITERATION.                             ADS91990
C     SCFO   - SCALE FACTOR OF OBJECTIVE.                               ADS92000
C     XSUB   - DESIGN VARIABLES OF SUB-PROBLEM.                         ADS92010
C     DFSUB  - GRADIENTS OF OBJECTIVE OF SUB-PROBLEM.                   ADS92020
C --- OUTPUT                                                            ADS92030
C     X     - OPTIMUM VALUE FOR DESIGN VARIABLES.                       ADS92040
C                                                                       ADS92050
      SUBROUTINE ADS320 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,ADS92060
     1A,NCOLA,NRA,X,IDG,WK,IWK,ICLP,SVLB,SVUB,S,GS,VLBST,VUBST,ZRO,IPRNTADS92070
     23,RMVLMZ,SCFX,XTMP,SCLMIN,DELSTR,DABSTR,ITRMST,JTMAX,SCFO,ITMP,CT,ADS92080
     3CTMIN,CTL,CTLMIN,NAC,NVC,XSUB,DFSUB,OBJR,OBJ1,RMVLMT,IDAB,IDEL,IMVADS92090
     4CM,JTER,ISIDC,SCFG)                                               ADS92100
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IC(1),DF(1),A(NRA,1),X(1),IDG(1),WK(1),IWK(1),ICLP(ADS92110
     11),SVLB(1),SVUB(1),S(1),GS(1),VLBST(1),VUBST(1),SCFX(1),XTMP(1),ITADS92120
     2MP(1),XSUB(1),DFSUB(1),ISIDC(1),SCFG(1)                           ADS92130
C                                                                       ADS92140
      IGOTO1=IGOTO+1                                                    ADS92150
      GO TO (10,30,80), IGOTO1                                          ADS92160
10    CONTINUE                                                          ADS92170
      IF (IPRNT3.GT.0) WRITE (6,280)                                    ADS92180
      DELST=5.*DELSTR                                                   ADS92190
      DABST=5.*DABSTR                                                   ADS92200
      OBJR=1.0E20                                                       ADS92210
      OBJ1=OBJ                                                          ADS92220
      RMVLMT=RMVLMZ                                                     ADS92230
      JTER=0                                                            ADS92240
      IMVCH=0                                                           ADS92250
C     INITIAL FUNCTION EVALUATION SHOULD BE FINISHED.                   ADS92260
      IF (IPRNT3.GT.0) WRITE (6,290)                                    ADS92270
      OBJW=OBJ/SCFO                                                     ADS92280
      IF (IPRNT3.GT.0) WRITE (6,300) OBJW                               ADS92290
      ICALC=1                                                           ADS92300
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS92310
      ICALC=4                                                           ADS92320
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS92330
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS92340
C     ------------------------------------------------------------------ADS92350
C     ITERATION START                                                   ADS92360
C     ------------------------------------------------------------------ADS92370
20    CONTINUE                                                          ADS92380
      JTER=JTER+1                                                       ADS92390
      IF (IPRNT3.GT.0) WRITE (6,310) JTER,RMVLMT                        ADS92400
      IF (JTER.EQ.1) GO TO 30                                           ADS92410
C     DETERMINE THE ACTIVE CONSTRAINTS SET.                             ADS92420
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS92430
C     CALCULATE GRADIENTS OF OBJECTIVE FUNCTION AND CONSTRAINTS.        ADS92440
      IGOTO=1                                                           ADS92450
      INFO=2                                                            ADS92460
      RETURN                                                            ADS92470
C                                                                       ADS92480
30    CONTINUE                                                          ADS92490
      IF (IPRNT3.GT.3) WRITE (6,330) NGT                                ADS92500
      IF (IPRNT3.GT.3) WRITE (6,340) (IC(I),I=1,NGT)                    ADS92510
      DO 40 J=1,NGT                                                     ADS92520
      ICJ=IC(J)                                                         ADS92530
      GS(J)=G(ICJ)                                                      ADS92540
      ITMP(J)=IDG(ICJ)                                                  ADS92550
40    CONTINUE                                                          ADS92560
      ICALC=3                                                           ADS92570
      IF (IPRNT3.LT.5) GO TO 60                                         ADS92580
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS92590
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS92600
      ISCAL=1                                                           ADS92610
      DO 50 J=1,NGT                                                     ADS92620
      WRITE (6,160) IC(J)                                               ADS92630
      CALL ADS401 (A(1,J),NDV,SCFG,ISCAL)                               ADS92640
50    CONTINUE                                                          ADS92650
60    CONTINUE                                                          ADS92660
C     ------------------------------------------------------------------ADS92670
C     SOLVE L.P. SUB-PROBLEM.                                           ADS92680
C     ------------------------------------------------------------------ADS92690
      CALL ADS503 (NDV,SVLB,VLBST)                                      ADS92700
      CALL ADS503 (NDV,SVUB,VUBST)                                      ADS92710
      CALL ADS321 (IOPT,IONED,NDV,NGT,X,XSUB,A,NRA,DF,GS,ICLP,G,DFSUB,SVADS92720
     1LB,SVUB,ITMP,NCOLA,WK,IWK,RMVLMT)                                 ADS92730
      CALL ADS503 (NDV,VLBST,SVLB)                                      ADS92740
      CALL ADS503 (NDV,VUBST,SVUB)                                      ADS92750
      DO 70 I=1,NDV                                                     ADS92760
      X(I)=X(I)+XSUB(I)                                                 ADS92770
70    CONTINUE                                                          ADS92780
C     CALCULATE THE VALUES OF OBJECTIVE FUNCTION AND CONSTRAINTS FOR NEWADS92790
C     DESIGN VARIABLES.                                                 ADS92800
      IGOTO=2                                                           ADS92810
      INFO=1                                                            ADS92820
      RETURN                                                            ADS92830
80    CONTINUE                                                          ADS92840
      OBJW=OBJ/SCFO                                                     ADS92850
      IF (IPRNT3.GT.0) WRITE (6,300) OBJW                               ADS92860
      ICALC=1                                                           ADS92870
      IF (IPRNT3.GT.1) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS92880
      ICALC=4                                                           ADS92890
      IF (IPRNT3.GT.2) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS92900
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS92910
      IF (IPRNT3.GT.3) WRITE (6,320) GMAX                               ADS92920
C     ------------------------------------------------------------------ADS92930
C     CHECK FOR CONVERGENCE.                                            ADS92940
C     ------------------------------------------------------------------ADS92950
      CALL ADS211 (OBJ,OBJ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)          ADS92960
      IF (JTER.GE.JTMAX) GO TO 90                                       ADS92970
      IF (IDEL.LT.ITRMST.AND.IDAB.LT.ITRMST) GO TO 140                  ADS92980
C     OPTIMIZATION IS COMPLETE.                                         ADS92990
90    CONTINUE                                                          ADS93000
      IF (IPRNT3.EQ.0) GO TO 130                                        ADS93010
      OBJW=OBJ/SCFO                                                     ADS93020
      WRITE (6,170) JTER,OBJW                                           ADS93030
C     X-VECTOR.                                                         ADS93040
      ICALC=1                                                           ADS93050
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS93060
C     G-VECTOR.                                                         ADS93070
      ICALC=4                                                           ADS93080
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS93090
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS93100
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS93110
C     SIDE CONSTRAINTS.                                                 ADS93120
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS93130
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS93140
      WRITE (6,180) CT,CTL                                              ADS93150
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS93160
      NAVC=NAC+NVC                                                      ADS93170
      WRITE (6,230) NAC,NVC                                             ADS93180
      IF (NAVC.GT.0) WRITE (6,240)                                      ADS93190
      IF (NAVC.GT.0) WRITE (6,250) (ITMP(I),I=1,NAVC)                   ADS93200
      WRITE (6,260) NACS                                                ADS93210
      IF (NACS.EQ.0) GO TO 120                                          ADS93220
      N=0                                                               ADS93230
      DO 110 I=1,NDV                                                    ADS93240
      IF (ISIDC(I).EQ.0) GO TO 110                                      ADS93250
      IF (ISIDC(I).EQ.1) GO TO 100                                      ADS93260
C     LOWER BOUND IS ACTIVE.                                            ADS93270
      N=N+1                                                             ADS93280
      ITMP(N)=-I                                                        ADS93290
100   IF (ISIDC(I).LT.0) GO TO 110                                      ADS93300
C     UPPER BOUND IS ACTIVE.                                            ADS93310
      N=N+1                                                             ADS93320
      ITMP(N)=I                                                         ADS93330
110   CONTINUE                                                          ADS93340
      WRITE (6,270)                                                     ADS93350
      WRITE (6,250) (ITMP(I),I=1,NACS)                                  ADS93360
120   WRITE (6,190)                                                     ADS93370
      IF (JTER.GE.JTMAX) WRITE (6,200)                                  ADS93380
      IF (IDEL.GE.ITRMST) WRITE (6,210) ITRMST                          ADS93390
      IF (IDAB.GE.ITRMST) WRITE (6,220) ITRMST                          ADS93400
130   CONTINUE                                                          ADS93410
      INFO=0                                                            ADS93420
      IGOTO=0                                                           ADS93430
      RETURN                                                            ADS93440
140   CONTINUE                                                          ADS93450
C     ------------------------------------------------------------------ADS93460
C     UPDATE MOVE LIMIT.                                                ADS93470
C     ------------------------------------------------------------------ADS93480
      RLOB=100.                                                         ADS93490
      ABOB=ABS(OBJ1)                                                    ADS93500
      IF (ABOB.LE.ZRO) GO TO 150                                        ADS93510
      RLOB=ABS(OBJ/OBJ1)                                                ADS93520
      RLOB=ABS(1.-RLOB)                                                 ADS93530
150   CONTINUE                                                          ADS93540
      ABOB=ABS(OBJ-OBJ1)                                                ADS93550
      CALL ADS323 (IMVCH,RLOB,ABOB,DELSTR,DABSTR,OBJ,OBJR,ZRO,RMVLMT)   ADS93560
      OBJR=OBJ1                                                         ADS93570
      OBJ1=OBJ                                                          ADS93580
      GO TO 20                                                          ADS93590
C     ------------------------------------------------------------------ADS93600
C     FORMATS.                                                          ADS93610
C     ------------------------------------------------------------------ADS93620
160   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS93630
170   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS93640
     1'NUMBER OF CONSTRAINED MINIMIZATIONS =',I5,//,5X,'OBJECTIVE =',   ADS93650
     21PE12.5)                                                          ADS93660
180   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,'  CTL =',      ADS93670
     11PE12.5)                                                          ADS93680
190   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS93690
200   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS93700
210   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS93710
     1' CONSECUTIVE ITERATIONS')                                        ADS93720
220   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS93730
     1' CONSECUTIVE ITERATIONS')                                        ADS93740
230   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS93750
     1' VIOLATED CONSTRAINTS')                                          ADS93760
240   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS93770
250   FORMAT (5X,10I5)                                                  ADS93780
260   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS93790
270   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS93800
280   FORMAT ('1',4X,75('-'),/,22X,                                     ADS93810
     1'ISTRAT = 6: SEQUENTIAL LINEAR PROGRAMMING',/,5X,75('-'))         ADS93820
290   FORMAT (/5X,17H-- INITIAL DESIGN)                                 ADS93830
300   FORMAT (/,5X,'OBJ =',1PE12.5)                                     ADS93840
310   FORMAT (/,5X,'-- BEGIN SLP ITERATION',I5,5X,'RMVLMT=',F8.5)       ADS93850
320   FORMAT (/,5X,'GMAX=',1PE12.4)                                     ADS93860
330   FORMAT (/,5X,'RETAINED CONSTRAINTS (',I4,' )')                    ADS93870
340   FORMAT (5X,15I5)                                                  ADS93880
      END                                                               ADS93890
C --- SUBROUTINE ADS321 (IOPT,IONED,NDV,NGT,X,D,A,NRA,DF,GS,ICLP,       ADS93900
C    1G,DFL,DLB,DUB,IDGLP,NCOLA,WK,IWK,RMVLMT)                          ADS93910
C     DIMENSION D(NDV),A(NRA,NCOLA),DF(NDV),GS(NCOLA),ICLP(NCOLA),      ADS93920
C    1G(NGT),DFL(NDV),IDGLP(NCOLA),WK(NRWK),IWK(NRIWK),DLB(NDV),        ADS93930
C    2DUB(NDV),S(NDV),ATMP(NATMP,NCOLA),X(1)                            ADS93940
C     ------------------------------------------------------------------ADS93950
C --- ROUTINE TO SOLVE THE LINEAR PROGRAMMING SUB-PROBLEM.              ADS93960
C     ------------------------------------------------------------------ADS93970
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                APRIL 1983   ADS93980
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS93990
C --- INPUT                                                             ADS94000
C     IOPT  - OPTIMIZATION IDENTIFIER.                                  ADS94010
C     IONED - ONE-DIMENSIONAL SEARCH IDENTIFIER.                        ADS94020
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS94030
C     NGT   - TOTAL NUMBER OF GRADIENTS CALCULATED.                     ADS94040
C     X     - CURRENT VALUE OF PRIMAL DESIGN VARIABLE.                  ADS94050
C     D     - DESIGN VARIABLES OF SUB-PROBLEM.                          ADS94060
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS94070
C             PROBLEM.                                                  ADS94080
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS94090
C     DF    - GRADIENT OF PRIMAL OBJECTIVE FUNCTION.                    ADS94100
C     GS    - CURRENT CONSTRAINT VALUE.                                 ADS94110
C     ICLP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR L.P. SUB-PROBLEM.ADS94120
C     G     - CONSTRAINT VALUE OF SUB-PROBLEM.                          ADS94130
C     DFL   - 'DF' FOR SUB-PROBLEM                                      ADS94140
C     DLB,DUB - UPPER AND LOWER BOUNDS ON DESIGN VARIABLES OF SUB-      ADS94150
C               PROBLEM.                                                ADS94160
C     IDGLP - CONSTRAINT TYPE IDENTIFIER FOR SUB-PROBLEM.               ADS94170
C     NCOLA - NUMBER OF COLUMNS OF ARRAY A.                             ADS94180
C     WK    - REAL WORK ARRAY.                                          ADS94190
C     IWK   - INTEGER WORK ARRAY.                                       ADS94200
C     RMVLMT- MOVE LIMIT.                                               ADS94210
C --- OUTPUT                                                            ADS94220
C     D     - OPTIMUM VALUE OF DESIGN VARIABLES OF SUB-PROBLEM.         ADS94230
C                                                                       ADS94240
      SUBROUTINE ADS321 (IOPT,IONED,NDV,NGT,X,D,A,NRA,DF,GS,ICLP,G,DFL,DADS94250
     1LB,DUB,IDGLP,NCOLA,WK,IWK,RMVLMT)                                 ADS94260
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION D(1),A(NRA,1),DF(1),GS(1),ICLP(1),G(1),DFL(1),IDGLP(1),WADS94270
     1K(1),IWK(1),DLB(1),DUB(1),X(1)                                    ADS94280
C                                                                       ADS94290
      INFO=0                                                            ADS94300
      ISTRAT=6                                                          ADS94310
      DO 10 I=1,NDV                                                     ADS94320
      ABSX=ABS(X(I))                                                    ADS94330
      IF (ABSX.LE.1.0) ABSX=1.0                                         ADS94340
      XLIM=RMVLMT*ABSX                                                  ADS94350
      DD=DLB(I)-X(I)                                                    ADS94360
      DLB(I)=DMAX1(-XLIM,DD)                                            ADS94370
      DD=DUB(I)-X(I)                                                    ADS94380
      DUB(I)=DMIN1(XLIM,DD)                                             ADS94390
10    CONTINUE                                                          ADS94400
      DO 20 I=1,NGT                                                     ADS94410
      ICLP(I)=I                                                         ADS94420
      IF (IDGLP(I).GE.0) IDGLP(I)=2                                     ADS94430
      IF (IDGLP(I).LT.0) IDGLP(I)=-2                                    ADS94440
20    CONTINUE                                                          ADS94450
C     CALCULATION OF  INITIAL VALUE.                                    ADS94460
      DO 30 I=1,NDV                                                     ADS94470
      D(I)=0.0                                                          ADS94480
      IF (D(I).LT.DLB(I)) D(I)=DLB(I)                                   ADS94490
      IF (D(I).GT.DUB(I)) D(I)=DUB(I)                                   ADS94500
30    CONTINUE                                                          ADS94510
C     BEGIN OPTIMIZATION OF LP-SUBPROBLEM                               ADS94520
40    CONTINUE                                                          ADS94530
      CALL ADS501 (NDV,D,DF,FL)                                         ADS94540
      DO 60 K=1,NGT                                                     ADS94550
      J=ICLP(K)                                                         ADS94560
      G(J)=GS(J)                                                        ADS94570
      DO 50 I=1,NDV                                                     ADS94580
      G(J)=G(J)+A(I,K)*D(I)                                             ADS94590
50    CONTINUE                                                          ADS94600
60    CONTINUE                                                          ADS94610
      IF (INFO.EQ.1) GO TO 80                                           ADS94620
      DO 70 I=1,NDV                                                     ADS94630
      DFL(I)=DF(I)                                                      ADS94640
70    CONTINUE                                                          ADS94650
80    CONTINUE                                                          ADS94660
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NGT,FL,D,DLB,DUB,G,IDGLP,NADS94670
     1GT,ICLP,DFL,A,NRA,NCOLA,WK,IWK)                                   ADS94680
      IF (INFO.EQ.0) RETURN                                             ADS94690
      GO TO 40                                                          ADS94700
      END                                                               ADS94710
C --- SUBROUTINE ADS323 (IMVCH,RLOB,ABOB,DELSTR,DABSTR,OBJ,O            ADS94720
C    1BJR,ZRO,RMVLMT)                                                   ADS94730
C     ------------------------------------------------------------------ADS94740
C --- ROUTINE FOR UPDATING THE MOVE-LIMIT FOR S.L.P..                   ADS94750
C     ------------------------------------------------------------------ADS94760
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                APRIL 1983   ADS94770
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS94780
C --- INPUT                                                             ADS94790
C     IMVCH  - SATISFIED NUMBER OF CONVERGENCE CRITERION.               ADS94800
C     RLOB   - RELATIVE OBJECTIVE CHANGE.                               ADS94810
C     ABOB   - ABSOLUTE OBJECTIVE CHANGE.                               ADS94820
C     DELSTR - CONVERGENCE CRITERION FOR RELATIVE OBJECTIVE CHANGE.     ADS94830
C     DABSTR - CONVERGENCE CRITERION FOR ABSOLUTE OBJECTIVE CHANGE.     ADS94840
C     OBJ    - CURRENT VALUE OF OBJECTIVE.                              ADS94850
C     OBJR   - OBJECTIVE VALUE OF TWO ITERATIONS BEFORE.                ADS94860
C     ZRO    - VALUE OF 0.                                              ADS94870
C --- OUTPUT                                                            ADS94880
C     RMVLMT - MOVELIMIT.                                               ADS94890
C                                                                       ADS94900
      SUBROUTINE ADS323 (IMVCH,RLOB,ABOB,DELSTR,DABSTR,OBJ,OBJR,ZRO,RMVLADS94910
     1MT)                                                               ADS94920
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C                                                                       ADS94930
      DELST=5.*DELSTR                                                   ADS94940
      DABST=5.*DABSTR                                                   ADS94950
      RLOSC=100.                                                        ADS94960
      ABOSC=ABS(OBJR)                                                   ADS94970
      IF (ABOSC.LE.ZRO) GO TO 10                                        ADS94980
      RLOSC=ABS(OBJ/OBJR)                                               ADS94990
      RLOSC=ABS(1.-RLOSC)                                               ADS95000
10    CONTINUE                                                          ADS95010
      ABOSC=ABS(OBJ-OBJR)                                               ADS95020
C                                                                       ADS95030
      JMV=IMVCH+1                                                       ADS95040
      GO TO (20,50,40), JMV                                             ADS95050
20    CONTINUE                                                          ADS95060
      IF (RLOB.LT.DELST.OR.ABOB.LT.DABSTR) GO TO 30                     ADS95070
      IF (RLOSC.LT.DELST.OR.ABOSC.LT.DABST) GO TO 30                    ADS95080
      RETURN                                                            ADS95090
30    CONTINUE                                                          ADS95100
      IMVCH=1                                                           ADS95110
40    CONTINUE                                                          ADS95120
      RMVLMT=0.5*RMVLMT                                                 ADS95130
      IF (RMVLMT.LT.1.0E-4) RMVLMT=1.0E-4                               ADS95140
      RETURN                                                            ADS95150
50    CONTINUE                                                          ADS95160
      IF (RLOB.LT.DELSTR.OR.ABOB.LT.DABSTR) GO TO 60                    ADS95170
      IF (RLOSC.LT.DELST.OR.ABOSC.LT.DABST) GO TO 60                    ADS95180
      RETURN                                                            ADS95190
60    CONTINUE                                                          ADS95200
      IMVCH=2                                                           ADS95210
      GO TO 40                                                          ADS95220
      END                                                               ADS95230
C --- SUBROUTINE ADS325                                                 ADS95240
C     ------------------------------------------------------------------ADS95250
C --- ROUTINE FOR PERFORM SEQUENTIAL CONVEX PROGRAMMING.                ADS95260
C     ------------------------------------------------------------------ADS95270
C     BY G. N. VANDERPLAATS                           SEPT. 1984        ADS95280
C     ENGINEERING DESIGN OPTIMIZATION, INC., SANTA BARBARA, CALIF.      ADS95290
C --- INPUT                                                             ADS95300
C     INFO  - CONTROL PARAMETER.                                        ADS95310
C     IOPT  - OPTIMIZER IDENTIFIER.                                     ADS95320
C     IONED - ONE DIMENSIONAL SEARCH IDENTIFIER.                        ADS95330
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS95340
C     NCON  - NUMBER OF CONSTRAINTS.                                    ADS95350
C     OBJ   - VALUE OF PRIMAL OBJECTIVE FUNCTION.                       ADS95360
C     G     - CURRENT VALUES OF CONSTRAINTS.                            ADS95370
C     NGT   - TOTAL NUMBER OF GRADIENTS CALCULATED.                     ADS95380
C     IC    - ARRAY IDENTIFYING CONSTRAINTS THAT ARE ACTIVE.            ADS95390
C     DF    - GRADIENT OF PRIMAL OBJECTIVE FUNCTION.                    ADS95400
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS95410
C             PROBLEM.                                                  ADS95420
C     NCOLA - DIMENSIONED COLUMNS OF A.                                 ADS95430
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS95440
C     X     - INITIAL VALLUE FOR DESIGN VARIABLES.                      ADS95450
C     IDG   - CONSTRAINT TYPE IDENTIFIER.                               ADS95460
C     WK    - REAL WORK ARRAY.                                          ADS95470
C     IWK   - INTEGER WORK ARRAY.                                       ADS95480
C     ICLP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR CONVEX           ADS95490
C             SUB-PROBLEM.                                              ADS95500
C     SVLB,SVUB - UPPER AND LOWER LIMITS FOR PRIMAL DESIGN VARIABLES.   ADS95510
C     S     - WORK ARRAY (NDV).                                         ADS95520
C     VUBST,VLBST - ARRAY TO TEMPORARILY STORE VLB,VUB.                 ADS95530
C     GS    - ARRAY TO STORE THE CURRENT VALUE OF CONSTRAINTS.          ADS95540
C     ZRO   - VALUE OF ZERO.                                            ADS95550
C     IPRNT3 - PRINT CONTROL FOR STRATEGY LEVEL.                        ADS95560
C     SCFX  - X-VECTOR SCALING FACTORS.                                 ADS95570
C     XTMP  - ARRAY FOR WORK.                                           ADS95580
C     SCLMIN - MINIMUM VALUE OF SCALING FACTOR.                         ADS95590
C     DELSTR - CONVERGENCE CRITERIOR OF RELATIVE CHANGE.                ADS95600
C     DABSTR - CONVERGENCE CRITERIOR OF ABSOLUTE CHANGE.                ADS95610
C     ITRMST - MINIMUM NUMBER OF CONVERGENCE CRITERIOR SATISFACTION.    ADS95620
C     JTMAX  - MAXIMUM NUMBER OF ITERATION.                             ADS95630
C     SCFO   - SCALE FACTOR OF OBJECTIVE.                               ADS95640
C     XSUB   - DESIGN VARIABLES OF SUB-PROBLEM.                         ADS95650
C     DFSUB  - GRADIENTS OF OBJECTIVE OF SUB-PROBLEM.                   ADS95660
C --- OUTPUT                                                            ADS95670
C     X     - OPTIMUM VALUE FOR DESIGN VARIABLES.                       ADS95680
C                                                                       ADS95690
      SUBROUTINE ADS325 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,ADS95700
     1A,NCOLA,NRA,X,IDG,WK,IWK,ICLP,SVLB,SVUB,S,GS,VLBST,VUBST,ZRO,IPRNTADS95710
     23,SCFX,XTMP,SCLMIN,DELSTR,DABSTR,ITRMST,JTMAX,SCFO,ITMP,CT,       ADS95720
     3CTMIN,CTL,CTLMIN,NAC,NVC,XSUB,DFSUB,OBJR,OBJ1,IDAB,IDEL,          ADS95730
     4JTER,ISIDC,B,NRB,C,THETA,GGL,AA,NRAA,RMVLMZ,RMVLMT,SCFG)          ADS95740
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IC(1),DF(1),A(NRA,1),X(1),IDG(1),WK(1),IWK(1),ICLP(ADS95750
     11),SVLB(1),SVUB(1),S(1),GS(1),VLBST(1),VUBST(1),SCFX(1),XTMP(1),ITADS95760
     2MP(1),XSUB(1),DFSUB(1),ISIDC(1),B(NRB,1),C(1),THETA(1),GGL(1),    ADS95770
     3AA(NRAA,1),SCFG(1)                                                ADS95780
C                                                                       ADS95790
      IGOTO1=IGOTO+1                                                    ADS95800
      GO TO (10,30,80), IGOTO1                                          ADS95810
10    CONTINUE                                                          ADS95820
      RMVLMT=RMVLMZ                                                     ADS95830
      IF (IPRNT3.GT.0) WRITE (6,280)                                    ADS95840
      DELST=5.*DELSTR                                                   ADS95850
      DABST=5.*DABSTR                                                   ADS95860
      OBJ1=OBJ                                                          ADS95870
      JTER=0                                                            ADS95880
C     INITIAL FUNCTION EVALUATION IS AVAILABLE.                         ADS95890
      IF (IPRNT3.GT.0) WRITE (6,290)                                    ADS95900
      OBJW=OBJ/SCFO                                                     ADS95910
      IF (IPRNT3.GT.0) WRITE (6,300) OBJW                               ADS95920
      ICALC=1                                                           ADS95930
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS95940
      ICALC=4                                                           ADS95950
      IF (IPRNT3.GT.0) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS95960
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS95970
C     ------------------------------------------------------------------ADS95980
C     ITERATION START                                                   ADS95990
C     ------------------------------------------------------------------ADS96000
20    CONTINUE                                                          ADS96010
      JTER=JTER+1                                                       ADS96020
      IF (IPRNT3.GT.0) WRITE (6,310) JTER                               ADS96030
      IF (JTER.EQ.1) GO TO 30                                           ADS96040
C     DETERMINE THE ACTIVE CONSTRAINTS SET.                             ADS96050
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS96060
C     CALCULATE GRADIENTS OF OBJECTIVE FUNCTION AND CONSTRAINTS.        ADS96070
      IGOTO=1                                                           ADS96080
      INFO=2                                                            ADS96090
      RETURN                                                            ADS96100
C                                                                       ADS96110
30    CONTINUE                                                          ADS96120
      IF(JTER.EQ.1) GO TO 35                                            ADS96130
      IF(GMAX.GT.GMAX1) RMVLMT=.5*RMVLMT                                ADS96140
      IF(RMVLMT.LT.0.01) RMVLMT=0.01                                    ADS96150
35    GMAX1=GMAX                                                        ADS96160
      IF (IPRNT3.GT.3) WRITE (6,330) NGT                                ADS96170
      IF (IPRNT3.GT.3) WRITE (6,340) (IC(I),I=1,NGT)                    ADS96180
      ICALC=3                                                           ADS96190
      IF (IPRNT3.LT.5) GO TO 60                                         ADS96200
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS96210
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS96220
      ISCAL=1                                                           ADS96230
      DO 50 J=1,NGT                                                     ADS96240
      WRITE (6,160) IC(J)                                               ADS96250
      CALL ADS401 (A(1,J),NDV,SCFG,ISCAL)                               ADS96260
50    CONTINUE                                                          ADS96270
60    CONTINUE                                                          ADS96280
C     ----------------------------------------------------------------- ADS96290
C     CHECK CONVERGENCE BY RESIDUAL OF KUHN TUCKER CONDITIONS.          ADS96300
C     ----------------------------------------------------------------- ADS96310
      RSDMX=1.                                                          ADS96320
      CTTMP=-CTMIN                                                      ADS96330
      CTLTMP=-CTLMIN                                                    ADS96340
C     DETERMINE ACTIVE/VIOLATED CONSTRAINTS SET.                        ADS96350
      CALL ADS009 (NDV,A,NRA,NGT,IC,IDG,G,CTTMP,CTMIN,CTLTMP,CTLMIN,NACTADS96360
     1,NVCT)                                                            ADS96370
      NAVCT=NACT+NVCT                                                   ADS96380
      IF (NVCT.GT.0) GO TO 70                                           ADS96390
C     IDENTIFY ACTIVE SIDE CONSTRAINTS.                                 ADS96400
      CALL ADS509 (NDV,X,SVLB,SVUB,ZRO,ISIDC,NACS)                      ADS96410
      IF(NAVCT.EQ.0) GO TO 61                                           ADS96420
      CALL ADS207 (NDV,NAVCT,A,NRA,ZRO,GGL)                             ADS96430
61    N2=1                                                              ADS96440
      CALL ADS207 (NDV,N2,DF,NDV,ZRO,SCALF)                             ADS96450
      IF(NAVCT.EQ.0.AND.NACS.EQ.0) GO TO 70                             ADS96460
C     CALCULATE THE RESIDUAL OF KUHN-TUCKER CONDITION.                  ADS96470
      THZ=0.                                                            ADS96480
      DUM=1.                                                            ADS96490
      NDUM=0                                                            ADS96500
      CALL ADS217 (NDV,NAVCT,NDUM,NACS,DUM,DUM,DUM,THZ,ZRO,DUM,GGL,G,DF,ADS96510
     1A,NRA,IDG,IC,ISIDC,ITMP,B,NRB,SCALF,S,THETA,C,BETA,RSDMX)         ADS96520
      IF (IPRNT3.GE.4) WRITE (6,350) RSDMX                              ADS96530
      IF (RSDMX.LT.0.001) GO TO 90                                      ADS96540
70    CONTINUE                                                          ADS96550
C     ------------------------------------------------------------------ADS96560
C     SOLVE CONVEX SUB-PROBLEM.                                         ADS96570
C     ------------------------------------------------------------------ADS96580
C     LOCAL FUNCTIONS.                                                  ADS96590
      DO 40 J=1,NGT                                                     ADS96600
      ICJ=IC(J)                                                         ADS96610
      GS(J)=G(ICJ)                                                      ADS96620
      ITMP(J)=IDG(ICJ)                                                  ADS96630
40    CONTINUE                                                          ADS96640
C     MOVE LIMITS                                                       ADS96650
      DO 71 I=1,NDV                                                     ADS96660
      XX=RMVLMT*ABS(X(I))                                               ADS96670
      XX=DMAX1(XX,RMVLMT)                                               ADS96680
      XI=X(I)-XX                                                        ADS96690
      VLBST(I)=DMAX1(XI,SVLB(I))                                        ADS96700
      XI=X(I)+XX                                                        ADS96710
      VUBST(I)=DMIN1(XI,SVUB(I))                                        ADS96720
71    CONTINUE                                                          ADS96730
      CALL ADS326 (IOPT,IONED,NDV,NGT,X,XSUB,A,NRA,DF,GS,ICLP,G,DFSUB,  ADS96740
     1ITMP,NCOLA,WK,IWK,OBJ,VLBST,VUBST,AA,NRAA,ZRO)                    ADS96750
C     CALCULATE THE VALUES OF OBJECTIVE FUNCTION AND CONSTRAINTS FOR NEWADS96760
C     DESIGN VARIABLES.                                                 ADS96770
      IGOTO=2                                                           ADS96780
      INFO=1                                                            ADS96790
      RETURN                                                            ADS96800
80    CONTINUE                                                          ADS96810
      OBJW=OBJ/SCFO                                                     ADS96820
      IF (IPRNT3.GT.0) WRITE (6,300) OBJW                               ADS96830
      ICALC=1                                                           ADS96840
      IF (IPRNT3.GT.1) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS96850
      ICALC=4                                                           ADS96860
      IF (IPRNT3.GT.2) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS96870
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS96880
      IF (IPRNT3.GT.3) WRITE (6,320) GMAX                               ADS96890
C     ------------------------------------------------------------------ADS96900
C     CHECK CONVERGENCE.                                                ADS96910
C     ------------------------------------------------------------------ADS96920
      CALL ADS211 (OBJ,OBJ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)          ADS96930
      OBJ1=OBJ                                                          ADS96940
      IF (JTER.GE.JTMAX) GO TO 90                                       ADS96950
      IF (IDEL.GT.0.OR.IDAB.GT.0) RMVLMT=.5*RMVLMT                      ADS96960
      IF (RMVLMT.LT.0.01) RMVLMT=.01                                    ADS96970
      IF (IDEL.LT.ITRMST.AND.IDAB.LT.ITRMST) GO TO 20                   ADS96980
C     ----------------------------------------------------------------- ADS96990
C     OPTIMIZATION IS COMPLETE.                                         ADS97000
C     ----------------------------------------------------------------- ADS97010
90    CONTINUE                                                          ADS97020
      IF (IPRNT3.EQ.0) GO TO 130                                        ADS97030
      OBJW=OBJ/SCFO                                                     ADS97040
      WRITE (6,170) JTER,OBJW                                           ADS97050
C     X-VECTOR.                                                         ADS97060
      ICALC=1                                                           ADS97070
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS97080
C     G-VECTOR.                                                         ADS97090
      ICALC=4                                                           ADS97100
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS97110
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS97120
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS97130
C     SIDE CONSTRAINTS.                                                 ADS97140
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS97150
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS97160
      WRITE (6,180) CT,CTL                                              ADS97170
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS97180
      NAVC=NAC+NVC                                                      ADS97190
      WRITE (6,230) NAC,NVC                                             ADS97200
      IF (NAVC.GT.0) WRITE (6,240)                                      ADS97210
      IF (NAVC.GT.0) WRITE (6,250) (ITMP(I),I=1,NAVC)                   ADS97220
      WRITE (6,260) NACS                                                ADS97230
      IF (NACS.EQ.0) GO TO 120                                          ADS97240
      N=0                                                               ADS97250
      DO 110 I=1,NDV                                                    ADS97260
      IF (ISIDC(I).EQ.0) GO TO 110                                      ADS97270
      IF (ISIDC(I).EQ.1) GO TO 100                                      ADS97280
C     LOWER BOUND IS ACTIVE.                                            ADS97290
      N=N+1                                                             ADS97300
      ITMP(N)=-I                                                        ADS97310
100   IF (ISIDC(I).LT.0) GO TO 110                                      ADS97320
C     UPPER BOUND IS ACTIVE.                                            ADS97330
      N=N+1                                                             ADS97340
      ITMP(N)=I                                                         ADS97350
110   CONTINUE                                                          ADS97360
      WRITE (6,270)                                                     ADS97370
      WRITE (6,250) (ITMP(I),I=1,NACS)                                  ADS97380
120   WRITE (6,190)                                                     ADS97390
      IF (JTER.GE.JTMAX) WRITE (6,200)                                  ADS97400
      IF (IDEL.GE.ITRMST) WRITE (6,210) ITRMST                          ADS97410
      IF (IDAB.GE.ITRMST) WRITE (6,220) ITRMST                          ADS97420
      IF(RSDMX.LT.0.001) WRITE(6,360) RSDMX                             ADS97430
130   CONTINUE                                                          ADS97440
      INFO=0                                                            ADS97450
      IGOTO=0                                                           ADS97460
      RETURN                                                            ADS97470
C     ------------------------------------------------------------------ADS97480
C     FORMATS.                                                          ADS97490
C     ------------------------------------------------------------------ADS97500
160   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS97510
170   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS97520
     1'NUMBER OF CONSTRAINED MINIMIZATIONS =',I5,//,5X,'OBJECTIVE =',   ADS97530
     21PE12.5)                                                          ADS97540
180   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,'  CTL =',      ADS97550
     11PE12.5)                                                          ADS97560
190   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS97570
200   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS97580
210   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS97590
     1' CONSECUTIVE ITERATIONS')                                        ADS97600
220   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS97610
     1' CONSECUTIVE ITERATIONS')                                        ADS97620
230   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS97630
     1' VIOLATED CONSTRAINTS')                                          ADS97640
240   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS97650
250   FORMAT (5X,10I5)                                                  ADS97660
260   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS97670
270   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS97680
280   FORMAT ('1',/,5X,75('-'),/,22X,                                   ADS97690
     1'ISTRAT = 9: SEQUENTIAL CONVEX PROGRAMMING',/,5X,75('-'))         ADS97700
290   FORMAT (/,5X,'-- INITIAL DESIGN')                                 ADS97710
300   FORMAT (/,5X,'OBJ =',1PE12.5)                                     ADS97720
310   FORMAT (/,5X,'-- BEGIN SCP ITERATION',I5)                         ADS97730
320   FORMAT (/,5X,'GMAX=',1PE12.4)                                     ADS97740
330   FORMAT (/,5X,'RETAINED CONSTRAINTS (',I4,' )')                    ADS97750
340   FORMAT (5X,15I5)                                                  ADS97760
350   FORMAT(/,5X,'MAXIMUM RESIDUAL OF KUHN-TUCKER CONDITIONS =',       ADS97770
     11PE12.5)                                                          ADS97780
360   FORMAT(/,5X,'MAXIMUM K-D RESIDUAL =',1PE12.5,                     ADS97790
     1' IS LESS THAN 0.001')                                            ADS97800
      END                                                               ADS97810
C --- SUBROUTINE ADS326                                                 ADS97820
C     ------------------------------------------------------------------ADS97830
C --- ROUTINE TO SOLVE THE CONVEX PROGRAMMING SUB-PROBLEM.              ADS97840
C     ------------------------------------------------------------------ADS97850
C     BY G. N. VANDERPLAATS                           SEPT. 1984        ADS97860
C     ENGINEERING DESIGN OPTIMIZATION, INC., SANTA BARBARA, CALIF.      ADS97870
C --- INPUT                                                             ADS97880
C     IOPT  - OPTIMIZATION IDENTIFIER.                                  ADS97890
C     IONED - ONE-DIMENSIONAL SEARCH IDENTIFIER.                        ADS97900
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS97910
C     NGT   - TOTAL NUMBER OF GRADIENTS CALCULATED.                     ADS97920
C     X     - CURRENT VALUE OF PRIMAL DESIGN VARIABLE.                  ADS97930
C     D     - DESIGN VARIABLES OF SUB-PROBLEM.                          ADS97940
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS97950
C             PROBLEM.                                                  ADS97960
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS97970
C     DF    - GRADIENT OF PRIMAL OBJECTIVE FUNCTION.                    ADS97980
C     GS    - CURRENT CONSTRAINT VALUE.                                 ADS97990
C     ICLP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR L.P. SUB-PROBLEM.ADS98000
C     G     - CONSTRAINT VALUE OF SUB-PROBLEM.                          ADS98010
C     DFL   - 'DF' FOR SUB-PROBLEM                                      ADS98020
C     DLB,DUB - UPPER AND LOWER BOUNDS ON DESIGN VARIABLES OF SUB-      ADS98030
C             PROBLEM.                                                  ADS98040
C     IDGLP - CONSTRAINT TYPE IDENTIFIER FOR SUB-PROBLEM.               ADS98050
C     NCOLA - NUMBER OF COLUMNS OF ARRAY A.                             ADS98060
C     WK    - REAL WORK ARRAY.                                          ADS98070
C     IWK   - INTEGER WORK ARRAY.                                       ADS98080
C --- OUTPUT                                                            ADS98090
C     D     - OPTIMUM VALUE OF DESIGN VARIABLES OF SUB-PROBLEM.         ADS98100
C                                                                       ADS98110
      SUBROUTINE ADS326 (IOPT,IONED,NDV,NGT,X,D,A,NRA,DF,GS,ICLP,G,DFL, ADS98120
     1IDGLP,NCOLA,WK,IWK,OBJ,VLBST,VUBST,AA,NRAA,ZRO)                   ADS98130
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION D(1),A(NRA,1),DF(1),GS(1),ICLP(1),G(1),DFL(1),IDGLP(1),WADS98140
     1K(1),IWK(1),X(1),VLBST(1),VUBST(1),AA(NRAA,1)                     ADS98150
      NCOLAA=NGT                                                        ADS98160
C                                                                       ADS98170
      INFO=0                                                            ADS98180
      ISTRAT=9                                                          ADS98190
      DO 20 I=1,NGT                                                     ADS98200
      ICLP(I)=I                                                         ADS98210
20    CONTINUE                                                          ADS98220
C     CONSTANT TERMS.                                                   ADS98230
      OBJZ=OBJ                                                          ADS98240
      DO 22 I=1,NDV                                                     ADS98250
      D(I)=X(I)                                                         ADS98260
22    CONTINUE                                                          ADS98270
C     BEGIN OPTIMIZATION OF CONVEX-SUBPROBLEM                           ADS98280
40    CONTINUE                                                          ADS98290
C     EVALUATE OBJECTIVE AND CONSTRAINTS.                               ADS98300
      FL=OBJZ                                                           ADS98310
      DO 42 J=1,NGT                                                     ADS98320
42    G(J)=GS(J)                                                        ADS98330
      DO 50 I=1,NDV                                                     ADS98340
      XI=D(I)                                                           ADS98350
      VLI=VLBST(I)                                                      ADS98360
      XIZ=X(I)                                                          ADS98370
      DFI=DF(I)                                                         ADS98380
      IF(VLI.GT.ZRO.AND.DFI.LT.0.) GO TO 43                             ADS98390
      FL=FL+DFI*(XI-XIZ)                                                ADS98400
      GO TO 44                                                          ADS98410
43    FL=FL+DFI*XIZ*(1.-XIZ/XI)                                         ADS98420
44    CONTINUE                                                          ADS98430
      DO 50 J=1,NGT                                                     ADS98440
      AJ=A(I,J)                                                         ADS98450
      K=IDGLP(J)                                                        ADS98460
      IF (K.LT.0.OR.K.GT.1) GO TO 51                                    ADS98470
      IF(VLI.GT.ZRO.AND.AJ.LT.0.) GO TO 52                              ADS98480
51    G(J)=G(J)+AJ*(XI-XIZ)                                             ADS98490
      GO TO 54                                                          ADS98500
52    G(J)=G(J)+AJ*XIZ*(1.-XIZ/XI)                                      ADS98510
54    CONTINUE                                                          ADS98520
50    CONTINUE                                                          ADS98530
      IF(INFO.EQ.1) GO TO 80                                            ADS98540
C     GRADIENTS.                                                        ADS98550
      DO 70 I=1,NDV                                                     ADS98560
      XI=D(I)                                                           ADS98570
      VLI=VLBST(I)                                                      ADS98580
      XIZ=X(I)                                                          ADS98590
      DFL(I)=DF(I)                                                      ADS98600
      IF(VLI.GT.ZRO.AND.DF(I).LT.0.) DFL(I)=DFL(I)*XIZ*XIZ/(XI*XI)      ADS98610
      DO 75 J=1,NGT                                                     ADS98620
      JJ=ICLP(J)                                                        ADS98630
      AA(I,J)=A(I,JJ)                                                   ADS98640
      K=IDGLP(J)                                                        ADS98650
      IF(K.LT.0.OR.K.GT.1) GO TO 75                                     ADS98660
      IF(VLI.GT.ZRO.AND.AA(I,J).LT.0.) AA(I,J)=AA(I,J)*XIZ*XIZ/(XI*XI)  ADS98670
75    CONTINUE                                                          ADS98680
70    CONTINUE                                                          ADS98690
80    CONTINUE                                                          ADS98700
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NGT,FL,D,VLBST,VUBST,G,   ADS98710
     1IDGLP,NGT,ICLP,DFL,AA,NRAA,NCOLAA,WK,IWK)                         ADS98720
      IF (INFO.GT.0) GO TO 40                                           ADS98730
      DO 85 I=1,NDV                                                     ADS98740
85    X(I)=D(I)                                                         ADS98750
      RETURN                                                            ADS98760
      END                                                               ADS98770
C --- SUBROUTINE ADS330 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF ADS98780
C    1 ,A,NCOLA,NRA,X,IDG,WK,IWK,ICCT,SVLB,SVUB,S,GS,SQRTG,VLBST,VUBST, ADS98790
C    2 ZRO,IPRNT3,RMVLMZ,SCFX,XTMP,SCLMIN,DELSTR,DABSTR,ITRMST,         ADS98800
C    3 JTMAX,BETA,SCFO,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC,XOLD,XSUB,DFSUB)ADS98810
C     DIMENSION G(NCON),IC(NCOLA),DF(NDV+1),A(NRA,NCOLA),X(NDV+1),      ADS98820
C    1IDG(NCON OR NGT+1),WK(NRWK),IWK(NRIWK),ICCT(NGT+1),SVLB(NDV+1),   ADS98830
C    2SVUB(NDV+1),S(NDV),GS(NGT+1),SQRTG(NGT+1),VLBST(NDV),VUBST(NDV),  ADS98840
C    3SCFX(NDV),XTMP(NDV),ITMP(NGT+1),XOLD(NDV),XSUB(NDV+1),DFSUB(NDV+1)ADS98850
C     ------------------------------------------------------------------ADS98860
C --- ROUTINE TO PERFORM THE METHOD OF CENTERS.                         ADS98870
C     ------------------------------------------------------------------ADS98880
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                APRIL 1983   ADS98890
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS98900
C --- INPUT                                                             ADS98910
C     INFO  - CONTROL PARAMETER.                                        ADS98920
C     IOPT  - OPTIMIZER IDENTIFIER.                                     ADS98930
C     IONED - ONE DIMENSIONAL SEARCH IDENTIFIER.                        ADS98940
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS98950
C     NCON  - NUMBER OF CONSTRAINTS.                                    ADS98960
C     OBJ   - VALUE OF PRIMAL OBJECTIVE FUNCTION.                       ADS98970
C     G     - CURRENT VALUES OF CONSTRAINTS.                            ADS98980
C     NGT   - TOTAL NUMBER OF GRADIENTS CALCULATED.                     ADS98990
C     IC    - ARRAY IDENTIFYING CONSTRAINTS THAT ARE ACTIVE.            ADS99000
C     DF    - GRADIENT OF PRIMAL OBJECTIVE FUNCTION.                    ADS99010
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS99020
C             PROBLEM.                                                  ADS99030
C     NCOLA - MAXIMUM NUMBER OF ACTIVE CONSTRAINTS.                     ADS99040
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS99050
C     X     - INITIAL VALLUE FOR DESIGN VARIABLES.                      ADS99060
C     IDG   - CONSTRAINT TYPE IDENTIFIER.                               ADS99070
C     WK    - REAL WORK ARRAY.                                          ADS99080
C     IWK   - INTEGER WORK ARRAY.                                       ADS99090
C     ICLP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR L.P. SUB-PROBLEM.ADS99100
C     SVLB,SVUB - UPPER AND LOWER LIMITS FOR PRIMAL DESIGN VARIABLES.   ADS99110
C     S     - WORK ARRAY (NDV).                                         ADS99120
C     VLBST,VUBST - ARRAY TO STORE VLB,VUB.                             ADS99130
C     GS    - ARRAY TO STORE THE CURRENT VALUE OF CONSTRAINTS.          ADS99140
C     SQRTG - ARRAY TO STORE THE VALUE OF SQUARE ROOT MEAN OF CONST-    ADS99150
C             RAINTS AND OBJECTIVE.                                     ADS99160
C     ZRO   - VALUE OF ZERO.                                            ADS99170
C     IPRNT3 - PRINT CONTROL FOR STRATEGY LEVEL.                        ADS99180
C     RMVLMZ - INITIAL VALUE OF MOVE LIMIT.                             ADS99190
C     SCFX  - SCALING FACTOR. ARRAY FOR QP SUB-PROBLEM.                 ADS99200
C     XTMP  - ARRAY FOR WORK.                                           ADS99210
C     SCLMIN - MINIMUM VALUE OF SCALING FACTOR.                         ADS99220
C     DELSTR - CONVERGENCE CRITERIOR OF RELATIVE CHANGE.                ADS99230
C              0.0010 IS RECOMMENDED.                                   ADS99240
C     DABSTR - CONVERGENCE CRITERIOR OF ABSOLUTE CHANGE.                ADS99250
C     ITRMST - MINIMUM NUMBER OF CONVERGENCE CRITERIOR SATISFACTION.    ADS99260
C     JTMAX  - MAXIMUM NUMBER OF ITERATION.                             ADS99270
C     BETA   - BETA. DEFAULT VALUE IS 0.                                ADS99280
C     SCFO   - SCALE FACTOR OF OBJECTIVE.                               ADS99290
C     XOLD  - ARRAY TO STORE OLD DESIGN VARIAVLES VALUE.                ADS99300
C     XSUB  - DESIGN VARIABLES OF SUB-PROBLEM                           ADS99310
C     DFSUB - GRADIENTS OF OBJECTIVE OF SUB-PROBLEM.                    ADS99320
C --- OUTPUT                                                            ADS99330
C     X     - OPTIMUM VALUE FOR DESIGN VARIABLES.                       ADS99340
C                                                                       ADS99350
      SUBROUTINE ADS330 (IGOTO,INFO,IOPT,IONED,NDV,NCON,OBJ,G,NGT,IC,DF,ADS99360
     1A,NCOLA,NRA,X,IDG,WK,IWK,ICCT,SVLB,SVUB,S,GS,SQRTG,VLBST,VUBST,ZROADS99370
     2,IPRNT3,RMVLMZ,SCFX,XTMP,SCLMIN,DELSTR,DABSTR,ITRMST,JTMAX,BETA,SCADS99380
     3FO,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC,XOLD,XSUB,DFSUB,GMAX,OBJR,OBJ1ADS99390
     4,IDAB,IDEL,IGMXC,IMVCH,JTER,ISIDC,SCFG)                           ADS99400
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IC(1),DF(1),A(NRA,1),X(1),IDG(1),WK(1),IWK(1),ICCT(ADS99410
     11),SVLB(1),SVUB(1),S(1),GS(1),SQRTG(1),VLBST(1),VUBST(1),SCFX(1),XADS99420
     2TMP(1),ITMP(1),XOLD(1),XSUB(1),DFSUB(1),ISIDC(1),SCFG(1)          ADS99430
C                                                                       ADS99440
      NDVR=NDV+1                                                        ADS99450
      IGOTO1=IGOTO+1                                                    ADS99460
      GO TO (10,40,150,180,210), IGOTO1                                 ADS99470
10    CONTINUE                                                          ADS99480
      IF (IPRNT3.GT.0) WRITE (6,450)                                    ADS99490
      OBJR=1.0E20                                                       ADS99500
      OBJ1=OBJ                                                          ADS99510
      RMVLMT=RMVLMZ                                                     ADS99520
      X(NDVR)=0.                                                        ADS99530
      SCFX(NDVR)=1.                                                     ADS99540
      JTER=0                                                            ADS99550
      IMVCH=0                                                           ADS99560
C     INITIAL FUNCTION EVALUATION IS ALREADY AVAILABLE.                 ADS99570
      IF (IPRNT3.EQ.0) GO TO 20                                         ADS99580
      WRITE (6,460)                                                     ADS99590
      OBJW=OBJ/SCFO                                                     ADS99600
      WRITE (6,470) OBJW                                                ADS99610
      ICALC=1                                                           ADS99620
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS99630
      ICALC=4                                                           ADS99640
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS99650
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS99660
20    CONTINUE                                                          ADS99670
C     IF INITIALLY VIOLATED , IGMXC=1.                                  ADS99680
      IGMXC=0                                                           ADS99690
      IF (GMAX.GT.0.) IGMXC=1                                           ADS99700
30    CONTINUE                                                          ADS99710
C     ------------------------------------------------------------------ADS99720
C     ITERATION START                                                   ADS99730
C     ------------------------------------------------------------------ADS99740
      JTER=JTER+1                                                       ADS99750
C     STORE THE FORMER VALUES.                                          ADS99760
      CALL ADS503 (NDV,X,XOLD)                                          ADS99770
      IF (IPRNT3.GT.0) WRITE (6,500) JTER,RMVLMT                        ADS99780
      IF (JTER.EQ.1) GO TO 40                                           ADS99790
C     DETERMINE THE ACTIVE CONSTRAINTS SET.                             ADS99800
      CALL ADS201 (NCON,NGT,G,IDG,IC,CT,CTMIN,CTL,CTLMIN,NAC,NVC)       ADS99810
C     CALCULATE GRADIENTS OF OBJECTIVE FUNCTION AND CONSTRAINTS.        ADS99820
      IGOTO=1                                                           ADS99830
      INFO=2                                                            ADS99840
      RETURN                                                            ADS99850
C                                                                       ADS99860
40    CONTINUE                                                          ADS99870
      NGT1=NGT                                                          ADS99880
      IF (GMAX.LE.0) NGT1=NGT+1                                         ADS99890
      IF (IPRNT3.GT.3) WRITE (6,530) NGT                                ADS99900
      IF (IPRNT3.GT.3) WRITE (6,480) (IC(I),I=1,NGT)                    ADS99910
      DO 50 J=1,NGT                                                     ADS99920
      ICJ=IC(J)                                                         ADS99930
      GS(J)=G(ICJ)                                                      ADS99940
      ITMP(J)=IDG(ICJ)                                                  ADS99950
50    CONTINUE                                                          ADS99960
      ICALC=3                                                           ADS99970
      IF (IPRNT3.LT.5) GO TO 70                                         ADS99980
      CALL ADS402 (ICALC,DF,NDV,SCFX)                                   ADS99990
C     WRITE GRADIENT OF ACTIVE/VIOLATED CONSTRAINTS.                    ADS00000
      ISCAL=0                                                           ADS00010
      DO 60 J=1,NGT                                                     ADS00020
      WRITE (6,320) IC(J)                                               ADS00030
      CALL ADS401 (A(1,J),NDV,SCFX,ISCAL)                               ADS00040
60    CONTINUE                                                          ADS00050
70    CONTINUE                                                          ADS00060
C     CALCULATE THE SQUARE ROOT MEAN OF GRADIENTS OF CONSTRAINTS AND    ADS00070
C     OBJECTIVE. AND NORMALIZE THE GRADIENTS.                           ADS00080
      DO 100 J=1,NGT                                                    ADS00090
      SQRTG(J)=0.                                                       ADS00100
      DO 80 K=1,NDV                                                     ADS00110
      SQRTG(J)=SQRTG(J)+A(K,J)*A(K,J)                                   ADS00120
80    CONTINUE                                                          ADS00130
      SQRTG(J)=SQRT(SQRTG(J))                                           ADS00140
      DO 90 K=1,NDV                                                     ADS00150
      A(K,J)=A(K,J)/SQRTG(J)                                            ADS00160
90    CONTINUE                                                          ADS00170
      GS(J)=GS(J)/SQRTG(J)                                              ADS00180
100   CONTINUE                                                          ADS00190
      NGTT=NGT+1                                                        ADS00200
      SQRTG(NGTT)=0.                                                    ADS00210
      DO 110 K=1,NDV                                                    ADS00220
      SQRTG(NGTT)=SQRTG(NGTT)+DF(K)*DF(K)                               ADS00230
110   CONTINUE                                                          ADS00240
      SQRTG(NGTT)=SQRT(SQRTG(NGTT))                                     ADS00250
      DO 120 K=1,NDV                                                    ADS00260
      DF(K)=DF(K)/SQRTG(NGTT)                                           ADS00270
      A(K,NGTT)=DF(K)                                                   ADS00280
120   CONTINUE                                                          ADS00290
      GS(NGTT)=0.                                                       ADS00300
      DO 130 J=1,NGT1                                                   ADS00310
      A(NDVR,J)=1.                                                      ADS00320
130   CONTINUE                                                          ADS00330
      CALL ADS503 (NDV,SVLB,VLBST)                                      ADS00340
      CALL ADS503 (NDV,SVUB,VUBST)                                      ADS00350
C     ------------------------------------------------------------------ADS00360
C     SOLVE L.P. SUB-PROBLEM.                                           ADS00370
C     ------------------------------------------------------------------ADS00380
      CALL ADS331 (IOPT,IONED,NDV,NDVR,NGT,NGT1,X,XSUB,A,NRA,GS,ICCT,G,DADS00390
     1F,SVLB,SVUB,ITMP,NCOLA,WK,IWK,RMVLMT,SQRTG)                       ADS00400
      CALL ADS503 (NDV,VLBST,SVLB)                                      ADS00410
      CALL ADS503 (NDV,VUBST,SVUB)                                      ADS00420
C     ------------------------------------------------------------------ADS00430
C     CALCULATE NEW DESIGN VARIABLES.                                   ADS00440
C     ------------------------------------------------------------------ADS00450
      DO 140 I=1,NDV                                                    ADS00460
      X(I)=X(I)+XSUB(I)-BETA*XSUB(NDVR)*DF(I)                           ADS00470
      IF (X(I).GT.SVUB(I)) X(I)=SVUB(I)                                 ADS00480
      IF (X(I).LT.SVLB(I)) X(I)=SVLB(I)                                 ADS00490
140   CONTINUE                                                          ADS00500
      IF (IPRNT3.GE.4) WRITE (6,510)                                    ADS00510
C     CALCULATE THE VALUES OF OBJECTIVE FUNCTION AND CONSTRAINTS FOR NEWADS00520
C     DESIGN VARIABLES.                                                 ADS00530
      IGOTO=2                                                           ADS00540
      INFO=1                                                            ADS00550
      RETURN                                                            ADS00560
150   CONTINUE                                                          ADS00570
      IWRITE=0                                                          ADS00580
      OBJW=OBJ/SCFO                                                     ADS00590
      IF (IPRNT3.GT.1) WRITE (6,470) OBJW                               ADS00600
      ICALC=1                                                           ADS00610
      IF (IPRNT3.GT.1) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS00620
      ICALC=4                                                           ADS00630
      IF (IPRNT3.GT.2) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS00640
      IF (IPRNT3.GT.3) WRITE (6,490) XSUB(NDVR)                         ADS00650
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS00660
      IF (IPRNT3.GT.3) WRITE (6,520) GMAX                               ADS00670
C     IF NEW DESIGN IS INFEASIBLE , PUSH BACK THEM INTO FEASIBLE REGION.ADS00680
      IF (GMAX.LE.0.) GO TO 220                                         ADS00690
      IF (IGMXC.EQ.1) GO TO 310                                         ADS00700
      IWRITE=1                                                          ADS00710
      IF (BETA.GT.0.0001) GO TO 190                                     ADS00720
160   CONTINUE                                                          ADS00730
      DO 170 I=1,NDV                                                    ADS00740
      X(I)=0.5*(X(I)+XOLD(I))                                           ADS00750
170   CONTINUE                                                          ADS00760
      IGOTO=3                                                           ADS00770
      INFO=1                                                            ADS00780
      RETURN                                                            ADS00790
180   CONTINUE                                                          ADS00800
      IWRITE=1                                                          ADS00810
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS00820
      IF (GMAX.GT.0.) GO TO 160                                         ADS00830
      RMVLMT=0.5*RMVLMT                                                 ADS00840
      GO TO 220                                                         ADS00850
190   CONTINUE                                                          ADS00860
      DO 200 I=1,NDV                                                    ADS00870
      X(I)=X(I)+BETA*XSUB(NDVR)*DF(I)                                   ADS00880
200   CONTINUE                                                          ADS00890
      IGOTO=4                                                           ADS00900
      INFO=1                                                            ADS00910
      RETURN                                                            ADS00920
210   CONTINUE                                                          ADS00930
      IWRITE=1                                                          ADS00940
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS00950
      IF (GMAX.GT.0.) GO TO 160                                         ADS00960
C     OUTPUT IF REQUIRED.                                               ADS00970
220   CONTINUE                                                          ADS00980
      IF (IWRITE.EQ.0) GO TO 230                                        ADS00990
      IF (IPRNT3.GT.1) WRITE (6,440)                                    ADS01000
      ICALC=1                                                           ADS01010
      IF (IPRNT3.GT.1) CALL ADS402 (ICALC,X,NDV,SCFX)                   ADS01020
      OBJW=OBJ/SCFO                                                     ADS01030
      IF (IPRNT3.GT.1) WRITE (6,470) OBJW                               ADS01040
      IF (IPRNT3.GT.3) WRITE (6,490) XSUB(NDVR)                         ADS01050
      ICALC=4                                                           ADS01060
      IF (IPRNT3.GT.2) CALL ADS402 (ICALC,G,NCON,SCFG)                  ADS01070
      CALL ADS333 (NCON,IDG,G,GMAX)                                     ADS01080
      IF (IPRNT3.GT.3) WRITE (6,520) GMAX                               ADS01090
C     ------------------------------------------------------------------ADS01100
C     CHECK FOR CONVERGENCE.                                            ADS01110
C     ------------------------------------------------------------------ADS01120
230   CONTINUE                                                          ADS01130
      CALL ADS211 (OBJ,OBJ1,DELSTR,DABSTR,DEL1,DAB1,IDEL,IDAB)          ADS01140
      IF (JTER.GE.JTMAX) GO TO 240                                      ADS01150
      IF (IDEL.LT.ITRMST.AND.IDAB.LT.ITRMST) GO TO 290                  ADS01160
C     ------------------------------------------------------------------ADS01170
C     OPTIMIZATION IS COMPLETE.                                         ADS01180
C     ------------------------------------------------------------------ADS01190
240   CONTINUE                                                          ADS01200
      IF (IPRNT3.EQ.0) GO TO 280                                        ADS01210
      OBJW=OBJ/SCFO                                                     ADS01220
      WRITE (6,330) JTER,OBJW                                           ADS01230
C     X-VECTOR.                                                         ADS01240
      ICALC=1                                                           ADS01250
      CALL ADS402 (ICALC,X,NDV,SCFX)                                    ADS01260
C     G-VECTOR.                                                         ADS01270
      ICALC=4                                                           ADS01280
      CALL ADS402 (ICALC,G,NCON,SCFG)                                   ADS01290
C     IDENTIFY ACTIVE/VIOLATED CONSTRAINTS.                             ADS01300
      CALL ADS201 (NCON,NGT,G,IDG,ITMP,CT,CTMIN,CTL,CTLMIN,NAC,NVC)     ADS01310
C     SIDE CONSTRAINTS.                                                 ADS01320
      CALL ADS223 (NDV,NACS,CTLMIN,X,SVLB,SVUB,ISIDC)                   ADS01330
C     WRITE FINAL CONSTRAINT TOLERANCES.                                ADS01340
      WRITE (6,340) CT,CTL                                              ADS01350
C     WRITE ACTIVE/VIOLATED CONSTRAINT NUMBERS.                         ADS01360
      NAVC=NAC+NVC                                                      ADS01370
      WRITE (6,390) NAC,NVC                                             ADS01380
      IF (NAVC.GT.0) WRITE (6,400)                                      ADS01390
      IF (NAVC.GT.0) WRITE (6,410) (ITMP(I),I=1,NAVC)                   ADS01400
      WRITE (6,420) NACS                                                ADS01410
      IF (NACS.EQ.0) GO TO 270                                          ADS01420
      N=0                                                               ADS01430
      DO 260 I=1,NDV                                                    ADS01440
      IF (ISIDC(I).EQ.0) GO TO 260                                      ADS01450
      IF (ISIDC(I).EQ.1) GO TO 250                                      ADS01460
C     LOWER BOUND IS ACTIVE.                                            ADS01470
      N=N+1                                                             ADS01480
      ITMP(N)=-I                                                        ADS01490
250   IF (ISIDC(I).LT.0) GO TO 260                                      ADS01500
C     UPPER BOUND IS ACTIVE.                                            ADS01510
      N=N+1                                                             ADS01520
      ITMP(N)=I                                                         ADS01530
260   CONTINUE                                                          ADS01540
      WRITE (6,430)                                                     ADS01550
      WRITE (6,410) (ITMP(I),I=1,NACS)                                  ADS01560
270   WRITE (6,350)                                                     ADS01570
      IF (JTER.GE.JTMAX) WRITE (6,360)                                  ADS01580
      IF (IDEL.GE.ITRMST) WRITE (6,370) ITRMST                          ADS01590
      IF (IDAB.GE.ITRMST) WRITE (6,380) ITRMST                          ADS01600
280   CONTINUE                                                          ADS01610
      INFO=0                                                            ADS01620
      IGOTO=0                                                           ADS01630
      RETURN                                                            ADS01640
C                                                                       ADS01650
290   CONTINUE                                                          ADS01660
C     ------------------------------------------------------------------ADS01670
C     UPDATE MOVE LIMIT.                                                ADS01680
C     ------------------------------------------------------------------ADS01690
      RLOB=100.                                                         ADS01700
      ABOB=ABS(OBJ1)                                                    ADS01710
      IF (ABOB.LE.ZRO) GO TO 300                                        ADS01720
      RLOB=ABS(OBJ/OBJ1)                                                ADS01730
      RLOB=ABS(1.-RLOB)                                                 ADS01740
300   CONTINUE                                                          ADS01750
      ABOB=ABS(OBJ-OBJ1)                                                ADS01760
      CALL ADS323 (IMVCH,RLOB,ABOB,DELSTR,DABSTR,OBJ,OBJR,ZRO,RMVLMT)   ADS01770
310   CONTINUE                                                          ADS01780
      OBJR=OBJ1                                                         ADS01790
      OBJ1=OBJ                                                          ADS01800
      GO TO 30                                                          ADS01810
C     ------------------------------------------------------------------ADS01820
C     FORMATS.                                                          ADS01830
C     ------------------------------------------------------------------ADS01840
320   FORMAT (/,5X,'GRADIENT OF CONSTRAINT NUMBER',I5)                  ADS01850
330   FORMAT (////,5X,'FINAL OPTIMIZATION RESULTS',//,5X,               ADS01860
     1'NUMBER OF CONSTRAINED MINIMIZATIONS =',I5,//,5X,'OBJECTIVE =',   ADS01870
     21PE12.5)                                                          ADS01880
340   FORMAT (/,5X,'CONSTRAINT TOLORANCE, CT =',1PE12.5,'  CTL =',      ADS01890
     11PE12.5)                                                          ADS01900
350   FORMAT (/,5X,'TERMINATION CRITERIA')                              ADS01910
360   FORMAT (/,5X,'NUMBER OF ITERATIONS = JTMAX')                      ADS01920
370   FORMAT (/,5X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS01930
     1' CONSECUTIVE ITERATIONS')                                        ADS01940
380   FORMAT (/,5X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,     ADS01950
     1' CONSECUTIVE ITERATIONS')                                        ADS01960
390   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS AND',I5,         ADS01970
     1' VIOLATED CONSTRAINTS')                                          ADS01980
400   FORMAT (5X,'CONSTRAINT NUMBERS')                                  ADS01990
410   FORMAT (5X,10I5)                                                  ADS02000
420   FORMAT (/,5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')           ADS02010
430   FORMAT (5X,'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)')      ADS02020
440   FORMAT (/,5X,'RESULT OF RE-CALCULATION')                          ADS02030
450   FORMAT ('1',4X,75('*'),/,28X,'ISTRAT = 7: METHOD OF CENTERS',/,   ADS02040
     15X,75('*'))                                                       ADS02050
460   FORMAT (/,5X,'-- INITIAL DESIGN')                                 ADS02060
470   FORMAT (/,5X,'OBJ =',1PE12.5)                                     ADS02070
480   FORMAT (5X,15I5)                                                  ADS02080
490   FORMAT (/,5X,'RADIUS OF HYPERSPHERE =',1PE12.5)                   ADS02090
500   FORMAT (/,5X,'-- BEGIN MC ITERATION',I5,5X,'RMVLMT=',F8.5)        ADS02100
510   FORMAT (/,5X,'L.P. IS FINISHED')                                  ADS02110
520   FORMAT (/,5X,'GMAX=',1PE12.4)                                     ADS02120
530   FORMAT (/,5X,'RETAINED CONSTRAINTS (',I4,' )')                    ADS02130
      END                                                               ADS02140
C --- SUBROUTINE ADS331 (IOPT,IONED,NDV,NDVR,NGT,NGT1,X,D,A,NRA,        ADS02150
C    1GS,ICCT,G,DFCT,DLB,DUB,IDGCT,NCOLA,WK,IWK,RMVLMT,SQRTG)           ADS02160
C     DIMENSION D(NDV+1),A(NRA,NCOLA),GS(NGT+1),ICCT(NGT+1),G(NGT),     ADS02170
C    1DFCT(NDV+1),IDGCT(NGT+1),WK(NRWK),IWK(NRIWK),DLB(NDV+1),X(NDV),   ADS02180
C    2DUB(NDV+1),SQRTG(NGT+1)                                           ADS02190
C     ------------------------------------------------------------------ADS02200
C --- ROUTINE TO SOLVE THE LINEAR PROGRAMMING SUB-PROBLEM FOR METHOD OF ADS02210
C     CENTERS.                                                          ADS02220
C     ------------------------------------------------------------------ADS02230
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                APRIL 1983   ADS02240
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS02250
C --- INPUT                                                             ADS02260
C     IOPT  - OPTIMIZATION IDENTIFIER.                                  ADS02270
C     IONED - ONE-DIMENSIONAL SEARCH IDENTIFIER.                        ADS02280
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS02290
C     NDVR  - NDV + 1.                                                  ADS02300
C     NGT   - TOTAL NUMBER OF GRADIENTS CALCULATED.                     ADS02310
C     NGT1  - NGT + 1. IF VIOLATED , NGT.                               ADS02320
C     X     - CURRENT VALUE OF PRIMAL DESIGN VARIABLE.                  ADS02330
C     D     - DESIGN VARIABLES OF SUB-PROBLEM.                          ADS02340
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS02350
C             PROBLEM.                                                  ADS02360
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS02370
C     GS    - CURRENT CONSTRAINT VALUE.                                 ADS02380
C     ICCT  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR L.P. SUB-PROBLEM.ADS02390
C     G     - CONSTRAINT VALUE OF SUB-PROBLEM.                          ADS02400
C     DFCT  - 'DF' FOR SUB=PROBLEM                                      ADS02410
C     DLB,DUB - UPPER AND LOWER BOUNDS ON DESIGN VARIABLES OF SUB-      ADS02420
C               PROBLEM.                                                ADS02430
C     IDGCT - CONSTRAINT TYPE IDENTIFIER FOR SUB-PROBLEM.               ADS02440
C     NCOLA - NUMBER OF COLUMNS OF ARRAY A.                             ADS02450
C     WK    - REAL WORK ARRAY.                                          ADS02460
C     IWK   - INTEGER WORK ARRAY.                                       ADS02470
C     RMVLMT- MOVE LIMIT.                                               ADS02480
C     SQRTG - SQUARE ROOT MEAN OF GRADIENTS OF CONSTRAINT AND OBJECTIVE.ADS02490
C --- OUTPUT                                                            ADS02500
C     D     - OPTIMUM VALUE OF DESIGN VARIABLES OF SUB-PROBLEM.         ADS02510
C                                                                       ADS02520
      SUBROUTINE ADS331 (IOPT,IONED,NDV,NDVR,NGT,NGT1,X,D,A,NRA,GS,ICCT,ADS02530
     1G,DFCT,DLB,DUB,IDGCT,NCOLA,WK,IWK,RMVLMT,SQRTG)                   ADS02540
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION D(1),A(NRA,1),GS(1),ICCT(1),G(1),DFCT(1),IDGCT(1),WK(1),ADS02550
     1IWK(1),DLB(1),DUB(1),SQRTG(1),X(1)                                ADS02560
C                                                                       ADS02570
      INFO=0                                                            ADS02580
      ISTRAT=7                                                          ADS02590
      DO 10 I=1,NDV                                                     ADS02600
      ABSX=ABS(X(I))                                                    ADS02610
      IF (ABSX.LE.1.0) ABSX=1.0                                         ADS02620
      XLIM=RMVLMT*ABSX                                                  ADS02630
      DD=DLB(I)-X(I)                                                    ADS02640
      DLB(I)=DMAX1(-XLIM,DD)                                            ADS02650
      DD=DUB(I)-X(I)                                                    ADS02660
      DUB(I)=DMIN1(XLIM,DD)                                             ADS02670
10    CONTINUE                                                          ADS02680
      DLB(NDVR)=0.0                                                     ADS02690
      DUB(NDVR)=1.0E20                                                  ADS02700
      DO 20 I=1,NGT                                                     ADS02710
      ICCT(I)=I                                                         ADS02720
      IF (IDGCT(I).GE.0) IDGCT(I)=2                                     ADS02730
      IF (IDGCT(I).LT.0) IDGCT(I)=-2                                    ADS02740
20    CONTINUE                                                          ADS02750
      IF (NGT1.GT.NGT) IDGCT(NGT1)=2                                    ADS02760
      IF (NGT1.GT.NGT) ICCT(NGT1)=NGT1                                  ADS02770
      DO 30 I=1,NDVR                                                    ADS02780
      D(I)=0.                                                           ADS02790
30    CONTINUE                                                          ADS02800
C     BEGIN OPTIMIZATION OF LP-SUBPROBLEM.                              ADS02810
40    CONTINUE                                                          ADS02820
      FCT=-D(NDVR)                                                      ADS02830
      DO 60 K=1,NGT1                                                    ADS02840
      J=ICCT(K)                                                         ADS02850
      G(J)=GS(J)+D(NDVR)                                                ADS02860
      DO 50 I=1,NDV                                                     ADS02870
      G(J)=G(J)+A(I,K)*D(I)                                             ADS02880
50    CONTINUE                                                          ADS02890
60    CONTINUE                                                          ADS02900
      IF (INFO.EQ.1) GO TO 80                                           ADS02910
      DO 70 I=1,NDV                                                     ADS02920
      DFCT(I)=0.                                                        ADS02930
70    CONTINUE                                                          ADS02940
      DFCT(NDVR)=-1.                                                    ADS02950
80    CONTINUE                                                          ADS02960
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDVR,NGT1,FCT,D,DLB,DUB,G,IDGCADS02970
     1T,NGT1,ICCT,DFCT,A,NRA,NCOLA,WK,IWK)                              ADS02980
      IF (INFO.EQ.0) RETURN                                             ADS02990
      GO TO 40                                                          ADS03000
      END                                                               ADS03010
C --- SUBROUTINE ADS333 (NCON,IDG,G,GMAX)                               ADS03020
C     DIMENSION G(NCON),IDG(NCON)                                       ADS03030
C     ------------------------------------------------------------------ADS03040
C --- ROUTINE TO CALCULATE THE MAXIMUM VALUE OF CONSTRAINTS.            ADS03050
C     ------------------------------------------------------------------ADS03060
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                APRIL 1983   ADS03070
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS03080
C --- INPUT                                                             ADS03090
C      NCON - NUMBER OF CONSTRAINTS.                                    ADS03100
C      IDG  - CONSTRAINT TYPE IDENTIFIER.                               ADS03110
C      G    - CONSTRAINT VALUES.                                        ADS03120
C --- OUTPUT                                                            ADS03130
C      GMAX - MAXIMUM VIOLATED VALUE IN CONSTRAINT VALUES.              ADS03140
C                                                                       ADS03150
      SUBROUTINE ADS333 (NCON,IDG,G,GMAX)                               ADS03160
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IDG(1)                                             ADS03170
      GMAX=-100.                                                        ADS03180
      DO 10 J=1,NCON                                                    ADS03190
      GJ=G(J)                                                           ADS03200
      IF (IDG(J).LT.0) GJ=ABS(GJ)                                       ADS03210
      IF (GJ.LE.GMAX) GO TO 10                                          ADS03220
      JMAX=J                                                            ADS03230
      GMAX=GJ                                                           ADS03240
10    CONTINUE                                                          ADS03250
      IF (IDG(JMAX).LT.0) GMAX=G(JMAX)                                  ADS03260
      RETURN                                                            ADS03270
      END                                                               ADS03280
C --- SUBROUTINE ADS371 (IOPT,IONED,NDV,NGT,X,D,A,NRA,B,DF,GZ,          ADS03290
C    1ICQP,GSUB,DQ,DLB,DUB,IDGQP,NCOLA,WK,IWK,P)                        ADS03300
C     DIMENSION D(NDV),A(NRA,NCOLA),B(NDV*)NDV+1)/2),DF(NDV),GZ(NCOLA)  ADS03310
C    1,ICQP(NCOLA),GSUB(NCOLA),DQ(NDV),IDGQP(NCOLA),WK(NRWK),IWK(NRIWK),ADS03320
C    2DLB(NDV),DUB(NDV),P(NDV),X(NDV)                                   ADS03330
C     ------------------------------------------------------------------ADS03340
C --- ROUTINE TO SOLVE THE QUADRATIC PROGRAMMING SUB-PROBLEM.           ADS03350
C     ------------------------------------------------------------------ADS03360
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 FEB. 1983   ADS03370
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS03380
C --- INPUT                                                             ADS03390
C     IOPT  - OPTIMIZATION IDENTIFIER.                                  ADS03400
C     IONED - ONE-DIMENSIONAL SEARCH IDENTIFIER.                        ADS03410
C     NDV   - NUMBER OF DESIGN VARIABLES.                               ADS03420
C     NGT   - NUMBER OF GRADIENTS CALCULATED.                           ADS03430
C     X     - CURRENT VALUE OF PROMAL DESIGN VARIABLE.                  ADS03440
C     D     - DESIGN VARIANLES OF SUB-PROBLEM.                          ADS03450
C     A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS OF PRIMALADS03460
C             PROBLEM.                                                  ADS03470
C     NRA   - NUMBER OF ROWS IN ARRAY A.                                ADS03480
C     B     - HESSIAN MATRIX.                                           ADS03490
C     DF    - GRADIENT OF PRIMAL OBJECTIVE FUNCTION.                    ADS03500
C     GZ    - PSAI TIMES CONSTRAINT VALUE.                              ADS03510
C     ICQP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS FOR Q.P. SUB-PROBLEM.ADS03520
C     GSUB- CONSTRAINT VALUE OF SUB-PROBLEM.                            ADS03530
C     DQ    - 'DF' FOR SUB-PROBLEM                                      ADS03540
C     DLB,DUB - UPPER AND LOWER LIMITS FOR DESIGN VARIABLES OF SUB-     ADS03550
C               PROBLEM.                                                ADS03560
C     IDGQP - CONSTRAINT TYPE IDENTIFIER FOR SUB-PROBLEM.               ADS03570
C     NCOLA - NUMBER OF COLUMNS OF ARRAY A.                             ADS03580
C     WK    - REAL WORK ARRAY.                                          ADS03590
C     IWK   - INTEGER WORK ARRAY.                                       ADS03600
C     P     - WORK ARRAY.                                               ADS03610
C --- OUTPUT                                                            ADS03620
C     D     - OPTIMUM VALUE OF DESIGN VARIABLES OF SUB-PROBLEM.         ADS03630
C                                                                       ADS03640
      SUBROUTINE ADS371 (IOPT,IONED,NDV,NGT,X,D,A,NRA,B,DF,GZ,ICQP,GSUB,ADS03650
     1DQ,DLB,DUB,IDGQP,NCOLA,WK,IWK,P)                                  ADS03660
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION D(1),A(NRA,1),B(1),DF(1),GZ(1),ICQP(1),GSUB(1),DQ(1),IDGADS03670
     1QP(1),WK(1),IWK(1),DLB(1),DUB(1),P(1),X(1)                        ADS03680
C                                                                       ADS03690
      INFO=0                                                            ADS03700
      ISTRAT=8                                                          ADS03710
      DO 10 I=1,NDV                                                     ADS03720
      DLB(I)=DLB(I)-X(I)                                                ADS03730
      DUB(I)=DUB(I)-X(I)                                                ADS03740
10    CONTINUE                                                          ADS03750
      DO 20 I=1,NGT                                                     ADS03760
      ICQP(I)=I                                                         ADS03770
      IF (IDGQP(I).GE.0) IDGQP(I)=2                                     ADS03780
      IF (IDGQP(I).LT.0) IDGQP(I)=-2                                    ADS03790
20    CONTINUE                                                          ADS03800
C     CALCULATION OF INITIAL VALUE.                                     ADS03810
      DO 30 I=1,NDV                                                     ADS03820
      D(I)=0.0                                                          ADS03830
      IF (D(I).LT.DLB(I)) D(I)=DLB(I)                                   ADS03840
      IF (D(I).GT.DUB(I)) D(I)=DUB(I)                                   ADS03850
30    CONTINUE                                                          ADS03860
C     BEGIN OPTIMIZATION OF QP-SUBPROBLEM.                              ADS03870
40    CONTINUE                                                          ADS03880
      CALL ADS504 (NDV,B,D,P)                                           ADS03890
      CALL ADS501 (NDV,D,P,Q)                                           ADS03900
      CALL ADS501 (NDV,D,DF,R)                                          ADS03910
      Q=R+0.5*Q                                                         ADS03920
      DO 60 K=1,NGT                                                     ADS03930
      J=ICQP(K)                                                         ADS03940
      GSUB(J)=GZ(J)                                                     ADS03950
      DO 50 I=1,NDV                                                     ADS03960
50    GSUB(J)=GSUB(J)+A(I,K)*D(I)                                       ADS03970
60    CONTINUE                                                          ADS03980
      IF (INFO.EQ.1) GO TO 80                                           ADS03990
      DO 70 I=1,NDV                                                     ADS04000
      DQ(I)=DF(I)+P(I)                                                  ADS04010
70    CONTINUE                                                          ADS04020
80    CONTINUE                                                          ADS04030
      CALL ADS002 (INFO,ISTRAT,IOPT,IONED,NDV,NGT,Q,D,DLB,DUB,GSUB,IDGQPADS04040
     1,NGT,ICQP,DQ,A,NRA,NCOLA,WK,IWK)                                  ADS04050
      IF (INFO.EQ.0) RETURN                                             ADS04060
      GO TO 40                                                          ADS04070
      END                                                               ADS04080
C --- SUBROUTINE ADS375 (NCON,IDG,G,GMULT,PENAL)                        ADS04090
C     DIMENSION IDG(NCON),G(NCON),GMULT(NCOLA)                          ADS04100
C     ------------------------------------------------------------------ADS04110
C --- ROUTINE TO CALCULATE THE TEMPORARY OBJECTIVE FUNCTION.            ADS04120
C     ------------------------------------------------------------------ADS04130
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 FEB. 1983   ADS04140
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS04150
C --- INPUT                                                             ADS04160
C     NCON  - NUMBER OF CONSTRAINTS.                                    ADS04170
C     IDG   - CONSTRAINT TYPE IDENTIFIER.                               ADS04180
C     G     - VALUES OF CONSTRAINTS.                                    ADS04190
C     GMULT - PENALTY PARAMETER.                                        ADS04200
C --- OUTPUT                                                            ADS04210
C     PENAL - VALUE OF PENALTY TERM.                                    ADS04220
C                                                                       ADS04230
      SUBROUTINE ADS375 (NCON,IDG,G,GMULT,PENAL)                        ADS04240
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION IDG(1),G(1),GMULT(1)                                    ADS04250
C                                                                       ADS04260
      PENAL=0.                                                          ADS04270
      DO 30 J=1,NCON                                                    ADS04280
      IF (IDG(J).LT.0) GO TO 10                                         ADS04290
      IF (G(J).LE.0.) GO TO 30                                          ADS04300
      GM=G(J)                                                           ADS04310
      GO TO 20                                                          ADS04320
10    CONTINUE                                                          ADS04330
      GM=ABS(G(J))                                                      ADS04340
20    CONTINUE                                                          ADS04350
      PENAL=PENAL+GMULT(J)*GM                                           ADS04360
30    CONTINUE                                                          ADS04370
      RETURN                                                            ADS04380
      END                                                               ADS04390
C --- SUBROUTINE ADS376 (NDV,NGT,G,IDG,IC,ICQP,DF,GMULT,A,NRA,S)        ADS04400
C     DIMENSION IC(NCOLA),DF(NDV),GMULT(NCOLA),A(NRA,NCOLA),S(NDV),     ADS04410
C    1G(NCON),IDG(NCON),ICQP(NCOLA)                                     ADS04420
C     ------------------------------------------------------------------ADS04430
C --- ROUTINE TO CALLCULATE THE GRADIENT OF PSEUDO-OBJECTIVE FOR ONE-D. ADS04440
C     SEARCH IN POWELL'S METHOD.                                        ADS04450
C     ------------------------------------------------------------------ADS04460
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 AUG. 1983   ADS04470
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS04480
C --- INPUT                                                             ADS04490
C      NDV   - NUMBER OF DESIGN VARIABLES.                              ADS04500
C      NGT  - NUMBER OF ACTIVE CONSTRAINTS.                             ADS04510
C      G     - CONSTRAINTS VALUES.                                      ADS04520
C      IDG   - CONSTRAINT TYPE IDENTIFIER.                              ADS04530
C      IC    - ARRAY IDENTIFYING ACTIVE CONSTRAINTS.                    ADS04540
C      ICQP  - ARRAY IDENTIFYING ACTIVE CONSTRAINTS OF QP-SUBPROBLEM.   ADS04550
C      DF    - GRADIENT OF PRIMAL OBJECTIVE.                            ADS04560
C      GMULT - PENALTY PARAMETER.                                       ADS04570
C      A     - ARRAY CONTAINING GRADIENTS OF CONSTRAINTS.               ADS04580
C      NRA   - NUMBER OF ROWS IN ARRAY A.                               ADS04590
C --- OUTPUT                                                            ADS04600
C      S     - GRADIENT OF PSEUDO OBJECTIVE.                            ADS04610
      SUBROUTINE ADS376 (NDV,NGT,G,IDG,IC,ICQP,DF,GMULT,A,NRA,S)        ADS04620
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IDG(1),IC(1),ICQP(1),DF(1),GMULT(1),A(NRA,1),S(1)  ADS04630
C                                                                       ADS04640
      DO 10 I=1,NDV                                                     ADS04650
      S(I)=DF(I)                                                        ADS04660
10    CONTINUE                                                          ADS04670
      IF (NGT.EQ.0) RETURN                                              ADS04680
      DO 20 K=1,NGT                                                     ADS04690
      KK=ICQP(K)                                                        ADS04700
      J=IC(KK)                                                          ADS04710
      IF (IDG(J).GE.0.AND.G(J).LT.0.) GO TO 20                          ADS04720
      V=GMULT(J)                                                        ADS04730
      IF (G(J).LT.0.) V=-V                                              ADS04740
      CALL ADS508 (NDV,K,V,A,NRA,S)                                     ADS04750
20    CONTINUE                                                          ADS04760
      RETURN                                                            ADS04770
      END                                                               ADS04780
C --- SUBROUTINE ADS377 (NDV,P,Y,B,S)                                   ADS04790
C     DIMENSION P(NDV),Y(NDV),B(NDV*(NDV+1)/2),S(NDV)                   ADS04800
C     ------------------------------------------------------------------ADS04810
C --- ROUTINE FOR CALCULATION OF TRANSFORMED DIFFERENCE BETWEEN FORMER  ADS04820
C     AND CURRENT GRADIENTS OF OBJECTIVE FUNCTION ('ETA' IN POWELL'S    ADS04830
C     PAPER).                                                           ADS04840
C     ------------------------------------------------------------------ADS04850
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 FEB. 1983   ADS04860
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS04870
C --- INPUT                                                             ADS04880
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS04890
C     P   - TRUE DIFFERENCE BETWEEN FORMER AND CURRENT GRADIENTS OF     ADS04900
C           OBJECTIVE FUNCTION.                                         ADS04910
C     Y    - DIFFERENCE BETWEEN FORMER AND CURRENT DESIGN VARIABLES.    ADS04920
C     B    - FORMER APPROXIMATION TO HESSIAN MATRIX. UPPER TRIANGLE OF  ADS04930
C            NDV X NDV MATRIX , STORED COLUMN BY COLUMN.                ADS04940
C     S    - ARRAY FOR WORK (NDV).                                      ADS04950
C --- OUTPUT                                                            ADS04960
C     P    - TRANSFORMED DIFFERENCE ('ETA').                            ADS04970
C                                                                       ADS04980
      SUBROUTINE ADS377 (NDV,P,Y,B,S)                                   ADS04990
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION P(1),Y(1),B(1),S(1)                                     ADS05000
C                                                                       ADS05010
      CALL ADS504 (NDV,B,Y,S)                                           ADS05020
      CALL ADS501 (NDV,Y,S,DBD)                                         ADS05030
      CALL ADS501 (NDV,Y,P,DG)                                          ADS05040
      DBD2=0.2*DBD                                                      ADS05050
      IF (DG.GE.DBD2) RETURN                                            ADS05060
      SITA=0.8*DBD/(DBD-DG)                                             ADS05070
      TITA=1.-SITA                                                      ADS05080
      DO 10 I=1,NDV                                                     ADS05090
      P(I)=SITA*P(I)+TITA*S(I)                                          ADS05100
10    CONTINUE                                                          ADS05110
      RETURN                                                            ADS05120
      END                                                               ADS05130
C --- SUBROUTINE ADS378 (NDV,P,Y,B,S,ZRO,ILL)                           ADS05140
C     DIMENSION P(NDV),Y(NDV),B(NDV*(NDV+1)/2),S(NDV)                   ADS05150
C     ------------------------------------------------------------------ADS05160
C --- ROUTINE FOR REVISION OF HESSIAN MATRIX IN PROJECTED LAGRANGIAN    ADS05170
C     METHOD.                                                           ADS05180
C                                                                       ADS05190
C   ------------------------------------------------------------------  ADS05200
C                                                                       ADS05210
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 FEB. 1983   ADS05220
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS05230
C --- INPUT                                                             ADS05240
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS05250
C     P   - ETA (SEE SUB. ADS377).                                      ADS05260
C     Y   - DIFFERENCE BETWEEN FORMER AND CURRENT DESIGN VARIABLES.     ADS05270
C     B   - FORMER APPROXIMATION TO HESSIAN MATRIX. UPPER TRIANGLE OF   ADS05280
C           NDV X NDV MATRIX , STORED COLUMN BY COLUMN.                 ADS05290
C     S   - ARRAY FOR WORK (NDV).                                       ADS05300
C --- OUTPUT                                                            ADS05310
C     B   - NEW APPROXIMATION TO HESSIAN MATRIX.                        ADS05320
C     ILL - RETURN CONDITION. =0 NORMAL , =1 ABNORMAL.                  ADS05330
C                                                                       ADS05340
      SUBROUTINE ADS378 (NDV,P,Y,B,S,ZRO,ILL)                           ADS05350
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION P(1),Y(1),B(1),S(1)                                     ADS05360
C                                                                       ADS05370
      ILL=1                                                             ADS05380
      CALL ADS504 (NDV,B,Y,S)                                           ADS05390
      CALL ADS501 (NDV,P,Y,SIGMA)                                       ADS05400
      CALL ADS501 (NDV,Y,S,TAU)                                         ADS05410
      SH=ABS(SIGMA)                                                     ADS05420
      ST=ABS(TAU)                                                       ADS05430
      IF (SH.LT.ZRO.OR.ST.LT.ZRO) RETURN                                ADS05440
      SIGMA=1./SIGMA                                                    ADS05450
      TAU=1./TAU                                                        ADS05460
C                                                                       ADS05470
      NA=0                                                              ADS05480
      DO 10 J=1,NDV                                                     ADS05490
      DO 10 I=1,J                                                       ADS05500
      NA=NA+1                                                           ADS05510
      B(NA)=B(NA)+SIGMA*P(I)*P(J)-TAU*S(I)*S(J)                         ADS05520
10    CONTINUE                                                          ADS05530
      ILL=0                                                             ADS05540
      RETURN                                                            ADS05550
      END                                                               ADS05560



C --- SUBROUTINE ADS401 (A,N,SCAL,ISCAL)                                ADS05570
C     DIMENSION A(N),SCAL(N),B(5)                                       ADS05580
C     ------------------------------------------------------------------ADS05590
C --- ROUTINE TO PRINT ARRAY A ON OUTPUT FILE 6.                        ADS05600
C     ------------------------------------------------------------------ADS05610
C     BY L. E. MADSEN                                      APRIL, 1981  ADS05620
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS05630
C                                                                       ADS05640
C --- INPUT                                                             ADS05650
C     A - VECTOR CONTAINING PARAMETERS TO PRINTED.  MINIMUM             ADS05660
C        DIMENSION IS N.                                                ADS05670
C     N - NUMBER OF PARAMETERS IN A WHICH ARE TO BE WRITTEN.            ADS05680
C     SCAL - SCAL FACTORS.                                              ADS05690
C     ISCAL - SCALING ID.  IF 0 DO NOT UNSCALE VECTOR A.                ADS05700
C                                                                       ADS05710
C --- OUTPUT - NONE TO CALLING PROGRAM.  ARGUMENTS REMAIN UNCHANGED.    ADS05720
C                                                                       ADS05730
      SUBROUTINE ADS401 (A,N,SCAL,ISCAL)                                ADS05740
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(1),SCAL(1),B(5)                                       ADS05750
C     WRITE THE ARRAY.                                                  ADS05760
C     WRITE FIVE VALUES PER ROW.                                        ADS05770
      IMAX=(N-1)/5+1                                                    ADS05780
      I1=1                                                              ADS05790
      DO 10 I=1,IMAX                                                    ADS05800
      I2=I1+4                                                           ADS05810
      I2=MIN0(I2,N)                                                     ADS05820
      J1=1                                                              ADS05830
      DO 20 J=I1,I2                                                     ADS05840
      B(J1)=A(J)                                                        ADS05850
      IF(ISCAL.GT.0) B(J1)=B(J1)/SCAL(J)                                ADS05860
20    J1=J1+1                                                           ADS05870
      NN=I2-I1+1                                                        ADS05880
      CALL ADS405(I1,B,NN)                                              ADS05890
      I1=I1+5                                                           ADS05900
10    CONTINUE                                                          ADS05910
      RETURN                                                            ADS05920
      END                                                               ADS05930
C --- SUBROUTINE ADS402(ICALC,A,N,SCAL)                                 ADS05940
C     DIMENSION A(N),SCAL(N)                                            ADS05950
C     ------------------------------------------------------------------ADS05960
C --- ROUTINE TO WRITE ARRAY A ON OUTPUT FILE 6.  ICALC DETERMINES      ADS05970
C     THE TITLE TO BE PRINTED.                                          ADS05980
C     ------------------------------------------------------------------ADS05990
C     BY G. N. VANDERPLAATS                            JAN 1981         ADS06000
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS06010
C                                                                       ADS06020
C --- INPUT                                                             ADS06030
C     ICALC - CALCULATION CONTROL.                                      ADS06040
C         1 - A IS THE X-VECTOR.                                        ADS06050
C         2 - A IS THE S-VECTOR.                                        ADS06060
C         3 - A IS THE DF-VECTOR.                                       ADS06070
C         4 - A IS THE G-VECTOR.                                        ADS06080
C         5 - A IS THE VLB-VECTOR.                                      ADS06090
C         6 - A IS THE VUB-VECTOR.                                      ADS06100
C         7 - A IS THE SCALING VECTOR, SCAL.                            ADS06110
C         8 - A IS THE DX-VECTOR.                                       ADS06120
C         9 - A IS THE DG-VECTOR.                                       ADS06130
C        10 - A IS THE ALAMDA-VECTOR.                                   ADS06140
C        11 - A IS THE GSCALE-VECTOR.                                   ADS06150
C        12 - A IS THE DF-VECTOR OF TEMPORARY OBJECTIVE.                ADS06160
C        13 - A IS THE DG/D-ALPHA VECTOR.                               ADS06170
C     A - VECTOR CONTAINING PARAMETERS TO BE WRITTEN ON UNIT 6.         ADS06180
C         MINIMUM DIMENSION IS N.                                       ADS06190
C     N - NUMBER OF PARAMETERS CONTAINED IN A WHICH ARE TO BE WRITTEN.  ADS06200
C     SCAL - SCALING VECTOR USED WITH X, VLB AND VUB PRINT. MINIMUM     ADS06210
C         DIMENSION IS N.                                               ADS06220
C                                                                       ADS06230
C --- OUTPUT - NONE TO CALLING PROGRAM. ICALC, A, AND N ARE UNCHANGED.  ADS06240
C                                                                       ADS06250
C --- EXTERNALS                                                         ADS06260
C     ADS401, ADS510                                                    ADS06270
C                                                                       ADS06280
      SUBROUTINE ADS402 (ICALC,A,N,SCAL)                                ADS06290
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(1),SCAL(1)                                            ADS06300
      JCALC=3                                                           ADS06310
      ISCAL=0                                                           ADS06320
      GO TO (10,20,30,40,50,60,70,80,90,110,120,100,115), ICALC         ADS06330
C     X-VECTOR                                                          ADS06340
10    WRITE (6,160)                                                     ADS06350
      ISCAL=1                                                           ADS06360
      GO TO 130                                                         ADS06370
C     S-VECTOR.                                                         ADS06380
20    WRITE (6,170)                                                     ADS06390
      GO TO 130                                                         ADS06400
C     DF-VECTOR.                                                        ADS06410
30    WRITE (6,190)                                                     ADS06420
      GO TO 130                                                         ADS06430
C     G-VECTOR.                                                         ADS06440
40    WRITE (6,180)                                                     ADS06450
      ISCAL=1                                                           ADS06460
      GO TO 130                                                         ADS06470
C     VLB-VECTOR.                                                       ADS06480
50    WRITE (6,140)                                                     ADS06490
      ISCAL=1                                                           ADS06500
      GO TO 130                                                         ADS06510
C     VUB-VECTOR.                                                       ADS06520
60    WRITE (6,150)                                                     ADS06530
      ISCAL=1                                                           ADS06540
      GO TO 130                                                         ADS06550
70    WRITE (6,200)                                                     ADS06560
      GO TO 130                                                         ADS06570
C     DX-VECTOR                                                         ADS06580
80    WRITE (6,210)                                                     ADS06590
      GO TO 130                                                         ADS06600
C     DG-VECTOR                                                         ADS06610
90    WRITE (6,220)                                                     ADS06620
      GO TO 130                                                         ADS06630
100   WRITE (6,250)                                                     ADS06640
      GO TO 130                                                         ADS06650
C     ALAMDA-VECTOR                                                     ADS06660
110   WRITE (6,230)                                                     ADS06670
      GO TO 130                                                         ADS06680
C     GSCALE-VECTOR                                                     ADS06690
120   WRITE (6,240)                                                     ADS06700
      GO TO 130                                                         ADS06710
115   WRITE (6,255)                                                     ADS06720
130   CONTINUE                                                          ADS06730
C     WRITE THE ARRAY.                                                  ADS06740
      CALL ADS401 (A,N,SCAL,ISCAL)                                      ADS06750
      RETURN                                                            ADS06760
140   FORMAT (/,5X,'LOWER BOUNDS ON THE DECISION VARIABLES (VLB-VECTOR)'ADS06770
     1)                                                                 ADS06780
150   FORMAT (/,5X,'UPPER BOUNDS ON THE DECISION VARIABLES (VUB-VECTOR)'ADS06790
     1)                                                                 ADS06800
160   FORMAT (/,5X,'DECISION VARIABLES (X-VECTOR)')                     ADS06810
170   FORMAT (/,5X,'SEARCH DIRECTION (S-VECTOR)')                       ADS06820
180   FORMAT (/,5X,'CONSTRAINT VALUES (G-VECTOR)')                      ADS06830
190   FORMAT (/,5X,'GRADIENT OF THE OBJECTIVE FUNCTION (DF-VECTOR)')    ADS06840
200   FORMAT (/,5X,'SCALING VECTOR (SCAL)')                             ADS06850
210   FORMAT (/,5X,'DIFFERENCE BETWEEN X AND XOLD')                     ADS06860
220   FORMAT (/,5X,'DIFFERENCE BETWEEN DF AND DFOLD')                   ADS06870
230   FORMAT (/,5X,'LAGRANGE MULTIPLIERS')                              ADS06880
240   FORMAT (/,5X,'PENALTY FACTOR OF EACH CONSTRAINT')                 ADS06890
250   FORMAT (/,5X,'GRADIENT OF THE TEMPORARY OBJECTIVE FUNCTION')      ADS06900
255   FORMAT (/,5X,'GRADIENTS OF CONSTRAINTS WITH RESPECT TO ALPHA',/,  ADS06910
     15X,'(IGNORE ZERO VALUES)')                                        ADS06920
      END                                                               ADS06930
C --- SUBROUTINE ADS403 (WK,NRWK,IWK,NRIWK)                             ADS06940
C     DIMENSION WK(NRWK),IWK(NRIWK)                                     ADS06950
C     ------------------------------------------------------------------ADS06960
C --- ROUTINE TO PRINT VALUES OF CONTROL PARAMETERS AND ARRAY STORAGE.  ADS06970
C     ------------------------------------------------------------------ADS06980
C     BY G. N. VANDERPLAATS                                AUG. 1982.   ADS06990
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS07000
C                                                                       ADS07010
C --- INPUT                                                             ADS07020
C     WK - REAL WORK ARRAY.  SCALAR PARAMETERS STORED AT THE BEGINNING  ADS07030
C         OF WK.  DIMENSION IS NRWK.                                    ADS07040
C  ETC.                                                                 ADS07050
C                                                                       ADS07060
C --- OUTPUT.                                                           ADS07070
C     NONE TO CALLING PROGRAM.  PRINTED OUTPUT TO FILE 6.               ADS07080
C                                                                       ADS07090
      SUBROUTINE ADS403 (WK,NRWK,IWK,NRIWK)                             ADS07100
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION WK(NRWK),IWK(NRIWK)                                     ADS07110
C     ------------------------------------------------------------------ADS07120
C     SCALAR PARAMETERS.                                                ADS07130
C     ------------------------------------------------------------------ADS07140
C     REAL PARAMETERS.                                                  ADS07150
      WRITE (6,10)                                                      ADS07160
      WRITE (6,20) (WK(I),WK(I+19),I=1,10)                              ADS07170
      WRITE (6,30) (WK(I),WK(I+19),I=11,19)                             ADS07180
C                                                                       ADS07190
C     INTEGER PARAMETERS.                                               ADS07200
      WRITE (6,40)                                                      ADS07210
      WRITE (6,50) (IWK(I),IWK(I+3),IWK(I+5),I=1,2),IWK(3)              ADS07220
      WRITE (6,60) NRWK,IWK(150),NRIWK,IWK(170)                         ADS07230
      RETURN                                                            ADS07240
10    FORMAT ('1',//,5X,'SCALAR PROGRAM PARAMETERS',/,                  ADS07250
     15X,'REAL PARAMETERS')                                             ADS07260
20    FORMAT (5X,'  1) ALAMDZ =',1PE13.5,4X,' 20) EXTRAP =',1PE13.5,/,  ADS07270
     15X,'  2) BETAMC =',1PE13.5,4X,' 21) FDCH   =',1PE13.5,/,5X,       ADS07280
     2'  3) CT     =',1PE13.5,4X,' 22) FDCHM  =',1PE13.5,/,5X,          ADS07290
     3'  4) CTL    =',1PE13.5,4X,' 23) GMULTZ =',1PE13.5,/,5X,          ADS07300
     4'  5) CTLMIN =',1PE13.5,4X,' 24) PSAIZ  =',1PE13.5,/,5X,          ADS07310
     5'  6) CTMIN  =',1PE13.5,4X,' 25) RMULT  =',1PE13.5,/,5X,          ADS07320
     6'  7) DABALP =',1PE13.5,4X,' 26) RMVLMZ =',1PE13.5,/,5X,          ADS07330
     7'  8) DABOBJ =',1PE13.5,4X,' 27) RP     =',1PE13.5,/,5X,          ADS07340
     8'  9) DABOBM =',1PE13.5,4X,' 28) RPMAX  =',1PE13.5,/,5X,          ADS07350
     9' 10) DABSTR =',1PE13.5,4X,' 29) RPMULT =',1PE13.5)               ADS07360
30    FORMAT (5X,' 11) DELALP =',1PE13.5,4X,' 30) RPPMIN =',1PE13.5,/,  ADS07370
     15X,' 12) DELOBJ =',1PE13.5,4X,' 31) RPPRIM =',1PE13.5,/,5X,       ADS07380
     2' 13) DELOBM =',1PE13.5,4X,' 32) SCFO   =',1PE13.5,/,5X,          ADS07390
     3' 14) DELSTR =',1PE13.5,4X,' 33) SCLMIN =',1PE13.5,/,5X,          ADS07400
     4' 15) DLOBJ1 =',1PE13.5,4X,' 34) STOL   =',1PE13.5,/,5X,          ADS07410
     5' 16) DLOBJ2 =',1PE13.5,4X,' 35) THETAZ =',1PE13.5,/,5X,          ADS07420
     6' 17) DX1    =',1PE13.5,4X,' 36) XMULT  =',1PE13.5,/,5X,          ADS07430
     7' 18) DX2    =',1PE13.5,4X,' 37) ZRO    =',1PE13.5,/,5X,          ADS07440
     8' 19) EPSPEN =',1PE13.5,4X,' 38) PMLT   =',1PE13.5)               ADS07450
40    FORMAT (/5X,18HINTEGER PARAMETERS)                                ADS07460
50    FORMAT (5X,'  1) ICNDIR =',I5,4X,'  4) ITRMOP =',I5,4X,           ADS07470
     1'  6) JONED  =',I5/5X,'  2) ISCAL  =',I5,4X,'  5) ITRMST =',I5,4X,ADS07480
     2'  7) JTMAX  =',I5/5X,'  3) ITMAX  =',I5)                         ADS07490
60    FORMAT (/,5X,'ARRAY STORAGE REQUIREMENTS',/,12X,                  ADS07500
     2'DIMENSIONED  REQUIRED',/,5X,'ARRAY     SIZE        SIZE',/,      ADS07510
     37X,'WK',4X,I6,6X,I6/6X,'IWK',4X,I6,6X,I6)                         ADS07520
      END                                                               ADS07530
C --- SUBROUTINE ADS404 (NDV,H)                                         ADS07540
C     DIMENSION H(NDV*(NDV+1)/2)                                        ADS07550
C     ------------------------------------------------------------------ADS07560
C --- ROUTINE TO PRINT HESSIAN.                                         ADS07570
C     ------------------------------------------------------------------ADS07580
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                NOV. 1982    ADS07590
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS07600
C                                                                       ADS07610
      SUBROUTINE ADS404 (NDV,H)                                         ADS07620
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION H(1)                                                    ADS07630
      WRITE (6,20)                                                      ADS07640
      NA=0                                                              ADS07650
      DO 10 I=1,NDV                                                     ADS07660
      NA=NA+I-1                                                         ADS07670
      I2=MIN0(5,I)                                                      ADS07680
      WRITE (6,30) I,(H(NA+J),J=1,I2)                                   ADS07690
      IF (I.GT.5) WRITE (6,40) (H(NA+J),J=6,I)                          ADS07700
10    CONTINUE                                                          ADS07710
      RETURN                                                            ADS07720
20    FORMAT (/,5X,'HESSIAN MATRIX (LOWER HALF)')                       ADS07730
30    FORMAT (I9,')',5(1PE14.5))                                        ADS07740
40    FORMAT (10X,1PE14.5,1PE14.5,1PE14.5,1PE14.5,1PE14.5)              ADS07750
      END                                                               ADS07760
C     SUBROUTINE ADS405(I1,B,NN)                                        ADS07770
C     DIMENSION B(NN)                                                   ADS07780
C     ------------------------------------------------------------------ADS07790
C --- ROUTINE TO PRINT ONE LINE OF AN ARRAY.                            ADS07800
C     ------------------------------------------------------------------ADS07810
C     BY G. N. VANDERPLAATS                               SEPT. 1986    ADS07820
C     ENGINEERING DESIGN OPTIMIZATION, INC., SANTA BARBARA, CALIF.      ADS07830
C --- INPUT                                                             ADS07840
C     I1 - STARTING NUMBER TO BE PRINTED.                               ADS07850
C     B - ARRAY OF NUMBERS.                                             ADS07860
C     NN - LENGTH OF B.                                                 ADS07870
C --- OUTPUT                                                            ADS07880
C     PRINTED ONLY.                                                     ADS07890
      SUBROUTINE ADS405(I1,B,NN)                                        ADS07900
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION B(NN)                                                   ADS07910
      WRITE(6,10)I1,B                                                   ADS07920
10    FORMAT(I9,')',5(1PE14.5))                                         ADS07930
      RETURN                                                            ADS07940
      END                                                               ADS07950
C --- SUBROUTINE ADS420 (NDV,NCON,X,VLB,VUB,G,OBJ,IFCALL,IGCALL,IGRAD)  ADS07960
C     DIMENSION X(NDV),VLB(NDV),VUB(NDV),G(NCON)                        ADS07970
C     ------------------------------------------------------------------ADS07980
C --- ROUTINE TO PRINT THE FINAL RESULTS.                               ADS07990
C     ------------------------------------------------------------------ADS08000
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                MARCH 1983   ADS08010
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS08020
C --- INPUT                                                             ADS08030
C      NDV    - NUMBER OF DESIGN VARIABLES.                             ADS08040
C      NCON   - NUMBER OF CONSTRAINTS.                                  ADS08050
C      X      - FINAL DESIGN VARIABLE VALUES.                           ADS08060
C      VLB    - LOWER BOUND ON DESIGN VARIABLE.                         ADS08070
C      VUB    - UPPER BOUND ON DESIGN VARIABLE.                         ADS08080
C      G      - FINAL CONSTRAINT VALUES.                                ADS08090
C      OBJ    - FINAL OBJECTIVE VALUE.                                  ADS08100
C      IFCALL - NUMBER OF FUNCTION EVALUATION CALL.                     ADS08110
C      IGCALL - NUMBER OF GRADIENT EVALUATION CALL.                     ADS08120
C      IGRAD  - GRADIENT CALCULATION CONTROL.                           ADS08130
C                                                                       ADS08140
      SUBROUTINE ADS420 (NDV,NCON,X,VLB,VUB,G,OBJ,IFCALL,IGCALL,IGRAD)  ADS08150
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),VLB(1),VUB(1),G(1)                                 ADS08160

	GOTO 35       !THIS LINE ADDED BY SONGSAK NOV2006 (NO PRINTING)	

      WRITE (6,40)                                                      ADS08170
      WRITE (6,50) OBJ                                                  ADS08180
      WRITE (6,60)                                                      ADS08190
      DO 10 I=1,NDV                                                     ADS08200
10    WRITE (6,70) I,VLB(I),X(I),VUB(I)                                 ADS08210
      IF (NCON.EQ.0) GO TO 30                                           ADS08220
      N=(NCON-1)/5+1                                                    ADS08230
      WRITE (6,80)                                                      ADS08240
      DO 20 I=1,N                                                       ADS08250
      I1=5*(I-1)+1                                                      ADS08260
      I2=5*I                                                            ADS08270
      IF (I2.GT.NCON) I2=NCON                                           ADS08280
20    WRITE (6,90) I1,(G(II),II=I1,I2)                                  ADS08290
30    CONTINUE                                                          ADS08300
      WRITE (6,100) IFCALL                                              ADS08310
      IF (IGRAD.GT.0) WRITE (6,110) IGCALL                              ADS08320

35	CONTINUE      !THIS LINE ADDED BY SONGSAK NOV2006 (NO PRINTING)	

      RETURN                                                            ADS08330
C     ------------------------------------------------------------------ADS08340
C     FORMATS.                                                          ADS08350
C     ------------------------------------------------------------------ADS08360
40    FORMAT ('1',//,23X,22('-'),/,24X,'OPTIMIZATION RESULTS',/,23X,    ADS08370
     122('-'))                                                          ADS08380
50    FORMAT (//,5X,'OBJECTIVE FUNCTION VALUE',1PE15.5)                 ADS08390
60    FORMAT (//,5X,'DESIGN VARIABLES',//,18X,'LOWER',23X,'UPPER',/,5X, ADS08400
     1'VARIABLE',5X,'BOUND',9X,'VALUE',9X,'BOUND')                      ADS08410
70    FORMAT (5X,I5,2X,3(1PE14.5))                                      ADS08420
80    FORMAT (//,5X,'DESIGN CONSTRAINTS',/)                             ADS08430
90    FORMAT (5X,I5,')',5(1PE12.4))                                     ADS08440
100   FORMAT (/,5X,'FUNCTION EVALUATIONS =',I6)                         ADS08450
110   FORMAT (/,5X,'GRADIENT EVALUATIONS =',I6)                         ADS08460
      END                                                               ADS08470
C --- SUBROUTINE ADS501(N,X,Y,PROD)                                     ADS08480
C     DIMENSION X(N),Y(N)                                               ADS08490
C     ------------------------------------------------------------------ADS08500
C --- ROUTINE TO EVALUATE THE VECTOR PRODUCT OF X AND Y.                ADS08510
C     ------------------------------------------------------------------ADS08520
C     BY G. N. VANDER PLAATS                            NOV. 1980       ADS08530
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS08540
C                                                                       ADS08550
C --- INPUT                                                             ADS08560
C     N - NUMBER OF PARAMETERS IN X AND Y (DIMENSION OF ACTUAL          ADS08570
C         ARRAY STORED IN X AND Y)                                      ADS08580
C     X - FIRST INPUT VECTOR.  MINIMUM DIMENSION IS N.                  ADS08590
C     Y - SECOND INPUT VECTOR.  MINIMUM DIMENSION IS N.                 ADS08600
C                                                                       ADS08610
C --- OUTPUT                                                            ADS08620
C     PROD - VECTOR PRODUCT, X TIMES Y.                                 ADS08630
C                                                                       ADS08640
      SUBROUTINE ADS501 (N,X,Y,PROD)                                    ADS08650
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),Y(1)                                               ADS08660
      PROD=0.                                                           ADS08670
      DO 10 I=1,N                                                       ADS08680
10    PROD=PROD+X(I)*Y(I)                                               ADS08690
      RETURN                                                            ADS08700
      END                                                               ADS08710
C --- SUBROUTINE ADS502(N,G,IDG,CTMIN,CTLMIN,GMCT,IG)                   ADS08720
C     DIMENSION G(N),IDG(N)                                             ADS08730
C     ------------------------------------------------------------------ADS08740
C --- ROUTINE TO FIND MAXIMUM CONSTRAINT TOLORANCE RELATIVE TO CTMIN    ADS08750
C     ------------------------------------------------------------------ADS08760
C     BY G. N. VANDERPLAATS                           NOV. 1986         ADS08770
C     ENGINEERING DESIGN OPTIMIZATION, INC., SANTA BARBARA, CA          ADS08780
C                                                                       ADS08790
C --- INPUT                                                             ADS08800
C     N - NUMBER OF PARAMETERS CONTAINED IN ARRAY G AND IDG.            ADS08810
C     G - ARRAY OF CONSTRAINTS.  MINIMUM DIMENSION IS N.                ADS08820
C     IDG - CONSTRAINT TYPE IDENTIFIER.                                 ADS08830
C     CTMIN - NONLINEAR CONSTRAINT TOLORANCE.                           ADS08840
C     CTLMIN - LINEAR CONSTRAINT TOLORANCE.                             ADS08850
C                                                                       ADS08860
C --- OUTPUT                                                            ADS08870
C     GMCT - VALUE OF MAXIMUM COMPONENT OF G MINUS CTMIN OR CTLMIN.     ADS08880
C     IG - LOCATION IN G OF MAXIMUM COMPONENT.                          ADS08890
C                                                                       ADS08900
      SUBROUTINE ADS502 (N,G,IDG,CTMIN,CTLMIN,GMCT,IG)                  ADS08910
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION G(1),IDG(1)                                             ADS08920
      IG=1                                                              ADS08930
      GMCT=G(1)-CTMIN                                                   ADS08940
      IF(IABS(IDG(1)).EQ.2) GMCT=G(1)-CTLMIN                            ADS08950
      IF (N.LE.1) RETURN                                                ADS08960
      DO 10 I=2,N                                                       ADS08970
      GG=G(I)-CTMIN                                                     ADS08980
      IF(IABS(IDG(I)).EQ.2) GG=G(I)-CTLMIN                              ADS08990
      IF (GG.LT.GMCT) GO TO 10                                          ADS09000
      GMCT=GG                                                           ADS09010
      IG=I                                                              ADS09020
10    CONTINUE                                                          ADS09030
      RETURN                                                            ADS09040
      END                                                               ADS09050
C --- SUBROUTINE ADS503(N,X,Y)                                          ADS09060
C     DIMENSION X(N),Y(N)                                               ADS09070
C     ------------------------------------------------------------------ADS09080
C --- ROUTINE TO EQUATE VECTOR Y TO VECTOR X                            ADS09090
C     ------------------------------------------------------------------ADS09100
C     BY G. N. VANDERPLAATS                        DEC. 1980            ADS09110
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS09120
C                                                                       ADS09130
C --- INPUT                                                             ADS09140
C     N - NUMBER OF PARAMETERS CONTAINED IN X.                          ADS09150
C     X - VECTOR OF PARAMETERS.  MINIMUM DIMENSION IS N.                ADS09160
C                                                                       ADS09170
C --- OUTPUT                                                            ADS09180
C     Y - VECTOR OF PARAMETERS, Y=X.  MINIMUM DIMENSION IS N.           ADS09190
C                                                                       ADS09200
      SUBROUTINE ADS503 (N,X,Y)                                         ADS09210
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),Y(1)                                               ADS09220
      DO 10 I=1,N                                                       ADS09230
10    Y(I)=X(I)                                                         ADS09240
      RETURN                                                            ADS09250
      END                                                               ADS09260
C --- SUBROUTINE ADS504 (N,H,Y,S)                                       ADS09270
C     DIMENSION H(1),Y(1),S(1)                                          ADS09280
C     ------------------------------------------------------------------ADS09290
C --- ROUTINE TO CALCULATE H(N,N)*Y(N)=S(N)                             ADS09300
C     ------------------------------------------------------------------ADS09310
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO.                 NIV. 1982  ADS09320
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS09330
C     H IS STORED AS .....                                              ADS09340
C                                                                       ADS09350
C        I 1 2 4 7...I                                                  ADS09360
C        I X 3 5 8...I                                                  ADS09370
C        I X X 6 9...I                                                  ADS09380
C        I X X X10...I = ( 1 2 3 4 5 6 7 8 9 10 .......)                ADS09390
C        I . . . ....I                                                  ADS09400
C        I . . . ....I                                                  ADS09410
C        I . . . ....I                                                  ADS09420
C                                                                       ADS09430
      SUBROUTINE ADS504 (N,H,Y,S)                                       ADS09440
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION H(1),Y(1),S(1)                                          ADS09450
      DO 10 I=1,N                                                       ADS09460
      S(I)=0.                                                           ADS09470
10    CONTINUE                                                          ADS09480
      NA=0                                                              ADS09490
      DO 30 I=1,N                                                       ADS09500
      DO 20 J=1,I                                                       ADS09510
      NA=NA+1                                                           ADS09520
      S(I)=S(I)+H(NA)*Y(J)                                              ADS09530
      IF (J.EQ.I) GO TO 30                                              ADS09540
      S(J)=S(J)+H(NA)*Y(I)                                              ADS09550
20    CONTINUE                                                          ADS09560
30    CONTINUE                                                          ADS09570
      RETURN                                                            ADS09580
      END                                                               ADS09590
C --- SUBROUTINE ADS506(N,A)                                            ADS09600
C     DIMENSION A(N*(N+1)/2)                                            ADS09610
C     ------------------------------------------------------------------ADS09620
C --- ROUTINE TO INITIALIZE A SYMMETRIC MATRIX OF DIMENSION N TO THE    ADS09630
C     IDENTITY MATRIX, WHERE THE MATRIX IS STORED COLUMN AFTER COLUMN   ADS09640
C     IN ARRAY A.                                                       ADS09650
C     ------------------------------------------------------------------ADS09660
C     BY G. N. VANDERPLAATS                           JAN 1981          ADS09670
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS09680
C                                                                       ADS09690
C --- INPUT                                                             ADS09700
C     N - DIMENSION OF ARRAY STORED IN A.                               ADS09710
C                                                                       ADS09720
C --- OUTPUT                                                            ADS09730
C     A - IDENTITY MATRIX, I.                                           ADS09740
C                                                                       ADS09750
      SUBROUTINE ADS506 (N,A)                                           ADS09760
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(1)                                                    ADS09770
C     INITIALIZE TO ZERO.                                               ADS09780
      NN=N*(N+1)/2                                                      ADS09790
      DO 10 I=1,NN                                                      ADS09800
10    A(I)=0.                                                           ADS09810
C     SET DIAGONAL TO 1.                                                ADS09820
      II=1                                                              ADS09830
      DO 20 I=1,N                                                       ADS09840
      A(II)=1.                                                          ADS09850
20    II=II+I+1                                                         ADS09860
      RETURN                                                            ADS09870
      END                                                               ADS09880
C --- SUBROUTINE ADS507 (N,A,ZRO,SCAL,IER)                              ADS09890
C     DIMENSION A(N)                                                    ADS09900
C     ------------------------------------------------------------------ADS09910
C --- ROUTINE TO NORMALIZE VECTOR A BY DIVIDING BY MAX ABS(A(I)).       ADS09920
C     BY G. N. VANDERPLAATS                           JAN 1981          ADS09930
C     ------------------------------------------------------------------ADS09940
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS09950
C                                                                       ADS09960
C --- INPUT                                                             ADS09970
C     N - NUMBER OF PARAMETERS CONTAINED IN A.                          ADS09980
C     A - VECTOR OF PARAMETERS.  MINIMUM DIMENSION IS N.                ADS09990
C     ZRO - TOLERANCE.  IF ABS(A(I)).LE.ZRO, IT IS ASSUMED TO BE ZERO.  ADS10000
C         THIS IS A MEASURE OF MACHINE ACCURACY.                        ADS10010
C                                                                       ADS10020
C --- OUTPUT                                                            ADS10030
C     A - NORMALIZED ARRAY.                                             ADS10040
C     SCAL - NORMALIZATION FACTOR.                                      ADS10050
C     IER - ERROR CODE.  IF IER.GT.0, ALL COMPONENTS OF A ARE LESS      ADS10060
C         IN MAGNITUDE THAN ZRO.                                        ADS10070
C                                                                       ADS10080
      SUBROUTINE ADS507 (N,A,ZRO,SCAL,IER)                              ADS10090
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(1)                                                    ADS10100
C     BEGIN NORMALIZATION.                                              ADS10110
      AMAX=0.                                                           ADS10120
      DO 10 I=1,N                                                       ADS10130
      AA=ABS(A(I))                                                      ADS10140
      IF (AA.GT.AMAX) AMAX=AA                                           ADS10150
10    CONTINUE                                                          ADS10160
C     CHECK FOR ZERO ARRAY.                                             ADS10170
      IER=1                                                             ADS10180
      SCAL=ZRO                                                          ADS10190
      IF (AMAX.LE.ZRO) RETURN                                           ADS10200
C     A HAS NON-ZERO COMPONENTS.                                        ADS10210
      IER=0                                                             ADS10220
C     SCALE FACTOR.                                                     ADS10230
      SCAL=AMAX                                                         ADS10240
C     PERFORM SCALING.                                                  ADS10250
C     GUARD AGAINST ZERO SCALE FACTOR (ONLY HAPPENS IF ZRO=0).          ADS10260
      SS=1./SCAL                                                        ADS10270
      DO 20 I=1,N                                                       ADS10280
20    A(I)=SS*A(I)                                                      ADS10290
      RETURN                                                            ADS10300
      END                                                               ADS10310
C --- SUBROUTINE ADS508 (NDV,NMG,GMULT,A,NRA,DF)                        ADS10320
C     DIMENSION A(NRA,NCOLA),DF(NDV)                                    ADS10330
C     ------------------------------------------------------------------ADS10340
C --- ROUTINE TO CALCULATE THE GRADIENTS OF PSEUDO-OBJECTIVE FOR SUMT.  ADS10350
C     ------------------------------------------------------------------ADS10360
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                JULY. 1983   ADS10370
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS10380
C --- INPUT                                                             ADS10390
C      NDV   - NUMBER OF DESIGN VARIABLES.                              ADS10400
C      NMG   - NUMBER OF COLUMN IN ARRAY A.                             ADS10410
C      GMULT - FACTOR OF GRADIENT OF SPECIFIED CONSTRAINT.              ADS10420
C      A     - ARRAY CONTAINING GRADIENTS OF ACTIVE CONSTRAINTS.        ADS10430
C     NRA    - NUMBER OF ROWS IN ARRAY A.                               ADS10440
C     DF     - GRADIENT OF PRIMAL OBJECTIVE.                            ADS10450
C --- OUTPUT                                                            ADS10460
C      DF    - GRADIENT OF PSEUDO-OBJECTIVE.                            ADS10470
C                                                                       ADS10480
      SUBROUTINE ADS508 (NDV,NMG,GMULT,A,NRA,DF)                        ADS10490
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(NRA,1),DF(1)                                          ADS10500
C                                                                       ADS10510
      DO 10 I=1,NDV                                                     ADS10520
      DF(I)=DF(I)+GMULT*A(I,NMG)                                        ADS10530
10    CONTINUE                                                          ADS10540
      RETURN                                                            ADS10550
      END                                                               ADS10560
C --- SUBROUTINE ADS509 (NDV,X,SVLB,SVUB,ZRO,ISIDC,NACS)                ADS10570
C     DIMENSION X(NDV),SVLB(NDV),SVUB(NDV),ISIDC(NDV)                   ADS10580
C     ------------------------------------------------------------------ADS10590
C --- ROUTINE TO IDENTIFY ACTIVE SIDE CONSTRAINTS.                      ADS10600
C     ------------------------------------------------------------------ADS10610
C     BY G. N. VANDERPLAATS AND H. SUGIMOTO                 JULY 1983   ADS10620
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS10630
C --- INPUT                                                             ADS10640
C      NDV  - NUMBER OF DESIGN VARIABLES.                               ADS10650
C      X    - DESIGN VARIABLES.                                         ADS10660
C      SVLB,SVUB - LOWER AND UPPER BOUNDS ON DESIGN VARIABLE.           ADS10670
C      ZRO  - NUMERICAL ZERO VALUE.                                     ADS10680
C --- OUTPUT                                                            ADS10690
C      ISIDC - = 0;SIDE CONSTRAINTS ARE NOT ACTIVE.                     ADS10700
C              =-1;LOWER SIDE CONSTRAINT IS ACTIVE.                     ADS10710
C              = 1;UPPER SIDE CONSTRAINT IS ACTIVE.                     ADS10720
C              = 2;BOTH SIDE CONSTRAINTS ARE ACTIVE.                    ADS10730
C      NACS  - NUMBER OF ACTIVE SIDE CONSTRAINTS.                       ADS10740
C                                                                       ADS10750
      SUBROUTINE ADS509 (NDV,X,SVLB,SVUB,ZRO,ISIDC,NACS)                ADS10760
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SVLB(1),SVUB(1),ISIDC(1)                           ADS10770
C                                                                       ADS10780
      NACS=0                                                            ADS10790
      DO 10 I=1,NDV                                                     ADS10800
      DXL=ABS(X(I)-SVLB(I))                                             ADS10810
      DXU=ABS(SVUB(I)-X(I))                                             ADS10820
      ISIDC(I)=0                                                        ADS10830
      IF (DXL.LE.ZRO) ISIDC(I)=-1                                       ADS10840
      IF (DXU.LE.ZRO) ISIDC(I)=1                                        ADS10850
      IF (DXL.LE.ZRO.AND.DXU.LE.ZRO) ISIDC(I)=2                         ADS10860
      IF (ISIDC(I).EQ.0) GO TO 10                                       ADS10870
      NACS=NACS+1                                                       ADS10880
10    CONTINUE                                                          ADS10890
      RETURN                                                            ADS10900
      END                                                               ADS10910
C --- SUBROUTINE ADS510 (ICALC,NDV,X,SX,SCAL,SCLMIN)                    ADS10920
C     DIMENSION X(NDV),SX(NDV),SCAL(NDV)                                ADS10930
C     ------------------------------------------------------------------ADS10940
C --- ROUTINE TO SET UP SCALING VECTOR, SCAL, SCALE THE X-VECTOR, OR    ADS10950
C     UN-SCALE THE X-VECTOR.                                            ADS10960
C     ------------------------------------------------------------------ADS10970
C     BY G. N. VANDERPLAATS                              MARCH, 1982.   ADS10980
C     NAVAL POSTGRADUATE SCHOOL, MONTEREY, CALIF.                       ADS10990
C --- INPUT                                                             ADS11000
C     ICALC - CALCULATION CONTROL.                                      ADS11010
C         0 - SET SCALING VECTOR TO UNITY (NO SCALING).                 ADS11020
C         1 - DEFINE SCALING VECTOR.                                    ADS11030
C         2 - SCALE THE X-VECTOR.                                       ADS11040
C         3 - UN-SCALE THE X-VECTOR.                                    ADS11050
C     NDV - NUMBER OF DESIGN VARIABLES.                                 ADS11060
C     X - VECTOR OF DESIGN VARIABLES, UNSCALED.                         ADS11070
C         MINIMUM DIMENSION IS NDV.                                     ADS11080
C     SX - VECTOR OF SCALED DESIGN VARIABLES.  MINIMUM DIMENSION IS NDV.ADS11090
C     SCAL - SCALING VECTOR.  SX(I)=SCAL(I)*X(I), I=1,NDV.              ADS11100
C         MINIMUM DIMENSION IS NDV.                                     ADS11110
C     SCLMIN - MINIMUM VALUE FOR SCALING PARAMETER.  USED WHEN X(I) IS  ADS11120
C         NEAR ZERO.                                                    ADS11130
C --- OUTPUT                                                            ADS11140
C     X, SX OR SCAL, DEPENDING ON THE VALUE OF ICALC.                   ADS11150
C                                                                       ADS11160
      SUBROUTINE ADS510 (ICALC,NDV,X,SX,SCAL,SCLMIN)                    ADS11170
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      DIMENSION X(1),SX(1),SCAL(1)                                      ADS11180
C                                                                       ADS11190
      ICALC1=ICALC+1                                                    ADS11200
      GO TO (10,30,50,70), ICALC1                                       ADS11210
10    CONTINUE                                                          ADS11220
C     SET SCALING VECTOR TO UNITY.                                      ADS11230
      DO 20 I=1,NDV                                                     ADS11240
20    SCAL(I)=1.                                                        ADS11250
      RETURN                                                            ADS11260
30    CONTINUE                                                          ADS11270
C     DEFINE THE SCALING VECTOR, SCAL, AS THE MAGNITUDE OF THE          ADS11280
C     COMPONENTS OF THE X-VECTOR OR SCLMIN, WHICHEVER IS LARGER.        ADS11290
      DO 40 I=1,NDV                                                     ADS11300
      XI=ABS(X(I))                                                      ADS11310
      IF (XI.LT.SCLMIN) XI=SCLMIN                                       ADS11320
40    SCAL(I)=1./XI                                                     ADS11330
      RETURN                                                            ADS11340
50    CONTINUE                                                          ADS11350
C     SCALE THE X-VECTOR.                                               ADS11360
      DO 60 I=1,NDV                                                     ADS11370
60    SX(I)=X(I)*SCAL(I)                                                ADS11380
      RETURN                                                            ADS11390
70    CONTINUE                                                          ADS11400
C     UN-SCALE THE X-VECTOR.                                            ADS11410
      DO 80 I=1,NDV                                                     ADS11420
80    X(I)=SX(I)/SCAL(I)                                                ADS11430
      RETURN                                                            ADS11440
      END                                                               ADS11450

