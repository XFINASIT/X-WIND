C	=======================================================================
C	===================== DAMAGE MECHANIC OF CONCRETE =====================
C	========== Design for Standard 8 node Solid Element XSolid81 ==========
C	=======================================================================
	SUBROUTINE	DAMAGE(STNCR,STSPS,STIFG,PROPM,WASOL)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
C	COMPUTE STRESS AND DAMAGE CONSTITUTIVE MATRIX OF PURE CONCRETE
C	PRODUCED BY SONGSAK MAR2006
C	=======================================================================
C	=======================================================================
	DIMENSION STSPS(6),STNPS(6),STNCR(6),STIFG(6,6),
	1		  DIREC(9),DAMPR(3),STNPT(3),DIRPT(3,3),TRMGL(6,6),
	2		  GRTST(3),PLSTN(3),ESTAR(3),PSTAR(3,3),GSTAR(3,3),
	3		  SGASH(6,6),STSCR(6),PROPM(*),WASOL(32)
C	=======================================================================
C	THE MEANING OF PARAMETERS
C	=======================================================================
C	YOUNG = YOUNG MODULUS	
C	POISN = POISON RATIO
C	GMODU = MODULUS OF RIGIDITY
C	GFACT = FRACTURE ENERGY PER UNIT LENGTH OF CRACK WIDTH
C	VLCH  = CHARACTERISTIC LENGTH
C	UNIAX = MAXIMUM COMPRESSIVE STRENGTH
C	TENST = MAXIMUM TENSILE STRENGTH
C	EMAXS = STRAIN AT MAXIMUM TENSILE STRESS
C	IDTCS = FLAG --->  0 = CONCRETE , 1 = STEEL
C	=======================================================================
C	--------------------------------------------------------
C	READ THE INPUT DATA 
C	--------------------------------------------------------
	YOUNG = PROPM(1)
	POISN = PROPM(2)
	GMODU = PROPM(3)
	GFACT = PROPM(4)
	VLCH  = PROPM(5)
	UNIAX = PROPM(6)
	TENST = PROPM(7)
	EMAXS = PROPM(8)
	IDTCS = PROPM(10)

C	WRITE(*,*) (PROPM(I),I=1,10)

C	--------------------------------------------------------
C	READ THE PROBLEM DATA INTO WORKING ARRAY
C	--------------------------------------------------------
	STSPS(1) = WASOL(1)
	STSPS(2) = WASOL(2)
	STSPS(3) = WASOL(3)
	STSPS(4) = WASOL(4)
	STSPS(5) = WASOL(5)
	STSPS(6) = WASOL(6)

	STNPS(1) = WASOL(7)
	STNPS(2) = WASOL(8)
	STNPS(3) = WASOL(9)
	STNPS(4) = WASOL(10)
	STNPS(5) = WASOL(11)
	STNPS(6) = WASOL(12)

	DIREC(1) = WASOL(13)
	DIREC(2) = WASOL(14)
	DIREC(3) = WASOL(15)
	DIREC(4) = WASOL(16)
	DIREC(5) = WASOL(17)
	DIREC(6) = WASOL(18)
	DIREC(7) = WASOL(19)
	DIREC(8) = WASOL(20)
	DIREC(9) = WASOL(21)

	DAMPR(1) = WASOL(22)
	DAMPR(2) = WASOL(23)
	DAMPR(3) = WASOL(24)

	GRTST(1) = WASOL(25)
	GRTST(2) = WASOL(26)
	GRTST(3) = WASOL(27)

	PLSTN(1) = WASOL(28)
	PLSTN(2) = WASOL(29)
	PLSTN(3) = WASOL(30)

	MSTAT    = WASOL(31)
	IPEL     = WASOL(32)

	IF(MSTAT.EQ.0.AND.IDTCS.EQ.0)     MSTAT = 1
	IF(MSTAT.EQ.0.AND.IDTCS.EQ.1)     MSTAT = 6
	IF(IPEL.EQ.0) IPEL = 1


C	-------------------------------------------------------
C	COMPUTE FRACTURE RESISTANT PARAMETER
C	-------------------------------------------------------
	EVALU = 2.718281828
	CON1  = (2.0*YOUNG*GFACT/VLCH/TENST/TENST) - 1.0
	AVALU = 3.0/EMAXS/CON1
	IF (AVALU.LT.0.0) AVALU = 0.1E-20

C	-------------------------------------------------------
C     COMPUTE THE CRITICAL STRAIN	
C	-------------------------------------------------------
	ECRIT = EMAXS + 4.6/AVALU

C	-------------------------------------------------------
C	MATERIAL IDENTIFICATION NUMBER
C	-------------------------------------------------------
C	NSTAT = 1 NO CRACKS OCCUR	
C	NSTAT = 2 CRACK IN ONE DIRECTION
C	NSTAT = 3 CRACK IN TWO DIRECTION
C	NSTAT = 4 CRACK IN THREE DIRECTION	
C	NSTAT = 5 CRUSHING OF CONCRETE	
C	-------------------------------------------------------	

	NSTAT = MSTAT
	
C	-------------------------------------------------------	
C     THIS LINE IS FOR CONCRETE THAT ALREADY CRACK
C	-------------------------------------------------------	
	IF(NSTAT.GT.1.AND.IDTCS.EQ.0) GO TO 500

	IF(IDTCS.EQ.1) GO TO 1000

C	-------------------------------------------------------	
C     DETERMINE THE PRINCIPAL STRESSES AND DIRECTION
C	-------------------------------------------------------		
	CALL	PRISTD(STNCR,STNPT,DIRPT)


C	-------------------------------------------------------	
C     TEST FOR NEW CRACKS
C	-------------------------------------------------------		
	NGASH = 1
	DO 5 I = 1,3
	IF(STNPT(I).GT.EMAXS) NGASH = NGASH + 1 !CRACK OCCURS
5	CONTINUE
	NSTAT = NGASH

C	-------------------------------------------------------	
C     IF THERE IS NO CRACK OCCUR --->ELASTIC BEHAVIOR
C	-------------------------------------------------------	
	IF(NSTAT.EQ.1) GO TO 1000

C	-------------------------------------------------------	
C     STORE THE PRINCIPAL DIRECTION 
C	-------------------------------------------------------	
	IGASH = 0
	DO 10 I = 1,3
	DO 10 J = 1,3
	IGASH = IGASH + 1
10	DIREC(IGASH) = DIRPT(I,J)

C	-------------------------------------------------------	
C     COMPUTE DAMAGE PARAMETERS
C	-------------------------------------------------------	
	DO 20 I = 1,3

	IF(STNPT(I).GT.EMAXS) THEN !!!

	CALL	DMPAR(AVALU,EMAXS,STNPT(I),FKO)
	DAMPR(I) = FKO

C	-------------------------------------------------------
C	STORE MAXIMUM STRAIN AND COMPUTE PLASTIC STRAIN
C	-------------------------------------------------------
	DELTA = 0.2
	GRTST(I) = STNPT(I)    
	PLSTN(I) = DELTA*STNPT(I)   
C	-------------------------------------------------------

	ELSE !!!

	DAMPR(I) = 0.0
	
	ENDIF !!!

20	CONTINUE

	
	GO TO 600

	
500	CONTINUE


C	--------------------------------------------------------
C	THIS PART IS FOR CRACKED CONCRETE
C	--------------------------------------------------------

C	--------------------------------------------------------
C	READ THE PRINCIPAL DIRECTION
C	--------------------------------------------------------	
	IGASH = 0
	DO 80 I = 1,3
	DO 80 J = 1,3
	IGASH = IGASH + 1
80	DIRPT(I,J) = DIREC(IGASH)

C	--------------------------------------------------------
C	TRANSFORM STRAIN INTO CRACK DIRECTION
C	--------------------------------------------------------
	CALL	TRNGLX(DIRPT,TRMGL)
	DO 90 I = 1,3
	STNPT(I) = 0.0
	DO 90 J = 1,6
90	STNPT(I) = STNPT(I) + TRMGL(I,J)*STNCR(J)

	DO 100 I = 1,3

	IF(GRTST(I).LT.EMAXS) GO TO 95
												     		   	
      IF(STNPT(I).GT.GRTST(I)) THEN  !    <-------------------! 
												     		!  
	CALL	DMPAR(AVALU,EMAXS,STNPT(I),FKO)                 ! 
	DAMPR(I) = FKO                                          ! 
												     		!  
C	------------------------------------------------------- ! 
C	STORE MAXIMUM STRAIN AND COMPUTE PLASTIC STRAIN         ! 
C	------------------------------------------------------- ! 
	DELTA = 0.2                                             ! 
	GRTST(I) = STNPT(I)                                     ! 
	PLSTN(I) = DELTA*STNPT(I)                               ! 
														    ! 
	ELSE  !  <----------------------------------------------! 
	DAMPR(I) = DAMPR(I)                                     ! 
	ENDIF !  <----------------------------------------------! 
												     		   
	IF(STNPT(I).LE.PLSTN(I)) GO TO 1000  ! CRACKS ARE CLOSED  
	
	GO TO 100
												     		  
95	CONTINUE 

												     		   
	IF(STNPT(I).GT.EMAXS) THEN !<---------------------------! 
												     		! 
	CALL	DMPAR(AVALU,EMAXS,STNPT(I),FKO)                 ! 
	DAMPR(I) = FKO                                          ! 
C	------------------------------------------------------- ! 
C	STORE MAXIMUM STRAIN AND COMPUTE PLASTIC STRAIN         ! 
C	------------------------------------------------------- ! 
	DELTA = 0.2                                             ! 
	GRTST(I) = STNPT(I)                                     ! 
	PLSTN(I) = DELTA*STNPT(I)                               ! 
														    ! 
C	------------------------------------------------------- ! 
C	UPDATE THE MATERIAL NUMBER IF NEW CRACKS OCCUR          ! 
C	------------------------------------------------------- ! 
	NSTAT = NSTAT + 1                                       ! 
												     		! 
	ENDIF   !<----------------------------------------------! 
												     		  
      
100	CONTINUE

	
600	CONTINUE	
C	-------------------------------------------------------
C	CONSTRUCT THE DAMAGE CONSTITUTIVE MATRIX
C	-------------------------------------------------------
	CALL ESTRFF (YOUNG,POISN,GMODU,DAMPR,STIFG,STNPT,ECRIT)
C	-------------------------------------------------------

C	-------------------------------------------------------
C	TRANSFORM THE CONSTITUTIVE MATRIX INTO GLOBAL SYSTEM
C	-------------------------------------------------------
	CALL	TRNGLX(DIRPT,TRMGL)

	

	DO 50 I = 1,6
	DO 50 J = 1,6
	SGASH(I,J) = 0.0
	DO 50 K = 1,6
50	SGASH(I,J) = SGASH(I,J) + STIFG(I,K)*TRMGL(K,J)

	DO 60 I = 1,6
	DO 60 J = 1,6
	STIFG(I,J) = 0.0
	DO 60 K = 1,6
60	STIFG(I,J) = STIFG(I,J) + SGASH(K,J)*TRMGL(K,I)


C	-------------------------------------------------------
C	FOR CRACKED BEHAVIOR, COMPUTE DAMAGE STRESSES
C	-------------------------------------------------------	
	DO 70 I  = 1,6
	STSPS(I) = 0.0
	DO 70 J  = 1,6
70	STSPS(I) = STSPS(I) + STIFG(I,J)*STNCR(J)
C	-------------------------------------------------------


	GO TO 2000


1000	CONTINUE
C	--------------------------------------------------------
C	FOR UNCRACK BEHAVIOR COMPUTE ELASTIC CONSTITUTIVE MATRIX
C	--------------------------------------------------------
c	CALL	DMMISE(STSPS,STNPS,UNIAX,IPEL,STNCR,STSCR,STIFG,
c	1			   YOUNG,POISN,IDTCS,MSTAT)
	CALL ELTIFF_DM(YOUNG,POISN,STIFG)


	DO 150 I  = 1,6
	STSPS(I) = 0.0
	DO 150 J  = 1,6
150	STSPS(I) = STSPS(I) + STIFG(I,J)*STNCR(J)


C	--------------------------------------------------------
C	FOR UNCRACK BEHAVIOR, COMPUTE ELASTIC STRESSES
C	--------------------------------------------------------
c	DO 160 I = 1,6
c	STSPS(I) = 0.0
c160	STSPS(I) = STSPS(I) + STSCR(I)


2000	CONTINUE

	
C	--------------------------------------------------------
C	UPDATE THE MATERIAL NUMBER
C	--------------------------------------------------------
	IF(IDTCS.EQ.1) GO TO 2500

	MSTAT = NSTAT

2500	CONTINUE

C	--------------------------------------------------------
C	STORE THE PROBLEM DATA INTO WORKING ARRAY
C	--------------------------------------------------------
	WASOL(1)  = STSPS(1)
	WASOL(2)  = STSPS(2)
	WASOL(3)  = STSPS(3)
	WASOL(4)  = STSPS(4)
	WASOL(5)  = STSPS(5)
	WASOL(6)  = STSPS(6)

	WASOL(7)  = STNCR(1)
	WASOL(8)  = STNCR(2)
	WASOL(9)  = STNCR(3)
	WASOL(10) = STNCR(4)
	WASOL(11) = STNCR(5)
	WASOL(12) = STNCR(6)

	WASOL(13) = DIREC(1)
	WASOL(14) = DIREC(2)
	WASOL(15) = DIREC(3)
	WASOL(16) = DIREC(4)
	WASOL(17) = DIREC(5)
	WASOL(18) = DIREC(6)
	WASOL(19) = DIREC(7)
	WASOL(20) = DIREC(8)
	WASOL(21) = DIREC(9)

	WASOL(22) = DAMPR(1)
	WASOL(23) = DAMPR(2)
	WASOL(24) = DAMPR(3)

	WASOL(25) = GRTST(1)
	WASOL(26) = GRTST(2)
	WASOL(27) = GRTST(3)

	WASOL(28) = PLSTN(1)
	WASOL(29) = PLSTN(2)
	WASOL(30) = PLSTN(3)

	WASOL(31) = MSTAT
	WASOL(32) = IPEL



	RETURN

	END


C	=======================================================================
C	=======================================================================
	SUBROUTINE	PRISTD(STNCR,STNPT,DIRPT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	DIMENSION STNCR(6),STNPT(3),DIRPT(3,3),
	1		  TRMGL(6,6),DIRPTT(3,3)	
C	=======================================================================


C	STNPT = PRINCIPAL STRAIN
C	DIRPT = PRINCIPAL DIRECTION

	SI1 = (STNCR(1) + STNCR(2) + STNCR(3))
	SI2 = STNCR(1)*STNCR(2)      + STNCR(2)*STNCR(3) +
	1	  STNCR(3)*STNCR(1)      - 0.25*STNCR(4)*STNCR(4) -
	2	  0.25*STNCR(5)*STNCR(5) - 0.25*STNCR(6)*STNCR(6)
	SI3 = STNCR(1)*STNCR(2)*STNCR(3)      -
	1	  0.25*STNCR(1)*STNCR(6)*STNCR(6) -
	2	  0.25*STNCR(2)*STNCR(5)*STNCR(5) -
	3	  0.25*STNCR(3)*STNCR(4)*STNCR(4) +
	4	  0.25*STNCR(4)*STNCR(5)*STNCR(6)

	SM = SI1/3.0

	DSI1 = 0.0
	DSI2 = SI1*SI1/3.0 - SI2
	DSI3 = SI3 - SM*SI2 + 2.0*SM*SM*SM 


	IF(DSI2.EQ.0.0) THEN
	STNPT(1) = 0.0
	STNPT(2) = 0.0
	STNPT(3) = 0.0
	ELSE
	CON1  = -1.5*SQRT(3.0)
	CON2  = 2.0*SQRT(DSI2/3.0)
	CON3  = SI1/3.0
	CON4  = DSI2**1.5



	IF(CON1*DSI3/CON4.GE.-1.000001.OR.CON1*DSI3/CON4.LE.-0.99999) THEN
	THETA = -1.57079
	ELSE
	THETA = ASIN(CON1*DSI3/CON4)
	ENDIF
c	write(*,*) theta

	


	THETA = THETA / 3.0
	PI    = 3.141592654
	STNPT(1) = CON2*SIN(THETA + 2.0*PI/3.0) + CON3
	STNPT(2) = CON2*SIN(THETA) + CON3
	STNPT(3) = CON2*SIN(THETA + 4.0*PI/3.0) + CON3
	ENDIF
C	--------------------------------------------------------
C	COMPUTE PRINCIPAL DIRECTION
C
C		          [L1  L2  L3]
C		 DIRPT =  [M1  M2  M3]
C		          [N1  N2  N3]
C
C	--------------------------------------------------------
	
	DO 20 I = 1,3

	PL = -2.0*STNCR(5)*(STNCR(2) - STNPT(I)) + (STNCR(4)*STNCR(6))	
	PM = -2.0*STNCR(6)*(STNCR(1) - STNPT(I)) + (STNCR(5)*STNCR(4))
	PN = 4.0*(STNCR(1) - STNPT(I))*(STNCR(2) - STNPT(I)) -
	1                                 (STNCR(4)*STNCR(4))
	FAC = SQRT(PL*PL+PM*PM+PN*PN)

	IF(FAC.EQ.0.0) THEN
	DIRPTT(1,I) = 0.0
	DIRPTT(2,I) = 0.0
	DIRPTT(3,I) = 0.0
	ELSE
	DIRPTT(1,I) = PL/FAC
	DIRPTT(2,I) = PM/FAC
	DIRPTT(3,I) = PN/FAC
	ENDIF

20	CONTINUE

	DO 25 I = 1,3
	DIRPT(I,1) = DIRPTT(I,1)
	DIRPT(I,2) = DIRPTT(I,2)
25	DIRPT(I,3) = DIRPTT(I,3)

	RETURN
	END


C     ============================================================================	
	SUBROUTINE ESTRFF (YOUNG,POISN,GMODU,DAMPR,ESTR,STNPT,ECRIT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	============================================================================
	DIMENSION DAMPR(3),ESTR(6,6),EI(3),GIJSR(3,3),VIJ(3,3),STNPT(3)
C	============================================================================
	
	DO 303 I = 1,3
	EI(I) = ((1.0-DAMPR(I))**2)*YOUNG
303	CONTINUE

	DO 403 I = 1,3
	DO 403 J = 1,3
	VIJ(I,J) = ((1.0-DAMPR(I))/(1.0-DAMPR(J)))*POISN
403	CONTINUE
	
	DELTA = 1.0-(VIJ(1,2)*VIJ(2,1))-(VIJ(2,3)*VIJ(3,2))-
	1            (VIJ(1,3)*VIJ(3,1))-2.0*(VIJ(2,1)*VIJ(3,2)*VIJ(1,3))
	
	DO 505 I   = 1,3
	DO 505 J   = 1,3
	GIJSR(I,J) = 2.0*(1.0-DAMPR(I))**2.0*(1.0-DAMPR(J))**2.0*GMODU/
	1                 ((1.0-DAMPR(I))**2.0+(1.0-DAMPR(J))**2.0)
505	CONTINUE

	DO 550  I = 1,6
	DO 550  J = 1,6
550	ESTR(I,J) = 0.0

	DO 610  I = 1,3
	DO 610  J = 1,3
	DO 610  K = 1,3
	IF(J.NE.K) THEN
	ESTR(I,I) = EI(I)*(1.0 - VIJ(J,K)*VIJ(K,J))/DELTA
	ENDIF
610	CONTINUE

	DO 650  I = 1,3
	DO 650  J = 1,3
	DO 650  K = 1,3
	IF(I.NE.J.AND.K.NE.I.AND.J.NE.K) THEN
	ESTR(I,J) = ESTR(I,J) + EI(I)*(VIJ(J,I) + VIJ(K,I)*VIJ(J,K))/DELTA
	ENDIF
650	CONTINUE

	ESTR(4,4) = GIJSR(1,2)
	ESTR(5,5) = GIJSR(2,3)
	ESTR(6,6) = GIJSR(3,1)

C	--------------------------------------------------------
C	IF THE CRACKED STRAIN IS VERY LARGE, STIFFNESS = 0.0
C	--------------------------------------------------------
	DO 700 I = 1,3
	DO 700 J = 1,3
	IF(STNPT(I).LT.ECRIT) GO TO 700
	ESTR(I,J) = 0.0
	ESTR(J,I) = 0.0
700	CONTINUE


	RETURN
	END 


C	=======================================================================
C	=======================================================================
	SUBROUTINE	TRNGLX(DIRPT,TRMGL)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	DIMENSION DIRPT(3,3),TRMGL(6,6)
C	=======================================================================

C						    EPSILON_L = TRMGL * EPSILON_G
C							SIGMA_G   = TRANSPOSE(TRMGL) * SIGMA_L
C
C				[L1L1   M1M1   N1N1   L1M1       L1N1       M1N1     ]
C				[L2L2   M2M2   N2N2   L2M2       L2N2       M2N2     ]
C				[L3L3   M3M3   N3N3   L3M3       L3N3       M3N3     ]
C				[2L1L2  2M1M2  2N1N2  L1M2+M1L2  L1N2+N1L2  M1N2+N1M2]
C				[2L1L3  2M1M3  2N1N3  L1M3+M1L3  L1N3+N1L3  M1N3+N1M3]
C				[2L2L3  2M2M3  2N2N3  L2M3+M2L3  L2N3+N2L3  M2N3+N2M3]

	DO 10 I = 1,3
	DO 10 J = 1,3
10	TRMGL(I,J)   = DIRPT(J,I)*DIRPT(J,I)
	
	DO 20 I = 1,3
	TRMGL(I,4)   = DIRPT(1,I)*DIRPT(2,I)
	TRMGL(I,5)   = DIRPT(1,I)*DIRPT(3,I)
20	TRMGL(I,6)   = DIRPT(2,I)*DIRPT(3,I)	

	DO 30 J = 1,3
	TRMGL(4,J)   = 2.0*DIRPT(J,1)*DIRPT(J,2)
	TRMGL(5,J)   = 2.0*DIRPT(J,1)*DIRPT(J,3)
30	TRMGL(6,J)   = 2.0*DIRPT(J,2)*DIRPT(J,3)

	TRMGL(4,4) = DIRPT(1,1)*DIRPT(2,2) + DIRPT(2,1)*DIRPT(1,2)
	TRMGL(4,5) = DIRPT(1,1)*DIRPT(3,2) + DIRPT(3,1)*DIRPT(1,2)
	TRMGL(4,6) = DIRPT(2,1)*DIRPT(3,2) + DIRPT(3,1)*DIRPT(2,2)
	TRMGL(5,4) = DIRPT(1,1)*DIRPT(2,3) + DIRPT(2,1)*DIRPT(1,3)
	TRMGL(5,5) = DIRPT(1,1)*DIRPT(3,3) + DIRPT(3,1)*DIRPT(1,3)
	TRMGL(5,6) = DIRPT(2,1)*DIRPT(3,3) + DIRPT(3,1)*DIRPT(2,3)
	TRMGL(6,4) = DIRPT(1,2)*DIRPT(2,3) + DIRPT(2,2)*DIRPT(1,3)
	TRMGL(6,5) = DIRPT(1,2)*DIRPT(3,3) + DIRPT(3,2)*DIRPT(1,3)
	TRMGL(6,6) = DIRPT(2,2)*DIRPT(3,3) + DIRPT(3,2)*DIRPT(2,3)
	

	RETURN
	END	


C	=======================================================================
C	=======================================================================
	SUBROUTINE	ELTIFF_DM(YOUNG,POISN,STIFG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	DIMENSION STIFG(6,6)
C	=======================================================================

	DO 110 I = 1,6
	DO 110 J = 1,6
110	STIFG(I,J) = 0.0

	STIFG(1,1) = 1.0 - POISN
	STIFG(2,2) = 1.0 - POISN
	STIFG(3,3) = 1.0 - POISN
	STIFG(1,2) = POISN
	STIFG(1,3) = POISN
	STIFG(2,3) = POISN
	STIFG(2,1) = POISN
	STIFG(3,1) = POISN
	STIFG(3,2) = POISN
	STIFG(4,4) = (1.0 - 2.0*POISN)/2.0
	STIFG(5,5) = (1.0 - 2.0*POISN)/2.0
	STIFG(6,6) = (1.0 - 2.0*POISN)/2.0
		 
	DO 150 I = 1,6
	DO 150 J = 1,6
150	STIFG(I,J) = STIFG(I,J)*YOUNG/(1.0+POISN)/(1.0-2.0*POISN)
	

	RETURN
	END	


C	=======================================================================
C	=======================================================================
	SUBROUTINE	DMPAR(AVALU,EMAXS,STRSS,FKO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================

	EVALU = 2.718281828
	CON2  = 2.0*EVALU**(-AVALU*(STRSS-EMAXS))
	CON3  = EVALU**(-2.0*AVALU*(STRSS-EMAXS))
	FKO   = 1.0 - ((EMAXS/STRSS)*(CON2-CON3))**0.5

	IF(FKO.LE.0.0)  FKO = 0.0D0
	IF(FKO.GE.1.0)  FKO = 0.99999
	

	RETURN
	END	


C	=======================================================================
C	==================== END OF DAMAGE MECHANIC MODULE ====================
C	=======================================================================