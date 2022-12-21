C	=========================================================================
C	=========================================================================
      SUBROUTINE SOLITHIN(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,MWG
     +				,FIN,WORKG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=========================================================================
C	=========================================================================
C	=========================================================================
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK

C	=========================================================================

C
C	STREAM INPUT AND OUTPUT FOR WORKING ARRAY
C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(*),COORD(3,NNO)
      DIMENSION EDIS(24),EDISI(24),RE(24),COORDO(3,NNO)
 
	DIMENSION SK(24,24)

	DIMENSION FIN(NEF)

	DIMENSION TRFNN(3,3),H(8),PR(4),PS(4),SHAP(3,24),COVR(3),COVS(3)
	DIMENSION VR(3),VS(3),VT(3),RV(3),SV(3),BMATX(3,24)
	DIMENSION STRESS(3),STRAIN(3),DEP(3,3),GLOTN(24,24),DEVET(3,24)

	DIMENSION WORKG(7,16)

	DO IRE = 1,24
	RE(IRE) = 0.0
	ENDDO

	KSK = 0
	DO ISK =1,24
	DO JSK =ISK,24
	KSK = KSK+1 
	S(KSK)= 0.0
	END DO
	END DO

	IJCOT = 0
	DO I = 1,8
	DO J = 1,3
	IJCOT = IJCOT + 1
	COORDO(J,I) = COORD(J,I) !- EDIS(IJCOT)
	ENDDO
	ENDDO

C     -----------------------------
C     SETTING INDEX FOR GAUSS POINT
C     -----------------------------

      IPT = 0

C	------------------------
C	 LOOP OVER GAUSS POINT
C	------------------------

      DO 900  IGR=1,4
      RI = GLOC(IGR,4)
      DO 900  IGS=1,4
      SI = GLOC(IGS,4)
	DO 900  IGT=1,1
      TI = GLOC(IGT,1)
	WT = GWT(IGR,4)*GWT(IGS,4)*GWT(IGT,1)
      IPT = IPT+1

C	TESTING OF WHICH SIDE IS ZERO THICKNESS
	XYZ1 = SQRT( (COORDO(1,2)-COORDO(1,1))**2.0 +
	1	         (COORDO(2,2)-COORDO(2,1))**2.0 +
	2	         (COORDO(3,2)-COORDO(3,1))**2.0 )
	XYZ2 = SQRT( (COORDO(1,4)-COORDO(1,1))**2.0 +
	1	         (COORDO(2,4)-COORDO(2,1))**2.0 +
	2	         (COORDO(3,4)-COORDO(3,1))**2.0 )
	XYZ3 = SQRT( (COORDO(1,5)-COORDO(1,1))**2.0 +
	1	         (COORDO(2,5)-COORDO(2,1))**2.0 +
	2	         (COORDO(3,5)-COORDO(3,1))**2.0 )

	NMX = 0
	IF(XYZ1.LT.XYZ2.AND.XYZ1.LT.XYZ3) NMX = 1
	IF(XYZ2.LT.XYZ1.AND.XYZ2.LT.XYZ3) NMX = 2
	IF(XYZ3.LT.XYZ1.AND.XYZ3.LT.XYZ2) NMX = 3

	IF(NMX.EQ.0) THEN
	IF(XYZ1.EQ.XYZ2) NMX = 1
	IF(XYZ1.EQ.XYZ3) NMX = 1
	IF(XYZ2.EQ.XYZ3) NMX = 2
	ENDIF

	IF(NMX.EQ.1) THEN
	THK = XYZ1
	ELSEIF(NMX.EQ.2) THEN
	THK = XYZ2
	ELSEIF(NMX.EQ.3) THEN
	THK = XYZ2
	ENDIF




	H(1)  = 0.25*(1.0+RI)*(1.0+SI)/THK
	H(2)  = 0.25*(1.0-RI)*(1.0+SI)/THK
	H(3)  = 0.25*(1.0-RI)*(1.0-SI)/THK
	H(4)  = 0.25*(1.0+RI)*(1.0-SI)/THK

	H(5)  = -0.25*(1.0+RI)*(1.0+SI)/THK
	H(6)  = -0.25*(1.0-RI)*(1.0+SI)/THK
	H(7)  = -0.25*(1.0-RI)*(1.0-SI)/THK
	H(8)  = -0.25*(1.0+RI)*(1.0-SI)/THK


	SHAP = 0.0
	IF(NMX.EQ.1) THEN
	SHAP(1,1 ) =  H(1)
	SHAP(1,4 ) = -H(1)
	SHAP(1,7 ) = -H(2)
	SHAP(1,10) =  H(2)
	SHAP(1,13) =  H(4)
	SHAP(1,16) = -H(4)
	SHAP(1,19) = -H(3)
	SHAP(1,22) =  H(3)
	ELSEIF(NMX.EQ.2) THEN
	SHAP(1,1 ) =  H(2)
	SHAP(1,4 ) =  H(1)
	SHAP(1,7 ) = -H(1)
	SHAP(1,10) = -H(2)
	SHAP(1,13) =  H(3)
	SHAP(1,16) =  H(4)
	SHAP(1,19) = -H(4)
	SHAP(1,22) = -H(3)
	ELSEIF(NMX.EQ.3) THEN
	SHAP(1,1 ) =  H(4)
	SHAP(1,4 ) =  H(1)
	SHAP(1,7 ) =  H(2)
	SHAP(1,10) =  H(3)
	SHAP(1,13) = -H(4)
	SHAP(1,16) = -H(1)
	SHAP(1,19) = -H(2)
	SHAP(1,22) = -H(3)
	ENDIF

	PR(1) =  0.25*(1.0+SI)
	PR(2) = -0.25*(1.0+SI)
	PR(3) = -0.25*(1.0-SI)
	PR(4) =  0.25*(1.0-SI)

	PS(1) =  0.25*(1.0+RI)
	PS(2) =  0.25*(1.0-RI)
	PS(3) = -0.25*(1.0-RI)
	PS(4) = -0.25*(1.0+RI)

	DO I = 1,8
	SHAP(2,3*I-1) = SHAP(1,3*I-2)
	SHAP(3,3*I-0) = SHAP(1,3*I-2)
	ENDDO


	COVR = 0.0
	COVS = 0.0
	IF(NMX.EQ.1) THEN
	COVR(1) = PR(1)*COORDO(1,1) + PR(2)*COORDO(1,4) +
	1		  PR(3)*COORDO(1,8) + PR(4)*COORDO(1,5)
	COVR(2) = PR(1)*COORDO(2,1) + PR(2)*COORDO(2,4) +
	1		  PR(3)*COORDO(2,8) + PR(4)*COORDO(2,5)
	COVR(3) = PR(1)*COORDO(3,1) + PR(2)*COORDO(3,4) +
	1		  PR(3)*COORDO(3,8) + PR(4)*COORDO(3,5)
	COVS(1) = PS(1)*COORDO(1,1) + PS(2)*COORDO(1,4) +
	1		  PS(3)*COORDO(1,8) + PS(4)*COORDO(1,5)
	COVS(2) = PS(1)*COORDO(2,1) + PS(2)*COORDO(2,4) +
	1		  PS(3)*COORDO(2,8) + PS(4)*COORDO(2,5)
	COVS(3) = PS(1)*COORDO(3,1) + PS(2)*COORDO(3,4) +
	1		  PS(3)*COORDO(3,8) + PS(4)*COORDO(3,5)
	ELSEIF(NMX.EQ.2) THEN
	COVR(1) = PR(1)*COORDO(1,2) + PR(2)*COORDO(1,1) +
	1		  PR(3)*COORDO(1,5) + PR(4)*COORDO(1,6)
	COVR(2) = PR(1)*COORDO(2,2) + PR(2)*COORDO(2,1) +
	1		  PR(3)*COORDO(2,5) + PR(4)*COORDO(2,6)
	COVR(3) = PR(1)*COORDO(3,2) + PR(2)*COORDO(3,1) +
	1		  PR(3)*COORDO(3,5) + PR(4)*COORDO(3,6)
	COVS(1) = PS(1)*COORDO(1,2) + PS(2)*COORDO(1,1) +
	1		  PS(3)*COORDO(1,5) + PS(4)*COORDO(1,6)
	COVS(2) = PS(1)*COORDO(2,2) + PS(2)*COORDO(2,1) +
	1		  PS(3)*COORDO(2,5) + PS(4)*COORDO(2,6)
	COVS(3) = PS(1)*COORDO(3,2) + PS(2)*COORDO(3,1) +
	1		  PS(3)*COORDO(3,5) + PS(4)*COORDO(3,6)
	ELSEIF(NMX.EQ.3) THEN
	COVR(1) = PR(1)*COORDO(1,2) + PR(2)*COORDO(1,3) +
	1		  PR(3)*COORDO(1,4) + PR(4)*COORDO(1,1)
	COVR(2) = PR(1)*COORDO(2,2) + PR(2)*COORDO(2,3) +
	1		  PR(3)*COORDO(2,4) + PR(4)*COORDO(2,1)
	COVR(3) = PR(1)*COORDO(3,2) + PR(2)*COORDO(3,3) +
	1		  PR(3)*COORDO(3,4) + PR(4)*COORDO(3,1)
	COVS(1) = PS(1)*COORDO(1,2) + PS(2)*COORDO(1,3) +
	1		  PS(3)*COORDO(1,4) + PS(4)*COORDO(1,1)
	COVS(2) = PS(1)*COORDO(2,2) + PS(2)*COORDO(2,3) +
	1		  PS(3)*COORDO(2,4) + PS(4)*COORDO(2,1)
	COVS(3) = PS(1)*COORDO(3,2) + PS(2)*COORDO(3,3) +
	1		  PS(3)*COORDO(3,4) + PS(4)*COORDO(3,1)
	ENDIF

	CALL VECPRD (COVR,COVS,VT)
	CALL SCALEN (VT,VT,DET,3)
      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)

      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)
	CALL SCALEN (VR,VR,DM,3)
	CALL SCALEN (VT,VT,DM,3)

	TRFNN(1,1) = VR(1)
	TRFNN(2,1) = VR(2)
	TRFNN(3,1) = VR(3)

	TRFNN(1,2) = VS(1)
	TRFNN(2,2) = VS(2)
	TRFNN(3,2) = VS(3)

	TRFNN(1,3) = VT(1)
	TRFNN(2,3) = VT(2)
	TRFNN(3,3) = VT(3)

C	WRITE(*,*) TRFNN(1,1),TRFNN(1,2),TRFNN(1,3)
C	WRITE(*,*) TRFNN(2,1),TRFNN(2,2),TRFNN(2,3)
C	WRITE(*,*) TRFNN(3,1),TRFNN(3,2),TRFNN(3,3)
C	PAUSE

	DO IT = 1,3
	DO JT = 1,3
	IF(ABS(TRFNN(IT,JT)).LE.1.0E-6) TRFNN(IT,JT) = 0.0
	ENDDO
	ENDDO


	GLOTN = 0.0
	DO IT = 1,8
	GLOTN(3*IT-2,3*IT-2) = TRFNN(1,1)
	GLOTN(3*IT-2,3*IT-1) = TRFNN(1,2)
	GLOTN(3*IT-2,3*IT-0) = TRFNN(1,3)
	GLOTN(3*IT-1,3*IT-2) = TRFNN(2,1)
	GLOTN(3*IT-1,3*IT-1) = TRFNN(2,2)
	GLOTN(3*IT-1,3*IT-0) = TRFNN(2,3)
	GLOTN(3*IT-0,3*IT-2) = TRFNN(3,1)
	GLOTN(3*IT-0,3*IT-1) = TRFNN(3,2)
	GLOTN(3*IT-0,3*IT-0) = TRFNN(3,3)
	ENDDO


	IF(DET.LE.0.0) THEN
	WRITE(*,*) 'Jacobian is not Poisitive Definition'
	WRITE(*,*) 'Pleas Check in the Thin Soid Element'
	STOP
	ENDIF

      DVOL = WT*DET*0.5*THK


	BMATX = MATMUL(SHAP,TRANSPOSE(GLOTN))


	UL = 0.0
	UM = 0.0
	UN = 0.0
	DO I = 1,24
	UL = UL + BMATX(1,I)*EDIS(I)
	UM = UM + BMATX(2,I)*EDIS(I)
	UN = UN + BMATX(3,I)*EDIS(I)
	ENDDO

      STRAIN(1) = UL
      STRAIN(2) = UM
      STRAIN(3) = UN


	CALL THIN3D(WORKG(1,IPT),WORKG(4,IPT),
	1			STRESS,STRAIN,DEP,PROPM,WORKG(7,IPT))


	SK = 0.0
	SK = MATMUL(TRANSPOSE(BMATX),MATMUL(DEP,BMATX))*DVOL


	DO IRE = 1,24
	DO JRE = 1,3
	RE(IRE) = RE(IRE) + BMATX(JRE,IRE)*STRESS(JRE)*DVOL
	ENDDO
	ENDDO

	
	KSK = 0
	DO ISK =1,24
	DO JSK =ISK,24
	KSK = KSK+1 
	S(KSK)=S(KSK)+SK(ISK,JSK)
	END DO
	END DO


900	CONTINUE



	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF



      RETURN
      END

C	=========================================================================
C	=========================================================================
C	=========================================================================
