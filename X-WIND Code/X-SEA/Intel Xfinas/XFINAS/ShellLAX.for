C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHSTORL1 (COORD,PROPG,WA,NWG,VR,VS,VT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION COORD(3,1),PROPG(1)
	DIMENSION WA(1),VR(3),VS(3),VT(3)
	DIMENSION VRN(3),VSN(3),VTN(3)


	THETA = PROPG(11)
	LXOPT = INT(PROPG(12))
	CALL SHNEWBV(COORD,VRN,VSN,VTN,THETA,LXOPT)


	N1 = NWG-17
	N2 = N1+2
	WA(N1:N2) = VR(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VS(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VT(1:3)

	N1 = NWG-8
	N2 = N1+2
	WA(N1:N2) = VRN(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VSN(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VTN(1:3)


	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHSTORL2 (WA,NWG,VR,VS,VT,VRN,VSN,VTN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION WA(1),VR(3),VS(3),VT(3)
	DIMENSION VRN(3),VSN(3),VTN(3)


	N1 = NWG-17
	N2 = N1+2
	VR(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VS(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VT(1:3) = WA(N1:N2)

	N1 = NWG-8
	N2 = N1+2
	VRN(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VSN(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VTN(1:3) = WA(N1:N2)



	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHSRLAX (WA,NWG,SIGR,SIGT,NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*6 NAME
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      DIMENSION COORD(3,1),PROPG(1),VRF(3),VSF(3),VTF(3)
      DIMENSION VRO(3),VSO(3),VTO(3),SIGR(1),SIGT(1),SIG(6)

	CALL SHSTORL2 (WA,NWG,VRO,VSO,VTO,VRF,VSF,VTF)

	IOPT = -1
	IF(NAME.EQ.'RESULT') IOPT = 0
	IF(NAME.EQ.'FIBER6') IOPT = 1
	IF(NAME.EQ.'FIBER3') IOPT = 2
	IF(IOPT.EQ.-1) RETURN

	SELECTCASE(IOPT)
	CASE(0)
	CALL STNRLAX (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	CASE(1)
	CALL STNFLAX (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	CASE(2)
	CALL STNLLAX (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	ENDSELECT


      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE STNRLAX (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      DIMENSION VRF(3),VSF(3),VTF(3)
      DIMENSION VRO(3),VSO(3),VTO(3)
	DIMENSION SIGR(1),SIGT(1),SIG(6)


C	RESULTANT MEMBRANE
	SIG(1) = SIGR(1)
	SIG(2) = SIGR(2)
	SIG(3) = 0.0D0
	SIG(4) = SIGR(3)
	SIG(5) = SIGR(7)
	SIG(6) = SIGR(8)

	CALL STSTRAN(SIG,VRO,VSO,VTO,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	CALL STSTRAN(SIG,VRF,VSF,VTF,1)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS

	SIGT(1) = SIG(1)
	SIGT(2) = SIG(2)
	SIGT(6) = SIG(5)
	SIGT(7) = SIG(6)
	SIGT(8) = SIG(4)

C	RESULTANT MOMENT
	SIG(1) = SIGR(4)
	SIG(2) = SIGR(5)
	SIG(3) = 0.0D0
	SIG(4) = SIGR(6)
	SIG(5) = 0.0D0
	SIG(6) = 0.0D0
	CALL STSTRAN(SIG,VRO,VSO,VTO,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	CALL STSTRAN(SIG,VRF,VSF,VTF,1)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	
	SIGT(3) = SIG(1)
	SIGT(4) = SIG(2)
	SIGT(5) = SIG(4)



      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================

      SUBROUTINE STNFLAX (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      DIMENSION VRF(3),VSF(3),VTF(3)
      DIMENSION VRO(3),VSO(3),VTO(3)
	DIMENSION SIGR(1),SIGT(1),SIG(6)


	SIG(1:6) = SIGR(1:6)

	CALL STSTRAN(SIG,VRO,VSO,VTO,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	CALL STSTRAN(SIG,VRF,VSF,VTF,1)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS

	SIGT(1) = SIG(1)
	SIGT(2) = SIG(2)
	SIGT(3) = SIG(3)
	SIGT(4) = SIG(4)
	SIGT(5) = SIG(5)
	SIGT(6) = SIG(6)


      RETURN
      END



C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE STNLLAX (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      DIMENSION VRF(3),VSF(3),VTF(3)
      DIMENSION VRO(3),VSO(3),VTO(3)
	DIMENSION SIGR(1),SIGT(1),SIG(6)


	SIG(1:2) = SIGR(1:2)
	SIG(3)   = 0.0D0
	SIG(4)   = SIGR(3)
	SIG(5:6) = 0.0D0

	CALL STSTRAN(SIG,VRO,VSO,VTO,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	CALL STSTRAN(SIG,VRF,VSF,VTF,1)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS

	SIGT(1) = SIG(1)
	SIGT(2) = SIG(2)
	SIGT(3) = SIG(4)


      RETURN
      END



C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHNEWBV (COORD,VR,VS,VT,THETA,LXOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
C
C     --------------------------------------------------------------
C     CALLING BASE VECTOR FROM PROPG
C     --------------------------------------------------------------
C
      DIMENSION COORD(3,1),VR(3),VS(3),VT(3),V1(3),V2(3),RV(3)

	VR(1:3) = COORD(1:3,2) - COORD(1:3,1)
	VS(1:3) = COORD(1:3,3) - COORD(1:3,1)

      CALL VECPRD (VR,VS,VT)
      CALL SCALEN (VT,VT,DET,3)

	RV(1:3) = [1.0D0,0.0D0,0.0D0]
      CALL VECPRD (VT,RV,V1)
	ELN = SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
	IF(ELN.LT.0.0001) RV(1:3) = [0.0D0,1.0D0,0.0D0]


	IF(LXOPT.EQ.1) THEN
      CALL SCALEN (VR,VR,DET,3)
	VS(1:3) = COORD(1:3,2) - COORD(1:3,3)
      CALL SCALEN (VS,VS,DET,3)
	CALL SCAPRD (VR,RV,DR,3)
	CALL SCAPRD (VS,RV,DS,3)
	RV(1:3) = VR(1:3)
	IF(ABS(DS).GT.ABS(DR)) RV(1:3) = VS(1:3)
	ENDIF


      CALL VECPRD (VT,RV,V1)
      CALL VECPRD (V1,VT,VR)
      CALL VECPRD (VT,VR,VS)
      CALL VECPRD (VR,VS,VT)
      CALL SCALEN (VR,VR,DET,3)
      CALL SCALEN (VS,VS,DET,3)
      CALL SCALEN (VT,VT,DET,3)

	THETR = THETA*0.01745329252
	IF(ABS(THETR).GT.0.0001) THEN
	V1(1:3) = VR(1:3)
	V2(1:3) = VS(1:3)
	V1(1:3) = V1(1:3)*COS(THETR) + V2(1:3)*SIN(THETR)
      CALL VECPRD (VT,V1,VS)
      CALL VECPRD (VS,VT,VR)
      CALL VECPRD (VR,VS,VT)
      CALL SCALEN (VR,VR,DET,3)
      CALL SCALEN (VS,VS,DET,3)
      CALL SCALEN (VT,VT,DET,3)
	ENDIF


      RETURN

      END

C	=====================================================================
C	=====================================================================
C	=====================================================================