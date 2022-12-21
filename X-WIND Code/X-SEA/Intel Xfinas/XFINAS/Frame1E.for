C	================================================================
C	================================================================
C	================================================================
      SUBROUTINE OFPROP (PROPO,IFSET,IGIDM,NOPS,NELE,OPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPT

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /sornOffset/ NELO,Nmemoffset(1000),NmemOffsetSelect(1000) ,Nmemtrueelement(1000),GloLoSw(1000) !total number of off set 
      COMMON /MGRAV/ NGRAV
      COMMON A(7000000),IA(6000000)
      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN
      
      DIMENSION PROPO(6,NOPS),IFSET(1),IGIDM(1),PROP(6)
      DIMENSION  SJoint1(3), SJoint2(3)
      DIMENSION VLX(3),VLY(3),VLZ(3)
      DIMENSION VLXn(3),VLYn(3),VLZn(3)
      DIMENSION soffinpi(3),soffinpj(3)
      DIMENSION GloOff(6)
      
      IF (OPT.EQ."READ") THEN
C     -----------------------------------------------
C     FRAME END OFFSET CONDITIONS
C     -----------------------------------------------

C	SET DEFAULT VALUES IF THERE IS NO OFFSET CONDITION
      IFSET(1:NELE) = 0
      
	IF (NOPS.EQ.0) RETURN 
	
	READ (ITI,*)
	READ (ITI,*) NELO
	
	PROP(1:6) = 0
      DO IOPS = 1,NOPS
        READ(ITI,*) II,PROP(1:6),GloLoSw(IOPS)
	  PROPO(1:6,IOPS) = PROP(1:6)
	ENDDO

      IF(NELO.LE.0) RETURN
C     ---------------------------------------
C     ASSIGN  TYPE OF HINGE  TO ELEMENTS
C     ---------------------------------------
      IELE = 0
	READ (ITI,*)
	
      Nsorn=0
      
300   READ (ITI,*)  MEMBA,ISET
      Nsorn=Nsorn+1
      Nmemtrueelement(Nsorn) =  MEMBA
	DO  IGID = 1,NELE
	    MEMGD = IGIDM(IGID) 
	    IF(MEMGD.EQ.MEMBA) THEN
	        MEMBA = IGID 
	        EXIT
	    ENDIF
	ENDDO
	IFSET(MEMBA) = ISET
	
      IELE = IELE+1
      
      Nmemoffset(Nsorn) = MEMBA
      NmemOffsetSelect(Nsorn) = ISET
 !     Nmemtrueelement(Nsorn) = 
      
      
      IF (IELE.LT.NELO)  GOTO 300

      ELSEIF   (OPT.EQ."WRIT") THEN    ! Loop 2 for Transform Local Cooridinate to Global


      !make loop here
      
      
      IF (NELO.EQ.0) RETURN 
      
      
      Do 50011 loop = 1,NELO 
          
          if (GloLoSw(loop) == 2 )  then
          
      iAnaOffCase= NmemOffsetSelect(loop)                !case of offset
      iAnaOffMem=  Nmemtrueelement(loop)                 !no. member
      
 !    NodeNum1
 !    NudeNum2
 !    
 !    Number
      
      
      ! find 2 Joint connect with selected member
      
      !NN = IELEMENT(184)
      !Node1 = N1(iAnaOffMem)
      !Node2 = N2(iAnaOffMem)
      CALL CALLNUMNODE_F (iAnaOffMem,Node1,Node2)
      
      
      ! 2 joints Coordinate member that will be make for Local X-AXIS
      
      !SJoint1(1) = AX(Node1)
      !SJoint1(2) = AY(Node1)
      !SJoint1(3) = AZ(Node1)
      CALL RELFILL('@XYZ',SJoint1(1),1,Node1,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',SJoint1(2),2,Node1,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',SJoint1(3),3,Node1,0)  !GETTING NODAL COORDINATE	
      
      !SJoint2(1) = AX(Node2)
      !SJoint2(2) = AY(Node2)
      !SJoint2(3) = AZ(Node2)
      CALL RELFILL('@XYZ',SJoint2(1),1,Node2,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',SJoint2(2),2,Node2,0)  !GETTING NODAL COORDINATE
	CALL RELFILL('@XYZ',SJoint2(3),3,Node2,0)  !GETTING NODAL COORDINATE	
      
      !  LOCAL X-AXIS VECTOR
     
      
      VLX(1) =  SJoint2(1) - SJoint1(1)      ! node VLX == Vector Local X-AXIS
      VLX(2) =  SJoint2(2) - SJoint1(2)
      VLX(3) =  SJoint2(3) - SJoint1(3)
      
      
      
      ! make Unit vector     
     
      CALL SCALEN(VLX,VLX,ELN,3)
      
      VLX(1) =  VLX(1)
      VLX(2) =  VLX(2)
      VLX(3) =  VLX(3)
      

      ! Find Local axis 
        
      CALL FMVEVR (VLX,VLY,VLZ)        
      
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
          
          ! bring offfset input data
      
       soffinpi(1) = A(LOP+(iAnaOffCase-1)*6)           
       soffinpi(2) = A(LOP+(iAnaOffCase-1)*6 +1)     
       soffinpi(3) = A(LOP+(iAnaOffCase-1)*6 +2)
                     
       soffinpj(1) = A(LOP+(iAnaOffCase-1)*6 +3)
       soffinpj(2) = A(LOP+(iAnaOffCase-1)*6 +4)
       soffinpj(3) = A(LOP+(iAnaOffCase-1)*6 +5)
          
          
       ! do 2 time i loop and j loop   
       iij = 0  
       ic = 0
       
       
       do iij = 1,2   
            
       
       !apply value (offset input) to vector
       
       ! make do iii = 1 for offset i 
       ! and iii = 2 for offset J
       
       if (iij == 1) then
       
       VLXn(1) = VLX(1)*soffinpi(1)
       VLXn(2) = VLX(2)*soffinpi(1)
       VLXn(3) = VLX(3)*soffinpi(1)
       
       VLYn(1) = VLY(1)*soffinpi(2)
       VLYn(2) = VLY(2)*soffinpi(2)
       VLYn(3) = VLY(3)*soffinpi(2)
       
       VLZn(1) = VLZ(1)*soffinpi(3)
       VLZn(2) = VLZ(2)*soffinpi(3)
       VLZn(3) = VLZ(3)*soffinpi(3)
       
       ic = 0
       
       else if (iij == 2) then
           
       VLXn(1) = VLX(1)*soffinpj(1)
       VLXn(2) = VLX(2)*soffinpj(1)
       VLXn(3) = VLX(3)*soffinpj(1)
       
       VLYn(1) = VLY(1)*soffinpj(2)
       VLYn(2) = VLY(2)*soffinpj(2)
       VLYn(3) = VLY(3)*soffinpj(2)
       
       VLZn(1) = VLZ(1)*soffinpj(3)
       VLZn(2) = VLZ(2)*soffinpj(3)
       VLZn(3) = VLZ(3)*soffinpj(3)    
           
       ic = 3   
       
       endif
           
        
       
       !global off set
       
       GloOff(ic+1)   = VLXn(1) + VLYn(1) + VLZn(1)
       GloOff(ic+2) = VLXn(2) + VLYn(2) + VLZn(2)
       GloOff(ic+3) = VLXn(3) + VLYn(3) + VLZn(3)
        
       enddo
      
      ! replace
       
       
       A(LOP+(iAnaOffCase-1)*6)        =  GloOff(1)
       A(LOP+(iAnaOffCase-1)*6 +1)     =  GloOff(2)
       A(LOP+(iAnaOffCase-1)*6 +2)     =  GloOff(3)
       
       A(LOP+(iAnaOffCase-1)*6 +3)     =  GloOff(4)
       A(LOP+(iAnaOffCase-1)*6 +4)     =  GloOff(5)
       A(LOP+(iAnaOffCase-1)*6 +5)     =  GloOff(6)
       
       
          else if (GloLoSw(loop) == 1 ) then
              
              
          endif
          
       
50011      continue 
       
      ENDIF
      
      RETURN
      END
C	================================================================
C	================================================================
C	================================================================
      SUBROUTINE HINGAS(IPIN,IRSET,IGIDM,NHIGE,NELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION IPIN(14,NHIGE),IRSET(1),IGIDM(1),LRPN(14)

C     -----------------------------------------------
C     FRAME END RELEASE CONDITIONS
C     -----------------------------------------------

C	SET DEFAULT VALUES IF THERE IS NO RELEASE CONDITION
      IRSET(1:NELE) = 0
      
	IF (NHIGE.EQ.0) RETURN 
	
	READ (ITI,*)
	READ (ITI,*) NELH
	
	LRPN(1:14) = 0
      DO IHIGE = 1,NHIGE
        READ(ITI,*) II,LRPN(1:6),LRPN(8:13)
	  IPIN(1:14,IHIGE) = LRPN(1:14)
	ENDDO

      IF(NELH.LE.0) RETURN
C     ---------------------------------------
C     ASSIGN  TYPE OF HINGE  TO ELEMENTS
C     ---------------------------------------
      IELE = 0
	READ (ITI,*)
	
300	READ (ITI,*)  MEMBA,ISET
	DO  IGID = 1,NELE
	    MEMGD = IGIDM(IGID) 
	    IF(MEMGD.EQ.MEMBA) THEN
	        MEMBA = IGID 
	        EXIT
	    ENDIF
	ENDDO
	IRSET(MEMBA) = ISET
	
      IELE = IELE+1
      IF (IELE.LT.NELH)  GOTO 300


      RETURN
      END
C	================================================================
C	================================================================
C	================================================================
      SUBROUTINE TT1A (VR,VS,VT,TT1)
	IMPLICIT REAL*8 (A-H,O-Z)
C	CALCULATE THE TRANSFORMATION MATRIX
      DIMENSION TT1(14,14),VR(3),VS(3),VT(3)

	CALL CLEAR (14,14,TT1)
      TT1(1,1)   = VR(1)
      TT1(1,2)   = VS(1)
      TT1(1,3)   = VT(1)
      TT1(2,1)   = VR(2)
      TT1(2,2)   = VS(2)
      TT1(2,3)   = VT(2)
      TT1(3,1)   = VR(3)
      TT1(3,2)   = VS(3)
      TT1(3,3)   = VT(3)
      TT1(4,4)   = VR(1)
      TT1(4,5)   = VS(1)
      TT1(4,6)   = VT(1)
      TT1(5,4)   = VR(2)
      TT1(5,5)   = VS(2)
      TT1(5,6)   = VT(2)
      TT1(6,4)   = VR(3)
      TT1(6,5)   = VS(3)
      TT1(6,6)   = VT(3)
      TT1(7,7)   = 1.
      TT1(8,8)   = VR(1)
      TT1(8,9)   = VS(1)
      TT1(8,10)  = VT(1)
      TT1(9,8)   = VR(2)
      TT1(9,9)   = VS(2)
      TT1(9,10)  = VT(2)
      TT1(10,8)  = VR(3)
      TT1(10,9)  = VS(3)
      TT1(10,10) = VT(3)
      TT1(11,11) = VR(1)
      TT1(11,12) = VS(1)
      TT1(11,13) = VT(1)
      TT1(12,11) = VR(2)
      TT1(12,12) = VS(2)
      TT1(12,13) = VT(2)
      TT1(13,11) = VR(3)
      TT1(13,12) = VS(3)
      TT1(13,13) = VT(3)
      TT1(14,14) = 1.


      RETURN
	END
C	================================================================
C	================================================================
C	================================================================
      SUBROUTINE TT1B (VR,VS,VT,TT1)
	IMPLICIT REAL*8 (A-H,O-Z)
C	CALCULATE THE TRANSFORMATION MATRIX
      DIMENSION TT1(6,6),VR(3),VS(3),VT(3)

	CALL CLEAR (6,6,TT1)
      TT1(1,1)   = VR(1)
      TT1(1,2)   = VS(1)
      TT1(1,3)   = VT(1)
      TT1(2,1)   = VR(2)
      TT1(2,2)   = VS(2)
      TT1(2,3)   = VT(2)
      TT1(3,1)   = VR(3)
      TT1(3,2)   = VS(3)
      TT1(3,3)   = VT(3)
      TT1(4,4)   = VR(1)
      TT1(4,5)   = VS(1)
      TT1(4,6)   = VT(1)
      TT1(5,4)   = VR(2)
      TT1(5,5)   = VS(2)
      TT1(5,6)   = VT(2)
      TT1(6,4)   = VR(3)
      TT1(6,5)   = VS(3)
      TT1(6,6)   = VT(3)
      
      RETURN
	END
C	================================================================
C	================================================================
C	================================================================
