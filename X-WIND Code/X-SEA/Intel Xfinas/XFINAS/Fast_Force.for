!     ==================================================================================================
      SUBROUTINE FASTTOPFORCE_DUMMY (IFAST,KFAST,NODEFAST,NTFASTPARA,NUMCASE,NFASTPARA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      

      Fast_Dummy = 0  
      
      !=====================================
      !  Fast_Dummy = 0  USE NORMAL TURBINE
      !  Fast_Dummy = 1  USE DUMMY TURBINE
      !=====================================
      
      if (Fast_Dummy == 0 ) return
      
C     FAST CONDITION      
!      IF (KFAST.EQ.2)THEN
!      READ (ITI,*) 
!      READ (ITI,*) NODEFAST
!      READ (ITI,*)
!      READ (ITI,*) NTFASTPARA
!          DO I = 1,NTFASTPARA
!          READ (ITI,*) NUMCASE(I),NFASTPARA(I)
!          ENDDO
!      ENDIF
      
      NODEFAST = 1
      
      IFAST = 2
      KFAST = 2 
      NTFASTPARA = 1
      NUMCASE = 1
      NFASTPARA = 1
      
      RETURN
      End SUBROUTINE
!     ==================================================================================================   
      SUBROUTINE FASTTOPFORCE (ID,MSF,NODEA,NCASE,NFPR,NTFPR,R,NEQ,ILC,OPTION)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION R(NEQ)
      DIMENSION NCASE(1000),NFPR(1000)
      CHARACTER*4 OPTION
      CHARACTER(LEN=200)::LINE
      DIMENSION ID(MSF,1)
      DIMENSION FFX(500),FFY(500),FFZ(500),FAMX(500),FAMY(500),FAMZ(500)
      COMMON NODEAFF(500),NODEAFX(500)
      CHARACTER(100) TurbStaFileName,TurbDynFileName,NumberofPara,ROOTFST
      COMMON / NODEFAST / NODEAF,NUMCASE(1000),NFASTPARA(1000)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /MGRAV/ NGRAV
      COMMON /NFAST/ FX,FY,FZ,AMX,AMY,AMZ,IILC
      COMMON /COUANA/ COUPLEANALYSIS
      
      IF (OPTION.EQ."WRIT")THEN ! WRITE 
      NODEAF             = NODEA
      NUMCASE(1:NTFPR)   = NCASE(1:NTFPR)
      NFASTPARA(1:NTFPR) = NFPR(1:NTFPR)
      
      ELSEIF (OPTION.EQ."READ")THEN ! LOAD CASE ( STATIC )
      OPEN(UNIT=553,FILE="INTEGER_TO_STRING"  ,STATUS='UNKNOWN')    
      WRITE (553,50) NFASTPARA(ILC)
50    FORMAT (I5)
      REWIND (553)
      READ (553,*) NumberofPara 
      CLOSE(553)
      ROOTFST         = '.dat'
      ! MODIFY BY TOEY 04/2018
      ! TurbStaFileName = "FASTLoadSta"//"/"//"TowerTopLoads"//"/"//"TurbStaFile"//"/"//TRIM(NumberofPara)//TRIM(ROOTFST)    
      ! TurbDynFileName = "FASTLoadDyn"//"/"//"TowerTopLoads"//"/"//"TurbDynFile"//"/"//TRIM(NumberofPara)//TRIM(ROOTFST)
      TurbStaFileName = "FASTLoadSta"//"/"//"TowerTopLoads"//"/"//"TurbStaFile"//TRIM(NumberofPara)//TRIM(ROOTFST)  
      TurbDynFileName = "FASTLoadDyn"//"/"//"TowerTopLoads"//"/"//"TurbDynFile"//TRIM(NumberofPara)//TRIM(ROOTFST)
      OPEN(UNIT=551,FILE=TurbStaFileName,STATUS='UNKNOWN')
      OPEN(UNIT=552,FILE=TurbDynFileName,STATUS='UNKNOWN')
      
      READ (551,'(A)',IOSTAT=ios) LINE
C      IF (LEN_TRIM(LINE).NE.22) RETURN
      READ (551,*) FX,FY,FZ,AMX,AMY,AMZ
      
      If( NGRAV .EQ. 2 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMZ*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+(-AMY*1000D0)
      ENDDO
      Elseif( NGRAV .EQ. 3 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMY*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+AMZ*1000D0
      ENDDO          
      Endif

      ELSEIF (OPTION.EQ."REDD") THEN ! LOAD CASE 1 ( DYNAMIC )
      IF (COUPLEANALYSIS.EQ.2.OR.COUPLEANALYSIS.EQ.3) RETURN
      OPEN(UNIT=553,FILE="INTEGER_TO_STRING"  ,STATUS='UNKNOWN')    
      WRITE (553,100) NFASTPARA(1) ! only single case for dynamic
100   FORMAT (I5)
      REWIND (553)
      READ (553,*) NumberofPara 
      CLOSE(553)
      IF (IILC.NE.ILC)THEN
      !NumberofPara = "1"
      !NumberofPara    = NFASTPARA(ILC)
      ROOTFST         = '.dat'
      TurbStaFileName = "FASTLoadSta"//"/"//"TowerTopLoads"//"/"//"TurbStaFile"//TRIM(NumberofPara)//TRIM(ROOTFST)    
      TurbDynFileName = "FASTLoadDyn"//"/"//"TowerTopLoads"//"/"//"TurbDynFile"//TRIM(NumberofPara)//TRIM(ROOTFST)
      OPEN(UNIT=551,FILE=TurbStaFileName,STATUS='UNKNOWN')
      OPEN(UNIT=552,FILE=TurbDynFileName,STATUS='UNKNOWN')
      
C      IF (ILC.EQ.1)THEN ! FIRST LOOP
C      READ (552,'(A)',IOSTAT=ios) LINE
C      IF (LEN_TRIM(LINE).NE.23) RETURN
C      READ (552,*) 
C      READ (552,*) 
C      ENDIF 
      
      READ (552,*) FX,FY,FZ,AMX,AMY,AMZ
      
C     KN TO N
      if( NGRAV .EQ. 2 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMZ*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+(-AMY*1000D0)
      ENDDO
      Elseif( NGRAV .EQ. 3 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMY*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+AMZ*1000D0
      ENDDO          
      Endif
      IILC = ILC
      ELSEIF (IILC.NE.ILC)THEN
          
      if( NGRAV .EQ. 2 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMZ*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+(-AMY*1000D0)
      ENDDO
      Elseif( NGRAV .EQ. 3 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMY*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+AMZ*1000D0
      ENDDO          
      Endif
      ENDIF
      
      ELSEIF (OPTION.EQ."REDT") THEN ! LOAD CASE 1 ( DYNAMIC )
      IF (COUPLEANALYSIS.EQ.2.OR.COUPLEANALYSIS.EQ.3) RETURN
      OPEN(UNIT=553,FILE="INTEGER_TO_STRING"  ,STATUS='UNKNOWN')    
      WRITE (553,200) NFASTPARA(1)
200   FORMAT (I5)
      REWIND (553)
      READ (553,*) NumberofPara 
      CLOSE(553)
c      NumberofPara    = "1"
      ROOTFST         = '.dat'
      TurbStaFileName = "FASTLoadSta"//"/"//"TowerTopLoads"//"/"//"TurbStaFile"//TRIM(NumberofPara)//TRIM(ROOTFST)    
      TurbDynFileName = "FASTLoadDyn"//"/"//"TowerTopLoads"//"/"//"TurbDynFile"//TRIM(NumberofPara)//TRIM(ROOTFST)
      OPEN(UNIT=551,FILE=TurbStaFileName,STATUS='UNKNOWN')
      OPEN(UNIT=552,FILE=TurbDynFileName,STATUS='UNKNOWN')
      
      READ (552,*) 
      READ (552,*) 
      
      ELSEIF (OPTION.EQ."REDA") THEN ! FOR TEST SUBDYN
      OPEN(UNIT=551,FILE="FORCE AND MOMENT FOR HD.TXT",STATUS='UNKNOWN')
      
      READ (551,*) DUMMY,NNODEFAST
      DO J=1,NNODEFAST
      READ (551,*) FFX(J),FFY(J),FFZ(J),FAMX(J),FAMY(J),FAMZ(J)
      ENDDO
      
      DO J=1,NNODEFAST
      DO I = 1,6
      IEQ = ID(I,NODEAFX(J))
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FFX(NODEAFF(J))
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FFY(NODEAFF(J))
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FFZ(NODEAFF(J))
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+FAMX(NODEAFF(J))
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+FAMY(NODEAFF(J))
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+FAMZ(NODEAFF(J))
      ENDDO  
      ENDDO
      
      ELSEIF (OPTION.EQ."REDB") THEN ! FOR TEST SUBDYN
      OPEN(UNIT=552,FILE="NODE ASSIGN.TXT",STATUS='UNKNOWN')
      READ (552,*) NELEMENTFAST
      DO J=1,NELEMENTFAST
      READ (552,*) NODEAFF(J),NODEAFX(J)
      ENDDO
      ENDIF
      
      return

!30    write(*,*) 'Over turbine load numbers'
 !     STOP
      
      END
!     ==================================================================================================
      SUBROUTINE COORDINATE_TRANSFER
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      ALLOCATABLE VALN(:,:)
      
C     CALLING NUMBER OF NODE NSN AND DIMENSION NSC ... SONGSAK OCT2019      
      CALL LOCATN('@XYZ',KXYZ,NSC,NSN,2) 
      
      WRITE(6000,1) 
1      FORMAT('---- STRUCTURE JOINTS: joints connect structure members (~Hydrodyn Input File)---')
      WRITE(6000,2)  NSN
2     FORMAT(I10,'  NJoints     - Number of joints (-)')      
      WRITE(6000,3)
3     FORMAT('JointID          JointXss               JointYss               JointZss 
     1 [Coordinates of Member joints in SS-Coordinate System]')
      WRITE(6000,4)
4     FORMAT('  (-)               (m)                    (m)                    (m)') 
      
      
C      DO I=1,NNO
C      WRITE(6000,10)NND(I),AX(I),AY(I),AZ(I)
C      ENDDO

C     ---------------------------------------      
C     SONGSAK MODIFY TO PRINT FASTER  OCT2019      
      ALLOCATE(VALN(4,NSN))
      
      DO ISN = 1,NSN
          VALN(1,ISN) = FLOAT(ISN)
          DO ISC = 1,NSC
	        CALL RELFILL('@XYZ',VALN(ISC+1,ISN),ISC,ISN,0)  !GETTING NODAL COORDINATE
          ENDDO
      ENDDO
      
      
      WRITE(6000,20) VALN
      
      DEALLOCATE(VALN)
C     ---------------------------------------  

10    FORMAT(I4,5X,F15.5,9X,F15.5,8X,F15.5)
20    FORMAT(F10.1,5X,F15.5,9X,F15.5,8X,F15.5)
      
      RETURN
      END SUBROUTINE
!     ==================================================================================================
      SUBROUTINE BUNDARY_CONDITION_TRANSFER
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      ALLOCATABLE IDRCT(:,:)
      
      
C     SONGSAK NEWLY IMPLEMENT ... OCT2019
          
	CALL LOCATN('@SUP',KRID,NDAT,NSUP,1)  !CALLING NUMBER OF SUPPORT NSUP
      
      ALLOCATE(IDRCT(7,NSUP))
      IDRCT = 0 !INITIALIZE

      DO ISUP = 1,NSUP
          DO ISAT = 1,7
	        CALL INTFILL('@SUP',IDRCT(ISAT,ISUP),ISAT,ISUP,0)  !CALLING SUPPORT DATA ... ONLY FIRST 6 DOF FOR FRAME ELEMENT
          ENDDO
      ENDDO
  
      
      WRITE(6000,1)
1     FORMAT('------------------- BASE REACTION JOINTS: 1/0 for Locked/Free DOF @ each Reaction Node ---------------------')
      WRITE(6000,2)  NSUP
2     FORMAT(I10,'   NReact      - Number of Joints with reaction forces; be sure to remove all rigid motion DOFs of the structure  
     1(else det([K])=[0])')   
      WRITE(6000,3)
3     FORMAT('RJointID   RctTDXss    RctTDYss    RctTDZss    RctRDXss    RctRDYss    RctRDZss     [Global Coordinate System]')
      WRITE(6000,4)
4     FORMAT('  (-)       (flag)      (flag)      (flag)      (flag)      (flag)      (flag)') 
           
      IF(NSUP.GT.0) THEN
          WRITE(6000,10) IDRCT(1:7,1:NSUP)
      ENDIF
      
      
      DEALLOCATE(IDRCT)
      
      
10    FORMAT(I10,2X,I10,2X,I10,2X,I10,2X,I10,2X,I10,2X,I10)
           
      RETURN
      END SUBROUTINE
!     ==================================================================================================
      SUBROUTINE CONNECTIVITY_FRAME_FAST
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      ALLOCATABLE LVAL(:,:)
      
C     SONGSAK NEWLY IMPLEMENT ... OCT2019
      CALL CALLNUMFRAME (NFRAM)
      
      WRITE(6000,1)
1     FORMAT('----------------------------------- MEMBERS --------------------------------------')
      WRITE(6000,2)  NFRAM
2     FORMAT(I10,'   NMembers    - Number of frame members')   
      WRITE(6000,3)
3     FORMAT('MemberID   MJointID1   MJointID2   MPropSetID1   MPropSetID2     COSMID')
      WRITE(6000,4)
4     FORMAT('  (-)         (-)         (-)          (-)           (-)           (-)') 
      
      
CC      DO I = 1,NELE
CC          CALL  MATERIALDATA_FOR_FAST(IELEMENT(I),IEGRP)
CC          WRITE(6000,10)IELEMENT(I),N1(I),N2(I),IEGRP,IEGRP 
CC      ENDDO
CC10    FORMAT(I8,3X,I10,2X,I10,2X,I10,3X,I10,3X,I10) 

      
C     ---------------------------------------      
C     SONGSAK MODIFY TO PRINT FASTER  OCT2019   
      
C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
      
      ALLOCATE(LVAL(5,NGIDM))
      
      II = 0
      DO 100 IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !IEG
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)  !IELE
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
      
      CALL INTFILL('@GMP',IEGRP,1 ,IGM,0)  !CALLING GID GROUP NUMBER ... SEE ALSO SUBROUTINE GIDMAP
      
      IF(ITYP.EQ.5) THEN
          CALL CALLNUMNODE_F (IGM,N1,N2)
          II = II+1
          LVAL(1,II) = IGM
          LVAL(2,II) = N1
          LVAL(3,II) = N2
          LVAL(4,II) = IEGRP
          LVAL(5,II) = IEGRP
      ENDIF
      
100   ENDDO
      
      IF(II.GT.0) WRITE(6000,20) LVAL(1:5,1:II)
      
      DEALLOCATE(LVAL)
C     ---------------------------------------    
      
      CALL MATPOPFAST
      
20    FORMAT(I8,3X,I10,2X,I10,2X,I10,3X,I10)   
      
      RETURN
      END SUBROUTINE
!     ==================================================================================================
      SUBROUTINE PROPERTY_TEMPORY_FAST(ISET,PROP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /MATERIAL_PROP_FAST/ ELASTICMODULOUSF(999),SHEARMODULOUSF(999),DENSITYF(999),NTOTALPROPERTY
      DIMENSION PROP(20)
      
      
      ELASTICMODULOUSF(ISET) = PROP(1)
      SHEARMODULOUSF(ISET) = ELASTICMODULOUSF(ISET)/(2.0D0*(1+PROP(2)))
      DENSITYF(ISET) = PROP(5)
      NTOTALPROPERTY = ISET

      
      END SUBROUTINE
!     ==================================================================================================
      SUBROUTINE MATPOPFAST
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /MATERIAL_PROP_FAST/ ELASTICMODULOUSF(999),SHEARMODULOUSF(999),DENSITYF(999),NTOTALPROPERTY
      COMMON /DISEC/ DIMENPROP(5000,11)
      
      WRITE(6000,1)
1     FORMAT('------------------ MEMBER X-SECTION PROPERTY data 1/2 [isotropic material for now: use this table for circular-tubular
     1 elements] ------------------------')
      WRITE(6000,2)  NTOTALPROPERTY
2     FORMAT(I10,'   NPropSets   - Number of structurally unique x-sections (i.e. how many groups of X-sectional properties are 
     1utilized throughout all of the members)')   
      WRITE(6000,3)
3     FORMAT('PropSetID     YoungE          ShearG          MatDens          XsecD           XsecT')
      WRITE(6000,4)
4     FORMAT('  (-)         (N/m2)          (N/m2)          (kg/m3)           (m)             (m)') 
      
      DO I=1,NTOTALPROPERTY
          
          OD = DIMENPROP(I,2)*2
          THICKNESS = DIMENPROP(I,2)-DIMENPROP(I,3)
      WRITE(6000,10)I,ELASTICMODULOUSF(I),SHEARMODULOUSF(I),DENSITYF(I),DIMENPROP(I,2)*2.0d0,THICKNESS
      ENDDO
10    FORMAT(I4,3X,E15.5,1X,E15.5,3X,F12.2,1X,F14.6,1X,F14.6)  
      
      
      CALL MATPOPSACS
      CALL MEMBER_SACS

      END SUBROUTINE
!     ==================================================================================================
      SUBROUTINE MEMBER_SACS
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)  

      ALLOCATABLE LVAL(:,:)
      
      
      WRITE(6005,1)
1     FORMAT('MEMBER')      

CC      DO I=1,NELE
CC          CALL  MATERIALDATA_FOR_FAST(IELEMENT(I),IEGRP)
CC          WRITE(6005,2)N1(I),N2(I),IEGRP
CC      ENDDO
CC2     FORMAT('MEMBER',1X,I4,I4,1X,I3)

      
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
      
      CALL INTFILL('@GMP',IEGRP,1 ,IGM,0)  !CALLING GID GROUP NUMBER ... SEE ALSO SUBROUTINE GIDMAP
      
      IF(ITYP.EQ.5) THEN
          CALL CALLNUMNODE_F (IGM,N1,N2)
          II = II+1
          LVAL(1,II) = N1
          LVAL(2,II) = N2
          LVAL(3,II) = IEGRP
      ENDIF
      
100   ENDDO
      
      IF(II.GT.0) WRITE(6005,3) LVAL(1:3,1:II)
      
      DEALLOCATE(LVAL)
C     ---------------------------------------          
      
      CALL SACS_JOINT
      
      
3     FORMAT('MEMBER',1X,I10,2X,I10,2X,I4)
      
      END SUBROUTINE
!     ==================================================================================================
      SUBROUTINE SACS_JOINT
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N) 

      ALLOCATABLE VALN(:,:)
            
      WRITE(6005,3)
3     FORMAT('JOINT')
      
      
C      DO I=1,NNO
C            WRITE(6005,4)NND(I),AX(I),AY(I),AZ(I)
C      ENDDO
C4     FORMAT('JOINT',1X,I4,1X,f14.3,X,f14.3,X,f14.3) 
      
      
C     CALLING NUMBER OF NODE NSN AND DIMENSION NSC ... SONGSAK OCT2019      
      CALL LOCATN('@XYZ',KXYZ,NSC,NSN,2)   
      
C     ---------------------------------------      
C     SONGSAK MODIFY TO PRINT FASTER  OCT2019      
      ALLOCATE(VALN(4,NSN))
      
      DO ISN = 1,NSN
          VALN(1,ISN) = FLOAT(ISN)
          DO ISC = 1,NSC
	        CALL RELFILL('@XYZ',VALN(ISC+1,ISN),ISC,ISN,0)  !GETTING NODAL COORDINATE
          ENDDO
      ENDDO
      
      
      WRITE(6005,5) VALN
      
      DEALLOCATE(VALN)
C     ---------------------------------------        
      
      
5     FORMAT('JOINT',1X,F10.1,5X,F15.5,9X,F15.5,8X,F15.5)
      
      END SUBROUTINE
!     ==================================================================================================     
!     ==================================================================================================
      SUBROUTINE MATPOPSACS
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /MATERIAL_PROP_FAST/ ELASTICMODULOUSF(999),SHEARMODULOUSF(999),DENSITYF(999),NTOTALPROPERTY
      COMMON /DISEC/ DIMENPROP(5000,11)
      
      
      WRITE(6005,1)
1     FORMAT('SECT')
      DO I=1,NTOTALPROPERTY
          THICKNESS = DIMENPROP(I,2)-DIMENPROP(I,3)
          WRITE(6005,2)I,DIMENPROP(I,2)*2.0d0*100.0D0,THICKNESS*100.0D0
2     FORMAT('SECT',1X,I7,3X,'TUB',31X,F6.1,F5.2)
      ENDDO
      
      WRITE(6005,3)
3     FORMAT('GRUP')
      DO I=1,NTOTALPROPERTY
          
      WRITE(6005,4)I,I,ELASTICMODULOUSF(I)/10.0D0**10.0D0,SHEARMODULOUSF(I)/10.0D0**10.0D0,24.8,1,1.0D0,1.0D0,0.5,DENSITYF(I)/1000.0D0
4         FORMAT('GRUP',1X,I3,1X,I7,14X,F5.2,F5.2,F5.2,1X,I1,4X,F4.2,F4.2,5X,F5.3,1X,F6.4)
      
      ENDDO
      
      
      !WRITE(6000,10)I,ELASTICMODULOUSF(I),SHEARMODULOUSF(I),DENSITYF(I),DIMENPROP(I,2)*2.0d0,THICKNESS
 
      END SUBROUTINE
!     ==================================================================================================   
      
!     ==================================================================================================
      SUBROUTINE FASTTOPFORCE_DUMMY_FATIGUE (IFAST,KFAST,NODEFAST,NTFASTPARA,NUMCASE,NFASTPARA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      

      Fast_Dummy = 0  
      
      !=====================================
      !  Fast_Dummy = 0  USE NORMAL TURBINE
      !  Fast_Dummy = 1  USE DUMMY TURBINE
      !=====================================
      
      if (Fast_Dummy == 0 ) return
      
C     FAST CONDITION      
!      IF (KFAST.EQ.2)THEN
!      READ (ITI,*) 
!      READ (ITI,*) NODEFAST
!      READ (ITI,*)
!      READ (ITI,*) NTFASTPARA
!          DO I = 1,NTFASTPARA
!          READ (ITI,*) NUMCASE(I),NFASTPARA(I)
!          ENDDO
!      ENDIF
      
      NODEFAST = 1
      
      IFAST = 2
      KFAST = 2 
      NTFASTPARA = 1
      NUMCASE = 1
      NFASTPARA = 1
      
      RETURN
      End SUBROUTINE
!     ==================================================================================================   
      SUBROUTINE FASTTOPFORCE_FATIGUE (ID,MSF,NODEA,NCASE,NFPR,NTFPR,R,NEQ,ILC,FILENUMBER,OPTION)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION R(NEQ)
      DIMENSION NCASE(1000),NFPR(1000)
      CHARACTER*4 OPTION
      CHARACTER*10 FILENUMBER
      CHARACTER(LEN=200)::LINE
      DIMENSION ID(MSF,1)
      DIMENSION FFX(500),FFY(500),FFZ(500),FAMX(500),FAMY(500),FAMZ(500)
      COMMON NODEAFF(500),NODEAFX(500)
      CHARACTER(100) TurbStaFileName,TurbDynFileName,NumberofPara,ROOTFST
      COMMON / NODEFAST / NODEAF,NUMCASE(1000),NFASTPARA(1000)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /MGRAV/ NGRAV
      COMMON /NFAST/ FX,FY,FZ,AMX,AMY,AMZ,IILC
      COMMON /COUANA/ COUPLEANALYSIS
      
      IF (OPTION.EQ."WRIT")THEN ! WRITE 
      NODEAF             = NODEA
      NUMCASE(1:NTFPR)   = NCASE(1:NTFPR)
      NFASTPARA(1:NTFPR) = NFPR(1:NTFPR)
      
      ELSEIF (OPTION.EQ."READ")THEN ! LOAD CASE ( STATIC )
      OPEN(UNIT=553,FILE="INTEGER_TO_STRING"  ,STATUS='UNKNOWN')    
      WRITE (553,50) NFASTPARA(ILC)
50    FORMAT (I5)
      REWIND (553)
      READ (553,*) NumberofPara 
      CLOSE(553)
      ROOTFST         = '.dat'
      TurbStaFileName = "FASTLoadSta"//"/"//"TowerTopLoads"//"/"//"TurbStaFile"//"/"//TRIM(NumberofPara)//TRIM(ROOTFST)    
      TurbDynFileName = "FASTLoadDyn"//"/"//"TowerTopLoads"//"/"//"TurbDynFile"//"/"//TRIM(NumberofPara)//TRIM(ROOTFST)
      OPEN(UNIT=551,FILE=TurbStaFileName,STATUS='UNKNOWN')
      OPEN(UNIT=552,FILE=TurbDynFileName,STATUS='UNKNOWN')
      
      READ (551,'(A)',IOSTAT=ios) LINE
C      IF (LEN_TRIM(LINE).NE.22) RETURN
      READ (551,*) FX,FY,FZ,AMX,AMY,AMZ
      
      If( NGRAV .EQ. 2 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMZ*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+(-AMY*1000D0)
      ENDDO
      Elseif( NGRAV .EQ. 3 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMY*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+AMZ*1000D0
      ENDDO          
      Endif

      ELSEIF (OPTION.EQ."REDD") THEN ! LOAD CASE 1 ( DYNAMIC )
      IF (COUPLEANALYSIS.EQ.2.OR.COUPLEANALYSIS.EQ.3) RETURN    
      IF (IILC.NE.ILC)THEN
      !NumberofPara = "1"
      !NumberofPara    = NFASTPARA(ILC)
      ROOTFST         = '.dat'
      TurbStaFileName = "FASTLoadSta"//"/"//"TowerTopLoads"//"/"//"TurbStaFile"//TRIM(FILENUMBER)//TRIM(ROOTFST)    
      TurbDynFileName = "FASTLoadDyn"//"/"//"TowerTopLoads"//"/"//"TurbDynFile"//TRIM(FILENUMBER)//TRIM(ROOTFST)
      OPEN(UNIT=551,FILE=TurbStaFileName,STATUS='UNKNOWN')
      OPEN(UNIT=552,FILE=TurbDynFileName,STATUS='UNKNOWN')
      
C      IF (ILC.EQ.1)THEN ! FIRST LOOP
C      READ (552,'(A)',IOSTAT=ios) LINE
C      IF (LEN_TRIM(LINE).NE.23) RETURN
C      READ (552,*) 
C      READ (552,*) 
C      ENDIF 
      
      READ (552,*) FX,FY,FZ,AMX,AMY,AMZ
      
C     KN TO N
      if( NGRAV .EQ. 2 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMZ*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+(-AMY*1000D0)
      ENDDO
      Elseif( NGRAV .EQ. 3 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMY*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+AMZ*1000D0
      ENDDO          
      Endif
      IILC = ILC
      ELSEIF (IILC.NE.ILC)THEN
          
      if( NGRAV .EQ. 2 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMZ*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+(-AMY*1000D0)
      ENDDO
      Elseif( NGRAV .EQ. 3 ) then
      DO I = 1,6
      IEQ = ID(I,NODEAF)
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FX*1000D0
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FY*1000D0
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FZ*1000D0
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+AMX*1000D0
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+AMY*1000D0
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+AMZ*1000D0
      ENDDO          
      Endif
      ENDIF
      
      ELSEIF (OPTION.EQ."REDT") THEN ! LOAD CASE 1 ( DYNAMIC )
      IF (COUPLEANALYSIS.EQ.2.OR.COUPLEANALYSIS.EQ.3) RETURN   
      ROOTFST         = '.dat'
      TurbStaFileName = "FASTLoadSta"//"/"//"TowerTopLoads"//"/"//"TurbStaFile"//TRIM(FILENUMBER)//TRIM(ROOTFST)    
      TurbDynFileName = "FASTLoadDyn"//"/"//"TowerTopLoads"//"/"//"TurbDynFile"//TRIM(FILENUMBER)//TRIM(ROOTFST)
      OPEN(UNIT=551,FILE=TurbStaFileName,STATUS='UNKNOWN')
      OPEN(UNIT=552,FILE=TurbDynFileName,STATUS='UNKNOWN')
      REWIND (551)
      REWIND (552)
      
      READ (552,*) 
      READ (552,*) 
      
      ELSEIF (OPTION.EQ."REDA") THEN ! FOR TEST SUBDYN
      OPEN(UNIT=551,FILE="FORCE AND MOMENT FOR HD.TXT",STATUS='UNKNOWN')
      
      READ (551,*) DUMMY,NNODEFAST
      DO J=1,NNODEFAST
      READ (551,*) FFX(J),FFY(J),FFZ(J),FAMX(J),FAMY(J),FAMZ(J)
      ENDDO
      
      DO J=1,NNODEFAST
      DO I = 1,6
      IEQ = ID(I,NODEAFX(J))
      IF (I.EQ.1)  R(IEQ) = R(IEQ)+FFX(NODEAFF(J))
      IF (I.EQ.2)  R(IEQ) = R(IEQ)+FFY(NODEAFF(J))
      IF (I.EQ.3)  R(IEQ) = R(IEQ)+FFZ(NODEAFF(J))
      IF (I.EQ.4)  R(IEQ) = R(IEQ)+FAMX(NODEAFF(J))
      IF (I.EQ.5)  R(IEQ) = R(IEQ)+FAMY(NODEAFF(J))
      IF (I.EQ.6)  R(IEQ) = R(IEQ)+FAMZ(NODEAFF(J))
      ENDDO  
      ENDDO
      
      ELSEIF (OPTION.EQ."REDB") THEN ! FOR TEST SUBDYN
      OPEN(UNIT=552,FILE="NODE ASSIGN.TXT",STATUS='UNKNOWN')
      READ (552,*) NELEMENTFAST
      DO J=1,NELEMENTFAST
      READ (552,*) NODEAFF(J),NODEAFX(J)
      ENDDO
      ENDIF
      
      return

!30    write(*,*) 'Over turbine load numbers'
 !     STOP
      
      END
!     ==================================================================================================
