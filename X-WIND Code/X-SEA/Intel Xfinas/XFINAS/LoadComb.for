C	=====================================================================
C	=====================================================================
C	=====================================================================
!      SUBROUTINE OFFSHORE_LOAD_CASE_NUMBER (ILC,ILOFF,OPT,OPT1)
!	IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER*4 (I-N)
!      CHARACTER*4 OPT,OPT1
!      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
!     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
!      COMMON/OFFSHORE_CASE/ LGEN,IOFFL,IORI_OFFSHORE
!      COMMON/OFFSHORE_CASE_START/ RUNNING_OFFSHORE_LOADCASE
!      COMMON /SOILDYNA/ NSOIL
!      COMMON / SOILRE/ NSOILRETURN
!      ! NUMBER OF STARING OFFSHORE LOAD CASE
!      IORI_OFFSHORE = LGEN+1  
!      
!      IF (ILC.LT.IORI_OFFSHORE)THEN ! NOT INCLUDE OFFSHORE LOAD CASE
!         OPT = "SKIP"
!      ELSEIF (ILC.GE.IORI_OFFSHORE)THEN ! INCLUDING OFFSHORE LOAD CASE
!         OPT = "CALC" 
!             IF (OPT1.EQ."CALL")THEN
!             ILOFF = RUNNING_OFFSHORE_LOADCASE
!             RETURN
!             ENDIF
!         IF (NSOIL.EQ.1.AND.NSOILRETURN.EQ.1)THEN ! PILE-SOIL STRUCTURE INTERACTION
!         RUNNING_OFFSHORE_LOADCASE = RUNNING_OFFSHORE_LOADCASE + 1
!         ILOFF = RUNNING_OFFSHORE_LOADCASE
!         ELSE
!         ILOFF = 1    
!         ENDIF
!         IF (INDEX_OFFSHORE_LOADCASE.GT.IOFFL)THEN
!             ! PROGRAM WILL BE ERROR
!         ENDIF
!      ENDIF
!      
!      END
      ! UPDATE 08/2020
      SUBROUTINE OFFSHORE_LOAD_CASE_NUMBER (ILC,ILOFF,OPT,OPT1)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPT,OPT1
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON/OFFSHORE_CASE/ LGEN,IOFFL,IORI_OFFSHORE
      COMMON/OFFSHORE_CASE_START/ RUNNING_OFFSHORE_LOADCASE
      ! NUMBER OF STARING OFFSHORE LOAD CASE
      IORI_OFFSHORE = LGEN+1  
      
      IF (ILC.LT.IORI_OFFSHORE)THEN ! NOT INCLUDE OFFSHORE LOAD CASE
         OPT = "SKIP"
      ELSEIF (ILC.GE.IORI_OFFSHORE)THEN ! INCLUDING OFFSHORE LOAD CASE
         OPT = "CALC" 
             IF (OPT1.EQ."CALL")THEN
             ILOFF = RUNNING_OFFSHORE_LOADCASE
             RETURN
             ENDIF
         RUNNING_OFFSHORE_LOADCASE = RUNNING_OFFSHORE_LOADCASE + 1
         ILOFF = RUNNING_OFFSHORE_LOADCASE
         IF (INDEX_OFFSHORE_LOADCASE.GT.IOFFL)THEN
             ! PROGRAM WILL BE ERROR
         ENDIF
      ENDIF
      END
C	=====================================================================
C	=====================================================================
      SUBROUTINE LDCMREAD
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON/OFFSHORE_CASE/ LGEN,IOFFL,IORI_OFFSHORE
C     ----------------------------------------------------------------
	CALL INTFILL('NOUT',LCSINP,1,2,0)   
	REWIND(LCSINP)

	READ(ITI,*) 
	READ(ITI,*) NUMCOM
	WRITE(LCSINP) NUMCOM

C	READ HEADER
	CALL LCHEADER(NUMCOM,'COMB')

	DO ICM = 1,NUMCOM
	READ(ITI,*) 
	READ(ITI,*) NCASE,IOPT
	WRITE(LCSINP) NCASE,IOPT
	DO ICS = 1,NCASE
	READ(ITI,*) ILCS,ILOFF,LTYP,FACT
      ! OFFSHORE LOAD COMBINATION
      IF (ILOFF.EQ.0)THEN
	WRITE(LCSINP) ILCS,LTYP,FACT
      ELSEIF (ILOFF.NE.0)THEN
          IF (LTYP.NE.1)THEN
          WRITE(LCSINP) LGEN+ILOFF,LTYP,FACT
          ELSEIF (LTYP.EQ.1)THEN
          WRITE(LCSINP) ILCS,LTYP,FACT
          ENDIF
      ENDIF
	ENDDO
	ENDDO

	
	
      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LDCMCALC
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------- 
	CALL INTFILL('NOUT',LCSOUT,1,1,0) 
	CALL INTFILL('NOUT',LCSINP,1,2,0)   
	CALL INTFILL('NOUT',LCSCOM,1,4,0)  

	REWIND(LCSINP)

	READ(LCSINP) NUMCOM

	DO 1000 IUMCOM = 1,NUMCOM

	LCSFIL = IUMCOM  
	FACTOR = 0.0D0
	CALL INIOPER(LCSCOM,LCSFIL,FACTOR,'CLEF','OLDF')   

	READ(LCSINP) NCASE,IOPT

	DO 900 ICASE = 1,NCASE

	READ(LCSINP) ILCS,LTYP,FACT
	IF(LTYP.EQ.0) LP = LCSOUT 
	IF(LTYP.EQ.1) LP = LCSCOM  

	CALL CLEROUT
	FACTOR = 1.0D0
	CALL INIOPER(LP,ILCS,FACTOR,'CALL','OLDF')  

	SELECTCASE(IOPT)
	CASE(0)
	FACTOR = FACT
	CALL INIOPER(LCSCOM,LCSFIL,FACTOR,'COMB','OLDF')  
	CASE(1)
	FACTOR = FACT
	CALL INIOPER(LCSCOM,LCSFIL,FACTOR,'SORT','OLDF')  !LOAD COMBINATION SORT
	ENDSELECT

900	CONTINUE


1000	CONTINUE


C	PRINT OUT ALL OF COMBINATION STEP
	CALL PRNFLAG('ELEM','LINK','GSUP','LSUP','DISP','GSPG','LSPG')
	CALL  PRNOUT('COMB','PALL','CALL',NUMCOM)
	


	
      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
