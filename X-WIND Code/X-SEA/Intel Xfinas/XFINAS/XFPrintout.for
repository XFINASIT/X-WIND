C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNFLAG(ELEM,LINK,GSUP,LSUP,DISP,GSPG,LSPG)  !PRINT OUTPUT FLAG
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 ELEM,LINK,GSUP,LSUP,DISP,GSPG,LSPG

	CALL INTZERO('NPRN')

	IF(ELEM.EQ.'ELEM') CALL INTFILL('NPRN',1,1,1,1)	
	IF(LINK.EQ.'LINK') CALL INTFILL('NPRN',1,1,2,1)	
	IF(GSUP.EQ.'GSUP') CALL INTFILL('NPRN',1,1,3,1)	
	IF(LSUP.EQ.'LSUP') CALL INTFILL('NPRN',1,1,4,1)	
	IF(DISP.EQ.'DISP') CALL INTFILL('NPRN',1,1,5,1)	
	IF(GSPG.EQ.'GSPG') CALL INTFILL('NPRN',1,1,6,1)	
	IF(LSPG.EQ.'LSPG') CALL INTFILL('NPRN',1,1,7,1)		

	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT(METH,OPER,FILE,NSTEP)  !PRINT ALL OUTPUT
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 METH,OPER,FILE

      COMMON /FTIM/ TIM(20),IDATE,ITIME
      
C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      
            COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,                      
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,                  
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
            
      CALL CPU_TIME (TIM1)
      
C     FIRST PRINT TO GIDOUT.FLAVIA.RES FOR GID OUTPUT      
      CALL INTFILL('%IOL',IGO,1,4 ,0)  !GIDOUT.FLAVIA
      LFPRN1 = IGO
      LFPRN2 = IGO
      LFPRN3 = IGO
      LFPRN4 = IGO
      KEXCEL = 0
      CALL PRNOUTP(METH,OPER,FILE,NSTEP)
       
C      NEXT PRINT TO OUTPUT FILE FOR EXCEL REPORT 
C      LFPRN1 = 106 !Displacement
C      LFPRN2 = 107 !Element
C      LFPRN3 = 108 !Support
C      LFPRN4 = 109 !Link
C      KEXCEL = 1
C      CALL PRNOUTP(METH,OPER,FILE,NSTEP)
           
      CALL CPU_TIME (TIM2)
      TIM(17) = TIM(17) + (TIM2-TIM1)
      
      RETURN
      END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUTP(METH,OPER,FILE,NSTEP)  !PRINT ALL OUTPUT
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 METH,OPER,FILE

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER
      
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,                      
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,                  
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON/WAVE_DISTRIBUTED/ ICASE
          
      ICASE = 1D0
	INM = -1
	IF(METH.EQ.'STND') INM = 0 
	IF(METH.EQ.'COMB') INM = 1 
	SELECTCASE(INM)
	CASE(0)
	CALL INTFILL('NOUT',NFIL,1,1,0)  !LOAD CASE OUTPUT FILE
	CASE(1)
	CALL INTFILL('NOUT',NFIL,1,4,0)  !LOAD COMBINATION OUTPUT FILE
	ENDSELECT
	CALL INTFILL('NPRN',INM,1,12,1)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)

	IND = -1
	IF(OPER.EQ.'PALL') IND = 0 
	IF(OPER.EQ.'PONE') IND = 1 
	CALL INTFILL('NPRN',IND,1,13,1)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
	SELECTCASE(IND)
	CASE(0)
	IFIRST = 1
	CASE(1)
	IFIRST = NSTEP
	ENDSELECT


	LLC = 4
	CALL INTFILL('NPRN',LLC,1,11,1)  !PRINT TOT OR POS OR NEG


	DO ISTEP = IFIRST,NSTEP

	IF(FILE.EQ.'CALL') THEN
	FACT = 1.0D0
	CALL INIOPER(NFIL,ISTEP,FACT,'CALL','OLDF')
	ENDIF

	CALL INTFILL('NPRN',LPRN,1,1,0)	
	LFPRIN = LFPRN2	
	IF(LPRN.NE.0) CALL PRNOUT1(ISTEP)  !ELEMENT DATA
	
	CALL INTFILL('NPRN',LPRN,1,2,0)	
	LFPRIN = LFPRN4		
	IF(LPRN.NE.0) CALL PRNOUT2(ISTEP)  !LINK ELEMENT DATA
	
	CALL INTFILL('NPRN',LPRN,1,3,0)
	LFPRIN = LFPRN3			
	IF(LPRN.NE.0) CALL PRNOUT3(ISTEP)  !GLOBAL SUPPORT DATA
      
      IF(KOFFL.EQ.1.AND.KSPEC.EQ.0)THEN
      IF(LPRN.NE.0) CALL PRNOUTWAVEFORCE(ISTEP) !OFFSHORE WAVE FORCE
      ENDIF
	
	CALL INTFILL('NPRN',LPRN,1,4,0)	
	LFPRIN = LFPRN3		
	IF(LPRN.NE.0) CALL PRNOUT4(ISTEP)  !LOCAL SUPPORT DATA
	
	CALL INTFILL('NPRN',LPRN,1,5,0)
	LFPRIN = LFPRN1			
	IF(LPRN.NE.0) CALL PRNOUT5(ISTEP)  !FREE NODE DATA
	
	CALL INTFILL('NPRN',LPRN,1,6,0)	
	LFPRIN = LFPRN3			
	IF(LPRN.NE.0) CALL PRNOUT6(ISTEP)  !GLOBAL SPRING DATA
	
	CALL INTFILL('NPRN',LPRN,1,7,0)	
	LFPRIN = LFPRN3			
	IF(LPRN.NE.0) CALL PRNOUT7(ISTEP)  !LOCAL SPRING DATA
	ENDDO


      RETURN
      END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT1(ISTEP)  !ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD
	COMMON /OSPEC/ ISTIME,IOUTSPEC,NFATIGUE
      COMMON /OUTFA/ CHANA,ICHANA
C     ----------------------------------------------------------------
	DIMENSION IDTEST(2,200)
	
      CHANA  = 0
      ICHANA = 1
      
C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 

C	-------------------------
	NUM = 0
	IDTEST(1:2,1) = 0
	DO 100 IGM = 1,NGIDM

	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)  !CALL IEG
	IF(IEG.EQ.0) GOTO 100
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)  !CALL IEL
	CALL INTFILL('OGRP',ITYPE,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTYP,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0)  !NNM

	ITT = 1
	DO II = 1,NUM
	ITEST1 = IDTEST(1,II)
	ITEST2 = IDTEST(2,II)
	IF(ITEST1.EQ.ITYPE.AND.ITEST2.EQ.ISTYP) ITT = 0
	ENDDO
	IF(ITT.EQ.1) THEN
	NUM = NUM + 1
	IDTEST(1,NUM) = ITYPE
	IDTEST(2,NUM) = ISTYP
	ENDIF

100   CONTINUE
C	-------------------------

	NUMP = NUM
	DO II = 1,NUMP
	ITYPE = IDTEST(1,II)
	ISTYP = IDTEST(2,II)

	DO IP = 1,2

	SELECTCASE(IP)
	CASE(1) 
	METHOD = 'MAX'
	CASE(2)
	METHOD = 'MIN'
	ENDSELECT

	SELECTCASE(ITYPE)
	CASE(2)
	CALL PRNTYP2(ISTYP,ISTEP,METHOD)
	CASE(5)
	CALL PRNTYP5(ISTYP,ISTEP,METHOD)
	CASE(6)
	CALL PRNTYP6(ISTYP,ISTEP,METHOD)
	CALL PRNTYP6N(ISTYP,ISTEP,METHOD)
	CASE(8)
	CALL PRNTYP8(ISTYP,ISTEP,METHOD)
	CALL PRNTYP8N(ISTYP,ISTEP,METHOD)
      CASE(9) 
          KAK = 1
	CALL PRNTYP9(ISTYP,ISTEP,METHOD)
	CALL PRNTYP9N(ISTYP,ISTEP,METHOD)
      CASE(10)
      	IF(ISTYP.EQ.4) THEN ! FOR SOLID-SHELL(LIU) BY BJ
	      !FOR RESULTANT STRESS
            CALL PRNTYP104(ISTYP,ISTEP,METHOD)
	      CALL PRNTYP104N(ISTYP,ISTEP,METHOD)
            !FOR STRESS
	      CALL PRNTYP104_ST(ISTYP,ISTEP,METHOD)
	      CALL PRNTYP104N_ST(ISTYP,ISTEP,METHOD)            
	  ELSE                ! FOR SOLID
	      CALL PRNTYP10(ISTYP,ISTEP,METHOD)
	      CALL PRNTYP10N(ISTYP,ISTEP,METHOD)
	  ENDIF                 
	CASE(11)
	CALL PRNTYP11(ISTYP,ISTEP,METHOD)
	CALL PRNTYP11N(ISTYP,ISTEP,METHOD)
	CASE(13)
	CALL PRNTYP13(ISTYP,ISTEP,METHOD,NNM)
	CASE(14)
	CALL PRNTYP14(ISTYP,ISTEP,METHOD)
	CASE(15)
	CALL PRNTYP15(ISTYP,ISTEP,METHOD,NNM)
	CASE(16)
	CALL PRNTYP16(ISTYP,ISTEP,METHOD)
	ENDSELECT

	ENDDO

	ENDDO

C
	CALL PRNTYP17(ISTEP)

C
      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT2(ISTEP)  !LINK ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO
C     ----------------------------------------------------------------
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
C     ----------------------------------------------------------------      
	ALLOCATABLE EDAT(:)

	IF(IFLOOR.EQ.3) THEN
	  NLINK = 0
	ELSE
	  CALL LOCATN('OLNK',KLNK,NOUT2,NLINK,2)  
      ENDIF
      
	IF(NLINK.LE.0) RETURN

	NOUT = NOUT2/2
	ALLOCATE(EDAT(NOUT2))

	DO IC = 1,2

	SELECTCASE(IC)
	CASE(1) 
	COMPO = 'FOC'
	CASE(2)
	COMPO = 'MOM'
	ENDSELECT


	DO IP = 1,2

	SELECTCASE(IP)
	CASE(1) 
	METHOD = 'MAX'
	CASE(2)
	METHOD = 'MIN'
	ENDSELECT

	CALL PRNHEAD2(ISTEP,METHOD,COMPO)
	DO IGM = 1,NLINK
	DO I = 1,NOUT2
	CALL RELFILL('OLNK',EDAT(I),I,IGM,0)  
	ENDDO
	CALL PRNVALV2(EDAT,IGM,NOUT,METHOD,COMPO)
	ENDDO
	CALL PRNEND(METHOD)

	ENDDO

	ENDDO


	DEALLOCATE(EDAT)


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT3(ISTEP)  !GLOBAL SUPPORT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO
C     OUTPUT INDEX      
      COMMON /OUTPUTCONTROL/ NLFORCE,NLSTRESS,NSTRESSP,NTORSION,NREATION,NDISPROTAWQ
C     ----------------------------------------------------------------
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
	DIMENSION IDOF(9),IPRN(4)

	ALLOCATABLE EDAT(:)

	CALL LOCATN('OGSP',KSUPG,NOUT2,NSN,2)  
	IF(NSN.LE.0) RETURN
      
C     FOR LINK CONSTRAINTS BY BJ	
	IF(IFLOOR.EQ.3) NSN = NSN0 

	NOUT = NOUT2/2
	ALLOCATE(EDAT(NOUT2))

	IDOF(1:9) = 0
	DO I = 1,9
	CALL INTFILL('%DOB',IK,1,I,0)
	IF(IK.EQ.0) IDOF(I) = 1
	ENDDO
	IPRN(1:4) = 0
	IF(IDOF(1).NE.0.OR.IDOF(2).NE.0.OR.IDOF(3).NE.0) IPRN(1) = 1
	IF(IDOF(4).NE.0.OR.IDOF(5).NE.0.OR.IDOF(6).NE.0) IPRN(2) = 1
	IF(IDOF(8).NE.0) IPRN(3) = 1
	IF(IDOF(9).NE.0) IPRN(4) = 1



	DO IC = 1,4


	IF(IPRN(IC).NE.0) THEN

	SELECTCASE(IC)
	CASE(1) 
	COMPO = 'FOC'
	CASE(2)
	COMPO = 'MOM'
	CASE(3)
	COMPO = 'TEM'
	CASE(4)
	COMPO = 'PRS'
	ENDSELECT

	DO IP = 1,2

	SELECTCASE(IP)
	CASE(1) 
	METHOD = 'MAX'
	CASE(2)
	METHOD = 'MIN'
      ENDSELECT
      
      IF (NREATION.NE.1) RETURN
	CALL PRNHEAD3(ISTEP,METHOD,COMPO)
	DO ISN = 1,NSN
	DO I = 1,NOUT2
	CALL RELFILL('OGSP',EDAT(I),I,ISN,0)  
	ENDDO
	CALL PRNVALV3(EDAT,ISN,NOUT,METHOD,COMPO)
	ENDDO
	CALL PRNEND(METHOD)

	ENDDO

	ENDIF

	ENDDO


	DEALLOCATE(EDAT)


C
      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT4(ISTEP)  !LOCAL SUPPORT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO
C     ----------------------------------------------------------------
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
C     ----------------------------------------------------------------
	DIMENSION IDOF(9),IPRN(4)
	ALLOCATABLE EDAT(:)


	CALL INTFILL('%NUB',NLS  ,1,11,0)  !NUMBER OF LOCAL SUPPORT
	IF(NLS.LE.0) RETURN
      
C     FOR LINK CONSTRAINTS BY BJ	
	IF(IFLOOR.EQ.3) NSN = NSN0 	
      
	CALL LOCATN('OLSP',KSUPL,NOUT2,NSN,2)  
	IF(NSN.LE.0) RETURN

	NOUT = NOUT2/2
	ALLOCATE(EDAT(NOUT2))


	IDOF(1:9) = 0
	DO I = 1,9
	CALL INTFILL('%DOB',IK,1,I,0)
	IF(IK.EQ.0) IDOF(I) = 1
	ENDDO
	IPRN(1:4) = 0
	IF(IDOF(1).NE.0.OR.IDOF(2).NE.0.OR.IDOF(3).NE.0) IPRN(1) = 1
	IF(IDOF(4).NE.0.OR.IDOF(5).NE.0.OR.IDOF(6).NE.0) IPRN(2) = 1
	IF(IDOF(8).NE.0) IPRN(3) = 1
	IF(IDOF(9).NE.0) IPRN(4) = 1


	DO IC = 1,4

	IF(IPRN(IC).NE.0) THEN

	SELECTCASE(IC)
	CASE(1) 
	COMPO = 'FOC'
	CASE(2)
	COMPO = 'MOM'
	CASE(3)
	COMPO = 'TEM'
	CASE(4)
	COMPO = 'PRS'
	ENDSELECT

	DO IP = 1,2

	SELECTCASE(IP)
	CASE(1) 
	METHOD = 'MAX'
	CASE(2)
	METHOD = 'MIN'
	ENDSELECT


	CALL PRNHEAD4(ISTEP,METHOD,COMPO)
	DO ISN = 1,NSN
	DO I = 1,NOUT2
	CALL RELFILL('OLSP',EDAT(I),I,ISN,0)  
      ENDDO
      
	CALL PRNVALV4(EDAT,ISN,NOUT,METHOD,COMPO)

      
	ENDDO
	CALL PRNEND(METHOD)

	ENDDO

	ENDIF

	ENDDO


	DEALLOCATE(EDAT)


C
      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT5(ISTEP)  !FREE NODE DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO
C     OUTPUT INDEX      
      COMMON /OUTPUTCONTROL/ NLFORCE,NLSTRESS,NSTRESSP,NTORSION,NREATION,NDISPROTAWQ
C     ----------------------------------------------------------------
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
C     ----------------------------------------------------------------
	DIMENSION IDOF(9),IPRN(5)
	ALLOCATABLE EDAT(:)

	CALL LOCATN('ONDS',KDISP,NOUT2,NSN,2)  
	IF(NSN.LE.0) RETURN
      
C     FOR LINK CONSTRAINTS BY BJ	
	IF(IFLOOR.EQ.3) NSN = NSN0 	
      
	NOUT = NOUT2/2
	ALLOCATE(EDAT(NOUT2))

	IDOF(1:9) = 0
	DO I = 1,9
	CALL INTFILL('%DOB',IK,1,I,0)
	IF(IK.EQ.0) IDOF(I) = 1
	ENDDO
	IPRN(1:5) = 0
	IF(IDOF(1).NE.0.OR.IDOF(2).NE.0.OR.IDOF(3).NE.0) IPRN(1) = 1
	IF(IDOF(4).NE.0.OR.IDOF(5).NE.0.OR.IDOF(6).NE.0) IPRN(2) = 1
	IF(IDOF(8).NE.0) IPRN(3) = 1
	IF(IDOF(9).NE.0) IPRN(4) = 1
	IF(IDOF(7).NE.0) IPRN(5) = 1

	DO IC = 1,5


	IF(IPRN(IC).NE.0) THEN

	SELECTCASE(IC)
	CASE(1) 
	COMPO = 'DIS'
	CASE(2)
	COMPO = 'ROT'
	CASE(3)
	COMPO = 'TEM'
	CASE(4)
	COMPO = 'PRS'
	CASE(5)
	COMPO = 'WRP'
	ENDSELECT

	DO IP = 1,2

	SELECTCASE(IP)
	CASE(1) 
	METHOD = 'MAX'
	CASE(2)
	METHOD = 'MIN'
	ENDSELECT

      IF (NDISPROTAWQ.NE.1) RETURN
	CALL PRNHEAD5(ISTEP,METHOD,COMPO)
	DO ISN = 1,NSN
	DO I = 1,NOUT2
	CALL RELFILL('ONDS',EDAT(I),I,ISN,0)  
	ENDDO
	CALL PRNVALV5(EDAT,ISN,NOUT,METHOD,COMPO)
	ENDDO
	CALL PRNEND(METHOD)

	ENDDO

	ENDIF

	ENDDO


	DEALLOCATE(EDAT)



      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT6(ISTEP)  !GLOBAL SPRING DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO
C     ----------------------------------------------------------------
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
C     ----------------------------------------------------------------
	ALLOCATABLE EDAT(:)

	CALL LOCATN('OGSS',KGSS,NOUT2,NSN,2)  
	IF(NSN.LE.0) RETURN

C     FOR LINK CONSTRAINTS BY BJ	
	IF(IFLOOR.EQ.3) NSN = NSN0 	
      
	NOUT = NOUT2/2
	ALLOCATE(EDAT(NOUT2))


	DO IC = 1,2


	SELECTCASE(IC)
	CASE(1) 
	COMPO = 'FOC'
	CASE(2)
	COMPO = 'MOM'
	ENDSELECT

	DO IP = 1,2

	SELECTCASE(IP)
	CASE(1) 
	METHOD = 'MAX'
	CASE(2)
	METHOD = 'MIN'
	ENDSELECT


	CALL PRNHEAD6(ISTEP,METHOD,COMPO)
	DO ISN = 1,NSN
	DO I = 1,NOUT2
	CALL RELFILL('OGSS',EDAT(I),I,ISN,0)  
	ENDDO
	CALL PRNVALV6(EDAT,ISN,NOUT,METHOD,COMPO)
	ENDDO
	CALL PRNEND(METHOD)

	ENDDO


	ENDDO


	DEALLOCATE(EDAT)



      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNOUT7(ISTEP)  !LOCAL SPRING DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO
C     ----------------------------------------------------------------
C     FOR LINK CONSTRAINTS BY BJ
      COMMON /FLOOR/ IFLOOR,NFLOOR,NSN0,NSNL
C     ----------------------------------------------------------------
	ALLOCATABLE EDAT(:)

	CALL LOCATN('OLSS',KLSS,NOUT2,NSN,2)  
	IF(NSN.LE.0) RETURN

C     FOR LINK CONSTRAINTS BY BJ	
	IF(IFLOOR.EQ.3) NSN = NSN0 	
      
	NOUT = NOUT2/2
	ALLOCATE(EDAT(NOUT2))


	DO IC = 1,2


	SELECTCASE(IC)
	CASE(1) 
	COMPO = 'FOC'
	CASE(2)
	COMPO = 'MOM'
	ENDSELECT

	DO IP = 1,2

	SELECTCASE(IP)
	CASE(1) 
	METHOD = 'MAX'
	CASE(2)
	METHOD = 'MIN'
	ENDSELECT


	CALL PRNHEAD7(ISTEP,METHOD,COMPO)
	DO ISN = 1,NSN
	DO I = 1,NOUT2
	CALL RELFILL('OLSS',EDAT(I),I,ISN,0)  
	ENDDO
	CALL PRNVALV7(EDAT,ISN,NOUT,METHOD,COMPO)
	ENDDO
	CALL PRNEND(METHOD)

	ENDDO


	ENDDO


	DEALLOCATE(EDAT)


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNGAUS(ITYPE,ISTYP,NPT,LGASPN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      
	COMMON /BF_SMOTH/ NP_SMH
C     ----------------------------------------------------------------
	DIMENSION LGASPN(20,20)

      CALL INTFILL('%IOL',IGO,1,4 ,0)  !GIDOUT.FLAVIA
      LFPRIN = IGO
      
	SELECTCASE(ITYPE)

	CASE(2)
	IF(ISTYP.EQ.3) THEN
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN
	
		
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,201) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,211)   !EXCEL  
	
	LGASPN(ITYPE,ISTYP) = 1
	ENDIF

	IF(ISTYP.EQ.4) THEN
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,202) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,212)   !EXCEL    
	
	LGASPN(ITYPE,ISTYP) = 1
	ENDIF

	IF(ISTYP.EQ.7) THEN
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,203) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,213)   !EXCEL  
	
	LGASPN(ITYPE,ISTYP) = 1
	ENDIF

	CASE(5)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,501) NP_SMH+2
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,511) NP_SMH+2  !EXCEL  
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,502) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,512)   !EXCEL  
      
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,503) 
	
	LGASPN(ITYPE,ISTYP) = 1

	CASE(6)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	SELECTCASE(ISTYP)
	CASE(0,1,2)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,601) NPT
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,611) NPT  !EXCEL  
	CASE(3,4,5)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,602) NPT
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,612) NPT  !EXCEL  
	ENDSELECT
	LGASPN(ITYPE,ISTYP) = 1
	RETURN

	CASE(8)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,801) NPT
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,811) NPT  !EXCEL  
	LGASPN(ITYPE,ISTYP) = 1
	RETURN

	CASE(9)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	SELECTCASE(ISTYP)
	CASE(11)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,901) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,911)   !EXCEL   
	CASE(1)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,902) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,912)   !EXCEL   
	CASE(4)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,903) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,913)   !EXCEL   
	CASE(12)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,904) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,914)   !EXCEL  	
	CASE(3)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,905) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,915)   !EXCEL  	
      CASE(13) ! FOR SHEAR WALL BY BJ
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,906)
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,916)   !EXCEL   
	ENDSELECT
	LGASPN(ITYPE,ISTYP) = 1
	RETURN

	CASE(10)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	SELECTCASE(ISTYP)
	CASE(1) ! FOR SOLID-SHELL(P) BY BJ
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1010) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1020)   !EXCEL  	          
	CASE(4) ! FOR SOLID-SHELL(LIU) BY BJ
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1011) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1021)   !EXCEL  	
	CASE(6)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1012) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1212)   !EXCEL  		
	CASE(13)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1013) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1213)   !EXCEL  		
	ENDSELECT
	RETURN

	CASE(11)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1111) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1211)   !EXCEL  		
	LGASPN(ITYPE,ISTYP) = 1
	RETURN

	CASE(13)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	SELECTCASE(ISTYP)
	CASE(3)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1301) NPT ! 4 NODE
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1311) NPT  !EXCEL  
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1302) NPT ! 8 NODE
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1312) NPT  !EXCEL  
	CASE(4)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1303)  ! 3 NODE
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1313)    !EXCEL  
	ENDSELECT
	RETURN

	CASE(14)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	SELECTCASE(ISTYP)
	CASE(3)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1401) ! 8 NODE HEX
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1411)   !EXCEL  
	CASE(4)
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,1402) ! 4 NODE TET
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1412)   !EXCEL  
	ENDSELECT
	RETURN

	CASE(15)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	SELECTCASE(ISTYP)
	CASE(3)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1501) NPT ! 4 NODE 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1511) NPT  !EXCEL  
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1502) NPT ! 8 NODE 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1512) NPT  !EXCEL  
	CASE(4)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1503)  ! 3 NODE 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1513)   !EXCEL  
	ENDSELECT
	RETURN

	CASE(16)
	IF(LGASPN(ITYPE,ISTYP).EQ.1) RETURN	
	SELECTCASE(ISTYP)
	CASE(3)
	WRITE (LFPRIN,1601)   ! 8 NODE HEX
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1601)  ! 8 NODE HEX
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1611)   !EXCEL  
	CASE(4)
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,1602)  ! 8 NODE TET
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,1612)   !EXCEL  
	ENDSELECT
	RETURN

	ENDSELECT


201	FORMAT('GaussPoints "XTruss" ElemType Linear'/,
	2'Number of Gauss Points: 1'/,
	3'Natural Coordinates: Internal'/,
	4'Nodes Not Included',/,
	5'End GaussPoints',/)

211	FORMAT('GaussPoints "XTruss" ElemType Linear'/,   !EXCEL  
	2"Number of Gauss Points:,","1,",/,
	3"Natural Coordinates:,","Internal,",/,
	4'Nodes Not Included',/,
	5'End GaussPoints',/)

202	FORMAT('GaussPoints "XCable" ElemType Linear'/,
	2'Number of Gauss Points: 1'/,
	3'Natural Coordinates: Internal'/,
	4'Nodes Not Included',/,
	5'End GaussPoints',/)

212	FORMAT('GaussPoints "XCable" ElemType Linear'/,   !EXCEL  
	2"Number of Gauss Points:,","1",/,
	3"Natural Coordinates:,","Internal,",/,
	4'Nodes Not Included',/,
	5'End GaussPoints',/)	

203   FORMAT('GaussPoints "XCableC" ElemType Linear'/,
	2'Number of Gauss Points: 1'/,
	3'Natural Coordinates: Internal'/,
	4'Nodes Not Included',/,
	5'End GaussPoints',/)
      
213	FORMAT('GaussPoints "XCableC" ElemType Linear'/,   !EXCEL  
	2"Number of Gauss Points:,","1,",/,
	3"Natural Coordinates:,","Internal,"/,
	4'Nodes Not Included',/,
	5'End GaussPoints',/)

501	FORMAT ('GaussPoints "XFrame" ElemType Linear'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

511	FORMAT ('GaussPoints "XFrame" ElemType Linear'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)	

502	FORMAT ('GaussPoints "XFrame2" ElemType Linear'/
     2        'Number of Gauss Points: 2'/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

512   FORMAT ('GaussPoints "XFrame2" ElemType Linear'/   !EXCEL  
     2        "Number of Gauss Points:,","2,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)	
      
503	FORMAT('GaussPoints "XWave" ElemType Linear'/,
	2        'Number of Gauss Points: 2'/,
	3        'Natural Coordinates: Internal'/,
	4        'Nodes Included',/,
	5        'End GaussPoints',/)

601	FORMAT ('GaussPoints "XPlane8" ElemType Quadrilateral'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

611	FORMAT ('GaussPoints "XPlane8" ElemType Quadrilateral'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)	

602	FORMAT ('GaussPoints "XPlane4" ElemType Quadrilateral'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

612	FORMAT ('GaussPoints "XPlane4" ElemType Quadrilateral'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)


801	FORMAT ('GaussPoints "XPlaneC" ElemType Quadrilateral'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

811	FORMAT ('GaussPoints "XPlaneC" ElemType Quadrilateral'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)	


901	FORMAT( 'GaussPoints "XShell3Q" Elemtype Triangle'/,
	1	    'Number of Gauss Points: 3'/,
	2        'Natural Coordinates: Given'/
	3        '   0.0  0.0'/
	4        '   1.0  0.0'/
	5        '   0.0  1.0'/
	7         'End GaussPoints'/)

911	FORMAT( 'GaussPoints "XShell3Q" Elemtype Triangle'/,   !EXCEL  
	1	    "Number of Gauss Points:,","3,",/,
	2      "Natural Coordinates:,","Given,",/
	3         "0.0,","0.0,",/
	4         "1.0,","0.0,",/
	5         "0.0,","1.0,",/
	7         'End GaussPoints'/)	

902	FORMAT('GaussPoints "XShell4" Elemtype Quadrilateral'/,
	1	   'Number of Gauss Points: 4'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints'/)

912   FORMAT('GaussPoints "XShell4" Elemtype Quadrilateral'/,   !EXCEL  
	1	   "Number of Gauss Points:,","4,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints'/)	
      
906   FORMAT('GaussPoints "XWALL4" Elemtype Quadrilateral'/,
	1	   'Number of Gauss Points: 4'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints'/)
      
916	FORMAT('GaussPoints "XWALL4" Elemtype Quadrilateral'/,   !EXCEL  
	1	   "Number of Gauss Points:,","4,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints'/)	


903	FORMAT ('GaussPoints "XShell4Q" ElemType Quadrilateral'/
     1        'Number of Gauss Points: 4'/
	2        'Natural Coordinates: Given'/
	3        '  -1.0 -1.0'/
	4        '   1.0 -1.0'/
	5        '   1.0  1.0'/
	6        '  -1.0  1.0'/
	7         'End GaussPoints'/)

913	FORMAT ('GaussPoints "XShell4Q" ElemType Quadrilateral'/   !EXCEL  
     1        "Number of Gauss Points:,","4,",/
	2        "Natural Coordinates:,","Given,",/
	3        "-1.0,","-1.0,",/
	4        "1.0,","-1.0,",/
	5        "1.0,","1.0,",/
	6        "-1.0,","1.0,",/
	7         'End GaussPoints'/)

904	FORMAT('GaussPoints "XShell9" Elemtype Quadrilateral'/,
	1	   'Number of Gauss Points: 9'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints'/)

914	FORMAT('GaussPoints "XShell9" Elemtype Quadrilateral'/,   !EXCEL  
	1	   "Number of Gauss Points:,","9,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints'/)	

905	FORMAT('GaussPoints "XShell8" Elemtype Quadrilateral'/,
	1	   'Number of Gauss Points: 9'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints'/)

915	FORMAT('GaussPoints "XShell8" Elemtype Quadrilateral'/,   !EXCEL  
	1	   "Number of Gauss Points:,","9,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints'/)	


cc1011	FORMAT('GaussPoints "XSolidANS" Elemtype Hexahedra'/,
cc     1        'Number of Gauss Points: 4'/
cc	2        'Natural Coordinates: Given'/
cc	3        '  -1.0 -1.0  0.0'/
cc	4        '   1.0 -1.0  0.0'/
cc	5        '   1.0  1.0  0.0'/
cc	6        '  -1.0  1.0  0.0'/
cc	7         'End GaussPoints'/)

CC1021	FORMAT('GaussPoints "XSolidANS" Elemtype Hexahedra'/,   !EXCEL  
CC     1        "Number of Gauss Points:,","4,",/
CC	2        "Natural Coordinates:,","Given,",/
CC	3        "-1.0,","-1.0,","0.0,",/
CC	4        "1.0,","-1.0,","0.0,",/
CC	5        "1.0,","1.0,","0.0,",/
CC	6        "-1.0,","1.0,","0.0,",/
CC	7         'End GaussPoints'/)	

1010	FORMAT('GaussPoints "XSolidSHL" Elemtype Hexahedra'/, !SOLID-SHELL(P) BY BJ
	1	   'Number of Gauss Points: 8'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1020	FORMAT('GaussPoints "XSolidSHL" Elemtype Hexahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","8,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')	

1011	FORMAT ('GaussPoints "XSolidSHL" ElemType Hexahedra'/ !SOLID-SHELL(LIU) BY BJ
     1        'Number of Gauss Points: 4'/
	2        'Natural Coordinates: Given'/
	3        '  -0.57735 0.0 -0.57735'/
	4        '   0.57735 0.0 -0.57735'/
	5        '   0.57735 0.0  0.57735'/
	6        '  -0.57735 0.0  0.57735'/
	7         'End GaussPoints'/)
	

1021	FORMAT('GaussPoints "XSolidSHL" Elemtype Hexahedra'/,   !EXCEL  
     1        "Number of Gauss Points:,","4,",/
	2        "Natural Coordinates:,","Given,",/
	3        "-1.0,","-1.0,","0.0,",/
	4        "1.0,","-1.0,","0.0,",/
	5        "1.0,","1.0,","0.0,",/
	6        "-1.0,","1.0,","0.0,",/
	7         'End GaussPoints'/)		
      

1012	FORMAT('GaussPoints "XSolidEAS" Elemtype Hexahedra'/,
	1	   'Number of Gauss Points: 8'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1212	FORMAT('GaussPoints "XSolidEAS" Elemtype Hexahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","8,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')	

1013	FORMAT('GaussPoints "XSolidTeT" Elemtype Tetrahedra'/,
	1	   'Number of Gauss Points: 4'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1213	FORMAT('GaussPoints "XSolidTeT" Elemtype Tetrahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","4,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')	


1111	FORMAT('GaussPoints "XSolidCon" Elemtype Hexahedra'/,
	1	   'Number of Gauss Points: 8'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1211	FORMAT('GaussPoints "XSolidCon" Elemtype Hexahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","8,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')	
	

1301	FORMAT ('GaussPoints "XSeepage4" ElemType Quadrilateral'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)


1311	FORMAT ('GaussPoints "XSeepage4" ElemType Quadrilateral'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)
		

1302	FORMAT ('GaussPoints "XSeepage8" ElemType Quadrilateral'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)


1312  FORMAT ('GaussPoints "XSeepage8" ElemType Quadrilateral'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)
		

1303	FORMAT ('GaussPoints "XSeepage3" ElemType Triangle'/
	1	   'Number of Gauss Points: 3'/,
	2        'Natural Coordinates: Given'/
	3        '   1.0  0.0'/
	4        '   0.0  0.0'/
	5        '   0.0  1.0'/
	7         'End GaussPoints'/)


1313	FORMAT ('GaussPoints "XSeepage3" ElemType Triangle'/   !EXCEL  
	1	    "Number of Gauss Points:,","3,",/,
	2      "Natural Coordinates:,","Given,",/
	3        "1.0,","0.0,",/
	4        "0.0,","0.0,",/
	5        "0.0,","1.0,",/
	7         'End GaussPoints'/)	

1401	FORMAT('GaussPoints "XSeepage8H" Elemtype Hexahedra'/,
	1	   'Number of Gauss Points: 8'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1411	FORMAT('GaussPoints "XSeepage8H" Elemtype Hexahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","8,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')	
	

1402	FORMAT('GaussPoints "XSeepage4T" Elemtype Tetrahedra'/,
	1	   'Number of Gauss Points: 4'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1412	FORMAT('GaussPoints "XSeepage4T" Elemtype Tetrahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","4,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')


1501	FORMAT ('GaussPoints "XHeat4" ElemType Quadrilateral'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

1511	FORMAT ('GaussPoints "XHeat4" ElemType Quadrilateral'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

1502	FORMAT ('GaussPoints "XHeat8" ElemType Quadrilateral'/
     2        'Number of Gauss Points:',I5,/
	3        'Natural Coordinates: Internal'/
	4        'Nodes Included'/
	5        'End GaussPoints'/)

1512	FORMAT ('GaussPoints "XHeat8" ElemType Quadrilateral'/   !EXCEL  
     2        "Number of Gauss Points:,","I5,",/
	3        "Natural Coordinates:,","Internal,",/
	4        'Nodes Included'/
	5        'End GaussPoints'/)


1503	FORMAT ('GaussPoints "XHeat3" ElemType Triangle'/
	1	   'Number of Gauss Points: 3'/,
	2        'Natural Coordinates: Given'/
	3        '   1.0  0.0'/
	4        '   0.0  0.0'/
	5        '   0.0  1.0'/
	7         'End GaussPoints'/)


1513	FORMAT ('GaussPoints "XHeat3" ElemType Triangle'/   !EXCEL  
	1	   "Number of Gauss Points,","3,",/,
	2     "Natural Coordinates:,","Given,",/
	3     "1.0,","0.0,"/
	4     "0.0,","0.0,",/
	5     "0.0,","1.0,",/
	7      'End GaussPoints'/)	
	

1601	FORMAT('GaussPoints "XHeat8H" Elemtype Hexahedra'/,
	1	   'Number of Gauss Points: 8'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1611	FORMAT('GaussPoints "XHeat8H" Elemtype Hexahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","8,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')	

1602	FORMAT('GaussPoints "XHeat4T" Elemtype Tetrahedra'/,
	1	   'Number of Gauss Points: 4'/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

1612	FORMAT('GaussPoints "XHeat4T" Elemtype Tetrahedra'/,   !EXCEL  
	1	   "Number of Gauss Points:,","4,",/,
	2	   "Natural Coordinates:,","internal,",/,
	3	   'end gausspoints')	


      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNEND(METHOD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------

	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	CASE('MIN')
	IP = 2
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) 
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,100) 

100	FORMAT('End Values'/)


      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================		
	SUBROUTINE PRNHEAD2(ISTEP,METHOD,COMPO)  !LINK ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD,COMPO
	CHARACTER(3) HED
	CHARACTER(30) HED1

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT	

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

	SELECTCASE(COMPO)
	CASE('FOC')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,101) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,201) HED,NAME(1:NAML),JSTEP   !EXCEL   	
	CASE('MOM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,102) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,202) HED,NAME(1:NAML),JSTEP   !EXCEL  
	ENDSELECT


101	FORMAT(/'Result "Link-Local-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "F-R" "F-S" "F-T"'/,
	3		'Values')

201	FORMAT(/'Result "Link-Local-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "F-R", "F-S", "F-T"'/,
	3		'Values')
	
102	FORMAT(/'Result "Link-Local-Moment ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "M-R" "M-S" "M-T"'/,
	3		'Values')

202	FORMAT(/'Result "Link-Local-Moment ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "M-R", "M-S", "M-T"'/,
	3		'Values')

      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================	
	SUBROUTINE PRNVALV2(EDATA,IGM,NOUT,METHOD,COMPO)  !Link ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10)
C     ----------------------------------------------------------------	
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	CASE('MIN')
	IP = 2
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN


	ND1 = INT(EDATA(1))
	ND2 = INT(EDATA(2))

	SELECTCASE(COMPO)
	CASE('FOC')
	NC  = 3
	NC1 = 1
	NC2 = 3
	CASE('MOM')
	NC  = 3
	NC1 = 4
	NC2 = 6
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	NP = 0
	N1 = 2+ NC1+NP
	N2 = 2+ NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	CASE('MIN')
	NP = NOUT
	N1 = 2+ NC1+NP
	N2 = 2+ NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	ENDSELECT


	IF(ABS(PDATA(1)).LT.1.0E-5) PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).LT.1.0E-5) PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).LT.1.0E-5) PDATA(3) = 0.0D0
	IF(ABS(PDATA(1)).EQ."NaN") PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).EQ."NaN") PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).EQ."NaN") PDATA(3) = 0.0D0

      IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) ND1,(PDATA(J),J=1,3) 
	IF(KEXCEL.EQ.1) THEN   !EXCEL  
	WRITE (LFPRIN,110) ND1,','
	DO J = 1,3
	IF(J.LT.3) WRITE(LFPRIN,120) PDATA(J),','
	IF(J.EQ.3) WRITE(LFPRIN,120) PDATA(J)
	ENDDO
	ENDIF
      
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) ND2,(PDATA(J),J=1,3) 
	IF(KEXCEL.EQ.1) THEN   !EXCEL  
	WRITE (LFPRIN,110) ND2,','
	DO J = 1,3
	IF(J.LT.3) WRITE(LFPRIN,120) PDATA(J),','
	IF(J.EQ.3) WRITE(LFPRIN,120) PDATA(J)
	ENDDO
	ENDIF


100	FORMAT(2X,I5,10E15.6)
110	FORMAT(I5,A,\)   !EXCEL  
120	FORMAT(E15.6,A,\)   !EXCEL  

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNHEAD3(ISTEP,METHOD,COMPO)  !GLOBAL SUPPORT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD,COMPO
	CHARACTER(3) HED
	CHARACTER(30) HED1

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON /OSPEC/ ISTIME,IOUTSPEC,NFATIGUE
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !

C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"GLobal Reaction"     '
	CASE(1)
	HED1 = '"Comb GLobal Reaction"'
      ENDSELECT	
      
      
	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)
      IF (IOUTSPEC.EQ.1) JSTEP = ISTIME
      IF (KFATL.NE.2) JSTEP = NFATIGUE-1 
	SELECTCASE(COMPO)
	CASE('FOC')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,101) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,201) HED,NAME(1:NAML),JSTEP  !EXCEL   
	CASE('MOM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,102) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,202) HED,NAME(1:NAML),JSTEP  !EXCEL   
	CASE('TEM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,103) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,203) HED,NAME(1:NAML),JSTEP  !EXCEL   
	CASE('PRS')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,104) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,204) HED,NAME(1:NAML),JSTEP  !EXCEL   
	ENDSELECT

C     ----------------------------------------------------------------

101	FORMAT(/'Result "Global Reaction-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "F-X" "F-Y" "F-Z"'/,
	3		'Values')
201	FORMAT(/'Result "Global Reaction-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "F-X", "F-Y", "F-Z"'/,
	3		'Values')


102	FORMAT(/'Result "Global Reaction-Moment ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "M-X" "M-Y" "M-Z"'/,
	3		'Values')
202	FORMAT(/'Result "Global Reaction-Moment ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "M-X", "M-Y", "M-Z"'/,
	3		'Values')


103	FORMAT(/'Result "Global Reaction-Temperature ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Temperature"'/,
	3		'Values')
203	FORMAT(/'Result "Global Reaction-Temperature ',A3,'"',2X,     !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames,  "Temperature"'/,
	3		'Values')


104	FORMAT(/'Result "Global Reaction-Flow-Potential ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Flow-Potential"'/,
	3		'Values')
204	FORMAT(/'Result "Global Reaction-Flow-Potential ',A3,'"',2X,     !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames,  "Flow-Potential"'/,
	3		'Values')

      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================	
	SUBROUTINE PRNVALV3(EDATA,ISN,NOUT,METHOD,COMPO)  !GLOBAL SUPPORT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON /INOU/  ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
      COMMON /PRINTREACTION/ NPRINTREACTION

C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	CASE('MIN')
	IP = 2
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN


	SELECTCASE(COMPO)
	CASE('FOC')
	NC  = 3
	NC1 = 1
	NC2 = 3
	CASE('MOM')
	NC  = 3
	NC1 = 4
	NC2 = 6
	CASE('TEM')
	NC  = 1
	NC1 = 8
	NC2 = 8
	CASE('PRS')
	NC  = 1
	NC1 = 9
	NC2 = 9
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	NP = 0
	N1 = NC1+NP
	N2 = NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	CASE('MIN')
	NP = NOUT
	N1 = NC1+NP
	N2 = NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	ENDSELECT


	IF(ABS(PDATA(1)).LT.1.0E-5) PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).LT.1.0E-5) PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).LT.1.0E-5) PDATA(3) = 0.0D0
      IF(ABS(PDATA(1)).EQ."NaN") PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).EQ."NaN") PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).EQ."NaN") PDATA(3) = 0.0D0

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) ISN,(PDATA(J),J=1,NC)
      
      SELECTCASE(COMPO)
	CASE('FOC')
	SELECTCASE(METHOD) ! CHANA
      CASE('MAX')
      
C     SONGSAK MODIFIED ... OCT2019
	CALL LOCATN('@SUP',KRID,NSAT,NSUP,1)  !CALLING NUMBER OF SUPPORT NSUP
      
      DO ISUP = 1,NSUP
          CALL INTFILL('@SUP',NODSP,1,ISUP,0)  !CALL NODE WHICH HAS BOUNDARY CONDITION 
          IF (NODSP.EQ.ISN)THEN
              !IF(NPRINTREACTION.EQ.1 .AND. KEXCEL.EQ.0) WRITE(ITO,*)' FORCE REACTION ( NODE, F-X , F-Y , F-Z )'
              !IF(NPRINTREACTION.EQ.1 .AND. KEXCEL.EQ.0) WRITE(ITO,100) ISN,(PDATA(J),J=1,NC)
              !IF(NPRINTREACTION.EQ.1 .AND. KEXCEL.EQ.0) WRITE(ITO,*)
              ! MODIFY BY TOEY
              WRITE (5203,200) ISN
              WRITE (5203,201) (PDATA(J),J=1,NC)
              EXIT
          ENDIF
      ENDDO
      
      ENDSELECT
      
      CASE('MOM')
      
	SELECTCASE(METHOD) ! CHANA
      CASE('MAX')
          
C     SONGSAK MODIFIED ... OCT2019
	CALL LOCATN('@SUP',KRID,NSAT,NSUP,1)  !CALLING NUMBER OF SUPPORT NSUP
      
      DO ISUP = 1,NSUP
          CALL INTFILL('@SUP',NODSP,1,ISUP,0)  !CALL NODE WHICH HAS BOUNDARY CONDITION 
          IF (NODSP.EQ.ISN)THEN
              !IF(NPRINTREACTION.EQ.1 .AND. KEXCEL.EQ.0) WRITE(ITO,*)' MOMENT REACTION ( NODE, M-X , M-Y , M-Z )' 
              !IF(NPRINTREACTION.EQ.1 .AND. KEXCEL.EQ.0) WRITE(ITO,100) ISN,(PDATA(J),J=1,NC)
              !IF(NPRINTREACTION.EQ.1 .AND. KEXCEL.EQ.0) WRITE(ITO,*)
              ! MODIFY BY TOEY
              WRITE (5204,200) ISN
              WRITE (5204,201) (PDATA(J),J=1,NC)
              EXIT
          ENDIF
      ENDDO
      
      ENDSELECT
      
      ENDSELECT
      
	IF(KEXCEL.EQ.1) THEN   !EXCEL  
	WRITE (LFPRIN,110) ISN,','
	DO J = 1,NC
	IF(J.LT.NC) WRITE (LFPRIN,120) PDATA(J),','
	IF(J.EQ.NC) WRITE (LFPRIN,120) PDATA(J)
      ENDDO
      ENDIF


100	FORMAT(2X,I5,10E15.6)
110	FORMAT(I5,A,\)   !EXCEL  
120   FORMAT(E15.6,A,\)   !EXCEL  
200   FORMAT (I5)
201   FORMAT (F20.7,2X,F20.7,2X,F20.7)
      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNHEAD4(ISTEP,METHOD,COMPO)  !LOCAL SUPPORT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD,COMPO
	CHARACTER(3) HED
	CHARACTER(30) HED1

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Local Reaction"     '
	CASE(1)
	HED1 = '"Comb Local Reaction"'
	ENDSELECT	

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

	SELECTCASE(COMPO)
	CASE('FOC')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,101) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,201) HED,NAME(1:NAML),JSTEP  !EXCEL  
	CASE('MOM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,102) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,202) HED,NAME(1:NAML),JSTEP  !EXCEL  
	CASE('TEM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,103) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,203) HED,NAME(1:NAML),JSTEP  !EXCEL  
	CASE('PRS')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,104) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,204) HED,NAME(1:NAML),JSTEP  !EXCEL  
	ENDSELECT

C     ----------------------------------------------------------------

101	FORMAT(/'Result "Local Reaction-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "F-X" "F-Y" "F-Z"'/,
	3		'Values')
201	FORMAT(/'Result "Local Reaction-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "F-X", "F-Y", "F-Z"'/,
	3		'Values')

102	FORMAT(/'Result "Local Reaction-Moment ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "M-X" "M-Y" "M-Z"'/,
	3		'Values')
202	FORMAT(/'Result "Local Reaction-Moment ',A3,'"',2X,   !EXCEL 
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "M-X", "M-Y", "M-Z"'/,
	3		'Values')

103	FORMAT(/'Result "Local Reaction-Temperature ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Temperature"'/,
	3		'Values')
203	FORMAT(/'Result "Local Reaction-Temperature ',A3,'"',2X,   !EXCEL 
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames , "Temperature"'/,
	3		'Values')


104	FORMAT(/'Result "Local Reaction-Flow-Potential ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Flow-Potential"'/,
	3		'Values')
204   FORMAT(/'Result "Local Reaction-Flow-Potential ',A3,'"',2X,   !EXCEL 
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames,  "Flow-Potential"'/,
	3		'Values')

      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================	
	SUBROUTINE PRNVALV4(EDATA,ISN,NOUT,METHOD,COMPO)  !Local SUPPORT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	CASE('MIN')
	IP = 2
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN


	SELECTCASE(COMPO)
	CASE('FOC')
	NC  = 3
	NC1 = 1
	NC2 = 3
	CASE('MOM')
	NC  = 3
	NC1 = 4
	NC2 = 6
	CASE('TEM')
	NC  = 1
	NC1 = 8
	NC2 = 8
	CASE('PRS')
	NC  = 1
	NC1 = 9
	NC2 = 9
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	NP = 0
	N1 = NC1+NP
	N2 = NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	CASE('MIN')
	NP = NOUT
	N1 = NC1+NP
	N2 = NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	ENDSELECT


	IF(ABS(PDATA(1)).LT.1.0E-5) PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).LT.1.0E-5) PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).LT.1.0E-5) PDATA(3) = 0.0D0
      IF(ABS(PDATA(1)).EQ."NaN") PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).EQ."NaN") PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).EQ."NaN") PDATA(3) = 0.0D0

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) ISN,(PDATA(J),J=1,NC)
	IF(KEXCEL.EQ.1) THEN   !EXCEL  
	WRITE (LFPRIN,110) ISN,','
      DO J = 1,NC
	IF(J.LT.NC) WRITE (LFPRIN,120) PDATA(J),','
	IF(J.EQ.NC) WRITE (LFPRIN,120) PDATA(J)
      ENDDO
      ENDIF

100	FORMAT(2X,I5,10E15.6)
110	FORMAT(I5,A,\)   !EXCEL  
120	FORMAT(E15.6,A,\)   !EXCEL  

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNHEAD5(ISTEP,METHOD,COMPO)  !FREE NODE DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD,COMPO
	CHARACTER(3) HED
	CHARACTER(30) HED1

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON /OSPEC/ ISTIME,IOUTSPEC,NFATIGUE
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !
C     ----------------------------------------------------------------	
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Displacement"     '
	CASE(1)
	HED1 = '"Comb Displacement"'
	ENDSELECT	

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)
	
      IF (IOUTSPEC.EQ.1) JSTEP = ISTIME
      IF (KFATL.NE.2) JSTEP = NFATIGUE-1
      
	SELECTCASE(COMPO)
	CASE('DIS')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,101) HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,111) HED,NAME(1:NAML),JSTEP  !EXCEL  
      
	CASE('ROT')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,102) HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,112) HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	CASE('TEM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,103) HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,113) HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	CASE('PRS')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,104) HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,114) HED,NAME(1:NAML),JSTEP  !EXCEL  
      
	CASE('WRP')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,105) HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,115) HED,NAME(1:NAML),JSTEP  !EXCEL  
      
	ENDSELECT

C     ----------------------------------------------------------------

101	FORMAT(/'Result "Displacement ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "U" "V" "W"'/,
	3		'Values')


111	FORMAT(/'Result "Displacement ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames, "U", "V", "W"'/,
	3		'Values')


102	FORMAT(/'Result "Disp & Rotation ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames  "U" "V" "W" "Rx" "Ry" "Rz"'/,
	3		'Values')

112	FORMAT(/'Result "Disp & Rotation ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames, "U", "V", "W", "Rx", "Ry", "Rz"'/,
	3		'Values')



103	FORMAT(/'Result "Temperature ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Temperature"'/,
	3		'Values')

113	FORMAT(/'Result "Temperature ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames, "Temperature"'/,
	3		'Values')


104	FORMAT(/'Result "Flow-Potential ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Flow-Potential"'/,
	3		'Values')

114	FORMAT(/'Result "Flow-Potential ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames, "Flow-Potential"'/,
	3		'Values')


105	FORMAT(/'Result "Warping ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Warping"'/,
	3		'Values')


115	FORMAT(/'Result "Warping ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames, "Warping"'/,
	3		'Values')	


      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================	
	SUBROUTINE PRNVALV5(EDATA,ISN,NOUT,METHOD,COMPO)  !FREE NODE DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	CASE('MIN')
	IP = 2
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN


	SELECTCASE(COMPO)
	CASE('DIS')
	NC  = 3
	NC1 = 1
	NC2 = 3
	CASE('ROT')
	NC  = 6
	NC1 = 1
	NC2 = 6
	CASE('TEM')
	NC  = 1
	NC1 = 8
	NC2 = 8
	CASE('PRS')
	NC  = 1
	NC1 = 9
	NC2 = 9
	CASE('WRP')
	NC  = 1
	NC1 = 7
	NC2 = 7
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	NP = 0
	N1 = NC1+NP
	N2 = NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	CASE('MIN')
	NP = NOUT
	N1 = NC1+NP
	N2 = NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	ENDSELECT


	IF(ABS(PDATA(1)).LT.1.0E-5) PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).LT.1.0E-5) PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).LT.1.0E-5) PDATA(3) = 0.0D0
	IF(ABS(PDATA(4)).LT.1.0E-5) PDATA(4) = 0.0D0
	IF(ABS(PDATA(5)).LT.1.0E-5) PDATA(5) = 0.0D0
	IF(ABS(PDATA(6)).LT.1.0E-5) PDATA(6) = 0.0D0
      IF(ABS(PDATA(1)).EQ."NaN") PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).EQ."NaN") PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).EQ."NaN") PDATA(3) = 0.0D0
      IF(ABS(PDATA(4)).EQ."NaN") PDATA(1) = 0.0D0
	IF(ABS(PDATA(5)).EQ."NaN") PDATA(2) = 0.0D0
	IF(ABS(PDATA(6)).EQ."NaN") PDATA(3) = 0.0D0

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) ISN,(PDATA(J),J=1,NC)
      
      IF(KEXCEL.EQ.1) THEN   !EXCEL  
         WRITE (LFPRIN,110) ISN,','
         DO J = 1,NC
         IF(J.LT.NC) WRITE (LFPRIN,120) PDATA(J),','
         IF(J.EQ.NC) WRITE (LFPRIN,120) PDATA(J)
         ENDDO
      ENDIF


100	FORMAT(2X,I5,10E15.6)
110	FORMAT(I5,A,\)   !EXCEL  
120	FORMAT(E15.3,A,\)   !EXCEL  



      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNHEAD6(ISTEP,METHOD,COMPO)  !GLOBAL SPRING DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD,COMPO
	CHARACTER(3) HED
	CHARACTER(30) HED1

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT	

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

	SELECTCASE(COMPO)
	CASE('FOC')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,101) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,201) HED,NAME(1:NAML),JSTEP  !EXCEL   	
	CASE('MOM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,102) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,202) HED,NAME(1:NAML),JSTEP  !EXCEL   
	ENDSELECT

C     ----------------------------------------------------------------

101	FORMAT(/'Result "Spring Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "F-X" "F-Y" "F-Z"'/,
	3		'Values')
201	FORMAT(/'Result "Spring Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "F-X", "F-Y", "F-Z"'/,
	3		'Values')

102	FORMAT(/'Result "Spring Moment ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "M-X" "M-Y" "M-Z"'/,
	3		'Values')
202	FORMAT(/'Result "Spring Moment ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "M-X", "M-Y", "M-Z"'/,
	3		'Values')



      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================	
	SUBROUTINE PRNVALV6(EDATA,ISN,NOUT,METHOD,COMPO)  !GLOBAL SPRING DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	CASE('MIN')
	IP = 2
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN


	SELECTCASE(COMPO)
	CASE('FOC')
	NC  = 3
	NC1 = 1
	NC2 = 3
	CASE('MOM')
	NC  = 3
	NC1 = 4
	NC2 = 6
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	NP = 0
	N1 = 1+ NC1+NP
	N2 = 1+ NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	CASE('MIN')
	NP = NOUT
	N1 = 1+ NC1+NP
	N2 = 1+ NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	ENDSELECT

	ND = INT(EDATA(1))
	IF(ABS(PDATA(1)).LT.1.0E-5) PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).LT.1.0E-5) PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).LT.1.0E-5) PDATA(3) = 0.0D0
      IF(ABS(PDATA(1)).EQ."NaN") PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).EQ."NaN") PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).EQ."NaN") PDATA(3) = 0.0D0

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) ND,(PDATA(J),J=1,NC)
	IF(KEXCEL.EQ.1) THEN   !EXCEL  
	WRITE (LFPRIN,110) ND,','
      DO J = 1,NC
	IF(J.LT.NC) WRITE (LFPRIN,120) PDATA(J),','
	IF(J.EQ.NC) WRITE (LFPRIN,120) PDATA(J)
      ENDDO
      ENDIF


100	FORMAT(2X,I5,10E15.6)
110	FORMAT(I5,A,\)   !EXCEL  
120	FORMAT(E15.6,A,\)   !EXCEL  


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNHEAD7(ISTEP,METHOD,COMPO)  !GLOBAL SPRING DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD,COMPO
	CHARACTER(3) HED
	CHARACTER(30) HED1

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT	

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

	SELECTCASE(COMPO)
	CASE('FOC')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,101) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,201) HED,NAME(1:NAML),JSTEP  !EXCEL   	
	CASE('MOM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,102) HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,202) HED,NAME(1:NAML),JSTEP  !EXCEL   
	ENDSELECT

C     ----------------------------------------------------------------

101	FORMAT(/'Result "Spring Force Local ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "F-X" "F-Y" "F-Z"'/,
	3		'Values')
201	FORMAT(/'Result "Spring Force Local ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "F-X", "F-Y", "F-Z"'/,
	3		'Values')

102	FORMAT(/'Result "Spring Moment Local ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "M-X" "M-Y" "M-Z"'/,
	3		'Values')
202	FORMAT(/'Result "Spring Moment Local ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "M-X", "M-Y", "M-Z"'/,
	3		'Values')


      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================	
	SUBROUTINE PRNVALV7(EDATA,ISN,NOUT,METHOD,COMPO)  !LOCAL SPRING DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
C     ----------------------------------------------------------------
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	CASE('MIN')
	IP = 2
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN


	SELECTCASE(COMPO)
	CASE('FOC')
	NC  = 3
	NC1 = 1
	NC2 = 3
	CASE('MOM')
	NC  = 3
	NC1 = 4
	NC2 = 6
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	NP = 0
	N1 = 1+ NC1+NP
	N2 = 1+ NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	CASE('MIN')
	NP = NOUT
	N1 = 1+ NC1+NP
	N2 = 1+ NC2+NP
	PDATA(1:NC) = EDATA(N1:N2)
	ENDSELECT

	ND = INT(EDATA(1))
	IF(ABS(PDATA(1)).LT.1.0E-5) PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).LT.1.0E-5) PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).LT.1.0E-5) PDATA(3) = 0.0D0
      IF(ABS(PDATA(1)).EQ."NaN") PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).EQ."NaN") PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).EQ."NaN") PDATA(3) = 0.0D0

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,100) ND,(PDATA(J),J=1,NC)
	IF(KEXCEL.EQ.1) THEN   !EXCEL  
	WRITE (LFPRIN,110) ND,','
      DO J = 1,NC
	IF(J.LT.NC) WRITE (LFPRIN,120) PDATA(J),','
	IF(J.EQ.NC) WRITE (LFPRIN,120) PDATA(J)
      ENDDO
      ENDIF


100	FORMAT(2X,I5,10E15.6)
110	FORMAT(I5,A,\)   !EXCEL  
120	FORMAT(E15.6,A,\)   !EXCEL  

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	
	SUBROUTINE PRNTYP2(ISTYP,ISTEP,METHOD)  !TRUSS&CABLE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------

      IEGO = 0
      
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 

	SELECTCASE(ISTYP)
	CASE(3)
	GOTO 100
	CASE(4)
	GOTO 200
	CASE(7)
	GOTO 300
      CASE(8)
      GOTO 400
	ENDSELECT

	RETURN

100	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HED,NAME(1:NAML),JSTEP	
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3100) HED,NAME(1:NAML),JSTEP	  !EXCEL  


	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 2
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 
      
	NV      = 1
	MAPP(1) = 1
	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1001	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)

	GOTO 9000


200   IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HED,NAME(1:NAML),JSTEP	
      IF(KEXCEL.EQ.0) WRITE (991,3001) HED,NAME(1:NAML),JSTEP         !CHANA
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3101) HED,NAME(1:NAML),JSTEP	  !EXCEL  



	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 2
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1002

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	NV      = 1
	MAPP(1) = 1
	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1002	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)


	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3002) HED,NAME(1:NAML),JSTEP	
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3102) HED,NAME(1:NAML),JSTEP	  !EXCEL  

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 2
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1003

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	NV      = 1
	MAPP(1) = 2
	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1003	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)

	GOTO 9000
C	-------------------------

300	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3003) HED,NAME(1:NAML),JSTEP	
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3103) HED,NAME(1:NAML),JSTEP	  !EXCEL  


	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 2
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1004

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	NV      = 1
	MAPP(1) = 1
	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1004	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)

	GOTO 9000
      
400	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3004) HED,NAME(1:NAML),JSTEP	
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3104) HED,NAME(1:NAML),JSTEP	  !EXCEL  


	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 2
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1005

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	NV      = 1
	MAPP(1) = 1
	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1005	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)

	GOTO 9000

C	-------------------------


9000	CONTINUE

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)

3000	FORMAT(/'Result "XTruss-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XTruss"'/,
	2		'ComponentNames  "Truss-Axial-Force"'/,
	3		'Values')

3100	FORMAT(/'Result "XTruss-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XTruss"'/,
	2		'ComponentNames, "Truss-Axial-Force"'/,
	3		'Values')	

3001	FORMAT(/'Result "XCable-Net-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCable"'/,
	2		'ComponentNames  "ParaCable-Net-Force"'/,
	3		'Values')

3101	FORMAT(/'Result "XCable-Net-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCable"'/,
	2		'ComponentNames, "ParaCable-Net-Force"'/,
	3		'Values')	

3002	FORMAT(/'Result "XCable-Inc-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCable"'/,
	2		'ComponentNames  "ParaCable-Inc-Force"'/,
	3		'Values')


3102	FORMAT(/'Result "XCable-Inc-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCable"'/,
	2		'ComponentNames, "ParaCable-Inc-Force"'/,
	3		'Values')	

3003	FORMAT(/'Result "XCable-Catenary-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCableC"'/,
	2		'ComponentNames  "Catenary-Net-Hor-Force"'/,
	3		'Values')

3103  FORMAT(/'Result "XCable-Catenary-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCableC"'/,
	2		'ComponentNames, "Catenary-Net-Hor-Force"'/,
	3		'Values')	
      
3004	FORMAT(/'Result "XMooring-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCableC"'/,
	2		'ComponentNames  "Mooring-Net-Hor-Force"'/,
	3		'Values')

3104	FORMAT(/'Result "XMooring-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XCableC"'/,
	2		'ComponentNames, "Mooring-Net-Hor-Force"'/,
	3		'Values')	

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP5(ISTYP,ISTEP,METHOD)  !FRAME ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME,NAMEG
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
      COMMON /OPTION/ NDESIGN
C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
	COMMON /OSPEC/ ISTIME,IOUTSPEC,NFATIGUE
C      COMMON /FATIGUE_STRESS/ CHANAS1(10000,7)
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /OUTFA/ CHANA,ICHANA,IFATIGUE

C     OUTPUT INDEX 
      COMMON /OUTPUTCONTROL/ NLFORCE,NLSTRESS,NSTRESSP,NTORSION,NREATION,NDISPROTAWQ
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10),ASTRESS(12,7)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO         = 0
        
	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 

      IF (IOUTSPEC.EQ.1) JSTEP = ISTIME
      
      IF (KFATL.NE.2.AND.METHOD.EQ.'MAX') THEN
      JSTEP    = IFATIGUE + 1
      IFATIGUE = IFATIGUE + 1
      ENDIF
      
      IF (NLFORCE.NE.1) GOTO 10
      
      IF (KFATL.NE.2.AND.METHOD.EQ.'MIN') JSTEP   =  IFATIGUE
      
      
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HED,NAME(1:NAML),JSTEP
C      WRITE (230,REC=1) HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3100) HED,NAME(1:NAML),JSTEP   !EXCEL  

	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]
	MTMOD = 0
      IF (METHOD.EQ."MAX") CALL DEFINDLOCALSTRESS_GID (IGM,"CLEA",HED,NAME,NAML,JSTEP)  
	DO IGM = 1,NGIDM
          
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
          
	ITEST = 0
	ITYPE = 5
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1
          
	IF(ITEST.EQ.1) GOTO 1001
	MTMOD = MTMD

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 
      
	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)
	

      IF (METHOD.EQ."MAX") CALL DEFINDLOCALSTRESS_GID (IGM,"BODY",HED,NAME,NAML,JSTEP)  
      
      
1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)
      
10    IF (NLSTRESS.NE.1) GOTO 20
      IF (METHOD.EQ."MAX") THEN
      CALL DEFINDLOCALSTRESS_GID (IGM,"HEAD",HED,NAME,NAML,JSTEP)  
	DO IGM = 1,NGIDM
          
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
          
	ITEST = 0
	ITYPE = 5
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1
          
	IF(ITEST.EQ.1) GOTO 1100
      
      CALL DEFINDLOCALSTRESS_GID (IGM,"WRIT",HED,NAME,NAML,JSTEP) 
      
1100	CONTINUE

      ENDDO
      
      CALL PRNEND(METHOD)
      
      ENDIF
	
20    IF (NTORSION.NE.1) GOTO 30
      

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3002) HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3102) HED,NAME(1:NAML),JSTEP  !EXCEL  
	NV = 1
	MAPP(1:1) = [7]
	MTMOD = 0

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	
	ITEST = 0
	ITYPE = 5
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1002
	MTMOD = MTMD

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)
	! ---------------------------------

1002	CONTINUE

      ENDDO
	
      
	CALL PRNEND(METHOD)
      
30    IF (NSTRESSP.NE.1) GOTO 40


	IF(MTMOD.NE.1) GOTO 9000

C
	CALL INTFILL('#GMP',NEGRP,1 ,1 ,0)	!CALL TOTAL GID GROUP NUMBER
	DO IEGP = 1,NEGRP

	CALL PRNGPHD(NAMEG,NAMLG,IEGP)

	
	NV = 6
	MAPP(1:6) = [8,9,10,11,12,13]


	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	CALL INTFILL('@GMP',IEGRP,1 ,IGM,0)	!GET GID GROUP NUMBER FOR A GIVEN GID ELEMENT

	ITEST = 0
	ITYPE = 5
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1
	IF(IEGP.NE.IEGRP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1005

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF


      IF(KEXCEL.EQ.0) 
	1WRITE (LFPRIN,3001) NAMEG(1:NAMLG),HED,NAME(1:NAML),JSTEP
      
C     ---------- PRINT FATIGUE HEAD --------- 
	IF(KEXCEL.EQ.0.AND.METHOD.EQ.'MAX'.AND.CHANA.EQ.0.AND.KFATL.NE.2)THEN
	CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL:') !READ TIME DATA FOR OFFSHORE LOAD
	WRITE(989,*) TIME , 'FREQUENCY'
C 	write(991,*)' ELEMENT STRESS_1  STRESS_2  STRESS_3  STRESS_4  STRESS_5  STRESS_6'
C      WRITE (991,*)  NAMEG(1:NAMLG),HED,NAME(1:NAML),JSTEP !JSTEP , NEGRP   !CHANA
      CHANA = 1
      ENDIF
C     ---------- PRINT FATIGUE HEAD --------- 
      
      
      IF(KEXCEL.EQ.1) 
	1WRITE (LFPRIN,3101) NAMEG(1:NAMLG),HED,NAME(1:NAML),JSTEP !EXCEL
      

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

      
1005	CONTINUE

	ENDDO

	IF(ICOUNT.NE.0) CALL PRNEND(METHOD)
	ENDDO
	

C
	CALL INTFILL('#GMP',NEGRP,1 ,1 ,0)	!CALL TOTAL GID GROUP NUMBER
	DO IEGP = 1,NEGRP

	CALL PRNGPHD(NAMEG,NAMLG,IEGP)


	NV = 6
	MAPP(1:6) = [14,15,16,17,18,19]


	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	CALL INTFILL('@GMP',IEGRP,1 ,IGM,0)	!GET GID GROUP NUMBER FOR A GIVEN GID ELEMENT

	ITEST = 0
	ITYPE = 5
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1
	IF(IEGP.NE.IEGRP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1006

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
	IF(KEXCEL.EQ.0) 
	1WRITE (LFPRIN,3003) NAMEG(1:NAMLG),HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) 
	1WRITE (LFPRIN,3103) NAMEG(1:NAMLG),HED,NAME(1:NAML),JSTEP  !EXCEL
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1006	CONTINUE

	ENDDO

	IF(ICOUNT.NE.0) CALL PRNEND(METHOD)
	ENDDO
	
C
9000	CONTINUE

40      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)

3000	FORMAT(/'Result "XFrame-Local-Force ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		'ComponentNames',
     3		' "Axial-Force" "Shear-S" "Shear-T"',
	4		' "Torsion" "Moment-S" "Moment-T"'/,
	5		'Values')


3100	FORMAT(/'Result "XFrame-Local-Force ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		"ComponentNames,",
     3		' "Axial-Force", "Shear-S", "Shear-T", 
     4"Torsion", "Moment-S", "Moment-T"'/,
	5		'Values')


3001	FORMAT(/'Result ','"Section-Stress-1-6 ',A,' ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		'ComponentNames',
     3		' "S1" "S2" "S3" "S4" "S5" "S6"'/,
	4		'Values')

3101	FORMAT(/'Result ','"Section-Stress-1-6 ',A,' ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		"ComponentNames,",
     3		' "S1", "S2", "S3", "S4", "S5", "S6"'/,
	4		'Values')	


3003	FORMAT(/'Result ','"Section-Stress-7-12 ',A,' ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		'ComponentNames',
     3		' "S7" "S8" "S9" "S10" "S11" "S12"'/,
	4		'Values')

	
3103	FORMAT(/'Result ','"Section-Stress-7-12 ',A,' ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		"ComponentNames,",
     3		' "S7", "S8", "S9", "S10", "S11", "S12"'/,
	4		'Values')

3002	FORMAT(/'Result "XFrame-Warping-Torsion ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		'ComponentNames',
     3		' "Warping-Torsion"'/,
	4		'Values')

3102	FORMAT(/'Result "XFrame-Warping-Torsion ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		"ComponentNames,",
     3		' "Warping-Torsion"'/,
	4		'Values')

3200	FORMAT(/'Result "Unity Check',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		'ComponentNames',
     3		' "Axial" "Shear-S" "Shear-T"',
	4		' "Bending-S" "Bending-T" "Max Combined"'/,
	5		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP6(ISTYP,ISTEP,METHOD)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------

      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(0,1,2)
	HSTP = 'XPlane8'
	CASE(3,4,5)
	HSTP = 'XPlane4'
	CASE DEFAULT
	RETURN
	ENDSELECT


	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
	

	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	 ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,3100) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL  
       ENDIF
      
	NV = 6
	MAPP(1:6) = [1,2,3,4,0,0]

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 6
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)
	
1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,3101) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL  
	ENDIF
	
	NV = 1
	MAPP(1) = 5


	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 6
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1002

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1002	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)


      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A7,'-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A7/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')

3100	FORMAT(/'Result "',A7,'-Stress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A7/,
	2		"ComponentNames,",
     3		'"S-X", "S-Y", "S-Z", "S-XY", "S-XZ", "S-YZ"'/,
	5		'Values')


3001	FORMAT(/'Result "',A7,'-VStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A7/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')

3101	FORMAT(/'Result "',A7,'-VStress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A7/,
	2		"ComponentNames,",
     3		'"VMise"'/,
	4		'Values')	


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP6N(ISTYP,ISTEP,METHOD)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE VALNOD(:,:),LM(:)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('%NUB',NSN,1,1 ,0)  !NUMBER OF NODES
	
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------

      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(0,1,2)
	HSTP = 'XPlane8'
	CASE(3,4,5)
	HSTP = 'XPlane4'
	CASE DEFAULT
	RETURN
	ENDSELECT


	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
	

C     -------------------------------
C     CALCULATE AVERAGE FACTOR
C     -------------------------------
      LEN = 1 + 5
	ALLOCATE(VALNOD(LEN,NSN))
	VALNOD(1:LEN,1:NSN) = 0.0D0
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	CALL INTFILL('OGRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !

      
	ITEST = 0
	ITYPE = 6
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 50

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

	ALLOCATE(LM(NNM))
	
      READ(NFL1,REC=IEL) AF1(1:N1) 
      
      CALL READCON(IGM,LM,NFL4) 
      DO INM = 1,NNM
      NOD = LM(INM)
      VALNOD(1,NOD) = VALNOD(1,NOD) + 1.0D0

          DO ILEN = 1,LENGTH
	        SELECTCASE(METHOD)
	        CASE('MAX')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1))
	        CASE('MIN')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1)+NOUT)
	        ENDSELECT
	        IF(ABS(VAL).LT.1.0E-12) VAL = 0.0D0    
	        VALNOD(1+ILEN,NOD) = VALNOD(1+ILEN,NOD) + VAL
          ENDDO
      
      ENDDO

	DEALLOCATE(LM)
	      
50	CONTINUE

	ENDDO
C     -------------------------------
C     -------------------------------


C     -------------------------------------------------
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3100) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	NV = 6
	MAPP(1:6) = [1,2,3,4,0,0]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1001
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1001	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C     -------------------------------------------------
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3101) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	NV = 1
	MAPP(1) = 5
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1002
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1002	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------


	DEALLOCATE(VALNOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
      
3000	FORMAT(/'Result "',A7,'-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')

3100	FORMAT(/'Result "',A7,'-NodalStress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		"ComponentNames,",
     3		'"S-X", "S-Y", "S-Z", "S-XY", "S-XZ", "S-YZ"'/,
	5		'Values')


3001	FORMAT(/'Result "',A7,'-NodalVStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')

3101	FORMAT(/'Result "',A7,'-NodalVStress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		"ComponentNames,",
     3		'"VMise"'/,
	4		'Values')	


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP8(ISTYP,ISTEP,METHOD)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------

      IEGO = 0
      
	HSTP = 'XPlaneC'

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


      IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
       ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,3100) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL  
	ENDIF
	
	NV = 6
	MAPP(1:6) = [1,2,3,4,0,0]


	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 8
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 
	
	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)


	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	 ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,3101) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL  
      ENDIF
      
	NV = 1
	MAPP(1) = 5


	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 8
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1002

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1002	CONTINUE

	ENDDO
	

	CALL PRNEND(METHOD)


      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)

3000	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')

3100	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
	2		"ComponentNames,",
     3		'"S-X", "S-Y", "S-Z", "S-XY", "S-XZ", "S-YZ"'/,
	5		'Values')

3001	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
	
3101	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP8N(ISTYP,ISTEP,METHOD)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE VALNOD(:,:),LM(:)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('%NUB',NSN,1,1 ,0)  !NUMBER OF NODES

	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------

      IEGO = 0
      
	HSTP = 'XPlaneC'

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 



C     -------------------------------
C     CALCULATE AVERAGE FACTOR
C     -------------------------------
      LEN = 1 + 5
	ALLOCATE(VALNOD(LEN,NSN))
	VALNOD(1:LEN,1:NSN) = 0.0D0
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	CALL INTFILL('OGRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !

      
	ITEST = 0
	ITYPE = 8
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 50

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
	ALLOCATE(LM(NNM))
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL READCON(IGM,LM,NFL4) 
      DO INM = 1,NNM
      NOD = LM(INM)
      VALNOD(1,NOD) = VALNOD(1,NOD) + 1.0D0

          DO ILEN = 1,LENGTH
	        SELECTCASE(METHOD)
	        CASE('MAX')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1))
	        CASE('MIN')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1)+NOUT)
	        ENDSELECT
	        IF(ABS(VAL).LT.1.0E-12) VAL = 0.0D0    
	        VALNOD(1+ILEN,NOD) = VALNOD(1+ILEN,NOD) + VAL
          ENDDO
      
      ENDDO

	DEALLOCATE(LM)
	      
50	CONTINUE

	ENDDO
C     -------------------------------
C     -------------------------------

C     -------------------------------------------------
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3100) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	NV = 6
	MAPP(1:6) = [1,2,3,4,0,0]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1001
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1001	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C     -------------------------------------------------
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3101) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
      
	NV = 1
	MAPP(1) = 5
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1002
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1002	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------


	DEALLOCATE(VALNOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')

3100	FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames, "S-X", "S-Y", "S-Z", "S-XY", "S-XZ", "S-YZ"'/,
	5		'Values')
	

3001	FORMAT(/'Result "',A9,'-NodalVStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')

3101	FORMAT(/'Result "',A9,'-NodalVStress ',A3,'"',2X,   !EXCEL  
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP9(ISTYP,ISTEP,METHOD)  !SHELL ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(8)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
     
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(11)
	HSTP = 'XShell3Q'
	CASE(1)
	HSTP = ' XShell4'
	CASE(4)
	HSTP = 'XShell4Q'
	CASE(12)
	HSTP = ' XShell9'
	CASE(3)
	HSTP = ' XShell8'
      CASE(13) ! FOR SHEAR WALL BY BJ
	HSTP = '  XWALL4'	
	CASE DEFAULT
	RETURN
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3010) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL 
	
	NV = 3
	MAPP(1:3) = [1,2,0]
	
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1001	CONTINUE

      ENDDO
      
	CALL PRNEND(METHOD)

      
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3011) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
	 
	NV = 3
	MAPP(1:3) = [3,4,5]
	
	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1002

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1002	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)


	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3002) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3012) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
		
	NV = 3
	MAPP(1:3) = [6,7,8]

	MTMOD = 0
	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1003

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
	MTMOD = MTMD

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1003	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)


	IF(MTMOD.EQ.4) THEN

	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3003) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3013) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
	
	NV = 6
	MAPP(1:6) = [9,10,0,11,0,0]

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1004

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1004	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)

C
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3004) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3014) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
	
	NV = 6
	MAPP(1:6) = [12,13,0,14,0,0]

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1005

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1005	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)

	ENDIF




	IF(MTMOD.EQ.1) THEN

	DO ILAYR = 1,3

	NV = 6

	IF(ILAYR.EQ.1) THEN
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,4001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,4011) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
	
	MAPP(1:6) = [9,10,11,12,13,14]
	ELSEIF(ILAYR.EQ.2) THEN
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,4002) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,4012) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
	
	MAPP(1:6) = [15,16,17,18,19,20]
	ELSEIF(ILAYR.EQ.3) THEN
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,4003) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,4013) HSTP,HED,NAME(1:NAML),JSTEP,HSTP	  !EXCEL
		
	MAPP(1:6) = [21,22,23,24,25,26]
	ENDIF
	
	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1006

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1006	CONTINUE

	ENDDO

	CALL PRNEND(METHOD)

	ENDDO

	ENDIF


      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A8,'-Membrane-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"'/,
	4		'Values')
3010	FORMAT(/'Result "',A8,'-Membrane-Stress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T"'/,
	4		   'Values')	

3001	FORMAT(/'Result "',A8,'-Bending-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-RR" "S-SS" "S-RS"'/,
	4		'Values')
3011	FORMAT(/'Result "',A8,'-Bending-Stress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-RR", "S-SS", "S-RS"'/,
	4		'Values')
	
3002	FORMAT(/'Result "',A8,'-Shear-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST" "S-RS"'/,
	4		'Values')
3012	FORMAT(/'Result "',A8,'-Shear-Stress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-RT", "S-ST", "S-RS"'/,
	4		'Values')

3003	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
3013	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S" ,"S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

3004	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
3014	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')


4001	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4011	FORMAT(/'Result "',A8,'-Fiber-Stress-Top ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

4002	FORMAT(/'Result "',A8,'-Fiber-Stress-Mid ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4012	FORMAT(/'Result "',A8,'-Fiber-Stress-Mid ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

4003	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4013	FORMAT(/'Result "',A8,'-Fiber-Stress-Bot ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A8/,
	2		'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')


      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP9N(ISTYP,ISTEP,METHOD)  !SHELL ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(8)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE VALNOD(:,:),LM(:)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('%NUB',NSN,1,1 ,0)  !NUMBER OF NODES

	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(11)
	HSTP = 'XShell3Q'
	CASE(1)
	HSTP = ' XShell4'
	CASE(4)
	HSTP = 'XShell4Q'
	CASE(12)
	HSTP = ' XShell9'
	CASE(3)
	HSTP = ' XShell8'
      	CASE(13) ! FOR SHEAR WALL BY BJ
	HSTP = '  XWALL4'
	CASE DEFAULT
	RETURN
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 



      MTMOD = 0
C     -------------------------------
C     CALCULATE AVERAGE FACTOR
C     -------------------------------
      LEN = 1 + (8+6*3)
	ALLOCATE(VALNOD(LEN,NSN))
	VALNOD(1:LEN,1:NSN) = 0.0D0
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	CALL INTFILL('OGRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !

      
	ITEST = 0
	ITYPE = 9
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 50

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
	ALLOCATE(LM(NNM))
	
	MTMOD = MTMD
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL READCON(IGM,LM,NFL4) 
      DO INM = 1,NNM
      NOD = LM(INM)
      VALNOD(1,NOD) = VALNOD(1,NOD) + 1.0D0

          DO ILEN = 1,LENGTH
	        SELECTCASE(METHOD)
	        CASE('MAX')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1))
	        CASE('MIN')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1)+NOUT)
	        ENDSELECT
	        IF(ABS(VAL).LT.1.0E-12) VAL = 0.0D0    
	        VALNOD(1+ILEN,NOD) = VALNOD(1+ILEN,NOD) + VAL
          ENDDO
      
      ENDDO

	DEALLOCATE(LM)
	      
50	CONTINUE

	ENDDO
C     -------------------------------
C     -------------------------------

C     -------------------------------------------------
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3010) HSTP,HED,NAME(1:NAML),JSTEP !EXCEL  
	
	NV = 3
	MAPP(1:3) = [1,2,0]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1001
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1001	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C     -------------------------------------------------
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3011) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	NV = 3
	MAPP(1:3) = [3,4,5]

	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1002
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1002	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C     -------------------------------------------------
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3002) HSTP,HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3012) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	NV = 3
	MAPP(1:3) = [6,7,8]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1003
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1003	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C     =================================================
	IF(MTMOD.EQ.1) THEN

	DO ILAYR = 1,3

C     -------------------------------------------------
	NV = 6

	IF(ILAYR.EQ.1) THEN
	    
	    IF(KEXCEL.EQ.0) WRITE (LFPRIN,4001) HSTP,HED,NAME(1:NAML),JSTEP
	    IF(KEXCEL.EQ.1) WRITE (LFPRIN,4011) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	    MAPP(1:6) = [9,10,11,12,13,14]
	ELSEIF(ILAYR.EQ.2) THEN
	   
	    IF(KEXCEL.EQ.0)  WRITE(LFPRIN,4002) HSTP,HED,NAME(1:NAML),JSTEP
	    IF(KEXCEL.EQ.1) WRITE (LFPRIN,4012) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	    
	    MAPP(1:6) = [15,16,17,18,19,20]
	ELSEIF(ILAYR.EQ.3) THEN
	   
	    IF(KEXCEL.EQ.0)  WRITE(LFPRIN,4003) HSTP,HED,NAME(1:NAML),JSTEP
	    IF(KEXCEL.EQ.1) WRITE (LFPRIN,4013) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	    
	    MAPP(1:6) = [21,22,23,24,25,26]
	ENDIF
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1004
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1004	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

	ENDDO

	ENDIF
C     =================================================

C     =================================================
	IF(MTMOD.EQ.4) THEN

      IF(KEXCEL.EQ.0) WRITE (LFPRIN,3003) HSTP,HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3013) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	NV = 6
	MAPP(1:6) = [9,10,0,11,0,0]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1005
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1005	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3004) HSTP,HED,NAME(1:NAML),JSTEP	
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3014) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL  
	
	NV = 6
	MAPP(1:6) = [12,13,0,14,0,0]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1006
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1006	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

	ENDIF
C     =================================================


	DEALLOCATE(VALNOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A8,'-Membrane-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"'/,
	4		'Values')
3010	FORMAT(/'Result "',A8,'-Membrane-NodalStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames, "S-R", "S-S", "S-T"'/,
	4		   'Values')
	
3001	FORMAT(/'Result "',A8,'-Bending-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-RR" "S-SS" "S-RS"'/,
	4		'Values')
3011	FORMAT(/'Result "',A8,'-Bending-NodalStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames, "S-RR", "S-SS", "S-RS"'/,
	4		'Values')	
	
3002	FORMAT(/'Result "',A8,'-Shear-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST" "S-RS"'/,
	4		'Values')
3012	FORMAT(/'Result "',A8,'-Shear-NodalStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames, "S-RT" ,"S-ST" ,"S-RS"'/,
	4		'Values')

3003	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Top ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
3013	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Top ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2	'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

3004	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Bot ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
3014	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Bot ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
     2   'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

4001	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Top ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4011	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Top ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2	'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

4002	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Mid ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4012	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Mid ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2	'ComponentNames, "S-R", "S-S", "S-T","S-RS", "S-RT", "S-ST"'/,
	5		'Values')

4003	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Bot ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"',
	4		' "S-RS" "S-RT" "S-ST"'/,
	5		'Values')
4013	FORMAT(/'Result "',A8,'-Fiber-NodalStress-Bot ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2  'ComponentNames, "S-R", "S-S", "S-T" ,"S-RS", "S-RT", "S-ST"'/,
	5		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP10(ISTYP,ISTEP,METHOD)  !SOLID ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(4)
CC	HSTP = 'XSolidANS'
      HSTP = 'XSolid-SH'
	NV = 9
	MAPP(1:9) = [1,2,3,4,5,6,7,8,9]
	CASE(6)
	HSTP = 'XSolidEAS'
	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]
	CASE(13)
	HSTP = 'XSolidTeT'
	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]
	CASE DEFAULT
	RETURN
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT


      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 



      IF(KEXCEL.EQ.0) THEN
       WRITE (LFPRIN,4010) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,3010) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	ENDIF
	

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

      
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)


	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	ENDIF
	IF(KEXCEL.EQ.1) THEN
	 WRITE (LFPRIN,3011) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	ENDIF
	
	
	IF(ISTYP.EQ.4) THEN ! FOR SOLID-SHELL
	  NV = 1
	  MAPP(1) = 10
	ELSE
	  NV = 1
	  MAPP(1) = 7
	ENDIF

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1002

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1002	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')
3010	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-X", "S-Y", "S-Z", "S-XY", "S-XZ", "S-YZ"'/,
	5    'Values')

3001	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
3011	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')

C    ====================== FOR SOLID-SHELL PRINT =========================

4000	FORMAT(/'Result "',A10,'-Membrane Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"'/,
	4		'Values')

4001	FORMAT(/'Result "',A10,'-Membrane Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-R", "S-S", "S-T"'/,
	5    'Values')


4010	FORMAT(/'Result "',A10,'-Bending Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-RR" "S-SS" "S-TT"'/,
	4		'Values')	

4011	FORMAT(/'Result "',A10,'-Bending Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-RR", "S-SS", "S-TT"'/,
	5    'Values')

4020	FORMAT(/'Result "',A10,'-Transverse Normal Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		'"S-TT"'/,
	4		'Values')
	
4021	FORMAT(/'Result "',A10,'-Transverse Normal Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-TT"'/,
	5    'Values')	
	
4030	FORMAT(/'Result "',A10,'-Transverse Shear Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST"'/,
	4		'Values')	

4031	FORMAT(/'Result "',A10,'-Transverse Shear Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-RT", "S-ST"'/,
	5    'Values')


4040	FORMAT(/'Result "',A10,'-VStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
4041	FORMAT(/'Result "',A10,'-VStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP10N(ISTYP,ISTEP,METHOD)  !SOLID ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE VALNOD(:,:),LM(:)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('%NUB',NSN,1,1 ,0)  !NUMBER OF NODES

	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(4)
	HSTP = 'XSolidANS'
	CASE(6)
	HSTP = 'XSolidEAS'
	CASE(13)
	HSTP = 'XSolidTeT'
	CASE DEFAULT
	RETURN
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT


      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


C     -------------------------------
C     CALCULATE AVERAGE FACTOR
C     -------------------------------
      LEN = 1 + 7
	ALLOCATE(VALNOD(LEN,NSN))
	VALNOD(1:LEN,1:NSN) = 0.0D0
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	CALL INTFILL('OGRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !

      
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 50

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

	ALLOCATE(LM(NNM))

      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL READCON(IGM,LM,NFL4) 
      DO INM = 1,NNM
      NOD = LM(INM)
      VALNOD(1,NOD) = VALNOD(1,NOD) + 1.0D0

          DO ILEN = 1,LENGTH
	        SELECTCASE(METHOD)
	        CASE('MAX')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1))
	        CASE('MIN')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1)+NOUT)
	        ENDSELECT
	        IF(ABS(VAL).LT.1.0E-12) VAL = 0.0D0    
	        VALNOD(1+ILEN,NOD) = VALNOD(1+ILEN,NOD) + VAL
          ENDDO
      
      ENDDO

	DEALLOCATE(LM)
	      
50	CONTINUE

	ENDDO
C     -------------------------------
C     -------------------------------

C     -------------------------------------------------
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3010) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL 
      
	
	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1001
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1001	CONTINUE
	ENDDO


      SELECTCASE(METHOD)
      CASE('MAX')
      IF(KEXCEL.EQ.0) THEN
         WRITE (120,3012) HSTP,HED,NAME(1:NAML),JSTEP
         CALL INTFILL('#GMP',NEGEL,1,2,0)  
         CALL INTFILL('#GMP',NEGRP,1,1,0)  
	   ENDIF
	ENDSELECT
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C     -------------------------------------------------
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3011) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL
	
	NV = 1
	MAPP(1) = 7
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1002
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1002	CONTINUE
	ENDDO
	
      SELECTCASE(METHOD)
      CASE('MAX')
      IF(KEXCEL.EQ.0) THEN
         WRITE (120,3014) HSTP,HED,NAME(1:NAML),JSTEP
         CALL INTFILL('#GMP',NEGEL,1,2,0)  
         CALL INTFILL('#GMP',NEGRP,1,1,0)  
	   ENDIF
	ENDSELECT
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------


	DEALLOCATE(VALNOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')
3010	FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
     2    'ComponentNames, "S-X", "S-Y", "S-Z","S-XY", "S-XZ", "S-YZ"'/,
	5		'Values')

3001	FORMAT(/'Result "',A9,'-NodalVStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
3011	FORMAT(/'Result "',A9,'-NodalVStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')
3012  FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes')	
3014  FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================      
      SUBROUTINE PRNTYP104(ISTYP,ISTEP,METHOD)  !SOLID-SHELL ELEMENT DATA (ON GAUSS POINTS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
      HSTP = 'XSolidSHL'


	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT


      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


CB      IF(KEXCEL.EQ.0) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON GAUSS POINTS","-",HED," ",">"
CB      ENDIF
CB      IF(KEXCEL.EQ.1) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON GAUSS POINTS","-",HED," ",">"   !EXCEL
CB	ENDIF
	

C ------------------ MEMBRANE STRESS PART --------------------------
      IF(KEXCEL.EQ.0) THEN
       WRITE (LFPRIN,4000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,4001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	ENDIF
	

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
     
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

	NV = 3
	MAPP(1:3) = [1,2,7]
		
      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD) ! "End Values" in flavia.res file
	
C ------------------ BENDING STRESS PART --------------------------
      IF(KEXCEL.EQ.0) THEN
       WRITE (LFPRIN,4010) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,4011) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	ENDIF
	
	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
     
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1002

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

	NV = 3
	!!MAPP(1:3) = [4,5,6]
      MAPP(1:3) = [3,4,5]
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1002	CONTINUE

	ENDDO

      CALL PRNEND(METHOD)
      
C ------------------ SHEAR STRESS PART --------------------------
      IF(KEXCEL.EQ.0) THEN
       WRITE (LFPRIN,4020) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,4021) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	ENDIF
	

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
     
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1003

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

	NV = 3
      !!MAPP(1:3) = [8,9,3]
      MAPP(1:3) = [6,8,9]

      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1003	CONTINUE

	ENDDO

      CALL PRNEND(METHOD)

C ------------------ TRANSVERSE STRESS PART --------------------------
C ------ THIS PART WAS INCLUDED IN UPPER PART, SO CANCELD ------------

CB      IF(KEXCEL.EQ.0) THEN
CB       WRITE (LFPRIN,4030) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
CB      ENDIF
CB      IF(KEXCEL.EQ.1) THEN
CB       WRITE (LFPRIN,4031) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
CB	ENDIF
	

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
     
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1004

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
	NV = 2
	MAPP(1:2) = [8,9]
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

CB      CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
CB     1              LFPRIN,1,NPT)

1004	CONTINUE

	ENDDO

CB      CALL PRNEND(METHOD)

C ---------------------- VONMISE STRESS PART ------------------------	
	
	!!IF(KEXCEL.EQ.0) THEN
	!!WRITE(LFPRIN,4040) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	!!ENDIF
	!!IF(KEXCEL.EQ.1) THEN
	!! WRITE (LFPRIN,4041) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	!!ENDIF

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1005

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

	NV = 1
	MAPP(1) = 10
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

      !!CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
     !!1              LFPRIN,1,NPT)

1005	CONTINUE

	ENDDO
	
	!!CALL PRNEND(METHOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


C    ====================== FOR SOLID-SHELL PRINT =========================

4000	FORMAT(/'Result "',A9,'-Membrane Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"'/,
	4		'Values')

4001	FORMAT(/'Result "',A9,'-Membrane Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-R", "S-S", "S-T"'/,
	5    'Values')


4010	FORMAT(/'Result "',A9,'-Bending Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-RR" "S-SS" "S-RS"'/,
	4		'Values')	

4011	FORMAT(/'Result "',A9,'-Bending Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-RR", "S-SS", "S-RS"'/,
	5    'Values')

4020	FORMAT(/'Result "',A9,'-Shear Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST" "S-RS"'/,
	4		'Values')
	
4021	FORMAT(/'Result "',A9,'-Shear Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-RT", "S-ST", "S-RS"'/,
	5    'Values')	
	
4030	FORMAT(/'Result "',A9,'-Transverse Shear Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST"'/,
	4		'Values')	

4031	FORMAT(/'Result "',A9,'-Transverse Shear Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-RT", "S-ST"'/,
	5    'Values')


4040	FORMAT(/'Result "',A9,'-VRStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
4041	FORMAT(/'Result "',A9,'-VRStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')


      RETURN
      END


C	=====================================================================
C	=====================================================================
C	=====================================================================           
	SUBROUTINE PRNTYP104N(ISTYP,ISTEP,METHOD)  !SOLID-SHELL ELEMENT DATA (ON NODAL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE VALNOD(:,:),LM(:)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('%NUB',NSN,1,1 ,0)  !NUMBER OF NODES

	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
 
      HSTP = 'XSolidSHL'


	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT


      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


C     -------------------------------
C     CALCULATE AVERAGE FACTOR
C     -------------------------------
      
      LEN = 1 + 10
      
	ALLOCATE(VALNOD(LEN,NSN))
	
	VALNOD(1:LEN,1:NSN) = 0.0D0
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	CALL INTFILL('OGRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !
    
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 50

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

	ALLOCATE(LM(NNM))
	

      READ(NFL1,REC=IEL) AF1(1:N1) 
		 	  
	CALL READCON(IGM,LM,NFL4) 
	 
	 DO INM = 1,NNM
	  NOD = LM(INM)
        VALNOD(1,NOD) = VALNOD(1,NOD) + 1.0D0
        
        DO ILEN = 1,LENGTH
           SELECTCASE(METHOD)
	        CASE('MAX')
CCB	        BUM = NOD*LENGTH+ILEN
CCB	        VAL = AF1(NOD*LENGTH+ILEN)
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1))                            
	        CASE('MIN')
CCB	        VAL = AF1(NOD*LENGTH+ILEN+NOUT)
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1)+NOUT)             
	        ENDSELECT
	        
	     IF(ABS(VAL).LT.1.0E-12) VAL = 0.0D0    
	        VALNOD(1+ILEN,NOD) = VALNOD(1+ILEN,NOD) + VAL
        
        ENDDO
	 
	 ENDDO
 
	DEALLOCATE(LM)
	      
50	CONTINUE

	ENDDO
C     -------------------------------
C     -------------------------------

CB      IF(KEXCEL.EQ.0) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON NODES","-",HED," ",">"
CB      ENDIF
CB      IF(KEXCEL.EQ.1) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON NODES","-",HED," ",">"  !EXCEL
CB	ENDIF
	

C -------------------------- MEMBRANE STRESS PART -----------------------------
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,4000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,4001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	
	
	NV = 3
	MAPP(1:3) = [1,2,7]

	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1001
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1001	CONTINUE
	ENDDO

      CALL PRNEND(METHOD)	

C --------------------------- BENDING STRESS PART -----------------------------
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,4010) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,4011) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	
	
	NV = 3
	!!MAPP(1:3) = [4,5,6]
      MAPP(1:3) = [3,4,5]

	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1002
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1002	CONTINUE
	ENDDO

      CALL PRNEND(METHOD)
      
C --------------------------- TRANSVERSE NORMAL STRESS PART -------------------
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,4020) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,4021) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	
	
	NV = 3
      !!MAPP(1:3) = [8,9,3]
      MAPP(1:3) = [6,8,9]
      
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1003
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1003	CONTINUE
	ENDDO

      CALL PRNEND(METHOD)
      
C ------------------------ TRANSVERSE SHEAR STRESS PART -----------------------
C -------------- THIS PART WAS INCLUDED IN UPPER PART, SO CANCELD -------------
CB      IF(KEXCEL.EQ.0) WRITE (LFPRIN,4030) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
CB      IF(KEXCEL.EQ.1) WRITE (LFPRIN,4031) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	
	
	!NV = 2
	!MAPP(1:2) = [8,9]

	!DO ISN = 1,NSN
      !FAC = VALNOD(1,ISN)
	!IF(FAC.EQ.0.0D0) GOTO 1004
	!FAC = 1.0D0/FAC
CB    CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
!1004	CONTINUE
!	ENDDO
	
CB	CALL PRNEND(METHOD)

C ----------------------- VONMISES PART --------------------------------------
	
	!!IF(KEXCEL.EQ.0) WRITE (LFPRIN,4040) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      !!IF(KEXCEL.EQ.1) WRITE (LFPRIN,4041) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL
		
      NV = 1
      MAPP(1) = 10
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1005
	FAC = 1.0D0/FAC
      !!CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1005	CONTINUE
	ENDDO
	
	!!CALL PRNEND(METHOD)
C     -------------------------------------------------


	DEALLOCATE(VALNOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


C    ====================== FOR SOLID-SHELL PRINT =========================

4000	FORMAT(/'Result "',A9,'-Membrane-Nodal Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes',2X,A10/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T"'/,
	4		'Values')

4001	FORMAT(/'Result "',A9,'-Membrane-Nodal Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnNodes',2X,A9/,
     2  'ComponentNames, "S-R", "S-S", "S-T"'/,
	5    'Values')


4010	FORMAT(/'Result "',A9,'-Bending-Nodal Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes',2X,A10/,
	2		'ComponentNames',
     3		' "S-RR" "S-SS" "S-RS"'/,
	4		'Values')	

4011	FORMAT(/'Result "',A9,'-Bending-Nodal Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnNodes',2X,A9/,
     2  'ComponentNames, "S-RR", "S-SS", "S-RS"'/,
	5    'Values')

4020	FORMAT(/'Result "',A9,'-Shear-Nodal Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes',2X,A10/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST" "S-RS"'/,
	4		'Values')
	
4021	FORMAT(/'Result "',A9,'-Shear-Nodal Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnNodes',2X,A9/,
     2  'ComponentNames, "S-RT", "S-ST", "S-RS"'/,
	5    'Values')	
	
4030	FORMAT(/'Result "',A9,'-Transverse Shear-Nodal Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnNodes',2X,A10/,
	2		'ComponentNames',
     3		' "S-RT" "S-ST"'/,
	4		'Values')	

4031	FORMAT(/'Result "',A9,'-Transverse Shear-Nodal Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Vector',2X,'OnNodes',2X,A9/,
     2  'ComponentNames, "S-RT", "S-ST"'/,
	5    'Values')


4040	FORMAT(/'Result "',A9,'-Nodal VRStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes',2X,A9/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
4041	FORMAT(/'Result "',A9,'-Nodal VRStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================        
      SUBROUTINE PRNTYP104_ST(ISTYP,ISTEP,METHOD)  !SOLID-SHELL ELEMENT DATA (ON GAUSS POINTS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF2(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
      HSTP = 'XSolidSHL'


	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT


      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


CB      IF(KEXCEL.EQ.0) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON GAUSS POINTS","-",HED," ",">"
CB      ENDIF
CB      IF(KEXCEL.EQ.1) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON GAUSS POINTS","-",HED," ",">"   !EXCEL
CB	ENDIF
	

C ------------------ STRESS PART --------------------------
      IF(KEXCEL.EQ.0) THEN
       WRITE (LFPRIN,4000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,4001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	ENDIF
	

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
      CALL INTFILL('OGRF',NFL5 ,15,IEG,0) !FOR SOLID-SHELL(LIU) STRESS OUTPUT BY BJ
     
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF2)) DEALLOCATE(AF2)
	    ALLOCATE(AF2(N1))
          READ(NFL5,REC=IEL) AF2(1:N1)  
	    IEGO = IEG
	ENDIF

	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]
		
      CALL PRNVALVL(IGM,AF2(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD) ! "End Values" in flavia.res file
	

C ---------------------- VONMISE STRESS PART ------------------------	
	
	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,4040) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	ENDIF
	IF(KEXCEL.EQ.1) THEN
	 WRITE (LFPRIN,4041) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL
	ENDIF

	DO IGM = 1,NGIDM

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
      CALL INTFILL('OGRF',NFL5 ,15,IEG,0) !FOR SOLID-SHELL(LIU) STRESS OUTPUT BY BJ

	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1005

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF2)) DEALLOCATE(AF2)
	    ALLOCATE(AF2(N1))
          READ(NFL5,REC=IEL) AF2(1:N1)  
	    IEGO = IEG
	ENDIF

	NV = 1
	MAPP(1) = 10
	
      CALL PRNVALVL(IGM,AF2(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1005	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

      IF(ALLOCATED(AF2)) DEALLOCATE(AF2)


C    ====================== FOR SOLID-SHELL PRINT =========================

4000	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T" "S-RS" "S-RT" "S-ST"'/,
	4		'Values')

4001	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
     2  'ComponentNames, "S-R", "S-S", "S-T","S-RS","S-RT","S-ST"'/,
	5    'Values')


4040	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "VMise-Stress"'/,
	4		'Values')
4041	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')


      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================           
	SUBROUTINE PRNTYP104N_ST(ISTYP,ISTEP,METHOD)  !SOLID-SHELL ELEMENT DATA (ON NODAL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE VALNOD(:,:),LM(:)
	ALLOCATABLE AF2(:)
C     ----------------------------------------------------------------
	CALL INTFILL('%NUB',NSN,1,1 ,0)  !NUMBER OF NODES

	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
 
      HSTP = 'XSolidSHL'


	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT


      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


C     -------------------------------
C     CALCULATE AVERAGE FACTOR
C     -------------------------------
      
      LEN = 1 + 10
      
	ALLOCATE(VALNOD(LEN,NSN))
	
	VALNOD(1:LEN,1:NSN) = 0.0D0
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	CALL INTFILL('OGRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !
      CALL INTFILL('OGRF',NFL5 ,15,IEG,0) !FOR SOLID-SHELL(LIU) STRESS OUTPUT BY BJ
    
	ITEST = 0
	ITYPE = 10
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 50

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF2)) DEALLOCATE(AF2)
	    ALLOCATE(AF2(N1))
          READ(NFL5,REC=IEL) AF2(1:N1)  
	    IEGO = IEG
	ENDIF

	ALLOCATE(LM(NNM))
	
		 	  
	CALL READCON(IGM,LM,NFL4) 
	 
	 DO INM = 1,NNM
	  NOD = LM(INM)
        VALNOD(1,NOD) = VALNOD(1,NOD) + 1.0D0
        
        DO ILEN = 1,LENGTH
           SELECTCASE(METHOD)
	        CASE('MAX')
CCB	        BUM = NOD*LENGTH+ILEN
CCB	        VAL = AF1(NOD*LENGTH+ILEN)
	        VAL = AF2(ILEN+LENGTH*(NPT)+LENGTH*(INM-1))
	        CASE('MIN')
CCB	        VAL = AF1(NOD*LENGTH+ILEN+NOUT)
	        VAL = AF2(ILEN+LENGTH*(NPT)+LENGTH*(INM-1)+NOUT)
	        ENDSELECT
	        
	     IF(ABS(VAL).LT.1.0E-12) VAL = 0.0D0    
	        VALNOD(1+ILEN,NOD) = VALNOD(1+ILEN,NOD) + VAL
        
        ENDDO
	 
	 ENDDO
 
	DEALLOCATE(LM)
	      
50	CONTINUE

	ENDDO
C     -------------------------------
C     -------------------------------

CB      IF(KEXCEL.EQ.0) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON NODES","-",HED," ",">"
CB      ENDIF
CB      IF(KEXCEL.EQ.1) THEN
CB       WRITE (LFPRIN,*) "< SOLID-SHELL STRESS ON NODES","-",HED," ",">"  !EXCEL
CB	ENDIF
	

C -------------------------- NODAL STRESS PART -----------------------------
      IF(KEXCEL.EQ.0) WRITE (LFPRIN,4000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,4001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	
	
	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]

	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1001
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1001	CONTINUE
	ENDDO

      CALL PRNEND(METHOD)	


C ----------------------- VONMISES PART --------------------------------------
	
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,4040) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,4041) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL
		
      NV = 1
      MAPP(1) = 10
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1005
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1005	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------


	DEALLOCATE(VALNOD)

      IF(ALLOCATED(AF2)) DEALLOCATE(AF2)


C    ====================== FOR SOLID-SHELL PRINT =========================

4000	FORMAT(/'Result "',A9,'-Nodal Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes',2X,A10/,
	2		'ComponentNames',
     3		' "S-R" "S-S" "S-T" "S-RS" "S-RT" "S-ST" '/,
	4		'Values')

4001	FORMAT(/'Result "',A9,'-Membrane-Nodal Stress ',A3,'"',2X,   !EXCEL
	1	 '"',A,'"',
	1   2X,I5,2X,'Matrix',2X,'OnNodes',2X,A9/,
     2  'ComponentNames, "S-R", "S-S", "S-T", "S-RS", "S-RT", "S-ST"'/,
	5    'Values')

4040	FORMAT(/'Result "',A9,'-Nodal VStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes',2X,A9/,
	2		'ComponentNames',
     3		' "VMise-Nodal Stress"'/,
	4		'Values')
4041	FORMAT(/'Result "',A9,'-Nodal VStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')

      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP11(ISTYP,ISTEP,METHOD)  !SOLID ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	HSTP = 'XSolidCon'

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,3010) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	ENDIF
	
	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]


	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 11
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN,1,NPT)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)


	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	ENDIF
	IF(KEXCEL.EQ.1) THEN
	 WRITE (LFPRIN,3011) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	ENDIF
	
	NV = 1
	MAPP(1) = 7



	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 11
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1002

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVL(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,
     1              LFPRIN,1,NPT)

1002	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')
3010	FORMAT(/'Result "',A9,'-Stress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,A9/,
	2	'ComponentNames, "S-X", "S-Y", "S-Z","S-XY", "S-XZ", "S-YZ"'/,
	5		'Values')

3001	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
3011	FORMAT(/'Result "',A9,'-VStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP11N(ISTYP,ISTEP,METHOD)  !SOLID ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE VALNOD(:,:),LM(:)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('%NUB',NSN,1,1 ,0)  !NUMBER OF NODES

	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	HSTP = 'XSolidCon'

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Element Stress"     '
	CASE(1)
	HED1 = '"Comb Element Stress"'
	ENDSELECT


	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)

C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 


C     -------------------------------
C     CALCULATE AVERAGE FACTOR
C     -------------------------------
      LEN = 1 + 7
	ALLOCATE(VALNOD(LEN,NSN))
	VALNOD(1:LEN,1:NSN) = 0.0D0
	DO IGM = 1,NGIDM

	ITYP = 0
	ISTY = 0
C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',MTMD ,3 ,IEG,0)  !
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NPT  ,5 ,IEG,0) !
	CALL INTFILL('OGRP',NNM  ,7 ,IEG,0) !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !
	CALL INTFILL('OGRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('OGRF',NFL4 ,14,IEG,0) !

      
	ITEST = 0
	ITYPE = 11
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 50

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF
	
	ALLOCATE(LM(NNM))
	
      READ(NFL1,REC=IEL) AF1(1:N1) 

      CALL READCON(IGM,LM,NFL4) 
      DO INM = 1,NNM
      NOD = LM(INM)
      VALNOD(1,NOD) = VALNOD(1,NOD) + 1.0D0

          DO ILEN = 1,LENGTH
	        SELECTCASE(METHOD)
	        CASE('MAX')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1))
	        CASE('MIN')
	        VAL = AF1(ILEN+LENGTH*(NPT)+LENGTH*(INM-1)+NOUT)
	        ENDSELECT
	        IF(ABS(VAL).LT.1.0E-12) VAL = 0.0D0    
	        VALNOD(1+ILEN,NOD) = VALNOD(1+ILEN,NOD) + VAL
          ENDDO
      
      ENDDO

	DEALLOCATE(LM)
	      
50	CONTINUE

	ENDDO
C     -------------------------------
C     -------------------------------

C     -------------------------------------------------
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP
      IF(KEXCEL.EQ.1) WRITE (LFPRIN,3010) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL 
	
	NV = 6
	MAPP(1:6) = [1,2,3,4,5,6]
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1001
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1001	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------

C     -------------------------------------------------
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,3011) HSTP,HED,NAME(1:NAML),JSTEP  !EXCEL 
	
	NV = 1
	MAPP(1) = 7
	
	DO ISN = 1,NSN
      FAC = VALNOD(1,ISN)
	IF(FAC.EQ.0.0D0) GOTO 1002
	FAC = 1.0D0/FAC
      CALL PRNVALVN(VALNOD(2,ISN),ISN,FAC,MAPP,NV)
1002	CONTINUE
	ENDDO
	
	CALL PRNEND(METHOD)
C     -------------------------------------------------


	DEALLOCATE(VALNOD)
	
      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)



3000	FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "S-X" "S-Y" "S-Z"',
	4		' "S-XY" "S-XZ" "S-YZ"'/,
	5		'Values')
3010	FORMAT(/'Result "',A9,'-NodalStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2	'ComponentNames, "S-X", "S-Y", "S-Z","S-XY", "S-XZ", "S-YZ"'/,
	5		'Values')
3001	FORMAT(/'Result "',A9,'-NodalVStress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames',
     3		' "VMise"'/,
	4		'Values')
3011	FORMAT(/'Result "',A9,'-NodalVStress ',A3,'"',2X,   !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames, "VMise"'/,
	4		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP13(ISTYP,ISTEP,METHOD,NNM)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(3)
	IF(NNM.EQ.4) HSTP = 'XSeepage4'
	IF(NNM.EQ.8) HSTP = 'XSeepage8'
	CASE(4)
	HSTP = 'XSeepage3'
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Velocity Field"     '
	CASE(1)
	HED1 = '"Comb Velocity Field"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
	


	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
      ENDIF
      IF(KEXCEL.EQ.1) THEN
       WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL 
	ENDIF
	
	NV = 3
	MAPP(1:3) = [1,2,0]

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 13
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 
      
	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A9,'-Velocity ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames',
     3		' "V-X" "V-Y" "V-Z"'/,
	4		'Values')
3001	FORMAT(/'Result "',A9,'-Velocity ',A3,'"',2X,      !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A9/,
	2		'ComponentNames, "V-X", "V-Y", "V-Z"'/,
	4		'Values')



      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP14(ISTYP,ISTEP,METHOD)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(3)
	HSTP = 'XSeepage8H'
	CASE(4)
	HSTP = 'XSeepage4T'
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Velocity Field"     '
	CASE(1)
	HED1 = '"Comb Velocity Field"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
	


      IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	ENDIF
	IF(KEXCEL.EQ.1) THEN
	 WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL 
	ENDIF
	
	NV = 3
	MAPP(1:3) = [1,2,3]

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !


	ITEST = 0
	ITYPE = 15
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1


	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
      

3000	FORMAT(/'Result "',A10,'-Velocity ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames',
     3		' "V-X" "V-Y" "V-Z"'/,
	4		'Values')
3001	FORMAT(/'Result "',A10,'-Velocity ',A3,'"',2X,     !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Vector',2X,'OnGaussPoints',2X,A10/,
	2		'ComponentNames, "V-X", "V-Y", "V-Z"'/,
	4		'Values')


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP15(ISTYP,ISTEP,METHOD,NNM)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(3)
	IF(NNM.EQ.4) HSTP = 'XHeat4'
	IF(NNM.EQ.8) HSTP = 'XHeat8'
	CASE(4)
	HSTP = 'XHeat3'
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Temperature Field"     '
	CASE(1)
	HED1 = '"Comb Temperature Field"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
	


	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	ENDIF
	IF(KEXCEL.EQ.1) THEN
	 WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP   !EXCEL 
	ENDIF
	
	NV = 3
	MAPP(1:3) = [1,2,0]

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 15
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)



3000	FORMAT(/'Result "',A6,'-Temperature ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A7/,
	2		'ComponentNames',
     3		' "T" '/,
	4		'Values')
3001	FORMAT(/'Result "',A6,'-Temperature ',A3,'"',2X,     !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A7/,
	2		'ComponentNames, "T" '/,
	4		'Values')




      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP16(ISTYP,ISTEP,METHOD)  !PLANE ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
	CHARACTER(30) HED1
	CHARACTER(9)  HSTP

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION NPM(10),NPI(10),MAPP(10)
	ALLOCATABLE AF1(:)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
      IEGO = 0
      
	SELECTCASE(ISTYP)
	CASE(3)
	HSTP = 'XHeat8H'
	CASE(4)
	HSTP = 'XHeat4T'
	ENDSELECT

	SELECTCASE(METHOD)
	CASE('MAX')
	IP = 1
	HED = 'MAX'
	CASE('MIN')
	IP = 2
	HED = 'MIN'
	ENDSELECT

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

	SELECTCASE(INM)
	CASE(0)
	HED1 = '"Temperature Field"     '
	CASE(1)
	HED1 = '"Comb Temperature Field"'
	ENDSELECT

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


C	CALL GID ELEMENT NUMBER  NGIDM
	CALL LOCATN('OGDM',KGDM,NN,NGIDM,1) 
	


	IF(KEXCEL.EQ.0) THEN
	WRITE(LFPRIN,3000) HSTP,HED,NAME(1:NAML),JSTEP,HSTP
	ENDIF
	IF(KEXCEL.EQ.1) THEN
	 WRITE (LFPRIN,3001) HSTP,HED,NAME(1:NAML),JSTEP,HSTP  !EXCEL
	ENDIF
	
	NV = 3
	MAPP(1:3) = [1,2,3]

	DO IGM = 1,NGIDM


C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('OGDM',IEG  ,1 ,IGM,0)
	CALL INTFILL('OGDM',IEL  ,2 ,IGM,0)
	CALL INTFILL('OGRP',ITYP ,1 ,IEG,0)  !ITYPE
	CALL INTFILL('OGRP',ISTY ,2 ,IEG,0)  !ISTYP
	CALL INTFILL('OGRP',NELE ,4 ,IEG,0)  !
	CALL INTFILL('OGRP',NOUT ,12,IEG,0)  !
	CALL INTFILL('OGRP',LENGTH,17,IEG,0) !
C	GROUP FILE
	CALL INTFILL('OGRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('OGRF',NFL1 ,11,IEG,0) !

	ITEST = 0
	ITYPE = 16
	IF(ITYP.NE.ITYPE) ITEST = 1
	IF(ISTY.NE.ISTYP) ITEST = 1

	IF(ITEST.EQ.1) GOTO 1001

      IF(IEG.NE.IEGO) THEN
          IF(ALLOCATED(AF1)) DEALLOCATE(AF1)
	    ALLOCATE(AF1(N1))
	    IEGO = IEG
	ENDIF

      READ(NFL1,REC=IEL) AF1(1:N1) 

	CALL PRNVALVE(IGM,AF1(1),NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)

1001	CONTINUE

	ENDDO
	
	CALL PRNEND(METHOD)

      IF(ALLOCATED(AF1)) DEALLOCATE(AF1)


3000	FORMAT(/'Result "',A7,'-Temperature ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A7/,
	2		'ComponentNames',
     3		' "T" '/,
	4		'Values')
3001	FORMAT(/'Result "',A7,'-Temperature ',A3,'"',2X,    !EXCEL
	1	    '"',A,'"',
	1		2X,I5,2X,'Scalar',2X,'OnGaussPoints',2X,A7/,
	2		'ComponentNames, "T" '/,
	4		'Values')

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE PRNTYP17(ISTEP)  !TENDON DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
	
C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN

      COMMON A(9000000),IA(9000000)

	DIMENSION LEST(2*NEG)


C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',INM,1,12,0)  !STORe PRINTING FLAG (0=STANDARD 1=COMBINATION)
	CALL INTFILL('NPRN',IND,1,13,0)  !STORE PRINTING FLAG (0=PRINT ALL 1=JUST ONE STEP)
C     ----------------------------------------------------------------
	IF(INM.NE.0) RETURN

	DO IEG = 1,2*NEG
	LEST(IEG) = IA(LEL+IEG-1)
	ENDDO

C     ----------------------------------------------
	DO 1000 IEG = 1,NEG
	KEG = IEG
	NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     
      IF (A(LMP).LE.0.001D0)THEN
      REWIND (NELEMA)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)  
      ENDIF
      CALL MOVLEV (2)

	IF(ITYPE.EQ.17) THEN
	CALL TDOUTP(IA(LTDN),ISTEP,INM,IND) 
	ENDIF

1000	CONTINUE

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================


	SUBROUTINE PRNVALVE(IGM,EDATA,NOUT,LENGTH,MAPP,NV,METHOD,LFPRIN)  !ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD
	CHARACTER*3 BLANK
	
C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRN1,LFPRN2,LFPRN3,LFPRN4	
	
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(20),MAPP(1)
C     ----------------------------------------------------------------


	NUM = NOUT/LENGTH

	DO IPT = 1,NUM

	DO J = 1,NV
	PDATA(J) = 0.0D0
	M = MAPP(J)

	IF(M.NE.0) THEN
	SELECTCASE(METHOD)
	CASE('MAX')
	PDATA(J) = EDATA(M+LENGTH*(IPT-1))
	CASE('MIN')
	PDATA(J) = EDATA(M+LENGTH*(IPT-1)+NOUT)
	ENDSELECT
	IF(ABS(PDATA(J)).LT.1.0E-5) PDATA(J) = 0.0D0
	ENDIF

	ENDDO


      IF(IPT.EQ.1) THEN
	  IF(KEXCEL.EQ.0) WRITE(LFPRIN,100) IGM,PDATA(1:NV)
C                        WRITE(230,REC=1) IGM,(PDATA(J),J=1,NV)
        IF(KEXCEL.EQ.1) THEN
         WRITE(LFPRIN,110) IGM,','
         DO J = 1,NV
         IF(J.LT.NV)  WRITE(LFPRIN,120) PDATA(J),','
         IF(J.EQ.NV)  WRITE(LFPRIN,120) PDATA(J)   !EXCEL
	ENDDO
        ENDIF
        
      ENDIF
         
      IF(IPT.NE.1) THEN
        IF(KEXCEL.EQ.0) WRITE(LFPRIN,101) PDATA(1:NV)
C                        WRITE(230,REC=1) (PDATA(J),J=1,NV)
        IF(KEXCEL.EQ.1) THEN
         WRITE(LFPRIN,111) BLANK,','
         DO J = 1,NV
         IF(J.LT.NV)  WRITE(LFPRIN,120) PDATA(J),','
         IF(J.EQ.NV)  WRITE(LFPRIN,120) PDATA(J)   !EXCEL
         ENDDO
        ENDIF
        
      ENDIF
                
      ENDDO


100	FORMAT(2X,I5,20E15.6)
 110	FORMAT(I5,A,\)     !EXCEL
 120	FORMAT(E15.6,A,\)  !EXCEL
   
101	FORMAT(2X,5X,20E15.6)
 111  FORMAT(A,\)        !EXCEL


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================	
      SUBROUTINE PRNVALVL(IGM,EDATA,NOUT,LENGTH,MAPP,NV,METHOD,
     1                      LFPRIN,NPR1,NPR2)  !ELEMENT DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD
	
C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRN1,LFPRN2,LFPRN3,LFPRN4	
	
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(20),MAPP(1)
C     ----------------------------------------------------------------


	NUM = NOUT/LENGTH

	DO IPT = NPR1,NPR2

	DO J = 1,NV
	PDATA(J) = 0.0D0
	M = MAPP(J)

	IF(M.NE.0) THEN
	SELECTCASE(METHOD)
	CASE('MAX')
	PDATA(J) = EDATA(M+LENGTH*(IPT-1))
	CASE('MIN')
	PDATA(J) = EDATA(M+LENGTH*(IPT-1)+NOUT)
	ENDSELECT
	IF(ABS(PDATA(J)).LT.1.0E-5) PDATA(J) = 0.0D0
	ENDIF

	ENDDO


      IF(IPT.EQ.1) THEN
      
	  IF(KEXCEL.EQ.0) WRITE(LFPRIN,100) IGM,PDATA(1:NV)
        IF(KEXCEL.EQ.1) THEN
         WRITE(LFPRIN,110) IGM,','
         DO J = 1,NV
         IF(J.LT.NV)  WRITE(LFPRIN,120) PDATA(J),','
         IF(J.EQ.NV)  WRITE(LFPRIN,120) PDATA(J)   !EXCEL
	ENDDO
        ENDIF 
              
      ENDIF 
             
      IF(IPT.NE.1) THEN
      
        IF(KEXCEL.EQ.0) WRITE(LFPRIN,101) PDATA(1:NV)
        IF(KEXCEL.EQ.1) THEN
         WRITE(LFPRIN,111) BLANK,','
         DO J = 1,NV
         IF(J.LT.NV)  WRITE(LFPRIN,120) PDATA(J),','
         IF(J.EQ.NV)  WRITE(LFPRIN,120) PDATA(J)   !EXCEL
         ENDDO
        ENDIF
        
      ENDIF
        
      ENDDO

100	FORMAT(2X,I5,20E15.6)
 110	FORMAT(I5,A,\)     !EXCEL
 120	FORMAT(E15.6,A,\)  !EXCEL
   
101	FORMAT(2X,5X,20E15.6)
 111  FORMAT(A,\)        !EXCEL

      RETURN
      END
C

C	=====================================================================
C	=====================================================================	

	SUBROUTINE PRNVALVN(EDATA,ISN,FAC,MAPP,NV)  !FREE NODE DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4

C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10),MAPP(1)
C     ----------------------------------------------------------------
	CALL INTFILL('NPRN',LLC,1,11,0)  !PRINT TOT OR POS OR NEG OR ALL
C     ----------------------------------------------------------------

      CALL INTFILL('OUTC',IPOS,1,1,0)
      CALL INTFILL('OUTC',INEG,1,2,0)
	IF(IP.EQ.1.AND.IPOS.NE.1) RETURN
	IF(IP.EQ.2.AND.INEG.NE.1) RETURN

      DO IV = 1,NV
      PDATA(IV) = 0.0D0
      MP = MAPP(IV)
      IF(MP.NE.0) PDATA(IV) = EDATA(MP)*FAC
      ENDDO
      CALL MAXVALUE(PDATA,ISN,NV)
      
	IF(KEXCEL.EQ.0) WRITE(LFPRIN,100) ISN,PDATA(1:NV)
	IF(KEXCEL.EQ.1) THEN   !EXCEL
	WRITE (LFPRIN,110) ISN,','
	DO J = 1,NV
	IF(J.LT.NV) WRITE(LFPRIN,120) PDATA(J),','
	IF(J.EQ.NV) WRITE(LFPRIN,120) PDATA(J)
	ENDDO
	ENDIF


100	FORMAT(2X,I5,10E15.6)
110   FORMAT(I5,A,\)       !EXCEL
120   FORMAT(E15.6,A,\)


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE MODPRN  !PRINT EIGENVALUE MODE SHAPE DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 COMPO
      
C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
C     ----------------------------------------------------------------
	DIMENSION IDOF(9),IPRN(4),EDAT(9)
	ALLOCATABLE REIG(:)

	CALL INTFILL('%IOL',ITO,1,2,0)  !SCREEN
	CALL INTFILL('%IOL',ISO,1,3,0)  !OUTPUT.OUT
	
      CALL INTFILL('%IOL',IGO,1,4 ,0)  !GIDOUT.FLAVIA
      LFPRIN = IGO
      
	CALL INTFILL('%NUB',NSN  ,1,1,0)

	CALL INTFILL('EIGM',ISOLV,1,1,0)
	CALL INTFILL('EIGM',IEIG ,1,2,0) !2-BUCKLING EIGEN SOLUTION,	3-NATUTAL FREQUANCEY
	CALL INTFILL('EIGM',NROOT,1,3,0)
	CALL INTFILL('EIGM',IVPRT,1,4,0)
	CALL INTFILL('EIGM',NFILM,1,5,0)
	IF(NROOT.LE.0) RETURN


	ALLOCATE(REIG(NROOT))
	

C	---------------------------------------------
	SELECTCASE(IEIG)
	CASE(2)
	WRITE(ISO,2000) 
	WRITE(ITO,4000)  
	WRITE(10,4000) 
	CASE(3)
	WRITE(ISO,2001)               !(PERIOD)
	WRITE(ITO,4001) 
      WRITE(10,4001) 
	ENDSELECT

	DO IROOT = 1,NROOT
	CALL RELFILL('EIVV',REIG(IROOT),1,IROOT,0)
	IF(IEIG.EQ.3) THEN
	REIG(IROOT) = SQRT(REIG(IROOT))/2.0/3.141526 !(FREQUENZY)
	REIG(IROOT) = 1.0/REIG(IROOT)                !(PERIOD)
C	EIGVPeriod(IROOT)=REIG(IROOT)
	ENDIF
	ENDDO

	DO IROOT = 1,NROOT
	WRITE(ISO,3000) IROOT,REIG(IROOT)
	WRITE(ITO,5000) IROOT,REIG(IROOT)
      WRITE(10,5000) IROOT,REIG(IROOT)
	ENDDO
	WRITE(ITO,5001)                         !BLANK
      WRITE(10,5001)                         !BLANK
C	---------------------------------------------


	IDOF(1:9) = 0
	DO I = 1,9
	CALL INTFILL('%DOB',IK,1,I,0)
	IF(IK.EQ.0) IDOF(I) = 1
	ENDDO
	IPRN(1:4) = 0
	IF(IDOF(1).NE.0.OR.IDOF(2).NE.0.OR.IDOF(3).NE.0) IPRN(1) = 1
	IF(IDOF(4).NE.0.OR.IDOF(5).NE.0.OR.IDOF(6).NE.0) IPRN(2) = 1
	IF(IDOF(8).NE.0) IPRN(3) = 1
	IF(IDOF(9).NE.0) IPRN(4) = 1

	DO 1000 IMOD= 1,NROOT

	DO IC = 1,4


	IF(IPRN(IC).NE.0) THEN

	SELECTCASE(IC)
	CASE(1) 
	COMPO = 'DIS'
	CASE(2)
	COMPO = 'ROT'
	CASE(3)
	COMPO = 'TEM'
	CASE(4)
	COMPO = 'PRS'
	ENDSELECT


      IF (KSPEC.NE.1)THEN
	CALL MODHEAD(IMOD,COMPO)
	DO ISN = 1,NSN
	NRC = ISN + NSN*(IMOD-1)  !RECORD NUMBER
	READ(NFILM,REC=NRC) EDAT(1:9)  
	CALL MODVALV(EDAT,ISN,COMPO)
	ENDDO
	CALL MODEND
	ENDIF


	ENDIF

	ENDDO

1000	CONTINUE

	DEALLOCATE(REIG)

2000	FORMAT(//20X,'EIGENVALUE ANALYSIS TYPE ---  BUCKLING ANALYSIS'/,
	1	     20X,'   MODE-NO.                      EIGENVALUE'/)
2001	FORMAT(//20X,'EIGENVALUE ANALYSIS TYPE --- FREQUENZY ANALYSIS'/,
	1	     20X,'   MODE-NO.                      PERIOD(S) '/)

3000	FORMAT(22X,I5,21X,10E15.6)


4000	FORMAT(//5X,'EIGENVALUE ANALYSIS TYPE ---  BUCKLING ANALYSIS'/,
	1	     5X,'   MODE-NO.                      EIGENVALUE'/)
4001	FORMAT(//5X,'EIGENVALUE ANALYSIS TYPE --- FREQUENZY ANALYSIS'/,
	1	     5X,'   MODE-NO.                      PERIOD(S) '/)

5000	FORMAT(7X,I5,21X,10E15.6)
5001	FORMAT(' ')

      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE MODHEAD(ISTEP,COMPO)  !FREE NODE DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------	
C     ----------------------------------------------------------------
	

	SELECTCASE(COMPO)
	CASE('DIS')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,101) ISTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,201) ISTEP  !EXCEL 
	CASE('ROT')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,102) ISTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,202) ISTEP  !EXCEL 
	CASE('TEM')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,103) ISTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,203) ISTEP  !EXCEL 
	CASE('PRS')
	IF(KEXCEL.EQ.0) WRITE (LFPRIN,104) ISTEP
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,204) ISTEP  !EXCEL 
	ENDSELECT

C     ----------------------------------------------------------------

101	FORMAT(/'Result "Mode Displacement "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames  "U" "V" "W"'/,
	3		'Values')
201	FORMAT(/'Result "Mode Displacement "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Vector',2X,'OnNodes'/,
	2		'ComponentNames,  "U", "V", "W"'/,
	3		'Values')                !Excel


102	FORMAT(/'Result "Mode Disp & Rotation "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames  "U" "V" "W" "Rx" "Ry" "Rz"'/,
	3		'Values')
202	FORMAT(/'Result "Mode Disp & Rotation "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Matrix',2X,'OnNodes'/,
	2		'ComponentNames,  "U", "V", "W", "Rx", "Ry", "Rz"'/,
	3		'Values')                !Excel


103	FORMAT(/'Result "Mode Temperature "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Temperature"'/,
	3		'Values')
203	FORMAT(/'Result "Mode Temperature "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames,  "Temperature"'/,
	3		'Values')                !Excel


104	FORMAT(/'Result "Mode Pressure "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames  "Pressure"'/,
	3		'Values')
204	FORMAT(/'Result "Mode Pressure "',2X,
	1	    '"Mode Shapes"',
	1		2X,I5,2X,'Scalar',2X,'OnNodes'/,
	2		'ComponentNames,  "Pressure"'/,
	3		'Values')                !Excel


      RETURN
      END
C

C	=====================================================================
C	=====================================================================
C	=====================================================================	
	SUBROUTINE MODVALV(EDATA,ISN,COMPO)  !FREE NODE DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 COMPO

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------
	DIMENSION EDATA(1),PDATA(10)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------


	SELECTCASE(COMPO)
	CASE('DIS')
	NC  = 3
	NC1 = 1
	NC2 = 3
	CASE('ROT')
	NC  = 6
	NC1 = 1
	NC2 = 6
	CASE('TEM')
	NC  = 1
	NC1 = 8
	NC2 = 8
	CASE('PRS')
	NC  = 1
	NC1 = 9
	NC2 = 9
	ENDSELECT


	N1 = NC1
	N2 = NC2
	PDATA(1:NC) = EDATA(N1:N2)

	IF(ABS(PDATA(1)).LT.1.0E-5) PDATA(1) = 0.0D0
	IF(ABS(PDATA(2)).LT.1.0E-5) PDATA(2) = 0.0D0
	IF(ABS(PDATA(3)).LT.1.0E-5) PDATA(3) = 0.0D0
	IF(ABS(PDATA(4)).LT.1.0E-5) PDATA(4) = 0.0D0
	IF(ABS(PDATA(5)).LT.1.0E-5) PDATA(5) = 0.0D0
	IF(ABS(PDATA(6)).LT.1.0E-5) PDATA(6) = 0.0D0

	IF(KEXCEL.EQ.0) WRITE(LFPRIN,100) ISN,(PDATA(J),J=1,NC)
	IF(KEXCEL.EQ.1) THEN                !Excel
	WRITE (LFPRIN,110) ISN,','
	DO J = 1,NC
	IF(J.LT.NC) WRITE(LFPRIN,120) PDATA(J),','
	IF(J.EQ.NC) WRITE(LFPRIN,120) PDATA(J)
	ENDDO
	ENDIF


100	FORMAT(2X,I5,10E15.6)
110   FORMAT(I5,A,\)       !EXCEL
120   FORMAT(E15.6,A,\)


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE MODEND
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     FILE INDEX FOR WRITING OUTPUT      
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
C     ----------------------------------------------------------------

	IF(KEXCEL.EQ.0) WRITE (LFPRIN,100)
	IF(KEXCEL.EQ.1) WRITE (LFPRIN,200)                 !Excel

100	FORMAT('End Values'/)
200	FORMAT(/'End Values'/)

      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================

	SUBROUTINE LINKKAK
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 METHOD,COMPO
C     ----------------------------------------------------------------
	ALLOCATABLE EDAT(:)

	CALL LOCATN('OLNK',KLNK,NOUT2,NLINK,2)  
	IF(NLINK.LE.0) RETURN

	NOUT = NOUT2/2
	ALLOCATE(EDAT(NOUT2))


	DO IGM = 1,1 !NLINK
	DO I = 1,NOUT2
	CALL RELFILL('OLNK',EDAT(I),I,IGM,0)  
	ENDDO
	WRITE(*,*) EDAT(1:2)
	ENDDO

	DEALLOCATE(EDAT)


      RETURN
      END
C


C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE PRINTOUTOFFSHORE (ILCAS,MLE,NDI,IWAVE,RK,CD,CM,WVIGHT,WNU,RAMDA,IGR,DIAM,LCS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*200 LHDSTD,LHDCOM,NAMEL
      COMMON /LHEADER/ LHDSTD(500),LHDCOM(500)
      COMMON /IPCON/ IPCONTROL
C      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
C     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
C     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C       COMMON /STFIN/ IWAVES(50,50000)
C       COMMON /STFRE/ RKS(50,50000),CDS(50,50000),CMS(50,50000),WVHIGHT(50,50000),WNUS(50,50000)
C     1                ,RAMDAS(50,50000),DIAMS(50,50000)

           
C       IWAVES(ILCAS,NUMBERELEMENT)    = IWAVE 
C       RKS(ILCAS,NUMBERELEMENT)       = RK
C       CDS(ILCAS,NUMBERELEMENT)       = CD
C       CMS(ILCAS,MLE)                 = CM
C       WVHIGHT(ILCAS,MLE)             = WVIGHT
C       WNUS(ILCAS,MLE)                = WNU 
C       RAMDAS(ILCAS,MLE)              = RAMDA
C       DIAMS(ILCAS,MLE)               = DIAM

      
      IF (IGR.EQ.1.0D0)THEN
      LENGTH = LEN_TRIM(LHDSTD(ILCAS))
      NAMEL  = LHDSTD(ILCAS)
         
           IF (IPCONTROL.EQ.0.0)THEN
           WRITE (50,1) ',',NAMEL(1:LENGTH)
1          FORMAT ('ELV',A,'"OFFSHORE DETAILS"',5X,'"',A,'"')
           WRITE  (50,2) ',',',',',',',',',',',',',',','
2          FORMAT ('NO.ELEMENT',A,'WAVE THEORY',A,'WAVE NUMBER',A,'WAVE LENGTH',A,'WAVE CREST',A,'CD',A,'CM',A,'CURRENT THEORY',A)
           IPCONTROL = 1.0D0
         ENDIF
       
       
      
      
C      IF (NDI.EQ.2.0D0)THEN  ! USER DEFINED
C      WRITE (50,1) MLE
C1     FORMAT ('   NUMBER OF ELEMENT    = ',I6)
C      IF (IWAVE.EQ.1.0)THEN
C      WRITE (50,2)
C2     FORMAT ('   WAVE THEORY          =   AIRY WAVE THEORY')
C      ELSEIF (IWAVE.EQ.2.0)THEN 
C      WRITE (50,3)
C3     FORMAT ('   WAVE THEORY          =   STOKE FIFTH ORDER THEORY')
C      ELSEIF (IWAVE.EQ.3.0)THEN
C      WRITE (50,4)
C4     FORMAT ('   WAVE THEORY          =   STREAM FUNCTION WAVE THEORY')
C      ELSEIF (IWAVE.EQ.4.0)THEN
C      WRITE (50,5)
C5     FORMAT ('   WAVE THEORY          =   CNOIDAL WAVE THEORY')
C      ELSEIF (IWAVE.EQ.5.0)THEN
C      WRITE (50,6)
C6     FORMAT ('   WAVE THEORY          =   SOLITATY WAVE THEORY')
C      ENDIF
C      WRITE (50,7) CD
C7     FORMAT ('   DRAG COEFFICIENT     = ',F8.4)
C      WRITE (50,8) CM
C8     FORMAT ('   INERTIA COEFFICIENT  = ',F8.4)
C      WRITE (50,9) RK
C9     FORMAT ('   WAVE NUMBER          = ',F8.4)
C      WRITE (50,10) RAMDA
C10    FORMAT ('   WAVE LENGTH          = ',F8.4)
C      WRITE (50,11) WNU
C11    FORMAT ('   WAVE CREST           = ',F8.4)
C      WRITE (50,100)
C100   FORMAT ('')

C      ELSEIF (IGR.EQ.6.0D0)THEN
      
C      ELSEIF (NDI.EQ.1.0D0)THEN ! AUTOMATIC
C      WRITE (50,101) MLE
C101   FORMAT ('   NUMBER OF ELEMENT    = ',I6)
C      IF (IWAVE.EQ.1.0)THEN
C      WRITE (50,102)
C102   FORMAT ('   WAVE THEORY          =   AIRY WAVE THEORY')
C      ELSEIF (IWAVE.EQ.2.0)THEN 
C      WRITE (50,103)
C103   FORMAT ('   WAVE THEORY          =   STOKE FIFTH ORDER THEORY')
C      ELSEIF (IWAVE.EQ.3.0)THEN
C      WRITE (50,104)
C104   FORMAT ('   WAVE THEORY          =   STREAM FUNCTION WAVE THEORY')
C      ELSEIF (IWAVE.EQ.4.0)THEN
C      WRITE (50,105)
C105   FORMAT ('   WAVE THEORY          =   CNOIDAL WAVE THEORY')
C      ELSEIF (IWAVE.EQ.5.0)THEN
C      WRITE (50,106)
C106   FORMAT ('   WAVE THEORY          =   SOLITATY WAVE THEORY')
C      ENDIF
C      WRITE (50,112) DIAM
C112   FORMAT ('   NORMINAL DIAMETER    = ',F8.4)
C      WRITE (50,107) CD
C107   FORMAT ('   DRAG COEFFICIENT     = ',F8.4)
C      WRITE (50,108) CM
C108   FORMAT ('   INERTIA COEFFICIENT  = ',F8.4)
C      WRITE (50,109) RK
C109   FORMAT ('   WAVE NUMBER          = ',F8.4)
C      WRITE (50,110) RAMDA
C110   FORMAT ('   WAVE LENGTH          = ',F8.4)
C      WRITE (50,111) WNU
C111   FORMAT ('   WAVE CREST           = ',F8.4)
C      WRITE (50,200)
C200   FORMAT ('')
      
      
C      ENDIF ! ENDIF NDI
      ENDIF ! ENDIF IGR
      END
!  ====================================================================================================
      SUBROUTINE PRINTOUTOFFSHORE1 (IGR,MLE,ILCAS,LWCASE,UZT,WFF,DEPTH,FX,FY)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION LWCASE(7)
      CHARACTER*200 LHDSTD,LHDCOM,NAMEL
      
      COMMON/ OFFSHOREFORCE / WAVEVELOCITYX,WAVEVELOCITYY,WAVEACCX,WAVEACCY,WAVEFORCEX,WAVEFORCEY
     1                        ,CURRENTVELOCITY,CURRENTFORCE,WAVEPLUSCURRENTX,WAVEPLUSCURRENTY
     1                        ,WAVEPLUSCURRENTZ,COORDINATEX,COORDINATEY,FDX
      COMMON /OUTCONTROL / NCONTROL,ICONTROL,JCONTROL,NCCONTROL,NWCONTROL
      COMMON /OUTCONTROL1/ XCONTROL,XCCONTROL,XWCONTROL
      
      COMMON /LHEADER/ LHDSTD(500),LHDCOM(500)
      
       LENGTH = LEN_TRIM(LHDSTD(ILCAS))
       NAMEL  = LHDSTD(ILCAS)
       
       ! BASE ON X,Y,Z AXIS
       IF (ICONTROL.EQ.0.0)THEN
       ICONTROL = ILCAS
       JCONTROL = 2.0D0
       ELSEIF (ILCAS-ICONTROL.NE.0.0)THEN
       ICONTROL = ILCAS
       JCONTROL = 0.0D0
       ELSEIF (ILCAS-ICONTROL.EQ.0.0)THEN
       ICONTROL = ILCAS
       JCONTROL = 1.0D0
       ENDIF
       
       IF (JCONTROL.EQ.1.0D0)THEN
       GOTO 1001
       ELSEIF (JCONTROL.EQ.0.0D0)THEN
       GOTO 1002
       ELSEIF (JCONTROL.EQ.2.0D0)THEN
       GOTO 1001
       ENDIF
       
       
1002   DO I = 1,4
       WRITE (51,1003)
       WRITE (52,1003)
       WRITE (53,1003)
1003   FORMAT ('')
       ENDDO
       
       IF (JCONTROL.EQ.0.0D0)THEN
       NCONTROL   = 0.0D0
       NCCONTROL = 0.0D0
       NWCONTROL  = 0.0D0
       XCONTROL   = 0.0D0
       ENDIF
       
       GOTO 1001
    
1001   IF (LWCASE(1).EQ.1.0D0)THEN ! WAVE FORCE
       
         IF (NCONTROL.EQ.0.0)THEN
           WRITE (51,1) ',',NAMEL(1:LENGTH)
1          FORMAT ('ELV',A,'"WAVE DETAILS"',5X,'"',A,'"')
           WRITE  (51,2) ',',',',',',',',',',',',',',','
2          FORMAT ('X',A,'Y',A,'UX',A,'UY',A,'AX',A,'AY',A,'FX',A,'FY',A)
           NCONTROL = 1.0D0
         ENDIF
         
      IF (COORDINATEY.EQ.0.0D0)THEN
      WRITE (51,101) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',WAVEVELOCITYY,',',WAVEACCX,',',WAVEACCY,',',WAVEFORCEX,',',
     1              WAVEFORCEY,','
101   FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A)
      ENDIF
      
      IF (COORDINATEY-XCONTROL.GT.0.0D0.OR.COORDINATEY-XCONTROL.LT.0.0D0)THEN
      XCONTROL = COORDINATEY
      
      WRITE (51,102) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',WAVEVELOCITYY,',',WAVEACCX,',',WAVEACCY,',',WAVEFORCEX,',',
     1              WAVEFORCEY,','
102   FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A)
      ENDIF
      ENDIF
        
      
      IF (LWCASE(2).EQ.1.0D0)THEN ! CURRENT
         IF (NCCONTROL.EQ.0.0)THEN
           WRITE (52,4) ',',NAMEL(1:LENGTH)
4          FORMAT ('ELV',A,'"CURRENT DETAILS"',5X,'"',A,'"')
           WRITE  (52,5) ',',',',',',','
5          FORMAT ('X',A,'Y',A,'UX',A,'FX',A)
           NCCONTROL = 1.0D0
         ENDIF
      
      IF (COORDINATEY.EQ.0.0D0)THEN
      WRITE (52,103) COORDINATEX,',',COORDINATEY,',',CURRENTVELOCITY,',',CURRENTFORCE
103   FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
      ENDIF
      
      IF (COORDINATEY-XCCONTROL.GT.0.0D0.OR.COORDINATEY-XCCONTROL.LT.0.0D0)THEN
      XCCONTROL = COORDINATEY
      
      WRITE (52,104) COORDINATEX,',',COORDINATEY,',',CURRENTVELOCITY,',',CURRENTFORCE
104   FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
      ENDIF
      
      ENDIF
      
      IF (LWCASE(4).EQ.1.0D0)THEN ! TIDAL LOAD
C     ========== UNNECESSARY BY TOEY ==========
      ENDIF
      
      IF (LWCASE(5).EQ.1.0D0)THEN ! WAVE BREAKING LOAD  PLUNGING  

         IF (NWCONTROL.EQ.0.0)THEN
           WRITE (54,9) ',',NAMEL(1:LENGTH)
9          FORMAT ('ELV',A,'"WAVE BREAKLING DETAILS ( PLUNGING )"',5X,'"',A,'"')
           WRITE  (54,10) ',',',',',',','
10          FORMAT ('X',A,'Y',A,'UX',A,'FX',A)
           NWCONTROL = 1.0D0
         ENDIF
         
      IF (COORDINATEY.GT.DEPTH)THEN
          
          IF (COORDINATEY.EQ.0.0D0)THEN
          WRITE (54,109) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',FDX
109       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
      
          IF (COORDINATEY-XCCONTROL.GT.0.0D0.OR.COORDINATEY-XCCONTROL.LT.0.0D0)THEN
          XCCONTROL = COORDINATEY
          WRITE (54,110) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',FDX
110       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
      
      ELSEIF (COORDINATEY.LE.DEPTH)THEN
      
          IF (COORDINATEY.EQ.0.0D0)THEN
          WRITE (54,107) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',WAVEFORCEX
107       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
      
          IF (COORDINATEY-XCCONTROL.GT.0.0D0.OR.COORDINATEY-XCCONTROL.LT.0.0D0)THEN
          XCCONTROL = COORDINATEY
          WRITE (54,108) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',WAVEFORCEX
108       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
          
      ENDIF
         
         

      ENDIF
      IF (LWCASE(6).EQ.1.0D0)THEN ! WAVE BREAKING LOAD  SURGING
C     ========== UNNECESSARY BY TOEY ==========      
         
         IF (NWCONTROL.EQ.0.0)THEN
           WRITE (55,11) ',',NAMEL(1:LENGTH)
11          FORMAT ('ELV',A,'"WAVE BREAKLING DETAILS ( SURGING )"',5X,'"',A,'"')
           WRITE  (55,12) ',',',',',',','
12          FORMAT ('X',A,'Y',A,'UX',A,'FX',A)
           NWCONTROL = 1.0D0
         ENDIF
         
      IF (COORDINATEY.GT.DEPTH)THEN
          
          IF (COORDINATEY.EQ.0.0D0)THEN
          WRITE (55,111) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',WAVEFORCEX
111       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
      
          IF (COORDINATEY-XCCONTROL.GT.0.0D0.OR.COORDINATEY-XCCONTROL.LT.0.0D0)THEN
          XCCONTROL = COORDINATEY
          WRITE (55,112) COORDINATEX,',',COORDINATEY,',',WAVEVELOCITYX,',',WAVEFORCEX
112       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
      
      ELSEIF (COORDINATEY.LE.DEPTH)THEN
      
          IF (COORDINATEY.EQ.0.0D0)THEN
          WRITE (55,113) COORDINATEX,',',COORDINATEY,',',UBW,',',FDXS
113       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
      
          IF (COORDINATEY-XCCONTROL.GT.0.0D0.OR.COORDINATEY-XCCONTROL.LT.0.0D0)THEN
          XCCONTROL = COORDINATEY
          WRITE (55,114) COORDINATEX,',',COORDINATEY,',',UBW,',',FDXS
114       FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
          ENDIF
          
      ENDIF


      ENDIF
      
      
      IF (LWCASE(7).EQ.1.0D0)THEN ! WIND LOAD 
 
         IF (NWCONTROL.EQ.0.0)THEN
           WRITE (53,7) ',',NAMEL(1:LENGTH)
7          FORMAT ('ELV',A,'"WIND DETAILS"',5X,'"',A,'"')
           WRITE  (53,8) ',',',',',',','
8          FORMAT ('X',A,'Y',A,'UX',A,'FX',A)
           NWCONTROL = 1.0D0
         ENDIF
      
      IF (COORDINATEY.EQ.0.0D0)THEN
      WRITE (53,105) COORDINATEX,',',COORDINATEY,',',UZT,',',WFF
105   FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
      ENDIF
      
      IF (COORDINATEY-XWCONTROL.GT.0.0D0.OR.COORDINATEY-XWCONTROL.LT.0.0D0)THEN
      XWCONTROL = COORDINATEY
      
      WRITE (53,106) COORDINATEX,',',COORDINATEY,',',UZT,',',WFF
106   FORMAT (F9.4,A,F9.4,A,E12.5,A,E12.5,A)
      ENDIF
      
      ENDIF
      END
!  ====================================================================================================
      SUBROUTINE PRINTOUTOFFSHORE2 (IGR,MLE,ILCAS,LWCASE,UZT,WFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      END
!  ====================================================================================================     
      SUBROUTINE MAXVALUE(PDATA,ISN,NV)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION PDATA(6)
      COMMON /SOD/ ANSX(500000),ANSY(500000),ANSZ(500000),ANSXY(500000)
     1            ,ANSXZ(500000),ANSYZ(500000) 
      
      ANSX(ISN)  = PDATA(1) 
      
      IF (NV.EQ.6)THEN
      ANSY(ISN)  = PDATA(2) 
      ANSZ(ISN)  = PDATA(3) 
      ANSXY(ISN) = PDATA(4) 
      ANSXZ(ISN) = PDATA(5) 
      ANSYZ(ISN) = PDATA(6)
      ENDIF 
      
      END
!  ====================================================================================================    
!  ====================================================================================================    
      SUBROUTINE DEFINDLOCALSTRESS_GID (NELEMENT,OPT,HED,NAMEHEAD,NAML,JSTEP) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPT
      CHARACTER*200 NAMEHEAD,NAMEG
	CHARACTER*3 METHOD
	CHARACTER(3)  HED
      COMMON /XFPOUT/ KEXCEL,LFPRIN,LFPRN1,LFPRN2,LFPRN3,LFPRN4
      COMMON / FRAMEFORCE / AR(30),NELEMENTT
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      COMMON /LOCAL_STRESS/ AXIAL(5,10000),SHEARS(5,10000),SHEART(5,10000),ATORSION(5,10000),AMOMENTS(5,10000),AMOMENTT(5,10000)
      
      
      IF (OPT.EQ."HEAD")THEN
      WRITE (LFPRIN,1000) HED,NAMEHEAD(1:NAML),JSTEP
      
      ELSEIF (OPT.EQ."CLEA")THEN
                
      AXIAL    = 0.
      SHEARS   = 0.
      SHEART   = 0.
      ATORSION = 0.
      AMOMENTS = 0.
      AMOMENTT = 0.
      
      ELSEIF (OPT.EQ."BODY")THEN
          
       CALL ARRAYELEMENT (NELEMENT,NSECTION)
       
       AXIAL(1,NELEMENT)    = AR(1)  /AREA
       AXIAL(2,NELEMENT)    = AR(7)  /AREA
       AXIAL(3,NELEMENT)    = AR(13) /AREA
       AXIAL(4,NELEMENT)    = AR(19) /AREA
       AXIAL(5,NELEMENT)    = AR(25) /AREA
       
       SHEARS(1,NELEMENT)   = 2D0*AR(2)  /AREA
       SHEARS(2,NELEMENT)   = 2D0*AR(8)  /AREA
       SHEARS(3,NELEMENT)   = 2D0*AR(14) /AREA
       SHEARS(4,NELEMENT)   = 2D0*AR(20) /AREA
       SHEARS(5,NELEMENT)   = 2D0*AR(26) /AREA
      
       SHEART(1,NELEMENT)   = 2D0*AR(3)  /AREA
       SHEART(2,NELEMENT)   = 2D0*AR(9)  /AREA
       SHEART(3,NELEMENT)   = 2D0*AR(15) /AREA
       SHEART(4,NELEMENT)   = 2D0*AR(21) /AREA
       SHEART(5,NELEMENT)   = 2D0*AR(27) /AREA
       
       IF (NSECTION.EQ.21) C = DPIPE  /2D0
       IF (NSECTION.EQ.22) C = DROUND /2D0
       IF (NSECTION.EQ.23) C = DROUND /2D0
      
       ATORSION(1,NELEMENT) = AR(4) *C /AJ
       ATORSION(2,NELEMENT) = AR(10)*C /AJ
       ATORSION(3,NELEMENT) = AR(16)*C /AJ
       ATORSION(4,NELEMENT) = AR(22)*C /AJ
       ATORSION(5,NELEMENT) = AR(28)*C /AJ
       
       AMOMENTS(1,NELEMENT) = AR(5)  /SECTIONT
       AMOMENTS(2,NELEMENT) = AR(11) /SECTIONT
       AMOMENTS(3,NELEMENT) = AR(17) /SECTIONT
       AMOMENTS(4,NELEMENT) = AR(23) /SECTIONT
       AMOMENTS(5,NELEMENT) = AR(29) /SECTIONT
      
       AMOMENTT(1,NELEMENT) = AR(6)  /SECTIONT
       AMOMENTT(2,NELEMENT) = AR(12) /SECTIONT
       AMOMENTT(3,NELEMENT) = AR(18) /SECTIONT
       AMOMENTT(4,NELEMENT) = AR(24) /SECTIONT
       AMOMENTT(5,NELEMENT) = AR(30) /SECTIONT
       
      ELSEIF (OPT.EQ."WRIT")THEN
          
      !DO NELEMENT = 1,NELE
       DO I = 1,5
         IF (I.EQ.1) WRITE (LFPRIN,101) NELEMENT,AXIAL(I,NELEMENT),SHEARS(I,NELEMENT),SHEART(I,NELEMENT)
     1                     ,ATORSION(I,NELEMENT),AMOMENTS(I,NELEMENT),AMOMENTT(I,NELEMENT)
         IF (I.NE.1) WRITE (LFPRIN,100) AXIAL(I,NELEMENT),SHEARS(I,NELEMENT),SHEART(I,NELEMENT)
     1                     ,ATORSION(I,NELEMENT),AMOMENTS(I,NELEMENT),AMOMENTT(I,NELEMENT)
       ENDDO
       !ENDDO
       
       !WRITE (LFPRIN,102)
      ENDIF
      
100   FORMAT (7X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)      
101   FORMAT (I5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)
!102   FORMAT ("End Values")
      
1000	FORMAT(/'Result "XFrame-Local-Stress ',A3,'"',2X,
	1	    '"',A,'"',
	1		2X,I5,2X,'Matrix',2X,'OnGaussPoints',2X,'"XFrame"'/,
	2		'ComponentNames',
     3		' "Axial-Stress" "Shear Stress-S" "Shear Stress-T"',
	4		' "Torsional Stress" "Bending Stress-S" "Bending Stress-T"'/,
	5		'Values')
      
      END
      