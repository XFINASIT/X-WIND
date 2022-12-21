      SUBROUTINE CLOSE_FILE
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      CLOSE(UNIT=810,IOSTAT= IOS,status='DELETE') ! 
      CLOSE(UNIT=811,IOSTAT= IOS,status='DELETE') !
      
      CLOSE(UNIT=70,IOSTAT= IOS,status='DELETE') ! MARINE GROWTH STORAGE
      CLOSE(UNIT=71,IOSTAT= IOS,status='DELETE') ! CURRENT DATA STORAGE
      
      CLOSE(UNIT=106,IOSTAT= IOS,STATUS='DELETE') !'Out Displacement.csv' 
      CLOSE(UNIT=107,IOSTAT= IOS,STATUS='DELETE') !'Out Element.csv'      
      CLOSE(UNIT=108,IOSTAT= IOS,STATUS='DELETE') !'Out Support.csv'      
      CLOSE(UNIT=109,IOSTAT= IOS,STATUS='DELETE') !'Out Link.csv'         
      CLOSE(UNIT=120,IOSTAT= IOS,STATUS='DELETE') !'Out Stress.csv'           
      CLOSE(UNIT=170,IOSTAT= IOS,STATUS='DELETE') !'Out So il.dat'         
      
      CLOSE(UNIT=201,IOSTAT= IOS,status='DELETE') ! DYNAMIC_MASS.DAT
      CLOSE(UNIT=220,IOSTAT= IOS,status='DELETE') ! DYNAMIC_NODAL.DAT
      CLOSE(UNIT=230,IOSTAT= IOS,status='DELETE') ! ASCII.BIN
      
      CLOSE(UNIT=4900,IOSTAT= IOS,status='DELETE') ! 4900 HYDRO SIFFNESS.dat
      CLOSE(UNIT=4901,IOSTAT= IOS,status='DELETE') ! 4901 HYDRO SIFFNESS ERROR.dat
      CLOSE(UNIT=4902,IOSTAT= IOS,status='DELETE') ! 4902 HYDRO SIFFNESS ORI.dat
        
      CLOSE(UNIT=5061,IOSTAT= IOS,STATUS='DELETE') ! 'GH_BLADED_JOINT_ANALYSIS.txt'      
      
      CLOSE(UNIT=5299,IOSTAT= IOS,status='DELETE') ! 5299 WIND.TXT
      CLOSE(UNIT=5300,IOSTAT= IOS,status='DELETE') ! 5300 WAVE.TXT
      CLOSE(UNIT=5301,IOSTAT= IOS,status='DELETE') ! 5300 FORCE.TXT
      
      CLOSE(UNIT=7100,IOSTAT= IOS,status='DELETE') ! Stress on element
      
      
      CLOSE(UNIT=57,IOSTAT= IOS,STATUS='DELETE')   ! 'Result of ocean analysis/Solid_WAVE_PRESSURE.out'       
      CLOSE(UNIT=58,IOSTAT= IOS,STATUS='DELETE')   ! 'Result of ocean analysis/Solid_Diffraction_PRESSURE.out'
      CLOSE(UNIT=59,IOSTAT= IOS,STATUS='DELETE')   ! 'Result of ocean analysis/Soild_Morison_PRESSURE.out'                  
      CLOSE(UNIT=5026,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5026 Stoke_wave_surface.out'            
      CLOSE(UNIT=5027,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5027 WAVE_FORCE_ELEVATION.out'          
      CLOSE(UNIT=5028,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5028 ELEMENT_COORDINATES.out'           
      CLOSE(UNIT=5029,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5029 TRAPIZOIDAL_RULE_FORCES.out'        
      CLOSE(UNIT=5031,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis5031 TRAPIZOIDAL_RULE_STRIP_FORCES.out'
      CLOSE(UNIT=5032,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5032 STRIP_FORCES.out'                 
      CLOSE(UNIT=5033,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5033 TRASFORMATION_STRIP_FORCES.out'                 
      CLOSE(UNIT=5038,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5038 JACOBIAN ELLIPTIC FUNCTION.out'   
      CLOSE(UNIT=5039,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5039 JACOBIAN ELLIPTIC INTEGRAL.out'    
      CLOSE(UNIT=5040,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5040 CNOIDAL SURFACE ELEVATION.out'     
      CLOSE(UNIT=5041,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5041 WAVE SPECTRUM.out'                
      CLOSE(UNIT=5042,IOSTAT= IOS,STATUS='DELETE') ! 'Result of ocean analysis/5042 STRIFF WAVE FORCE.out'    
      
      CLOSE(UNIT=7110,IOSTAT= IOS,STATUS='DELETE')  
      
      ! CLOSE INTERNAL STORAGE
      !OPEN(UNIT=999,FILE='fort.112',STATUS='UNKNOWN',IOSTAT= IOS)
      CLOSE(UNIT=112,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=299,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=791,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=989,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=990,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=4001,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=4004,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=5203,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=5204,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=7101,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=7102,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=8001,IOSTAT= IOS,STATUS='DELETE') 
      CLOSE(UNIT=8002,IOSTAT= IOS,STATUS='DELETE') 
      END