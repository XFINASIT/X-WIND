C	===================================================================================       
      SUBROUTINE Solitary_Wave (H,waterdepth,g,Atime,XX,z,u,v,au,av,Surface_Elevation,q)

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
                  
      IOPTION = 2 ! IOPTION = 1 ; Original
                  ! IOPTION = 2 ; Modify      
      

      d=waterdepth
      E=H/d
      
      ! Wave celerity,=> c
      c=((g*d)**0.5d0)*(1d0+0.5d0*E-0.15d0*(E**2d0))
    
      ! The particle velocity at the crest relative to the crest,=> q
      q=(sqrt(3d0*E))/d*(1d0-0.625d0*E)*(XX-c*Atime)
    
      ! Surface_Elevation
      Surface_Elevation= d * ( E * (1/((cosh(q))**2d0)) - 0.75d0* (E**2d0) * (1/((cosh(q))**2d0)) * tanh(q) )
      
      If( IOPTION .EQ. 1 ) THEN
          
     
      ! Horizontal particle velocity,=> u
      u =(sqrt(g*d))*(  E*(1/((cosh(q))**2d0)) + (E**2d0)*(1/((cosh(q))**2d0))*(0.25d0- (1/((cosh(q))**2d0))
     1-0.75d0*((z/d)**2d0)*(2d0-3d0*(1/((cosh(q))**2d0)) ) ) )
    
      ! Vertical particle velocity,=> v
      v =(sqrt(g*d))*E*(sqrt(3d0*E))*(z/d)*(1/((cosh(q))**2d0))* tanh(q)*(1d0-E*(0.375d0-2d0*(1/((cosh(q))**2d0))
     1     -0.5d0*((z/d)**2d0)*(1d0+3d0*(1/((cosh(q))**2d0))   )))
    
      ! Horizontal particle acceleration ,=> au
      au=g*(sqrt(3d0*E))*(1/((cosh(q))**2d0))* tanh(q) *(1d0+E*(0.125d0-2d0*(1/((cosh(q))**2d0))
     1    -1.5d0*((z/d)**2d0)*(1d0+3d0*(1/((cosh(q))**2d0)) )))
     
      ! Vertical particle acceleration ,=> av
      av=g*(sqrt(3d0*E))*(1.5d0*(E**2d0))*(z/d)*(1/((cosh(q))**2d0))*(2d0-3d0*(1/((cosh(q))**2d0)))
      
      
      
      ELSEIf( IOPTION .EQ. 2 ) THEN
                    
      SHAPEFUNCTION = 0.57D0
          
      ! Horizontal particle velocity,=> u
      u =1.031D0 * (sqrt(g*d))*(  E*(1/((cosh(q))**2d0)) + (E**2d0)*(1/((cosh(q))**2d0))*(0.25d0- (1/((cosh(q))**2d0))
     1-  SHAPEFUNCTION *0.75d0*((z/d)**2d0)*(2d0-3d0*(1/((cosh(q))**2d0)) ) ) )
    
      ! Vertical particle velocity,=> v
      v =(sqrt(g*d))*E*(sqrt(3d0*E))*(z/d)*(1/((cosh(q))**2d0))* tanh(q)*(1d0-E*(0.375d0-2d0*(1/((cosh(q))**2d0))
     1     -0.5d0*((z/d)**2d0)*(1d0+3d0*(1/((cosh(q))**2d0))   )))
    
      ! Horizontal particle acceleration ,=> au
      au=g*(sqrt(3d0*E))*(1/((cosh(q))**2d0))* tanh(q) *(1d0+E*(0.125d0-2d0*(1/((cosh(q))**2d0))
     1    -1.5d0*((z/d)**2d0)*(1d0+3d0*(1/((cosh(q))**2d0)) )))
     
      ! Vertical particle acceleration ,=> av
      av=g*(sqrt(3d0*E))*(1.5d0*(E**2d0))*(z/d)*(1/((cosh(q))**2d0))*(2d0-3d0*(1/((cosh(q))**2d0)))    
      
      
          
      ENDIF

      RETURN
      END SUBROUTINE
C	=========================================================================================
C	========================================================================================= 
