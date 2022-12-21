      Subroutine  Wave_Diameter(MLE,RadiousOUTCAL)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      
      ! OUTPUT FOR CALCULATION
     	COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREAxx,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS

      Call ARRAYELEMENT (MLE,NSECTION)
      
      RadiousOUTCAL = 0.5d0 * DPIPE
      
      Return
      End Subroutine
      !===============================================================================================