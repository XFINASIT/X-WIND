
C=====================================================================
      SUBROUTINE TRICB3_C(RS,SLN,RN,SN,BL,BLR,BLS,FLR,FLS,FLRR,FLSS,
	1                    DWR,DWRS,DWS,DWSR,CB)
	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cb FOR BENDING STRAINS
C
C	INPUT VARIABLES
C	RS(2,3)			= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)			= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	    = ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	BL(3),BLR,BLS	= SHEAR DEFORMATION PARAMETERS ALONG THE BOUNDARY,
C						R, S, RESP.
C	FLR,FLS,
C	DWR,DWS			= COMMON FACTORS DEFINED IN SUBROUTINE COMFACT
C
C	OUTPUT VARIAVBLES
C	CB(9,18)		= BENDING STRAIN PARAMETER
C	--------------------------------------------------------------
      DIMENSION RS(2,3),SLN(3),RN(3),SN(3),BL(3)
	DIMENSION FLR(3),FLS(3),FLRR(3),FLSS(3)
      DIMENSION DWR(9),DWS(9),DWRS(9),DWSR(9)
	DIMENSION CB(9,18)

C	--------------------------------
C	SOLVE FOR Cb FOR BENDING STRAINS
C	--------------------------------
	CB = 0.0D0

C	FIRST ROW
C	------------------
	CB(1,4) =
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2)
	CB(1,10)=
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(2)*RN(2)**2*SN(2)*BL(2)/2)
	CB(1,16)=
     #(-SLN(2)*RN(2)**2*SN(2)*BL(2)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2)

	CB(1,5) =
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))*SLN(1)+
     #(RN(3)**3/2+(1.E0/2.E0-BL(3)/2)*SN(3)**2*RN(3))*SLN(3))
	CB(1,11)=
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))*SLN(1)+
     #(RN(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2))
	CB(1,17)=
     #((RN(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2)+
     #(RN(3)**3/2+(-BL(3)/2+1.E0/2.E0)*SN(3)**2*RN(3))*SLN(3))

	CB(1,3) =(RN(3)*SN(3)*BL(3)-RN(1)*SN(1)*BL(1))
	CB(1,9) =(RN(1)*SN(1)*BL(1)-RN(2)*SN(2)*BL(2))
	CB(1,15)=(RN(2)*SN(2)*BL(2)-RN(3)*SN(3)*BL(3))

C	SECOND ROW
C	------------------
	CB(2,4) =
     #((-RS(1,1)/4-RS(1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,1)
     #/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)
     #-DWR(1)
	CB(2,10)=
     #((-RS(1,1)/4-RS(
     #1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,2)/4-RS(1,3)/4)*BL(2)*S
     #LN(2)*SN(2)*RN(2)**2)
     #-DWR(2)
	CB(2,16)=
     #((-RS(1,2)/4-RS(1,3)/4)*BL(2)*SLN(2)*
     #SN(2)*RN(2)**2+(-RS(1,1)/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)
     #-DWR(3)

	CB(2,5) =
     #((RS(1,2)/6+RS(1,1)/3)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2
     #)/4)*BL(1)+RS(1,2)/6+RS(1,1)/3)*SLN(1)*SN(1)**2*RN(1)+(RS(1,1)/3+R
     #S(1,3)/6)*SLN(3)*RN(3)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/3+
     #RS(1,3)/6)*SLN(3)*SN(3)**2*RN(3))
     #-DWR(4)-FLS(1)
	CB(2,11)=
     #((RS(1,1)/6+RS(1,2)/3)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2
     #)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(1)**2*RN(1)+(RS(1,2)/3+R
     #S(1,3)/6)*SLN(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/3+
     #RS(1,3)/6)*SLN(2)*SN(2)**2*RN(2))
     #-DWR(5)-FLS(2)
	CB(2,17)=
     #((RS(1,2)/6+RS(1,3)/3)*SL
     #N(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/3)*S
     #LN(2)*SN(2)**2*RN(2)+(RS(1,1)/6+RS(1,3)/3)*SLN(3)*RN(3)**3+((-RS(1
     #,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/6+RS(1,3)/3)*SLN(3)*SN(3)**2*RN(3))
     #-DWR(6)-FLS(3)

	CB(2,3) =((-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)*RN(1)+(RS(1,1)/2+RS(1,
     #3)/2)*BL(3)*SN(3)*RN(3))
     #-DWR(7)
	CB(2,9) =((RS(1,1)/2+RS(1,2)/2)*BL(1)*SN(1)*
     #RN(1)+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)*RN(2))
     #-DWR(8)
	CB(2,15)=((RS(1,3)/2+
     #RS(1,2)/2)*BL(2)*SN(2)*RN(2)+(-RS(1,3)/2-RS(1,1)/2)*BL(3)*SN(3)*RN
     #(3))
     #-DWR(9)

C	THIRD ROW
C	------------------
	CB(3,4) =
     #((-RS(2,2)/4-RS(2,1)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,1)
     #/4-RS(2,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)
	CB(3,10)=
     #((-RS(2,2)/4-RS(
     #2,1)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,2)/4-RS(2,3)/4)*BL(2)*S
     #LN(2)*SN(2)*RN(2)**2)
	CB(3,16)=
     #((-RS(2,2)/4-RS(2,3)/4)*BL(2)*SLN(2)*
     #SN(2)*RN(2)**2+(-RS(2,1)/4-RS(2,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)

	CB(3,5) =
     #((RS(2,1)/3+RS(2,2)/6)*SLN(1)*RN(1)**3+((-RS(2,2)/4-RS(2,1
     #)/4)*BL(1)+RS(2,1)/3+RS(2,2)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/6+R
     #S(2,1)/3)*SLN(3)*RN(3)**3+((-RS(2,1)/4-RS(2,3)/4)*BL(3)+RS(2,3)/6+
     #RS(2,1)/3)*SLN(3)*SN(3)**2*RN(3))
 	CB(3,11)=
     #((RS(2,2)/3+RS(2,1)/6)*SLN(1)*RN(1)**3+((-RS(2,2)/4-RS(2,1
     #)/4)*BL(1)+RS(2,2)/3+RS(2,1)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/6+R
     #S(2,2)/3)*SLN(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,3)/6+
     #RS(2,2)/3)*SLN(2)*SN(2)**2*RN(2))
	CB(3,17)=
     #((RS(2,2)/6+RS(2,3)/3)*SL
     #N(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,2)/6+RS(2,3)/3)*S
     #LN(2)*SN(2)**2*RN(2)+(RS(2,1)/6+RS(2,3)/3)*SLN(3)*RN(3)**3+((-RS(2
     #,1)/4-RS(2,3)/4)*BL(3)+RS(2,1)/6+RS(2,3)/3)*SLN(3)*SN(3)**2*RN(3))

	CB(3,3) =
     #((-RS(2,1)/2-RS(2,2)/2)*BL(1)*SN(1)*RN(1)+(RS(2,1)/2+RS(2,
     #3)/2)*BL(3)*SN(3)*RN(3))
	CB(3,9) =
     #((RS(2,1)/2+RS(2,2)/2)*BL(1)*SN(1)*
     #RN(1)+(-RS(2,2)/2-RS(2,3)/2)*BL(2)*SN(2)*RN(2))
	CB(3,15)=
     #((RS(2,2)/2+
     #RS(2,3)/2)*BL(2)*SN(2)*RN(2)+(-RS(2,3)/2-RS(2,1)/2)*BL(3)*SN(3)*RN
     #(3))

C	FOURTH ROW
C	------------------
	CB(4,4) =
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(BL(3)/2-1.E0/2.E0
     #)*SLN(3)*SN(3)*RN(3)**2-SLN(1)*SN(1)**3/2-SLN(3)*SN(3)**3/2)
	CB(4,10)=
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(-1.E0/2.E0+BL(2)/2)*
     #SLN(2)*SN(2)*RN(2)**2-SLN(1)*SN(1)**3/2-SLN(2)*SN(2)**3/2)
	CB(4,16)=
     #((-1.E0/2.E0+BL(2)/2)*SLN(2)*SN(2)*RN(2)**2+(BL(3)/2-1.E0/2.E0)*SL
     #N(3)*SN(3)*RN(3)**2-SLN(2)*SN(2)**3/2-SLN(3)*SN(3)**3/2)

	CB(4,5) =
     #(SLN(1)*SN(1)**2*RN(1)*BL(1)/2+SLN(3)*SN(3)**2*RN(3)*BL(3)/2)
	CB(4,11)=
     #(SLN(2)*SN(2)**2*RN(2)*BL(2)/2+SLN(1)*SN(1)**2*RN(1)*BL(1)/2)
	CB(4,17)=
     #(SLN(3)*SN(3)**2*RN(3)*BL(3)/2+SLN(2)*SN(2)**2*RN(2)*BL(2)/2)

	CB(4,3) =(SN(1)*RN(1)*BL(1)-SN(3)*RN(3)*BL(3))
	CB(4,9) =(SN(2)*RN(2)*BL(2)-SN(1)*RN(1)*BL(1))
	CB(4,15)=(SN(3)*RN(3)*BL(3)-SN(2)*RN(2)*BL(2))

C	FIFTH ROW
C	------------------
	CB(5,4) =
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1)-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,3)/4+RS(1,1)/4)*BL(3)-RS(1,1)/3-RS(1,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN(1)**3+(-RS(1,1)
     #/3-RS(1,3)/6)*SLN(3)*SN(3)**3)
	CB(5,10)=
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1
     #)-RS(1,2)/3-RS(1,1)/6)*SLN(1)*SN(1)*RN(1)**2+((RS(1,3)/4+RS(1,2)/4
     #)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(1,2)/3-RS(
     #1,1)/6)*SLN(1)*SN(1)**3+(-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)**3)
	CB(5,16)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(1,3)/4+RS(1,1)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(1,3)/3-RS(1,2)/6)*SLN(2)*SN(2)**3+(-RS(1
     #,1)/6-RS(1,3)/3)*SLN(3)*SN(3)**3)

	CB(5,5) =
     #((RS(1,2)/4+RS(1,1)/4)*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,3)/4+RS(1,1)/4)*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))
	CB(5,11)=
     #((RS(1,2)/4+RS(1,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,3
     #)/4+RS(1,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))
	CB(5,17)=
     #((RS(1,3)/4+RS(
     #1,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(1,3)/4+RS(1,1)/4)*BL(3)*SL
     #N(3)*SN(3)**2*RN(3))

	CB(5,3) =
     #((RS(1,1)/2+RS(1,2)/2)*BL(1)*SN(1)*RN(
     #1)+(-RS(1,1)/2-RS(1,3)/2)*BL(3)*SN(3)*RN(3))
	CB(5,9) =
     #((-RS(1,2)/2-RS
     #(1,1)/2)*BL(1)*SN(1)*RN(1)+(RS(1,2)/2+RS(1,3)/2)*BL(2)*SN(2)*RN(2)
     #)
	CB(5,15)=
     #((-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)*RN(2)+(RS(1,3)/2+RS(1,
     #1)/2)*BL(3)*SN(3)*RN(3))

C	SIXTH ROW
C	------------------
	CB(6,4) =
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,1)/4+RS(2,3)/4)*BL(3)-RS(2,3)/6-RS(2,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(1)**3+(-RS(2,3)
     #/6-RS(2,1)/3)*SLN(3)*SN(3)**3)
     #-DWS(1)-FLR(1)
	CB(6,10)=
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1
     #)-RS(2,2)/3-RS(2,1)/6)*SLN(1)*SN(1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4
     #)*BL(2)-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(2,2)/3-RS(
     #2,1)/6)*SLN(1)*SN(1)**3+(-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)**3)
     #-DWS(2)-FLR(2)
	CB(6,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(2,1)/4+RS(2,3)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN(2)**3+(-RS(2
     #,1)/6-RS(2,3)/3)*SLN(3)*SN(3)**3)
     #-DWS(3)-FLR(3)

	CB(6,5) =
     #((RS(2,2)/4+RS(2,1)/4)*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,1)/4+RS(2,3)/4)*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))
     #-DWS(4)
	CB(6,11)=
     #((RS(2,2)/4+RS(2,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3
     #)/4+RS(2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))
     #-DWS(5)
	CB(6,17)=
     #((RS(2,3)/4+RS(
     #2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(2,1)/4+RS(2,3)/4)*BL(3)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWS(6)

	CB(6,3) =
     #((RS(2,1)/2+RS(2,2)/2)*BL(1)*SN(1)*RN(
     #1)+(-RS(2,1)/2-RS(2,3)/2)*BL(3)*SN(3)*RN(3))
     #-DWS(7)
	CB(6,9) =
     #((-RS(2,1)/2-RS
     #(2,2)/2)*BL(1)*SN(1)*RN(1)+(RS(2,3)/2+RS(2,2)/2)*BL(2)*SN(2)*RN(2)
     #)
     #-DWS(8)
	CB(6,15)=
     #((-RS(2,3)/2-RS(2,2)/2)*BL(2)*SN(2)*RN(2)+(RS(2,3)/2+RS(2,
     #1)/2)*BL(3)*SN(3)*RN(3))
     #-DWS(9)

C	SEVENTH ROW
C	------------------
	CB(7,4) =
     #((-1.E0/2.E0+BL(1)/2)*SLN(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SL
     #N(1)*SN(1)**2*RN(1)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-BL(3)/2-
     #1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))
	CB(7,10)=((-1.E0/2.E0+BL(1)/2)*SLN
     #(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SLN(1)*SN(1)**2*RN(1)+(-1.E0/2.E
     #0+BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2.E0)*SLN(2)*SN(2)**2*RN
     #(2))
	CB(7,16)=((-1.E0/2.E0+BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2
     #.E0)*SLN(2)*SN(2)**2*RN(2)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-B
     #L(3)/2-1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))

	CB(7,5) =((BL(1)/2+1.E0/2.E0)*SL
     #N(1)*SN(1)*RN(1)**2+(1.E0/2.E0+BL(3)/2)*SLN(3)*SN(3)*RN(3)**2+(-BL
     #(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3+(1.E0/2.E0-BL(3)/2)*SLN(3)*SN(3)*
     #*3)
	CB(7,11)=((BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(1.E0/2.E0+B
     #L(2)/2)*SLN(2)*SN(2)*RN(2)**2+(-BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3
     #+(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3)
	CB(7,17)=
     #((1.E0/2.E0+BL(2)/2)*SLN(2)*SN(2)*RN(2)**2+(BL(3)/2+1.E0/2
     #.E0)*SLN(3)*SN(3)*RN(3)**2+(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3+(1.
     #E0/2.E0-BL(3)/2)*SLN(3)*SN(3)**3)

	CB(7,3) =
     #(RN(1)**2*BL(1)-SN(1)**2*BL(1)+SN(3)**2*BL(3)-RN(3)**2*BL(3))
	CB(7,9) =
     #(-RN(1)**2*BL(1)-SN(2)**2*BL(2)+SN(1)**2*BL(1)+RN(2)**2*BL(2))
	CB(7,15)=
     #(-RN(2)**2*BL(2)-SN(3)**2*BL(3)+SN(2)**2*BL(2)+RN(3)**2*BL(3))

C	EIGHTH ROW
C	------------------
	CB(8,4) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SLN
     #(3)*RN(3)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWS(1)-FLR(1)
	CB(8,10)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN
     #(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SL
     #N(2)*SN(2)**2*RN(2))
     #-DWS(2)-FLR(2)
	CB(8,16)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,2)/6-RS(1,3)/3)*SLN(2)*RN(
     #2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,2)/6-RS(1,3)/3)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*SLN
     #(3)*RN(3)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWS(3)-FLR(3)

	CB(8,5) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,2)/6+RS(1,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,3)/6+RS(1,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,2)/6+RS(1,1)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)+RS(1,3)/6+RS(1,1)
     #/3)*SLN(3)*SN(3)**3)
     #-DWS(4)
	CB(8,11)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,3)/6+RS(1,2)/3)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,1)/6+RS(1,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,3)/6+RS(1,2)
     #/3)*SLN(2)*SN(2)**3)
     #-DWS(5)
	CB(8,17)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,2)/6+RS(1,3)/3)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,1)/6+RS(1,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/
     #3)*SLN(2)*SN(2)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)+RS(1,1)/6+RS(1,3)
     #/3)*SLN(3)*SN(3)**3)
     #-DWS(6)

	CB(8,3) =
     #((RS(1,2)/2+RS(1,1)/2)*BL(1)*RN(1)**2+(-RS(1,1)/2-RS(1,3)/
     #2)*BL(3)*RN(3)**2+(-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)**2+(RS(1,1)/2
     #+RS(1,3)/2)*BL(3)*SN(3)**2)
     #-DWS(7)
	CB(8,9) =
     #((-RS(1,1)/2-RS(1,2)/2)*BL(1)*RN
     #(1)**2+(RS(1,2)/2+RS(1,3)/2)*BL(2)*RN(2)**2+(RS(1,2)/2+RS(1,1)/2)*
     #BL(1)*SN(1)**2+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)**2)
     #-DWS(8)
	CB(8,15)=
     #((-RS(
     #1,2)/2-RS(1,3)/2)*BL(2)*RN(2)**2+(RS(1,1)/2+RS(1,3)/2)*BL(3)*RN(3)
     #**2+(RS(1,2)/2+RS(1,3)/2)*BL(2)*SN(2)**2+(-RS(1,1)/2-RS(1,3)/2)*BL
     #(3)*SN(3)**2)
     #-DWS(9)

C	NINTH ROW
C	------------------
	CB(9,4) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*RN(
     #1)**3+((-RS(2,1)/4-RS(2,2)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SLN
     #(3)*RN(3)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWR(1)
	CB(9,10)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(2,1)/4-RS(2,2)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/3-RS(2,3)/6)*SLN
     #(2)*RN(2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,2)/3-RS(2,3)/6)*SL
     #N(2)*SN(2)**2*RN(2))
     #-DWR(2)
	CB(9,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*RN(
     #2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*SLN
     #(3)*RN(3)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWR(3)

	CB(9,5) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,2)/6+RS(2,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)+RS(2,3)/6+RS(2,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,1)/4-RS(2,2)/4)*BL(1)+RS(2,2)/6+RS(2,1)/
     #3)*SLN(1)*SN(1)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,3)/6+RS(2,1)
     #/3)*SLN(3)*SN(3)**3)
     #-DWR(4)-FLS(1)
	CB(9,11)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,1)/6+RS(2,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)/6)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(2,1)/4-RS(2,2)/4)*BL(1)+RS(2,1)/6+RS(2,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)
     #/6)*SLN(2)*SN(2)**3)
     #-DWR(5)-FLS(2)
	CB(9,17)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/6)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)+RS(2,1)/6+RS(2,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/
     #6)*SLN(2)*SN(2)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,1)/6+RS(2,3)
     #/3)*SLN(3)*SN(3)**3)
     #-DWR(6)-FLS(3)

	CB(9,3) =
     #((RS(2,2)/2+RS(2,1)/2)*BL(1)*RN(1)**2+(-RS(2,3)/2-RS(2,1)/
     #2)*BL(3)*RN(3)**2+(-RS(2,1)/2-RS(2,2)/2)*BL(1)*SN(1)**2+(RS(2,1)/2
     #+RS(2,3)/2)*BL(3)*SN(3)**2)
     #-DWR(7)
	CB(9,9) =
     #((-RS(2,1)/2-RS(2,2)/2)*BL(1)*RN
     #(1)**2+(RS(2,2)/2+RS(2,3)/2)*BL(2)*RN(2)**2+(RS(2,2)/2+RS(2,1)/2)*
     #BL(1)*SN(1)**2+(-RS(2,2)/2-RS(2,3)/2)*BL(2)*SN(2)**2)
     #-DWR(8)
	CB(9,15)=
     #((-RS(
     #2,2)/2-RS(2,3)/2)*BL(2)*RN(2)**2+(RS(2,1)/2+RS(2,3)/2)*BL(3)*RN(3)
     #**2+(RS(2,2)/2+RS(2,3)/2)*BL(2)*SN(2)**2+(-RS(2,3)/2-RS(2,1)/2)*BL
     #(3)*SN(3)**2)
     #-DWR(9)

	RETURN
	END
C
C=====================================================================


C
C=====================================================================
      SUBROUTINE PSIGRT2(IPT,SIGR,A,AR,AS,PAPM,PAPB,PAPSH)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTES FOR INTEGRAL(TRANS(P)*SIGR)dA
C
C	INPUT VARIABLES
C	SIGR(8)  = NODAL STRESS RESULTANTS
C	A   = FACTOR FOR INT(SIGR)dA
C     AR  = FACTOR FOR INT(R*SIGR)dA
C     AS  = FACTOR FOR INT(S*SIGR)dA
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C	PAPM(7)  = INTEGRAL(TRANS(Pm)*N)
C	PAPB(9) = INTEGRAL(TRANS(Pb)*M)
C	PAPSH(2) = INTEGRAL(TRANS(Ps)*Q)
C     ----------------------------------------------------------------
	DIMENSION SIGR(8)
	DIMENSION PAPM(7),PAPB(9),PAPSH(2)

	IF (IPT.EQ.1) THEN
C	  ----------------------------
C	  COMPUTE PAPM - MEMBRANE PART
C	  ----------------------------
	  PAPM(1) = SIGR(1)*A
	  PAPM(2) = SIGR(1)*AR-SIGR(3)*AS
	  PAPM(3) = SIGR(1)*AS
	  PAPM(4) = SIGR(2)*A
	  PAPM(5) = SIGR(2)*AR
	  PAPM(6) = SIGR(2)*AS-SIGR(3)*AR
	  PAPM(7) = SIGR(3)*A
C	  ---------------------------
C	  COMPUTE PAPB - BENDING PART
C	  ---------------------------
	  PAPB(1) = SIGR(4)*A
	  PAPB(2) = SIGR(4)*AR
	  PAPB(3) = SIGR(4)*AS
	  PAPB(4) = SIGR(5)*A
	  PAPB(5) = SIGR(5)*AR
	  PAPB(6) = SIGR(5)*AS
	  PAPB(7) = SIGR(6)*A
	  PAPB(8) = SIGR(6)*AR
	  PAPB(9) = SIGR(6)*AS
C	  -------------------------------------
C	  COMPUTE PAPSH - TRANSVERSE SHEAR PART
C	  -------------------------------------
	  PAPSH(1) = SIGR(7)*A
	  PAPSH(2) = SIGR(8)*A
	ELSE
C	  ----------------------------
C	  COMPUTE PAPM - MEMBRANE PART
C	  ----------------------------
	  PAPM(1) = PAPM(1)+SIGR(1)*A
	  PAPM(2) = PAPM(2)+SIGR(1)*AR-SIGR(3)*AS
	  PAPM(3) = PAPM(3)+SIGR(1)*AS
	  PAPM(4) = PAPM(4)+SIGR(2)*A
	  PAPM(5) = PAPM(5)+SIGR(2)*AR
	  PAPM(6) = PAPM(6)+SIGR(2)*AS-SIGR(3)*AR
	  PAPM(7) = PAPM(7)+SIGR(3)*A
C	  ---------------------------
C	  COMPUTE PAPB - BENDING PART
C	  ---------------------------
	  PAPB(1) = PAPB(1)+SIGR(4)*A
	  PAPB(2) = PAPB(2)+SIGR(4)*AR
	  PAPB(3) = PAPB(3)+SIGR(4)*AS
	  PAPB(4) = PAPB(4)+SIGR(5)*A
	  PAPB(5) = PAPB(5)+SIGR(5)*AR
	  PAPB(6) = PAPB(6)+SIGR(5)*AS
	  PAPB(7) = PAPB(7)+SIGR(6)*A
	  PAPB(8) = PAPB(8)+SIGR(6)*AR
	  PAPB(9) = PAPB(9)+SIGR(6)*AS
C	  -------------------------------------
C	  COMPUTE PAPSH - TRANSVERSE SHEAR PART
C	  -------------------------------------
	  PAPSH(1) = PAPSH(1)+SIGR(7)*A
	  PAPSH(2) = PAPSH(2)+SIGR(8)*A
	ENDIF

      RETURN
      END
C
C=====================================================================
      SUBROUTINE PDRPTRI(IPT,DR,A,AR,AS,ARS,ARR,ASS,PMPM,PMPB,PBPB,PTP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTES FOR INTEGRAL(TRANS(P)*DR*P)dA
C
C	INPUT VARIABLES
C	DR(64)  = NODAL STRESS RESULTANTS
C	A       = FACTOR FOR INT(DR)dA
C     AR      = FACTOR FOR INT(R*DR)dA
C     AS      = FACTOR FOR INT(S*DR)dA
C     ARS     = FACTOR FOR INT(R*S*DR)dA
C     ARR     = FACTOR FOR INT(R*R*DR)dA
C     ASS     = FACTOR FOR INT(S*S*DR)dA
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C	PMPM(7,7)   = INTEGRAL(TRANS(Pm)*DR*Pm) - MEMBRANE PART
C	PMPB(7,9)  = INTEGRAL(TRANS(Pm)*DR*Pb) - MEMBRANE-BENDING PART
C	PBPB(9,9) = INTEGRAL(TRANS(Pb)*DR*Pb) - BENDING PART
C	PTP(2,2)	= INTEGRAL(TRANS(Ps)*DR*Ps) - TRANSVERSE SHEAR PART
C     ----------------------------------------------------------------
	DIMENSION DR(64)
	DIMENSION PMPM(7,7),PMPB(7,9),PBPB(9,9), PTP(2,2)

	IF (IPT.EQ.1) THEN
C	  ----------------------------
C	  COMPUTE PMPM - MEMBRANE PART
C	  UPPER TRIANGULAR ONLY
C	  ----------------------------
        PMPM(1,1) = A*DR(1)
        PMPM(1,2) = AR*DR(1)-AS*DR(3)
        PMPM(1,3) = AS*DR(1)
        PMPM(1,4) = A*DR(2)
        PMPM(1,5) = AR*DR(2)
        PMPM(1,6) = AS*DR(2)-AR*DR(3)
        PMPM(1,7) = A*DR(3)

        PMPM(2,2) = ARR*DR(1)-2.0D0*ARS*DR(3)+ASS*DR(19)
        PMPM(2,3) = ARS*DR(1)-ASS*DR(3)
        PMPM(2,4) = AR*DR(2)-AS*DR(11)
        PMPM(2,5) = ARR*DR(2)-ARS*DR(11)
        PMPM(2,6) = ARS*DR(2)-ASS*DR(11)-ARR*DR(3)+ARS*DR(19)
        PMPM(2,7) = AR*DR(3)-AS*DR(19)

        PMPM(3,3) = ASS*DR(1)
        PMPM(3,4) = AS*DR(2)
        PMPM(3,5) = ARS*DR(2)
        PMPM(3,6) = ASS*DR(2)-ARS*DR(3)
        PMPM(3,7) = AS*DR(3)

        PMPM(4,4) = A*DR(10)
        PMPM(4,5) = AR*DR(10)
        PMPM(4,6) = AS*DR(10)-AR*DR(11)
        PMPM(4,7) = A*DR(11)

        PMPM(5,5) = ARR*DR(10)
        PMPM(5,6) = ARS*DR(10)-ARR*DR(11)
        PMPM(5,7) = AR*DR(11)

        PMPM(6,6) = ASS*DR(10)-2.0D0*ARS*DR(11)+ARR*DR(19)
        PMPM(6,7) = AS*DR(11)-AR*DR(19)

        PMPM(7,7) = A*DR(19)

C	  ---------------------------------------------
C	  COMPUTE PMPB - MEMBRANE-BENDING COUPLING PART
C	  FULL MATRIX
C	  ---------------------------------------------
        PMPB(1,1) = A*DR(4)
        PMPB(1,2) = AR*DR(4)
        PMPB(1,3) = AS*DR(4)
        PMPB(1,4) = A*DR(5)
        PMPB(1,5) = AR*DR(5)
        PMPB(1,6) = AS*DR(5)
        PMPB(1,7) = A*DR(6)
        PMPB(1,8) = AR*DR(6)
        PMPB(1,9) = AS*DR(6)

        PMPB(2,1) = AR*DR(4)-AS*DR(6)
        PMPB(2,2) = ARR*DR(4)-ARS*DR(6)
        PMPB(2,3) = ARS*DR(4)-ASS*DR(6)
        PMPB(2,4) = AR*DR(5)-AS*DR(14)
        PMPB(2,5) = ARR*DR(5)-ARS*DR(14)
        PMPB(2,6) = ARS*DR(5)-ASS*DR(14)
        PMPB(2,7) = AR*DR(6)-AS*DR(22)
        PMPB(2,8) = ARR*DR(6)-ARS*DR(22)
        PMPB(2,9) = ARS*DR(6)-ASS*DR(22)

        PMPB(3,1) = AS*DR(4)
        PMPB(3,2) = ARS*DR(4)
        PMPB(3,3) = ASS*DR(4)
        PMPB(3,4) = AS*DR(5)
        PMPB(3,5) = ARS*DR(5)
        PMPB(3,6) = ASS*DR(5)
        PMPB(3,7) = AS*DR(6)
        PMPB(3,8) = ARS*DR(6)
        PMPB(3,9) = ASS*DR(6)

        PMPB(4,1) = A*DR(5)
        PMPB(4,2) = AR*DR(5)
        PMPB(4,3) = AS*DR(5)
        PMPB(4,4) = A*DR(13)
        PMPB(4,5) = AR*DR(13)
        PMPB(4,6) = AS*DR(13)
        PMPB(4,7) = A*DR(14)
        PMPB(4,8) = AR*DR(14)
        PMPB(4,9) = AS*DR(14)

        PMPB(5,1) = AR*DR(5)
        PMPB(5,2) = ARR*DR(5)
        PMPB(5,3) = ARS*DR(5)
        PMPB(5,4) = AR*DR(13)
        PMPB(5,5) = ARR*DR(13)
        PMPB(5,6) = ARS*DR(13)
        PMPB(5,7) = AR*DR(14)
        PMPB(5,8) = ARR*DR(14)
        PMPB(5,9) = ARS*DR(14)

        PMPB(6,1) = AS*DR(5)-AR*DR(6)
        PMPB(6,2) = ARS*DR(5)-ARR*DR(6)
        PMPB(6,3) = ASS*DR(5)-ARS*DR(6)
        PMPB(6,4) = AS*DR(13)-AR*DR(14)
        PMPB(6,5) = ARS*DR(13)-ARR*DR(14)
        PMPB(6,6) = ASS*DR(13)-ARS*DR(14)
        PMPB(6,7) = AS*DR(14)-AR*DR(22)
        PMPB(6,8) = ARS*DR(14)-ARR*DR(22)
        PMPB(6,9) = ASS*DR(14)-ARS*DR(22)

        PMPB(7,1) = A*DR(6)
        PMPB(7,2) = AR*DR(6)
        PMPB(7,3) = AS*DR(6)
        PMPB(7,4) = A*DR(14)
        PMPB(7,5) = AR*DR(14)
        PMPB(7,6) = AS*DR(14)
        PMPB(7,7) = A*DR(22)
        PMPB(7,8) = AR*DR(22)
        PMPB(7,9) = AS*DR(22)

C	  ---------------------------
C	  COMPUTE PBPB - BENDING PART
C	  UPPER TRIANGULAR ONLY
C	  ---------------------------
        PBPB(1,1) = A*DR(28)
        PBPB(1,2) = AR*DR(28)
        PBPB(1,3) = AS*DR(28)
        PBPB(1,4) = A*DR(29)
        PBPB(1,5) = AR*DR(29)
        PBPB(1,6) = AS*DR(29)
        PBPB(1,7) = A*DR(30)
        PBPB(1,8) = AR*DR(30)
        PBPB(1,9) = AS*DR(30)

        PBPB(2,2) = ARR*DR(28)
        PBPB(2,3) = ARS*DR(28)
        PBPB(2,4) = AR*DR(29)
        PBPB(2,5) = ARR*DR(29)
        PBPB(2,6) = ARS*DR(29)
        PBPB(2,7) = AR*DR(30)
        PBPB(2,8) = ARR*DR(30)
        PBPB(2,9) = ARS*DR(30)

        PBPB(3,3) = ASS*DR(28)
        PBPB(3,4) = AS*DR(29)
        PBPB(3,5) = ARS*DR(29)
        PBPB(3,6) = ASS*DR(29)
        PBPB(3,7) = AS*DR(30)
        PBPB(3,8) = ARS*DR(30)
        PBPB(3,9) = ASS*DR(30)

        PBPB(4,4) = A*DR(37)
        PBPB(4,5) = AR*DR(37)
        PBPB(4,6) = AS*DR(37)
        PBPB(4,7) = A*DR(38)
        PBPB(4,8) = AR*DR(38)
        PBPB(4,9) = AS*DR(38)

        PBPB(5,5) = ARR*DR(37)
        PBPB(5,6) = ARS*DR(37)
        PBPB(5,7) = AR*DR(38)
        PBPB(5,8) = ARR*DR(38)
        PBPB(5,9) = ARS*DR(38)

        PBPB(6,6) = ASS*DR(37)
        PBPB(6,7) = AS*DR(38)
        PBPB(6,8) = ARS*DR(38)
        PBPB(6,9) = ASS*DR(38)

        PBPB(7,7) = A*DR(46)
        PBPB(7,8) = AR*DR(46)
        PBPB(7,9) = AS*DR(46)

        PBPB(8,8) = ARR*DR(46)
        PBPB(8,9) = ARS*DR(46)

        PBPB(9,9) = ASS*DR(46)

C	  -----------------------------------
C	  COMPUTE PTP - TRANSVERSE SHEAR PART
C	  UPPER TRIANGULAR ONLY
C	  -----------------------------------
        PTP(1,1) = A*DR(55)
        PTP(1,2) = A*DR(56)

        PTP(2,2) = A*DR(64)
	ELSE
C	  ----------------------------
C	  COMPUTE PMPM - MEMBRANE PART
C	  UPPER TRIANGULAR ONLY
C	  ----------------------------
        PMPM(1,1) = PMPM(1,1)+A*DR(1)
        PMPM(1,2) = PMPM(1,2)+AR*DR(1)-AS*DR(3)
        PMPM(1,3) = PMPM(1,3)+AS*DR(1)
        PMPM(1,4) = PMPM(1,4)+A*DR(2)
        PMPM(1,5) = PMPM(1,5)+AR*DR(2)
        PMPM(1,6) = PMPM(1,6)+AS*DR(2)-AR*DR(3)
        PMPM(1,7) = PMPM(1,7)+A*DR(3)

        PMPM(2,2) = PMPM(2,2)+ARR*DR(1)-2.0D0*ARS*DR(3)+ASS*DR(19)
        PMPM(2,3) = PMPM(2,3)+ARS*DR(1)-ASS*DR(3)
        PMPM(2,4) = PMPM(2,4)+AR*DR(2)-AS*DR(11)
        PMPM(2,5) = PMPM(2,5)+ARR*DR(2)-ARS*DR(11)
        PMPM(2,6) = PMPM(2,6)+ARS*DR(2)-ASS*DR(11)-ARR*DR(3)+ARS*DR(19)
        PMPM(2,7) = PMPM(2,7)+AR*DR(3)-AS*DR(19)

        PMPM(3,3) = PMPM(3,3)+ASS*DR(1)
        PMPM(3,4) = PMPM(3,4)+AS*DR(2)
        PMPM(3,5) = PMPM(3,5)+ARS*DR(2)
        PMPM(3,6) = PMPM(3,6)+ASS*DR(2)-ARS*DR(3)
        PMPM(3,7) = PMPM(3,7)+AS*DR(3)

        PMPM(4,4) = PMPM(4,4)+A*DR(10)
        PMPM(4,5) = PMPM(4,5)+AR*DR(10)
        PMPM(4,6) = PMPM(4,6)+AS*DR(10)-AR*DR(11)
        PMPM(4,7) = PMPM(4,7)+A*DR(11)

        PMPM(5,5) = PMPM(5,5)+ARR*DR(10)
        PMPM(5,6) = PMPM(5,6)+ARS*DR(10)-ARR*DR(11)
        PMPM(5,7) = PMPM(5,7)+AR*DR(11)

        PMPM(6,6) = PMPM(6,6)+ASS*DR(10)-2.0D0*ARS*DR(11)+ARR*DR(19)
        PMPM(6,7) = PMPM(6,7)+AS*DR(11)-AR*DR(19)

        PMPM(7,7) = PMPM(7,7)+A*DR(19)

C	  ---------------------------------------------
C	  COMPUTE PMPB - MEMBRANE-BENDING COUPLING PART
C	  FULL MATRIX
C	  ---------------------------------------------
        PMPB(1,1) = PMPB(1,1)+A*DR(4)
        PMPB(1,2) = PMPB(1,2)+AR*DR(4)
        PMPB(1,3) = PMPB(1,3)+AS*DR(4)
        PMPB(1,4) = PMPB(1,4)+A*DR(5)
        PMPB(1,5) = PMPB(1,5)+AR*DR(5)
        PMPB(1,6) = PMPB(1,6)+AS*DR(5)
        PMPB(1,7) = PMPB(1,7)+A*DR(6)
        PMPB(1,8) = PMPB(1,8)+AR*DR(6)
        PMPB(1,9) = PMPB(1,9)+AS*DR(6)

        PMPB(2,1) = PMPB(2,1)+AR*DR(4)-AS*DR(6)
        PMPB(2,2) = PMPB(2,2)+ARR*DR(4)-ARS*DR(6)
        PMPB(2,3) = PMPB(2,3)+ARS*DR(4)-ASS*DR(6)
        PMPB(2,4) = PMPB(2,4)+AR*DR(5)-AS*DR(14)
        PMPB(2,5) = PMPB(2,5)+ARR*DR(5)-ARS*DR(14)
        PMPB(2,6) = PMPB(2,6)+ARS*DR(5)-ASS*DR(14)
        PMPB(2,7) = PMPB(2,7)+AR*DR(6)-AS*DR(22)
        PMPB(2,8) = PMPB(2,8)+ARR*DR(6)-ARS*DR(22)
        PMPB(2,9) = PMPB(2,9)+ARS*DR(6)-ASS*DR(22)

        PMPB(3,1) = PMPB(3,1)+AS*DR(4)
        PMPB(3,2) = PMPB(3,2)+ARS*DR(4)
        PMPB(3,3) = PMPB(3,3)+ASS*DR(4)
        PMPB(3,4) = PMPB(3,4)+AS*DR(5)
        PMPB(3,5) = PMPB(3,5)+ARS*DR(5)
        PMPB(3,6) = PMPB(3,6)+ASS*DR(5)
        PMPB(3,7) = PMPB(3,7)+AS*DR(6)
        PMPB(3,8) = PMPB(3,8)+ARS*DR(6)
        PMPB(3,9) = PMPB(3,9)+ASS*DR(6)

        PMPB(4,1) = PMPB(4,1)+A*DR(5)
        PMPB(4,2) = PMPB(4,2)+AR*DR(5)
        PMPB(4,3) = PMPB(4,3)+AS*DR(5)
        PMPB(4,4) = PMPB(4,4)+A*DR(13)
        PMPB(4,5) = PMPB(4,5)+AR*DR(13)
        PMPB(4,6) = PMPB(4,6)+AS*DR(13)
        PMPB(4,7) = PMPB(4,7)+A*DR(14)
        PMPB(4,8) = PMPB(4,8)+AR*DR(14)
        PMPB(4,9) = PMPB(4,9)+AS*DR(14)

        PMPB(5,1) = PMPB(5,1)+AR*DR(5)
        PMPB(5,2) = PMPB(5,2)+ARR*DR(5)
        PMPB(5,3) = PMPB(5,3)+ARS*DR(5)
        PMPB(5,4) = PMPB(5,4)+AR*DR(13)
        PMPB(5,5) = PMPB(5,5)+ARR*DR(13)
        PMPB(5,6) = PMPB(5,6)+ARS*DR(13)
        PMPB(5,7) = PMPB(5,7)+AR*DR(14)
        PMPB(5,8) = PMPB(5,8)+ARR*DR(14)
        PMPB(5,9) = PMPB(5,9)+ARS*DR(14)

        PMPB(6,1) = PMPB(6,1)+AS*DR(5)-AR*DR(6)
        PMPB(6,2) = PMPB(6,2)+ARS*DR(5)-ARR*DR(6)
        PMPB(6,3) = PMPB(6,3)+ASS*DR(5)-ARS*DR(6)
        PMPB(6,4) = PMPB(6,4)+AS*DR(13)-AR*DR(14)
        PMPB(6,5) = PMPB(6,5)+ARS*DR(13)-ARR*DR(14)
        PMPB(6,6) = PMPB(6,6)+ASS*DR(13)-ARS*DR(14)
        PMPB(6,7) = PMPB(6,7)+AS*DR(14)-AR*DR(22)
        PMPB(6,8) = PMPB(6,8)+ARS*DR(14)-ARR*DR(22)
        PMPB(6,9) = PMPB(6,9)+ASS*DR(14)-ARS*DR(22)

        PMPB(7,1) = PMPB(7,1)+A*DR(6)
        PMPB(7,2) = PMPB(7,2)+AR*DR(6)
        PMPB(7,3) = PMPB(7,3)+AS*DR(6)
        PMPB(7,4) = PMPB(7,4)+A*DR(14)
        PMPB(7,5) = PMPB(7,5)+AR*DR(14)
        PMPB(7,6) = PMPB(7,6)+AS*DR(14)
        PMPB(7,7) = PMPB(7,7)+A*DR(22)
        PMPB(7,8) = PMPB(7,8)+AR*DR(22)
        PMPB(7,9) = PMPB(7,9)+AS*DR(22)

C	  ---------------------------
C	  COMPUTE PBPB - BENDING PART
C	  UPPER TRIANGULAR ONLY
C	  ---------------------------
        PBPB(1,1) = PBPB(1,1)+A*DR(28)
        PBPB(1,2) = PBPB(1,2)+AR*DR(28)
        PBPB(1,3) = PBPB(1,3)+AS*DR(28)
        PBPB(1,4) = PBPB(1,4)+A*DR(29)
        PBPB(1,5) = PBPB(1,5)+AR*DR(29)
        PBPB(1,6) = PBPB(1,6)+AS*DR(29)
        PBPB(1,7) = PBPB(1,7)+A*DR(30)
        PBPB(1,8) = PBPB(1,8)+AR*DR(30)
        PBPB(1,9) = PBPB(1,9)+AS*DR(30)

        PBPB(2,2) = PBPB(2,2)+ARR*DR(28)
        PBPB(2,3) = PBPB(2,3)+ARS*DR(28)
        PBPB(2,4) = PBPB(2,4)+AR*DR(29)
        PBPB(2,5) = PBPB(2,5)+ARR*DR(29)
        PBPB(2,6) = PBPB(2,6)+ARS*DR(29)
        PBPB(2,7) = PBPB(2,7)+AR*DR(30)
        PBPB(2,8) = PBPB(2,8)+ARR*DR(30)
        PBPB(2,9) = PBPB(2,9)+ARS*DR(30)

        PBPB(3,3) = PBPB(3,3)+ASS*DR(28)
        PBPB(3,4) = PBPB(3,4)+AS*DR(29)
        PBPB(3,5) = PBPB(3,5)+ARS*DR(29)
        PBPB(3,6) = PBPB(3,6)+ASS*DR(29)
        PBPB(3,7) = PBPB(3,7)+AS*DR(30)
        PBPB(3,8) = PBPB(3,8)+ARS*DR(30)
        PBPB(3,9) = PBPB(3,9)+ASS*DR(30)

        PBPB(4,4) = PBPB(4,4)+A*DR(37)
        PBPB(4,5) = PBPB(4,5)+AR*DR(37)
        PBPB(4,6) = PBPB(4,6)+AS*DR(37)
        PBPB(4,7) = PBPB(4,7)+A*DR(38)
        PBPB(4,8) = PBPB(4,8)+AR*DR(38)
        PBPB(4,9) = PBPB(4,9)+AS*DR(38)

        PBPB(5,5) = PBPB(5,5)+ARR*DR(37)
        PBPB(5,6) = PBPB(5,6)+ARS*DR(37)
        PBPB(5,7) = PBPB(5,7)+AR*DR(38)
        PBPB(5,8) = PBPB(5,8)+ARR*DR(38)
        PBPB(5,9) = PBPB(5,9)+ARS*DR(38)

        PBPB(6,6) = PBPB(6,6)+ASS*DR(37)
        PBPB(6,7) = PBPB(6,7)+AS*DR(38)
        PBPB(6,8) = PBPB(6,8)+ARS*DR(38)
        PBPB(6,9) = PBPB(6,9)+ASS*DR(38)

        PBPB(7,7) = PBPB(7,7)+A*DR(46)
        PBPB(7,8) = PBPB(7,8)+AR*DR(46)
        PBPB(7,9) = PBPB(7,9)+AS*DR(46)

        PBPB(8,8) = PBPB(8,8)+ARR*DR(46)
        PBPB(8,9) = PBPB(8,9)+ARS*DR(46)

        PBPB(9,9) = PBPB(9,9)+ASS*DR(46)

C	  -----------------------------------
C	  COMPUTE PTP - TRANSVERSE SHEAR PART
C	  UPPER TRIANGULAR ONLY
C	  -----------------------------------
        PTP(1,1) = PTP(1,1)+A*DR(55)
        PTP(1,2) = PTP(1,2)+A*DR(56)

        PTP(2,2) = PTP(2,2)+A*DR(64)
	ENDIF

      RETURN
      END
C
C=====================================================================
      SUBROUTINE KMEBE_TRI(ACM,ACB,PMPB,ST)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-----------------------------------------------------------
C	PURPOSE:	TO ADD THE FLEXURAL CONTRIBUTION TO THE
C					ELEMENT STIFFNESS MATRIX
C
C	INPUT VARIABLES
C	ACM(7,18)	= INVERSE(Am)*(Cm)
C	ACB(9,18)	= INVERSE(Ab)*(Cb)
C	PMPB(7,9)	= Nm*B*Pb
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIAVBLES
C	ST(171)		= STIFFNESS MATRIX IN ROW FORMAT
C	-----------------------------------------------------------
      DIMENSION ACM(7,18),ACB(9,18),PMPB(7,9)
	DIMENSION ST(171)

      t8 = ACM(1,1)*PMPB(1,1)+ACM(2,1)*PMPB(2,1)+ACM(3,1)*PMPB(3,1)+ACM(
     #4,1)*PMPB(4,1)+ACM(5,1)*PMPB(5,1)+ACM(6,1)*PMPB(6,1)+ACM(7,1)*PMPB
     #(7,1)
      t17 = ACM(1,1)*PMPB(1,2)+ACM(2,1)*PMPB(2,2)+ACM(3,1)*PMPB(3,2)+ACM
     #(4,1)*PMPB(4,2)+ACM(5,1)*PMPB(5,2)+ACM(6,1)*PMPB(6,2)+ACM(7,1)*PMP
     #B(7,2)
      t26 = ACM(1,1)*PMPB(1,3)+ACM(2,1)*PMPB(2,3)+ACM(3,1)*PMPB(3,3)+ACM
     #(4,1)*PMPB(4,3)+ACM(5,1)*PMPB(5,3)+ACM(6,1)*PMPB(6,3)+ACM(7,1)*PMP
     #B(7,3)
      t35 = ACM(1,1)*PMPB(1,4)+ACM(2,1)*PMPB(2,4)+ACM(3,1)*PMPB(3,4)+ACM
     #(4,1)*PMPB(4,4)+ACM(5,1)*PMPB(5,4)+ACM(6,1)*PMPB(6,4)+ACM(7,1)*PMP
     #B(7,4)
      t44 = ACM(1,1)*PMPB(1,5)+ACM(2,1)*PMPB(2,5)+ACM(3,1)*PMPB(3,5)+ACM
     #(4,1)*PMPB(4,5)+ACM(5,1)*PMPB(5,5)+ACM(6,1)*PMPB(6,5)+ACM(7,1)*PMP
     #B(7,5)
      t53 = ACM(1,1)*PMPB(1,6)+ACM(2,1)*PMPB(2,6)+ACM(3,1)*PMPB(3,6)+ACM
     #(4,1)*PMPB(4,6)+ACM(5,1)*PMPB(5,6)+ACM(6,1)*PMPB(6,6)+ACM(7,1)*PMP
     #B(7,6)
      t62 = ACM(1,1)*PMPB(1,7)+ACM(2,1)*PMPB(2,7)+ACM(3,1)*PMPB(3,7)+ACM
     #(4,1)*PMPB(4,7)+ACM(5,1)*PMPB(5,7)+ACM(6,1)*PMPB(6,7)+ACM(7,1)*PMP
     #B(7,7)
      t71 = ACM(1,1)*PMPB(1,8)+ACM(2,1)*PMPB(2,8)+ACM(3,1)*PMPB(3,8)+ACM
     #(4,1)*PMPB(4,8)+ACM(5,1)*PMPB(5,8)+ACM(6,1)*PMPB(6,8)+ACM(7,1)*PMP
     #B(7,8)
      t80 = ACM(1,1)*PMPB(1,9)+ACM(2,1)*PMPB(2,9)+ACM(3,1)*PMPB(3,9)+ACM
     #(4,1)*PMPB(4,9)+ACM(5,1)*PMPB(5,9)+ACM(6,1)*PMPB(6,9)+ACM(7,1)*PMP
     #B(7,9)
      t170 = ACM(1,2)*PMPB(1,1)+ACM(2,2)*PMPB(2,1)+ACM(3,2)*PMPB(3,1)+AC
     #M(4,2)*PMPB(4,1)+ACM(5,2)*PMPB(5,1)+ACM(6,2)*PMPB(6,1)+ACM(7,2)*PM
     #PB(7,1)
      t179 = ACM(1,2)*PMPB(1,2)+ACM(2,2)*PMPB(2,2)+ACM(3,2)*PMPB(3,2)+AC
     #M(4,2)*PMPB(4,2)+ACM(5,2)*PMPB(5,2)+ACM(6,2)*PMPB(6,2)+ACM(7,2)*PM
     #PB(7,2)
      t188 = ACM(1,2)*PMPB(1,3)+ACM(2,2)*PMPB(2,3)+ACM(3,2)*PMPB(3,3)+AC
     #M(4,2)*PMPB(4,3)+ACM(5,2)*PMPB(5,3)+ACM(6,2)*PMPB(6,3)+ACM(7,2)*PM
     #PB(7,3)
      t197 = ACM(1,2)*PMPB(1,4)+ACM(2,2)*PMPB(2,4)+ACM(3,2)*PMPB(3,4)+AC
     #M(4,2)*PMPB(4,4)+ACM(5,2)*PMPB(5,4)+ACM(6,2)*PMPB(6,4)+ACM(7,2)*PM
     #PB(7,4)
      t206 = ACM(1,2)*PMPB(1,5)+ACM(2,2)*PMPB(2,5)+ACM(3,2)*PMPB(3,5)+AC
     #M(4,2)*PMPB(4,5)+ACM(5,2)*PMPB(5,5)+ACM(6,2)*PMPB(6,5)+ACM(7,2)*PM
     #PB(7,5)
      t215 = ACM(1,2)*PMPB(1,6)+ACM(2,2)*PMPB(2,6)+ACM(3,2)*PMPB(3,6)+AC
     #M(4,2)*PMPB(4,6)+ACM(5,2)*PMPB(5,6)+ACM(6,2)*PMPB(6,6)+ACM(7,2)*PM
     #PB(7,6)
      t224 = ACM(1,2)*PMPB(1,7)+ACM(2,2)*PMPB(2,7)+ACM(3,2)*PMPB(3,7)+AC
     #M(4,2)*PMPB(4,7)+ACM(5,2)*PMPB(5,7)+ACM(6,2)*PMPB(6,7)+ACM(7,2)*PM
     #PB(7,7)
      t233 = ACM(1,2)*PMPB(1,8)+ACM(2,2)*PMPB(2,8)+ACM(3,2)*PMPB(3,8)+AC
     #M(4,2)*PMPB(4,8)+ACM(5,2)*PMPB(5,8)+ACM(6,2)*PMPB(6,8)+ACM(7,2)*PM
     #PB(7,8)
      t242 = ACM(1,2)*PMPB(1,9)+ACM(2,2)*PMPB(2,9)+ACM(3,2)*PMPB(3,9)+AC
     #M(4,2)*PMPB(4,9)+ACM(5,2)*PMPB(5,9)+ACM(6,2)*PMPB(6,9)+ACM(7,2)*PM
     #PB(7,9)
      t332 = ACM(1,6)*PMPB(1,1)+ACM(2,6)*PMPB(2,1)+ACM(3,6)*PMPB(3,1)+AC
     #M(4,6)*PMPB(4,1)+ACM(5,6)*PMPB(5,1)+ACM(6,6)*PMPB(6,1)+ACM(7,6)*PM
     #PB(7,1)
      t341 = ACM(1,6)*PMPB(1,2)+ACM(2,6)*PMPB(2,2)+ACM(3,6)*PMPB(3,2)+AC
     #M(4,6)*PMPB(4,2)+ACM(5,6)*PMPB(5,2)+ACM(6,6)*PMPB(6,2)+ACM(7,6)*PM
     #PB(7,2)
      t350 = ACM(1,6)*PMPB(1,3)+ACM(2,6)*PMPB(2,3)+ACM(3,6)*PMPB(3,3)+AC
     #M(4,6)*PMPB(4,3)+ACM(5,6)*PMPB(5,3)+ACM(6,6)*PMPB(6,3)+ACM(7,6)*PM
     #PB(7,3)
      t359 = ACM(1,6)*PMPB(1,4)+ACM(2,6)*PMPB(2,4)+ACM(3,6)*PMPB(3,4)+AC
     #M(4,6)*PMPB(4,4)+ACM(5,6)*PMPB(5,4)+ACM(6,6)*PMPB(6,4)+ACM(7,6)*PM
     #PB(7,4)
      t368 = ACM(1,6)*PMPB(1,5)+ACM(2,6)*PMPB(2,5)+ACM(3,6)*PMPB(3,5)+AC
     #M(4,6)*PMPB(4,5)+ACM(5,6)*PMPB(5,5)+ACM(6,6)*PMPB(6,5)+ACM(7,6)*PM
     #PB(7,5)
      t377 = ACM(1,6)*PMPB(1,6)+ACM(2,6)*PMPB(2,6)+ACM(3,6)*PMPB(3,6)+AC
     #M(4,6)*PMPB(4,6)+ACM(5,6)*PMPB(5,6)+ACM(6,6)*PMPB(6,6)+ACM(7,6)*PM
     #PB(7,6)
      t386 = ACM(1,6)*PMPB(1,7)+ACM(2,6)*PMPB(2,7)+ACM(3,6)*PMPB(3,7)+AC
     #M(4,6)*PMPB(4,7)+ACM(5,6)*PMPB(5,7)+ACM(6,6)*PMPB(6,7)+ACM(7,6)*PM
     #PB(7,7)
      t395 = ACM(1,6)*PMPB(1,8)+ACM(2,6)*PMPB(2,8)+ACM(3,6)*PMPB(3,8)+AC
     #M(4,6)*PMPB(4,8)+ACM(5,6)*PMPB(5,8)+ACM(6,6)*PMPB(6,8)+ACM(7,6)*PM
     #PB(7,8)
      t404 = ACM(1,6)*PMPB(1,9)+ACM(2,6)*PMPB(2,9)+ACM(3,6)*PMPB(3,9)+AC
     #M(4,6)*PMPB(4,9)+ACM(5,6)*PMPB(5,9)+ACM(6,6)*PMPB(6,9)+ACM(7,6)*PM
     #PB(7,9)
      t414 = ACM(1,7)*PMPB(1,1)+ACM(2,7)*PMPB(2,1)+ACM(3,7)*PMPB(3,1)+AC
     #M(4,7)*PMPB(4,1)+ACM(5,7)*PMPB(5,1)+ACM(6,7)*PMPB(6,1)+ACM(7,7)*PM
     #PB(7,1)
      t423 = ACM(1,7)*PMPB(1,2)+ACM(2,7)*PMPB(2,2)+ACM(3,7)*PMPB(3,2)+AC
     #M(4,7)*PMPB(4,2)+ACM(5,7)*PMPB(5,2)+ACM(6,7)*PMPB(6,2)+ACM(7,7)*PM
     #PB(7,2)
      t432 = ACM(1,7)*PMPB(1,3)+ACM(2,7)*PMPB(2,3)+ACM(3,7)*PMPB(3,3)+AC
     #M(4,7)*PMPB(4,3)+ACM(5,7)*PMPB(5,3)+ACM(6,7)*PMPB(6,3)+ACM(7,7)*PM
     #PB(7,3)
      t441 = ACM(1,7)*PMPB(1,4)+ACM(2,7)*PMPB(2,4)+ACM(3,7)*PMPB(3,4)+AC
     #M(4,7)*PMPB(4,4)+ACM(5,7)*PMPB(5,4)+ACM(6,7)*PMPB(6,4)+ACM(7,7)*PM
     #PB(7,4)
      t450 = ACM(1,7)*PMPB(1,5)+ACM(2,7)*PMPB(2,5)+ACM(3,7)*PMPB(3,5)+AC
     #M(4,7)*PMPB(4,5)+ACM(5,7)*PMPB(5,5)+ACM(6,7)*PMPB(6,5)+ACM(7,7)*PM
     #PB(7,5)
      t459 = ACM(1,7)*PMPB(1,6)+ACM(2,7)*PMPB(2,6)+ACM(3,7)*PMPB(3,6)+AC
     #M(4,7)*PMPB(4,6)+ACM(5,7)*PMPB(5,6)+ACM(6,7)*PMPB(6,6)+ACM(7,7)*PM
     #PB(7,6)
      t468 = ACM(1,7)*PMPB(1,7)+ACM(2,7)*PMPB(2,7)+ACM(3,7)*PMPB(3,7)+AC
     #M(4,7)*PMPB(4,7)+ACM(5,7)*PMPB(5,7)+ACM(6,7)*PMPB(6,7)+ACM(7,7)*PM
     #PB(7,7)
      t477 = ACM(1,7)*PMPB(1,8)+ACM(2,7)*PMPB(2,8)+ACM(3,7)*PMPB(3,8)+AC
     #M(4,7)*PMPB(4,8)+ACM(5,7)*PMPB(5,8)+ACM(6,7)*PMPB(6,8)+ACM(7,7)*PM
     #PB(7,8)
      t486 = ACM(1,7)*PMPB(1,9)+ACM(2,7)*PMPB(2,9)+ACM(3,7)*PMPB(3,9)+AC
     #M(4,7)*PMPB(4,9)+ACM(5,7)*PMPB(5,9)+ACM(6,7)*PMPB(6,9)+ACM(7,7)*PM
     #PB(7,9)
      t496 = ACM(1,8)*PMPB(1,1)+ACM(2,8)*PMPB(2,1)+ACM(3,8)*PMPB(3,1)+AC
     #M(4,8)*PMPB(4,1)+ACM(5,8)*PMPB(5,1)+ACM(6,8)*PMPB(6,1)+ACM(7,8)*PM
     #PB(7,1)
      t505 = ACM(1,8)*PMPB(1,2)+ACM(2,8)*PMPB(2,2)+ACM(3,8)*PMPB(3,2)+AC
     #M(4,8)*PMPB(4,2)+ACM(5,8)*PMPB(5,2)+ACM(6,8)*PMPB(6,2)+ACM(7,8)*PM
     #PB(7,2)
      t514 = ACM(1,8)*PMPB(1,3)+ACM(2,8)*PMPB(2,3)+ACM(3,8)*PMPB(3,3)+AC
     #M(4,8)*PMPB(4,3)+ACM(5,8)*PMPB(5,3)+ACM(6,8)*PMPB(6,3)+ACM(7,8)*PM
     #PB(7,3)
      t523 = ACM(1,8)*PMPB(1,4)+ACM(2,8)*PMPB(2,4)+ACM(3,8)*PMPB(3,4)+AC
     #M(4,8)*PMPB(4,4)+ACM(5,8)*PMPB(5,4)+ACM(6,8)*PMPB(6,4)+ACM(7,8)*PM
     #PB(7,4)
      t532 = ACM(1,8)*PMPB(1,5)+ACM(2,8)*PMPB(2,5)+ACM(3,8)*PMPB(3,5)+AC
     #M(4,8)*PMPB(4,5)+ACM(5,8)*PMPB(5,5)+ACM(6,8)*PMPB(6,5)+ACM(7,8)*PM
     #PB(7,5)
      t541 = ACM(1,8)*PMPB(1,6)+ACM(2,8)*PMPB(2,6)+ACM(3,8)*PMPB(3,6)+AC
     #M(4,8)*PMPB(4,6)+ACM(5,8)*PMPB(5,6)+ACM(6,8)*PMPB(6,6)+ACM(7,8)*PM
     #PB(7,6)
      t550 = ACM(1,8)*PMPB(1,7)+ACM(2,8)*PMPB(2,7)+ACM(3,8)*PMPB(3,7)+AC
     #M(4,8)*PMPB(4,7)+ACM(5,8)*PMPB(5,7)+ACM(6,8)*PMPB(6,7)+ACM(7,8)*PM
     #PB(7,7)
      t559 = ACM(1,8)*PMPB(1,8)+ACM(2,8)*PMPB(2,8)+ACM(3,8)*PMPB(3,8)+AC
     #M(4,8)*PMPB(4,8)+ACM(5,8)*PMPB(5,8)+ACM(6,8)*PMPB(6,8)+ACM(7,8)*PM
     #PB(7,8)
      t568 = ACM(1,8)*PMPB(1,9)+ACM(2,8)*PMPB(2,9)+ACM(3,8)*PMPB(3,9)+AC
     #M(4,8)*PMPB(4,9)+ACM(5,8)*PMPB(5,9)+ACM(6,8)*PMPB(6,9)+ACM(7,8)*PM
     #PB(7,9)
      t578 = ACM(1,12)*PMPB(1,1)+ACM(2,12)*PMPB(2,1)+ACM(3,12)*PMPB(3,1)
     #+ACM(4,12)*PMPB(4,1)+ACM(5,12)*PMPB(5,1)+ACM(6,12)*PMPB(6,1)+ACM(7
     #,12)*PMPB(7,1)
      t587 = ACM(1,12)*PMPB(1,2)+ACM(2,12)*PMPB(2,2)+ACM(3,12)*PMPB(3,2)
     #+ACM(4,12)*PMPB(4,2)+ACM(5,12)*PMPB(5,2)+ACM(6,12)*PMPB(6,2)+ACM(7
     #,12)*PMPB(7,2)
      t596 = ACM(1,12)*PMPB(1,3)+ACM(2,12)*PMPB(2,3)+ACM(3,12)*PMPB(3,3)
     #+ACM(4,12)*PMPB(4,3)+ACM(5,12)*PMPB(5,3)+ACM(6,12)*PMPB(6,3)+ACM(7
     #,12)*PMPB(7,3)
      t605 = ACM(1,12)*PMPB(1,4)+ACM(2,12)*PMPB(2,4)+ACM(3,12)*PMPB(3,4)
     #+ACM(4,12)*PMPB(4,4)+ACM(5,12)*PMPB(5,4)+ACM(6,12)*PMPB(6,4)+ACM(7
     #,12)*PMPB(7,4)
      t614 = ACM(1,12)*PMPB(1,5)+ACM(2,12)*PMPB(2,5)+ACM(3,12)*PMPB(3,5)
     #+ACM(4,12)*PMPB(4,5)+ACM(5,12)*PMPB(5,5)+ACM(6,12)*PMPB(6,5)+ACM(7
     #,12)*PMPB(7,5)
      t623 = ACM(1,12)*PMPB(1,6)+ACM(2,12)*PMPB(2,6)+ACM(3,12)*PMPB(3,6)
     #+ACM(4,12)*PMPB(4,6)+ACM(5,12)*PMPB(5,6)+ACM(6,12)*PMPB(6,6)+ACM(7
     #,12)*PMPB(7,6)
      t632 = ACM(1,12)*PMPB(1,7)+ACM(2,12)*PMPB(2,7)+ACM(3,12)*PMPB(3,7)
     #+ACM(4,12)*PMPB(4,7)+ACM(5,12)*PMPB(5,7)+ACM(6,12)*PMPB(6,7)+ACM(7
     #,12)*PMPB(7,7)
      t641 = ACM(1,12)*PMPB(1,8)+ACM(2,12)*PMPB(2,8)+ACM(3,12)*PMPB(3,8)
     #+ACM(4,12)*PMPB(4,8)+ACM(5,12)*PMPB(5,8)+ACM(6,12)*PMPB(6,8)+ACM(7
     #,12)*PMPB(7,8)
      t650 = ACM(1,12)*PMPB(1,9)+ACM(2,12)*PMPB(2,9)+ACM(3,12)*PMPB(3,9)
     #+ACM(4,12)*PMPB(4,9)+ACM(5,12)*PMPB(5,9)+ACM(6,12)*PMPB(6,9)+ACM(7
     #,12)*PMPB(7,9)
      t660 = ACM(1,13)*PMPB(1,1)+ACM(2,13)*PMPB(2,1)+ACM(3,13)*PMPB(3,1)
     #+ACM(4,13)*PMPB(4,1)+ACM(5,13)*PMPB(5,1)+ACM(6,13)*PMPB(6,1)+ACM(7
     #,13)*PMPB(7,1)
      t669 = ACM(1,13)*PMPB(1,2)+ACM(2,13)*PMPB(2,2)+ACM(3,13)*PMPB(3,2)
     #+ACM(4,13)*PMPB(4,2)+ACM(5,13)*PMPB(5,2)+ACM(6,13)*PMPB(6,2)+ACM(7
     #,13)*PMPB(7,2)
      t678 = ACM(1,13)*PMPB(1,3)+ACM(2,13)*PMPB(2,3)+ACM(3,13)*PMPB(3,3)
     #+ACM(4,13)*PMPB(4,3)+ACM(5,13)*PMPB(5,3)+ACM(6,13)*PMPB(6,3)+ACM(7
     #,13)*PMPB(7,3)
      t687 = ACM(1,13)*PMPB(1,4)+ACM(2,13)*PMPB(2,4)+ACM(3,13)*PMPB(3,4)
     #+ACM(4,13)*PMPB(4,4)+ACM(5,13)*PMPB(5,4)+ACM(6,13)*PMPB(6,4)+ACM(7
     #,13)*PMPB(7,4)
      t696 = ACM(1,13)*PMPB(1,5)+ACM(2,13)*PMPB(2,5)+ACM(3,13)*PMPB(3,5)
     #+ACM(4,13)*PMPB(4,5)+ACM(5,13)*PMPB(5,5)+ACM(6,13)*PMPB(6,5)+ACM(7
     #,13)*PMPB(7,5)
      t705 = ACM(1,13)*PMPB(1,6)+ACM(2,13)*PMPB(2,6)+ACM(3,13)*PMPB(3,6)
     #+ACM(4,13)*PMPB(4,6)+ACM(5,13)*PMPB(5,6)+ACM(6,13)*PMPB(6,6)+ACM(7
     #,13)*PMPB(7,6)
      t714 = ACM(1,13)*PMPB(1,7)+ACM(2,13)*PMPB(2,7)+ACM(3,13)*PMPB(3,7)
     #+ACM(4,13)*PMPB(4,7)+ACM(5,13)*PMPB(5,7)+ACM(6,13)*PMPB(6,7)+ACM(7
     #,13)*PMPB(7,7)
      t723 = ACM(1,13)*PMPB(1,8)+ACM(2,13)*PMPB(2,8)+ACM(3,13)*PMPB(3,8)
     #+ACM(4,13)*PMPB(4,8)+ACM(5,13)*PMPB(5,8)+ACM(6,13)*PMPB(6,8)+ACM(7
     #,13)*PMPB(7,8)
      t732 = ACM(1,13)*PMPB(1,9)+ACM(2,13)*PMPB(2,9)+ACM(3,13)*PMPB(3,9)
     #+ACM(4,13)*PMPB(4,9)+ACM(5,13)*PMPB(5,9)+ACM(6,13)*PMPB(6,9)+ACM(7
     #,13)*PMPB(7,9)
      t742 = ACM(1,14)*PMPB(1,1)+ACM(2,14)*PMPB(2,1)+ACM(3,14)*PMPB(3,1)
     #+ACM(4,14)*PMPB(4,1)+ACM(5,14)*PMPB(5,1)+ACM(6,14)*PMPB(6,1)+ACM(7
     #,14)*PMPB(7,1)
      t751 = ACM(1,14)*PMPB(1,2)+ACM(2,14)*PMPB(2,2)+ACM(3,14)*PMPB(3,2)
     #+ACM(4,14)*PMPB(4,2)+ACM(5,14)*PMPB(5,2)+ACM(6,14)*PMPB(6,2)+ACM(7
     #,14)*PMPB(7,2)
      t760 = ACM(1,14)*PMPB(1,3)+ACM(2,14)*PMPB(2,3)+ACM(3,14)*PMPB(3,3)
     #+ACM(4,14)*PMPB(4,3)+ACM(5,14)*PMPB(5,3)+ACM(6,14)*PMPB(6,3)+ACM(7
     #,14)*PMPB(7,3)
      t769 = ACM(1,14)*PMPB(1,4)+ACM(2,14)*PMPB(2,4)+ACM(3,14)*PMPB(3,4)
     #+ACM(4,14)*PMPB(4,4)+ACM(5,14)*PMPB(5,4)+ACM(6,14)*PMPB(6,4)+ACM(7
     #,14)*PMPB(7,4)
      t778 = ACM(1,14)*PMPB(1,5)+ACM(2,14)*PMPB(2,5)+ACM(3,14)*PMPB(3,5)
     #+ACM(4,14)*PMPB(4,5)+ACM(5,14)*PMPB(5,5)+ACM(6,14)*PMPB(6,5)+ACM(7
     #,14)*PMPB(7,5)
      t787 = ACM(1,14)*PMPB(1,6)+ACM(2,14)*PMPB(2,6)+ACM(3,14)*PMPB(3,6)
     #+ACM(4,14)*PMPB(4,6)+ACM(5,14)*PMPB(5,6)+ACM(6,14)*PMPB(6,6)+ACM(7
     #,14)*PMPB(7,6)
      t796 = ACM(1,14)*PMPB(1,7)+ACM(2,14)*PMPB(2,7)+ACM(3,14)*PMPB(3,7)
     #+ACM(4,14)*PMPB(4,7)+ACM(5,14)*PMPB(5,7)+ACM(6,14)*PMPB(6,7)+ACM(7
     #,14)*PMPB(7,7)
      t805 = ACM(1,14)*PMPB(1,8)+ACM(2,14)*PMPB(2,8)+ACM(3,14)*PMPB(3,8)
     #+ACM(4,14)*PMPB(4,8)+ACM(5,14)*PMPB(5,8)+ACM(6,14)*PMPB(6,8)+ACM(7
     #,14)*PMPB(7,8)
      t814 = ACM(1,14)*PMPB(1,9)+ACM(2,14)*PMPB(2,9)+ACM(3,14)*PMPB(3,9)
     #+ACM(4,14)*PMPB(4,9)+ACM(5,14)*PMPB(5,9)+ACM(6,14)*PMPB(6,9)+ACM(7
     #,14)*PMPB(7,9)
      t824 = ACM(1,18)*PMPB(1,1)+ACM(2,18)*PMPB(2,1)+ACM(3,18)*PMPB(3,1)
     #+ACM(4,18)*PMPB(4,1)+ACM(5,18)*PMPB(5,1)+ACM(6,18)*PMPB(6,1)+ACM(7
     #,18)*PMPB(7,1)
      t833 = ACM(1,18)*PMPB(1,2)+ACM(2,18)*PMPB(2,2)+ACM(3,18)*PMPB(3,2)
     #+ACM(4,18)*PMPB(4,2)+ACM(5,18)*PMPB(5,2)+ACM(6,18)*PMPB(6,2)+ACM(7
     #,18)*PMPB(7,2)
      t842 = ACM(1,18)*PMPB(1,3)+ACM(2,18)*PMPB(2,3)+ACM(3,18)*PMPB(3,3)
     #+ACM(4,18)*PMPB(4,3)+ACM(5,18)*PMPB(5,3)+ACM(6,18)*PMPB(6,3)+ACM(7
     #,18)*PMPB(7,3)
      t851 = ACM(1,18)*PMPB(1,4)+ACM(2,18)*PMPB(2,4)+ACM(3,18)*PMPB(3,4)
     #+ACM(4,18)*PMPB(4,4)+ACM(5,18)*PMPB(5,4)+ACM(6,18)*PMPB(6,4)+ACM(7
     #,18)*PMPB(7,4)
      t860 = ACM(1,18)*PMPB(1,5)+ACM(2,18)*PMPB(2,5)+ACM(3,18)*PMPB(3,5)
     #+ACM(4,18)*PMPB(4,5)+ACM(5,18)*PMPB(5,5)+ACM(6,18)*PMPB(6,5)+ACM(7
     #,18)*PMPB(7,5)
      t869 = ACM(1,18)*PMPB(1,6)+ACM(2,18)*PMPB(2,6)+ACM(3,18)*PMPB(3,6)
     #+ACM(4,18)*PMPB(4,6)+ACM(5,18)*PMPB(5,6)+ACM(6,18)*PMPB(6,6)+ACM(7
     #,18)*PMPB(7,6)
      t878 = ACM(1,18)*PMPB(1,7)+ACM(2,18)*PMPB(2,7)+ACM(3,18)*PMPB(3,7)
     #+ACM(4,18)*PMPB(4,7)+ACM(5,18)*PMPB(5,7)+ACM(6,18)*PMPB(6,7)+ACM(7
     #,18)*PMPB(7,7)
      t887 = ACM(1,18)*PMPB(1,8)+ACM(2,18)*PMPB(2,8)+ACM(3,18)*PMPB(3,8)
     #+ACM(4,18)*PMPB(4,8)+ACM(5,18)*PMPB(5,8)+ACM(6,18)*PMPB(6,8)+ACM(7
     #,18)*PMPB(7,8)
      t896 = ACM(1,18)*PMPB(1,9)+ACM(2,18)*PMPB(2,9)+ACM(3,18)*PMPB(3,9)
     #+ACM(4,18)*PMPB(4,9)+ACM(5,18)*PMPB(5,9)+ACM(6,18)*PMPB(6,9)+ACM(7
     #,18)*PMPB(7,9)

      ST(3) = ST(3)+t8*ACB(1,3)+t17*ACB(2,3)+t26*ACB(3,3)+t35*ACB(4,3)+t
     #44*ACB(5,3)+t53*ACB(6,3)+t62*ACB(7,3)+t71*ACB(8,3)+t80*ACB(9,3)
      ST(4) = ST(4)+t8*ACB(1,4)+t17*ACB(2,4)+t26*ACB(3,4)+t35*ACB(4,4)+t
     #44*ACB(5,4)+t53*ACB(6,4)+t62*ACB(7,4)+t71*ACB(8,4)+t80*ACB(9,4)
      ST(5) = ST(5)+t8*ACB(1,5)+t17*ACB(2,5)+t26*ACB(3,5)+t35*ACB(4,5)+t
     #44*ACB(5,5)+t53*ACB(6,5)+t62*ACB(7,5)+t71*ACB(8,5)+t80*ACB(9,5)
      ST(9) = ST(9)+t8*ACB(1,9)+t17*ACB(2,9)+t26*ACB(3,9)+t35*ACB(4,9)+t
     #44*ACB(5,9)+t53*ACB(6,9)+t62*ACB(7,9)+t71*ACB(8,9)+t80*ACB(9,9)
      ST(10) = ST(10)+t8*ACB(1,10)+t17*ACB(2,10)+t26*ACB(3,10)+t35*ACB(4
     #,10)+t44*ACB(5,10)+t53*ACB(6,10)+t62*ACB(7,10)+t71*ACB(8,10)+t80*A
     #CB(9,10)
      ST(11) = ST(11)+t8*ACB(1,11)+t17*ACB(2,11)+t26*ACB(3,11)+t35*ACB(4
     #,11)+t44*ACB(5,11)+t53*ACB(6,11)+t62*ACB(7,11)+t71*ACB(8,11)+t80*A
     #CB(9,11)
      ST(15) = ST(15)+t8*ACB(1,15)+t17*ACB(2,15)+t26*ACB(3,15)+t35*ACB(4
     #,15)+t44*ACB(5,15)+t53*ACB(6,15)+t62*ACB(7,15)+t71*ACB(8,15)+t80*A
     #CB(9,15)
      ST(16) = ST(16)+t8*ACB(1,16)+t17*ACB(2,16)+t26*ACB(3,16)+t35*ACB(4
     #,16)+t44*ACB(5,16)+t53*ACB(6,16)+t62*ACB(7,16)+t71*ACB(8,16)+t80*A
     #CB(9,16)
      ST(17) = ST(17)+t8*ACB(1,17)+t17*ACB(2,17)+t26*ACB(3,17)+t35*ACB(4
     #,17)+t44*ACB(5,17)+t53*ACB(6,17)+t62*ACB(7,17)+t71*ACB(8,17)+t80*A
     #CB(9,17)
      ST(20) = ST(20)+t170*ACB(1,3)+t179*ACB(2,3)+t188*ACB(3,3)+t197*ACB
     #(4,3)+t206*ACB(5,3)+t215*ACB(6,3)+t224*ACB(7,3)+t233*ACB(8,3)+t242
     #*ACB(9,3)
      ST(21) = ST(21)+t170*ACB(1,4)+t179*ACB(2,4)+t188*ACB(3,4)+t197*ACB
     #(4,4)+t206*ACB(5,4)+t215*ACB(6,4)+t224*ACB(7,4)+t233*ACB(8,4)+t242
     #*ACB(9,4)
      ST(22) = ST(22)+t170*ACB(1,5)+t179*ACB(2,5)+t188*ACB(3,5)+t197*ACB
     #(4,5)+t206*ACB(5,5)+t215*ACB(6,5)+t224*ACB(7,5)+t233*ACB(8,5)+t242
     #*ACB(9,5)
      ST(26) = ST(26)+t170*ACB(1,9)+t179*ACB(2,9)+t188*ACB(3,9)+t197*ACB
     #(4,9)+t206*ACB(5,9)+t215*ACB(6,9)+t224*ACB(7,9)+t233*ACB(8,9)+t242
     #*ACB(9,9)
      ST(27) = ST(27)+t170*ACB(1,10)+t179*ACB(2,10)+t188*ACB(3,10)+t197*
     #ACB(4,10)+t206*ACB(5,10)+t215*ACB(6,10)+t224*ACB(7,10)+t233*ACB(8,
     #10)+t242*ACB(9,10)
      ST(28) = ST(28)+t170*ACB(1,11)+t179*ACB(2,11)+t188*ACB(3,11)+t197*
     #ACB(4,11)+t206*ACB(5,11)+t215*ACB(6,11)+t224*ACB(7,11)+t233*ACB(8,
     #11)+t242*ACB(9,11)
      ST(32) = ST(32)+t170*ACB(1,15)+t179*ACB(2,15)+t188*ACB(3,15)+t197*
     #ACB(4,15)+t206*ACB(5,15)+t215*ACB(6,15)+t224*ACB(7,15)+t233*ACB(8,
     #15)+t242*ACB(9,15)
      ST(33) = ST(33)+t170*ACB(1,16)+t179*ACB(2,16)+t188*ACB(3,16)+t197*
     #ACB(4,16)+t206*ACB(5,16)+t215*ACB(6,16)+t224*ACB(7,16)+t233*ACB(8,
     #16)+t242*ACB(9,16)
      ST(34) = ST(34)+t170*ACB(1,17)+t179*ACB(2,17)+t188*ACB(3,17)+t197*
     #ACB(4,17)+t206*ACB(5,17)+t215*ACB(6,17)+t224*ACB(7,17)+t233*ACB(8,
     #17)+t242*ACB(9,17)
      ST(39) = ST(39)+t332*ACB(1,3)+t341*ACB(2,3)+t350*ACB(3,3)+t359*ACB
     #(4,3)+t368*ACB(5,3)+t377*ACB(6,3)+t386*ACB(7,3)+t395*ACB(8,3)+t404
     #*ACB(9,3)
      ST(40) = ST(40)+t414*ACB(1,3)+t423*ACB(2,3)+t432*ACB(3,3)+t441*ACB
     #(4,3)+t450*ACB(5,3)+t459*ACB(6,3)+t468*ACB(7,3)+t477*ACB(8,3)+t486
     #*ACB(9,3)
      ST(41) = ST(41)+t496*ACB(1,3)+t505*ACB(2,3)+t514*ACB(3,3)+t523*ACB
     #(4,3)+t532*ACB(5,3)+t541*ACB(6,3)+t550*ACB(7,3)+t559*ACB(8,3)+t568
     #*ACB(9,3)
      ST(45) = ST(45)+t578*ACB(1,3)+t587*ACB(2,3)+t596*ACB(3,3)+t605*ACB
     #(4,3)+t614*ACB(5,3)+t623*ACB(6,3)+t632*ACB(7,3)+t641*ACB(8,3)+t650
     #*ACB(9,3)
      ST(46) = ST(46)+t660*ACB(1,3)+t669*ACB(2,3)+t678*ACB(3,3)+t687*ACB
     #(4,3)+t696*ACB(5,3)+t705*ACB(6,3)+t714*ACB(7,3)+t723*ACB(8,3)+t732
     #*ACB(9,3)
      ST(47) = ST(47)+t742*ACB(1,3)+t751*ACB(2,3)+t760*ACB(3,3)+t769*ACB
     #(4,3)+t778*ACB(5,3)+t787*ACB(6,3)+t796*ACB(7,3)+t805*ACB(8,3)+t814
     #*ACB(9,3)
      ST(51) = ST(51)+t824*ACB(1,3)+t833*ACB(2,3)+t842*ACB(3,3)+t851*ACB
     #(4,3)+t860*ACB(5,3)+t869*ACB(6,3)+t878*ACB(7,3)+t887*ACB(8,3)+t896
     #*ACB(9,3)
      ST(54) = ST(54)+t332*ACB(1,4)+t341*ACB(2,4)+t350*ACB(3,4)+t359*ACB
     #(4,4)+t368*ACB(5,4)+t377*ACB(6,4)+t386*ACB(7,4)+t395*ACB(8,4)+t404
     #*ACB(9,4)
      ST(55) = ST(55)+t414*ACB(1,4)+t423*ACB(2,4)+t432*ACB(3,4)+t441*ACB
     #(4,4)+t450*ACB(5,4)+t459*ACB(6,4)+t468*ACB(7,4)+t477*ACB(8,4)+t486
     #*ACB(9,4)
      ST(56) = ST(56)+t496*ACB(1,4)+t505*ACB(2,4)+t514*ACB(3,4)+t523*ACB
     #(4,4)+t532*ACB(5,4)+t541*ACB(6,4)+t550*ACB(7,4)+t559*ACB(8,4)+t568
     #*ACB(9,4)
      ST(60) = ST(60)+t578*ACB(1,4)+t587*ACB(2,4)+t596*ACB(3,4)+t605*ACB
     #(4,4)+t614*ACB(5,4)+t623*ACB(6,4)+t632*ACB(7,4)+t641*ACB(8,4)+t650
     #*ACB(9,4)
      ST(61) = ST(61)+t660*ACB(1,4)+t669*ACB(2,4)+t678*ACB(3,4)+t687*ACB
     #(4,4)+t696*ACB(5,4)+t705*ACB(6,4)+t714*ACB(7,4)+t723*ACB(8,4)+t732
     #*ACB(9,4)
      ST(62) = ST(62)+t742*ACB(1,4)+t751*ACB(2,4)+t760*ACB(3,4)+t769*ACB
     #(4,4)+t778*ACB(5,4)+t787*ACB(6,4)+t796*ACB(7,4)+t805*ACB(8,4)+t814
     #*ACB(9,4)
      ST(66) = ST(66)+t824*ACB(1,4)+t833*ACB(2,4)+t842*ACB(3,4)+t851*ACB
     #(4,4)+t860*ACB(5,4)+t869*ACB(6,4)+t878*ACB(7,4)+t887*ACB(8,4)+t896
     #*ACB(9,4)
      ST(68) = ST(68)+t332*ACB(1,5)+t341*ACB(2,5)+t350*ACB(3,5)+t359*ACB
     #(4,5)+t368*ACB(5,5)+t377*ACB(6,5)+t386*ACB(7,5)+t395*ACB(8,5)+t404
     #*ACB(9,5)
      ST(69) = ST(69)+t414*ACB(1,5)+t423*ACB(2,5)+t432*ACB(3,5)+t441*ACB
     #(4,5)+t450*ACB(5,5)+t459*ACB(6,5)+t468*ACB(7,5)+t477*ACB(8,5)+t486
     #*ACB(9,5)
      ST(70) = ST(70)+t496*ACB(1,5)+t505*ACB(2,5)+t514*ACB(3,5)+t523*ACB
     #(4,5)+t532*ACB(5,5)+t541*ACB(6,5)+t550*ACB(7,5)+t559*ACB(8,5)+t568
     #*ACB(9,5)
      ST(74) = ST(74)+t578*ACB(1,5)+t587*ACB(2,5)+t596*ACB(3,5)+t605*ACB
     #(4,5)+t614*ACB(5,5)+t623*ACB(6,5)+t632*ACB(7,5)+t641*ACB(8,5)+t650
     #*ACB(9,5)
      ST(75) = ST(75)+t660*ACB(1,5)+t669*ACB(2,5)+t678*ACB(3,5)+t687*ACB
     #(4,5)+t696*ACB(5,5)+t705*ACB(6,5)+t714*ACB(7,5)+t723*ACB(8,5)+t732
     #*ACB(9,5)
      ST(76) = ST(76)+t742*ACB(1,5)+t751*ACB(2,5)+t760*ACB(3,5)+t769*ACB
     #(4,5)+t778*ACB(5,5)+t787*ACB(6,5)+t796*ACB(7,5)+t805*ACB(8,5)+t814
     #*ACB(9,5)
      ST(80) = ST(80)+t824*ACB(1,5)+t833*ACB(2,5)+t842*ACB(3,5)+t851*ACB
     #(4,5)+t860*ACB(5,5)+t869*ACB(6,5)+t878*ACB(7,5)+t887*ACB(8,5)+t896
     #*ACB(9,5)
      ST(84) = ST(84)+t332*ACB(1,9)+t341*ACB(2,9)+t350*ACB(3,9)+t359*ACB
     #(4,9)+t368*ACB(5,9)+t377*ACB(6,9)+t386*ACB(7,9)+t395*ACB(8,9)+t404
     #*ACB(9,9)
      ST(85) = ST(85)+t332*ACB(1,10)+t341*ACB(2,10)+t350*ACB(3,10)+t359*
     #ACB(4,10)+t368*ACB(5,10)+t377*ACB(6,10)+t386*ACB(7,10)+t395*ACB(8,
     #10)+t404*ACB(9,10)
      ST(86) = ST(86)+t332*ACB(1,11)+t341*ACB(2,11)+t350*ACB(3,11)+t359*
     #ACB(4,11)+t368*ACB(5,11)+t377*ACB(6,11)+t386*ACB(7,11)+t395*ACB(8,
     #11)+t404*ACB(9,11)
      ST(90) = ST(90)+t332*ACB(1,15)+t341*ACB(2,15)+t350*ACB(3,15)+t359*
     #ACB(4,15)+t368*ACB(5,15)+t377*ACB(6,15)+t386*ACB(7,15)+t395*ACB(8,
     #15)+t404*ACB(9,15)
      ST(91) = ST(91)+t332*ACB(1,16)+t341*ACB(2,16)+t350*ACB(3,16)+t359*
     #ACB(4,16)+t368*ACB(5,16)+t377*ACB(6,16)+t386*ACB(7,16)+t395*ACB(8,
     #16)+t404*ACB(9,16)
      ST(92) = ST(92)+t332*ACB(1,17)+t341*ACB(2,17)+t350*ACB(3,17)+t359*
     #ACB(4,17)+t368*ACB(5,17)+t377*ACB(6,17)+t386*ACB(7,17)+t395*ACB(8,
     #17)+t404*ACB(9,17)
      ST(96) = ST(96)+t414*ACB(1,9)+t423*ACB(2,9)+t432*ACB(3,9)+t441*ACB
     #(4,9)+t450*ACB(5,9)+t459*ACB(6,9)+t468*ACB(7,9)+t477*ACB(8,9)+t486
     #*ACB(9,9)
      ST(97) = ST(97)+t414*ACB(1,10)+t423*ACB(2,10)+t432*ACB(3,10)+t441*
     #ACB(4,10)+t450*ACB(5,10)+t459*ACB(6,10)+t468*ACB(7,10)+t477*ACB(8,
     #10)+t486*ACB(9,10)
      ST(98) = ST(98)+t414*ACB(1,11)+t423*ACB(2,11)+t432*ACB(3,11)+t441*
     #ACB(4,11)+t450*ACB(5,11)+t459*ACB(6,11)+t468*ACB(7,11)+t477*ACB(8,
     #11)+t486*ACB(9,11)
      ST(102) = ST(102)+t414*ACB(1,15)+t423*ACB(2,15)+t432*ACB(3,15)+t44
     #1*ACB(4,15)+t450*ACB(5,15)+t459*ACB(6,15)+t468*ACB(7,15)+t477*ACB(
     #8,15)+t486*ACB(9,15)
      ST(103) = ST(103)+t414*ACB(1,16)+t423*ACB(2,16)+t432*ACB(3,16)+t44
     #1*ACB(4,16)+t450*ACB(5,16)+t459*ACB(6,16)+t468*ACB(7,16)+t477*ACB(
     #8,16)+t486*ACB(9,16)
      ST(104) = ST(104)+t414*ACB(1,17)+t423*ACB(2,17)+t432*ACB(3,17)+t44
     #1*ACB(4,17)+t450*ACB(5,17)+t459*ACB(6,17)+t468*ACB(7,17)+t477*ACB(
     #8,17)+t486*ACB(9,17)
      ST(107) = ST(107)+t496*ACB(1,9)+t505*ACB(2,9)+t514*ACB(3,9)+t523*A
     #CB(4,9)+t532*ACB(5,9)+t541*ACB(6,9)+t550*ACB(7,9)+t559*ACB(8,9)+t5
     #68*ACB(9,9)
      ST(108) = ST(108)+t496*ACB(1,10)+t505*ACB(2,10)+t514*ACB(3,10)+t52
     #3*ACB(4,10)+t532*ACB(5,10)+t541*ACB(6,10)+t550*ACB(7,10)+t559*ACB(
     #8,10)+t568*ACB(9,10)
      ST(109) = ST(109)+t496*ACB(1,11)+t505*ACB(2,11)+t514*ACB(3,11)+t52
     #3*ACB(4,11)+t532*ACB(5,11)+t541*ACB(6,11)+t550*ACB(7,11)+t559*ACB(
     #8,11)+t568*ACB(9,11)
      ST(113) = ST(113)+t496*ACB(1,15)+t505*ACB(2,15)+t514*ACB(3,15)+t52
     #3*ACB(4,15)+t532*ACB(5,15)+t541*ACB(6,15)+t550*ACB(7,15)+t559*ACB(
     #8,15)+t568*ACB(9,15)
      ST(114) = ST(114)+t496*ACB(1,16)+t505*ACB(2,16)+t514*ACB(3,16)+t52
     #3*ACB(4,16)+t532*ACB(5,16)+t541*ACB(6,16)+t550*ACB(7,16)+t559*ACB(
     #8,16)+t568*ACB(9,16)
      ST(115) = ST(115)+t496*ACB(1,17)+t505*ACB(2,17)+t514*ACB(3,17)+t52
     #3*ACB(4,17)+t532*ACB(5,17)+t541*ACB(6,17)+t550*ACB(7,17)+t559*ACB(
     #8,17)+t568*ACB(9,17)
      ST(120) = ST(120)+t578*ACB(1,9)+t587*ACB(2,9)+t596*ACB(3,9)+t605*A
     #CB(4,9)+t614*ACB(5,9)+t623*ACB(6,9)+t632*ACB(7,9)+t641*ACB(8,9)+t6
     #50*ACB(9,9)
      ST(121) = ST(121)+t660*ACB(1,9)+t669*ACB(2,9)+t678*ACB(3,9)+t687*A
     #CB(4,9)+t696*ACB(5,9)+t705*ACB(6,9)+t714*ACB(7,9)+t723*ACB(8,9)+t7
     #32*ACB(9,9)
      ST(122) = ST(122)+t742*ACB(1,9)+t751*ACB(2,9)+t760*ACB(3,9)+t769*A
     #CB(4,9)+t778*ACB(5,9)+t787*ACB(6,9)+t796*ACB(7,9)+t805*ACB(8,9)+t8
     #14*ACB(9,9)
      ST(126) = ST(126)+t824*ACB(1,9)+t833*ACB(2,9)+t842*ACB(3,9)+t851*A
     #CB(4,9)+t860*ACB(5,9)+t869*ACB(6,9)+t878*ACB(7,9)+t887*ACB(8,9)+t8
     #96*ACB(9,9)
      ST(129) = ST(129)+t578*ACB(1,10)+t587*ACB(2,10)+t596*ACB(3,10)+t60
     #5*ACB(4,10)+t614*ACB(5,10)+t623*ACB(6,10)+t632*ACB(7,10)+t641*ACB(
     #8,10)+t650*ACB(9,10)
      ST(130) = ST(130)+t660*ACB(1,10)+t669*ACB(2,10)+t678*ACB(3,10)+t68
     #7*ACB(4,10)+t696*ACB(5,10)+t705*ACB(6,10)+t714*ACB(7,10)+t723*ACB(
     #8,10)+t732*ACB(9,10)
      ST(131) = ST(131)+t742*ACB(1,10)+t751*ACB(2,10)+t760*ACB(3,10)+t76
     #9*ACB(4,10)+t778*ACB(5,10)+t787*ACB(6,10)+t796*ACB(7,10)+t805*ACB(
     #8,10)+t814*ACB(9,10)
      ST(135) = ST(135)+t824*ACB(1,10)+t833*ACB(2,10)+t842*ACB(3,10)+t85
     #1*ACB(4,10)+t860*ACB(5,10)+t869*ACB(6,10)+t878*ACB(7,10)+t887*ACB(
     #8,10)+t896*ACB(9,10)
      ST(137) = ST(137)+t578*ACB(1,11)+t587*ACB(2,11)+t596*ACB(3,11)+t60
     #5*ACB(4,11)+t614*ACB(5,11)+t623*ACB(6,11)+t632*ACB(7,11)+t641*ACB(
     #8,11)+t650*ACB(9,11)
      ST(138) = ST(138)+t660*ACB(1,11)+t669*ACB(2,11)+t678*ACB(3,11)+t68
     #7*ACB(4,11)+t696*ACB(5,11)+t705*ACB(6,11)+t714*ACB(7,11)+t723*ACB(
     #8,11)+t732*ACB(9,11)
      ST(139) = ST(139)+t742*ACB(1,11)+t751*ACB(2,11)+t760*ACB(3,11)+t76
     #9*ACB(4,11)+t778*ACB(5,11)+t787*ACB(6,11)+t796*ACB(7,11)+t805*ACB(
     #8,11)+t814*ACB(9,11)
      ST(143) = ST(143)+t824*ACB(1,11)+t833*ACB(2,11)+t842*ACB(3,11)+t85
     #1*ACB(4,11)+t860*ACB(5,11)+t869*ACB(6,11)+t878*ACB(7,11)+t887*ACB(
     #8,11)+t896*ACB(9,11)
      ST(147) = ST(147)+t578*ACB(1,15)+t587*ACB(2,15)+t596*ACB(3,15)+t60
     #5*ACB(4,15)+t614*ACB(5,15)+t623*ACB(6,15)+t632*ACB(7,15)+t641*ACB(
     #8,15)+t650*ACB(9,15)
      ST(148) = ST(148)+t578*ACB(1,16)+t587*ACB(2,16)+t596*ACB(3,16)+t60
     #5*ACB(4,16)+t614*ACB(5,16)+t623*ACB(6,16)+t632*ACB(7,16)+t641*ACB(
     #8,16)+t650*ACB(9,16)
      ST(149) = ST(149)+t578*ACB(1,17)+t587*ACB(2,17)+t596*ACB(3,17)+t60
     #5*ACB(4,17)+t614*ACB(5,17)+t623*ACB(6,17)+t632*ACB(7,17)+t641*ACB(
     #8,17)+t650*ACB(9,17)
      ST(153) = ST(153)+t660*ACB(1,15)+t669*ACB(2,15)+t678*ACB(3,15)+t68
     #7*ACB(4,15)+t696*ACB(5,15)+t705*ACB(6,15)+t714*ACB(7,15)+t723*ACB(
     #8,15)+t732*ACB(9,15)
      ST(154) = ST(154)+t660*ACB(1,16)+t669*ACB(2,16)+t678*ACB(3,16)+t68
     #7*ACB(4,16)+t696*ACB(5,16)+t705*ACB(6,16)+t714*ACB(7,16)+t723*ACB(
     #8,16)+t732*ACB(9,16)
      ST(155) = ST(155)+t660*ACB(1,17)+t669*ACB(2,17)+t678*ACB(3,17)+t68
     #7*ACB(4,17)+t696*ACB(5,17)+t705*ACB(6,17)+t714*ACB(7,17)+t723*ACB(
     #8,17)+t732*ACB(9,17)
      ST(158) = ST(158)+t742*ACB(1,15)+t751*ACB(2,15)+t760*ACB(3,15)+t76
     #9*ACB(4,15)+t778*ACB(5,15)+t787*ACB(6,15)+t796*ACB(7,15)+t805*ACB(
     #8,15)+t814*ACB(9,15)
      ST(159) = ST(159)+t742*ACB(1,16)+t751*ACB(2,16)+t760*ACB(3,16)+t76
     #9*ACB(4,16)+t778*ACB(5,16)+t787*ACB(6,16)+t796*ACB(7,16)+t805*ACB(
     #8,16)+t814*ACB(9,16)
      ST(160) = ST(160)+t742*ACB(1,17)+t751*ACB(2,17)+t760*ACB(3,17)+t76
     #9*ACB(4,17)+t778*ACB(5,17)+t787*ACB(6,17)+t796*ACB(7,17)+t805*ACB(
     #8,17)+t814*ACB(9,17)
      ST(165) = ST(165)+t824*ACB(1,15)+t833*ACB(2,15)+t842*ACB(3,15)+t85
     #1*ACB(4,15)+t860*ACB(5,15)+t869*ACB(6,15)+t878*ACB(7,15)+t887*ACB(
     #8,15)+t896*ACB(9,15)
      ST(168) = ST(168)+t824*ACB(1,16)+t833*ACB(2,16)+t842*ACB(3,16)+t85
     #1*ACB(4,16)+t860*ACB(5,16)+t869*ACB(6,16)+t878*ACB(7,16)+t887*ACB(
     #8,16)+t896*ACB(9,16)
      ST(170) = ST(170)+t824*ACB(1,17)+t833*ACB(2,17)+t842*ACB(3,17)+t85
     #1*ACB(4,17)+t860*ACB(5,17)+t869*ACB(6,17)+t878*ACB(7,17)+t887*ACB(
     #8,17)+t896*ACB(9,17)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE NFPPL_TRI(PAPM,PAPB,PAPSH,PFP)

C=============================================================
C	PURPOSE: TO COMPUTE FOR THE INTEGRAL(TRANSPOSE(Pg)*F*Pg)
C
C	INPUT VARIABLES
C		PAPM(7)   = INTEGRAL(TRANS(Pm)*N)
C		PAPB(9)  = INTEGRAL(TRANS(Pb)*M)
C		PAPSH(2)  = INTEGRAL(TRANS(Ps)*Q)
C
C	LOCAL VARIABLES
C		NR,NS,NRS = INTEGRAL(N)dA
C		MR,MS,MRS = INTEGRAL(M)dA
C		QR,QS	  = INTEGRAL(Q)dA
C
C	OUTPUT VARIABLE
C		PFP(15,15) = INTEGRAL(TRANSPOSE(Pg)*F*Pg)
C=============================================================
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION PAPM(7),PAPB(9),PAPSH(2)
	DIMENSION PFP(15,15)

C	-------------------------------------
C	ASSEMBLE INTEGRAL(TRANSPOSE(Pg)*F*Pg)
C	-------------------------------------
      PFP(1,1) = PAPM(1)
      PFP(1,2) = PAPM(7)
      PFP(1,9) = PAPB(1)
      PFP(1,10) = PAPB(7)
      PFP(1,14) = PAPSH(1)

      PFP(2,1) = PAPM(7)
      PFP(2,2) = PAPM(4)
      PFP(2,9) = PAPB(7)
      PFP(2,10) = PAPB(4)
      PFP(2,14) = PAPSH(2)

      PFP(3,3) = PAPM(1)
      PFP(3,4) = PAPM(7)
      PFP(3,7) = -PAPB(1)
      PFP(3,8) = -PAPB(7)
      PFP(3,13) = -PAPSH(1)

      PFP(4,3) = PAPM(7)
      PFP(4,4) = PAPM(4)
      PFP(4,7) = -PAPB(7)
      PFP(4,8) = -PAPB(4)
      PFP(4,13) = -PAPSH(2)

      PFP(5,5) = PAPM(1)
      PFP(5,6) = PAPM(7)

      PFP(6,5) = PAPM(7)
      PFP(6,6) = PAPM(4)

      PFP(7,3) = -PAPB(1)
      PFP(7,4) = -PAPB(7)
      PFP(7,15) = PAPB(1)/2.0D0

      PFP(8,3) = -PAPB(7)
      PFP(8,4) = -PAPB(4)
      PFP(8,15) = PAPB(7)/2.0D0

      PFP(9,1) = PAPB(1)
      PFP(9,2) = PAPB(7)
      PFP(9,15) = PAPB(7)/2.0D0

      PFP(10,1) = PAPB(7)
      PFP(10,2) = PAPB(4)
      PFP(10,15) = PAPB(4)/2.0D0

      PFP(11,13) = PAPB(1)/2.0D0
      PFP(11,14) = PAPB(7)/2.0D0

      PFP(12,13) = PAPB(7)/2.0D0
      PFP(12,14) = PAPB(4)/2.0D0

      PFP(13,3) = -PAPSH(1)
      PFP(13,4) = -PAPSH(2)
      PFP(13,11) = PAPB(1)/2.0D0
      PFP(13,12) = PAPB(7)/2.0D0
      PFP(13,15) = PAPSH(1)/2.0D0

      PFP(14,1) = PAPSH(1)
      PFP(14,2) = PAPSH(2)
      PFP(14,11) = PAPB(7)/2.0D0
      PFP(14,12) = PAPB(4)/2.0D0
      PFP(14,15) = PAPSH(2)/2.0D0

      PFP(15,7) = PAPB(1)/2.0D0
      PFP(15,8) = PAPB(7)/2.0D0
      PFP(15,9) = PAPB(7)/2.0D0
      PFP(15,10) = PAPB(4)/2.0D0
      PFP(15,13) = PAPSH(1)/2.0D0
      PFP(15,14) = PAPSH(2)/2.0D0

	RETURN
	END
C
C=============================================================
      SUBROUTINE COMFACT_TRI(RS,ANTR,ANTS,ANTRS,ANTRR,ANTSS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------
C     PURPOSE:	COMPUTES FOR INTEGRATION FACTORS FOR INTEGRALS
C				USING ISOPARAMETRIC MAPPING
C
C	INPUT VARIABLES
C	RS(2,3)  = NODAL COORDINATES
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C     ANTR(3)  = FACTORS FOR INT(R*PHI)
C     ANTS(3)  = FACTORS FOR INT(S*PHI)
C     ANTRS(3) = FACTORS FOR INT(R*S*PHI)
C     ANTRR(3) = FACTORS FOR INT(R*R*PHI)
C     ANTSS(3) = FACTORS FOR INT(S*S*PHI)
C     ----------------------------------------------------------
	DIMENSION RS(2,3)
	DIMENSION ANTR(3),ANTS(3),ANTRS(3),ANTRR(3),ANTSS(3)

C	-----------------------
C	COMPUTE FOR FACTOR ANTR
C	-----------------------
	ANTR(1) = 
     #((RS(2,2)/12-RS(2,3)/12)*RS(1,1)**2+((RS(2,2)/24-RS(2,1)/12+R
     #S(2,3)/24)*RS(1,2)+(-RS(2,2)/24+RS(2,1)/12-RS(2,3)/24)*RS(1,3))*RS
     #(1,1)+(RS(2,3)/24-RS(2,1)/24)*RS(1,2)**2+(-RS(2,2)/24+RS(2,3)/24)*
     #RS(1,3)*RS(1,2)+(-RS(2,2)/24+RS(2,1)/24)*RS(1,3)**2)
	ANTR(2) = 
     #((-RS(2,3)/24+RS(2,2)/24)*RS(1,1)**2+((RS(2,2)/12-RS(2,1)/24-
     #RS(2,3)/24)*RS(1,2)+(RS(2,1)/24-RS(2,3)/24)*RS(1,3))*RS(1,1)+(-RS(
     #2,1)/12+RS(2,3)/12)*RS(1,2)**2+(RS(2,1)/24-RS(2,2)/12+RS(2,3)/24)*
     #RS(1,3)*RS(1,2)+(-RS(2,2)/24+RS(2,1)/24)*RS(1,3)**2)
	ANTR(3) = 
     #((-RS(
     #2,3)/24+RS(2,2)/24)*RS(1,1)**2+((RS(2,2)/24-RS(2,1)/24)*RS(1,2)+(R
     #S(2,2)/24-RS(2,3)/12+RS(2,1)/24)*RS(1,3))*RS(1,1)+(RS(2,3)/24-RS(2
     #,1)/24)*RS(1,2)**2+(-RS(2,1)/24-RS(2,2)/24+RS(2,3)/12)*RS(1,3)*RS(
     #1,2)+(RS(2,1)/12-RS(2,2)/12)*RS(1,3)**2)

C	-----------------------
C	COMPUTE FOR FACTOR ANTS
C	-----------------------
	ANTS(1) = 
     #((RS(2,2)**2/24+RS(2,1)*RS(2,2)/12-RS(2,3)**2/24-RS(2,1)*RS(2
     #,3)/12)*RS(1,1)+(-RS(2,1)**2/12-RS(2,1)*RS(2,2)/24+RS(2,1)*RS(2,3)
     #/24+RS(2,2)*RS(2,3)/24+RS(2,3)**2/24)*RS(1,2)+(-RS(2,2)*RS(2,3)/24
     #+RS(2,1)**2/12+RS(2,1)*RS(2,3)/24-RS(2,2)**2/24-RS(2,1)*RS(2,2)/24
     #)*RS(1,3))
	ANTS(2) = 
     #((RS(2,2)**2/12-RS(2,1)*RS(2,3)/24-RS(2,2)*RS(2,3)/24+RS(2,1)
     #*RS(2,2)/24-RS(2,3)**2/24)*RS(1,1)+(-RS(2,1)*RS(2,2)/12+RS(2,3)**2
     #/24-RS(2,1)**2/24+RS(2,2)*RS(2,3)/12)*RS(1,2)+(RS(2,1)*RS(2,3)/24-
     #RS(2,2)**2/12-RS(2,2)*RS(2,3)/24+RS(2,1)*RS(2,2)/24+RS(2,1)**2/24)
     #*RS(1,3))
	ANTS(3) = 
     #((-RS(2,1)*RS(2,3)/24+RS(2,2)**2/24+RS(2,1)*RS(2,
     #2)/24+RS(2,2)*RS(2,3)/24-RS(2,3)**2/12)*RS(1,1)+(RS(2,3)**2/12-RS(
     #2,1)*RS(2,3)/24-RS(2,1)**2/24-RS(2,1)*RS(2,2)/24+RS(2,2)*RS(2,3)/2
     #4)*RS(1,2)+(RS(2,1)*RS(2,3)/12-RS(2,2)**2/24+RS(2,1)**2/24-RS(2,2)
     #*RS(2,3)/12)*RS(1,3))

C	------------------------
C	COMPUTE FOR FACTOR ANTRS
C	------------------------
      ANTRS(1) = 
     #((RS(2,1)*RS(2,2)/20-RS(2,3)**2/60+RS(2,2)**2/60-RS(2,1)*RS(2
     #,3)/20)*RS(1,1)**2+((RS(2,2)*RS(2,3)/120+RS(2,3)**2/120+RS(2,2)**2
     #/60-RS(2,1)**2/20+RS(2,1)*RS(2,3)/60)*RS(1,2)+(-RS(2,1)*RS(2,2)/60
     #-RS(2,2)**2/120+RS(2,1)**2/20-RS(2,3)**2/60-RS(2,2)*RS(2,3)/120)*R
     #S(1,3))*RS(1,1)+(-RS(2,1)**2/60+RS(2,1)*RS(2,3)/120-RS(2,1)*RS(2,2
     #)/60+RS(2,2)*RS(2,3)/60+RS(2,3)**2/120)*RS(1,2)**2+(RS(2,1)*RS(2,3
     #)/120-RS(2,2)**2/60+RS(2,3)**2/60-RS(2,1)*RS(2,2)/120)*RS(1,3)*RS(
     #1,2)+(-RS(2,2)**2/120+RS(2,1)**2/60-RS(2,1)*RS(2,2)/120-RS(2,2)*RS
     #(2,3)/60+RS(2,1)*RS(2,3)/60)*RS(1,3)**2)
      ANTRS(2) = 
     #((RS(2,1)*RS(2,2)/60-RS(2,3)**2/120+RS(2,2)**2/60-RS(2,1)*RS(
     #2,3)/60-RS(2,2)*RS(2,3)/120)*RS(1,1)**2+((-RS(2,1)*RS(2,3)/120-RS(
     #2,2)*RS(2,3)/60-RS(2,3)**2/120+RS(2,2)**2/20-RS(2,1)**2/60)*RS(1,2
     #)+(RS(2,1)*RS(2,2)/120+RS(2,1)**2/60-RS(2,3)**2/60-RS(2,2)*RS(2,3)
     #/120)*RS(1,3))*RS(1,1)+(-RS(2,1)**2/60+RS(2,3)**2/60-RS(2,1)*RS(2,
     #2)/20+RS(2,2)*RS(2,3)/20)*RS(1,2)**2+(-RS(2,2)**2/20+RS(2,1)*RS(2,
     #2)/60+RS(2,1)*RS(2,3)/120+RS(2,3)**2/60+RS(2,1)**2/120)*RS(1,3)*RS
     #(1,2)+(-RS(2,2)**2/60+RS(2,1)**2/120-RS(2,2)*RS(2,3)/60+RS(2,1)*RS
     #(2,2)/120+RS(2,1)*RS(2,3)/60)*RS(1,3)**2)
      ANTRS(3) = 
     #((RS(2,2)*RS(2,3)/120+RS(2,2)**2/120-RS(2,3)**2/60+RS(2,1)*RS
     #(2,2)/60-RS(2,1)*RS(2,3)/60)*RS(1,1)**2+((-RS(2,1)**2/60-RS(2,1)*R
     #S(2,3)/120+RS(2,2)**2/60+RS(2,2)*RS(2,3)/120)*RS(1,2)+(RS(2,2)*RS(
     #2,3)/60+RS(2,1)*RS(2,2)/120+RS(2,1)**2/60-RS(2,3)**2/20+RS(2,2)**2
     #/120)*RS(1,3))*RS(1,1)+(-RS(2,1)**2/120+RS(2,3)**2/60-RS(2,1)*RS(2
     #,3)/120-RS(2,1)*RS(2,2)/60+RS(2,2)*RS(2,3)/60)*RS(1,2)**2+(-RS(2,2
     #)**2/60-RS(2,1)*RS(2,3)/60-RS(2,1)*RS(2,2)/120-RS(2,1)**2/120+RS(2
     #,3)**2/20)*RS(1,3)*RS(1,2)+(-RS(2,2)**2/60+RS(2,1)**2/60+RS(2,1)*R
     #S(2,3)/20-RS(2,2)*RS(2,3)/20)*RS(1,3)**2)

C	------------------------
C	COMPUTE FOR FACTOR ANTRR
C	------------------------
      ANTRR(1) = 
     #((-RS(2,3)/20+RS(2,2)/20)*RS(1,1)**3+((RS(2,2)/30+RS(2,3)/60-
     #RS(2,1)/20)*RS(1,2)+(RS(2,1)/20-RS(2,3)/30-RS(2,2)/60)*RS(1,3))*RS
     #(1,1)**2+((RS(2,2)/60+RS(2,3)/60-RS(2,1)/30)*RS(1,2)**2+(RS(2,3)/6
     #0-RS(2,2)/60)*RS(1,3)*RS(1,2)+(RS(2,1)/30-RS(2,2)/60-RS(2,3)/60)*R
     #S(1,3)**2)*RS(1,1)+(-RS(2,1)/60+RS(2,3)/60)*RS(1,2)**3+(RS(2,3)/60
     #-RS(2,2)/60)*RS(1,3)*RS(1,2)**2+(RS(2,3)/60-RS(2,2)/60)*RS(1,3)**2
     #*RS(1,2)+(-RS(2,2)/60+RS(2,1)/60)*RS(1,3)**3)
      ANTRR(2) = 
     #((-RS(2,3)/60+RS(2,2)/60)*RS(1,1)**3+((RS(2,2)/30-RS(2,3)/60-
     #RS(2,1)/60)*RS(1,2)+(RS(2,1)/60-RS(2,3)/60)*RS(1,3))*RS(1,1)**2+((
     #-RS(2,3)/60+RS(2,2)/20-RS(2,1)/30)*RS(1,2)**2+(RS(2,1)/60-RS(2,3)/
     #60)*RS(1,3)*RS(1,2)+(RS(2,1)/60-RS(2,3)/60)*RS(1,3)**2)*RS(1,1)+(R
     #S(2,3)/20-RS(2,1)/20)*RS(1,2)**3+(-RS(2,2)/20+RS(2,1)/60+RS(2,3)/3
     #0)*RS(1,3)*RS(1,2)**2+(RS(2,1)/60-RS(2,2)/30+RS(2,3)/60)*RS(1,3)**
     #2*RS(1,2)+(-RS(2,2)/60+RS(2,1)/60)*RS(1,3)**3)
      ANTRR(3) = 
     #((-RS(2,3)/60+RS(2,2)/60)*RS(1,1)**3+((RS(2,2)/60-RS(2,1)/60)
     #*RS(1,2)+(RS(2,2)/60+RS(2,1)/60-RS(2,3)/30)*RS(1,3))*RS(1,1)**2+((
     #RS(2,2)/60-RS(2,1)/60)*RS(1,2)**2+(RS(2,2)/60-RS(2,1)/60)*RS(1,3)*
     #RS(1,2)+(-RS(2,3)/20+RS(2,1)/30+RS(2,2)/60)*RS(1,3)**2)*RS(1,1)+(-
     #RS(2,1)/60+RS(2,3)/60)*RS(1,2)**3+(RS(2,3)/30-RS(2,2)/60-RS(2,1)/6
     #0)*RS(1,3)*RS(1,2)**2+(-RS(2,2)/30+RS(2,3)/20-RS(2,1)/60)*RS(1,3)*
     #*2*RS(1,2)+(RS(2,1)/20-RS(2,2)/20)*RS(1,3)**3)

C	------------------------
C	COMPUTE FOR FACTOR ANTSS
C	------------------------
      ANTSS(1) = 
     #((RS(2,1)**2*RS(2,2)/20+RS(2,2)**3/60-RS(2,1)*RS(2,3)**2/30-R
     #S(2,1)**2*RS(2,3)/20-RS(2,3)**3/60+RS(2,1)*RS(2,2)**2/30)*RS(1,1)+
     #(RS(2,1)**2*RS(2,3)/60+RS(2,1)*RS(2,2)*RS(2,3)/60+RS(2,2)**2*RS(2,
     #3)/60+RS(2,2)*RS(2,3)**2/60+RS(2,3)**3/60-RS(2,1)**3/20+RS(2,1)*RS
     #(2,3)**2/60-RS(2,1)**2*RS(2,2)/30-RS(2,1)*RS(2,2)**2/60)*RS(1,2)+(
     #-RS(2,2)**2*RS(2,3)/60+RS(2,1)**3/20-RS(2,1)*RS(2,2)**2/60-RS(2,1)
     #*RS(2,2)*RS(2,3)/60+RS(2,1)**2*RS(2,3)/30-RS(2,2)*RS(2,3)**2/60-RS
     #(2,1)**2*RS(2,2)/60-RS(2,2)**3/60+RS(2,1)*RS(2,3)**2/60)*RS(1,3))
      ANTSS(2) = 
     #((-RS(2,1)**2*RS(2,3)/60-RS(2,2)*RS(2,3)**2/60+RS(2,1)**2*RS(
     #2,2)/60-RS(2,2)**2*RS(2,3)/60+RS(2,1)*RS(2,2)**2/30-RS(2,1)*RS(2,2
     #)*RS(2,3)/60+RS(2,2)**3/20-RS(2,3)**3/60-RS(2,1)*RS(2,3)**2/60)*RS
     #(1,1)+(RS(2,2)**2*RS(2,3)/20+RS(2,3)**3/60-RS(2,1)**2*RS(2,2)/30-R
     #S(2,1)**3/60+RS(2,2)*RS(2,3)**2/30-RS(2,1)*RS(2,2)**2/20)*RS(1,2)+
     #(RS(2,1)*RS(2,2)*RS(2,3)/60-RS(2,2)*RS(2,3)**2/60-RS(2,2)**2*RS(2,
     #3)/30-RS(2,2)**3/20+RS(2,1)**2*RS(2,3)/60+RS(2,1)**2*RS(2,2)/60+RS
     #(2,1)*RS(2,3)**2/60+RS(2,1)*RS(2,2)**2/60+RS(2,1)**3/60)*RS(1,3))
      ANTSS(3) = 
     #((RS(2,2)**3/60+RS(2,1)**2*RS(2,2)/60+RS(2,1)*RS(2,2)**2/60+R
     #S(2,2)*RS(2,3)**2/60-RS(2,1)*RS(2,3)**2/30+RS(2,2)**2*RS(2,3)/60-R
     #S(2,3)**3/20+RS(2,1)*RS(2,2)*RS(2,3)/60-RS(2,1)**2*RS(2,3)/60)*RS(
     #1,1)+(RS(2,2)**2*RS(2,3)/60-RS(2,1)**3/60-RS(2,1)*RS(2,2)*RS(2,3)/
     #60-RS(2,1)*RS(2,3)**2/60+RS(2,2)*RS(2,3)**2/30-RS(2,1)**2*RS(2,2)/
     #60-RS(2,1)**2*RS(2,3)/60+RS(2,3)**3/20-RS(2,1)*RS(2,2)**2/60)*RS(1
     #,2)+(-RS(2,2)**2*RS(2,3)/30-RS(2,2)*RS(2,3)**2/20-RS(2,2)**3/60+RS
     #(2,1)**2*RS(2,3)/30+RS(2,1)**3/60+RS(2,1)*RS(2,3)**2/20)*RS(1,3))

      RETURN
      END
C
C=====================================================================
      SUBROUTINE TRIVRS(COORD,VR,VS,VT,RS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------
C     PROGRAM TO DETERMINE LOCAL CORDINATE VECTORS - VR,VS,VT
C     -------------------------------------------------------
      DIMENSION COORD(3,3) ! COORD(XYZ,NNO) - INPUT
	DIMENSION V13(3),DUMMY(3) ! LOCAL VARIABLES
	DIMENSION VR(3),VS(3),VT(3),RS(2,3) !RS(RS,NNO)-OUTPUT VARIABLES
	DIMENSION EC(3),DUM1(3),DUM2(3)

C	----------------------------------------
C	COMPUTE FOR LOCAL COORDINATE VECTOR - VR
C	----------------------------------------
	VR(1) = COORD(1,2)-COORD(1,1)
	VR(2) = COORD(2,2)-COORD(2,1)
	VR(3) = COORD(3,2)-COORD(3,1)
	CALL SCALEN(VR,VR,DUM,3)

C	---------------------------------
C	COMPUTE FOR SIDE 1-3 VECTOR - V13
C	---------------------------------
	V13(1) = COORD(1,3)-COORD(1,1)
	V13(2) = COORD(2,3)-COORD(2,1)
	V13(3) = COORD(3,3)-COORD(3,1)
	CALL SCALEN(V13,V13,DUM,3)

C	-----------------------------------------------
C	COMPUTE FOR LOCAL COORDINATE VECTORS - VR,VS,VT
C	-----------------------------------------------
	CALL VECPRD(VR,V13,VT)
	CALL SCALEN(VT,VT,DUM,3)

	CALL VECPRD(VT,VR,VS)
	CALL SCALEN(VS,VS,DUM,3)

C	------------------------------
C	COMPUTE FOR THE ELEMENT CENTER
C	------------------------------
	EC(1) = COORD(1,1)
	EC(2) = COORD(2,1)
	EC(3) = COORD(3,1)

C	-----------------------------------------------
C	DETERMINE NODAL COORDINATES IN LOCAL SPACE - RS
C	-----------------------------------------------
	DO 20 I=1,3
	 DO 30 J=1,3
	  DUMMY(J)=COORD(J,I)-EC(J)
30	 CONTINUE
	 RS(1,I)=DUMMY(1)*VR(1)+DUMMY(2)*VR(2)+DUMMY(3)*VR(3)
	 RS(2,I)=DUMMY(1)*VS(1)+DUMMY(2)*VS(2)+DUMMY(3)*VS(3)
20	CONTINUE

      RETURN
      END
C
C=====================================================================
      SUBROUTINE TRIVRS2(COORD,VR,VS,VT,RS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------
C     PROGRAM TO DETERMINE LOCAL CORDINATE VECTORS - VR,VS,VT
C     -------------------------------------------------------
      DIMENSION COORD(3,3)
	DIMENSION XB(3,3),EC(3),VL1(3),VL2(3),DUMMY(3)
	DIMENSION VR(3),VS(3),VT(3),RS(2,3),T(4)
	DIMENSION V15(3)

C	------------------------------------------------------------------
C	COMPUTE FOR VECTOR PASSING THROUGH NODE 1 AND MIDPOINT OF SIDE 2-3
C	------------------------------------------------------------------
	V15(1) = (COORD(1,2)+COORD(1,3))/2.0d0 - COORD(1,1)
	V15(2) = (COORD(2,2)+COORD(2,3))/2.0d0 - COORD(2,1)
	V15(3) = (COORD(3,2)+COORD(3,3))/2.0d0 - COORD(3,1)
	CALL SCALEN(V15,V15,DUM,3)

C	-------------------------------------------
C	COMPUTE FOR ELEMENT CENTER COORDINATES - EC
C	-------------------------------------------
	EC(1) = COORD(1,1) + V15(1)*DUM*(2.0D0/3.0D0)
	EC(2) = COORD(2,1) + V15(2)*DUM*(2.0D0/3.0D0)
	EC(3) = COORD(3,1) + V15(3)*DUM*(2.0D0/3.0D0)

C	--------------------------------------
C	COMPUTE FOR MID-POINT COORDINATES - XB
C	--------------------------------------
	XB(1,1)=(COORD(1,1)+COORD(1,3))/2.0d0
	XB(2,1)=(COORD(2,1)+COORD(2,3))/2.0d0
	XB(3,1)=(COORD(3,1)+COORD(3,3))/2.0d0
	DO 10 I=2,NNO
	 XB(1,I)=(COORD(1,I)+COORD(1,I-1))/2.0d0
	 XB(2,I)=(COORD(2,I)+COORD(2,I-1))/2.0d0
	 XB(3,I)=(COORD(3,I)+COORD(3,I-1))/2.0d0
10	CONTINUE

C	----------------------------------------------
C	COMPUTE FOR OPPOSITE MID-POINT VECTORS - L1,L2
C	----------------------------------------------
	VL1(1)=COORD(1,2)-EC(1)
	VL1(2)=COORD(2,2)-EC(2)
	VL1(3)=COORD(3,2)-EC(3)
	CALL SCALEN(VL1,VL1,DUM,3)

	VL2(1)=COORD(1,3)-EC(1)
	VL2(2)=COORD(2,3)-EC(2)
	VL2(3)=COORD(3,3)-EC(3)
	CALL SCALEN(VL2,VL2,DUM,3)

C	-----------------------------------------------
C	COMPUTE FOR LOCAL COORDINATE VECTORS - VR,VS,VT
C	-----------------------------------------------
	CALL VECPRD(VL1,VL2,VT)
	CALL SCALEN(VT,VT,DUM,3)

	CALL VECPRD(VT,VL1,VS)
	CALL ADDVEC(VS,VL2,VS)
	CALL SCALEN(VS,VS,DUM,3)

	CALL VECPRD(VS,VT,VR)
	CALL SCALEN(VR,VR,DUM,3)
C	-----------------------------------------------
C	DETERMINE NODAL COORDINATES IN LOCAL SPACE - RS
C	-----------------------------------------------
	DO 20 I=1,3 
	 DO 30 J=1,3
	  DUMMY(J)=COORD(J,I)-EC(J)
30	 CONTINUE
	 RS(1,I)=DUMMY(1)*VR(1)+DUMMY(2)*VR(2)+DUMMY(3)*VR(3)
	 RS(2,I)=DUMMY(1)*VS(1)+DUMMY(2)*VS(2)+DUMMY(3)*VS(3)
	 T(I)   =DUMMY(1)*VT(1)+DUMMY(2)*VT(2)+DUMMY(3)*VT(3)
20	CONTINUE

      RETURN
      END
C
C=====================================================================
      SUBROUTINE TRILNT(RS,THK,PR,SLN,RN,SN,BL,BLR,BLS)
	IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION RS(2,3)
	DIMENSION SLN(3),RN(3),SN(3),BL(3)

C	------------------------------------------
C	COMPUTE FOR ELEMENT BOUNDARY LENGTHS - SLN
C	------------------------------------------
	SLN(1)=DSQRT((RS(1,2)-RS(1,1))**2.0D0+(RS(2,2)-RS(2,1))**2.0D0)
	SLN(2)=DSQRT((RS(1,3)-RS(1,2))**2.0D0+(RS(2,3)-RS(2,2))**2.0D0)
	SLN(3)=DSQRT((RS(1,1)-RS(1,3))**2.0D0+(RS(2,1)-RS(2,3))**2.0D0)
C	-----------------------------------------
C	COMPUTE FOR ELEMENT BOUNDARY NORMALS - RN
C	-----------------------------------------
	RN(1)=(RS(2,2)-RS(2,1))/SLN(1)
	RN(2)=(RS(2,3)-RS(2,2))/SLN(2)
	RN(3)=(RS(2,1)-RS(2,3))/SLN(3)
C	------------------------------------------
C	COMPUTE FOR ELEMENT BOUNDARY TANGENTS - SN
C	------------------------------------------
	SN(1)=(RS(1,1)-RS(1,2))/SLN(1)
	SN(2)=(RS(1,2)-RS(1,3))/SLN(2)
	SN(3)=(RS(1,3)-RS(1,1))/SLN(3)
C	-------------------------------------------------------
C	COMPUTE FOR ELEMENT SHEAR DEFORMATION PARAMETER, LAMBDA
C	-------------------------------------------------------
	FACT=(12.0D0*THK*THK)/(5.0D0*(1.0D0-PR))

C	PARAMETERS ALONG THE BOUNDARY
	BL(1)=1.0D0/(1.0D0+FACT/(SLN(1)*SLN(1)))
	BL(2)=1.0D0/(1.0D0+FACT/(SLN(2)*SLN(2)))
	BL(3)=1.0D0/(1.0D0+FACT/(SLN(3)*SLN(3)))

C	EFFECTIVE LENGTHS ALONG S AND R
	EFR = (SLN(1)+SLN(2))/2.0D0
	EFS = (SLN(2)+SLN(3))/2.0D0

C	ELEMENT PARAMETERS
	BLR=1.0D0/(1.0D0+FACT/(EFR*EFR))
	BLS=1.0D0/(1.0D0+FACT/(EFS*EFS))

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRICF(RS,BLR,BLS,SLN,RN,SN,FLR,FLS,DWR,DWS,ANT)
	IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION RS(2,3),FLR(3),FLS(3)
	DIMENSION SLN(3),SN(3),RN(3)
	DIMENSION DWR(9),DWS(9),ANT(3)

C	-------------------
C	AREA OF THE ELEMENT
C	-------------------
	AREA = 0.50*( (RS(1,1)*RS(2,2)+RS(1,2)*RS(2,3)+RS(1,3)*RS(2,1))-
     #              (RS(2,1)*RS(1,2)+RS(2,2)*RS(1,3)+RS(2,3)*RS(1,1))  )

C	------------------------------------------------
C	COMPUTE FOR COMMON FACTORS FOR AREA INTEGRATIONS
C	------------------------------------------------
	F1 = 1.0D0-BLR
	F2 = BLS-1.0D0

C	------------------------------
C	AREA INTEGRATION FACTORS - ANT
C	------------------------------
	ANT(1) = AREA/3.0D0
	ANT(2) = AREA/3.0D0
	ANT(3) = AREA/3.0D0

C	----------------------------
C	FACTORS FOR ROTATION ABOUT S
C	----------------------------
	FLS(1) = F1*ANT(1)
	FLS(2) = F1*ANT(2)
	FLS(3) = F1*ANT(3)

C	----------------------------
C	FACTORS FOR ROTATION ABOUT R
C	----------------------------
	FLR(1) = F2*ANT(1)
	FLR(2) = F2*ANT(2)
	FLR(3) = F2*ANT(3)

C	---------------------
C	FACTORS FOR W WITH RN
C	---------------------
	DWR(1) = (RN(3)**2*SLN(3)**2-RN(1)**2*SLN(1)**2)*BLR/12
	DWR(2) = (RN(1)**2*SLN(1)**2-RN(2)**2*SLN(2)**2)*BLR/12
	DWR(3) = (RN(2)**2*SLN(2)**2-RN(3)**2*SLN(3)**2)*BLR/12

	DWR(4) = (RN(3)*SLN(3)**2*SN(3)-RN(1)*SLN(1)**2*SN(1))*BLR/12
	DWR(5) = (RN(1)*SLN(1)**2*SN(1)-RN(2)*SLN(2)**2*SN(2))*BLR/12
	DWR(6) = (RN(2)*SLN(2)**2*SN(2)-RN(3)*SLN(3)**2*SN(3))*BLR/12

	DWR(7) = (-SLN(1)*RN(1)-SLN(3)*RN(3))*BLR/2
	DWR(8) = (-SLN(1)*RN(1)-SLN(2)*RN(2))*BLR/2
	DWR(9) = (-SLN(3)*RN(3)-SLN(2)*RN(2))*BLR/2

C	---------------------
C	FACTORS FOR W WITH NS
C	---------------------
	DWS(1) = (SN(3)*SLN(3)**2*RN(3)-SN(1)*SLN(1)**2*RN(1))*BLS/12
	DWS(2) = (SN(1)*SLN(1)**2*RN(1)-SN(2)*SLN(2)**2*RN(2))*BLS/12
	DWS(3) = (SN(2)*SLN(2)**2*RN(2)-SN(3)*SLN(3)**2*RN(3))*BLS/12

	DWS(4) = (SN(3)**2*SLN(3)**2-SN(1)**2*SLN(1)**2)*BLS/12
	DWS(5) = (SN(1)**2*SLN(1)**2-SN(2)**2*SLN(2)**2)*BLS/12
	DWS(6) = (SN(2)**2*SLN(2)**2-SN(3)**2*SLN(3)**2)*BLS/12

	DWS(7) = (-SLN(1)*SN(1)-SLN(3)*SN(3))*BLS/2
	DWS(8) = (-SLN(1)*SN(1)-SLN(2)*SN(2))*BLS/2
	DWS(9) = (-SLN(3)*SN(3)-SLN(2)*SN(2))*BLS/2

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRITRE(TT,REDIS,TEDIS)
	IMPLICIT REAL*8 (A-H,O-Z)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTE FOR COROTATIONAL DISPLACEMENTS IN LOCAL
C				COORDINATES
C
C	INPUT VARIABLES
C	TT(3,3) = TRANSPOSE OF TRANSFORMATION MATRIX
C	REDIS(18) = VECTOR OF COROTATIONAL DISPLACEMENTS IN GLOBAL SPACE
C
C	OUTPUT VARIABLES
C	TEDIS(18) = TT*REDIS
C     ----------------------------------------------------------------
	DIMENSION TT(3,3),REDIS(18)
	DIMENSION TEDIS(18)

      TEDIS(1) = TT(1,1)*REDIS(1)+TT(1,2)*REDIS(2)+TT(1,3)*REDIS(3)
      TEDIS(2) = TT(2,1)*REDIS(1)+TT(2,2)*REDIS(2)+TT(2,3)*REDIS(3)
      TEDIS(3) = TT(3,1)*REDIS(1)+TT(3,2)*REDIS(2)+TT(3,3)*REDIS(3)
      TEDIS(4) = TT(1,1)*REDIS(4)+TT(1,2)*REDIS(5)+TT(1,3)*REDIS(6)
      TEDIS(5) = TT(2,1)*REDIS(4)+TT(2,2)*REDIS(5)+TT(2,3)*REDIS(6)
      TEDIS(6) = TT(3,1)*REDIS(4)+TT(3,2)*REDIS(5)+TT(3,3)*REDIS(6)
      TEDIS(7) = TT(1,1)*REDIS(7)+TT(1,2)*REDIS(8)+TT(1,3)*REDIS(9)
      TEDIS(8) = TT(2,1)*REDIS(7)+TT(2,2)*REDIS(8)+TT(2,3)*REDIS(9)
      TEDIS(9) = TT(3,1)*REDIS(7)+TT(3,2)*REDIS(8)+TT(3,3)*REDIS(9)
      TEDIS(10) = TT(1,1)*REDIS(10)+TT(1,2)*REDIS(11)+TT(1,3)*REDIS(12)
      TEDIS(11) = TT(2,1)*REDIS(10)+TT(2,2)*REDIS(11)+TT(2,3)*REDIS(12)
      TEDIS(12) = TT(3,1)*REDIS(10)+TT(3,2)*REDIS(11)+TT(3,3)*REDIS(12)
      TEDIS(13) = TT(1,1)*REDIS(13)+TT(1,2)*REDIS(14)+TT(1,3)*REDIS(15)
      TEDIS(14) = TT(2,1)*REDIS(13)+TT(2,2)*REDIS(14)+TT(2,3)*REDIS(15)
      TEDIS(15) = TT(3,1)*REDIS(13)+TT(3,2)*REDIS(14)+TT(3,3)*REDIS(15)
      TEDIS(16) = TT(1,1)*REDIS(16)+TT(1,2)*REDIS(17)+TT(1,3)*REDIS(18)
      TEDIS(17) = TT(2,1)*REDIS(16)+TT(2,2)*REDIS(17)+TT(2,3)*REDIS(18)
      TEDIS(18) = TT(3,1)*REDIS(16)+TT(3,2)*REDIS(17)+TT(3,3)*REDIS(18)

      RETURN
      END
C
C=====================================================================
      SUBROUTINE KGLB18(ST,TT,S)

	IMPLICIT REAL*8 (A-H,O-Z)
C	------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR THE GLOBAL STIFFNESS MATRIX T*KL*TT
C
C	INPUT VARIABLES
C	ST(171) = ELEMENT STIFFNESS MATRIX IN LOCAL COORDINATES
C	TT(3,3) = TRANSPOSE OF THE TRANSFORMATION MATRIX
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIAVBLES
C	S(*)		= GLOBAL STIFFNESS MATRIX IN ROW FORMAT
C	------------------------------------------------------------
	DIMENSION ST(171),TT(3,3)
	DIMENSION S(171)

      t4 = TT(1,1)*ST(1)+TT(2,1)*ST(2)+TT(3,1)*ST(3)
      t9 = TT(1,1)*ST(2)+TT(2,1)*ST(19)+TT(3,1)*ST(20)
      t14 = TT(1,1)*ST(3)+TT(2,1)*ST(20)+TT(3,1)*ST(36)
      t28 = TT(1,1)*ST(4)+TT(2,1)*ST(21)+TT(3,1)*ST(37)
      t33 = TT(1,1)*ST(5)+TT(2,1)*ST(22)+TT(3,1)*ST(38)
      t38 = TT(1,1)*ST(6)+TT(2,1)*ST(23)+TT(3,1)*ST(39)
      t52 = TT(1,1)*ST(7)+TT(2,1)*ST(24)+TT(3,1)*ST(40)
      t57 = TT(1,1)*ST(8)+TT(2,1)*ST(25)+TT(3,1)*ST(41)
      t62 = TT(1,1)*ST(9)+TT(2,1)*ST(26)+TT(3,1)*ST(42)
      t76 = TT(1,1)*ST(10)+TT(2,1)*ST(27)+TT(3,1)*ST(43)
      t81 = TT(1,1)*ST(11)+TT(2,1)*ST(28)+TT(3,1)*ST(44)
      t86 = TT(1,1)*ST(12)+TT(2,1)*ST(29)+TT(3,1)*ST(45)
      t100 = TT(1,1)*ST(13)+TT(2,1)*ST(30)+TT(3,1)*ST(46)
      t105 = TT(1,1)*ST(14)+TT(2,1)*ST(31)+TT(3,1)*ST(47)
      t110 = TT(1,1)*ST(15)+TT(2,1)*ST(32)+TT(3,1)*ST(48)
      t124 = TT(1,1)*ST(16)+TT(2,1)*ST(33)+TT(3,1)*ST(49)
      t129 = TT(1,1)*ST(17)+TT(2,1)*ST(34)+TT(3,1)*ST(50)
      t134 = TT(1,1)*ST(18)+TT(2,1)*ST(35)+TT(3,1)*ST(51)
      t148 = TT(1,2)*ST(1)+TT(2,2)*ST(2)+TT(3,2)*ST(3)
      t153 = TT(1,2)*ST(2)+TT(2,2)*ST(19)+TT(3,2)*ST(20)
      t158 = TT(1,2)*ST(3)+TT(2,2)*ST(20)+TT(3,2)*ST(36)
      t168 = TT(1,2)*ST(4)+TT(2,2)*ST(21)+TT(3,2)*ST(37)
      t173 = TT(1,2)*ST(5)+TT(2,2)*ST(22)+TT(3,2)*ST(38)
      t178 = TT(1,2)*ST(6)+TT(2,2)*ST(23)+TT(3,2)*ST(39)
      t192 = TT(1,2)*ST(7)+TT(2,2)*ST(24)+TT(3,2)*ST(40)
      t197 = TT(1,2)*ST(8)+TT(2,2)*ST(25)+TT(3,2)*ST(41)
      t202 = TT(1,2)*ST(9)+TT(2,2)*ST(26)+TT(3,2)*ST(42)
      t216 = TT(1,2)*ST(10)+TT(2,2)*ST(27)+TT(3,2)*ST(43)
      t221 = TT(1,2)*ST(11)+TT(2,2)*ST(28)+TT(3,2)*ST(44)
      t226 = TT(1,2)*ST(12)+TT(2,2)*ST(29)+TT(3,2)*ST(45)
      t240 = TT(1,2)*ST(13)+TT(2,2)*ST(30)+TT(3,2)*ST(46)
      t245 = TT(1,2)*ST(14)+TT(2,2)*ST(31)+TT(3,2)*ST(47)
      t250 = TT(1,2)*ST(15)+TT(2,2)*ST(32)+TT(3,2)*ST(48)
      t264 = TT(1,2)*ST(16)+TT(2,2)*ST(33)+TT(3,2)*ST(49)
      t269 = TT(1,2)*ST(17)+TT(2,2)*ST(34)+TT(3,2)*ST(50)
      t274 = TT(1,2)*ST(18)+TT(2,2)*ST(35)+TT(3,2)*ST(51)
      t304 = TT(1,3)*ST(4)+TT(2,3)*ST(21)+TT(3,3)*ST(37)
      t309 = TT(1,3)*ST(5)+TT(2,3)*ST(22)+TT(3,3)*ST(38)
      t314 = TT(1,3)*ST(6)+TT(2,3)*ST(23)+TT(3,3)*ST(39)
      t328 = TT(1,3)*ST(7)+TT(2,3)*ST(24)+TT(3,3)*ST(40)
      t333 = TT(1,3)*ST(8)+TT(2,3)*ST(25)+TT(3,3)*ST(41)
      t338 = TT(1,3)*ST(9)+TT(2,3)*ST(26)+TT(3,3)*ST(42)
      t352 = TT(1,3)*ST(10)+TT(2,3)*ST(27)+TT(3,3)*ST(43)
      t357 = TT(1,3)*ST(11)+TT(2,3)*ST(28)+TT(3,3)*ST(44)
      t362 = TT(1,3)*ST(12)+TT(2,3)*ST(29)+TT(3,3)*ST(45)
      t376 = TT(1,3)*ST(13)+TT(2,3)*ST(30)+TT(3,3)*ST(46)
      t381 = TT(1,3)*ST(14)+TT(2,3)*ST(31)+TT(3,3)*ST(47)
      t386 = TT(1,3)*ST(15)+TT(2,3)*ST(32)+TT(3,3)*ST(48)
      t400 = TT(1,3)*ST(16)+TT(2,3)*ST(33)+TT(3,3)*ST(49)
      t405 = TT(1,3)*ST(17)+TT(2,3)*ST(34)+TT(3,3)*ST(50)
      t410 = TT(1,3)*ST(18)+TT(2,3)*ST(35)+TT(3,3)*ST(51)
      t424 = TT(1,1)*ST(52)+TT(2,1)*ST(53)+TT(3,1)*ST(54)
      t429 = TT(1,1)*ST(53)+TT(2,1)*ST(67)+TT(3,1)*ST(68)
      t434 = TT(1,1)*ST(54)+TT(2,1)*ST(68)+TT(3,1)*ST(81)
      t448 = TT(1,1)*ST(55)+TT(2,1)*ST(69)+TT(3,1)*ST(82)
      t453 = TT(1,1)*ST(56)+TT(2,1)*ST(70)+TT(3,1)*ST(83)
      t458 = TT(1,1)*ST(57)+TT(2,1)*ST(71)+TT(3,1)*ST(84)
      t472 = TT(1,1)*ST(58)+TT(2,1)*ST(72)+TT(3,1)*ST(85)
      t477 = TT(1,1)*ST(59)+TT(2,1)*ST(73)+TT(3,1)*ST(86)
      t482 = TT(1,1)*ST(60)+TT(2,1)*ST(74)+TT(3,1)*ST(87)
      t496 = TT(1,1)*ST(61)+TT(2,1)*ST(75)+TT(3,1)*ST(88)
      t501 = TT(1,1)*ST(62)+TT(2,1)*ST(76)+TT(3,1)*ST(89)
      t506 = TT(1,1)*ST(63)+TT(2,1)*ST(77)+TT(3,1)*ST(90)
      t520 = TT(1,1)*ST(64)+TT(2,1)*ST(78)+TT(3,1)*ST(91)
      t525 = TT(1,1)*ST(65)+TT(2,1)*ST(79)+TT(3,1)*ST(92)
      t530 = TT(1,1)*ST(66)+TT(2,1)*ST(80)+TT(3,1)*ST(93)
      t544 = TT(1,2)*ST(52)+TT(2,2)*ST(53)+TT(3,2)*ST(54)
      t549 = TT(1,2)*ST(53)+TT(2,2)*ST(67)+TT(3,2)*ST(68)
      t554 = TT(1,2)*ST(54)+TT(2,2)*ST(68)+TT(3,2)*ST(81)
      t564 = TT(1,2)*ST(55)+TT(2,2)*ST(69)+TT(3,2)*ST(82)
      t569 = TT(1,2)*ST(56)+TT(2,2)*ST(70)+TT(3,2)*ST(83)
      t574 = TT(1,2)*ST(57)+TT(2,2)*ST(71)+TT(3,2)*ST(84)
      t588 = TT(1,2)*ST(58)+TT(2,2)*ST(72)+TT(3,2)*ST(85)
      t593 = TT(1,2)*ST(59)+TT(2,2)*ST(73)+TT(3,2)*ST(86)
      t598 = TT(1,2)*ST(60)+TT(2,2)*ST(74)+TT(3,2)*ST(87)
      t612 = TT(1,2)*ST(61)+TT(2,2)*ST(75)+TT(3,2)*ST(88)
      t617 = TT(1,2)*ST(62)+TT(2,2)*ST(76)+TT(3,2)*ST(89)
      t622 = TT(1,2)*ST(63)+TT(2,2)*ST(77)+TT(3,2)*ST(90)
      t636 = TT(1,2)*ST(64)+TT(2,2)*ST(78)+TT(3,2)*ST(91)
      t641 = TT(1,2)*ST(65)+TT(2,2)*ST(79)+TT(3,2)*ST(92)
      t646 = TT(1,2)*ST(66)+TT(2,2)*ST(80)+TT(3,2)*ST(93)
      t676 = TT(1,3)*ST(55)+TT(2,3)*ST(69)+TT(3,3)*ST(82)
      t681 = TT(1,3)*ST(56)+TT(2,3)*ST(70)+TT(3,3)*ST(83)
      t686 = TT(1,3)*ST(57)+TT(2,3)*ST(71)+TT(3,3)*ST(84)
      t700 = TT(1,3)*ST(58)+TT(2,3)*ST(72)+TT(3,3)*ST(85)
      t705 = TT(1,3)*ST(59)+TT(2,3)*ST(73)+TT(3,3)*ST(86)
      t710 = TT(1,3)*ST(60)+TT(2,3)*ST(74)+TT(3,3)*ST(87)
      t724 = TT(1,3)*ST(61)+TT(2,3)*ST(75)+TT(3,3)*ST(88)
      t729 = TT(1,3)*ST(62)+TT(2,3)*ST(76)+TT(3,3)*ST(89)
      t734 = TT(1,3)*ST(63)+TT(2,3)*ST(77)+TT(3,3)*ST(90)
      t748 = TT(1,3)*ST(64)+TT(2,3)*ST(78)+TT(3,3)*ST(91)
      t753 = TT(1,3)*ST(65)+TT(2,3)*ST(79)+TT(3,3)*ST(92)
      t758 = TT(1,3)*ST(66)+TT(2,3)*ST(80)+TT(3,3)*ST(93)
      t772 = TT(1,1)*ST(94)+TT(2,1)*ST(95)+TT(3,1)*ST(96)
      t777 = TT(1,1)*ST(95)+TT(2,1)*ST(106)+TT(3,1)*ST(107)
      t782 = TT(1,1)*ST(96)+TT(2,1)*ST(107)+TT(3,1)*ST(117)
      t796 = TT(1,1)*ST(97)+TT(2,1)*ST(108)+TT(3,1)*ST(118)
      t801 = TT(1,1)*ST(98)+TT(2,1)*ST(109)+TT(3,1)*ST(119)
      t806 = TT(1,1)*ST(99)+TT(2,1)*ST(110)+TT(3,1)*ST(120)
      t820 = TT(1,1)*ST(100)+TT(2,1)*ST(111)+TT(3,1)*ST(121)
      t825 = TT(1,1)*ST(101)+TT(2,1)*ST(112)+TT(3,1)*ST(122)
      t830 = TT(1,1)*ST(102)+TT(2,1)*ST(113)+TT(3,1)*ST(123)
      t844 = TT(1,1)*ST(103)+TT(2,1)*ST(114)+TT(3,1)*ST(124)
      t849 = TT(1,1)*ST(104)+TT(2,1)*ST(115)+TT(3,1)*ST(125)
      t854 = TT(1,1)*ST(105)+TT(2,1)*ST(116)+TT(3,1)*ST(126)
      t868 = TT(1,2)*ST(94)+TT(2,2)*ST(95)+TT(3,2)*ST(96)
      t873 = TT(1,2)*ST(95)+TT(2,2)*ST(106)+TT(3,2)*ST(107)
      t878 = TT(1,2)*ST(96)+TT(2,2)*ST(107)+TT(3,2)*ST(117)
      t888 = TT(1,2)*ST(97)+TT(2,2)*ST(108)+TT(3,2)*ST(118)
      t893 = TT(1,2)*ST(98)+TT(2,2)*ST(109)+TT(3,2)*ST(119)
      t898 = TT(1,2)*ST(99)+TT(2,2)*ST(110)+TT(3,2)*ST(120)
      t912 = TT(1,2)*ST(100)+TT(2,2)*ST(111)+TT(3,2)*ST(121)
      t917 = TT(1,2)*ST(101)+TT(2,2)*ST(112)+TT(3,2)*ST(122)
      t922 = TT(1,2)*ST(102)+TT(2,2)*ST(113)+TT(3,2)*ST(123)
      t936 = TT(1,2)*ST(103)+TT(2,2)*ST(114)+TT(3,2)*ST(124)
      t941 = TT(1,2)*ST(104)+TT(2,2)*ST(115)+TT(3,2)*ST(125)
      t946 = TT(1,2)*ST(105)+TT(2,2)*ST(116)+TT(3,2)*ST(126)
      t976 = TT(1,3)*ST(97)+TT(2,3)*ST(108)+TT(3,3)*ST(118)
      t981 = TT(1,3)*ST(98)+TT(2,3)*ST(109)+TT(3,3)*ST(119)
      t986 = TT(1,3)*ST(99)+TT(2,3)*ST(110)+TT(3,3)*ST(120)
      t1000 = TT(1,3)*ST(100)+TT(2,3)*ST(111)+TT(3,3)*ST(121)
      t1005 = TT(1,3)*ST(101)+TT(2,3)*ST(112)+TT(3,3)*ST(122)
      t1010 = TT(1,3)*ST(102)+TT(2,3)*ST(113)+TT(3,3)*ST(123)
      t1024 = TT(1,3)*ST(103)+TT(2,3)*ST(114)+TT(3,3)*ST(124)
      t1029 = TT(1,3)*ST(104)+TT(2,3)*ST(115)+TT(3,3)*ST(125)
      t1034 = TT(1,3)*ST(105)+TT(2,3)*ST(116)+TT(3,3)*ST(126)
      t1048 = TT(1,1)*ST(127)+TT(2,1)*ST(128)+TT(3,1)*ST(129)
      t1053 = TT(1,1)*ST(128)+TT(2,1)*ST(136)+TT(3,1)*ST(137)
      t1058 = TT(1,1)*ST(129)+TT(2,1)*ST(137)+TT(3,1)*ST(144)
      t1072 = TT(1,1)*ST(130)+TT(2,1)*ST(138)+TT(3,1)*ST(145)
      t1077 = TT(1,1)*ST(131)+TT(2,1)*ST(139)+TT(3,1)*ST(146)
      t1082 = TT(1,1)*ST(132)+TT(2,1)*ST(140)+TT(3,1)*ST(147)
      t1096 = TT(1,1)*ST(133)+TT(2,1)*ST(141)+TT(3,1)*ST(148)
      t1101 = TT(1,1)*ST(134)+TT(2,1)*ST(142)+TT(3,1)*ST(149)
      t1106 = TT(1,1)*ST(135)+TT(2,1)*ST(143)+TT(3,1)*ST(150)
      t1120 = TT(1,2)*ST(127)+TT(2,2)*ST(128)+TT(3,2)*ST(129)
      t1125 = TT(1,2)*ST(128)+TT(2,2)*ST(136)+TT(3,2)*ST(137)
      t1130 = TT(1,2)*ST(129)+TT(2,2)*ST(137)+TT(3,2)*ST(144)
      t1140 = TT(1,2)*ST(130)+TT(2,2)*ST(138)+TT(3,2)*ST(145)
      t1145 = TT(1,2)*ST(131)+TT(2,2)*ST(139)+TT(3,2)*ST(146)
      t1150 = TT(1,2)*ST(132)+TT(2,2)*ST(140)+TT(3,2)*ST(147)
      t1164 = TT(1,2)*ST(133)+TT(2,2)*ST(141)+TT(3,2)*ST(148)
      t1169 = TT(1,2)*ST(134)+TT(2,2)*ST(142)+TT(3,2)*ST(149)
      t1174 = TT(1,2)*ST(135)+TT(2,2)*ST(143)+TT(3,2)*ST(150)
      t1204 = TT(1,3)*ST(130)+TT(2,3)*ST(138)+TT(3,3)*ST(145)
      t1209 = TT(1,3)*ST(131)+TT(2,3)*ST(139)+TT(3,3)*ST(146)
      t1214 = TT(1,3)*ST(132)+TT(2,3)*ST(140)+TT(3,3)*ST(147)
      t1228 = TT(1,3)*ST(133)+TT(2,3)*ST(141)+TT(3,3)*ST(148)
      t1233 = TT(1,3)*ST(134)+TT(2,3)*ST(142)+TT(3,3)*ST(149)
      t1238 = TT(1,3)*ST(135)+TT(2,3)*ST(143)+TT(3,3)*ST(150)
      t1252 = TT(1,1)*ST(151)+TT(2,1)*ST(152)+TT(3,1)*ST(153)
      t1257 = TT(1,1)*ST(152)+TT(2,1)*ST(157)+TT(3,1)*ST(158)
      t1262 = TT(1,1)*ST(153)+TT(2,1)*ST(158)+TT(3,1)*ST(162)
      t1276 = TT(1,1)*ST(154)+TT(2,1)*ST(159)+TT(3,1)*ST(163)
      t1281 = TT(1,1)*ST(155)+TT(2,1)*ST(160)+TT(3,1)*ST(164)
      t1286 = TT(1,1)*ST(156)+TT(2,1)*ST(161)+TT(3,1)*ST(165)
      t1300 = TT(1,2)*ST(151)+TT(2,2)*ST(152)+TT(3,2)*ST(153)
      t1305 = TT(1,2)*ST(152)+TT(2,2)*ST(157)+TT(3,2)*ST(158)
      t1310 = TT(1,2)*ST(153)+TT(2,2)*ST(158)+TT(3,2)*ST(162)
      t1320 = TT(1,2)*ST(154)+TT(2,2)*ST(159)+TT(3,2)*ST(163)
      t1325 = TT(1,2)*ST(155)+TT(2,2)*ST(160)+TT(3,2)*ST(164)
      t1330 = TT(1,2)*ST(156)+TT(2,2)*ST(161)+TT(3,2)*ST(165)
      t1360 = TT(1,3)*ST(154)+TT(2,3)*ST(159)+TT(3,3)*ST(163)
      t1365 = TT(1,3)*ST(155)+TT(2,3)*ST(160)+TT(3,3)*ST(164)
      t1370 = TT(1,3)*ST(156)+TT(2,3)*ST(161)+TT(3,3)*ST(165)
      t1384 = TT(1,1)*ST(166)+TT(2,1)*ST(167)+TT(3,1)*ST(168)
      t1389 = TT(1,1)*ST(167)+TT(2,1)*ST(169)+TT(3,1)*ST(170)
      t1394 = TT(1,1)*ST(168)+TT(2,1)*ST(170)+TT(3,1)*ST(171)
      t1408 = TT(1,2)*ST(166)+TT(2,2)*ST(167)+TT(3,2)*ST(168)
      t1413 = TT(1,2)*ST(167)+TT(2,2)*ST(169)+TT(3,2)*ST(170)
      t1418 = TT(1,2)*ST(168)+TT(2,2)*ST(170)+TT(3,2)*ST(171)
      S(1) = t4*TT(1,1)+t9*TT(2,1)+t14*TT(3,1)
      S(2) = t4*TT(1,2)+t9*TT(2,2)+t14*TT(3,2)
      S(3) = t4*TT(1,3)+t9*TT(2,3)+t14*TT(3,3)
      S(4) = t28*TT(1,1)+t33*TT(2,1)+t38*TT(3,1)
      S(5) = t28*TT(1,2)+t33*TT(2,2)+t38*TT(3,2)
      S(6) = t28*TT(1,3)+t33*TT(2,3)+t38*TT(3,3)
      S(7) = t52*TT(1,1)+t57*TT(2,1)+t62*TT(3,1)
      S(8) = t52*TT(1,2)+t57*TT(2,2)+t62*TT(3,2)
      S(9) = t52*TT(1,3)+t57*TT(2,3)+t62*TT(3,3)
      S(10) = t76*TT(1,1)+t81*TT(2,1)+t86*TT(3,1)
      S(11) = t76*TT(1,2)+t81*TT(2,2)+t86*TT(3,2)
      S(12) = t76*TT(1,3)+t81*TT(2,3)+t86*TT(3,3)
      S(13) = t100*TT(1,1)+t105*TT(2,1)+t110*TT(3,1)
      S(14) = t100*TT(1,2)+t105*TT(2,2)+t110*TT(3,2)
      S(15) = t100*TT(1,3)+t105*TT(2,3)+t110*TT(3,3)
      S(16) = t124*TT(1,1)+t129*TT(2,1)+t134*TT(3,1)
      S(17) = t124*TT(1,2)+t129*TT(2,2)+t134*TT(3,2)
      S(18) = t124*TT(1,3)+t129*TT(2,3)+t134*TT(3,3)
      S(19) = t148*TT(1,2)+t153*TT(2,2)+t158*TT(3,2)
      S(20) = t148*TT(1,3)+t153*TT(2,3)+t158*TT(3,3)
      S(21) = t168*TT(1,1)+t173*TT(2,1)+t178*TT(3,1)
      S(22) = t168*TT(1,2)+t173*TT(2,2)+t178*TT(3,2)
      S(23) = t168*TT(1,3)+t173*TT(2,3)+t178*TT(3,3)
      S(24) = t192*TT(1,1)+t197*TT(2,1)+t202*TT(3,1)
      S(25) = t192*TT(1,2)+t197*TT(2,2)+t202*TT(3,2)
      S(26) = t192*TT(1,3)+t197*TT(2,3)+t202*TT(3,3)
      S(27) = t216*TT(1,1)+t221*TT(2,1)+t226*TT(3,1)
      S(28) = t216*TT(1,2)+t221*TT(2,2)+t226*TT(3,2)
      S(29) = t216*TT(1,3)+t221*TT(2,3)+t226*TT(3,3)
      S(30) = t240*TT(1,1)+t245*TT(2,1)+t250*TT(3,1)
      S(31) = t240*TT(1,2)+t245*TT(2,2)+t250*TT(3,2)
      S(32) = t240*TT(1,3)+t245*TT(2,3)+t250*TT(3,3)
      S(33) = t264*TT(1,1)+t269*TT(2,1)+t274*TT(3,1)
      S(34) = t264*TT(1,2)+t269*TT(2,2)+t274*TT(3,2)
      S(35) = t264*TT(1,3)+t269*TT(2,3)+t274*TT(3,3)
      S(36) = (TT(1,3)*ST(1)+TT(2,3)*ST(2)+TT(3,3)*ST(3))*TT(1,3)+(TT(1,
     #3)*ST(2)+TT(2,3)*ST(19)+TT(3,3)*ST(20))*TT(2,3)+(TT(1,3)*ST(3)+TT(
     #2,3)*ST(20)+TT(3,3)*ST(36))*TT(3,3)
      S(37) = t304*TT(1,1)+t309*TT(2,1)+t314*TT(3,1)
      S(38) = t304*TT(1,2)+t309*TT(2,2)+t314*TT(3,2)
      S(39) = t304*TT(1,3)+t309*TT(2,3)+t314*TT(3,3)
      S(40) = t328*TT(1,1)+t333*TT(2,1)+t338*TT(3,1)
      S(41) = t328*TT(1,2)+t333*TT(2,2)+t338*TT(3,2)
      S(42) = t328*TT(1,3)+t333*TT(2,3)+t338*TT(3,3)
      S(43) = t352*TT(1,1)+t357*TT(2,1)+t362*TT(3,1)
      S(44) = t352*TT(1,2)+t357*TT(2,2)+t362*TT(3,2)
      S(45) = t352*TT(1,3)+t357*TT(2,3)+t362*TT(3,3)
      S(46) = t376*TT(1,1)+t381*TT(2,1)+t386*TT(3,1)
      S(47) = t376*TT(1,2)+t381*TT(2,2)+t386*TT(3,2)
      S(48) = t376*TT(1,3)+t381*TT(2,3)+t386*TT(3,3)
      S(49) = t400*TT(1,1)+t405*TT(2,1)+t410*TT(3,1)
      S(50) = t400*TT(1,2)+t405*TT(2,2)+t410*TT(3,2)
      S(51) = t400*TT(1,3)+t405*TT(2,3)+t410*TT(3,3)
      S(52) = t424*TT(1,1)+t429*TT(2,1)+t434*TT(3,1)
      S(53) = t424*TT(1,2)+t429*TT(2,2)+t434*TT(3,2)
      S(54) = t424*TT(1,3)+t429*TT(2,3)+t434*TT(3,3)
      S(55) = t448*TT(1,1)+t453*TT(2,1)+t458*TT(3,1)
      S(56) = t448*TT(1,2)+t453*TT(2,2)+t458*TT(3,2)
      S(57) = t448*TT(1,3)+t453*TT(2,3)+t458*TT(3,3)
      S(58) = t472*TT(1,1)+t477*TT(2,1)+t482*TT(3,1)
      S(59) = t472*TT(1,2)+t477*TT(2,2)+t482*TT(3,2)
      S(60) = t472*TT(1,3)+t477*TT(2,3)+t482*TT(3,3)
      S(61) = t496*TT(1,1)+t501*TT(2,1)+t506*TT(3,1)
      S(62) = t496*TT(1,2)+t501*TT(2,2)+t506*TT(3,2)
      S(63) = t496*TT(1,3)+t501*TT(2,3)+t506*TT(3,3)
      S(64) = t520*TT(1,1)+t525*TT(2,1)+t530*TT(3,1)
      S(65) = t520*TT(1,2)+t525*TT(2,2)+t530*TT(3,2)
      S(66) = t520*TT(1,3)+t525*TT(2,3)+t530*TT(3,3)
      S(67) = t544*TT(1,2)+t549*TT(2,2)+t554*TT(3,2)
      S(68) = t544*TT(1,3)+t549*TT(2,3)+t554*TT(3,3)
      S(69) = t564*TT(1,1)+t569*TT(2,1)+t574*TT(3,1)
      S(70) = t564*TT(1,2)+t569*TT(2,2)+t574*TT(3,2)
      S(71) = t564*TT(1,3)+t569*TT(2,3)+t574*TT(3,3)
      S(72) = t588*TT(1,1)+t593*TT(2,1)+t598*TT(3,1)
      S(73) = t588*TT(1,2)+t593*TT(2,2)+t598*TT(3,2)
      S(74) = t588*TT(1,3)+t593*TT(2,3)+t598*TT(3,3)
      S(75) = t612*TT(1,1)+t617*TT(2,1)+t622*TT(3,1)
      S(76) = t612*TT(1,2)+t617*TT(2,2)+t622*TT(3,2)
      S(77) = t612*TT(1,3)+t617*TT(2,3)+t622*TT(3,3)
      S(78) = t636*TT(1,1)+t641*TT(2,1)+t646*TT(3,1)
      S(79) = t636*TT(1,2)+t641*TT(2,2)+t646*TT(3,2)
      S(80) = t636*TT(1,3)+t641*TT(2,3)+t646*TT(3,3)
      S(81) = (TT(1,3)*ST(52)+TT(2,3)*ST(53)+TT(3,3)*ST(54))*TT(1,3)+(TT
     #(1,3)*ST(53)+TT(2,3)*ST(67)+TT(3,3)*ST(68))*TT(2,3)+(TT(1,3)*ST(54
     #)+TT(2,3)*ST(68)+TT(3,3)*ST(81))*TT(3,3)
      S(82) = t676*TT(1,1)+t681*TT(2,1)+t686*TT(3,1)
      S(83) = t676*TT(1,2)+t681*TT(2,2)+t686*TT(3,2)
      S(84) = t676*TT(1,3)+t681*TT(2,3)+t686*TT(3,3)
      S(85) = t700*TT(1,1)+t705*TT(2,1)+t710*TT(3,1)
      S(86) = t700*TT(1,2)+t705*TT(2,2)+t710*TT(3,2)
      S(87) = t700*TT(1,3)+t705*TT(2,3)+t710*TT(3,3)
      S(88) = t724*TT(1,1)+t729*TT(2,1)+t734*TT(3,1)
      S(89) = t724*TT(1,2)+t729*TT(2,2)+t734*TT(3,2)
      S(90) = t724*TT(1,3)+t729*TT(2,3)+t734*TT(3,3)
      S(91) = t748*TT(1,1)+t753*TT(2,1)+t758*TT(3,1)
      S(92) = t748*TT(1,2)+t753*TT(2,2)+t758*TT(3,2)
      S(93) = t748*TT(1,3)+t753*TT(2,3)+t758*TT(3,3)
      S(94) = t772*TT(1,1)+t777*TT(2,1)+t782*TT(3,1)
      S(95) = t772*TT(1,2)+t777*TT(2,2)+t782*TT(3,2)
      S(96) = t772*TT(1,3)+t777*TT(2,3)+t782*TT(3,3)
      S(97) = t796*TT(1,1)+t801*TT(2,1)+t806*TT(3,1)
      S(98) = t796*TT(1,2)+t801*TT(2,2)+t806*TT(3,2)
      S(99) = t796*TT(1,3)+t801*TT(2,3)+t806*TT(3,3)
      S(100) = t820*TT(1,1)+t825*TT(2,1)+t830*TT(3,1)
      S(101) = t820*TT(1,2)+t825*TT(2,2)+t830*TT(3,2)
      S(102) = t820*TT(1,3)+t825*TT(2,3)+t830*TT(3,3)
      S(103) = t844*TT(1,1)+t849*TT(2,1)+t854*TT(3,1)
      S(104) = t844*TT(1,2)+t849*TT(2,2)+t854*TT(3,2)
      S(105) = t844*TT(1,3)+t849*TT(2,3)+t854*TT(3,3)
      S(106) = t868*TT(1,2)+t873*TT(2,2)+t878*TT(3,2)
      S(107) = t868*TT(1,3)+t873*TT(2,3)+t878*TT(3,3)
      S(108) = t888*TT(1,1)+t893*TT(2,1)+t898*TT(3,1)
      S(109) = t888*TT(1,2)+t893*TT(2,2)+t898*TT(3,2)
      S(110) = t888*TT(1,3)+t893*TT(2,3)+t898*TT(3,3)
      S(111) = t912*TT(1,1)+t917*TT(2,1)+t922*TT(3,1)
      S(112) = t912*TT(1,2)+t917*TT(2,2)+t922*TT(3,2)
      S(113) = t912*TT(1,3)+t917*TT(2,3)+t922*TT(3,3)
      S(114) = t936*TT(1,1)+t941*TT(2,1)+t946*TT(3,1)
      S(115) = t936*TT(1,2)+t941*TT(2,2)+t946*TT(3,2)
      S(116) = t936*TT(1,3)+t941*TT(2,3)+t946*TT(3,3)
      S(117) = (TT(1,3)*ST(94)+TT(2,3)*ST(95)+TT(3,3)*ST(96))*TT(1,3)+(T
     #T(1,3)*ST(95)+TT(2,3)*ST(106)+TT(3,3)*ST(107))*TT(2,3)+(TT(1,3)*ST
     #(96)+TT(2,3)*ST(107)+TT(3,3)*ST(117))*TT(3,3)
      S(118) = t976*TT(1,1)+t981*TT(2,1)+t986*TT(3,1)
      S(119) = t976*TT(1,2)+t981*TT(2,2)+t986*TT(3,2)
      S(120) = t976*TT(1,3)+t981*TT(2,3)+t986*TT(3,3)
      S(121) = t1000*TT(1,1)+t1005*TT(2,1)+t1010*TT(3,1)
      S(122) = t1000*TT(1,2)+t1005*TT(2,2)+t1010*TT(3,2)
      S(123) = t1000*TT(1,3)+t1005*TT(2,3)+t1010*TT(3,3)
      S(124) = t1024*TT(1,1)+t1029*TT(2,1)+t1034*TT(3,1)
      S(125) = t1024*TT(1,2)+t1029*TT(2,2)+t1034*TT(3,2)
      S(126) = t1024*TT(1,3)+t1029*TT(2,3)+t1034*TT(3,3)
      S(127) = t1048*TT(1,1)+t1053*TT(2,1)+t1058*TT(3,1)
      S(128) = t1048*TT(1,2)+t1053*TT(2,2)+t1058*TT(3,2)
      S(129) = t1048*TT(1,3)+t1053*TT(2,3)+t1058*TT(3,3)
      S(130) = t1072*TT(1,1)+t1077*TT(2,1)+t1082*TT(3,1)
      S(131) = t1072*TT(1,2)+t1077*TT(2,2)+t1082*TT(3,2)
      S(132) = t1072*TT(1,3)+t1077*TT(2,3)+t1082*TT(3,3)
      S(133) = t1096*TT(1,1)+t1101*TT(2,1)+t1106*TT(3,1)
      S(134) = t1096*TT(1,2)+t1101*TT(2,2)+t1106*TT(3,2)
      S(135) = t1096*TT(1,3)+t1101*TT(2,3)+t1106*TT(3,3)
      S(136) = t1120*TT(1,2)+t1125*TT(2,2)+t1130*TT(3,2)
      S(137) = t1120*TT(1,3)+t1125*TT(2,3)+t1130*TT(3,3)
      S(138) = t1140*TT(1,1)+t1145*TT(2,1)+t1150*TT(3,1)
      S(139) = t1140*TT(1,2)+t1145*TT(2,2)+t1150*TT(3,2)
      S(140) = t1140*TT(1,3)+t1145*TT(2,3)+t1150*TT(3,3)
      S(141) = t1164*TT(1,1)+t1169*TT(2,1)+t1174*TT(3,1)
      S(142) = t1164*TT(1,2)+t1169*TT(2,2)+t1174*TT(3,2)
      S(143) = t1164*TT(1,3)+t1169*TT(2,3)+t1174*TT(3,3)
      S(144) = (TT(1,3)*ST(127)+TT(2,3)*ST(128)+TT(3,3)*ST(129))*TT(1,3)
     #+(TT(1,3)*ST(128)+TT(2,3)*ST(136)+TT(3,3)*ST(137))*TT(2,3)+(TT(1,3
     #)*ST(129)+TT(2,3)*ST(137)+TT(3,3)*ST(144))*TT(3,3)
      S(145) = t1204*TT(1,1)+t1209*TT(2,1)+t1214*TT(3,1)
      S(146) = t1204*TT(1,2)+t1209*TT(2,2)+t1214*TT(3,2)
      S(147) = t1204*TT(1,3)+t1209*TT(2,3)+t1214*TT(3,3)
      S(148) = t1228*TT(1,1)+t1233*TT(2,1)+t1238*TT(3,1)
      S(149) = t1228*TT(1,2)+t1233*TT(2,2)+t1238*TT(3,2)
      S(150) = t1228*TT(1,3)+t1233*TT(2,3)+t1238*TT(3,3)
      S(151) = t1252*TT(1,1)+t1257*TT(2,1)+t1262*TT(3,1)
      S(152) = t1252*TT(1,2)+t1257*TT(2,2)+t1262*TT(3,2)
      S(153) = t1252*TT(1,3)+t1257*TT(2,3)+t1262*TT(3,3)
      S(154) = t1276*TT(1,1)+t1281*TT(2,1)+t1286*TT(3,1)
      S(155) = t1276*TT(1,2)+t1281*TT(2,2)+t1286*TT(3,2)
      S(156) = t1276*TT(1,3)+t1281*TT(2,3)+t1286*TT(3,3)
      S(157) = t1300*TT(1,2)+t1305*TT(2,2)+t1310*TT(3,2)
      S(158) = t1300*TT(1,3)+t1305*TT(2,3)+t1310*TT(3,3)
      S(159) = t1320*TT(1,1)+t1325*TT(2,1)+t1330*TT(3,1)
      S(160) = t1320*TT(1,2)+t1325*TT(2,2)+t1330*TT(3,2)
      S(161) = t1320*TT(1,3)+t1325*TT(2,3)+t1330*TT(3,3)
      S(162) = (TT(1,3)*ST(151)+TT(2,3)*ST(152)+TT(3,3)*ST(153))*TT(1,3)
     #+(TT(1,3)*ST(152)+TT(2,3)*ST(157)+TT(3,3)*ST(158))*TT(2,3)+(TT(1,3
     #)*ST(153)+TT(2,3)*ST(158)+TT(3,3)*ST(162))*TT(3,3)
      S(163) = t1360*TT(1,1)+t1365*TT(2,1)+t1370*TT(3,1)
      S(164) = t1360*TT(1,2)+t1365*TT(2,2)+t1370*TT(3,2)
      S(165) = t1360*TT(1,3)+t1365*TT(2,3)+t1370*TT(3,3)
      S(166) = t1384*TT(1,1)+t1389*TT(2,1)+t1394*TT(3,1)
      S(167) = t1384*TT(1,2)+t1389*TT(2,2)+t1394*TT(3,2)
      S(168) = t1384*TT(1,3)+t1389*TT(2,3)+t1394*TT(3,3)
      S(169) = t1408*TT(1,2)+t1413*TT(2,2)+t1418*TT(3,2)
      S(170) = t1408*TT(1,3)+t1413*TT(2,3)+t1418*TT(3,3)
      S(171) = (TT(1,3)*ST(166)+TT(2,3)*ST(167)+TT(3,3)*ST(168))*TT(1,3)
     #+(TT(1,3)*ST(167)+TT(2,3)*ST(169)+TT(3,3)*ST(170))*TT(2,3)+(TT(1,3
     #)*ST(168)+TT(2,3)*ST(170)+TT(3,3)*ST(171))*TT(3,3)

      RETURN
      END
C
C=====================================================================





C======================================================================
C	START 3 NODE SHELL ELEMENT - BENDING SUBROUTINES
C======================================================================
C=====================================================================
      SUBROUTINE TRICB(RS,SLN,RN,SN,BL,BLR,BLS,FLR,FLS,DWR,DWS,CB)
	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cb FOR BENDING STRAINS
C
C	INPUT VARIABLES
C	RS(2,3)			= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)			= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	    = ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	BL(3),BLR,BLS	= SHEAR DEFORMATION PARAMETERS ALONG THE BOUNDARY,
C						R, S, RESP.
C	FLR,FLS,
C	DWR,DWS			= COMMON FACTORS DEFINED IN SUBROUTINE COMFACT
C
C	OUTPUT VARIAVBLES
C	CB(9,18)		= BENDING STRAIN PARAMETER
C	--------------------------------------------------------------
      DIMENSION RS(2,3),SLN(3),RN(3),SN(3),BL(3)
	DIMENSION DWR(9),DWS(9),FLR(3),FLS(3)
	DIMENSION CB(9,18)

C	--------------------------------
C	SOLVE FOR Cb FOR BENDING STRAINS
C	--------------------------------

C	FIRST ROW
C	------------------
	CB(1,4) =
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2)
	CB(1,10)=
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(2)*RN(2)**2*SN(2)*BL(2)/2)
	CB(1,16)=
     #(-SLN(2)*RN(2)**2*SN(2)*BL(2)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2)

	CB(1,5) =
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))*SLN(1)+
     #(RN(3)**3/2+(1.E0/2.E0-BL(3)/2)*SN(3)**2*RN(3))*SLN(3))
	CB(1,11)=
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))*SLN(1)+
     #(RN(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2))
	CB(1,17)=
     #((RN(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2)+
     #(RN(3)**3/2+(-BL(3)/2+1.E0/2.E0)*SN(3)**2*RN(3))*SLN(3))

	CB(1,3) =(RN(3)*SN(3)*BL(3)-RN(1)*SN(1)*BL(1))
	CB(1,9) =(RN(1)*SN(1)*BL(1)-RN(2)*SN(2)*BL(2))
	CB(1,15)=(RN(2)*SN(2)*BL(2)-RN(3)*SN(3)*BL(3))

C	SECOND ROW
C	------------------
	CB(2,4) =
     #((-SLN(1)*RS(1,1)/4-SLN(1)*RS(1,2)/4)*BL(1)*SN(1)*RN(1)**2+(-
     #SLN(3)*RS(1,1)/4-SLN(3)*RS(1,3)/4)*BL(3)*SN(3)*RN(3)**2)-DWR(1)
	CB(2,10)=((
     #-SLN(1)*RS(1,1)/4-SLN(1)*RS(1,2)/4)*BL(1)*SN(1)*RN(1)**2+(-SLN(2)*
     #RS(1,2)/4-SLN(2)*RS(1,3)/4)*BL(2)*SN(2)*RN(2)**2)-DWR(2)
	CB(2,16)=((-SLN(2)
     #*RS(1,2)/4-SLN(2)*RS(1,3)/4)*BL(2)*SN(2)*RN(2)**2+(-SLN(3)*RS(1,3)
     #/4-SLN(3)*RS(1,3)/4)*BL(3)*SN(3)*RN(3)**2)-DWR(3)

	CB(2,5) =
     #((SLN(1)*RS(1,2)/6+SLN(1)*RS(1,1)/3)*RN(1)**3+((-SLN(1)*RS(1,
     #1)/4-SLN(1)*RS(1,2)/4)*BL(1)+SLN(1)*RS(1,2)/6+SLN(1)*RS(1,1)/3)*SN
     #(1)**2*RN(1)+(SLN(3)*RS(1,3)/6+SLN(3)*RS(1,1)/3)*RN(3)**3+((-SLN(3
     #)*RS(1,1)/4-SLN(3)*RS(1,3)/4)*BL(3)+SLN(3)*RS(1,3)/6+SLN(3)*RS(1,1
     #)/3)*SN(3)**2*RN(3))-DWR(4)-FLS(1)
	CB(2,11)=
     #((SLN(1)*RS(1,2)/3+SLN(1)*RS(1,1)/6)*RN(1)**3+((-SLN(1)*RS(1,
     #1)/4-SLN(1)*RS(1,2)/4)*BL(1)+SLN(1)*RS(1,2)/3+SLN(1)*RS(1,1)/6)*SN
     #(1)**2*RN(1)+(SLN(2)*RS(1,3)/6+SLN(2)*RS(1,2)/3)*RN(2)**3+((-SLN(2
     #)*RS(1,2)/4-SLN(2)*RS(1,3)/4)*BL(2)+SLN(2)*RS(1,3)/6+SLN(2)*RS(1,2
     #)/3)*SN(2)**2*RN(2))-DWR(5)-FLS(2)
	CB(2,17)=((SLN(2)*RS(1,2)/6+SLN(2)*RS(1,3)/3)*RN(2)**3+
     #((-SLN(2)*RS(1,
     #2)/4-SLN(2)*RS(1,3)/4)*BL(2)+SLN(2)*RS(1,2)/6+SLN(2)*RS(1,3)/3)*SN
     #(2)**2*RN(2)+(SLN(3)*RS(1,3)/6+SLN(3)*RS(1,3)/3)*RN(3)**3+((-SLN(3
     #)*RS(1,3)/4-SLN(3)*RS(1,3)/4)*BL(3)+SLN(3)*RS(1,3)/6+SLN(3)*RS(1,3
     #)/3)*SN(3)**2*RN(3))-DWR(6)-FLS(3)

	CB(2,3) =((-RS(1,1)-RS(1,2))*BL(1)*SN(1)*RN(1)+(RS(1,3)+RS(1,
     #1))*BL(3)*SN(3)*RN(3))/2-DWR(7)
	CB(2,9) =((RS(1,2)+RS(1,1))*BL(1)*SN(1)*
     #RN(1)+(-RS(1,3)-RS(1,2))*BL(2)*SN(2)*RN(2))/2-DWR(8)
	CB(2,15)=((RS(1,3)+RS(1,2))*BL(2)*SN(2)*RN(2)+
     #(-RS(1,3)-RS(1,3))*BL(3)*SN(3)*RN(3))/2-DWR(9)

C	THIRD ROW
C	------------------
      t1 = -RS(2,2)-RS(2,1)
      t2 = t1*BL(1)/4
      t4 = RN(1)**2
      t6 = t2*SLN(1)*SN(1)*t4
      t7 = -RS(2,1)-RS(2,3)
      t8 = t7*BL(3)/4
      t10 = RN(3)**2
      t12 = t8*SLN(3)*SN(3)*t10
      t15 = -RS(2,2)-RS(2,3)
      t16 = t15*BL(2)/4
      t18 = RN(2)**2
      t20 = t16*SLN(2)*SN(2)*t18
      t23 = -RS(2,3)-RS(2,3)
      t24 = t23*BL(3)/4
      t26 = RN(3)**2
      t28 = t24*SLN(3)*SN(3)*t26
      t33 = RS(2,1)/3
      t34 = RS(2,2)/6
      t37 = t4*RN(1)
      t41 = SN(1)**2
      t42 = t41*RN(1)
      t44 = RS(2,3)/6
      t47 = t10*RN(3)
      t51 = SN(3)**2
      t52 = t51*RN(3)
      t56 = RS(2,2)/3
      t57 = RS(2,1)/6
      t64 = RS(2,3)/6
      t67 = t18*RN(2)
      t71 = SN(2)**2
      t72 = t71*RN(2)
      t76 = RS(2,3)/3
      t85 = t26*RN(3)
      t89 = SN(3)**2
      t90 = t89*RN(3)
      t94 = RS(2,3)/3
      t110 = SN(1)*RN(1)
      t113 = SN(3)*RN(3)
      t120 = SN(2)*RN(2)
      t127 = SN(3)*RN(3)

	CB(3,4) =(t6+t12)
	CB(3,10)=(t6+t20)
	CB(3,16)=(t20+t28)

	CB(3,5) =((t33+t34)*SLN(1)*t37+(t2+t33+t34)*SLN(1)*t42+
     #(t33+t44)*SLN(3
     #)*t47+(t8+t33+t44)*SLN(3)*t52)
 	CB(3,11)=((t56+t57)*SLN(1)*t37+(t2+t5
     #6+t57)*SLN(1)*t42+(t56+t64)*SLN(2)*t67+(t16+t56+t64)*SLN(2)*t72)
	CB(3,17)=((t76+t34)*SLN(2)*t67+(t16+t76+t34)*SLN(2)*t72+(t44+t76)
     #*SLN(3)*t85+(t24+t44+t76)*SLN(3)*t90)

	CB(3,3) =(t1*BL(1)*t110-t7*BL(3)*t113)/2
	CB(3,9) =(t15*BL(2)*t120-t1*BL(1)*t110)/2
	CB(3,15)=(t23*BL(3)*t127-t15*BL(2)*t120)/2

C	FOURTH ROW
C	------------------
	CB(4,4) =
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(BL(3)/2-1.E0/2.E0
     #)*SLN(3)*SN(3)*RN(3)**2-SLN(1)*SN(1)**3/2-SLN(3)*SN(3)**3/2)
	CB(4,10)=
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(-1.E0/2.E0+BL(2)/2)*
     #SLN(2)*SN(2)*RN(2)**2-SLN(1)*SN(1)**3/2-SLN(2)*SN(2)**3/2)
	CB(4,16)=
     #((-1.E0/2.E0+BL(2)/2)*SLN(2)*SN(2)*RN(2)**2+(BL(3)/2-1.E0/2.E0)*SL
     #N(3)*SN(3)*RN(3)**2-SLN(2)*SN(2)**3/2-SLN(3)*SN(3)**3/2)

	CB(4,5) =
     #(SLN(1)*SN(1)**2*RN(1)*BL(1)/2+SLN(3)*SN(3)**2*RN(3)*BL(3)/2)
	CB(4,11)=
     #(SLN(2)*SN(2)**2*RN(2)*BL(2)/2+SLN(1)*SN(1)**2*RN(1)*BL(1)/2)
	CB(4,17)=
     #(SLN(3)*SN(3)**2*RN(3)*BL(3)/2+SLN(2)*SN(2)**2*RN(2)*BL(2)/2)

	CB(4,3) =(SN(1)*RN(1)*BL(1)-SN(3)*RN(3)*BL(3))
	CB(4,9) =(SN(2)*RN(2)*BL(2)-SN(1)*RN(1)*BL(1))
	CB(4,15)=(SN(3)*RN(3)*BL(3)-SN(2)*RN(2)*BL(2))

C	FIFTH ROW
C	------------------
	CB(5,4) =
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,1)/3-RS(1,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(1,1)/3-RS(1,2)/6)*SLN(1)*SN(1)**3+(-RS(1,1)
     #/3-RS(1,3)/6)*SLN(3)*SN(3)**3)
	CB(5,10)=
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1)-RS(1,2)/3-RS(1,1)/6)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,2)/4+RS(1,3)/4)*BL(2)-RS(1,2)/3-RS(1,3)/6)*SLN(
     #2)*SN(2)*RN(2)**2+(-RS(1,2)/3-RS(1,1)/6)*SLN(1)*SN(1)**3+(-RS(1,2)
     #/3-RS(1,3)/6)*SLN(2)*SN(2)**3)
	CB(5,16)=(((RS(1,2)/4+RS(1,3)/4)*BL(2
     #)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*SN(2)*RN(2)**2+((RS(1,3)/4+RS(1,3)/4
     #)*BL(3)-RS(1,3)/6-RS(1,3)/3)*SLN(3)*SN(3)*RN(3)**2+(-RS(1,3)/3-RS(
     #1,2)/6)*SLN(2)*SN(2)**3+(-RS(1,3)/6-RS(1,3)/3)*SLN(3)*SN(3)**3)

	CB(5,5) =((RS(1,2)/4+RS(1,1)/4)*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,1)/4+RS(1,3)/4)*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))
	CB(5,11)=((RS(1,2)/4+RS(1,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(
     #1)+(RS(1,2)/4+RS(1,3)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))
	CB(5,17)=
     #((RS(1,2)/4+RS(1,3)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(1,3
     #)/4+RS(1,3)/4)*BL(3)*SLN(3)*SN(3)**2*RN(3))

	CB(5,3) =
     #((RS(1,2)/2+RS(1,1)/2)*BL(1)*SN(1)*RN(1)+(-RS(1,1)/2-RS(1,
     #3)/2)*BL(3)*SN(3)*RN(3))
	CB(5,9) =((-RS(1,2)/2-RS(1,1)/2)*BL(1)*SN(1)
     #*RN(1)+(RS(1,3)/2+RS(1,2)/2)*BL(2)*SN(2)*RN(2))
	CB(5,15)=((-RS(1,2)/2
     #-RS(1,3)/2)*BL(2)*SN(2)*RN(2)+(RS(1,3)/2+RS(1,3)/2)*BL(3)*SN(3)*RN
     #(3))

C	SIXTH ROW
C	------------------
	CB(6,4) =
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(1)**3+(-RS(2,1)
     #/3-RS(2,3)/6)*SLN(3)*SN(3)**3)-DWS(1)-FLR(1)
	CB(6,10)=
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/3-RS(2,3)/6)*SLN(
     #2)*SN(2)*RN(2)**2+(-RS(2,1)/6-RS(2,2)/3)*SLN(1)*SN(1)**3+(-RS(2,3)
     #/6-RS(2,2)/3)*SLN(2)*SN(2)**3)-DWS(2)-FLR(2)
	CB(6,16)=(((RS(2,3)/4+RS(2,2)/4)*BL(2
     #)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN(2)*RN(2)**2+((RS(2,3)/4+RS(2,3)/4
     #)*BL(3)-RS(2,3)/3-RS(2,3)/6)*SLN(3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(
     #2,3)/3)*SLN(2)*SN(2)**3+(-RS(2,3)/3-RS(2,3)/6)*SLN(3)*SN(3)**3)-
     #DWS(3)-FLR(3)

	CB(6,5) =((RS(2,2)+RS(2,1))*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)+RS(2,1))*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))/4-DWS(4)
	CB(6,11)=((RS(2,2)+RS(2,1))*BL(1)*SLN(1)*SN(1)**2*RN(
     #1)+(RS(2,3)+RS(2,2))*BL(2)*SLN(2)*SN(2)**2*RN(2))/4-DWS(5)
	CB(6,17)=
     #((RS(2,3)+RS(2,2))*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(2,3
     #)+RS(2,3))*BL(3)*SLN(3)*SN(3)**2*RN(3))/4-DWS(6)

	CB(6,3) =((RS(2,2)+RS(2,1))*BL(1)*SN(1)*RN(1)+(-RS(2,1)-RS(2,
     #3))*BL(3)*SN(3)*RN(3))/2-DWS(7)
	CB(6,9) =((-RS(2,2)-RS(2,1))*BL(1)*SN(1)
     #*RN(1)+(RS(2,2)+RS(2,3))*BL(2)*SN(2)*RN(2))/2-DWS(8)
	CB(6,15)=((-RS(2,3)
     #-RS(2,2))*BL(2)*SN(2)*RN(2)+(RS(2,3)+RS(2,3))*BL(3)*SN(3)*RN
     #(3))/2-DWS(9)

C	SEVENTH ROW
C	------------------
	CB(7,4) =
     #((-1.E0/2.E0+BL(1)/2)*SLN(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SL
     #N(1)*SN(1)**2*RN(1)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-BL(3)/2-
     #1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))
	CB(7,10)=((-1.E0/2.E0+BL(1)/2)*SLN
     #(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SLN(1)*SN(1)**2*RN(1)+(-1.E0/2.E
     #0+BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2.E0)*SLN(2)*SN(2)**2*RN
     #(2))
	CB(7,16)=((-1.E0/2.E0+BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2
     #.E0)*SLN(2)*SN(2)**2*RN(2)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-B
     #L(3)/2-1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))

	CB(7,5) =((BL(1)/2+1.E0/2.E0)*SL
     #N(1)*SN(1)*RN(1)**2+(1.E0/2.E0+BL(3)/2)*SLN(3)*SN(3)*RN(3)**2+(-BL
     #(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3+(1.E0/2.E0-BL(3)/2)*SLN(3)*SN(3)*
     #*3)
	CB(7,11)=((BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(1.E0/2.E0+B
     #L(2)/2)*SLN(2)*SN(2)*RN(2)**2+(-BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3
     #+(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3)
	CB(7,17)=
     #((1.E0/2.E0+BL(2)/2)*SLN(2)*SN(2)*RN(2)**2+(BL(3)/2+1.E0/2
     #.E0)*SLN(3)*SN(3)*RN(3)**2+(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3+(1.
     #E0/2.E0-BL(3)/2)*SLN(3)*SN(3)**3)

	CB(7,3) =
     #(RN(1)**2*BL(1)-SN(1)**2*BL(1)+SN(3)**2*BL(3)-RN(3)**2*BL(3))
	CB(7,9) =
     #(-RN(1)**2*BL(1)-SN(2)**2*BL(2)+SN(1)**2*BL(1)+RN(2)**2*BL(2))
	CB(7,15)=
     #(-RN(2)**2*BL(2)-SN(3)**2*BL(3)+SN(2)**2*BL(2)+RN(3)**2*BL(3))

C	EIGHTH ROW
C	------------------
	CB(8,4) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,2)/6-RS(1,1)/3)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,1)/3-RS(1,3)/6)*SLN
     #(3)*RN(3)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)-RS(1,1)/3-RS(1,3)/6)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWS(1)-FLR(1)
	CB(8,10) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,2)/3-RS(1,3)/6)*SLN
     #(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,2)/3-RS(1,3)/6)*SL
     #N(2)*SN(2)**2*RN(2))
     #-DWS(2)-FLR(2)
	CB(8,16)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*RN(
     #2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(1,3)/4+RS(1,3)/4)*BL(3)-RS(1,3)/3-RS(1,3)/6)*SLN
     #(3)*RN(3)**3+((-RS(1,3)/4-RS(1,3)/4)*BL(3)-RS(1,3)/3-RS(1,3)/6)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWS(3)-FLR(3)

	CB(8,5) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,1)/3+RS(1,2)/6)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,1)/3+RS(1,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,1)/3+RS(1,2)/
     #6)*SLN(1)*SN(1)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/3+RS(1,3)
     #/6)*SLN(3)*SN(3)**3)-DWS(4)
	CB(8,11)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,3)/6+RS(1,2)/3)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,1)/6+RS(1,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,3)/6+RS(1,2)
     #/3)*SLN(2)*SN(2)**3)-DWS(5)
	CB(8,17)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,2)/6+RS(1,3)/3)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(1,3)/4+RS(1,3)/4)*BL(3)+RS(1,3)/6+RS(1,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/
     #3)*SLN(2)*SN(2)**3+((-RS(1,3)/4-RS(1,3)/4)*BL(3)+RS(1,3)/6+RS(1,3)
     #/3)*SLN(3)*SN(3)**3)-DWS(6)

	CB(8,3) =
     #((RS(1,2)/2+RS(1,1)/2)*BL(1)*RN(1)**2+(-RS(1,3)/2-RS(1,1)/
     #2)*BL(3)*RN(3)**2+(-RS(1,2)/2-RS(1,1)/2)*BL(1)*SN(1)**2+(RS(1,3)/2
     #+RS(1,1)/2)*BL(3)*SN(3)**2)-DWS(7)
	CB(8,9) =
     #((-RS(1,2)/2-RS(1,1)/2)*BL(1)*RN(1)**2+(RS(1,3)/2+RS(1,2)/
     #2)*BL(2)*RN(2)**2+(RS(1,2)/2+RS(1,1)/2)*BL(1)*SN(1)**2+(-RS(1,3)/2
     #-RS(1,2)/2)*BL(2)*SN(2)**2)-DWS(8)
	CB(8,15)=((-RS(1,3)/2-RS(1,2)/2)*BL(2)*RN
     #(2)**2+(RS(1,3)/2+RS(1,3)/2)*BL(3)*RN(3)**2+(RS(1,3)/2+RS(1,2)/2)*
     #BL(2)*SN(2)**2+(-RS(1,3)/2-RS(1,3)/2)*BL(3)*SN(3)**2)-DWS(9)

C	NINTH ROW
C	------------------
	CB(9,4) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,1)/3-RS(2,2)/6)*SLN(1)*RN(
     #1)**3+((-RS(2,2)/4-RS(2,1)/4)*BL(1)-RS(2,1)/3-RS(2,2)/6)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,1)/4+RS(2,3)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SLN
     #(3)*RN(3)**3+((-RS(2,1)/4-RS(2,3)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SL
     #N(3)*SN(3)**2*RN(3))-DWR(1)
	CB(9,10)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(2,2)/4-RS(2,1)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,3)/6-RS(2,2)/3)*SLN
     #(2)*RN(2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,3)/6-RS(2,2)/3)*SL
     #N(2)*SN(2)**2*RN(2))-DWR(2)
	CB(9,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,3)/3-RS(2,2)/6)*SLN(2)*RN(
     #2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,3)/3-RS(2,2)/6)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(2,3)/4+RS(2,3)/4)*BL(3)-RS(2,3)/6-RS(2,3)/3)*SLN
     #(3)*RN(3)**3+((-RS(2,3)/4-RS(2,3)/4)*BL(3)-RS(2,3)/6-RS(2,3)/3)*SL
     #N(3)*SN(3)**2*RN(3))-DWR(3)

	CB(9,5) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,1)/3+RS(2,2)/6)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,1)/4+RS(2,3)/4)*BL(3)+RS(2,1)/3+RS(2,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,2)/4-RS(2,1)/4)*BL(1)+RS(2,1)/3+RS(2,2)/
     #6)*SLN(1)*SN(1)**3+((-RS(2,1)/4-RS(2,3)/4)*BL(3)+RS(2,1)/3+RS(2,3)
     #/6)*SLN(3)*SN(3)**3)-DWR(4)-FLS(1)
	CB(9,11)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,1)/6+RS(2,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)/6)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(2,2)/4-RS(2,1)/4)*BL(1)+RS(2,1)/6+RS(2,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)
     #/6)*SLN(2)*SN(2)**3)-DWR(5)-FLS(2)
	CB(9,17)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/6)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(2,3)/4+RS(2,3)/4)*BL(3)+RS(2,3)/6+RS(2,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/
     #6)*SLN(2)*SN(2)**3+((-RS(2,3)/4-RS(2,3)/4)*BL(3)+RS(2,3)/3+RS(2,3)
     #/6)*SLN(3)*SN(3)**3)-DWR(6)-FLS(3)

	CB(9,3) =
     #((RS(2,1)/2+RS(2,2)/2)*BL(1)*RN(1)**2+(-RS(2,1)/2-RS(2,3)/
     #2)*BL(3)*RN(3)**2+(-RS(2,1)/2-RS(2,2)/2)*BL(1)*SN(1)**2+(RS(2,1)/2
     #+RS(2,3)/2)*BL(3)*SN(3)**2)-DWR(7)
	CB(9,9) =
     #((-RS(2,1)/2-RS(2,2)/2)*BL(1)*RN(1)**2+(RS(2,2)/2+RS(2,3)/
     #2)*BL(2)*RN(2)**2+(RS(2,1)/2+RS(2,2)/2)*BL(1)*SN(1)**2+(-RS(2,3)/2
     #-RS(2,2)/2)*BL(2)*SN(2)**2)-DWR(8)
	CB(9,15)=((-RS(2,3)/2-RS(2,2)/2)*BL(2)*RN
     #(2)**2+(RS(2,3)/2+RS(2,3)/2)*BL(3)*RN(3)**2+(RS(2,2)/2+RS(2,3)/2)*
     #BL(2)*SN(2)**2+(-RS(2,3)/2-RS(2,3)/2)*BL(3)*SN(3)**2)-DWR(9)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRICB3(ANT,RS,SLN,RN,SN,BL,BLR,BLS,CB)
	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cb FOR BENDING STRAINS
C
C	INPUT VARIABLES
C	RS(2,3)			= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)			= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	    = ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	BL(3),BLR,BLS	= SHEAR DEFORMATION PARAMETERS ALONG THE BOUNDARY,
C						R, S, RESP.
C	FLR,FLS,
C	DWR,DWS			= COMMON FACTORS DEFINED IN SUBROUTINE COMFACT
C
C	OUTPUT VARIAVBLES
C	CB(9,18)		= BENDING STRAIN PARAMETER
C	--------------------------------------------------------------
      DIMENSION ANT(3),RS(2,3),SLN(3),RN(3),SN(3),BL(3)
	DIMENSION CB(9,18)

C	--------------------------------
C	SOLVE FOR Cb FOR BENDING STRAINS
C	--------------------------------

C	FIRST ROW
C	------------------
	CB(1,4) =
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2)
	CB(1,10)=
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(2)*RN(2)**2*SN(2)*BL(2)/2)
	CB(1,16)=
     #(-SLN(2)*RN(2)**2*SN(2)*BL(2)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2)

	CB(1,5) =
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))*SLN(1)+
     #(RN(3)**3/2+(1.E0/2.E0-BL(3)/2)*SN(3)**2*RN(3))*SLN(3))
	CB(1,11)=
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))*SLN(1)+
     #(RN(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2))
	CB(1,17)=
     #((RN(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2)+
     #(RN(3)**3/2+(-BL(3)/2+1.E0/2.E0)*SN(3)**2*RN(3))*SLN(3))

	CB(1,3) =(RN(3)*SN(3)*BL(3)-RN(1)*SN(1)*BL(1))
	CB(1,9) =(RN(1)*SN(1)*BL(1)-RN(2)*SN(2)*BL(2))
	CB(1,15)=(RN(2)*SN(2)*BL(2)-RN(3)*SN(3)*BL(3))

C	SECOND ROW
C	------------------
	CB(2,4) =
     #((-RS(1,1)/4-RS(1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,1)
     #/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)
	CB(2,10)=
     #((-RS(1,1)/4-RS(
     #1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,2)/4-RS(1,3)/4)*BL(2)*S
     #LN(2)*SN(2)*RN(2)**2)
	CB(2,16)=
     #((-RS(1,2)/4-RS(1,3)/4)*BL(2)*SLN(2)*
     #SN(2)*RN(2)**2+(-RS(1,1)/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)

	CB(2,5) =
     #((RS(1,2)/6+RS(1,1)/3)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2
     #)/4)*BL(1)+RS(1,2)/6+RS(1,1)/3)*SLN(1)*SN(1)**2*RN(1)+(RS(1,1)/3+R
     #S(1,3)/6)*SLN(3)*RN(3)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/3+
     #RS(1,3)/6)*SLN(3)*SN(3)**2*RN(3))-ANT(1)
	CB(2,11)=
     #((RS(1,1)/6+RS(1,2)/3)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2
     #)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(1)**2*RN(1)+(RS(1,2)/3+R
     #S(1,3)/6)*SLN(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/3+
     #RS(1,3)/6)*SLN(2)*SN(2)**2*RN(2))-ANT(1)
	CB(2,17)=
     #((RS(1,2)/6+RS(1,3)/3)*SL
     #N(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/3)*S
     #LN(2)*SN(2)**2*RN(2)+(RS(1,1)/6+RS(1,3)/3)*SLN(3)*RN(3)**3+((-RS(1
     #,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/6+RS(1,3)/3)*SLN(3)*SN(3)**2*RN(3))
     #-ANT(1)

	CB(2,3) =((-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)*RN(1)+(RS(1,1)/2+RS(1,
     #3)/2)*BL(3)*SN(3)*RN(3))
	CB(2,9) =((RS(1,1)/2+RS(1,2)/2)*BL(1)*SN(1)*
     #RN(1)+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)*RN(2))
	CB(2,15)=((RS(1,3)/2+
     #RS(1,2)/2)*BL(2)*SN(2)*RN(2)+(-RS(1,3)/2-RS(1,1)/2)*BL(3)*SN(3)*RN
     #(3))

C	THIRD ROW
C	------------------
	CB(3,4) =
     #((-RS(2,2)/4-RS(2,1)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,1)
     #/4-RS(2,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)
	CB(3,10)=
     #((-RS(2,2)/4-RS(
     #2,1)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,2)/4-RS(2,3)/4)*BL(2)*S
     #LN(2)*SN(2)*RN(2)**2)
	CB(3,16)=
     #((-RS(2,2)/4-RS(2,3)/4)*BL(2)*SLN(2)*
     #SN(2)*RN(2)**2+(-RS(2,1)/4-RS(2,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)

	CB(3,5) =
     #((RS(2,1)/3+RS(2,2)/6)*SLN(1)*RN(1)**3+((-RS(2,2)/4-RS(2,1
     #)/4)*BL(1)+RS(2,1)/3+RS(2,2)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/6+R
     #S(2,1)/3)*SLN(3)*RN(3)**3+((-RS(2,1)/4-RS(2,3)/4)*BL(3)+RS(2,3)/6+
     #RS(2,1)/3)*SLN(3)*SN(3)**2*RN(3))
 	CB(3,11)=
     #((RS(2,2)/3+RS(2,1)/6)*SLN(1)*RN(1)**3+((-RS(2,2)/4-RS(2,1
     #)/4)*BL(1)+RS(2,2)/3+RS(2,1)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/6+R
     #S(2,2)/3)*SLN(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,3)/6+
     #RS(2,2)/3)*SLN(2)*SN(2)**2*RN(2))
	CB(3,17)=
     #((RS(2,2)/6+RS(2,3)/3)*SL
     #N(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,2)/6+RS(2,3)/3)*S
     #LN(2)*SN(2)**2*RN(2)+(RS(2,1)/6+RS(2,3)/3)*SLN(3)*RN(3)**3+((-RS(2
     #,1)/4-RS(2,3)/4)*BL(3)+RS(2,1)/6+RS(2,3)/3)*SLN(3)*SN(3)**2*RN(3))

	CB(3,3) =
     #((-RS(2,1)/2-RS(2,2)/2)*BL(1)*SN(1)*RN(1)+(RS(2,1)/2+RS(2,
     #3)/2)*BL(3)*SN(3)*RN(3))
	CB(3,9) =
     #((RS(2,1)/2+RS(2,2)/2)*BL(1)*SN(1)*
     #RN(1)+(-RS(2,2)/2-RS(2,3)/2)*BL(2)*SN(2)*RN(2))
	CB(3,15)=
     #((RS(2,2)/2+
     #RS(2,3)/2)*BL(2)*SN(2)*RN(2)+(-RS(2,3)/2-RS(2,1)/2)*BL(3)*SN(3)*RN
     #(3))

C	FOURTH ROW
C	------------------
	CB(4,4) =
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(BL(3)/2-1.E0/2.E0
     #)*SLN(3)*SN(3)*RN(3)**2-SLN(1)*SN(1)**3/2-SLN(3)*SN(3)**3/2)
	CB(4,10)=
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(-1.E0/2.E0+BL(2)/2)*
     #SLN(2)*SN(2)*RN(2)**2-SLN(1)*SN(1)**3/2-SLN(2)*SN(2)**3/2)
	CB(4,16)=
     #((-1.E0/2.E0+BL(2)/2)*SLN(2)*SN(2)*RN(2)**2+(BL(3)/2-1.E0/2.E0)*SL
     #N(3)*SN(3)*RN(3)**2-SLN(2)*SN(2)**3/2-SLN(3)*SN(3)**3/2)

	CB(4,5) =
     #(SLN(1)*SN(1)**2*RN(1)*BL(1)/2+SLN(3)*SN(3)**2*RN(3)*BL(3)/2)
	CB(4,11)=
     #(SLN(2)*SN(2)**2*RN(2)*BL(2)/2+SLN(1)*SN(1)**2*RN(1)*BL(1)/2)
	CB(4,17)=
     #(SLN(3)*SN(3)**2*RN(3)*BL(3)/2+SLN(2)*SN(2)**2*RN(2)*BL(2)/2)

	CB(4,3) =(SN(1)*RN(1)*BL(1)-SN(3)*RN(3)*BL(3))
	CB(4,9) =(SN(2)*RN(2)*BL(2)-SN(1)*RN(1)*BL(1))
	CB(4,15)=(SN(3)*RN(3)*BL(3)-SN(2)*RN(2)*BL(2))

C	FIFTH ROW
C	------------------
	CB(5,4) =
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1)-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,3)/4+RS(1,1)/4)*BL(3)-RS(1,1)/3-RS(1,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN(1)**3+(-RS(1,1)
     #/3-RS(1,3)/6)*SLN(3)*SN(3)**3)
	CB(5,10)=
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1
     #)-RS(1,2)/3-RS(1,1)/6)*SLN(1)*SN(1)*RN(1)**2+((RS(1,3)/4+RS(1,2)/4
     #)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(1,2)/3-RS(
     #1,1)/6)*SLN(1)*SN(1)**3+(-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)**3)
	CB(5,16)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(1,3)/4+RS(1,1)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(1,3)/3-RS(1,2)/6)*SLN(2)*SN(2)**3+(-RS(1
     #,1)/6-RS(1,3)/3)*SLN(3)*SN(3)**3)

	CB(5,5) =
     #((RS(1,2)/4+RS(1,1)/4)*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,3)/4+RS(1,1)/4)*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))
	CB(5,11)=
     #((RS(1,2)/4+RS(1,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,3
     #)/4+RS(1,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))
	CB(5,17)=
     #((RS(1,3)/4+RS(
     #1,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(1,3)/4+RS(1,1)/4)*BL(3)*SL
     #N(3)*SN(3)**2*RN(3))

	CB(5,3) =
     #((RS(1,1)/2+RS(1,2)/2)*BL(1)*SN(1)*RN(
     #1)+(-RS(1,1)/2-RS(1,3)/2)*BL(3)*SN(3)*RN(3))
	CB(5,9) =
     #((-RS(1,2)/2-RS
     #(1,1)/2)*BL(1)*SN(1)*RN(1)+(RS(1,2)/2+RS(1,3)/2)*BL(2)*SN(2)*RN(2)
     #)
	CB(5,15)=
     #((-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)*RN(2)+(RS(1,3)/2+RS(1,
     #1)/2)*BL(3)*SN(3)*RN(3))

C	SIXTH ROW
C	------------------
	CB(6,4) =
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,1)/4+RS(2,3)/4)*BL(3)-RS(2,3)/6-RS(2,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(1)**3+(-RS(2,3)
     #/6-RS(2,1)/3)*SLN(3)*SN(3)**3)
     #+ANT(1)
	CB(6,10)=
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1
     #)-RS(2,2)/3-RS(2,1)/6)*SLN(1)*SN(1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4
     #)*BL(2)-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(2,2)/3-RS(
     #2,1)/6)*SLN(1)*SN(1)**3+(-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)**3)
     #+ANT(1)
	CB(6,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(2,1)/4+RS(2,3)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN(2)**3+(-RS(2
     #,1)/6-RS(2,3)/3)*SLN(3)*SN(3)**3)
     #+ANT(1)

	CB(6,5) =
     #((RS(2,2)/4+RS(2,1)/4)*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,1)/4+RS(2,3)/4)*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))
	CB(6,11)=
     #((RS(2,2)/4+RS(2,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3
     #)/4+RS(2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))
	CB(6,17)=
     #((RS(2,3)/4+RS(
     #2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(2,1)/4+RS(2,3)/4)*BL(3)*SL
     #N(3)*SN(3)**2*RN(3))

	CB(6,3) =
     #((RS(2,1)/2+RS(2,2)/2)*BL(1)*SN(1)*RN(
     #1)+(-RS(2,1)/2-RS(2,3)/2)*BL(3)*SN(3)*RN(3))
	CB(6,9) =
     #((-RS(2,1)/2-RS
     #(2,2)/2)*BL(1)*SN(1)*RN(1)+(RS(2,3)/2+RS(2,2)/2)*BL(2)*SN(2)*RN(2)
     #)
	CB(6,15)=
     #((-RS(2,3)/2-RS(2,2)/2)*BL(2)*SN(2)*RN(2)+(RS(2,3)/2+RS(2,
     #1)/2)*BL(3)*SN(3)*RN(3))

C	SEVENTH ROW
C	------------------
	CB(7,4) =
     #((-1.E0/2.E0+BL(1)/2)*SLN(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SL
     #N(1)*SN(1)**2*RN(1)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-BL(3)/2-
     #1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))
	CB(7,10)=((-1.E0/2.E0+BL(1)/2)*SLN
     #(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SLN(1)*SN(1)**2*RN(1)+(-1.E0/2.E
     #0+BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2.E0)*SLN(2)*SN(2)**2*RN
     #(2))
	CB(7,16)=((-1.E0/2.E0+BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2
     #.E0)*SLN(2)*SN(2)**2*RN(2)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-B
     #L(3)/2-1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))

	CB(7,5) =((BL(1)/2+1.E0/2.E0)*SL
     #N(1)*SN(1)*RN(1)**2+(1.E0/2.E0+BL(3)/2)*SLN(3)*SN(3)*RN(3)**2+(-BL
     #(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3+(1.E0/2.E0-BL(3)/2)*SLN(3)*SN(3)*
     #*3)
	CB(7,11)=((BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(1.E0/2.E0+B
     #L(2)/2)*SLN(2)*SN(2)*RN(2)**2+(-BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3
     #+(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3)
	CB(7,17)=
     #((1.E0/2.E0+BL(2)/2)*SLN(2)*SN(2)*RN(2)**2+(BL(3)/2+1.E0/2
     #.E0)*SLN(3)*SN(3)*RN(3)**2+(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3+(1.
     #E0/2.E0-BL(3)/2)*SLN(3)*SN(3)**3)

	CB(7,3) =
     #(RN(1)**2*BL(1)-SN(1)**2*BL(1)+SN(3)**2*BL(3)-RN(3)**2*BL(3))
	CB(7,9) =
     #(-RN(1)**2*BL(1)-SN(2)**2*BL(2)+SN(1)**2*BL(1)+RN(2)**2*BL(2))
	CB(7,15)=
     #(-RN(2)**2*BL(2)-SN(3)**2*BL(3)+SN(2)**2*BL(2)+RN(3)**2*BL(3))

C	EIGHTH ROW
C	------------------
	CB(8,4) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SLN
     #(3)*RN(3)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SL
     #N(3)*SN(3)**2*RN(3))+ANT(1)
	CB(8,10)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN
     #(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SL
     #N(2)*SN(2)**2*RN(2))+ANT(1)
	CB(8,16)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,2)/6-RS(1,3)/3)*SLN(2)*RN(
     #2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,2)/6-RS(1,3)/3)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*SLN
     #(3)*RN(3)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*SL
     #N(3)*SN(3)**2*RN(3))+ANT(1)

	CB(8,5) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,2)/6+RS(1,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,3)/6+RS(1,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,2)/6+RS(1,1)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)+RS(1,3)/6+RS(1,1)
     #/3)*SLN(3)*SN(3)**3)
	CB(8,11)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,3)/6+RS(1,2)/3)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,1)/6+RS(1,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,3)/6+RS(1,2)
     #/3)*SLN(2)*SN(2)**3)
	CB(8,17)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,2)/6+RS(1,3)/3)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,1)/6+RS(1,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/
     #3)*SLN(2)*SN(2)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)+RS(1,1)/6+RS(1,3)
     #/3)*SLN(3)*SN(3)**3)

	CB(8,3) =
     #((RS(1,2)/2+RS(1,1)/2)*BL(1)*RN(1)**2+(-RS(1,1)/2-RS(1,3)/
     #2)*BL(3)*RN(3)**2+(-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)**2+(RS(1,1)/2
     #+RS(1,3)/2)*BL(3)*SN(3)**2)
	CB(8,9) =
     #((-RS(1,1)/2-RS(1,2)/2)*BL(1)*RN
     #(1)**2+(RS(1,2)/2+RS(1,3)/2)*BL(2)*RN(2)**2+(RS(1,2)/2+RS(1,1)/2)*
     #BL(1)*SN(1)**2+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)**2)
	CB(8,15)=
     #((-RS(
     #1,2)/2-RS(1,3)/2)*BL(2)*RN(2)**2+(RS(1,1)/2+RS(1,3)/2)*BL(3)*RN(3)
     #**2+(RS(1,2)/2+RS(1,3)/2)*BL(2)*SN(2)**2+(-RS(1,1)/2-RS(1,3)/2)*BL
     #(3)*SN(3)**2)

C	NINTH ROW
C	------------------
	CB(9,4) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*RN(
     #1)**3+((-RS(2,1)/4-RS(2,2)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SLN
     #(3)*RN(3)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SL
     #N(3)*SN(3)**2*RN(3))
	CB(9,10)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(2,1)/4-RS(2,2)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/3-RS(2,3)/6)*SLN
     #(2)*RN(2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,2)/3-RS(2,3)/6)*SL
     #N(2)*SN(2)**2*RN(2))
	CB(9,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*RN(
     #2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*SLN
     #(3)*RN(3)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*SL
     #N(3)*SN(3)**2*RN(3))

	CB(9,5) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,2)/6+RS(2,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)+RS(2,3)/6+RS(2,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,1)/4-RS(2,2)/4)*BL(1)+RS(2,2)/6+RS(2,1)/
     #3)*SLN(1)*SN(1)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,3)/6+RS(2,1)
     #/3)*SLN(3)*SN(3)**3)-ANT(1)
	CB(9,11)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,1)/6+RS(2,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)/6)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(2,1)/4-RS(2,2)/4)*BL(1)+RS(2,1)/6+RS(2,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)
     #/6)*SLN(2)*SN(2)**3)-ANT(1)
	CB(9,17)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/6)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)+RS(2,1)/6+RS(2,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/
     #6)*SLN(2)*SN(2)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,1)/6+RS(2,3)
     #/3)*SLN(3)*SN(3)**3)-ANT(1)

	CB(9,3) =
     #((RS(2,2)/2+RS(2,1)/2)*BL(1)*RN(1)**2+(-RS(2,3)/2-RS(2,1)/
     #2)*BL(3)*RN(3)**2+(-RS(2,1)/2-RS(2,2)/2)*BL(1)*SN(1)**2+(RS(2,1)/2
     #+RS(2,3)/2)*BL(3)*SN(3)**2)
	CB(9,9) =
     #((-RS(2,1)/2-RS(2,2)/2)*BL(1)*RN
     #(1)**2+(RS(2,2)/2+RS(2,3)/2)*BL(2)*RN(2)**2+(RS(2,2)/2+RS(2,1)/2)*
     #BL(1)*SN(1)**2+(-RS(2,2)/2-RS(2,3)/2)*BL(2)*SN(2)**2)
	CB(9,15)=
     #((-RS(
     #2,2)/2-RS(2,3)/2)*BL(2)*RN(2)**2+(RS(2,1)/2+RS(2,3)/2)*BL(3)*RN(3)
     #**2+(RS(2,2)/2+RS(2,3)/2)*BL(2)*SN(2)**2+(-RS(2,3)/2-RS(2,1)/2)*BL
     #(3)*SN(3)**2)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRICB4(RS,SLN,RN,SN,BL,BLR,BLS,FLR,FLS,DWR,DWS,CB)
	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cb FOR BENDING STRAINS
C
C	INPUT VARIABLES
C	RS(2,3)			= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)			= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	    = ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	BL(3),BLR,BLS	= SHEAR DEFORMATION PARAMETERS ALONG THE BOUNDARY,
C						R, S, RESP.
C	FLR,FLS,
C	DWR,DWS			= COMMON FACTORS DEFINED IN SUBROUTINE COMFACT
C
C	OUTPUT VARIAVBLES
C	CB(9,18)		= BENDING STRAIN PARAMETER
C	--------------------------------------------------------------
      DIMENSION RS(2,3),SLN(3),RN(3),SN(3),BL(3)
	DIMENSION DWR(9),DWS(9),FLR(3),FLS(3)
	DIMENSION CB(9,18)

C	--------------------------------
C	SOLVE FOR Cb FOR BENDING STRAINS
C	--------------------------------

C	FIRST ROW
C	------------------
	CB(1,4) =
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2)
	CB(1,10)=
     #(-SLN(2)*RN(2)**2*SN(2)*BL(2)/2-SLN(1)*RN(1)**2*SN(1)*BL(1)/2)
	CB(1,16)=
     #(-SLN(3)*RN(3)**2*SN(3)*BL(3)/2-SLN(2)*RN(2)**2*SN(2)*BL(2)/2)

	CB(1,5) =
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))
     #*SLN(1)+(RN(3)**3/2+(1.E0/2.E0-BL(3)/2)*SN(3)**2*RN(3))*SLN(3))
	CB(1,11)=
     #((RN(1)**3/2+(-BL(1)/2+1.E0/2.E0)*SN(1)**2*RN(1))*SLN(1)+(
     #RN(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2))
	CB(1,17)=
     #((RN
     #(2)**3/2+(1.E0/2.E0-BL(2)/2)*SN(2)**2*RN(2))*SLN(2)+(RN(3)**3/2+(1
     #.E0/2.E0-BL(3)/2)*SN(3)**2*RN(3))*SLN(3))

	CB(1,3) =(RN(3)*SN(3)*BL(3)-RN(1)*SN(1)*BL(1))
	CB(1,9) =(-RN(2)*SN(2)*BL(2)+RN(1)*SN(1)*BL(1))
	CB(1,15)=(-RN(3)*SN(3)*BL(3)+RN(2)*SN(2)*BL(2))

C	SECOND ROW
C	------------------
	CB(2,4) =
     #((-RS(1,1)/4-RS(1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,1)
     #/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)-DWR(1)
	CB(2,10)=
     #((-RS(1,1)/4-RS(
     #1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,2)/4-RS(1,3)/4)*BL(2)*S
     #LN(2)*SN(2)*RN(2)**2)-DWR(2)
	CB(2,16)=
     #((-RS(1,2)/4-RS(1,3)/4)*BL(2)*SLN(2)*
     #SN(2)*RN(2)**2+(-RS(1,1)/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)
     #-DWR(3)

	CB(2,5) =
     #((RS(1,2)/6+RS(1,1)/3)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2
     #)/4)*BL(1)+RS(1,2)/6+RS(1,1)/3)*SLN(1)*SN(1)**2*RN(1)+(RS(1,1)/3+R
     #S(1,3)/6)*SLN(3)*RN(3)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/3+
     #RS(1,3)/6)*SLN(3)*SN(3)**2*RN(3))-DWR(4)-FLS(1)
	CB(2,11)=
     #((RS(1,1)/6+RS(1,2)/3)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2
     #)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(1)**2*RN(1)+(RS(1,2)/3+R
     #S(1,3)/6)*SLN(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/3+
     #RS(1,3)/6)*SLN(2)*SN(2)**2*RN(2))-DWR(5)-FLS(2)
	CB(2,17)=
     #((RS(1,2)/6+RS(1,3)/3)*SL
     #N(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/3)*S
     #LN(2)*SN(2)**2*RN(2)+(RS(1,1)/6+RS(1,3)/3)*SLN(3)*RN(3)**3+((-RS(1
     #,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/6+RS(1,3)/3)*SLN(3)*SN(3)**2*RN(3))
     #-DWR(6)-FLS(3)

	CB(2,3) =((-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)*RN(1)+(RS(1,1)/2+RS(1,
     #3)/2)*BL(3)*SN(3)*RN(3))-DWR(7)
	CB(2,9) =((RS(1,1)/2+RS(1,2)/2)*BL(1)*SN(1)*
     #RN(1)+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)*RN(2))-DWR(8)
	CB(2,15)=((RS(1,3)/2+
     #RS(1,2)/2)*BL(2)*SN(2)*RN(2)+(-RS(1,3)/2-RS(1,1)/2)*BL(3)*SN(3)*RN
     #(3))-DWR(9)

C	THIRD ROW
C	------------------
	CB(3,4) =
     #((-RS(2,2)/4-RS(2,1)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,1)
     #/4-RS(2,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)
	CB(3,10)=
     #((-RS(2,2)/4-RS(
     #2,1)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,2)/4-RS(2,3)/4)*BL(2)*S
     #LN(2)*SN(2)*RN(2)**2)
	CB(3,16)=
     #((-RS(2,2)/4-RS(2,3)/4)*BL(2)*SLN(2)*
     #SN(2)*RN(2)**2+(-RS(2,1)/4-RS(2,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)

	CB(3,5) =
     #((RS(2,1)/3+RS(2,2)/6)*SLN(1)*RN(1)**3+((-RS(2,2)/4-RS(2,1
     #)/4)*BL(1)+RS(2,1)/3+RS(2,2)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/6+R
     #S(2,1)/3)*SLN(3)*RN(3)**3+((-RS(2,1)/4-RS(2,3)/4)*BL(3)+RS(2,3)/6+
     #RS(2,1)/3)*SLN(3)*SN(3)**2*RN(3))
 	CB(3,11)=
     #((RS(2,2)/3+RS(2,1)/6)*SLN(1)*RN(1)**3+((-RS(2,2)/4-RS(2,1
     #)/4)*BL(1)+RS(2,2)/3+RS(2,1)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/6+R
     #S(2,2)/3)*SLN(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,3)/6+
     #RS(2,2)/3)*SLN(2)*SN(2)**2*RN(2))
	CB(3,17)=
     #((RS(2,2)/6+RS(2,3)/3)*SL
     #N(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,2)/6+RS(2,3)/3)*S
     #LN(2)*SN(2)**2*RN(2)+(RS(2,1)/6+RS(2,3)/3)*SLN(3)*RN(3)**3+((-RS(2
     #,1)/4-RS(2,3)/4)*BL(3)+RS(2,1)/6+RS(2,3)/3)*SLN(3)*SN(3)**2*RN(3))

	CB(3,3) =
     #((-RS(2,1)/2-RS(2,2)/2)*BL(1)*SN(1)*RN(1)+(RS(2,1)/2+RS(2,
     #3)/2)*BL(3)*SN(3)*RN(3))
	CB(3,9) =
     #((RS(2,1)/2+RS(2,2)/2)*BL(1)*SN(1)*
     #RN(1)+(-RS(2,2)/2-RS(2,3)/2)*BL(2)*SN(2)*RN(2))
	CB(3,15)=
     #((RS(2,2)/2+
     #RS(2,3)/2)*BL(2)*SN(2)*RN(2)+(-RS(2,3)/2-RS(2,1)/2)*BL(3)*SN(3)*RN
     #(3))

C	FOURTH ROW
C	------------------
	CB(4,4) =
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(BL(3)/2-1.E0/2.E0
     #)*SLN(3)*SN(3)*RN(3)**2-SLN(3)*SN(3)**3/2-SLN(1)*SN(1)**3/2)
	CB(4,10)=
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)*RN(1)**2+(BL(2)/2-1.E0/2.E0)*S
     #LN(2)*SN(2)*RN(2)**2-SLN(1)*SN(1)**3/2-SLN(2)*SN(2)**3/2)
	CB(4,16)=
     #(
     #(BL(2)/2-1.E0/2.E0)*SLN(2)*SN(2)*RN(2)**2+(BL(3)/2-1.E0/2.E0)*SLN(
     #3)*SN(3)*RN(3)**2-SLN(3)*SN(3)**3/2-SLN(2)*SN(2)**3/2)

	CB(4,5) =
     #(SLN(1)*SN(1)**2*RN(1)*BL(1)/2+SLN(3)*SN(3)**2*RN(3)*BL(3)/2)
	CB(4,11)=
     #(SLN(1)*SN(1)**2*RN(1)*BL(1)/2+SLN(2)*SN(2)**2*RN(2)*BL(2)/2)
	CB(4,17)=
     #(SLN(2)*SN(2)**2*RN(2)*BL(2)/2+SLN(3)*SN(3)**2*RN(3)*BL(3)/2)

	CB(4,3) =(SN(1)*RN(1)*BL(1)-SN(3)*RN(3)*BL(3))
	CB(4,9) =(-SN(1)*RN(1)*BL(1)+SN(2)*RN(2)*BL(2))
	CB(4,15)=(SN(3)*RN(3)*BL(3)-SN(2)*RN(2)*BL(2))

C	FIFTH ROW
C	------------------
	CB(5,4) =
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1)-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,3)/4+RS(1,1)/4)*BL(3)-RS(1,1)/3-RS(1,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN(1)**3+(-RS(1,1)
     #/3-RS(1,3)/6)*SLN(3)*SN(3)**3)
	CB(5,10)=
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1
     #)-RS(1,2)/3-RS(1,1)/6)*SLN(1)*SN(1)*RN(1)**2+((RS(1,3)/4+RS(1,2)/4
     #)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(1,2)/3-RS(
     #1,1)/6)*SLN(1)*SN(1)**3+(-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)**3)
	CB(5,16)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(1,3)/4+RS(1,1)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(1,3)/3-RS(1,2)/6)*SLN(2)*SN(2)**3+(-RS(1
     #,1)/6-RS(1,3)/3)*SLN(3)*SN(3)**3)

	CB(5,5) =
     #((RS(1,2)/4+RS(1,1)/4)*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,3)/4+RS(1,1)/4)*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))
	CB(5,11)=
     #((RS(1,2)/4+RS(1,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,3
     #)/4+RS(1,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))
	CB(5,17)=
     #((RS(1,3)/4+RS(
     #1,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(1,3)/4+RS(1,1)/4)*BL(3)*SL
     #N(3)*SN(3)**2*RN(3))

	CB(5,3) =
     #((RS(1,1)/2+RS(1,2)/2)*BL(1)*SN(1)*RN(
     #1)+(-RS(1,1)/2-RS(1,3)/2)*BL(3)*SN(3)*RN(3))
	CB(5,9) =
     #((-RS(1,2)/2-RS
     #(1,1)/2)*BL(1)*SN(1)*RN(1)+(RS(1,2)/2+RS(1,3)/2)*BL(2)*SN(2)*RN(2)
     #)
	CB(5,15)=
     #((-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)*RN(2)+(RS(1,3)/2+RS(1,
     #1)/2)*BL(3)*SN(3)*RN(3))

C	SIXTH ROW
C	------------------
	CB(6,4) =
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,1)/4+RS(2,3)/4)*BL(3)-RS(2,3)/6-RS(2,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN(1)**3+(-RS(2,3)
     #/6-RS(2,1)/3)*SLN(3)*SN(3)**3)
     #-DWS(1)-FLR(1)
	CB(6,10)=
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1
     #)-RS(2,2)/3-RS(2,1)/6)*SLN(1)*SN(1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4
     #)*BL(2)-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(2,2)/3-RS(
     #2,1)/6)*SLN(1)*SN(1)**3+(-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)**3)
     #-DWS(2)-FLR(2)
	CB(6,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(2,1)/4+RS(2,3)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN(2)**3+(-RS(2
     #,1)/6-RS(2,3)/3)*SLN(3)*SN(3)**3)
     #-DWS(3)-FLR(3)

	CB(6,5) =
     #((RS(2,2)/4+RS(2,1)/4)*BL
     #(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,1)/4+RS(2,3)/4)*BL(3)*SLN(3)*SN(3)
     #**2*RN(3))
     #-DWS(4)
	CB(6,11)=
     #((RS(2,2)/4+RS(2,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3
     #)/4+RS(2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))
     #-DWS(5)
	CB(6,17)=
     #((RS(2,3)/4+RS(
     #2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(2,1)/4+RS(2,3)/4)*BL(3)*SL
     #N(3)*SN(3)**2*RN(3))
     #-DWS(6)

	CB(6,3) =
     #((RS(2,1)/2+RS(2,2)/2)*BL(1)*SN(1)*RN(
     #1)+(-RS(2,1)/2-RS(2,3)/2)*BL(3)*SN(3)*RN(3))
     #-DWS(7)
	CB(6,9) =
     #((-RS(2,1)/2-RS
     #(2,2)/2)*BL(1)*SN(1)*RN(1)+(RS(2,3)/2+RS(2,2)/2)*BL(2)*SN(2)*RN(2)
     #)
     #-DWS(8)
	CB(6,15)=
     #((-RS(2,3)/2-RS(2,2)/2)*BL(2)*SN(2)*RN(2)+(RS(2,3)/2+RS(2,
     #1)/2)*BL(3)*SN(3)*RN(3))
     #-DWS(9)

C	SEVENTH ROW
C	------------------
	CB(7,4) =
     #((BL(1)/2-1.E0/2.E0)*SLN(1)*RN(1)**3+(-BL(1)/2-1.E0/2.E0)*SLN
     #(1)*SN(1)**2*RN(1)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-1.E0/2.E0
     #-BL(3)/2)*SLN(3)*SN(3)**2*RN(3))
	CB(7,10)=
     #((BL(1)/2-1.E0/2.E0)*SLN(1
     #)*RN(1)**3+(-BL(1)/2-1.E0/2.E0)*SLN(1)*SN(1)**2*RN(1)+(-1.E0/2.E0+
     #BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2.E0)*SLN(2)*SN(2)**2*RN(2
     #))
	CB(7,16)=
     #((-1.E0/2.E0+BL(2)/2)*SLN(2)*RN(2)**3+(-BL(2)/2-1.E0/2.E0)
     #*SLN(2)*SN(2)**2*RN(2)+(BL(3)/2-1.E0/2.E0)*SLN(3)*RN(3)**3+(-1.E0/
     #2.E0-BL(3)/2)*SLN(3)*SN(3)**2*RN(3))

	CB(7,5) =
     #((1.E0/2.E0+BL(1)/2)*S
     #LN(1)*SN(1)*RN(1)**2+(1.E0/2.E0+BL(3)/2)*SLN(3)*SN(3)*RN(3)**2+(-B
     #L(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3+(1.E0/2.E0-BL(3)/2)*SLN(3)*SN(3)
     #**3)
	CB(7,11)=
     #((1.E0/2.E0+BL(1)/2)*SLN(1)*SN(1)*RN(1)**2+(1.E0/2.E0+BL(2
     #)/2)*SLN(2)*SN(2)*RN(2)**2+(-BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3+(-
     #BL(2)/2+1.E0/2.E0)*SLN(2)*SN(2)**3)
	CB(7,17)=
     #((1.E0/2.E0+BL(2)/2)*SL
     #N(2)*SN(2)*RN(2)**2+(1.E0/2.E0+BL(3)/2)*SLN(3)*SN(3)*RN(3)**2+(-BL
     #(2)/2+1.E0/2.E0)*SLN(2)*SN(2)**3+(1.E0/2.E0-BL(3)/2)*SLN(3)*SN(3)*
     #*3)

	CB(7,3) =
     #(-RN(3)**2*BL(3)-SN(1)**2*BL(1)+RN(1)**2*BL(1)+SN(3)**2*BL(3))
	CB(7,9) =
     #(-RN(1)**2*BL(1)+RN(2)**2*BL(2)+SN(1)**2*BL(1)-SN(2)**2*BL(2))
	CB(7,15)=
     #(RN(3)**2*BL(3)+SN(2)**2*BL(2)-RN(2)**2*BL(2)-SN(3)**2*BL(3))

C	EIGHTH ROW
C	------------------
	CB(8,4) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SLN
     #(3)*RN(3)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SL
     #N(3)*SN(3)**2*RN(3))+0.1/3.0
     #-DWS(1)-FLR(1)
	CB(8,10)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN
     #(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SL
     #N(2)*SN(2)**2*RN(2))+0.1/3.0
     #-DWS(2)-FLR(2)
	CB(8,16)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)-RS(1,2)/6-RS(1,3)/3)*SLN(2)*RN(
     #2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)-RS(1,2)/6-RS(1,3)/3)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*SLN
     #(3)*RN(3)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)-RS(1,1)/6-RS(1,3)/3)*SL
     #N(3)*SN(3)**2*RN(3))+0.1/3.0
     #-DWS(3)-FLR(3)

	CB(8,5) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,2)/6+RS(1,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,3)/6+RS(1,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,2)/6+RS(1,1)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)+RS(1,3)/6+RS(1,1)
     #/3)*SLN(3)*SN(3)**3)-DWS(4)
	CB(8,11)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,3)/6+RS(1,2)/3)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,1)/6+RS(1,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,3)/6+RS(1,2)
     #/3)*SLN(2)*SN(2)**3)-DWS(5)
	CB(8,17)=
     #(((RS(1,3)/4+RS(1,2)/4)*BL(2)+RS(1,2)/6+RS(1,3)/3)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,1)/6+RS(1,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/
     #3)*SLN(2)*SN(2)**3+((-RS(1,3)/4-RS(1,1)/4)*BL(3)+RS(1,1)/6+RS(1,3)
     #/3)*SLN(3)*SN(3)**3)-DWS(6)

	CB(8,3) =
     #((RS(1,2)/2+RS(1,1)/2)*BL(1)*RN(1)**2+(-RS(1,1)/2-RS(1,3)/
     #2)*BL(3)*RN(3)**2+(-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)**2+(RS(1,1)/2
     #+RS(1,3)/2)*BL(3)*SN(3)**2)-DWS(7)
	CB(8,9) =
     #((-RS(1,1)/2-RS(1,2)/2)*BL(1)*RN
     #(1)**2+(RS(1,2)/2+RS(1,3)/2)*BL(2)*RN(2)**2+(RS(1,2)/2+RS(1,1)/2)*
     #BL(1)*SN(1)**2+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)**2)-DWS(8)
	CB(8,15)=
     #((-RS(
     #1,2)/2-RS(1,3)/2)*BL(2)*RN(2)**2+(RS(1,1)/2+RS(1,3)/2)*BL(3)*RN(3)
     #**2+(RS(1,2)/2+RS(1,3)/2)*BL(2)*SN(2)**2+(-RS(1,1)/2-RS(1,3)/2)*BL
     #(3)*SN(3)**2)-DWS(9)

C	NINTH ROW
C	------------------
	CB(9,4) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*RN(
     #1)**3+((-RS(2,1)/4-RS(2,2)/4)*BL(1)-RS(2,2)/6-RS(2,1)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SLN
     #(3)*RN(3)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SL
     #N(3)*SN(3)**2*RN(3))-DWR(1)
	CB(9,10)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(2,1)/4-RS(2,2)/4)*BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/3-RS(2,3)/6)*SLN
     #(2)*RN(2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,2)/3-RS(2,3)/6)*SL
     #N(2)*SN(2)**2*RN(2))-DWR(2)
	CB(9,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*RN(
     #2)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*SLN
     #(3)*RN(3)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)-RS(2,1)/6-RS(2,3)/3)*SL
     #N(3)*SN(3)**2*RN(3))-DWR(3)

	CB(9,5) =
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,2)/6+RS(2,1)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)+RS(2,3)/6+RS(2,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,1)/4-RS(2,2)/4)*BL(1)+RS(2,2)/6+RS(2,1)/
     #3)*SLN(1)*SN(1)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,3)/6+RS(2,1)
     #/3)*SLN(3)*SN(3)**3)-DWR(4)-FLS(1)
	CB(9,11)=
     #(((RS(2,1)/4+RS(2,2)/4)*BL(1)+RS(2,1)/6+RS(2,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)/6)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(2,1)/4-RS(2,2)/4)*BL(1)+RS(2,1)/6+RS(2,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,2)/3+RS(2,3)
     #/6)*SLN(2)*SN(2)**3)-DWR(5)-FLS(2)
	CB(9,17)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/6)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)+RS(2,1)/6+RS(2,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(2,3)/4-RS(2,2)/4)*BL(2)+RS(2,3)/3+RS(2,2)/
     #6)*SLN(2)*SN(2)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,1)/6+RS(2,3)
     #/3)*SLN(3)*SN(3)**3)-DWR(6)-FLS(3)

	CB(9,3) =
     #((RS(2,2)/2+RS(2,1)/2)*BL(1)*RN(1)**2+(-RS(2,3)/2-RS(2,1)/
     #2)*BL(3)*RN(3)**2+(-RS(2,1)/2-RS(2,2)/2)*BL(1)*SN(1)**2+(RS(2,1)/2
     #+RS(2,3)/2)*BL(3)*SN(3)**2)-DWR(7)
	CB(9,9) =
     #((-RS(2,1)/2-RS(2,2)/2)*BL(1)*RN
     #(1)**2+(RS(2,2)/2+RS(2,3)/2)*BL(2)*RN(2)**2+(RS(2,2)/2+RS(2,1)/2)*
     #BL(1)*SN(1)**2+(-RS(2,2)/2-RS(2,3)/2)*BL(2)*SN(2)**2)-DWR(8)
	CB(9,15)=
     #((-RS(
     #2,2)/2-RS(2,3)/2)*BL(2)*RN(2)**2+(RS(2,1)/2+RS(2,3)/2)*BL(3)*RN(3)
     #**2+(RS(2,2)/2+RS(2,3)/2)*BL(2)*SN(2)**2+(-RS(2,3)/2-RS(2,1)/2)*BL
     #(3)*SN(3)**2)-DWR(9)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRICB5(ANT,ARS,SL,BL,BLR,BLS,CB)
	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cb FOR BENDING STRAINS
C
C	INPUT VARIABLES
C	ARS(2,3)		= ELEMENT COORDINATES IN LOCAL SPACE
C	SL(3)			= ELEMENT BOUNDARY LENGTHS
C	BL(3),BLR,BLS	= SHEAR DEFORMATION PARAMETERS ALONG THE BOUNDARY,
C						R, S, RESP.
C	FLR,FLS,
C	DWR,DWS			= COMMON FACTORS DEFINED IN SUBROUTINE COMFACT
C
C	LOCAL VARIABLES
C	RS(2,3)			= NATURAL NODAL COORDINATES
C	SLN(3)			= NATURAL SIDE LENGTHS
C	RN(3),SN(3)	    = ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C
C	OUTPUT VARIAVBLES
C	CB(9,18)		= BENDING STRAIN PARAMETER
C	--------------------------------------------------------------
      DIMENSION ANT(3),ARS(2,3),SL(3),BL(3)
	DIMENSION RS(2,3),SLN(3),RN(3),SN(3)
	DIMENSION CB(9,18)

C	---------------------------------------
C	COMPUTE THE DETERMINANT OF THE JACOBIAN
C	---------------------------------------
	AJD = 
     #-ARS(1,1)*ARS(2,3)-ARS(1,2)*ARS(2,1)+ARS(1,2)*ARS(2,3)+ARS(2,
     #1)*ARS(1,3)+ARS(2,2)*ARS(1,1)-ARS(2,2)*ARS(1,3)

C	--------------------------------
C	DEFINE NATURAL NODAL COORDINATES
C	--------------------------------
	RS(1,1) = 0.0
	RS(2,1) = 0.0

	RS(1,2) = 1.0
	RS(2,2) = 0.0

	RS(1,3) = 0.0
	RS(2,3) = 1.0

C	----------------------
C	DEFINE NATURAL NORMALS
C	----------------------
	RN(1) =  0.0
	RN(2) =  DSQRT(1.0D0/2.0D0)
	RN(3) = -1.0

	SN(1) = -1.0
	SN(2) =  DSQRT(1.0D0/2.0D0)
	SN(3) =  0.0

C	---------------------------
C	DEFINE NATURAL SIDE LENGTHS
C	---------------------------
	SLN(1) = 1.0
	SLN(2) = DSQRT(2.0D0)
	SLN(3) = 1.0

C	--------------------------------
C	SOLVE FOR Cb FOR BENDING STRAINS
C	--------------------------------

C	FIRST ROW
C	------------------
	CB(1,4) =
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(3)*RN(3)**2*SN(3)*BL(3)/2
     #)*AJD
	CB(1,10)=
     #(-SLN(1)*RN(1)**2*SN(1)*BL(1)/2-SLN(2)*RN(2)**2*SN(2)
     #*BL(2)/2)*AJD
	CB(1,16)=
     #(-SLN(2)*RN(2)**2*SN(2)*BL(2)/2-SLN(3)*RN(3)*
     #*2*SN(3)*BL(3)/2)*AJD

	CB(1,5) =
     #(SLN(1)*RN(1)**3/2+SLN(1)*RN(1)*SN(1)
     #**2/2+SLN(3)*RN(3)**3/2-SLN(3)*RN(3)*SN(3)**2*BL(3)/2+SLN(3)*RN(3)
     #*SN(3)**2/2-SLN(1)*RN(1)*SN(1)**2*BL(1)/2)*AJD
	CB(1,11)=
     #(SLN(1)*RN(1)*SN(1)**2/2-SLN(1)*RN(1)*SN(1)**2*BL(1)/2+SLN
     #(2)*RN(2)*SN(2)**2/2-SLN(2)*RN(2)*SN(2)**2*BL(2)/2+SLN(2)*RN(2)**3
     #/2+SLN(1)*RN(1)**3/2)*AJD
	CB(1,17)=
     #(-SLN(2)*RN(2)*SN(2)**2*BL(2)/2+S
     #LN(2)*RN(2)**3/2+SLN(3)*RN(3)**3/2-SLN(3)*RN(3)*SN(3)**2*BL(3)/2+S
     #LN(2)*RN(2)*SN(2)**2/2+SLN(3)*RN(3)*SN(3)**2/2)*AJD

	CB(1,3) =(RN(3)*SN(3)*BL(3)-RN(1)*SN(1)*BL(1))*AJD
	CB(1,9) =(RN(1)*SN(1)*BL(1)-RN(2)*SN(2)*BL(2))*AJD
	CB(1,15)=(-RN(3)*SN(3)*BL(3)+RN(2)*SN(2)*BL(2))*AJD

C	SECOND ROW
C	------------------
	CB(2,4) =
     #((-RS(1,1)/4-RS(1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,1)
     #/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)*AJD
     	CB(2,10)=
     #((-RS(1,1)/4
     #-RS(1,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(1,2)/4-RS(1,3)/4)*BL(
     #2)*SLN(2)*SN(2)*RN(2)**2)*AJD
     	CB(2,16)=
     #((-RS(1,2)/4-RS(1,3)/4)*BL(2)*SLN(2)*SN(2)*RN(2)**2+(-RS(1
     #,1)/4-RS(1,3)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)*AJD

	CB(2,5) =
     #((RS(1,1)
     #/3+RS(1,2)/6)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2)/4)*BL(1)+RS(1,1
     #)/3+RS(1,2)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(1,1)/3+RS(1,3)/6)*SLN(3)*
     #RN(3)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,1)/3+RS(1,3)/6)*SLN(3)
     #*SN(3)**2*RN(3))*AJD
     #-ANT(1)
	CB(2,11)=
     #((RS(1,2)/3+RS(1,1)/6)*SLN(1)*RN(1)**3+((-RS(1,1)/4-RS(1,2
     #)/4)*BL(1)+RS(1,2)/3+RS(1,1)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(1,2)/3+R
     #S(1,3)/6)*SLN(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/3+
     #RS(1,3)/6)*SLN(2)*SN(2)**2*RN(2))*AJD
     #-ANT(1)
	CB(2,17)=
     #((RS(1,2)/6+RS(1,3)/3
     #)*SLN(2)*RN(2)**3+((-RS(1,2)/4-RS(1,3)/4)*BL(2)+RS(1,2)/6+RS(1,3)/
     #3)*SLN(2)*SN(2)**2*RN(2)+(RS(1,3)/3+RS(1,1)/6)*SLN(3)*RN(3)**3+((-
     #RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,3)/3+RS(1,1)/6)*SLN(3)*SN(3)**2*RN
     #(3))*AJD
     #-ANT(1)

	CB(2,3) =
     #((-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)*RN(1)+(RS(1,1)/2+RS(1,
     #3)/2)*BL(3)*SN(3)*RN(3))*AJD
	CB(2,9) =
     #((RS(1,1)/2+RS(1,2)/2)*BL(1)*SN
     #(1)*RN(1)+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)*RN(2))*AJD
	CB(2,15)=
     #((RS
     #(1,2)/2+RS(1,3)/2)*BL(2)*SN(2)*RN(2)+(-RS(1,1)/2-RS(1,3)/2)*BL(3)*
     #SN(3)*RN(3))*AJD

C	THIRD ROW
C	------------------
	CB(3,4) =
     #((-RS(2,1)/4-RS(2,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,3)
     #/4-RS(2,1)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)*AJD
	CB(3,10)=
     #((-RS(2,1)/4
     #-RS(2,2)/4)*BL(1)*SLN(1)*SN(1)*RN(1)**2+(-RS(2,2)/4-RS(2,3)/4)*BL(
     #2)*SLN(2)*SN(2)*RN(2)**2)*AJD
	CB(3,16)=
     #((-RS(2,2)/4-RS(2,3)/4)*BL(2)*SLN(2)*SN(2)*RN(2)**2+(-RS(2
     #,3)/4-RS(2,1)/4)*BL(3)*SLN(3)*SN(3)*RN(3)**2)*AJD

	CB(3,5) =
     #((RS(2,2)
     #/6+RS(2,1)/3)*SLN(1)*RN(1)**3+((-RS(2,1)/4-RS(2,2)/4)*BL(1)+RS(2,2
     #)/6+RS(2,1)/3)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/6+RS(2,1)/3)*SLN(3)*
     #RN(3)**3+((-RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,3)/6+RS(2,1)/3)*SLN(3)
     #*SN(3)**2*RN(3))*AJD
 	CB(3,11)=
     #((RS(2,2)/3+RS(2,1)/6)*SLN(1)*RN(1)**3+((-RS(2,1)/4-RS(2,2
     #)/4)*BL(1)+RS(2,2)/3+RS(2,1)/6)*SLN(1)*SN(1)**2*RN(1)+(RS(2,2)/3+R
     #S(2,3)/6)*SLN(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,2)/3+
     #RS(2,3)/6)*SLN(2)*SN(2)**2*RN(2))*AJD
	CB(3,17)=
     #((RS(2,3)/3+RS(2,2)/6
     #)*SLN(2)*RN(2)**3+((-RS(2,2)/4-RS(2,3)/4)*BL(2)+RS(2,3)/3+RS(2,2)/
     #6)*SLN(2)*SN(2)**2*RN(2)+(RS(2,1)/6+RS(2,3)/3)*SLN(3)*RN(3)**3+((-
     #RS(2,3)/4-RS(2,1)/4)*BL(3)+RS(2,1)/6+RS(2,3)/3)*SLN(3)*SN(3)**2*RN
     #(3))*AJD

	CB(3,3) =
     #((-RS(2,2)/2-RS(2,1)/2)*BL(1)*SN(1)*RN(1)+(RS(2,1)/2+RS(2,
     #3)/2)*BL(3)*SN(3)*RN(3))*AJD
	CB(3,9) =
     #((RS(2,2)/2+RS(2,1)/2)*BL(1)*SN
     #(1)*RN(1)+(-RS(2,2)/2-RS(2,3)/2)*BL(2)*SN(2)*RN(2))*AJD
	CB(3,15)=
     #((RS
     #(2,3)/2+RS(2,2)/2)*BL(2)*SN(2)*RN(2)+(-RS(2,3)/2-RS(2,1)/2)*BL(3)*
     #SN(3)*RN(3))*AJD


C	FOURTH ROW
C	------------------
	CB(4,4) =
     #((-1.E0/2.E0+BL(1)/2)*SLN(1)*SN(1)*RN(1)**2+(-1.E0/2.E0+BL(3)
     #/2)*SLN(3)*SN(3)*RN(3)**2-SLN(1)*SN(1)**3/2-SLN(3)*SN(3)**3/2)*AJD
	CB(4,10)=
     #((-1.E0/2.E0+BL(1)/2)*SLN(1)*SN(1)*RN(1)**2+(-1.E0/2.E0+BL
     #(2)/2)*SLN(2)*SN(2)*RN(2)**2-SLN(1)*SN(1)**3/2-SLN(2)*SN(2)**3/2)*
     #AJD
	CB(4,16)=
     #((-1.E0/2.E0+BL(2)/2)*SLN(2)*SN(2)*RN(2)**2+(-1.E0/2.E0
     #+BL(3)/2)*SLN(3)*SN(3)*RN(3)**2-SLN(2)*SN(2)**3/2-SLN(3)*SN(3)**3/
     #2)*AJD

	CB(4,5) =
     #(SLN(1)*SN(1)**2*RN(1)*BL(1)/2+SLN(3)*SN(3)**2*RN(3)
     #*BL(3)/2)*AJD
	CB(4,11)=
     #(SLN(2)*SN(2)**2*RN(2)*BL(2)/2+SLN(1)*SN(1)**2*RN(1)*BL(1)
     #/2)*AJD
	CB(4,17)=
     #(SLN(3)*SN(3)**2*RN(3)*BL(3)/2+SLN(2)*SN(2)**2*RN(2
     #)*BL(2)/2)*AJD

	CB(4,3) =(SN(1)*RN(1)*BL(1)-SN(3)*RN(3)*BL(3))*AJD
	CB(4,9) =(-SN(1)*RN(1)*BL(1)+SN(2)*RN(2)*BL(2))*AJD
	CB(4,15)=(-SN(2)*RN(2)*BL(2)+SN(3)*RN(3)*BL(3))*AJD


C	FIFTH ROW
C	------------------
	CB(5,4) =
     #(((RS(1,2)/4+RS(1,1)/4)*BL(1)-RS(1,1)/3-RS(1,2)/6)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,1)/3-RS(1,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(1,1)/3-RS(1,2)/6)*SLN(1)*SN(1)**3+(-RS(1,1)
     #/3-RS(1,3)/6)*SLN(3)*SN(3)**3)*AJD
	CB(5,10)=
     #(((RS(1,2)/4+RS(1,1)/4)*
     #BL(1)-RS(1,2)/3-RS(1,1)/6)*SLN(1)*SN(1)*RN(1)**2+((RS(1,2)/4+RS(1,
     #3)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(1,2)/3
     #-RS(1,1)/6)*SLN(1)*SN(1)**3+(-RS(1,3)/6-RS(1,2)/3)*SLN(2)*SN(2)**3
     #)*AJD
	CB(5,16)=
     #(((RS(1,2)/4+RS(1,3)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,3)/3-RS(1,1)/6)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(1,3)/3-RS(1,2)/6)*SLN(2)*SN(2)**3+(-RS(1
     #,3)/3-RS(1,1)/6)*SLN(3)*SN(3)**3)*AJD

	CB(5,5) =
     #((RS(1,2)/4+RS(1,1)/4
     #)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,1)/4+RS(1,3)/4)*BL(3)*SLN(3)*S
     #N(3)**2*RN(3))*AJD
	CB(5,11)=
     #((RS(1,2)/4+RS(1,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(1,2
     #)/4+RS(1,3)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))*AJD
	CB(5,17)=
     #((RS(1,2)/4
     #+RS(1,3)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(1,1)/4+RS(1,3)/4)*BL(3
     #)*SLN(3)*SN(3)**2*RN(3))*AJD

	CB(5,3) =
     #((RS(1,2)/2+RS(1,1)/2)*BL(1)*S
     #N(1)*RN(1)+(-RS(1,3)/2-RS(1,1)/2)*BL(3)*SN(3)*RN(3))*AJD
	CB(5,9) =
     #((-
     #RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)*RN(1)+(RS(1,3)/2+RS(1,2)/2)*BL(2)
     #*SN(2)*RN(2))*AJD
	CB(5,15)=
     #((-RS(1,3)/2-RS(1,2)/2)*BL(2)*SN(2)*RN(2)+
     #(RS(1,1)/2+RS(1,3)/2)*BL(3)*SN(3)*RN(3))*AJD

C	SIXTH ROW
C	------------------
	CB(6,4) =
     #(((RS(2,2)/4+RS(2,1)/4)*BL(1)-RS(2,1)/3-RS(2,2)/6)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,1)/3-RS(2,3)/6)*SLN(
     #3)*SN(3)*RN(3)**2+(-RS(2,1)/3-RS(2,2)/6)*SLN(1)*SN(1)**3+(-RS(2,1)
     #/3-RS(2,3)/6)*SLN(3)*SN(3)**3)*AJD
     #+ANT(1)
	CB(6,10)=
     #(((RS(2,2)/4+RS(2,1)/4)*
     #BL(1)-RS(2,1)/6-RS(2,2)/3)*SLN(1)*SN(1)*RN(1)**2+((RS(2,3)/4+RS(2,
     #2)/4)*BL(2)-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)*RN(2)**2+(-RS(2,1)/6
     #-RS(2,2)/3)*SLN(1)*SN(1)**3+(-RS(2,3)/6-RS(2,2)/3)*SLN(2)*SN(2)**3
     #)*AJD
     #+ANT(1)
	CB(6,16)=
     #(((RS(2,3)/4+RS(2,2)/4)*BL(2)-RS(2,2)/6-RS(2,3)/3)*SLN(2)*
     #SN(2)*RN(2)**2+((RS(2,3)/4+RS(2,1)/4)*BL(3)-RS(2,3)/3-RS(2,1)/6)*S
     #LN(3)*SN(3)*RN(3)**2+(-RS(2,2)/6-RS(2,3)/3)*SLN(2)*SN(2)**3+(-RS(2
     #,3)/3-RS(2,1)/6)*SLN(3)*SN(3)**3)*AJD
     #+ANT(1)

	CB(6,5) =
     #((RS(2,2)/4+RS(2,1)/4
     #)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3)/4+RS(2,1)/4)*BL(3)*SLN(3)*S
     #N(3)**2*RN(3))*AJD
	CB(6,11)=
     #((RS(2,2)/4+RS(2,1)/4)*BL(1)*SLN(1)*SN(1)**2*RN(1)+(RS(2,3
     #)/4+RS(2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2))*AJD
	CB(6,17)=
     #((RS(2,3)/4
     #+RS(2,2)/4)*BL(2)*SLN(2)*SN(2)**2*RN(2)+(RS(2,3)/4+RS(2,1)/4)*BL(3
     #)*SLN(3)*SN(3)**2*RN(3))*AJD

	CB(6,3) =
     #((RS(2,2)/2+RS(2,1)/2)*BL(1)*S
     #N(1)*RN(1)+(-RS(2,1)/2-RS(2,3)/2)*BL(3)*SN(3)*RN(3))*AJD
	CB(6,9) =
     #((-
     #RS(2,2)/2-RS(2,1)/2)*BL(1)*SN(1)*RN(1)+(RS(2,2)/2+RS(2,3)/2)*BL(2)
     #*SN(2)*RN(2))*AJD
	CB(6,15)=
     #((-RS(2,3)/2-RS(2,2)/2)*BL(2)*SN(2)*RN(2)+
     #(RS(2,1)/2+RS(2,3)/2)*BL(3)*SN(3)*RN(3))*AJD

C	SEVENTH ROW
C	------------------
	CB(7,4) =
     #((-1.E0/2.E0+BL(1)/2)*SLN(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SL
     #N(1)*SN(1)**2*RN(1)+(-1.E0/2.E0+BL(3)/2)*SLN(3)*RN(3)**3+(-BL(3)/2
     #-1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))*AJD
	CB(7,10)=
     #((-1.E0/2.E0+BL(1)/2
     #)*SLN(1)*RN(1)**3+(-1.E0/2.E0-BL(1)/2)*SLN(1)*SN(1)**2*RN(1)+(-1.E
     #0/2.E0+BL(2)/2)*SLN(2)*RN(2)**3+(-1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)*
     #*2*RN(2))*AJD
	CB(7,16)=
     #((-1.E0/2.E0+BL(2)/2)*SLN(2)*RN(2)**3+(-1.E0/2.E0-BL(2)/2)
     #*SLN(2)*SN(2)**2*RN(2)+(-1.E0/2.E0+BL(3)/2)*SLN(3)*RN(3)**3+(-BL(3
     #)/2-1.E0/2.E0)*SLN(3)*SN(3)**2*RN(3))*AJD

	CB(7,5) =
     #((1.E0/2.E0+BL(1)
     #/2)*SLN(1)*SN(1)*RN(1)**2+(1.E0/2.E0+BL(3)/2)*SLN(3)*SN(3)*RN(3)**
     #2+(1.E0/2.E0-BL(1)/2)*SLN(1)*SN(1)**3+(-BL(3)/2+1.E0/2.E0)*SLN(3)*
     #SN(3)**3)*AJD
	CB(7,11)=
     #((1.E0/2.E0+BL(1)/2)*SLN(1)*SN(1)*RN(1)**2+(1.E0/2.E0+BL(2
     #)/2)*SLN(2)*SN(2)*RN(2)**2+(1.E0/2.E0-BL(1)/2)*SLN(1)*SN(1)**3+(-B
     #L(2)/2+1.E0/2.E0)*SLN(2)*SN(2)**3)*AJD
	CB(7,17)=
     #((1.E0/2.E0+BL(2)/2)
     #*SLN(2)*SN(2)*RN(2)**2+(1.E0/2.E0+BL(3)/2)*SLN(3)*SN(3)*RN(3)**2+(
     #-BL(2)/2+1.E0/2.E0)*SLN(2)*SN(2)**3+(-BL(3)/2+1.E0/2.E0)*SLN(3)*SN
     #(3)**3)*AJD

	CB(7,3) =
     #(RN(1)**2*BL(1)-RN(3)**2*BL(3)-SN(1)**2*BL(1)+S
     #N(3)**2*BL(3))*AJD
	CB(7,9) =
     #(-SN(2)**2*BL(2)+SN(1)**2*BL(1)-RN(1)**2*
     #BL(1)+RN(2)**2*BL(2))*AJD
	CB(7,15)=
     #(-SN(3)**2*BL(3)-RN(2)**2*BL(2)+RN
     #(3)**2*BL(3)+SN(2)**2*BL(2))*AJD

C	EIGHTH ROW
C	------------------
	CB(8,4) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,2)/6-RS(1,1)/3)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,2)/6-RS(1,1)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SLN
     #(3)*RN(3)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)-RS(1,3)/6-RS(1,1)/3)*SL
     #N(3)*SN(3)**2*RN(3))*AJD
     #+ANT(1)
	CB(8,10)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*RN(
     #1)**3+((-RS(1,2)/4-RS(1,1)/4)*BL(1)-RS(1,1)/6-RS(1,2)/3)*SLN(1)*SN
     #(1)**2*RN(1)+((RS(1,2)/4+RS(1,3)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SLN
     #(2)*RN(2)**3+((-RS(1,3)/4-RS(1,2)/4)*BL(2)-RS(1,3)/6-RS(1,2)/3)*SL
     #N(2)*SN(2)**2*RN(2))*AJD
     #+ANT(1)
	CB(8,16)=
     #(((RS(1,2)/4+RS(1,3)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*RN(
     #2)**3+((-RS(1,3)/4-RS(1,2)/4)*BL(2)-RS(1,3)/3-RS(1,2)/6)*SLN(2)*SN
     #(2)**2*RN(2)+((RS(1,1)/4+RS(1,3)/4)*BL(3)-RS(1,3)/3-RS(1,1)/6)*SLN
     #(3)*RN(3)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)-RS(1,3)/3-RS(1,1)/6)*SL
     #N(3)*SN(3)**2*RN(3))*AJD
     #+ANT(1)

	CB(8,5) =
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,1)/3+RS(1,2)/6)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,3)/6+RS(1,1)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,1)/3+RS(1,2)/
     #6)*SLN(1)*SN(1)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,3)/6+RS(1,1)
     #/3)*SLN(3)*SN(3)**3)*AJD
	CB(8,11)=
     #(((RS(1,1)/4+RS(1,2)/4)*BL(1)+RS(1,1)/6+RS(1,2)/3)*SLN(1)*SN(
     #1)*RN(1)**2+((RS(1,2)/4+RS(1,3)/4)*BL(2)+RS(1,2)/3+RS(1,3)/6)*SLN(
     #2)*SN(2)*RN(2)**2+((-RS(1,2)/4-RS(1,1)/4)*BL(1)+RS(1,1)/6+RS(1,2)/
     #3)*SLN(1)*SN(1)**3+((-RS(1,3)/4-RS(1,2)/4)*BL(2)+RS(1,2)/3+RS(1,3)
     #/6)*SLN(2)*SN(2)**3)*AJD
	CB(8,17)=
     #(((RS(1,2)/4+RS(1,3)/4)*BL(2)+RS(1,3)/3+RS(1,2)/6)*SLN(2)*SN(
     #2)*RN(2)**2+((RS(1,1)/4+RS(1,3)/4)*BL(3)+RS(1,1)/6+RS(1,3)/3)*SLN(
     #3)*SN(3)*RN(3)**2+((-RS(1,3)/4-RS(1,2)/4)*BL(2)+RS(1,3)/3+RS(1,2)/
     #6)*SLN(2)*SN(2)**3+((-RS(1,1)/4-RS(1,3)/4)*BL(3)+RS(1,3)/3+RS(1,1)
     #/6)*SLN(3)*SN(3)**3)*AJD

	CB(8,3) =
     #((RS(1,2)/2+RS(1,1)/2)*BL(1)*RN(1)**2+(-RS(1,3)/2-RS(1,1)/
     #2)*BL(3)*RN(3)**2+(-RS(1,1)/2-RS(1,2)/2)*BL(1)*SN(1)**2+(RS(1,3)/2
     #+RS(1,1)/2)*BL(3)*SN(3)**2)*AJD
	CB(8,9) =
     #((-RS(1,1)/2-RS(1,2)/2)*BL(1
     #)*RN(1)**2+(RS(1,3)/2+RS(1,2)/2)*BL(2)*RN(2)**2+(RS(1,2)/2+RS(1,1)
     #/2)*BL(1)*SN(1)**2+(-RS(1,2)/2-RS(1,3)/2)*BL(2)*SN(2)**2)*AJD
	CB(8,15)=
     #((-RS(1,2)/2-RS(1,3)/2)*BL(2)*RN(2)**2+(RS(1,3)/2+RS(1,1)/2)*BL(
     #3)*RN(3)**2+(RS(1,3)/2+RS(1,2)/2)*BL(2)*SN(2)**2+(-RS(1,3)/2-RS(1,
     #1)/2)*BL(3)*SN(3)**2)*AJD

C	NINTH ROW
C	------------------
	CB(9,4) =
     #(((AJD*RS(2,1)/4+AJD*RS(2,2)/4)*BL(1)-AJD*RS(2,2)/6-AJD*RS(2,
     #1)/3)*SLN(1)*RN(1)**3+((-AJD*RS(2,2)/4-AJD*RS(2,1)/4)*BL(1)-AJD*RS
     #(2,2)/6-AJD*RS(2,1)/3)*SLN(1)*SN(1)**2*RN(1)+((AJD*RS(2,3)/4+AJD*R
     #S(2,1)/4)*BL(3)-AJD*RS(2,1)/3-AJD*RS(2,3)/6)*SLN(3)*RN(3)**3+((-AJ
     #D*RS(2,1)/4-AJD*RS(2,3)/4)*BL(3)-AJD*RS(2,1)/3-AJD*RS(2,3)/6)*SLN(
     #3)*SN(3)**2*RN(3))
	CB(9,10)=
     #(((AJD*RS(2,1)/4+AJD*RS(2,2)/4)*BL(1)-AJD*RS(2,2)/3-AJD*RS(2,
     #1)/6)*SLN(1)*RN(1)**3+((-AJD*RS(2,2)/4-AJD*RS(2,1)/4)*BL(1)-AJD*RS
     #(2,2)/3-AJD*RS(2,1)/6)*SLN(1)*SN(1)**2*RN(1)+((AJD*RS(2,3)/4+AJD*R
     #S(2,2)/4)*BL(2)-AJD*RS(2,2)/3-AJD*RS(2,3)/6)*SLN(2)*RN(2)**3+((-AJ
     #D*RS(2,2)/4-AJD*RS(2,3)/4)*BL(2)-AJD*RS(2,2)/3-AJD*RS(2,3)/6)*SLN(
     #2)*SN(2)**2*RN(2))
	CB(9,16)=
     #(((AJD*RS(2,3)/4+AJD*RS(2,2)/4)*BL(2)-AJD*RS(2,3)/3-AJD*RS(2,
     #2)/6)*SLN(2)*RN(2)**3+((-AJD*RS(2,2)/4-AJD*RS(2,3)/4)*BL(2)-AJD*RS
     #(2,3)/3-AJD*RS(2,2)/6)*SLN(2)*SN(2)**2*RN(2)+((AJD*RS(2,3)/4+AJD*R
     #S(2,1)/4)*BL(3)-AJD*RS(2,1)/6-AJD*RS(2,3)/3)*SLN(3)*RN(3)**3+((-AJ
     #D*RS(2,1)/4-AJD*RS(2,3)/4)*BL(3)-AJD*RS(2,1)/6-AJD*RS(2,3)/3)*SLN(
     #3)*SN(3)**2*RN(3))

	CB(9,5) =
     #(((AJD*RS(2,1)/4+AJD*RS(2,2)/4)*BL(1)+AJD*RS(2,2)/6+AJD*RS(2,
     #1)/3)*SLN(1)*SN(1)*RN(1)**2+((AJD*RS(2,3)/4+AJD*RS(2,1)/4)*BL(3)+A
     #JD*RS(2,1)/3+AJD*RS(2,3)/6)*SLN(3)*SN(3)*RN(3)**2+((-AJD*RS(2,2)/4
     #-AJD*RS(2,1)/4)*BL(1)+AJD*RS(2,2)/6+AJD*RS(2,1)/3)*SLN(1)*SN(1)**3
     #+((-AJD*RS(2,1)/4-AJD*RS(2,3)/4)*BL(3)+AJD*RS(2,1)/3+AJD*RS(2,3)/6
     #)*SLN(3)*SN(3)**3)
     #-ANT(1)
	CB(9,11)=
     #(((AJD*RS(2,1)/4+AJD*RS(2,2)/4)*BL(1)+AJD*RS(2,1)/6+AJD*RS(2,
     #2)/3)*SLN(1)*SN(1)*RN(1)**2+((AJD*RS(2,3)/4+AJD*RS(2,2)/4)*BL(2)+A
     #JD*RS(2,2)/3+AJD*RS(2,3)/6)*SLN(2)*SN(2)*RN(2)**2+((-AJD*RS(2,2)/4
     #-AJD*RS(2,1)/4)*BL(1)+AJD*RS(2,1)/6+AJD*RS(2,2)/3)*SLN(1)*SN(1)**3
     #+((-AJD*RS(2,2)/4-AJD*RS(2,3)/4)*BL(2)+AJD*RS(2,2)/3+AJD*RS(2,3)/6
     #)*SLN(2)*SN(2)**3)
     #-ANT(1)
	CB(9,17)=
     #(((AJD*RS(2,3)/4+AJD*RS(2,2)/4)*BL(2)+AJD*RS(2,3)/3+AJD*RS(2,
     #2)/6)*SLN(2)*SN(2)*RN(2)**2+((AJD*RS(2,3)/4+AJD*RS(2,1)/4)*BL(3)+A
     #JD*RS(2,3)/3+AJD*RS(2,1)/6)*SLN(3)*SN(3)*RN(3)**2+((-AJD*RS(2,2)/4
     #-AJD*RS(2,3)/4)*BL(2)+AJD*RS(2,3)/3+AJD*RS(2,2)/6)*SLN(2)*SN(2)**3
     #+((-AJD*RS(2,1)/4-AJD*RS(2,3)/4)*BL(3)+AJD*RS(2,3)/3+AJD*RS(2,1)/6
     #)*SLN(3)*SN(3)**3)
     #-ANT(1)

	CB(9,3) =
     #((AJD*RS(2,2)/2+AJD*RS(2,1)/2)*BL(1)*RN(1)**2+(-AJD*RS(2,1
     #)/2-AJD*RS(2,3)/2)*BL(3)*RN(3)**2+(-AJD*RS(2,1)/2-AJD*RS(2,2)/2)*B
     #L(1)*SN(1)**2+(AJD*RS(2,1)/2+AJD*RS(2,3)/2)*BL(3)*SN(3)**2)
	CB(9,9) =
     #((-AJD*RS(2,1)/2-AJD*RS(2,2)/2)*BL(1)*RN(1)**2+(AJD*RS(2,2
     #)/2+AJD*RS(2,3)/2)*BL(2)*RN(2)**2+(AJD*RS(2,2)/2+AJD*RS(2,1)/2)*BL
     #(1)*SN(1)**2+(-AJD*RS(2,3)/2-AJD*RS(2,2)/2)*BL(2)*SN(2)**2)
	CB(9,15)=
     #((-AJD*RS(2,3)/2-AJD*RS(2,2)/2)*BL(2)*RN(2)**2+(AJD*RS(2,1)/2+AJD*
     #RS(2,3)/2)*BL(3)*RN(3)**2+(AJD*RS(2,2)/2+AJD*RS(2,3)/2)*BL(2)*SN(2
     #)**2+(-AJD*RS(2,1)/2-AJD*RS(2,3)/2)*BL(3)*SN(3)**2)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRIAB(RS,SLN,RN,SN,AB1)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR STRAIN PARAMETERS AB1
C
C	INPUT VARIABLES
C	RS(2,3)		= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)		= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C
C	OUTPUT VARIAVBLES
C	AB1(3,3)	= STRAIN PARAMETERS
C	-------------------------------------------------------
      DIMENSION RS(2,3),SLN(3),RN(3),SN(3)
	DIMENSION AB1(3,3)

C	-------------
C	SOLVE FOR AB1
C	-------------
C	UPPER TRIANGLE
      AB1(1,1) = SLN(1)*RN(1)*RS(1,1)/2+SLN(1)*RN(1)*RS(1,2)/2+SLN(2)*RN
     #(2)*RS(1,2)/2+SLN(2)*RN(2)*RS(1,3)/2+SLN(3)*RN(3)*RS(1,3)/2+SLN(3)
     #*RN(3)*RS(1,1)/2

      AB1(1,2) = SLN(1)*RN(1)*RS(1,1)**2/6+SLN(1)*RN(1)*RS(1,1)*RS(1,2)/
     #6+SLN(1)*RN(1)*RS(1,2)**2/6+SLN(2)*RN(2)*RS(1,2)**2/6+SLN(2)*RN(2)
     #*RS(1,2)*RS(1,3)/6+SLN(2)*RN(2)*RS(1,3)**2/6+SLN(3)*RN(3)*RS(1,3)*
     #*2/6+SLN(3)*RN(3)*RS(1,3)*RS(1,1)/6+SLN(3)*RN(3)*RS(1,1)**2/6

      AB1(1,3) = SLN(1)*SN(1)*RS(2,1)**2/6+SLN(1)*SN(1)*RS(2,1)*RS(2,2)/
     #6+SLN(1)*SN(1)*RS(2,2)**2/6+SLN(2)*SN(2)*RS(2,2)**2/6+SLN(2)*SN(2)
     #*RS(2,2)*RS(2,3)/6+SLN(2)*SN(2)*RS(2,3)**2/6+SLN(3)*SN(3)*RS(2,3)*
     #*2/6+SLN(3)*SN(3)*RS(2,3)*RS(2,1)/6+SLN(3)*SN(3)*RS(2,1)**2/6

      AB1(2,2) = SLN(1)*RN(1)*RS(1,1)**3/12+SLN(1)*RN(1)*RS(1,1)**2*RS(1
     #,2)/12+SLN(1)*RN(1)*RS(1,1)*RS(1,2)**2/12+SLN(1)*RN(1)*RS(1,2)**3/
     #12+SLN(2)*RN(2)*RS(1,2)**3/12+SLN(2)*RN(2)*RS(1,2)**2*RS(1,3)/12+S
     #LN(2)*RN(2)*RS(1,2)*RS(1,3)**2/12+SLN(2)*RN(2)*RS(1,3)**3/12+SLN(3
     #)*RN(3)*RS(1,3)**3/12+SLN(3)*RN(3)*RS(1,3)**2*RS(1,1)/12+SLN(3)*RN
     #(3)*RS(1,3)*RS(1,1)**2/12+SLN(3)*RN(3)*RS(1,1)**3/12

      s1 = SLN(1)*RN(1)*RS(1,1)*RS(1,2)*RS(2,1)/12+SLN(1)*RN(1)*RS(1,1)*
     #RS(1,2)*RS(2,2)/12+SLN(2)*RN(2)*RS(1,2)*RS(1,3)*RS(2,2)/12+SLN(2)*
     #RN(2)*RS(1,2)*RS(1,3)*RS(2,3)/12+SLN(3)*RN(3)*RS(1,3)*RS(1,1)*RS(2
     #,3)/12+SLN(3)*RN(3)*RS(1,3)*RS(1,1)*RS(2,1)/12+SLN(1)*RN(1)*RS(1,1
     #)**2*RS(2,1)/8+SLN(1)*RN(1)*RS(1,1)**2*RS(2,2)/24+SLN(1)*RN(1)*RS(
     #1,2)**2*RS(2,1)/24

      AB1(2,3) = s1+SLN(1)*RN(1)*RS(1,2)**2*RS(2,2)/8+SLN(2)*RN(2)*RS(1,
     #2)**2*RS(2,2)/8+SLN(2)*RN(2)*RS(1,2)**2*RS(2,3)/24+SLN(2)*RN(2)*RS
     #(1,3)**2*RS(2,2)/24+SLN(2)*RN(2)*RS(1,3)**2*RS(2,3)/8+SLN(3)*RN(3)
     #*RS(1,3)**2*RS(2,3)/8+SLN(3)*RN(3)*RS(1,3)**2*RS(2,1)/24+SLN(3)*RN
     #(3)*RS(1,1)**2*RS(2,3)/24+SLN(3)*RN(3)*RS(1,1)**2*RS(2,1)/8

      AB1(3,3) = SLN(1)*SN(1)*RS(2,1)**3/12+SLN(1)*SN(1)*RS(2,1)**2*RS(2
     #,2)/12+SLN(1)*SN(1)*RS(2,1)*RS(2,2)**2/12+SLN(1)*SN(1)*RS(2,2)**3/
     #12+SLN(2)*SN(2)*RS(2,2)**3/12+SLN(2)*SN(2)*RS(2,2)**2*RS(2,3)/12+S
     #LN(2)*SN(2)*RS(2,2)*RS(2,3)**2/12+SLN(2)*SN(2)*RS(2,3)**3/12+SLN(3
     #)*SN(3)*RS(2,3)**3/12+SLN(3)*SN(3)*RS(2,3)**2*RS(2,1)/12+SLN(3)*SN
     #(3)*RS(2,3)*RS(2,1)**2/12+SLN(3)*SN(3)*RS(2,1)**3/12


	AMAX=0.0
	DO 10 I=1,3
	  DO 10 J=I,3
	    IF (DABS(AB1(I,J)).GE.AMAX) AMAX=DABS(AB1(I,J))
10	CONTINUE
	DO 20 I=1,3
	  DO 20 J=I,3
	    RES=DABS(AB1(I,J)/AMAX)
	    IF (I.NE.J.AND.RES.LE.1.0E-8) AB1(I,J)=0.0
20	CONTINUE


C	LOWER TRIANGLE
	AB1(2,1)=AB1(1,2)

	AB1(3,1)=AB1(1,3)
	AB1(3,2)=AB1(2,3)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRIAB2(RS,SLN,RN,SN,AB1)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR STRAIN PARAMETERS AB1
C
C	INPUT VARIABLES
C	RS(2,3)		= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)		= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C
C	OUTPUT VARIAVBLES
C	AB1(3,3)	= STRAIN PARAMETERS
C	-------------------------------------------------------
      DIMENSION RS(2,3),SLN(3),RN(3),SN(3)
	DIMENSION AB1(3,3)

C	-------------
C	SOLVE FOR AB1
C	-------------
      AB1(1,1) = SLN(1)*RN(1)*RS(1,1)/2+SLN(1)*RN(1)*RS(1,2)/2+SLN(2)*RN
     #(2)*RS(1,2)/2+SLN(2)*RN(2)*RS(1,3)/2+SLN(3)*RN(3)*RS(1,3)/2+SLN(3)
     #*RN(3)*RS(1,1)/2

      AB1(1,2) = 0.D0

      AB1(1,3) = 0.D0

      AB1(2,1) = 0.D0

      AB1(2,2) = 0.D0

      AB1(2,3) = 0.D0

      AB1(3,1) = 0.D0

      AB1(3,2) = 0.D0

      AB1(3,3) = 0.D0

	AMAX=0.0
	DO 315 I=1,3
	  DO 315 J=I,3
	    IF (DABS(AB1(I,J)).GE.AMAX) AMAX=DABS(AB1(I,J))
315	CONTINUE
	DO 305 I=1,3
	  DO 305 J=I,3
	    RES=DABS(AB1(I,J)/AMAX)
	    IF (I.NE.J.AND.RES.LE.1.0E-8) AB1(I,J)=0.0
305	CONTINUE

	AB1(2,1)=AB1(1,2)
	AB1(3,1)=AB1(1,3)
	AB1(3,2)=AB1(2,3)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TRIAB7(RS,SLN,RN,SN,AB1)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR STRAIN PARAMETERS AB1
C
C	INPUT VARIABLES
C	RS(2,3)		= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)		= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C
C	OUTPUT VARIAVBLES
C	AB1(3,3)	= STRAIN PARAMETERS
C	-------------------------------------------------------
      DIMENSION RS(2,3),SLN(3),RN(3),SN(3)
	DIMENSION AB1(3,3)

C	----------------------
C	SOLVE FOR ELEMENT AREA
C	----------------------
	AREA = 0.5*( (RS(1,1)*RS(2,2)+RS(1,2)*RS(2,3)+RS(1,3)*RS(2,1))-
     #             (RS(2,1)*RS(1,2)+RS(2,2)*RS(1,3)+RS(2,3)*RS(1,1))  )

	FACT = AREA/12.0

C	-------------
C	SOLVE FOR AB1
C	-------------
      AB1(1,1) = 2.0*FACT
      AB1(1,2) = 1.0*FACT
      AB1(1,3) = 1.0*FACT

      AB1(2,1) = 1.0*FACT
      AB1(2,2) = 2.0*FACT
      AB1(2,3) = 1.0*FACT

      AB1(3,1) = 1.0*FACT
      AB1(3,2) = 1.0*FACT
      AB1(3,3) = 2.0*FACT

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TABCB(AB1,CB,ACB)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-----------------------------------------------
C	PURPOSE:	TO SOLVE FOR INVERSE(Ab)*(Cb)=ACB
C
C	INPUT VARIABLES
C	AB1(3,3)	    	= STRAIN PARAMETERS
C	CB(9,18)			= BENDING STRAIN PARAMETERS
C
C	LOCAL VARIABLES
C	AI1(3,3)			= INVERSE OF AB1
C
C	OUTPUT VARIAVBLES
C	ACB(9,18)			= INVERSE(Ab)*(CB)=ACB
C	-----------------------------------------------
      DIMENSION AB1(3,3),CB(9,18)
	DIMENSION ABI(3,3)
	DIMENSION ACB(9,18)

C	---------------------------------
C	DETERMINE THE INVERSE OF Ab - AI1
C	---------------------------------
	CALL  INVMAT(AB1,ABI,3)

C	----------------------
C	INVERSE(Ab)*(CB) - ACB
C	----------------------
      ACB(1,1) = 0.E0
      ACB(1,2) = 0.E0
      ACB(1,3) = ABI(1,1)*CB(1,3)+ABI(1,2)*CB(2,3)+ABI(1,3)*CB(3,3)
      ACB(1,4) = ABI(1,1)*CB(1,4)+ABI(1,2)*CB(2,4)+ABI(1,3)*CB(3,4)
      ACB(1,5) = ABI(1,1)*CB(1,5)+ABI(1,2)*CB(2,5)+ABI(1,3)*CB(3,5)
      ACB(1,6) = 0.E0
      ACB(1,7) = 0.E0
      ACB(1,8) = 0.E0
      ACB(1,9) = ABI(1,1)*CB(1,9)+ABI(1,2)*CB(2,9)+ABI(1,3)*CB(3,9)
      ACB(1,10) = ABI(1,1)*CB(1,10)+ABI(1,2)*CB(2,10)+ABI(1,3)*CB(3,10)
      ACB(1,11) = ABI(1,1)*CB(1,11)+ABI(1,2)*CB(2,11)+ABI(1,3)*CB(3,11)
      ACB(1,12) = 0.E0
      ACB(1,13) = 0.E0
      ACB(1,14) = 0.E0
      ACB(1,15) = ABI(1,1)*CB(1,15)+ABI(1,2)*CB(2,15)+ABI(1,3)*CB(3,15)
      ACB(1,16) = ABI(1,1)*CB(1,16)+ABI(1,2)*CB(2,16)+ABI(1,3)*CB(3,16)
      ACB(1,17) = ABI(1,1)*CB(1,17)+ABI(1,2)*CB(2,17)+ABI(1,3)*CB(3,17)
      ACB(1,18) = 0.E0
      ACB(2,1) = 0.E0
      ACB(2,2) = 0.E0
      ACB(2,3) = ABI(1,2)*CB(1,3)+ABI(2,2)*CB(2,3)+ABI(2,3)*CB(3,3)
      ACB(2,4) = ABI(1,2)*CB(1,4)+ABI(2,2)*CB(2,4)+ABI(2,3)*CB(3,4)
      ACB(2,5) = ABI(1,2)*CB(1,5)+ABI(2,2)*CB(2,5)+ABI(2,3)*CB(3,5)
      ACB(2,6) = 0.E0
      ACB(2,7) = 0.E0
      ACB(2,8) = 0.E0
      ACB(2,9) = ABI(1,2)*CB(1,9)+ABI(2,2)*CB(2,9)+ABI(2,3)*CB(3,9)
      ACB(2,10) = ABI(1,2)*CB(1,10)+ABI(2,2)*CB(2,10)+ABI(2,3)*CB(3,10)
      ACB(2,11) = ABI(1,2)*CB(1,11)+ABI(2,2)*CB(2,11)+ABI(2,3)*CB(3,11)
      ACB(2,12) = 0.E0
      ACB(2,13) = 0.E0
      ACB(2,14) = 0.E0
      ACB(2,15) = ABI(1,2)*CB(1,15)+ABI(2,2)*CB(2,15)+ABI(2,3)*CB(3,15)
      ACB(2,16) = ABI(1,2)*CB(1,16)+ABI(2,2)*CB(2,16)+ABI(2,3)*CB(3,16)
      ACB(2,17) = ABI(1,2)*CB(1,17)+ABI(2,2)*CB(2,17)+ABI(2,3)*CB(3,17)
      ACB(2,18) = 0.E0
      ACB(3,1) = 0.E0
      ACB(3,2) = 0.E0
      ACB(3,3) = ABI(1,3)*CB(1,3)+ABI(2,3)*CB(2,3)+ABI(3,3)*CB(3,3)
      ACB(3,4) = ABI(1,3)*CB(1,4)+ABI(2,3)*CB(2,4)+ABI(3,3)*CB(3,4)
      ACB(3,5) = ABI(1,3)*CB(1,5)+ABI(2,3)*CB(2,5)+ABI(3,3)*CB(3,5)
      ACB(3,6) = 0.E0
      ACB(3,7) = 0.E0
      ACB(3,8) = 0.E0
      ACB(3,9) = ABI(1,3)*CB(1,9)+ABI(2,3)*CB(2,9)+ABI(3,3)*CB(3,9)
      ACB(3,10) = ABI(1,3)*CB(1,10)+ABI(2,3)*CB(2,10)+ABI(3,3)*CB(3,10)
      ACB(3,11) = ABI(1,3)*CB(1,11)+ABI(2,3)*CB(2,11)+ABI(3,3)*CB(3,11)
      ACB(3,12) = 0.E0
      ACB(3,13) = 0.E0
      ACB(3,14) = 0.E0
      ACB(3,15) = ABI(1,3)*CB(1,15)+ABI(2,3)*CB(2,15)+ABI(3,3)*CB(3,15)
      ACB(3,16) = ABI(1,3)*CB(1,16)+ABI(2,3)*CB(2,16)+ABI(3,3)*CB(3,16)
      ACB(3,17) = ABI(1,3)*CB(1,17)+ABI(2,3)*CB(2,17)+ABI(3,3)*CB(3,17)
      ACB(3,18) = 0.E0
      ACB(4,1) = 0.E0
      ACB(4,2) = 0.E0
      ACB(4,3) = ABI(1,1)*CB(4,3)+ABI(1,2)*CB(5,3)+ABI(1,3)*CB(6,3)
      ACB(4,4) = ABI(1,1)*CB(4,4)+ABI(1,2)*CB(5,4)+ABI(1,3)*CB(6,4)
      ACB(4,5) = ABI(1,1)*CB(4,5)+ABI(1,2)*CB(5,5)+ABI(1,3)*CB(6,5)
      ACB(4,6) = 0.E0
      ACB(4,7) = 0.E0
      ACB(4,8) = 0.E0
      ACB(4,9) = ABI(1,1)*CB(4,9)+ABI(1,2)*CB(5,9)+ABI(1,3)*CB(6,9)
      ACB(4,10) = ABI(1,1)*CB(4,10)+ABI(1,2)*CB(5,10)+ABI(1,3)*CB(6,10)
      ACB(4,11) = ABI(1,1)*CB(4,11)+ABI(1,2)*CB(5,11)+ABI(1,3)*CB(6,11)
      ACB(4,12) = 0.E0
      ACB(4,13) = 0.E0
      ACB(4,14) = 0.E0
      ACB(4,15) = ABI(1,1)*CB(4,15)+ABI(1,2)*CB(5,15)+ABI(1,3)*CB(6,15)
      ACB(4,16) = ABI(1,1)*CB(4,16)+ABI(1,2)*CB(5,16)+ABI(1,3)*CB(6,16)
      ACB(4,17) = ABI(1,1)*CB(4,17)+ABI(1,2)*CB(5,17)+ABI(1,3)*CB(6,17)
      ACB(4,18) = 0.E0
      ACB(5,1) = 0.E0
      ACB(5,2) = 0.E0
      ACB(5,3) = ABI(1,2)*CB(4,3)+ABI(2,2)*CB(5,3)+ABI(2,3)*CB(6,3)
      ACB(5,4) = ABI(1,2)*CB(4,4)+ABI(2,2)*CB(5,4)+ABI(2,3)*CB(6,4)
      ACB(5,5) = ABI(1,2)*CB(4,5)+ABI(2,2)*CB(5,5)+ABI(2,3)*CB(6,5)
      ACB(5,6) = 0.E0
      ACB(5,7) = 0.E0
      ACB(5,8) = 0.E0
      ACB(5,9) = ABI(1,2)*CB(4,9)+ABI(2,2)*CB(5,9)+ABI(2,3)*CB(6,9)
      ACB(5,10) = ABI(1,2)*CB(4,10)+ABI(2,2)*CB(5,10)+ABI(2,3)*CB(6,10)
      ACB(5,11) = ABI(1,2)*CB(4,11)+ABI(2,2)*CB(5,11)+ABI(2,3)*CB(6,11)
      ACB(5,12) = 0.E0
      ACB(5,13) = 0.E0
      ACB(5,14) = 0.E0
      ACB(5,15) = ABI(1,2)*CB(4,15)+ABI(2,2)*CB(5,15)+ABI(2,3)*CB(6,15)
      ACB(5,16) = ABI(1,2)*CB(4,16)+ABI(2,2)*CB(5,16)+ABI(2,3)*CB(6,16)
      ACB(5,17) = ABI(1,2)*CB(4,17)+ABI(2,2)*CB(5,17)+ABI(2,3)*CB(6,17)
      ACB(5,18) = 0.E0
      ACB(6,1) = 0.E0
      ACB(6,2) = 0.E0
      ACB(6,3) = ABI(1,3)*CB(4,3)+ABI(2,3)*CB(5,3)+ABI(3,3)*CB(6,3)
      ACB(6,4) = ABI(1,3)*CB(4,4)+ABI(2,3)*CB(5,4)+ABI(3,3)*CB(6,4)
      ACB(6,5) = ABI(1,3)*CB(4,5)+ABI(2,3)*CB(5,5)+ABI(3,3)*CB(6,5)
      ACB(6,6) = 0.E0
      ACB(6,7) = 0.E0
      ACB(6,8) = 0.E0
      ACB(6,9) = ABI(1,3)*CB(4,9)+ABI(2,3)*CB(5,9)+ABI(3,3)*CB(6,9)
      ACB(6,10) = ABI(1,3)*CB(4,10)+ABI(2,3)*CB(5,10)+ABI(3,3)*CB(6,10)
      ACB(6,11) = ABI(1,3)*CB(4,11)+ABI(2,3)*CB(5,11)+ABI(3,3)*CB(6,11)
      ACB(6,12) = 0.E0
      ACB(6,13) = 0.E0
      ACB(6,14) = 0.E0
      ACB(6,15) = ABI(1,3)*CB(4,15)+ABI(2,3)*CB(5,15)+ABI(3,3)*CB(6,15)
      ACB(6,16) = ABI(1,3)*CB(4,16)+ABI(2,3)*CB(5,16)+ABI(3,3)*CB(6,16)
      ACB(6,17) = ABI(1,3)*CB(4,17)+ABI(2,3)*CB(5,17)+ABI(3,3)*CB(6,17)
      ACB(6,18) = 0.E0
      ACB(7,1) = 0.E0
      ACB(7,2) = 0.E0
      ACB(7,3) = ABI(1,1)*CB(7,3)+ABI(1,2)*CB(8,3)+ABI(1,3)*CB(9,3)
      ACB(7,4) = ABI(1,1)*CB(7,4)+ABI(1,2)*CB(8,4)+ABI(1,3)*CB(9,4)
      ACB(7,5) = ABI(1,1)*CB(7,5)+ABI(1,2)*CB(8,5)+ABI(1,3)*CB(9,5)
      ACB(7,6) = 0.E0
      ACB(7,7) = 0.E0
      ACB(7,8) = 0.E0
      ACB(7,9) = ABI(1,1)*CB(7,9)+ABI(1,2)*CB(8,9)+ABI(1,3)*CB(9,9)
      ACB(7,10) = ABI(1,1)*CB(7,10)+ABI(1,2)*CB(8,10)+ABI(1,3)*CB(9,10)
      ACB(7,11) = ABI(1,1)*CB(7,11)+ABI(1,2)*CB(8,11)+ABI(1,3)*CB(9,11)
      ACB(7,12) = 0.E0
      ACB(7,13) = 0.E0
      ACB(7,14) = 0.E0
      ACB(7,15) = ABI(1,1)*CB(7,15)+ABI(1,2)*CB(8,15)+ABI(1,3)*CB(9,15)
      ACB(7,16) = ABI(1,1)*CB(7,16)+ABI(1,2)*CB(8,16)+ABI(1,3)*CB(9,16)
      ACB(7,17) = ABI(1,1)*CB(7,17)+ABI(1,2)*CB(8,17)+ABI(1,3)*CB(9,17)
      ACB(7,18) = 0.E0
      ACB(8,1) = 0.E0
      ACB(8,2) = 0.E0
      ACB(8,3) = ABI(1,2)*CB(7,3)+ABI(2,2)*CB(8,3)+ABI(2,3)*CB(9,3)
      ACB(8,4) = ABI(1,2)*CB(7,4)+ABI(2,2)*CB(8,4)+ABI(2,3)*CB(9,4)
      ACB(8,5) = ABI(1,2)*CB(7,5)+ABI(2,2)*CB(8,5)+ABI(2,3)*CB(9,5)
      ACB(8,6) = 0.E0
      ACB(8,7) = 0.E0
      ACB(8,8) = 0.E0
      ACB(8,9) = ABI(1,2)*CB(7,9)+ABI(2,2)*CB(8,9)+ABI(2,3)*CB(9,9)
      ACB(8,10) = ABI(1,2)*CB(7,10)+ABI(2,2)*CB(8,10)+ABI(2,3)*CB(9,10)
      ACB(8,11) = ABI(1,2)*CB(7,11)+ABI(2,2)*CB(8,11)+ABI(2,3)*CB(9,11)
      ACB(8,12) = 0.E0
      ACB(8,13) = 0.E0
      ACB(8,14) = 0.E0
      ACB(8,15) = ABI(1,2)*CB(7,15)+ABI(2,2)*CB(8,15)+ABI(2,3)*CB(9,15)
      ACB(8,16) = ABI(1,2)*CB(7,16)+ABI(2,2)*CB(8,16)+ABI(2,3)*CB(9,16)
      ACB(8,17) = ABI(1,2)*CB(7,17)+ABI(2,2)*CB(8,17)+ABI(2,3)*CB(9,17)
      ACB(8,18) = 0.E0
      ACB(9,1) = 0.E0
      ACB(9,2) = 0.E0
      ACB(9,3) = ABI(1,3)*CB(7,3)+ABI(2,3)*CB(8,3)+ABI(3,3)*CB(9,3)
      ACB(9,4) = ABI(1,3)*CB(7,4)+ABI(2,3)*CB(8,4)+ABI(3,3)*CB(9,4)
      ACB(9,5) = ABI(1,3)*CB(7,5)+ABI(2,3)*CB(8,5)+ABI(3,3)*CB(9,5)
      ACB(9,6) = 0.E0
      ACB(9,7) = 0.E0
      ACB(9,8) = 0.E0
      ACB(9,9) = ABI(1,3)*CB(7,9)+ABI(2,3)*CB(8,9)+ABI(3,3)*CB(9,9)
      ACB(9,10) = ABI(1,3)*CB(7,10)+ABI(2,3)*CB(8,10)+ABI(3,3)*CB(9,10)
      ACB(9,11) = ABI(1,3)*CB(7,11)+ABI(2,3)*CB(8,11)+ABI(3,3)*CB(9,11)
      ACB(9,12) = 0.E0
      ACB(9,13) = 0.E0
      ACB(9,14) = 0.E0
      ACB(9,15) = ABI(1,3)*CB(7,15)+ABI(2,3)*CB(8,15)+ABI(3,3)*CB(9,15)
      ACB(9,16) = ABI(1,3)*CB(7,16)+ABI(2,3)*CB(8,16)+ABI(3,3)*CB(9,16)
      ACB(9,17) = ABI(1,3)*CB(7,17)+ABI(2,3)*CB(8,17)+ABI(3,3)*CB(9,17)
      ACB(9,18) = 0.E0

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TABCB2(AB1,CB,ACB)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-----------------------------------------------
C	PURPOSE:	TO SOLVE FOR INVERSE(Ab)*(Cb)=ACB
C
C	INPUT VARIABLES
C	AB1(3,3)	    	= STRAIN PARAMETERS
C	CB(9,18)			= BENDING STRAIN PARAMETERS
C
C	LOCAL VARIABLES
C	AI1(3,3)			= INVERSE OF AB1
C
C	OUTPUT VARIAVBLES
C	ACB(9,18)			= INVERSE(Ab)*(CB)=ACB
C	-----------------------------------------------
      DIMENSION AB1(3,3),CB(9,18)
	DIMENSION ABI(3,3)
	DIMENSION ACB(9,18)

C	---------------------------------
C	DETERMINE THE INVERSE OF Ab - AI1
C	---------------------------------
c	CALL  INVMAT(AB1,ABI,3)
	abi = 0.d0
	abi(1,1) = 1.0/ab1(1,1)
C	----------------------
C	INVERSE(Ab)*(CB) - ACB
C	----------------------
      ACB(1,1) = 0.E0
      ACB(1,2) = 0.E0
      ACB(1,3) = ABI(1,1)*CB(1,3)+ABI(1,2)*CB(2,3)+ABI(1,3)*CB(3,3)
      ACB(1,4) = ABI(1,1)*CB(1,4)+ABI(1,2)*CB(2,4)+ABI(1,3)*CB(3,4)
      ACB(1,5) = ABI(1,1)*CB(1,5)+ABI(1,2)*CB(2,5)+ABI(1,3)*CB(3,5)
      ACB(1,6) = 0.E0
      ACB(1,7) = 0.E0
      ACB(1,8) = 0.E0
      ACB(1,9) = ABI(1,1)*CB(1,9)+ABI(1,2)*CB(2,9)+ABI(1,3)*CB(3,9)
      ACB(1,10) = ABI(1,1)*CB(1,10)+ABI(1,2)*CB(2,10)+ABI(1,3)*CB(3,10)
      ACB(1,11) = ABI(1,1)*CB(1,11)+ABI(1,2)*CB(2,11)+ABI(1,3)*CB(3,11)
      ACB(1,12) = 0.E0
      ACB(1,13) = 0.E0
      ACB(1,14) = 0.E0
      ACB(1,15) = ABI(1,1)*CB(1,15)+ABI(1,2)*CB(2,15)+ABI(1,3)*CB(3,15)
      ACB(1,16) = ABI(1,1)*CB(1,16)+ABI(1,2)*CB(2,16)+ABI(1,3)*CB(3,16)
      ACB(1,17) = ABI(1,1)*CB(1,17)+ABI(1,2)*CB(2,17)+ABI(1,3)*CB(3,17)
      ACB(1,18) = 0.E0
      ACB(2,1) = 0.E0
      ACB(2,2) = 0.E0
      ACB(2,3) = ABI(1,2)*CB(1,3)+ABI(2,2)*CB(2,3)+ABI(2,3)*CB(3,3)
      ACB(2,4) = ABI(1,2)*CB(1,4)+ABI(2,2)*CB(2,4)+ABI(2,3)*CB(3,4)
      ACB(2,5) = ABI(1,2)*CB(1,5)+ABI(2,2)*CB(2,5)+ABI(2,3)*CB(3,5)
      ACB(2,6) = 0.E0
      ACB(2,7) = 0.E0
      ACB(2,8) = 0.E0
      ACB(2,9) = ABI(1,2)*CB(1,9)+ABI(2,2)*CB(2,9)+ABI(2,3)*CB(3,9)
      ACB(2,10) = ABI(1,2)*CB(1,10)+ABI(2,2)*CB(2,10)+ABI(2,3)*CB(3,10)
      ACB(2,11) = ABI(1,2)*CB(1,11)+ABI(2,2)*CB(2,11)+ABI(2,3)*CB(3,11)
      ACB(2,12) = 0.E0
      ACB(2,13) = 0.E0
      ACB(2,14) = 0.E0
      ACB(2,15) = ABI(1,2)*CB(1,15)+ABI(2,2)*CB(2,15)+ABI(2,3)*CB(3,15)
      ACB(2,16) = ABI(1,2)*CB(1,16)+ABI(2,2)*CB(2,16)+ABI(2,3)*CB(3,16)
      ACB(2,17) = ABI(1,2)*CB(1,17)+ABI(2,2)*CB(2,17)+ABI(2,3)*CB(3,17)
      ACB(2,18) = 0.E0
      ACB(3,1) = 0.E0
      ACB(3,2) = 0.E0
      ACB(3,3) = ABI(1,3)*CB(1,3)+ABI(2,3)*CB(2,3)+ABI(3,3)*CB(3,3)
      ACB(3,4) = ABI(1,3)*CB(1,4)+ABI(2,3)*CB(2,4)+ABI(3,3)*CB(3,4)
      ACB(3,5) = ABI(1,3)*CB(1,5)+ABI(2,3)*CB(2,5)+ABI(3,3)*CB(3,5)
      ACB(3,6) = 0.E0
      ACB(3,7) = 0.E0
      ACB(3,8) = 0.E0
      ACB(3,9) = ABI(1,3)*CB(1,9)+ABI(2,3)*CB(2,9)+ABI(3,3)*CB(3,9)
      ACB(3,10) = ABI(1,3)*CB(1,10)+ABI(2,3)*CB(2,10)+ABI(3,3)*CB(3,10)
      ACB(3,11) = ABI(1,3)*CB(1,11)+ABI(2,3)*CB(2,11)+ABI(3,3)*CB(3,11)
      ACB(3,12) = 0.E0
      ACB(3,13) = 0.E0
      ACB(3,14) = 0.E0
      ACB(3,15) = ABI(1,3)*CB(1,15)+ABI(2,3)*CB(2,15)+ABI(3,3)*CB(3,15)
      ACB(3,16) = ABI(1,3)*CB(1,16)+ABI(2,3)*CB(2,16)+ABI(3,3)*CB(3,16)
      ACB(3,17) = ABI(1,3)*CB(1,17)+ABI(2,3)*CB(2,17)+ABI(3,3)*CB(3,17)
      ACB(3,18) = 0.E0
      ACB(4,1) = 0.E0
      ACB(4,2) = 0.E0
      ACB(4,3) = ABI(1,1)*CB(4,3)+ABI(1,2)*CB(5,3)+ABI(1,3)*CB(6,3)
      ACB(4,4) = ABI(1,1)*CB(4,4)+ABI(1,2)*CB(5,4)+ABI(1,3)*CB(6,4)
      ACB(4,5) = ABI(1,1)*CB(4,5)+ABI(1,2)*CB(5,5)+ABI(1,3)*CB(6,5)
      ACB(4,6) = 0.E0
      ACB(4,7) = 0.E0
      ACB(4,8) = 0.E0
      ACB(4,9) = ABI(1,1)*CB(4,9)+ABI(1,2)*CB(5,9)+ABI(1,3)*CB(6,9)
      ACB(4,10) = ABI(1,1)*CB(4,10)+ABI(1,2)*CB(5,10)+ABI(1,3)*CB(6,10)
      ACB(4,11) = ABI(1,1)*CB(4,11)+ABI(1,2)*CB(5,11)+ABI(1,3)*CB(6,11)
      ACB(4,12) = 0.E0
      ACB(4,13) = 0.E0
      ACB(4,14) = 0.E0
      ACB(4,15) = ABI(1,1)*CB(4,15)+ABI(1,2)*CB(5,15)+ABI(1,3)*CB(6,15)
      ACB(4,16) = ABI(1,1)*CB(4,16)+ABI(1,2)*CB(5,16)+ABI(1,3)*CB(6,16)
      ACB(4,17) = ABI(1,1)*CB(4,17)+ABI(1,2)*CB(5,17)+ABI(1,3)*CB(6,17)
      ACB(4,18) = 0.E0
      ACB(5,1) = 0.E0
      ACB(5,2) = 0.E0
      ACB(5,3) = ABI(1,2)*CB(4,3)+ABI(2,2)*CB(5,3)+ABI(2,3)*CB(6,3)
      ACB(5,4) = ABI(1,2)*CB(4,4)+ABI(2,2)*CB(5,4)+ABI(2,3)*CB(6,4)
      ACB(5,5) = ABI(1,2)*CB(4,5)+ABI(2,2)*CB(5,5)+ABI(2,3)*CB(6,5)
      ACB(5,6) = 0.E0
      ACB(5,7) = 0.E0
      ACB(5,8) = 0.E0
      ACB(5,9) = ABI(1,2)*CB(4,9)+ABI(2,2)*CB(5,9)+ABI(2,3)*CB(6,9)
      ACB(5,10) = ABI(1,2)*CB(4,10)+ABI(2,2)*CB(5,10)+ABI(2,3)*CB(6,10)
      ACB(5,11) = ABI(1,2)*CB(4,11)+ABI(2,2)*CB(5,11)+ABI(2,3)*CB(6,11)
      ACB(5,12) = 0.E0
      ACB(5,13) = 0.E0
      ACB(5,14) = 0.E0
      ACB(5,15) = ABI(1,2)*CB(4,15)+ABI(2,2)*CB(5,15)+ABI(2,3)*CB(6,15)
      ACB(5,16) = ABI(1,2)*CB(4,16)+ABI(2,2)*CB(5,16)+ABI(2,3)*CB(6,16)
      ACB(5,17) = ABI(1,2)*CB(4,17)+ABI(2,2)*CB(5,17)+ABI(2,3)*CB(6,17)
      ACB(5,18) = 0.E0
      ACB(6,1) = 0.E0
      ACB(6,2) = 0.E0
      ACB(6,3) = ABI(1,3)*CB(4,3)+ABI(2,3)*CB(5,3)+ABI(3,3)*CB(6,3)
      ACB(6,4) = ABI(1,3)*CB(4,4)+ABI(2,3)*CB(5,4)+ABI(3,3)*CB(6,4)
      ACB(6,5) = ABI(1,3)*CB(4,5)+ABI(2,3)*CB(5,5)+ABI(3,3)*CB(6,5)
      ACB(6,6) = 0.E0
      ACB(6,7) = 0.E0
      ACB(6,8) = 0.E0
      ACB(6,9) = ABI(1,3)*CB(4,9)+ABI(2,3)*CB(5,9)+ABI(3,3)*CB(6,9)
      ACB(6,10) = ABI(1,3)*CB(4,10)+ABI(2,3)*CB(5,10)+ABI(3,3)*CB(6,10)
      ACB(6,11) = ABI(1,3)*CB(4,11)+ABI(2,3)*CB(5,11)+ABI(3,3)*CB(6,11)
      ACB(6,12) = 0.E0
      ACB(6,13) = 0.E0
      ACB(6,14) = 0.E0
      ACB(6,15) = ABI(1,3)*CB(4,15)+ABI(2,3)*CB(5,15)+ABI(3,3)*CB(6,15)
      ACB(6,16) = ABI(1,3)*CB(4,16)+ABI(2,3)*CB(5,16)+ABI(3,3)*CB(6,16)
      ACB(6,17) = ABI(1,3)*CB(4,17)+ABI(2,3)*CB(5,17)+ABI(3,3)*CB(6,17)
      ACB(6,18) = 0.E0
      ACB(7,1) = 0.E0
      ACB(7,2) = 0.E0
      ACB(7,3) = ABI(1,1)*CB(7,3)+ABI(1,2)*CB(8,3)+ABI(1,3)*CB(9,3)
      ACB(7,4) = ABI(1,1)*CB(7,4)+ABI(1,2)*CB(8,4)+ABI(1,3)*CB(9,4)
      ACB(7,5) = ABI(1,1)*CB(7,5)+ABI(1,2)*CB(8,5)+ABI(1,3)*CB(9,5)
      ACB(7,6) = 0.E0
      ACB(7,7) = 0.E0
      ACB(7,8) = 0.E0
      ACB(7,9) = ABI(1,1)*CB(7,9)+ABI(1,2)*CB(8,9)+ABI(1,3)*CB(9,9)
      ACB(7,10) = ABI(1,1)*CB(7,10)+ABI(1,2)*CB(8,10)+ABI(1,3)*CB(9,10)
      ACB(7,11) = ABI(1,1)*CB(7,11)+ABI(1,2)*CB(8,11)+ABI(1,3)*CB(9,11)
      ACB(7,12) = 0.E0
      ACB(7,13) = 0.E0
      ACB(7,14) = 0.E0
      ACB(7,15) = ABI(1,1)*CB(7,15)+ABI(1,2)*CB(8,15)+ABI(1,3)*CB(9,15)
      ACB(7,16) = ABI(1,1)*CB(7,16)+ABI(1,2)*CB(8,16)+ABI(1,3)*CB(9,16)
      ACB(7,17) = ABI(1,1)*CB(7,17)+ABI(1,2)*CB(8,17)+ABI(1,3)*CB(9,17)
      ACB(7,18) = 0.E0
      ACB(8,1) = 0.E0
      ACB(8,2) = 0.E0
      ACB(8,3) = ABI(1,2)*CB(7,3)+ABI(2,2)*CB(8,3)+ABI(2,3)*CB(9,3)
      ACB(8,4) = ABI(1,2)*CB(7,4)+ABI(2,2)*CB(8,4)+ABI(2,3)*CB(9,4)
      ACB(8,5) = ABI(1,2)*CB(7,5)+ABI(2,2)*CB(8,5)+ABI(2,3)*CB(9,5)
      ACB(8,6) = 0.E0
      ACB(8,7) = 0.E0
      ACB(8,8) = 0.E0
      ACB(8,9) = ABI(1,2)*CB(7,9)+ABI(2,2)*CB(8,9)+ABI(2,3)*CB(9,9)
      ACB(8,10) = ABI(1,2)*CB(7,10)+ABI(2,2)*CB(8,10)+ABI(2,3)*CB(9,10)
      ACB(8,11) = ABI(1,2)*CB(7,11)+ABI(2,2)*CB(8,11)+ABI(2,3)*CB(9,11)
      ACB(8,12) = 0.E0
      ACB(8,13) = 0.E0
      ACB(8,14) = 0.E0
      ACB(8,15) = ABI(1,2)*CB(7,15)+ABI(2,2)*CB(8,15)+ABI(2,3)*CB(9,15)
      ACB(8,16) = ABI(1,2)*CB(7,16)+ABI(2,2)*CB(8,16)+ABI(2,3)*CB(9,16)
      ACB(8,17) = ABI(1,2)*CB(7,17)+ABI(2,2)*CB(8,17)+ABI(2,3)*CB(9,17)
      ACB(8,18) = 0.E0
      ACB(9,1) = 0.E0
      ACB(9,2) = 0.E0
      ACB(9,3) = ABI(1,3)*CB(7,3)+ABI(2,3)*CB(8,3)+ABI(3,3)*CB(9,3)
      ACB(9,4) = ABI(1,3)*CB(7,4)+ABI(2,3)*CB(8,4)+ABI(3,3)*CB(9,4)
      ACB(9,5) = ABI(1,3)*CB(7,5)+ABI(2,3)*CB(8,5)+ABI(3,3)*CB(9,5)
      ACB(9,6) = 0.E0
      ACB(9,7) = 0.E0
      ACB(9,8) = 0.E0
      ACB(9,9) = ABI(1,3)*CB(7,9)+ABI(2,3)*CB(8,9)+ABI(3,3)*CB(9,9)
      ACB(9,10) = ABI(1,3)*CB(7,10)+ABI(2,3)*CB(8,10)+ABI(3,3)*CB(9,10)
      ACB(9,11) = ABI(1,3)*CB(7,11)+ABI(2,3)*CB(8,11)+ABI(3,3)*CB(9,11)
      ACB(9,12) = 0.E0
      ACB(9,13) = 0.E0
      ACB(9,14) = 0.E0
      ACB(9,15) = ABI(1,3)*CB(7,15)+ABI(2,3)*CB(8,15)+ABI(3,3)*CB(9,15)
      ACB(9,16) = ABI(1,3)*CB(7,16)+ABI(2,3)*CB(8,16)+ABI(3,3)*CB(9,16)
      ACB(9,17) = ABI(1,3)*CB(7,17)+ABI(2,3)*CB(8,17)+ABI(3,3)*CB(9,17)
      ACB(9,18) = 0.E0

	RETURN
	END
C
C=====================================================================
      SUBROUTINE NBDPB18(AB1,DR,PDPB)

	IMPLICIT REAL*8 (A-H,O-Z)
C	---------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Nb*D*Pb - PDPB
C
C	INPUT VARIABLES
C	AB1(3,3)	= STRAIN PARAMETERS
C	DR(64)		= MEMBERANE,BENDING AND SHEAR STIFFNESS
C
C	OUTPUT VARIAVBLES
C	PDPB(9,9)	= Nb*D*Pb
C	---------------------------------------------------
      DIMENSION AB1(3,3),DR(64)
	DIMENSION PDPB(9,9)

	PDPB = 0.0D0

      PDPB(1,1) = DR(28)*AB1(1,1)
      PDPB(1,2) = DR(28)*AB1(1,2)
      PDPB(1,3) = DR(28)*AB1(1,3)
      PDPB(1,4) = DR(29)*AB1(1,1)
      PDPB(1,5) = DR(29)*AB1(1,2)
      PDPB(1,6) = DR(29)*AB1(1,3)

      PDPB(2,1) = PDPB(1,2)
      PDPB(2,2) = DR(28)*AB1(2,2)
      PDPB(2,3) = DR(28)*AB1(2,3)
      PDPB(2,4) = DR(29)*AB1(1,2)
      PDPB(2,5) = DR(29)*AB1(2,2)
      PDPB(2,6) = DR(29)*AB1(2,3)

      PDPB(3,1) = PDPB(1,3)
      PDPB(3,2) = PDPB(2,3)
      PDPB(3,3) = DR(28)*AB1(3,3)
      PDPB(3,4) = DR(29)*AB1(1,3)
      PDPB(3,5) = DR(29)*AB1(2,3)
      PDPB(3,6) = DR(29)*AB1(3,3)

      PDPB(4,1) = PDPB(1,4)
      PDPB(4,2) = PDPB(2,4)
      PDPB(4,3) = PDPB(3,4)
      PDPB(4,4) = DR(28)*AB1(1,1)
      PDPB(4,5) = DR(28)*AB1(1,2)
      PDPB(4,6) = DR(28)*AB1(1,3)

      PDPB(5,1) = PDPB(1,5)
      PDPB(5,2) = PDPB(2,5)
      PDPB(5,3) = PDPB(3,5)
      PDPB(5,4) = PDPB(4,5)
      PDPB(5,5) = DR(28)*AB1(2,2)
      PDPB(5,6) = DR(28)*AB1(2,3)

      PDPB(6,1) = PDPB(1,6)
      PDPB(6,2) = PDPB(2,6)
      PDPB(6,3) = PDPB(3,6)
      PDPB(6,4) = PDPB(4,6)
      PDPB(6,5) = PDPB(5,6)
      PDPB(6,6) = DR(28)*AB1(3,3)

      PDPB(7,7) = DR(46)*AB1(1,1)
      PDPB(7,8) = DR(46)*AB1(1,2)
      PDPB(7,9) = DR(46)*AB1(1,3)

      PDPB(8,7) = PDPB(7,8)
      PDPB(8,8) = DR(46)*AB1(2,2)
      PDPB(8,9) = DR(46)*AB1(2,3)

      PDPB(9,7) = PDPB(7,9)
      PDPB(9,8) = PDPB(8,9)
      PDPB(9,9) = DR(46)*AB1(3,3)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE APHAB18(ACB,TEDIS,APB)

	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR PARAMETER ALPHA
C
C	INPUT VARIABLES
C	ACB(9,14)	= INVERSE(Ab)*(Cb)
C	TEDIS(18)	= TRANSPOSE(T)*REDIS
C
C	OUTPUT VARIAVBLES
C	APB(9)		= MEMBRANE STRAIN PARAMETER
C	--------------------------------------------------------
      DIMENSION ACB(9,18),TEDIS(18)
	DIMENSION APB(9)

      APB(1) = ACB(1,3)*TEDIS(3)+ACB(1,4)*TEDIS(4)+ACB(1,5)*TEDIS(5)+ACB
     #(1,9)*TEDIS(9)+ACB(1,10)*TEDIS(10)+ACB(1,11)*TEDIS(11)+ACB(1,15)*T
     #EDIS(15)+ACB(1,16)*TEDIS(16)+ACB(1,17)*TEDIS(17)

      APB(2) = ACB(2,3)*TEDIS(3)+ACB(2,4)*TEDIS(4)+ACB(2,5)*TEDIS(5)+ACB
     #(2,9)*TEDIS(9)+ACB(2,10)*TEDIS(10)+ACB(2,11)*TEDIS(11)+ACB(2,15)*T
     #EDIS(15)+ACB(2,16)*TEDIS(16)+ACB(2,17)*TEDIS(17)

      APB(3) = ACB(3,3)*TEDIS(3)+ACB(3,4)*TEDIS(4)+ACB(3,5)*TEDIS(5)+ACB
     #(3,9)*TEDIS(9)+ACB(3,10)*TEDIS(10)+ACB(3,11)*TEDIS(11)+ACB(3,15)*T
     #EDIS(15)+ACB(3,16)*TEDIS(16)+ACB(3,17)*TEDIS(17)

      APB(4) = ACB(4,3)*TEDIS(3)+ACB(4,4)*TEDIS(4)+ACB(4,5)*TEDIS(5)+ACB
     #(4,9)*TEDIS(9)+ACB(4,10)*TEDIS(10)+ACB(4,11)*TEDIS(11)+ACB(4,15)*T
     #EDIS(15)+ACB(4,16)*TEDIS(16)+ACB(4,17)*TEDIS(17)
 
      APB(5) = ACB(5,3)*TEDIS(3)+ACB(5,4)*TEDIS(4)+ACB(5,5)*TEDIS(5)+ACB
     #(5,9)*TEDIS(9)+ACB(5,10)*TEDIS(10)+ACB(5,11)*TEDIS(11)+ACB(5,15)*T
     #EDIS(15)+ACB(5,16)*TEDIS(16)+ACB(5,17)*TEDIS(17)

      APB(6) = ACB(6,3)*TEDIS(3)+ACB(6,4)*TEDIS(4)+ACB(6,5)*TEDIS(5)+ACB
     #(6,9)*TEDIS(9)+ACB(6,10)*TEDIS(10)+ACB(6,11)*TEDIS(11)+ACB(6,15)*T
     #EDIS(15)+ACB(6,16)*TEDIS(16)+ACB(6,17)*TEDIS(17)

      APB(7) = ACB(7,3)*TEDIS(3)+ACB(7,4)*TEDIS(4)+ACB(7,5)*TEDIS(5)+ACB
     #(7,9)*TEDIS(9)+ACB(7,10)*TEDIS(10)+ACB(7,11)*TEDIS(11)+ACB(7,15)*T
     #EDIS(15)+ACB(7,16)*TEDIS(16)+ACB(7,17)*TEDIS(17)

      APB(8) = ACB(8,3)*TEDIS(3)+ACB(8,4)*TEDIS(4)+ACB(8,5)*TEDIS(5)+ACB
     #(8,9)*TEDIS(9)+ACB(8,10)*TEDIS(10)+ACB(8,11)*TEDIS(11)+ACB(8,15)*T
     #EDIS(15)+ACB(8,16)*TEDIS(16)+ACB(8,17)*TEDIS(17)

      APB(9) = ACB(9,3)*TEDIS(3)+ACB(9,4)*TEDIS(4)+ACB(9,5)*TEDIS(5)+ACB
     #(9,9)*TEDIS(9)+ACB(9,10)*TEDIS(10)+ACB(9,11)*TEDIS(11)+ACB(9,15)*T
     #EDIS(15)+ACB(9,16)*TEDIS(16)+ACB(9,17)*TEDIS(17)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE KBND18(ACB,PDPB,ST)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-----------------------------------------------------------
C	PURPOSE:	TO ADD THE FLEXURAL CONTRIBUTION TO THE
C					ELEMENT STIFFNESS MATRIX
C
C	INPUT VARIABLES
C	ACB(9,18)		= INVERSE(Ab)*(Cb)
C	PDPB(9,9)		= Nb*D*Pb
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIAVBLES
C	ST(171)			= STIFFNESS MATRIX IN ROW FORMAT
C	-----------------------------------------------------------
      DIMENSION ACB(9,18),PDPB(9,9)
	DIMENSION ST(171)

      t7 = ACB(1,3)*PDPB(1,1)+ACB(2,3)*PDPB(1,2)+ACB(3,3)*PDPB(1,3)+ACB(
     #4,3)*PDPB(1,4)+ACB(5,3)*PDPB(1,5)+ACB(6,3)*PDPB(1,6)
      t15 = ACB(1,3)*PDPB(1,2)+ACB(2,3)*PDPB(2,2)+ACB(3,3)*PDPB(2,3)+ACB
     #(4,3)*PDPB(2,4)+ACB(5,3)*PDPB(2,5)+ACB(6,3)*PDPB(2,6)
      t23 = ACB(1,3)*PDPB(1,3)+ACB(2,3)*PDPB(2,3)+ACB(3,3)*PDPB(3,3)+ACB
     #(4,3)*PDPB(3,4)+ACB(5,3)*PDPB(3,5)+ACB(6,3)*PDPB(3,6)
      t31 = ACB(1,3)*PDPB(1,4)+ACB(2,3)*PDPB(2,4)+ACB(3,3)*PDPB(3,4)+ACB
     #(4,3)*PDPB(4,4)+ACB(5,3)*PDPB(4,5)+ACB(6,3)*PDPB(4,6)
      t39 = ACB(1,3)*PDPB(1,5)+ACB(2,3)*PDPB(2,5)+ACB(3,3)*PDPB(3,5)+ACB
     #(4,3)*PDPB(4,5)+ACB(5,3)*PDPB(5,5)+ACB(6,3)*PDPB(5,6)
      t47 = ACB(1,3)*PDPB(1,6)+ACB(2,3)*PDPB(2,6)+ACB(3,3)*PDPB(3,6)+ACB
     #(4,3)*PDPB(4,6)+ACB(5,3)*PDPB(5,6)+ACB(6,3)*PDPB(6,6)
      t52 = ACB(7,3)*PDPB(7,7)+ACB(8,3)*PDPB(7,8)+ACB(9,3)*PDPB(7,9)
      t57 = ACB(7,3)*PDPB(7,8)+ACB(8,3)*PDPB(8,8)+ACB(9,3)*PDPB(8,9)
      t62 = ACB(7,3)*PDPB(7,9)+ACB(8,3)*PDPB(8,9)+ACB(9,3)*PDPB(9,9)
      t151 = ACB(1,4)*PDPB(1,1)+ACB(2,4)*PDPB(1,2)+ACB(3,4)*PDPB(1,3)+AC
     #B(4,4)*PDPB(1,4)+ACB(5,4)*PDPB(1,5)+ACB(6,4)*PDPB(1,6)
      t159 = ACB(1,4)*PDPB(1,2)+ACB(2,4)*PDPB(2,2)+ACB(3,4)*PDPB(2,3)+AC
     #B(4,4)*PDPB(2,4)+ACB(5,4)*PDPB(2,5)+ACB(6,4)*PDPB(2,6)
      t167 = ACB(1,4)*PDPB(1,3)+ACB(2,4)*PDPB(2,3)+ACB(3,4)*PDPB(3,3)+AC
     #B(4,4)*PDPB(3,4)+ACB(5,4)*PDPB(3,5)+ACB(6,4)*PDPB(3,6)
      t175 = ACB(1,4)*PDPB(1,4)+ACB(2,4)*PDPB(2,4)+ACB(3,4)*PDPB(3,4)+AC
     #B(4,4)*PDPB(4,4)+ACB(5,4)*PDPB(4,5)+ACB(6,4)*PDPB(4,6)
      t183 = ACB(1,4)*PDPB(1,5)+ACB(2,4)*PDPB(2,5)+ACB(3,4)*PDPB(3,5)+AC
     #B(4,4)*PDPB(4,5)+ACB(5,4)*PDPB(5,5)+ACB(6,4)*PDPB(5,6)
      t191 = ACB(1,4)*PDPB(1,6)+ACB(2,4)*PDPB(2,6)+ACB(3,4)*PDPB(3,6)+AC
     #B(4,4)*PDPB(4,6)+ACB(5,4)*PDPB(5,6)+ACB(6,4)*PDPB(6,6)
      t196 = ACB(7,4)*PDPB(7,7)+ACB(8,4)*PDPB(7,8)+ACB(9,4)*PDPB(7,9)
      t201 = ACB(7,4)*PDPB(7,8)+ACB(8,4)*PDPB(8,8)+ACB(9,4)*PDPB(8,9)
      t206 = ACB(7,4)*PDPB(7,9)+ACB(8,4)*PDPB(8,9)+ACB(9,4)*PDPB(9,9)
      t285 = ACB(1,5)*PDPB(1,1)+ACB(2,5)*PDPB(1,2)+ACB(3,5)*PDPB(1,3)+AC
     #B(4,5)*PDPB(1,4)+ACB(5,5)*PDPB(1,5)+ACB(6,5)*PDPB(1,6)
      t293 = ACB(1,5)*PDPB(1,2)+ACB(2,5)*PDPB(2,2)+ACB(3,5)*PDPB(2,3)+AC
     #B(4,5)*PDPB(2,4)+ACB(5,5)*PDPB(2,5)+ACB(6,5)*PDPB(2,6)
      t301 = ACB(1,5)*PDPB(1,3)+ACB(2,5)*PDPB(2,3)+ACB(3,5)*PDPB(3,3)+AC
     #B(4,5)*PDPB(3,4)+ACB(5,5)*PDPB(3,5)+ACB(6,5)*PDPB(3,6)
      t309 = ACB(1,5)*PDPB(1,4)+ACB(2,5)*PDPB(2,4)+ACB(3,5)*PDPB(3,4)+AC
     #B(4,5)*PDPB(4,4)+ACB(5,5)*PDPB(4,5)+ACB(6,5)*PDPB(4,6)
      t317 = ACB(1,5)*PDPB(1,5)+ACB(2,5)*PDPB(2,5)+ACB(3,5)*PDPB(3,5)+AC
     #B(4,5)*PDPB(4,5)+ACB(5,5)*PDPB(5,5)+ACB(6,5)*PDPB(5,6)
      t325 = ACB(1,5)*PDPB(1,6)+ACB(2,5)*PDPB(2,6)+ACB(3,5)*PDPB(3,6)+AC
     #B(4,5)*PDPB(4,6)+ACB(5,5)*PDPB(5,6)+ACB(6,5)*PDPB(6,6)
      t330 = ACB(7,5)*PDPB(7,7)+ACB(8,5)*PDPB(7,8)+ACB(9,5)*PDPB(7,9)
      t335 = ACB(7,5)*PDPB(7,8)+ACB(8,5)*PDPB(8,8)+ACB(9,5)*PDPB(8,9)
      t340 = ACB(7,5)*PDPB(7,9)+ACB(8,5)*PDPB(8,9)+ACB(9,5)*PDPB(9,9)
      t409 = ACB(1,9)*PDPB(1,1)+ACB(2,9)*PDPB(1,2)+ACB(3,9)*PDPB(1,3)+AC
     #B(4,9)*PDPB(1,4)+ACB(5,9)*PDPB(1,5)+ACB(6,9)*PDPB(1,6)
      t417 = ACB(1,9)*PDPB(1,2)+ACB(2,9)*PDPB(2,2)+ACB(3,9)*PDPB(2,3)+AC
     #B(4,9)*PDPB(2,4)+ACB(5,9)*PDPB(2,5)+ACB(6,9)*PDPB(2,6)
      t425 = ACB(1,9)*PDPB(1,3)+ACB(2,9)*PDPB(2,3)+ACB(3,9)*PDPB(3,3)+AC
     #B(4,9)*PDPB(3,4)+ACB(5,9)*PDPB(3,5)+ACB(6,9)*PDPB(3,6)
      t433 = ACB(1,9)*PDPB(1,4)+ACB(2,9)*PDPB(2,4)+ACB(3,9)*PDPB(3,4)+AC
     #B(4,9)*PDPB(4,4)+ACB(5,9)*PDPB(4,5)+ACB(6,9)*PDPB(4,6)
      t441 = ACB(1,9)*PDPB(1,5)+ACB(2,9)*PDPB(2,5)+ACB(3,9)*PDPB(3,5)+AC
     #B(4,9)*PDPB(4,5)+ACB(5,9)*PDPB(5,5)+ACB(6,9)*PDPB(5,6)
      t449 = ACB(1,9)*PDPB(1,6)+ACB(2,9)*PDPB(2,6)+ACB(3,9)*PDPB(3,6)+AC
     #B(4,9)*PDPB(4,6)+ACB(5,9)*PDPB(5,6)+ACB(6,9)*PDPB(6,6)
      t454 = ACB(7,9)*PDPB(7,7)+ACB(8,9)*PDPB(7,8)+ACB(9,9)*PDPB(7,9)
      t459 = ACB(7,9)*PDPB(7,8)+ACB(8,9)*PDPB(8,8)+ACB(9,9)*PDPB(8,9)
      t464 = ACB(7,9)*PDPB(7,9)+ACB(8,9)*PDPB(8,9)+ACB(9,9)*PDPB(9,9)
      t523 = ACB(1,10)*PDPB(1,1)+ACB(2,10)*PDPB(1,2)+ACB(3,10)*PDPB(1,3)
     #+ACB(4,10)*PDPB(1,4)+ACB(5,10)*PDPB(1,5)+ACB(6,10)*PDPB(1,6)
      t531 = ACB(1,10)*PDPB(1,2)+ACB(2,10)*PDPB(2,2)+ACB(3,10)*PDPB(2,3)
     #+ACB(4,10)*PDPB(2,4)+ACB(5,10)*PDPB(2,5)+ACB(6,10)*PDPB(2,6)
      t539 = ACB(1,10)*PDPB(1,3)+ACB(2,10)*PDPB(2,3)+ACB(3,10)*PDPB(3,3)
     #+ACB(4,10)*PDPB(3,4)+ACB(5,10)*PDPB(3,5)+ACB(6,10)*PDPB(3,6)
      t547 = ACB(1,10)*PDPB(1,4)+ACB(2,10)*PDPB(2,4)+ACB(3,10)*PDPB(3,4)
     #+ACB(4,10)*PDPB(4,4)+ACB(5,10)*PDPB(4,5)+ACB(6,10)*PDPB(4,6)
      t555 = ACB(1,10)*PDPB(1,5)+ACB(2,10)*PDPB(2,5)+ACB(3,10)*PDPB(3,5)
     #+ACB(4,10)*PDPB(4,5)+ACB(5,10)*PDPB(5,5)+ACB(6,10)*PDPB(5,6)
      t563 = ACB(1,10)*PDPB(1,6)+ACB(2,10)*PDPB(2,6)+ACB(3,10)*PDPB(3,6)
     #+ACB(4,10)*PDPB(4,6)+ACB(5,10)*PDPB(5,6)+ACB(6,10)*PDPB(6,6)
      t568 = ACB(7,10)*PDPB(7,7)+ACB(8,10)*PDPB(7,8)+ACB(9,10)*PDPB(7,9)
      t573 = ACB(7,10)*PDPB(7,8)+ACB(8,10)*PDPB(8,8)+ACB(9,10)*PDPB(8,9)
      t578 = ACB(7,10)*PDPB(7,9)+ACB(8,10)*PDPB(8,9)+ACB(9,10)*PDPB(9,9)
      t627 = ACB(1,11)*PDPB(1,1)+ACB(2,11)*PDPB(1,2)+ACB(3,11)*PDPB(1,3)
     #+ACB(4,11)*PDPB(1,4)+ACB(5,11)*PDPB(1,5)+ACB(6,11)*PDPB(1,6)
      t635 = ACB(1,11)*PDPB(1,2)+ACB(2,11)*PDPB(2,2)+ACB(3,11)*PDPB(2,3)
     #+ACB(4,11)*PDPB(2,4)+ACB(5,11)*PDPB(2,5)+ACB(6,11)*PDPB(2,6)
      t643 = ACB(1,11)*PDPB(1,3)+ACB(2,11)*PDPB(2,3)+ACB(3,11)*PDPB(3,3)
     #+ACB(4,11)*PDPB(3,4)+ACB(5,11)*PDPB(3,5)+ACB(6,11)*PDPB(3,6)
      t651 = ACB(1,11)*PDPB(1,4)+ACB(2,11)*PDPB(2,4)+ACB(3,11)*PDPB(3,4)
     #+ACB(4,11)*PDPB(4,4)+ACB(5,11)*PDPB(4,5)+ACB(6,11)*PDPB(4,6)
      t659 = ACB(1,11)*PDPB(1,5)+ACB(2,11)*PDPB(2,5)+ACB(3,11)*PDPB(3,5)
     #+ACB(4,11)*PDPB(4,5)+ACB(5,11)*PDPB(5,5)+ACB(6,11)*PDPB(5,6)
      t667 = ACB(1,11)*PDPB(1,6)+ACB(2,11)*PDPB(2,6)+ACB(3,11)*PDPB(3,6)
     #+ACB(4,11)*PDPB(4,6)+ACB(5,11)*PDPB(5,6)+ACB(6,11)*PDPB(6,6)
      t672 = ACB(7,11)*PDPB(7,7)+ACB(8,11)*PDPB(7,8)+ACB(9,11)*PDPB(7,9)
      t677 = ACB(7,11)*PDPB(7,8)+ACB(8,11)*PDPB(8,8)+ACB(9,11)*PDPB(8,9)
      t682 = ACB(7,11)*PDPB(7,9)+ACB(8,11)*PDPB(8,9)+ACB(9,11)*PDPB(9,9)
      t721 = ACB(1,15)*PDPB(1,1)+ACB(2,15)*PDPB(1,2)+ACB(3,15)*PDPB(1,3)
     #+ACB(4,15)*PDPB(1,4)+ACB(5,15)*PDPB(1,5)+ACB(6,15)*PDPB(1,6)
      t729 = ACB(1,15)*PDPB(1,2)+ACB(2,15)*PDPB(2,2)+ACB(3,15)*PDPB(2,3)
     #+ACB(4,15)*PDPB(2,4)+ACB(5,15)*PDPB(2,5)+ACB(6,15)*PDPB(2,6)
      t737 = ACB(1,15)*PDPB(1,3)+ACB(2,15)*PDPB(2,3)+ACB(3,15)*PDPB(3,3)
     #+ACB(4,15)*PDPB(3,4)+ACB(5,15)*PDPB(3,5)+ACB(6,15)*PDPB(3,6)
      t745 = ACB(1,15)*PDPB(1,4)+ACB(2,15)*PDPB(2,4)+ACB(3,15)*PDPB(3,4)
     #+ACB(4,15)*PDPB(4,4)+ACB(5,15)*PDPB(4,5)+ACB(6,15)*PDPB(4,6)
      t753 = ACB(1,15)*PDPB(1,5)+ACB(2,15)*PDPB(2,5)+ACB(3,15)*PDPB(3,5)
     #+ACB(4,15)*PDPB(4,5)+ACB(5,15)*PDPB(5,5)+ACB(6,15)*PDPB(5,6)
      t761 = ACB(1,15)*PDPB(1,6)+ACB(2,15)*PDPB(2,6)+ACB(3,15)*PDPB(3,6)
     #+ACB(4,15)*PDPB(4,6)+ACB(5,15)*PDPB(5,6)+ACB(6,15)*PDPB(6,6)
      t766 = ACB(7,15)*PDPB(7,7)+ACB(8,15)*PDPB(7,8)+ACB(9,15)*PDPB(7,9)
      t771 = ACB(7,15)*PDPB(7,8)+ACB(8,15)*PDPB(8,8)+ACB(9,15)*PDPB(8,9)
      t776 = ACB(7,15)*PDPB(7,9)+ACB(8,15)*PDPB(8,9)+ACB(9,15)*PDPB(9,9)
      t805 = ACB(1,16)*PDPB(1,1)+ACB(2,16)*PDPB(1,2)+ACB(3,16)*PDPB(1,3)
     #+ACB(4,16)*PDPB(1,4)+ACB(5,16)*PDPB(1,5)+ACB(6,16)*PDPB(1,6)
      t813 = ACB(1,16)*PDPB(1,2)+ACB(2,16)*PDPB(2,2)+ACB(3,16)*PDPB(2,3)
     #+ACB(4,16)*PDPB(2,4)+ACB(5,16)*PDPB(2,5)+ACB(6,16)*PDPB(2,6)
      t821 = ACB(1,16)*PDPB(1,3)+ACB(2,16)*PDPB(2,3)+ACB(3,16)*PDPB(3,3)
     #+ACB(4,16)*PDPB(3,4)+ACB(5,16)*PDPB(3,5)+ACB(6,16)*PDPB(3,6)
      t829 = ACB(1,16)*PDPB(1,4)+ACB(2,16)*PDPB(2,4)+ACB(3,16)*PDPB(3,4)
     #+ACB(4,16)*PDPB(4,4)+ACB(5,16)*PDPB(4,5)+ACB(6,16)*PDPB(4,6)
      t837 = ACB(1,16)*PDPB(1,5)+ACB(2,16)*PDPB(2,5)+ACB(3,16)*PDPB(3,5)
     #+ACB(4,16)*PDPB(4,5)+ACB(5,16)*PDPB(5,5)+ACB(6,16)*PDPB(5,6)
      t845 = ACB(1,16)*PDPB(1,6)+ACB(2,16)*PDPB(2,6)+ACB(3,16)*PDPB(3,6)
     #+ACB(4,16)*PDPB(4,6)+ACB(5,16)*PDPB(5,6)+ACB(6,16)*PDPB(6,6)
      t850 = ACB(7,16)*PDPB(7,7)+ACB(8,16)*PDPB(7,8)+ACB(9,16)*PDPB(7,9)
      t855 = ACB(7,16)*PDPB(7,8)+ACB(8,16)*PDPB(8,8)+ACB(9,16)*PDPB(8,9)
      t860 = ACB(7,16)*PDPB(7,9)+ACB(8,16)*PDPB(8,9)+ACB(9,16)*PDPB(9,9)
      ST(1) = 0
      ST(2) = 0
      ST(3) = 0
      ST(4) = 0
      ST(5) = 0
      ST(6) = 0
      ST(7) = 0
      ST(8) = 0
      ST(9) = 0
      ST(10) = 0
      ST(11) = 0
      ST(12) = 0
      ST(13) = 0
      ST(14) = 0
      ST(15) = 0
      ST(16) = 0
      ST(17) = 0
      ST(18) = 0
      ST(19) = 0
      ST(20) = 0
      ST(21) = 0
      ST(22) = 0
      ST(23) = 0
      ST(24) = 0
      ST(25) = 0
      ST(26) = 0
      ST(27) = 0
      ST(28) = 0
      ST(29) = 0
      ST(30) = 0
      ST(31) = 0
      ST(32) = 0
      ST(33) = 0
      ST(34) = 0
      ST(35) = 0
      ST(36) = t7*ACB(1,3)+t15*ACB(2,3)+t23*ACB(3,3)+t31*ACB(4,3)+t39*AC
     #B(5,3)+t47*ACB(6,3)+t52*ACB(7,3)+t57*ACB(8,3)+t62*ACB(9,3)
      ST(37) = t7*ACB(1,4)+t15*ACB(2,4)+t23*ACB(3,4)+t31*ACB(4,4)+t39*AC
     #B(5,4)+t47*ACB(6,4)+t52*ACB(7,4)+t57*ACB(8,4)+t62*ACB(9,4)
      ST(38) = t7*ACB(1,5)+t15*ACB(2,5)+t23*ACB(3,5)+t31*ACB(4,5)+t39*AC
     #B(5,5)+t47*ACB(6,5)+t52*ACB(7,5)+t57*ACB(8,5)+t62*ACB(9,5)
      ST(39) = 0
      ST(40) = 0
      ST(41) = 0
      ST(42) = t7*ACB(1,9)+t15*ACB(2,9)+t23*ACB(3,9)+t31*ACB(4,9)+t39*AC
     #B(5,9)+t47*ACB(6,9)+t52*ACB(7,9)+t57*ACB(8,9)+t62*ACB(9,9)
      ST(43) = t7*ACB(1,10)+t15*ACB(2,10)+t23*ACB(3,10)+t31*ACB(4,10)+t3
     #9*ACB(5,10)+t47*ACB(6,10)+t52*ACB(7,10)+t57*ACB(8,10)+t62*ACB(9,10
     #)
      ST(44) = t7*ACB(1,11)+t15*ACB(2,11)+t23*ACB(3,11)+t31*ACB(4,11)+t3
     #9*ACB(5,11)+t47*ACB(6,11)+t52*ACB(7,11)+t57*ACB(8,11)+t62*ACB(9,11
     #)
      ST(45) = 0
      ST(46) = 0
      ST(47) = 0
      ST(48) = t7*ACB(1,15)+t15*ACB(2,15)+t23*ACB(3,15)+t31*ACB(4,15)+t3
     #9*ACB(5,15)+t47*ACB(6,15)+t52*ACB(7,15)+t57*ACB(8,15)+t62*ACB(9,15
     #)
      ST(49) = t7*ACB(1,16)+t15*ACB(2,16)+t23*ACB(3,16)+t31*ACB(4,16)+t3
     #9*ACB(5,16)+t47*ACB(6,16)+t52*ACB(7,16)+t57*ACB(8,16)+t62*ACB(9,16
     #)
      ST(50) = t7*ACB(1,17)+t15*ACB(2,17)+t23*ACB(3,17)+t31*ACB(4,17)+t3
     #9*ACB(5,17)+t47*ACB(6,17)+t52*ACB(7,17)+t57*ACB(8,17)+t62*ACB(9,17
     #)
      ST(51) = 0
      ST(52) = t151*ACB(1,4)+t159*ACB(2,4)+t167*ACB(3,4)+t175*ACB(4,4)+t
     #183*ACB(5,4)+t191*ACB(6,4)+t196*ACB(7,4)+t201*ACB(8,4)+t206*ACB(9,
     #4)
      ST(53) = t151*ACB(1,5)+t159*ACB(2,5)+t167*ACB(3,5)+t175*ACB(4,5)+t
     #183*ACB(5,5)+t191*ACB(6,5)+t196*ACB(7,5)+t201*ACB(8,5)+t206*ACB(9,
     #5)
      ST(54) = 0
      ST(55) = 0
      ST(56) = 0
      ST(57) = t151*ACB(1,9)+t159*ACB(2,9)+t167*ACB(3,9)+t175*ACB(4,9)+t
     #183*ACB(5,9)+t191*ACB(6,9)+t196*ACB(7,9)+t201*ACB(8,9)+t206*ACB(9,
     #9)
      ST(58) = t151*ACB(1,10)+t159*ACB(2,10)+t167*ACB(3,10)+t175*ACB(4,1
     #0)+t183*ACB(5,10)+t191*ACB(6,10)+t196*ACB(7,10)+t201*ACB(8,10)+t20
     #6*ACB(9,10)
      ST(59) = t151*ACB(1,11)+t159*ACB(2,11)+t167*ACB(3,11)+t175*ACB(4,1
     #1)+t183*ACB(5,11)+t191*ACB(6,11)+t196*ACB(7,11)+t201*ACB(8,11)+t20
     #6*ACB(9,11)
      ST(60) = 0
      ST(61) = 0
      ST(62) = 0
      ST(63) = t151*ACB(1,15)+t159*ACB(2,15)+t167*ACB(3,15)+t175*ACB(4,1
     #5)+t183*ACB(5,15)+t191*ACB(6,15)+t196*ACB(7,15)+t201*ACB(8,15)+t20
     #6*ACB(9,15)
      ST(64) = t151*ACB(1,16)+t159*ACB(2,16)+t167*ACB(3,16)+t175*ACB(4,1
     #6)+t183*ACB(5,16)+t191*ACB(6,16)+t196*ACB(7,16)+t201*ACB(8,16)+t20
     #6*ACB(9,16)
      ST(65) = t151*ACB(1,17)+t159*ACB(2,17)+t167*ACB(3,17)+t175*ACB(4,1
     #7)+t183*ACB(5,17)+t191*ACB(6,17)+t196*ACB(7,17)+t201*ACB(8,17)+t20
     #6*ACB(9,17)
      ST(66) = 0
      ST(67) = t285*ACB(1,5)+t293*ACB(2,5)+t301*ACB(3,5)+t309*ACB(4,5)+t
     #317*ACB(5,5)+t325*ACB(6,5)+t330*ACB(7,5)+t335*ACB(8,5)+t340*ACB(9,
     #5)
      ST(68) = 0
      ST(69) = 0
      ST(70) = 0
      ST(71) = t285*ACB(1,9)+t293*ACB(2,9)+t301*ACB(3,9)+t309*ACB(4,9)+t
     #317*ACB(5,9)+t325*ACB(6,9)+t330*ACB(7,9)+t335*ACB(8,9)+t340*ACB(9,
     #9)
      ST(72) = t285*ACB(1,10)+t293*ACB(2,10)+t301*ACB(3,10)+t309*ACB(4,1
     #0)+t317*ACB(5,10)+t325*ACB(6,10)+t330*ACB(7,10)+t335*ACB(8,10)+t34
     #0*ACB(9,10)
      ST(73) = t285*ACB(1,11)+t293*ACB(2,11)+t301*ACB(3,11)+t309*ACB(4,1
     #1)+t317*ACB(5,11)+t325*ACB(6,11)+t330*ACB(7,11)+t335*ACB(8,11)+t34
     #0*ACB(9,11)
      ST(74) = 0
      ST(75) = 0
      ST(76) = 0
      ST(77) = t285*ACB(1,15)+t293*ACB(2,15)+t301*ACB(3,15)+t309*ACB(4,1
     #5)+t317*ACB(5,15)+t325*ACB(6,15)+t330*ACB(7,15)+t335*ACB(8,15)+t34
     #0*ACB(9,15)
      ST(78) = t285*ACB(1,16)+t293*ACB(2,16)+t301*ACB(3,16)+t309*ACB(4,1
     #6)+t317*ACB(5,16)+t325*ACB(6,16)+t330*ACB(7,16)+t335*ACB(8,16)+t34
     #0*ACB(9,16)
      ST(79) = t285*ACB(1,17)+t293*ACB(2,17)+t301*ACB(3,17)+t309*ACB(4,1
     #7)+t317*ACB(5,17)+t325*ACB(6,17)+t330*ACB(7,17)+t335*ACB(8,17)+t34
     #0*ACB(9,17)
      ST(80) = 0
      ST(81) = 0
      ST(82) = 0
      ST(83) = 0
      ST(84) = 0
      ST(85) = 0
      ST(86) = 0
      ST(87) = 0
      ST(88) = 0
      ST(89) = 0
      ST(90) = 0
      ST(91) = 0
      ST(92) = 0
      ST(93) = 0
      ST(94) = 0
      ST(95) = 0
      ST(96) = 0
      ST(97) = 0
      ST(98) = 0
      ST(99) = 0
      ST(100) = 0
      ST(101) = 0
      ST(102) = 0
      ST(103) = 0
      ST(104) = 0
      ST(105) = 0
      ST(106) = 0
      ST(107) = 0
      ST(108) = 0
      ST(109) = 0
      ST(110) = 0
      ST(111) = 0
      ST(112) = 0
      ST(113) = 0
      ST(114) = 0
      ST(115) = 0
      ST(116) = 0
      ST(117) = t409*ACB(1,9)+t417*ACB(2,9)+t425*ACB(3,9)+t433*ACB(4,9)+
     #t441*ACB(5,9)+t449*ACB(6,9)+t454*ACB(7,9)+t459*ACB(8,9)+t464*ACB(9
     #,9)
      ST(118) = t409*ACB(1,10)+t417*ACB(2,10)+t425*ACB(3,10)+t433*ACB(4,
     #10)+t441*ACB(5,10)+t449*ACB(6,10)+t454*ACB(7,10)+t459*ACB(8,10)+t4
     #64*ACB(9,10)
      ST(119) = t409*ACB(1,11)+t417*ACB(2,11)+t425*ACB(3,11)+t433*ACB(4,
     #11)+t441*ACB(5,11)+t449*ACB(6,11)+t454*ACB(7,11)+t459*ACB(8,11)+t4
     #64*ACB(9,11)
      ST(120) = 0
      ST(121) = 0
      ST(122) = 0
      ST(123) = t409*ACB(1,15)+t417*ACB(2,15)+t425*ACB(3,15)+t433*ACB(4,
     #15)+t441*ACB(5,15)+t449*ACB(6,15)+t454*ACB(7,15)+t459*ACB(8,15)+t4
     #64*ACB(9,15)
      ST(124) = t409*ACB(1,16)+t417*ACB(2,16)+t425*ACB(3,16)+t433*ACB(4,
     #16)+t441*ACB(5,16)+t449*ACB(6,16)+t454*ACB(7,16)+t459*ACB(8,16)+t4
     #64*ACB(9,16)
      ST(125) = t409*ACB(1,17)+t417*ACB(2,17)+t425*ACB(3,17)+t433*ACB(4,
     #17)+t441*ACB(5,17)+t449*ACB(6,17)+t454*ACB(7,17)+t459*ACB(8,17)+t4
     #64*ACB(9,17)
      ST(126) = 0
      ST(127) = t523*ACB(1,10)+t531*ACB(2,10)+t539*ACB(3,10)+t547*ACB(4,
     #10)+t555*ACB(5,10)+t563*ACB(6,10)+t568*ACB(7,10)+t573*ACB(8,10)+t5
     #78*ACB(9,10)
      ST(128) = t523*ACB(1,11)+t531*ACB(2,11)+t539*ACB(3,11)+t547*ACB(4,
     #11)+t555*ACB(5,11)+t563*ACB(6,11)+t568*ACB(7,11)+t573*ACB(8,11)+t5
     #78*ACB(9,11)
      ST(129) = 0
      ST(130) = 0
      ST(131) = 0
      ST(132) = t523*ACB(1,15)+t531*ACB(2,15)+t539*ACB(3,15)+t547*ACB(4,
     #15)+t555*ACB(5,15)+t563*ACB(6,15)+t568*ACB(7,15)+t573*ACB(8,15)+t5
     #78*ACB(9,15)
      ST(133) = t523*ACB(1,16)+t531*ACB(2,16)+t539*ACB(3,16)+t547*ACB(4,
     #16)+t555*ACB(5,16)+t563*ACB(6,16)+t568*ACB(7,16)+t573*ACB(8,16)+t5
     #78*ACB(9,16)
      ST(134) = t523*ACB(1,17)+t531*ACB(2,17)+t539*ACB(3,17)+t547*ACB(4,
     #17)+t555*ACB(5,17)+t563*ACB(6,17)+t568*ACB(7,17)+t573*ACB(8,17)+t5
     #78*ACB(9,17)
      ST(135) = 0
      ST(136) = t627*ACB(1,11)+t635*ACB(2,11)+t643*ACB(3,11)+t651*ACB(4,
     #11)+t659*ACB(5,11)+t667*ACB(6,11)+t672*ACB(7,11)+t677*ACB(8,11)+t6
     #82*ACB(9,11)
      ST(137) = 0
      ST(138) = 0
      ST(139) = 0
      ST(140) = t627*ACB(1,15)+t635*ACB(2,15)+t643*ACB(3,15)+t651*ACB(4,
     #15)+t659*ACB(5,15)+t667*ACB(6,15)+t672*ACB(7,15)+t677*ACB(8,15)+t6
     #82*ACB(9,15)
      ST(141) = t627*ACB(1,16)+t635*ACB(2,16)+t643*ACB(3,16)+t651*ACB(4,
     #16)+t659*ACB(5,16)+t667*ACB(6,16)+t672*ACB(7,16)+t677*ACB(8,16)+t6
     #82*ACB(9,16)
      ST(142) = t627*ACB(1,17)+t635*ACB(2,17)+t643*ACB(3,17)+t651*ACB(4,
     #17)+t659*ACB(5,17)+t667*ACB(6,17)+t672*ACB(7,17)+t677*ACB(8,17)+t6
     #82*ACB(9,17)
      ST(143) = 0
      ST(144) = 0
      ST(145) = 0
      ST(146) = 0
      ST(147) = 0
      ST(148) = 0
      ST(149) = 0
      ST(150) = 0
      ST(151) = 0
      ST(152) = 0
      ST(153) = 0
      ST(154) = 0
      ST(155) = 0
      ST(156) = 0
      ST(157) = 0
      ST(158) = 0
      ST(159) = 0
      ST(160) = 0
      ST(161) = 0
      ST(162) = t721*ACB(1,15)+t729*ACB(2,15)+t737*ACB(3,15)+t745*ACB(4,
     #15)+t753*ACB(5,15)+t761*ACB(6,15)+t766*ACB(7,15)+t771*ACB(8,15)+t7
     #76*ACB(9,15)
      ST(163) = t721*ACB(1,16)+t729*ACB(2,16)+t737*ACB(3,16)+t745*ACB(4,
     #16)+t753*ACB(5,16)+t761*ACB(6,16)+t766*ACB(7,16)+t771*ACB(8,16)+t7
     #76*ACB(9,16)
      ST(164) = t721*ACB(1,17)+t729*ACB(2,17)+t737*ACB(3,17)+t745*ACB(4,
     #17)+t753*ACB(5,17)+t761*ACB(6,17)+t766*ACB(7,17)+t771*ACB(8,17)+t7
     #76*ACB(9,17)
      ST(165) = 0
      ST(166) = t805*ACB(1,16)+t813*ACB(2,16)+t821*ACB(3,16)+t829*ACB(4,
     #16)+t837*ACB(5,16)+t845*ACB(6,16)+t850*ACB(7,16)+t855*ACB(8,16)+t8
     #60*ACB(9,16)
      ST(167) = t805*ACB(1,17)+t813*ACB(2,17)+t821*ACB(3,17)+t829*ACB(4,
     #17)+t837*ACB(5,17)+t845*ACB(6,17)+t850*ACB(7,17)+t855*ACB(8,17)+t8
     #60*ACB(9,17)
      ST(168) = 0
      s1 = (ACB(1,17)*PDPB(1,1)+ACB(2,17)*PDPB(1,2)+ACB(3,17)*PDPB(1,3)+
     #ACB(4,17)*PDPB(1,4)+ACB(5,17)*PDPB(1,5)+ACB(6,17)*PDPB(1,6))*ACB(1
     #,17)+(ACB(1,17)*PDPB(1,2)+ACB(2,17)*PDPB(2,2)+ACB(3,17)*PDPB(2,3)+
     #ACB(4,17)*PDPB(2,4)+ACB(5,17)*PDPB(2,5)+ACB(6,17)*PDPB(2,6))*ACB(2
     #,17)+(ACB(1,17)*PDPB(1,3)+ACB(2,17)*PDPB(2,3)+ACB(3,17)*PDPB(3,3)+
     #ACB(4,17)*PDPB(3,4)+ACB(5,17)*PDPB(3,5)+ACB(6,17)*PDPB(3,6))*ACB(3
     #,17)+(ACB(1,17)*PDPB(1,4)+ACB(2,17)*PDPB(2,4)+ACB(3,17)*PDPB(3,4)+
     #ACB(4,17)*PDPB(4,4)+ACB(5,17)*PDPB(4,5)+ACB(6,17)*PDPB(4,6))*ACB(4
     #,17)
      ST(169) = s1+(ACB(1,17)*PDPB(1,5)+ACB(2,17)*PDPB(2,5)+ACB(3,17)*PD
     #PB(3,5)+ACB(4,17)*PDPB(4,5)+ACB(5,17)*PDPB(5,5)+ACB(6,17)*PDPB(5,6
     #))*ACB(5,17)+(ACB(1,17)*PDPB(1,6)+ACB(2,17)*PDPB(2,6)+ACB(3,17)*PD
     #PB(3,6)+ACB(4,17)*PDPB(4,6)+ACB(5,17)*PDPB(5,6)+ACB(6,17)*PDPB(6,6
     #))*ACB(6,17)+(ACB(7,17)*PDPB(7,7)+ACB(8,17)*PDPB(7,8)+ACB(9,17)*PD
     #PB(7,9))*ACB(7,17)+(ACB(7,17)*PDPB(7,8)+ACB(8,17)*PDPB(8,8)+ACB(9,
     #17)*PDPB(8,9))*ACB(8,17)+(ACB(7,17)*PDPB(7,9)+ACB(8,17)*PDPB(8,9)+
     #ACB(9,17)*PDPB(9,9))*ACB(9,17)
      ST(170) = 0
      ST(171) = 0

	RETURN
	END
C
C=====================================================================





C======================================================================
C	START 3 NODE SHELL ELEMENT - MEMBRANE SUBROUTINES
C======================================================================
C=====================================================================
      SUBROUTINE TRICM(RS,SLN,RN,SN,CM)
	IMPLICIT REAL*8 (A-H,O-Z)

C	-------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cm FOR MEMBRANE STRAINS
C
C	INPUT VARIABLES
C	RS(2,3)		= ELEMENT COORDINATES IN LOCAL SPACE
C	SLN(3)		= ELEMENT BOUNDARY LENGTHS
C	RN(3),SN(3)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C
C	OUTPUT VARIAVBLES
C	CM(5,18)	= MEMBRANE AND DRILLING STRAIN PARAMETER
C	-------------------------------------------------------
      DIMENSION RS(2,3),SLN(3),RN(3),SN(3)
	DIMENSION CM(7,18)

C	INITIALIZE
	CM = 0.0D0

C	FIRST ROW
C	---------
C	disp u
	CM(1,1)  = (SLN(1)*RN(1)/2+SLN(3)*RN(3)/2)
	CM(1,7)  = (SLN(2)*RN(2)/2+SLN(1)*RN(1)/2)
	CM(1,13) = (SLN(2)*RN(2)/2+SLN(3)*RN(3)/2)
C	rot t
	CM(1,6)  = ((-RS(2,2)/12+RS(2,1)/12)*RN(1)*SLN(1)+
     #            (RS(2,1)/12-RS(2,3)/12)*RN(3)*SLN(3))
	CM(1,12) = ((-RS(2,1)/12+RS(2,2)/12)*RN(1)*SLN(1)+
     #            (RS(2,2)/12-RS(2,3)/12)*RN(2)*SLN(2))
	CM(1,18) = ((-RS(2,2)/12+RS(2,3)/12)*RN(2)*SLN(2)+
     #            (RS(2,3)/12-RS(2,1)/12)*RN(3)*SLN(3))

C	SECOND ROW
C	----------
C	disp u
	CM(2,1)  = 
     #(((RS(1,1)/3+RS(1,2)/6)*RN(1)+(-RS(2,2)/6-RS(2,1)/3)*SN(1))*S
     #LN(1)+((RS(1,1)/3+RS(1,3)/6)*RN(3)+(-RS(2,1)/3-RS(2,3)/6)*SN(3))*S
     #LN(3))
	CM(2,7)  = 
     #(((RS(1,1)/6+RS(1,2)/3)*RN(1)+(-RS(2,2)/3-RS(2,1)/6)*
     #SN(1))*SLN(1)+((RS(1,2)/3+RS(1,3)/6)*RN(2)+(-RS(2,2)/3-RS(2,3)/6)*
     #SN(2))*SLN(2))
	CM(2,13) = 
     #(((RS(1,3)/3+RS(1,2)/6)*RN(2)+(-RS(2,2)/6-RS(
     #2,3)/3)*SN(2))*SLN(2)+((RS(1,3)/3+RS(1,1)/6)*RN(3)+(-RS(2,3)/3-RS(
     #2,1)/6)*SN(3))*SLN(3))
C	disp v
	CM(2,2)  = 
     #((-RS(2,2)/6-RS(2,1)/3)*RN(1)*SLN(1)+
     #(-RS(2,1)/3-RS(2,3)/6)*RN(3)*SLN(3))
	CM(2,8) = 
     #((-RS(2,2)/3-RS(2,1)/6)*RN(1)*SLN(1)+(-RS(2,2)/3-RS(2,3)/6
     #)*RN(2)*SLN(2))
	CM(2,14) = 
     #((-RS(2,2)/6-RS(2,3)/3)*RN(2)*SLN(2)+(-RS(2,
     #3)/3-RS(2,1)/6)*RN(3)*SLN(3))
C	rot t
	CM(2,6)  = 
     #(((-RS(2,2)*RS(1,2)/12+RS(2,1)*RS(1,1)/12)*RN(1)+(RS(2,2)*
     #*2/24-RS(2,1)**2/24)*SN(1))*SLN(1)+((-RS(2,3)*RS(1,3)/12+RS(2,1)*R
     #S(1,1)/12)*RN(3)+(-RS(2,1)**2/24+RS(2,3)**2/24)*SN(3))*SLN(3))
	CM(2,12) = 
     #(((RS(2,2)*RS(1,2)/12-RS(2,1)*RS(1,1)/12)*RN(1)+(-RS(2,2)*
     #*2/24+RS(2,1)**2/24)*SN(1))*SLN(1)+((RS(2,2)*RS(1,2)/12-RS(2,3)*RS
     #(1,3)/12)*RN(2)+(-RS(2,2)**2/24+RS(2,3)**2/24)*SN(2))*SLN(2))
	CM(2,18) = 
     #(((RS(2,3)*RS(1,3)/12-RS(2,2)*RS(1,2)/12)*RN(2)+(RS(2,2)**2/24-
     #RS(2,3)**2/24)*SN(2))*SLN(2)+((RS(2,3)*RS(1,3)/12-RS(2,1)*RS(1,1)/
     #12)*RN(3)+(RS(2,1)**2/24-RS(2,3)**2/24)*SN(3))*SLN(3))

C	THIRD ROW
C	---------
C	disp u
	CM(3,1)  = 
     #((RS(2,2)/6+RS(2,1)/3)*RN(1)*SLN(1)+(RS(2,3)/6+RS(2,1)/3)*RN(
     #3)*SLN(3))
	CM(3,7)  = 
     #((RS(2,2)/3+RS(2,1)/6)*RN(1)*SLN(1)+(RS(2,2)/3+RS
     #(2,3)/6)*RN(2)*SLN(2))
	CM(3,13) = 
     #((RS(2,2)/6+RS(2,3)/3)*RN(2)*SLN(2)+(
     #RS(2,3)/3+RS(2,1)/6)*RN(3)*SLN(3))
C	rot t
	CM(3,6)  = 
     #((-RS(2,2)**2/24+RS(2,1)*
     #*2/24)*RN(1)*SLN(1)+(RS(2,1)**2/24-RS(2,3)**2/24)*RN(3)*SLN(3))
	CM(3,12) = 
     #((RS(2,2)**2/24-RS(2,1)**2/24)*RN(1)*SLN(1)+(RS(2,2)**2/24-RS
     #(2,3)**2/24)*RN(2)*SLN(2))
	CM(3,18) = 
     #((-RS(2,2)**2/24+RS(2,3)**2/24)*
     #RN(2)*SLN(2)+(-RS(2,1)**2/24+RS(2,3)**2/24)*RN(3)*SLN(3))

C	FOURTH ROW
C	----------
C	disp v
	CM(4,2)  = (SLN(1)*SN(1)/2+SLN(3)*SN(3)/2)
	CM(4,8)  = (SLN(2)*SN(2)/2+SLN(1)*SN(1)/2)
	CM(4,14) = (SLN(2)*SN(2)/2+SLN(3)*SN(3)/2)
C	rot t
	CM(4,6)  = 
     #((RS(1,2)/12-
     #RS(1,1)/12)*SN(1)*SLN(1)+(-RS(1,1)/12+RS(1,3)/12)*SN(3)*SLN(3))
	CM(4,12) = 
     #((RS(1,1)/12-RS(1,2)/12)*SN(1)*SLN(1)+(-RS(1,2)/12+RS(1,3)/12
     #)*SN(2)*SLN(2))
	CM(4,18) = 
     #((RS(1,2)/12-RS(1,3)/12)*SN(2)*SLN(2)+(-RS(
     #1,3)/12+RS(1,1)/12)*SN(3)*SLN(3))

C	FIFTH ROW
C	---------
C	disp v
	CM(5,2)  = 
     #((RS(1,1)/3+RS(1,2)/6)*SN(1)*SLN(1)+(RS(1,1)/3+RS(1,3)/6)*SN(
     #3)*SLN(3))
	CM(5,8)  = 
     #((RS(1,1)/6+RS(1,2)/3)*SN(1)*SLN(1)+(RS(1,2)/3+RS
     #(1,3)/6)*SN(2)*SLN(2))
	CM(5,14) = 
     #((RS(1,3)/3+RS(1,2)/6)*SN(2)*SLN(2)+(
     #RS(1,3)/3+RS(1,1)/6)*SN(3)*SLN(3))
C	rot t
	CM(5,6)  = 
     #((RS(1,2)**2/24-RS(1,1)**
     #2/24)*SN(1)*SLN(1)+(RS(1,3)**2/24-RS(1,1)**2/24)*SN(3)*SLN(3))
	CM(5,12) = 
     #((RS(1,1)**2/24-RS(1,2)**2/24)*SN(1)*SLN(1)+(RS(1,3)**2/24-RS(
     #1,2)**2/24)*SN(2)*SLN(2))
	CM(5,18) = 
     #((RS(1,2)**2/24-RS(1,3)**2/24)*SN
     #(2)*SLN(2)+(-RS(1,3)**2/24+RS(1,1)**2/24)*SN(3)*SLN(3))

C	SIXTH ROW
C	---------
C	disp u
	CM(6,1)  = 
     #((-RS(1,2)/6-RS(1,1)/3)*SN(1)*SLN(1)+(-RS(1,1)/3-RS(1,3)/6)*S
     #N(3)*SLN(3))
	CM(6,7)  = 
     #((-RS(1,1)/6-RS(1,2)/3)*SN(1)*SLN(1)+(-RS(1,3)/
     #6-RS(1,2)/3)*SN(2)*SLN(2))
	CM(6,13) = 
     #((-RS(1,3)/3-RS(1,2)/6)*SN(2)*SLN
     #(2)+(-RS(1,3)/3-RS(1,1)/6)*SN(3)*SLN(3))
C	disp v
	CM(6,2)  = 
     #(((-RS(1,2)/6-RS(1,
     #1)/3)*RN(1)+(RS(2,2)/6+RS(2,1)/3)*SN(1))*SLN(1)+((-RS(1,1)/3-RS(1,
     #3)/6)*RN(3)+(RS(2,3)/6+RS(2,1)/3)*SN(3))*SLN(3))
	CM(6,8)  = 
     #(((-RS(1,1)/6-RS(1,2)/3)*RN(1)+(RS(2,2)/3+RS(2,1)/6)*SN(1)
     #)*SLN(1)+((-RS(1,3)/6-RS(1,2)/3)*RN(2)+(RS(2,2)/3+RS(2,3)/6)*SN(2)
     #)*SLN(2))
	CM(6,14) = 
     #(((-RS(1,3)/3-RS(1,2)/6)*RN(2)+(RS(2,2)/6+RS(2,3)/
     #3)*SN(2))*SLN(2)+((-RS(1,3)/3-RS(1,1)/6)*RN(3)+(RS(2,3)/3+RS(2,1)/
     #6)*SN(3))*SLN(3))
C	rot t
	CM(6,6)  = 
     #(((RS(1,1)**2/24-RS(1,2)**2/24)*RN(1)+(RS(2,2)*RS(1,2)/12-
     #RS(2,1)*RS(1,1)/12)*SN(1))*SLN(1)+((-RS(1,3)**2/24+RS(1,1)**2/24)*
     #RN(3)+(RS(2,3)*RS(1,3)/12-RS(2,1)*RS(1,1)/12)*SN(3))*SLN(3))
	CM(6,12) = 
     #(((RS(1,2)**2/24-RS(1,1)**2/24)*RN(1)+(-RS(2,2)*RS(1,2)/12
     #+RS(2,1)*RS(1,1)/12)*SN(1))*SLN(1)+((RS(1,2)**2/24-RS(1,3)**2/24)*
     #RN(2)+(RS(2,3)*RS(1,3)/12-RS(2,2)*RS(1,2)/12)*SN(2))*SLN(2))
	CM(6,18) = 
     #(((RS(1,3)**2/24-RS(1,2)**2/24)*RN(2)+(RS(2,2)*RS(1,2)/12-RS(2,3
     #)*RS(1,3)/12)*SN(2))*SLN(2)+((RS(1,3)**2/24-RS(1,1)**2/24)*RN(3)+(
     #-RS(2,3)*RS(1,3)/12+RS(2,1)*RS(1,1)/12)*SN(3))*SLN(3))

C	SEVENTH ROW
C	-----------
C	disp u
	CM(7,1)  = (SLN(1)*SN(1)/2+SLN(3)*SN(3)/2)
	CM(7,7)  = (SLN(2)*SN(2)/2+SLN(1)*SN(1)/2)
	CM(7,13) = (SLN(2)*SN(2)/2+SLN(3)*SN(3)/2)
C	disp v
	CM(7,2)  = (SLN(1)*RN(1)/2+SLN(3)*RN(3)/2)
	CM(7,8)  = (SLN(2)*RN(2)/2+SLN(1)*RN(1)/2)
	CM(7,14) = (SLN(2)*RN(2)/2+SLN(3)*RN(3)/2)
C	rot t
	CM(7,6)  = 
     #(((RS(1,2)/12-RS(1,1)/12)*RN(1)+(-RS(2,2)/12+RS(2,1)/12)*S
     #N(1))*SLN(1)+((-RS(1,1)/12+RS(1,3)/12)*RN(3)+(RS(2,1)/12-RS(2,3)/1
     #2)*SN(3))*SLN(3))
	CM(7,12) = 
     #(((RS(1,1)/12-RS(1,2)/12)*RN(1)+(-RS(2,1)
     #/12+RS(2,2)/12)*SN(1))*SLN(1)+((-RS(1,2)/12+RS(1,3)/12)*RN(2)+(RS(
     #2,2)/12-RS(2,3)/12)*SN(2))*SLN(2))
	CM(7,18) = 
     #(((RS(1,2)/12-RS(1,3)/12
     #)*RN(2)+(-RS(2,2)/12+RS(2,3)/12)*SN(2))*SLN(2)+((-RS(1,3)/12+RS(1,
     #1)/12)*RN(3)+(RS(2,3)/12-RS(2,1)/12)*SN(3))*SLN(3))

c	ak = 0.0d0
c	do i = 1,7
c		cm(i,6)  = cm(i,6)*ak
c		cm(i,12) = cm(i,12)*ak
c		cm(i,18) = cm(i,18)*ak
c	enddo

	RETURN
	END
C
C=====================================================================
      SUBROUTINE TABCM(AB1,CM,ACM)

	IMPLICIT REAL*8 (A-H,O-Z)
C	------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR INVERSE(Ab)*(Cb)=ACB
C
C	INPUT VARIABLES
C	AB1(3,3)	= STRAIN PARAMETERS
C	CM(3,18)	= MEMBRANE AND DRILLING STRAIN PARAMETER
C
C	LOCAL VARIABLE
C	AMI(5,5)	= INVERSE OF Am
C
C	OUTPUT VARIAVBLES
C	ACM(5,24)	= INVERSE(Am)*(Cm)
C	------------------------------------------------------------
	DIMENSION AB1(3,3),CM(7,18)
	DIMENSION AM(7,7),AMI(7,7)
	DIMENSION ACM(7,18)

C	------------
C	SOLVE FOR AM
C	------------
	AM = 0.0D0

	AM(1,1) =  AB1(1,1)
	AM(1,2) =  AB1(1,2)
	AM(1,3) =  AB1(1,3)

	AM(2,2) =  AB1(2,2)+AB1(3,3)
	AM(2,3) =  AB1(2,3)
	AM(2,6) =  AB1(2,3)
	AM(2,7) = -AB1(1,3)

	AM(3,3) =  AB1(3,3)

	AM(4,4) =  AB1(1,1)
	AM(4,5) =  AB1(1,2)
	AM(4,6) =  AB1(1,3)

	AM(5,5) =  AB1(2,2)
	AM(5,6) =  AB1(2,3)

	AM(6,6) =  AB1(2,2)+AB1(3,3)
	AM(6,7) = -AB1(1,2)

	AM(7,7) =  AB1(1,1)

	DO 10 I = 1,7
	  DO 20 J = I,7
	    IF (J.NE.I) AM(J,I) = AM(I,J)
20	  CONTINUE
10	CONTINUE

C	----------------------
C	SOLVE FOR INVERSE(AMI)
C	----------------------
	CALL  INVMAT(AM,AMI,7)

C	--------------------------------
C	SOLVE FOR INVERSE(Am)*(Cm) - ACM
C	--------------------------------
      ACM(1,1) = AMI(1,1)*CM(1,1)+AMI(1,2)*CM(2,1)+AMI(1,3)*CM(3,1)+AMI(
     #1,6)*CM(6,1)+AMI(1,7)*CM(7,1)
      ACM(1,2) = AMI(1,2)*CM(2,2)+AMI(1,4)*CM(4,2)+AMI(1,5)*CM(5,2)+AMI(
     #1,6)*CM(6,2)+AMI(1,7)*CM(7,2)
      ACM(1,3) = 0.0D0
      ACM(1,4) = 0.0D0
      ACM(1,5) = 0.0D0
      ACM(1,6) = AMI(1,1)*CM(1,6)+AMI(1,2)*CM(2,6)+AMI(1,3)*CM(3,6)+AMI(
     #1,4)*CM(4,6)+AMI(1,5)*CM(5,6)+AMI(1,6)*CM(6,6)+AMI(1,7)*CM(7,6)
      ACM(1,7) = AMI(1,1)*CM(1,7)+AMI(1,2)*CM(2,7)+AMI(1,3)*CM(3,7)+AMI(
     #1,6)*CM(6,7)+AMI(1,7)*CM(7,7)
      ACM(1,8) = AMI(1,2)*CM(2,8)+AMI(1,4)*CM(4,8)+AMI(1,5)*CM(5,8)+AMI(
     #1,6)*CM(6,8)+AMI(1,7)*CM(7,8)
      ACM(1,9) = 0.0D0
      ACM(1,10) = 0.0D0
      ACM(1,11) = 0.0D0
      ACM(1,12) = AMI(1,1)*CM(1,12)+AMI(1,2)*CM(2,12)+AMI(1,3)*CM(3,12)+
     #AMI(1,4)*CM(4,12)+AMI(1,5)*CM(5,12)+AMI(1,6)*CM(6,12)+AMI(1,7)*CM(
     #7,12)
      ACM(1,13) = AMI(1,1)*CM(1,13)+AMI(1,2)*CM(2,13)+AMI(1,3)*CM(3,13)+
     #AMI(1,6)*CM(6,13)+AMI(1,7)*CM(7,13)
      ACM(1,14) = AMI(1,2)*CM(2,14)+AMI(1,4)*CM(4,14)+AMI(1,5)*CM(5,14)+
     #AMI(1,6)*CM(6,14)+AMI(1,7)*CM(7,14)
      ACM(1,15) = 0.0D0
      ACM(1,16) = 0.0D0
      ACM(1,17) = 0.0D0
      ACM(1,18) = AMI(1,1)*CM(1,18)+AMI(1,2)*CM(2,18)+AMI(1,3)*CM(3,18)+
     #AMI(1,4)*CM(4,18)+AMI(1,5)*CM(5,18)+AMI(1,6)*CM(6,18)+AMI(1,7)*CM(
     #7,18)
      ACM(2,1) = AMI(2,1)*CM(1,1)+AMI(2,2)*CM(2,1)+AMI(2,3)*CM(3,1)+AMI(
     #2,6)*CM(6,1)+AMI(2,7)*CM(7,1)
      ACM(2,2) = AMI(2,2)*CM(2,2)+AMI(2,4)*CM(4,2)+AMI(2,5)*CM(5,2)+AMI(
     #2,6)*CM(6,2)+AMI(2,7)*CM(7,2)
      ACM(2,3) = 0.0D0
      ACM(2,4) = 0.0D0
      ACM(2,5) = 0.0D0
      ACM(2,6) = AMI(2,1)*CM(1,6)+AMI(2,2)*CM(2,6)+AMI(2,3)*CM(3,6)+AMI(
     #2,4)*CM(4,6)+AMI(2,5)*CM(5,6)+AMI(2,6)*CM(6,6)+AMI(2,7)*CM(7,6)
      ACM(2,7) = AMI(2,1)*CM(1,7)+AMI(2,2)*CM(2,7)+AMI(2,3)*CM(3,7)+AMI(
     #2,6)*CM(6,7)+AMI(2,7)*CM(7,7)
      ACM(2,8) = AMI(2,2)*CM(2,8)+AMI(2,4)*CM(4,8)+AMI(2,5)*CM(5,8)+AMI(
     #2,6)*CM(6,8)+AMI(2,7)*CM(7,8)
      ACM(2,9) = 0.0D0
      ACM(2,10) = 0.0D0
      ACM(2,11) = 0.0D0
      ACM(2,12) = AMI(2,1)*CM(1,12)+AMI(2,2)*CM(2,12)+AMI(2,3)*CM(3,12)+
     #AMI(2,4)*CM(4,12)+AMI(2,5)*CM(5,12)+AMI(2,6)*CM(6,12)+AMI(2,7)*CM(
     #7,12)
      ACM(2,13) = AMI(2,1)*CM(1,13)+AMI(2,2)*CM(2,13)+AMI(2,3)*CM(3,13)+
     #AMI(2,6)*CM(6,13)+AMI(2,7)*CM(7,13)
      ACM(2,14) = AMI(2,2)*CM(2,14)+AMI(2,4)*CM(4,14)+AMI(2,5)*CM(5,14)+
     #AMI(2,6)*CM(6,14)+AMI(2,7)*CM(7,14)
      ACM(2,15) = 0.0D0
      ACM(2,16) = 0.0D0
      ACM(2,17) = 0.0D0
      ACM(2,18) = AMI(2,1)*CM(1,18)+AMI(2,2)*CM(2,18)+AMI(2,3)*CM(3,18)+
     #AMI(2,4)*CM(4,18)+AMI(2,5)*CM(5,18)+AMI(2,6)*CM(6,18)+AMI(2,7)*CM(
     #7,18)
      ACM(3,1) = AMI(3,1)*CM(1,1)+AMI(3,2)*CM(2,1)+AMI(3,3)*CM(3,1)+AMI(
     #3,6)*CM(6,1)+AMI(3,7)*CM(7,1)
      ACM(3,2) = AMI(3,2)*CM(2,2)+AMI(3,4)*CM(4,2)+AMI(3,5)*CM(5,2)+AMI(
     #3,6)*CM(6,2)+AMI(3,7)*CM(7,2)
      ACM(3,3) = 0.0D0
      ACM(3,4) = 0.0D0
      ACM(3,5) = 0.0D0
      ACM(3,6) = AMI(3,1)*CM(1,6)+AMI(3,2)*CM(2,6)+AMI(3,3)*CM(3,6)+AMI(
     #3,4)*CM(4,6)+AMI(3,5)*CM(5,6)+AMI(3,6)*CM(6,6)+AMI(3,7)*CM(7,6)
      ACM(3,7) = AMI(3,1)*CM(1,7)+AMI(3,2)*CM(2,7)+AMI(3,3)*CM(3,7)+AMI(
     #3,6)*CM(6,7)+AMI(3,7)*CM(7,7)
      ACM(3,8) = AMI(3,2)*CM(2,8)+AMI(3,4)*CM(4,8)+AMI(3,5)*CM(5,8)+AMI(
     #3,6)*CM(6,8)+AMI(3,7)*CM(7,8)
      ACM(3,9) = 0.0D0
      ACM(3,10) = 0.0D0
      ACM(3,11) = 0.0D0
      ACM(3,12) = AMI(3,1)*CM(1,12)+AMI(3,2)*CM(2,12)+AMI(3,3)*CM(3,12)+
     #AMI(3,4)*CM(4,12)+AMI(3,5)*CM(5,12)+AMI(3,6)*CM(6,12)+AMI(3,7)*CM(
     #7,12)
      ACM(3,13) = AMI(3,1)*CM(1,13)+AMI(3,2)*CM(2,13)+AMI(3,3)*CM(3,13)+
     #AMI(3,6)*CM(6,13)+AMI(3,7)*CM(7,13)
      ACM(3,14) = AMI(3,2)*CM(2,14)+AMI(3,4)*CM(4,14)+AMI(3,5)*CM(5,14)+
     #AMI(3,6)*CM(6,14)+AMI(3,7)*CM(7,14)
      ACM(3,15) = 0.0D0
      ACM(3,16) = 0.0D0
      ACM(3,17) = 0.0D0
      ACM(3,18) = AMI(3,1)*CM(1,18)+AMI(3,2)*CM(2,18)+AMI(3,3)*CM(3,18)+
     #AMI(3,4)*CM(4,18)+AMI(3,5)*CM(5,18)+AMI(3,6)*CM(6,18)+AMI(3,7)*CM(
     #7,18)
      ACM(4,1) = AMI(4,1)*CM(1,1)+AMI(4,2)*CM(2,1)+AMI(4,3)*CM(3,1)+AMI(
     #4,6)*CM(6,1)+AMI(4,7)*CM(7,1)
      ACM(4,2) = AMI(4,2)*CM(2,2)+AMI(4,4)*CM(4,2)+AMI(4,5)*CM(5,2)+AMI(
     #4,6)*CM(6,2)+AMI(4,7)*CM(7,2)
      ACM(4,3) = 0.0D0
      ACM(4,4) = 0.0D0
      ACM(4,5) = 0.0D0
      ACM(4,6) = AMI(4,1)*CM(1,6)+AMI(4,2)*CM(2,6)+AMI(4,3)*CM(3,6)+AMI(
     #4,4)*CM(4,6)+AMI(4,5)*CM(5,6)+AMI(4,6)*CM(6,6)+AMI(4,7)*CM(7,6)
      ACM(4,7) = AMI(4,1)*CM(1,7)+AMI(4,2)*CM(2,7)+AMI(4,3)*CM(3,7)+AMI(
     #4,6)*CM(6,7)+AMI(4,7)*CM(7,7)
      ACM(4,8) = AMI(4,2)*CM(2,8)+AMI(4,4)*CM(4,8)+AMI(4,5)*CM(5,8)+AMI(
     #4,6)*CM(6,8)+AMI(4,7)*CM(7,8)
      ACM(4,9) = 0.0D0
      ACM(4,10) = 0.0D0
      ACM(4,11) = 0.0D0
      ACM(4,12) = AMI(4,1)*CM(1,12)+AMI(4,2)*CM(2,12)+AMI(4,3)*CM(3,12)+
     #AMI(4,4)*CM(4,12)+AMI(4,5)*CM(5,12)+AMI(4,6)*CM(6,12)+AMI(4,7)*CM(
     #7,12)
      ACM(4,13) = AMI(4,1)*CM(1,13)+AMI(4,2)*CM(2,13)+AMI(4,3)*CM(3,13)+
     #AMI(4,6)*CM(6,13)+AMI(4,7)*CM(7,13)
      ACM(4,14) = AMI(4,2)*CM(2,14)+AMI(4,4)*CM(4,14)+AMI(4,5)*CM(5,14)+
     #AMI(4,6)*CM(6,14)+AMI(4,7)*CM(7,14)
      ACM(4,15) = 0.0D0
      ACM(4,16) = 0.0D0
      ACM(4,17) = 0.0D0
      ACM(4,18) = AMI(4,1)*CM(1,18)+AMI(4,2)*CM(2,18)+AMI(4,3)*CM(3,18)+
     #AMI(4,4)*CM(4,18)+AMI(4,5)*CM(5,18)+AMI(4,6)*CM(6,18)+AMI(4,7)*CM(
     #7,18)
      ACM(5,1) = AMI(5,1)*CM(1,1)+AMI(5,2)*CM(2,1)+AMI(5,3)*CM(3,1)+AMI(
     #5,6)*CM(6,1)+AMI(5,7)*CM(7,1)
      ACM(5,2) = AMI(5,2)*CM(2,2)+AMI(5,4)*CM(4,2)+AMI(5,5)*CM(5,2)+AMI(
     #5,6)*CM(6,2)+AMI(5,7)*CM(7,2)
      ACM(5,3) = 0.0D0
      ACM(5,4) = 0.0D0
      ACM(5,5) = 0.0D0
      ACM(5,6) = AMI(5,1)*CM(1,6)+AMI(5,2)*CM(2,6)+AMI(5,3)*CM(3,6)+AMI(
     #5,4)*CM(4,6)+AMI(5,5)*CM(5,6)+AMI(5,6)*CM(6,6)+AMI(5,7)*CM(7,6)
      ACM(5,7) = AMI(5,1)*CM(1,7)+AMI(5,2)*CM(2,7)+AMI(5,3)*CM(3,7)+AMI(
     #5,6)*CM(6,7)+AMI(5,7)*CM(7,7)
      ACM(5,8) = AMI(5,2)*CM(2,8)+AMI(5,4)*CM(4,8)+AMI(5,5)*CM(5,8)+AMI(
     #5,6)*CM(6,8)+AMI(5,7)*CM(7,8)
      ACM(5,9) = 0.0D0
      ACM(5,10) = 0.0D0
      ACM(5,11) = 0.0D0
      ACM(5,12) = AMI(5,1)*CM(1,12)+AMI(5,2)*CM(2,12)+AMI(5,3)*CM(3,12)+
     #AMI(5,4)*CM(4,12)+AMI(5,5)*CM(5,12)+AMI(5,6)*CM(6,12)+AMI(5,7)*CM(
     #7,12)
      ACM(5,13) = AMI(5,1)*CM(1,13)+AMI(5,2)*CM(2,13)+AMI(5,3)*CM(3,13)+
     #AMI(5,6)*CM(6,13)+AMI(5,7)*CM(7,13)
      ACM(5,14) = AMI(5,2)*CM(2,14)+AMI(5,4)*CM(4,14)+AMI(5,5)*CM(5,14)+
     #AMI(5,6)*CM(6,14)+AMI(5,7)*CM(7,14)
      ACM(5,15) = 0.0D0
      ACM(5,16) = 0.0D0
      ACM(5,17) = 0.0D0
      ACM(5,18) = AMI(5,1)*CM(1,18)+AMI(5,2)*CM(2,18)+AMI(5,3)*CM(3,18)+
     #AMI(5,4)*CM(4,18)+AMI(5,5)*CM(5,18)+AMI(5,6)*CM(6,18)+AMI(5,7)*CM(
     #7,18)
      ACM(6,1) = AMI(6,1)*CM(1,1)+AMI(6,2)*CM(2,1)+AMI(6,3)*CM(3,1)+AMI(
     #6,6)*CM(6,1)+AMI(6,7)*CM(7,1)
      ACM(6,2) = AMI(6,2)*CM(2,2)+AMI(6,4)*CM(4,2)+AMI(6,5)*CM(5,2)+AMI(
     #6,6)*CM(6,2)+AMI(6,7)*CM(7,2)
      ACM(6,3) = 0.0D0
      ACM(6,4) = 0.0D0
      ACM(6,5) = 0.0D0
      ACM(6,6) = AMI(6,1)*CM(1,6)+AMI(6,2)*CM(2,6)+AMI(6,3)*CM(3,6)+AMI(
     #6,4)*CM(4,6)+AMI(6,5)*CM(5,6)+AMI(6,6)*CM(6,6)+AMI(6,7)*CM(7,6)
      ACM(6,7) = AMI(6,1)*CM(1,7)+AMI(6,2)*CM(2,7)+AMI(6,3)*CM(3,7)+AMI(
     #6,6)*CM(6,7)+AMI(6,7)*CM(7,7)
      ACM(6,8) = AMI(6,2)*CM(2,8)+AMI(6,4)*CM(4,8)+AMI(6,5)*CM(5,8)+AMI(
     #6,6)*CM(6,8)+AMI(6,7)*CM(7,8)
      ACM(6,9) = 0.0D0
      ACM(6,10) = 0.0D0
      ACM(6,11) = 0.0D0
      ACM(6,12) = AMI(6,1)*CM(1,12)+AMI(6,2)*CM(2,12)+AMI(6,3)*CM(3,12)+
     #AMI(6,4)*CM(4,12)+AMI(6,5)*CM(5,12)+AMI(6,6)*CM(6,12)+AMI(6,7)*CM(
     #7,12)
      ACM(6,13) = AMI(6,1)*CM(1,13)+AMI(6,2)*CM(2,13)+AMI(6,3)*CM(3,13)+
     #AMI(6,6)*CM(6,13)+AMI(6,7)*CM(7,13)
      ACM(6,14) = AMI(6,2)*CM(2,14)+AMI(6,4)*CM(4,14)+AMI(6,5)*CM(5,14)+
     #AMI(6,6)*CM(6,14)+AMI(6,7)*CM(7,14)
      ACM(6,15) = 0.0D0
      ACM(6,16) = 0.0D0
      ACM(6,17) = 0.0D0
      ACM(6,18) = AMI(6,1)*CM(1,18)+AMI(6,2)*CM(2,18)+AMI(6,3)*CM(3,18)+
     #AMI(6,4)*CM(4,18)+AMI(6,5)*CM(5,18)+AMI(6,6)*CM(6,18)+AMI(6,7)*CM(
     #7,18)
      ACM(7,1) = AMI(7,1)*CM(1,1)+AMI(7,2)*CM(2,1)+AMI(7,3)*CM(3,1)+AMI(
     #7,6)*CM(6,1)+AMI(7,7)*CM(7,1)
      ACM(7,2) = AMI(7,2)*CM(2,2)+AMI(7,4)*CM(4,2)+AMI(7,5)*CM(5,2)+AMI(
     #7,6)*CM(6,2)+AMI(7,7)*CM(7,2)
      ACM(7,3) = 0.0D0
      ACM(7,4) = 0.0D0
      ACM(7,5) = 0.0D0
      ACM(7,6) = AMI(7,1)*CM(1,6)+AMI(7,2)*CM(2,6)+AMI(7,3)*CM(3,6)+AMI(
     #7,4)*CM(4,6)+AMI(7,5)*CM(5,6)+AMI(7,6)*CM(6,6)+AMI(7,7)*CM(7,6)
      ACM(7,7) = AMI(7,1)*CM(1,7)+AMI(7,2)*CM(2,7)+AMI(7,3)*CM(3,7)+AMI(
     #7,6)*CM(6,7)+AMI(7,7)*CM(7,7)
      ACM(7,8) = AMI(7,2)*CM(2,8)+AMI(7,4)*CM(4,8)+AMI(7,5)*CM(5,8)+AMI(
     #7,6)*CM(6,8)+AMI(7,7)*CM(7,8)
      ACM(7,9) = 0.0D0
      ACM(7,10) = 0.0D0
      ACM(7,11) = 0.0D0
      ACM(7,12) = AMI(7,1)*CM(1,12)+AMI(7,2)*CM(2,12)+AMI(7,3)*CM(3,12)+
     #AMI(7,4)*CM(4,12)+AMI(7,5)*CM(5,12)+AMI(7,6)*CM(6,12)+AMI(7,7)*CM(
     #7,12)
      ACM(7,13) = AMI(7,1)*CM(1,13)+AMI(7,2)*CM(2,13)+AMI(7,3)*CM(3,13)+
     #AMI(7,6)*CM(6,13)+AMI(7,7)*CM(7,13)
      ACM(7,14) = AMI(7,2)*CM(2,14)+AMI(7,4)*CM(4,14)+AMI(7,5)*CM(5,14)+
     #AMI(7,6)*CM(6,14)+AMI(7,7)*CM(7,14)
      ACM(7,15) = 0.0D0
      ACM(7,16) = 0.0D0
      ACM(7,17) = 0.0D0
      ACM(7,18) = AMI(7,1)*CM(1,18)+AMI(7,2)*CM(2,18)+AMI(7,3)*CM(3,18)+
     #AMI(7,4)*CM(4,18)+AMI(7,5)*CM(5,18)+AMI(7,6)*CM(6,18)+AMI(7,7)*CM(
     #7,18)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE NMSPM18(AB1,DR,PSP)
	IMPLICIT REAL*8 (A-H,O-Z)
C	---------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Nm*S*Pm - PSP
C
C	INPUT VARIABLES
C	AREA		= AREA OF TRIANGULAR ELEMENT
C	DR(64)		= MEMBERANE,BENDING AND SHEAR STIFFNESS
C
C	OUTPUT VARIAVBLES
C	PSP(3,3)	= Nm*S*Pm
C	---------------------------------------------------
      DIMENSION DR(64),AB1(3,3)
	DIMENSION PSP(7,7)

C	AREA INTEGRATION RESULTS
C	------------------------
	A = AB1(1,1) !1
	R = AB1(1,2) !R
	S = AB1(1,3) !S
	RS = AB1(2,3) !R*S
	R2 = AB1(2,2) !R**2
	S2 = AB1(3,3) !S**2
C	INT(TRANS(P),DR,P)
C	------------------
      PSP(1,1) = DR(1)*A
      PSP(1,2) = R*DR(1)
      PSP(1,3) = S*DR(1)
      PSP(1,4) = DR(2)*A
      PSP(1,5) = R*DR(2)
      PSP(1,6) = S*DR(2)
      PSP(1,7) = 0.0D0

      PSP(2,1) = R*DR(1)
      PSP(2,2) = R2*DR(1)+S2*DR(19)
      PSP(2,3) = RS*DR(1)
      PSP(2,4) = R*DR(2)
      PSP(2,5) = R2*DR(2)
      PSP(2,6) = RS*(DR(2)+DR(19))
      PSP(2,7) = -S*DR(19)

      PSP(3,1) = S*DR(1)
      PSP(3,2) = RS*DR(1)
      PSP(3,3) = S2*DR(1)
      PSP(3,4) = S*DR(2)
      PSP(3,5) = RS*DR(2)
      PSP(3,6) = S2*DR(2)
      PSP(3,7) = 0.0D0

      PSP(4,1) = DR(2)*A
      PSP(4,2) = R*DR(2)
      PSP(4,3) = S*DR(2)
      PSP(4,4) = DR(10)*A
      PSP(4,5) = R*DR(10)
      PSP(4,6) = S*DR(10)
      PSP(4,7) = 0.0D0

      PSP(5,1) = R*DR(2)
      PSP(5,2) = R2*DR(2)
      PSP(5,3) = RS*DR(2)
      PSP(5,4) = R*DR(10)
      PSP(5,5) = R2*DR(10)
      PSP(5,6) = RS*DR(10)
      PSP(5,7) = 0.0D0

      PSP(6,1) = S*DR(2)
      PSP(6,2) = RS*(DR(2)+DR(19))
      PSP(6,3) = S2*DR(2)
      PSP(6,4) = S*DR(10)
      PSP(6,5) = RS*DR(10)
      PSP(6,6) = S2*DR(10)+R2*DR(19)
      PSP(6,7) = -R*DR(19)

      PSP(7,1) = 0.0D0
      PSP(7,2) = -S*DR(19)
      PSP(7,3) = 0.0D0
      PSP(7,4) = 0.0D0
      PSP(7,5) = 0.0D0
      PSP(7,6) = -R*DR(19)
      PSP(7,7) = DR(19)*A

	RETURN
	END
C
C=====================================================================
      SUBROUTINE APHAM18(ACM,TEDIS,APM)
	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR PARAMETER ALPHA
C
C	INPUT VARIABLES
C	ACM(5,18)	= INVERSE(Am)*(Cm)
C	TEDIS(18)	= TRANSPOSE(T)*REDIS
C
C	OUTPUT VARIAVBLES
C	APM(5)		= MEMBRANE STRAIN PARAMETER
C	--------------------------------------------------------
      DIMENSION ACM(7,18),TEDIS(18)
	DIMENSION APM(7)

      APM(1) = ACM(1,1)*TEDIS(1)+ACM(1,2)*TEDIS(2)+ACM(1,6)*TEDIS(6)+ACM
     #(1,7)*TEDIS(7)+ACM(1,8)*TEDIS(8)+ACM(1,12)*TEDIS(12)+ACM(1,13)*TED
     #IS(13)+ACM(1,14)*TEDIS(14)+ACM(1,18)*TEDIS(18)

      APM(2) = ACM(2,1)*TEDIS(1)+ACM(2,2)*TEDIS(2)+ACM(2,6)*TEDIS(6)+ACM
     #(2,7)*TEDIS(7)+ACM(2,8)*TEDIS(8)+ACM(2,12)*TEDIS(12)+ACM(2,13)*TED
     #IS(13)+ACM(2,14)*TEDIS(14)+ACM(2,18)*TEDIS(18)

      APM(3) = ACM(3,1)*TEDIS(1)+ACM(3,2)*TEDIS(2)+ACM(3,6)*TEDIS(6)+ACM
     #(3,7)*TEDIS(7)+ACM(3,8)*TEDIS(8)+ACM(3,12)*TEDIS(12)+ACM(3,13)*TED
     #IS(13)+ACM(3,14)*TEDIS(14)+ACM(3,18)*TEDIS(18)

      APM(4) = ACM(4,1)*TEDIS(1)+ACM(4,2)*TEDIS(2)+ACM(4,6)*TEDIS(6)+ACM
     #(4,7)*TEDIS(7)+ACM(4,8)*TEDIS(8)+ACM(4,12)*TEDIS(12)+ACM(4,13)*TED
     #IS(13)+ACM(4,14)*TEDIS(14)+ACM(4,18)*TEDIS(18)

      APM(5) = ACM(5,1)*TEDIS(1)+ACM(5,2)*TEDIS(2)+ACM(5,6)*TEDIS(6)+ACM
     #(5,7)*TEDIS(7)+ACM(5,8)*TEDIS(8)+ACM(5,12)*TEDIS(12)+ACM(5,13)*TED
     #IS(13)+ACM(5,14)*TEDIS(14)+ACM(5,18)*TEDIS(18)

      APM(6) = ACM(6,1)*TEDIS(1)+ACM(6,2)*TEDIS(2)+ACM(6,6)*TEDIS(6)+ACM
     #(6,7)*TEDIS(7)+ACM(6,8)*TEDIS(8)+ACM(6,12)*TEDIS(12)+ACM(6,13)*TED
     #IS(13)+ACM(6,14)*TEDIS(14)+ACM(6,18)*TEDIS(18)

      APM(7) = ACM(7,1)*TEDIS(1)+ACM(7,2)*TEDIS(2)+ACM(7,6)*TEDIS(6)+ACM
     #(7,7)*TEDIS(7)+ACM(7,8)*TEDIS(8)+ACM(7,12)*TEDIS(12)+ACM(7,13)*TED
     #IS(13)+ACM(7,14)*TEDIS(14)+ACM(7,18)*TEDIS(18)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE KMEM18(ACM,PSP,ST)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-----------------------------------------------------------
C	PURPOSE:	TO ADD THE FLEXURAL CONTRIBUTION TO THE
C					ELEMENT STIFFNESS MATRIX
C
C	INPUT VARIABLES
C	ACM(5,18)	= INVERSE(Am)*(Cm)
C	PSP(5,5)	= Nm*S*Pm
C
C	OUTPUT VARIAVBLES
C	S(171)		= STIFFNESS MATRIX IN ROW FORMAT
C	-----------------------------------------------------------
      DIMENSION ACM(7,18),PSP(7,7)
	DIMENSION ST(171)

      t8 = ACM(1,1)*PSP(1,1)+ACM(2,1)*PSP(1,2)+ACM(3,1)*PSP(1,3)+ACM(4,1
     #)*PSP(1,4)+ACM(5,1)*PSP(1,5)+ACM(6,1)*PSP(1,6)+ACM(7,1)*PSP(1,7)
      t17 = ACM(1,1)*PSP(1,2)+ACM(2,1)*PSP(2,2)+ACM(3,1)*PSP(2,3)+ACM(4,
     #1)*PSP(2,4)+ACM(5,1)*PSP(2,5)+ACM(6,1)*PSP(2,6)+ACM(7,1)*PSP(2,7)
      t26 = ACM(1,1)*PSP(1,3)+ACM(2,1)*PSP(2,3)+ACM(3,1)*PSP(3,3)+ACM(4,
     #1)*PSP(3,4)+ACM(5,1)*PSP(3,5)+ACM(6,1)*PSP(3,6)+ACM(7,1)*PSP(3,7)
      t35 = ACM(1,1)*PSP(1,4)+ACM(2,1)*PSP(2,4)+ACM(3,1)*PSP(3,4)+ACM(4,
     #1)*PSP(4,4)+ACM(5,1)*PSP(4,5)+ACM(6,1)*PSP(4,6)+ACM(7,1)*PSP(4,7)
      t44 = ACM(1,1)*PSP(1,5)+ACM(2,1)*PSP(2,5)+ACM(3,1)*PSP(3,5)+ACM(4,
     #1)*PSP(4,5)+ACM(5,1)*PSP(5,5)+ACM(6,1)*PSP(5,6)+ACM(7,1)*PSP(5,7)
      t53 = ACM(1,1)*PSP(1,6)+ACM(2,1)*PSP(2,6)+ACM(3,1)*PSP(3,6)+ACM(4,
     #1)*PSP(4,6)+ACM(5,1)*PSP(5,6)+ACM(6,1)*PSP(6,6)+ACM(7,1)*PSP(6,7)
      t62 = ACM(1,1)*PSP(1,7)+ACM(2,1)*PSP(2,7)+ACM(3,1)*PSP(3,7)+ACM(4,
     #1)*PSP(4,7)+ACM(5,1)*PSP(5,7)+ACM(6,1)*PSP(6,7)+ACM(7,1)*PSP(7,7)
      t136 = ACM(1,2)*PSP(1,1)+ACM(2,2)*PSP(1,2)+ACM(3,2)*PSP(1,3)+ACM(4
     #,2)*PSP(1,4)+ACM(5,2)*PSP(1,5)+ACM(6,2)*PSP(1,6)+ACM(7,2)*PSP(1,7)
      t145 = ACM(1,2)*PSP(1,2)+ACM(2,2)*PSP(2,2)+ACM(3,2)*PSP(2,3)+ACM(4
     #,2)*PSP(2,4)+ACM(5,2)*PSP(2,5)+ACM(6,2)*PSP(2,6)+ACM(7,2)*PSP(2,7)
      t154 = ACM(1,2)*PSP(1,3)+ACM(2,2)*PSP(2,3)+ACM(3,2)*PSP(3,3)+ACM(4
     #,2)*PSP(3,4)+ACM(5,2)*PSP(3,5)+ACM(6,2)*PSP(3,6)+ACM(7,2)*PSP(3,7)
      t163 = ACM(1,2)*PSP(1,4)+ACM(2,2)*PSP(2,4)+ACM(3,2)*PSP(3,4)+ACM(4
     #,2)*PSP(4,4)+ACM(5,2)*PSP(4,5)+ACM(6,2)*PSP(4,6)+ACM(7,2)*PSP(4,7)
      t172 = ACM(1,2)*PSP(1,5)+ACM(2,2)*PSP(2,5)+ACM(3,2)*PSP(3,5)+ACM(4
     #,2)*PSP(4,5)+ACM(5,2)*PSP(5,5)+ACM(6,2)*PSP(5,6)+ACM(7,2)*PSP(5,7)
      t181 = ACM(1,2)*PSP(1,6)+ACM(2,2)*PSP(2,6)+ACM(3,2)*PSP(3,6)+ACM(4
     #,2)*PSP(4,6)+ACM(5,2)*PSP(5,6)+ACM(6,2)*PSP(6,6)+ACM(7,2)*PSP(6,7)
      t190 = ACM(1,2)*PSP(1,7)+ACM(2,2)*PSP(2,7)+ACM(3,2)*PSP(3,7)+ACM(4
     #,2)*PSP(4,7)+ACM(5,2)*PSP(5,7)+ACM(6,2)*PSP(6,7)+ACM(7,2)*PSP(7,7)
      t256 = ACM(1,6)*PSP(1,1)+ACM(2,6)*PSP(1,2)+ACM(3,6)*PSP(1,3)+ACM(4
     #,6)*PSP(1,4)+ACM(5,6)*PSP(1,5)+ACM(6,6)*PSP(1,6)+ACM(7,6)*PSP(1,7)
      t265 = ACM(1,6)*PSP(1,2)+ACM(2,6)*PSP(2,2)+ACM(3,6)*PSP(2,3)+ACM(4
     #,6)*PSP(2,4)+ACM(5,6)*PSP(2,5)+ACM(6,6)*PSP(2,6)+ACM(7,6)*PSP(2,7)
      t274 = ACM(1,6)*PSP(1,3)+ACM(2,6)*PSP(2,3)+ACM(3,6)*PSP(3,3)+ACM(4
     #,6)*PSP(3,4)+ACM(5,6)*PSP(3,5)+ACM(6,6)*PSP(3,6)+ACM(7,6)*PSP(3,7)
      t283 = ACM(1,6)*PSP(1,4)+ACM(2,6)*PSP(2,4)+ACM(3,6)*PSP(3,4)+ACM(4
     #,6)*PSP(4,4)+ACM(5,6)*PSP(4,5)+ACM(6,6)*PSP(4,6)+ACM(7,6)*PSP(4,7)
      t292 = ACM(1,6)*PSP(1,5)+ACM(2,6)*PSP(2,5)+ACM(3,6)*PSP(3,5)+ACM(4
     #,6)*PSP(4,5)+ACM(5,6)*PSP(5,5)+ACM(6,6)*PSP(5,6)+ACM(7,6)*PSP(5,7)
      t301 = ACM(1,6)*PSP(1,6)+ACM(2,6)*PSP(2,6)+ACM(3,6)*PSP(3,6)+ACM(4
     #,6)*PSP(4,6)+ACM(5,6)*PSP(5,6)+ACM(6,6)*PSP(6,6)+ACM(7,6)*PSP(6,7)
      t310 = ACM(1,6)*PSP(1,7)+ACM(2,6)*PSP(2,7)+ACM(3,6)*PSP(3,7)+ACM(4
     #,6)*PSP(4,7)+ACM(5,6)*PSP(5,7)+ACM(6,6)*PSP(6,7)+ACM(7,6)*PSP(7,7)
      t368 = ACM(1,7)*PSP(1,1)+ACM(2,7)*PSP(1,2)+ACM(3,7)*PSP(1,3)+ACM(4
     #,7)*PSP(1,4)+ACM(5,7)*PSP(1,5)+ACM(6,7)*PSP(1,6)+ACM(7,7)*PSP(1,7)
      t377 = ACM(1,7)*PSP(1,2)+ACM(2,7)*PSP(2,2)+ACM(3,7)*PSP(2,3)+ACM(4
     #,7)*PSP(2,4)+ACM(5,7)*PSP(2,5)+ACM(6,7)*PSP(2,6)+ACM(7,7)*PSP(2,7)
      t386 = ACM(1,7)*PSP(1,3)+ACM(2,7)*PSP(2,3)+ACM(3,7)*PSP(3,3)+ACM(4
     #,7)*PSP(3,4)+ACM(5,7)*PSP(3,5)+ACM(6,7)*PSP(3,6)+ACM(7,7)*PSP(3,7)
      t395 = ACM(1,7)*PSP(1,4)+ACM(2,7)*PSP(2,4)+ACM(3,7)*PSP(3,4)+ACM(4
     #,7)*PSP(4,4)+ACM(5,7)*PSP(4,5)+ACM(6,7)*PSP(4,6)+ACM(7,7)*PSP(4,7)
      t404 = ACM(1,7)*PSP(1,5)+ACM(2,7)*PSP(2,5)+ACM(3,7)*PSP(3,5)+ACM(4
     #,7)*PSP(4,5)+ACM(5,7)*PSP(5,5)+ACM(6,7)*PSP(5,6)+ACM(7,7)*PSP(5,7)
      t413 = ACM(1,7)*PSP(1,6)+ACM(2,7)*PSP(2,6)+ACM(3,7)*PSP(3,6)+ACM(4
     #,7)*PSP(4,6)+ACM(5,7)*PSP(5,6)+ACM(6,7)*PSP(6,6)+ACM(7,7)*PSP(6,7)
      t422 = ACM(1,7)*PSP(1,7)+ACM(2,7)*PSP(2,7)+ACM(3,7)*PSP(3,7)+ACM(4
     #,7)*PSP(4,7)+ACM(5,7)*PSP(5,7)+ACM(6,7)*PSP(6,7)+ACM(7,7)*PSP(7,7)
      t472 = ACM(1,8)*PSP(1,1)+ACM(2,8)*PSP(1,2)+ACM(3,8)*PSP(1,3)+ACM(4
     #,8)*PSP(1,4)+ACM(5,8)*PSP(1,5)+ACM(6,8)*PSP(1,6)+ACM(7,8)*PSP(1,7)
      t481 = ACM(1,8)*PSP(1,2)+ACM(2,8)*PSP(2,2)+ACM(3,8)*PSP(2,3)+ACM(4
     #,8)*PSP(2,4)+ACM(5,8)*PSP(2,5)+ACM(6,8)*PSP(2,6)+ACM(7,8)*PSP(2,7)
      t490 = ACM(1,8)*PSP(1,3)+ACM(2,8)*PSP(2,3)+ACM(3,8)*PSP(3,3)+ACM(4
     #,8)*PSP(3,4)+ACM(5,8)*PSP(3,5)+ACM(6,8)*PSP(3,6)+ACM(7,8)*PSP(3,7)
      t499 = ACM(1,8)*PSP(1,4)+ACM(2,8)*PSP(2,4)+ACM(3,8)*PSP(3,4)+ACM(4
     #,8)*PSP(4,4)+ACM(5,8)*PSP(4,5)+ACM(6,8)*PSP(4,6)+ACM(7,8)*PSP(4,7)
      t508 = ACM(1,8)*PSP(1,5)+ACM(2,8)*PSP(2,5)+ACM(3,8)*PSP(3,5)+ACM(4
     #,8)*PSP(4,5)+ACM(5,8)*PSP(5,5)+ACM(6,8)*PSP(5,6)+ACM(7,8)*PSP(5,7)
      t517 = ACM(1,8)*PSP(1,6)+ACM(2,8)*PSP(2,6)+ACM(3,8)*PSP(3,6)+ACM(4
     #,8)*PSP(4,6)+ACM(5,8)*PSP(5,6)+ACM(6,8)*PSP(6,6)+ACM(7,8)*PSP(6,7)
      t526 = ACM(1,8)*PSP(1,7)+ACM(2,8)*PSP(2,7)+ACM(3,8)*PSP(3,7)+ACM(4
     #,8)*PSP(4,7)+ACM(5,8)*PSP(5,7)+ACM(6,8)*PSP(6,7)+ACM(7,8)*PSP(7,7)
      t568 = ACM(1,12)*PSP(1,1)+ACM(2,12)*PSP(1,2)+ACM(3,12)*PSP(1,3)+AC
     #M(4,12)*PSP(1,4)+ACM(5,12)*PSP(1,5)+ACM(6,12)*PSP(1,6)+ACM(7,12)*P
     #SP(1,7)
      t577 = ACM(1,12)*PSP(1,2)+ACM(2,12)*PSP(2,2)+ACM(3,12)*PSP(2,3)+AC
     #M(4,12)*PSP(2,4)+ACM(5,12)*PSP(2,5)+ACM(6,12)*PSP(2,6)+ACM(7,12)*P
     #SP(2,7)
      t586 = ACM(1,12)*PSP(1,3)+ACM(2,12)*PSP(2,3)+ACM(3,12)*PSP(3,3)+AC
     #M(4,12)*PSP(3,4)+ACM(5,12)*PSP(3,5)+ACM(6,12)*PSP(3,6)+ACM(7,12)*P
     #SP(3,7)
      t595 = ACM(1,12)*PSP(1,4)+ACM(2,12)*PSP(2,4)+ACM(3,12)*PSP(3,4)+AC
     #M(4,12)*PSP(4,4)+ACM(5,12)*PSP(4,5)+ACM(6,12)*PSP(4,6)+ACM(7,12)*P
     #SP(4,7)
      t604 = ACM(1,12)*PSP(1,5)+ACM(2,12)*PSP(2,5)+ACM(3,12)*PSP(3,5)+AC
     #M(4,12)*PSP(4,5)+ACM(5,12)*PSP(5,5)+ACM(6,12)*PSP(5,6)+ACM(7,12)*P
     #SP(5,7)
      t613 = ACM(1,12)*PSP(1,6)+ACM(2,12)*PSP(2,6)+ACM(3,12)*PSP(3,6)+AC
     #M(4,12)*PSP(4,6)+ACM(5,12)*PSP(5,6)+ACM(6,12)*PSP(6,6)+ACM(7,12)*P
     #SP(6,7)
      t622 = ACM(1,12)*PSP(1,7)+ACM(2,12)*PSP(2,7)+ACM(3,12)*PSP(3,7)+AC
     #M(4,12)*PSP(4,7)+ACM(5,12)*PSP(5,7)+ACM(6,12)*PSP(6,7)+ACM(7,12)*P
     #SP(7,7)
      t656 = ACM(1,13)*PSP(1,1)+ACM(2,13)*PSP(1,2)+ACM(3,13)*PSP(1,3)+AC
     #M(4,13)*PSP(1,4)+ACM(5,13)*PSP(1,5)+ACM(6,13)*PSP(1,6)+ACM(7,13)*P
     #SP(1,7)
      t665 = ACM(1,13)*PSP(1,2)+ACM(2,13)*PSP(2,2)+ACM(3,13)*PSP(2,3)+AC
     #M(4,13)*PSP(2,4)+ACM(5,13)*PSP(2,5)+ACM(6,13)*PSP(2,6)+ACM(7,13)*P
     #SP(2,7)
      t674 = ACM(1,13)*PSP(1,3)+ACM(2,13)*PSP(2,3)+ACM(3,13)*PSP(3,3)+AC
     #M(4,13)*PSP(3,4)+ACM(5,13)*PSP(3,5)+ACM(6,13)*PSP(3,6)+ACM(7,13)*P
     #SP(3,7)
      t683 = ACM(1,13)*PSP(1,4)+ACM(2,13)*PSP(2,4)+ACM(3,13)*PSP(3,4)+AC
     #M(4,13)*PSP(4,4)+ACM(5,13)*PSP(4,5)+ACM(6,13)*PSP(4,6)+ACM(7,13)*P
     #SP(4,7)
      t692 = ACM(1,13)*PSP(1,5)+ACM(2,13)*PSP(2,5)+ACM(3,13)*PSP(3,5)+AC
     #M(4,13)*PSP(4,5)+ACM(5,13)*PSP(5,5)+ACM(6,13)*PSP(5,6)+ACM(7,13)*P
     #SP(5,7)
      t701 = ACM(1,13)*PSP(1,6)+ACM(2,13)*PSP(2,6)+ACM(3,13)*PSP(3,6)+AC
     #M(4,13)*PSP(4,6)+ACM(5,13)*PSP(5,6)+ACM(6,13)*PSP(6,6)+ACM(7,13)*P
     #SP(6,7)
      t710 = ACM(1,13)*PSP(1,7)+ACM(2,13)*PSP(2,7)+ACM(3,13)*PSP(3,7)+AC
     #M(4,13)*PSP(4,7)+ACM(5,13)*PSP(5,7)+ACM(6,13)*PSP(6,7)+ACM(7,13)*P
     #SP(7,7)
      t736 = ACM(1,14)*PSP(1,1)+ACM(2,14)*PSP(1,2)+ACM(3,14)*PSP(1,3)+AC
     #M(4,14)*PSP(1,4)+ACM(5,14)*PSP(1,5)+ACM(6,14)*PSP(1,6)+ACM(7,14)*P
     #SP(1,7)
      t745 = ACM(1,14)*PSP(1,2)+ACM(2,14)*PSP(2,2)+ACM(3,14)*PSP(2,3)+AC
     #M(4,14)*PSP(2,4)+ACM(5,14)*PSP(2,5)+ACM(6,14)*PSP(2,6)+ACM(7,14)*P
     #SP(2,7)
      t754 = ACM(1,14)*PSP(1,3)+ACM(2,14)*PSP(2,3)+ACM(3,14)*PSP(3,3)+AC
     #M(4,14)*PSP(3,4)+ACM(5,14)*PSP(3,5)+ACM(6,14)*PSP(3,6)+ACM(7,14)*P
     #SP(3,7)
      t763 = ACM(1,14)*PSP(1,4)+ACM(2,14)*PSP(2,4)+ACM(3,14)*PSP(3,4)+AC
     #M(4,14)*PSP(4,4)+ACM(5,14)*PSP(4,5)+ACM(6,14)*PSP(4,6)+ACM(7,14)*P
     #SP(4,7)
      t772 = ACM(1,14)*PSP(1,5)+ACM(2,14)*PSP(2,5)+ACM(3,14)*PSP(3,5)+AC
     #M(4,14)*PSP(4,5)+ACM(5,14)*PSP(5,5)+ACM(6,14)*PSP(5,6)+ACM(7,14)*P
     #SP(5,7)
      t781 = ACM(1,14)*PSP(1,6)+ACM(2,14)*PSP(2,6)+ACM(3,14)*PSP(3,6)+AC
     #M(4,14)*PSP(4,6)+ACM(5,14)*PSP(5,6)+ACM(6,14)*PSP(6,6)+ACM(7,14)*P
     #SP(6,7)
      t790 = ACM(1,14)*PSP(1,7)+ACM(2,14)*PSP(2,7)+ACM(3,14)*PSP(3,7)+AC
     #M(4,14)*PSP(4,7)+ACM(5,14)*PSP(5,7)+ACM(6,14)*PSP(6,7)+ACM(7,14)*P
     #SP(7,7)
      ST(1) = ST(1)+t8*ACM(1,1)+t17*ACM(2,1)+t26*ACM(3,1)+t35*ACM(4,1)+t
     #44*ACM(5,1)+t53*ACM(6,1)+t62*ACM(7,1)
      ST(2) = ST(2)+t8*ACM(1,2)+t17*ACM(2,2)+t26*ACM(3,2)+t35*ACM(4,2)+t
     #44*ACM(5,2)+t53*ACM(6,2)+t62*ACM(7,2)

      ST(6) = ST(6)+t8*ACM(1,6)+t17*ACM(2,6)+t26*ACM(3,6)+t35*ACM(4,6)+t
     #44*ACM(5,6)+t53*ACM(6,6)+t62*ACM(7,6)
      ST(7) = ST(7)+t8*ACM(1,7)+t17*ACM(2,7)+t26*ACM(3,7)+t35*ACM(4,7)+t
     #44*ACM(5,7)+t53*ACM(6,7)+t62*ACM(7,7)
      ST(8) = ST(8)+t8*ACM(1,8)+t17*ACM(2,8)+t26*ACM(3,8)+t35*ACM(4,8)+t
     #44*ACM(5,8)+t53*ACM(6,8)+t62*ACM(7,8)

      ST(12) = ST(12)+t8*ACM(1,12)+t17*ACM(2,12)+t26*ACM(3,12)+t35*ACM(4
     #,12)+t44*ACM(5,12)+t53*ACM(6,12)+t62*ACM(7,12)
      ST(13) = ST(13)+t8*ACM(1,13)+t17*ACM(2,13)+t26*ACM(3,13)+t35*ACM(4
     #,13)+t44*ACM(5,13)+t53*ACM(6,13)+t62*ACM(7,13)
      ST(14) = ST(14)+t8*ACM(1,14)+t17*ACM(2,14)+t26*ACM(3,14)+t35*ACM(4
     #,14)+t44*ACM(5,14)+t53*ACM(6,14)+t62*ACM(7,14)

      ST(18) = ST(18)+t8*ACM(1,18)+t17*ACM(2,18)+t26*ACM(3,18)+t35*ACM(4
     #,18)+t44*ACM(5,18)+t53*ACM(6,18)+t62*ACM(7,18)
      ST(19) = ST(19)+t136*ACM(1,2)+t145*ACM(2,2)+t154*ACM(3,2)+t163*ACM
     #(4,2)+t172*ACM(5,2)+t181*ACM(6,2)+t190*ACM(7,2)

      ST(23) = ST(23)+t136*ACM(1,6)+t145*ACM(2,6)+t154*ACM(3,6)+t163*ACM
     #(4,6)+t172*ACM(5,6)+t181*ACM(6,6)+t190*ACM(7,6)
      ST(24) = ST(24)+t136*ACM(1,7)+t145*ACM(2,7)+t154*ACM(3,7)+t163*ACM
     #(4,7)+t172*ACM(5,7)+t181*ACM(6,7)+t190*ACM(7,7)
      ST(25) = ST(25)+t136*ACM(1,8)+t145*ACM(2,8)+t154*ACM(3,8)+t163*ACM
     #(4,8)+t172*ACM(5,8)+t181*ACM(6,8)+t190*ACM(7,8)

      ST(29) = ST(29)+t136*ACM(1,12)+t145*ACM(2,12)+t154*ACM(3,12)+t163*
     #ACM(4,12)+t172*ACM(5,12)+t181*ACM(6,12)+t190*ACM(7,12)
      ST(30) = ST(30)+t136*ACM(1,13)+t145*ACM(2,13)+t154*ACM(3,13)+t163*
     #ACM(4,13)+t172*ACM(5,13)+t181*ACM(6,13)+t190*ACM(7,13)
      ST(31) = ST(31)+t136*ACM(1,14)+t145*ACM(2,14)+t154*ACM(3,14)+t163*
     #ACM(4,14)+t172*ACM(5,14)+t181*ACM(6,14)+t190*ACM(7,14)

      ST(35) = ST(35)+t136*ACM(1,18)+t145*ACM(2,18)+t154*ACM(3,18)+t163*
     #ACM(4,18)+t172*ACM(5,18)+t181*ACM(6,18)+t190*ACM(7,18)

      ST(81) = ST(81)+t256*ACM(1,6)+t265*ACM(2,6)+t274*ACM(3,6)+t283*ACM
     #(4,6)+t292*ACM(5,6)+t301*ACM(6,6)+t310*ACM(7,6)
      ST(82) = ST(82)+t256*ACM(1,7)+t265*ACM(2,7)+t274*ACM(3,7)+t283*ACM
     #(4,7)+t292*ACM(5,7)+t301*ACM(6,7)+t310*ACM(7,7)
      ST(83) = ST(83)+t256*ACM(1,8)+t265*ACM(2,8)+t274*ACM(3,8)+t283*ACM
     #(4,8)+t292*ACM(5,8)+t301*ACM(6,8)+t310*ACM(7,8)

      ST(87) = ST(87)+t256*ACM(1,12)+t265*ACM(2,12)+t274*ACM(3,12)+t283*
     #ACM(4,12)+t292*ACM(5,12)+t301*ACM(6,12)+t310*ACM(7,12)
      ST(88) = ST(88)+t256*ACM(1,13)+t265*ACM(2,13)+t274*ACM(3,13)+t283*
     #ACM(4,13)+t292*ACM(5,13)+t301*ACM(6,13)+t310*ACM(7,13)
      ST(89) = ST(89)+t256*ACM(1,14)+t265*ACM(2,14)+t274*ACM(3,14)+t283*
     #ACM(4,14)+t292*ACM(5,14)+t301*ACM(6,14)+t310*ACM(7,14)

      ST(93) = ST(93)+t256*ACM(1,18)+t265*ACM(2,18)+t274*ACM(3,18)+t283*
     #ACM(4,18)+t292*ACM(5,18)+t301*ACM(6,18)+t310*ACM(7,18)
      ST(94) = ST(94)+t368*ACM(1,7)+t377*ACM(2,7)+t386*ACM(3,7)+t395*ACM
     #(4,7)+t404*ACM(5,7)+t413*ACM(6,7)+t422*ACM(7,7)
      ST(95) = ST(95)+t368*ACM(1,8)+t377*ACM(2,8)+t386*ACM(3,8)+t395*ACM
     #(4,8)+t404*ACM(5,8)+t413*ACM(6,8)+t422*ACM(7,8)

      ST(99) = ST(99)+t368*ACM(1,12)+t377*ACM(2,12)+t386*ACM(3,12)+t395*
     #ACM(4,12)+t404*ACM(5,12)+t413*ACM(6,12)+t422*ACM(7,12)
      ST(100) = ST(100)+t368*ACM(1,13)+t377*ACM(2,13)+t386*ACM(3,13)+t39
     #5*ACM(4,13)+t404*ACM(5,13)+t413*ACM(6,13)+t422*ACM(7,13)
      ST(101) = ST(101)+t368*ACM(1,14)+t377*ACM(2,14)+t386*ACM(3,14)+t39
     #5*ACM(4,14)+t404*ACM(5,14)+t413*ACM(6,14)+t422*ACM(7,14)

      ST(105) = ST(105)+t368*ACM(1,18)+t377*ACM(2,18)+t386*ACM(3,18)+t39
     #5*ACM(4,18)+t404*ACM(5,18)+t413*ACM(6,18)+t422*ACM(7,18)
      ST(106) = ST(106)+t472*ACM(1,8)+t481*ACM(2,8)+t490*ACM(3,8)+t499*A
     #CM(4,8)+t508*ACM(5,8)+t517*ACM(6,8)+t526*ACM(7,8)

      ST(110) = ST(110)+t472*ACM(1,12)+t481*ACM(2,12)+t490*ACM(3,12)+t49
     #9*ACM(4,12)+t508*ACM(5,12)+t517*ACM(6,12)+t526*ACM(7,12)
      ST(111) = ST(111)+t472*ACM(1,13)+t481*ACM(2,13)+t490*ACM(3,13)+t49
     #9*ACM(4,13)+t508*ACM(5,13)+t517*ACM(6,13)+t526*ACM(7,13)
      ST(112) = ST(112)+t472*ACM(1,14)+t481*ACM(2,14)+t490*ACM(3,14)+t49
     #9*ACM(4,14)+t508*ACM(5,14)+t517*ACM(6,14)+t526*ACM(7,14)

      ST(116) = ST(116)+t472*ACM(1,18)+t481*ACM(2,18)+t490*ACM(3,18)+t49
     #9*ACM(4,18)+t508*ACM(5,18)+t517*ACM(6,18)+t526*ACM(7,18)

      ST(144) = ST(144)+t568*ACM(1,12)+t577*ACM(2,12)+t586*ACM(3,12)+t59
     #5*ACM(4,12)+t604*ACM(5,12)+t613*ACM(6,12)+t622*ACM(7,12)
      ST(145) = ST(145)+t568*ACM(1,13)+t577*ACM(2,13)+t586*ACM(3,13)+t59
     #5*ACM(4,13)+t604*ACM(5,13)+t613*ACM(6,13)+t622*ACM(7,13)
      ST(146) = ST(146)+t568*ACM(1,14)+t577*ACM(2,14)+t586*ACM(3,14)+t59
     #5*ACM(4,14)+t604*ACM(5,14)+t613*ACM(6,14)+t622*ACM(7,14)

      ST(150) = ST(150)+t568*ACM(1,18)+t577*ACM(2,18)+t586*ACM(3,18)+t59
     #5*ACM(4,18)+t604*ACM(5,18)+t613*ACM(6,18)+t622*ACM(7,18)
      ST(151) = ST(151)+t656*ACM(1,13)+t665*ACM(2,13)+t674*ACM(3,13)+t68
     #3*ACM(4,13)+t692*ACM(5,13)+t701*ACM(6,13)+t710*ACM(7,13)
      ST(152) = ST(152)+t656*ACM(1,14)+t665*ACM(2,14)+t674*ACM(3,14)+t68
     #3*ACM(4,14)+t692*ACM(5,14)+t701*ACM(6,14)+t710*ACM(7,14)

      ST(156) = ST(156)+t656*ACM(1,18)+t665*ACM(2,18)+t674*ACM(3,18)+t68
     #3*ACM(4,18)+t692*ACM(5,18)+t701*ACM(6,18)+t710*ACM(7,18)
      ST(157) = ST(157)+t736*ACM(1,14)+t745*ACM(2,14)+t754*ACM(3,14)+t76
     #3*ACM(4,14)+t772*ACM(5,14)+t781*ACM(6,14)+t790*ACM(7,14)

      ST(161) = ST(161)+t736*ACM(1,18)+t745*ACM(2,18)+t754*ACM(3,18)+t76
     #3*ACM(4,18)+t772*ACM(5,18)+t781*ACM(6,18)+t790*ACM(7,18)

      s1 = ST(171)+(ACM(1,18)*PSP(1,1)+ACM(2,18)*PSP(1,2)+ACM(3,18)*PSP(
     #1,3)+ACM(4,18)*PSP(1,4)+ACM(5,18)*PSP(1,5)+ACM(6,18)*PSP(1,6)+ACM(
     #7,18)*PSP(1,7))*ACM(1,18)+(ACM(1,18)*PSP(1,2)+ACM(2,18)*PSP(2,2)+A
     #CM(3,18)*PSP(2,3)+ACM(4,18)*PSP(2,4)+ACM(5,18)*PSP(2,5)+ACM(6,18)*
     #PSP(2,6)+ACM(7,18)*PSP(2,7))*ACM(2,18)+(ACM(1,18)*PSP(1,3)+ACM(2,1
     #8)*PSP(2,3)+ACM(3,18)*PSP(3,3)+ACM(4,18)*PSP(3,4)+ACM(5,18)*PSP(3,
     #5)+ACM(6,18)*PSP(3,6)+ACM(7,18)*PSP(3,7))*ACM(3,18)
      s2 = s1+(ACM(1,18)*PSP(1,4)+ACM(2,18)*PSP(2,4)+ACM(3,18)*PSP(3,4)+
     #ACM(4,18)*PSP(4,4)+ACM(5,18)*PSP(4,5)+ACM(6,18)*PSP(4,6)+ACM(7,18)
     #*PSP(4,7))*ACM(4,18)
      ST(171) = s2+(ACM(1,18)*PSP(1,5)+ACM(2,18)*PSP(2,5)+ACM(3,18)*PSP(
     #3,5)+ACM(4,18)*PSP(4,5)+ACM(5,18)*PSP(5,5)+ACM(6,18)*PSP(5,6)+ACM(
     #7,18)*PSP(5,7))*ACM(5,18)+(ACM(1,18)*PSP(1,6)+ACM(2,18)*PSP(2,6)+A
     #CM(3,18)*PSP(3,6)+ACM(4,18)*PSP(4,6)+ACM(5,18)*PSP(5,6)+ACM(6,18)*
     #PSP(6,6)+ACM(7,18)*PSP(6,7))*ACM(6,18)+(ACM(1,18)*PSP(1,7)+ACM(2,1
     #8)*PSP(2,7)+ACM(3,18)*PSP(3,7)+ACM(4,18)*PSP(4,7)+ACM(5,18)*PSP(5,
     #7)+ACM(6,18)*PSP(6,7)+ACM(7,18)*PSP(7,7))*ACM(7,18)

	RETURN
	END
C
C=====================================================================





C======================================================================
C	START 3 NODE SHELL ELEMENT - TRANSVERSE SHEAR SUBROUTINES
C======================================================================
C=====================================================================
      SUBROUTINE TRICS(SLN,RN,SN,FLR,FLS,DWR,DWS,CSC,BLR,BLS)

	IMPLICIT REAL*8 (A-H,O-Z)
C	----------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cs FOR TRANSVERSE SHEAR STRAINS
C
C	INPUT VARIABLES
C	BR,BS		= SHEAR DEFORMATION PARAMETER
C	SLN(3)		= ELEMENT BOUNDARY LENGTHS
C	RN(4),SN(3)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	FLR,FLS
C	DWR,DWS		= COMMON FACTORS DEFINED IN SUBROUTINE COMFACT
C
C	OUTPUT VARIAVBLES
C	CS(2,18)	= MEMBRANE AND DRILLING STRAIN PARAMETER
C	CSC(4,18)	= COMPONENTS OF CS
C				= CSC(1,*) -> LINE INTEGRAL (DW/Dr)
C				= CSC(2,*) -> AREA INTEGRAL D(PHI)s
C				= CSC(3,*) -> INT(DW/Ds)
C				= CSC(4,*) -> AREA INTEGRAL -D(PHI)r
C	----------------------------------------------------------
      DIMENSION SLN(3),RN(3),SN(3),FLR(3),FLS(3)
	DIMENSION DWR(9),DWS(9)
	DIMENSION CSC(4,18)

	BR=-BLR
	BS=-BLS

	CSC(1,4) =DWR(1)/BR
	CSC(1,10)=DWR(2)/BR
	CSC(1,16)=DWR(3)/BR
	CSC(1,5) =DWR(4)/BR
	CSC(1,11)=DWR(5)/BR
	CSC(1,17)=DWR(6)/BR
	CSC(1,3) =DWR(7)/BR
	CSC(1,9) =DWR(8)/BR
	CSC(1,15)=DWR(9)/BR

	CSC(2,4) =DWR(1)
	CSC(2,10)=DWR(2)
	CSC(2,16)=DWR(3)
	CSC(2,5) =DWR(4)+FLS(1)
	CSC(2,11)=DWR(5)+FLS(2)
	CSC(2,17)=DWR(6)+FLS(3)
	CSC(2,3) =DWR(7)
	CSC(2,9) =DWR(8)
	CSC(2,15)=DWR(9)

	CSC(3,4) =DWS(1)/BS
	CSC(3,10)=DWS(2)/BS
	CSC(3,16)=DWS(3)/BS
	CSC(3,5) =DWS(4)/BS
	CSC(3,11)=DWS(5)/BS
	CSC(3,17)=DWS(6)/BS
	CSC(3,3) =DWS(7)/BS
	CSC(3,9) =DWS(8)/BS
	CSC(3,15)=DWS(9)/BS

	CSC(4,4) =DWS(1)+FLR(1)
	CSC(4,10)=DWS(2)+FLR(2)
	CSC(4,16)=DWS(3)+FLR(3)
	CSC(4,5) =DWS(4)
	CSC(4,11)=DWS(5)
	CSC(4,17)=DWS(6)
	CSC(4,3) =DWS(7)
	CSC(4,9) =DWS(8)
	CSC(4,15)=DWS(9)

	RETURN
	END
C=====================================================================
      SUBROUTINE TRICS3(ANT,SLN,RN,SN,FLR,FLS,DWR,DWS,CSC,BLR,BLS)

	IMPLICIT REAL*8 (A-H,O-Z)
C	----------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cs FOR TRANSVERSE SHEAR STRAINS
C
C	INPUT VARIABLES
C	BR,BS		= SHEAR DEFORMATION PARAMETER
C	SLN(3)		= ELEMENT BOUNDARY LENGTHS
C	RN(4),SN(3)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	FLR,FLS
C	DWR,DWS		= COMMON FACTORS DEFINED IN SUBROUTINE COMFACT
C
C	OUTPUT VARIAVBLES
C	CS(2,18)	= MEMBRANE AND DRILLING STRAIN PARAMETER
C	CSC(4,18)	= COMPONENTS OF CS
C				= CSC(1,*) -> LINE INTEGRAL (DW/Dr)
C				= CSC(2,*) -> AREA INTEGRAL D(PHI)s
C				= CSC(3,*) -> INT(DW/Ds)
C				= CSC(4,*) -> AREA INTEGRAL -D(PHI)r
C	----------------------------------------------------------
      DIMENSION SLN(3),RN(3),SN(3),FLR(3),FLS(3)
	DIMENSION DWR(9),DWS(9),ANT(3)
	DIMENSION CSC(4,18)

	BR=-BLR
	BS=-BLS

	CSC(1,4) =DWR(1)/BR*(1+BS)
	CSC(1,10)=DWR(2)/BR*(1+BS)
	CSC(1,16)=DWR(3)/BR*(1+BS)
	CSC(1,5) =DWR(4)/BR*(1+BS)
	CSC(1,11)=DWR(5)/BR*(1+BS)
	CSC(1,17)=DWR(6)/BR*(1+BS)
	CSC(1,3) =DWR(7)/BR*(1+BS)
	CSC(1,9) =DWR(8)/BR*(1+BS)
	CSC(1,15)=DWR(9)/BR*(1+BS)

	CSC(2,4) =0.0
	CSC(2,10)=0.0
	CSC(2,16)=0.0
	CSC(2,5) =ANT(1)*(1+BS)
	CSC(2,11)=ANT(1)*(1+BS)
	CSC(2,17)=ANT(1)*(1+BS)
	CSC(2,3) =0.0
	CSC(2,9) =0.0
	CSC(2,15)=0.0

	CSC(3,4) =DWS(1)/BS*(1+BR)
	CSC(3,10)=DWS(2)/BS*(1+BR)
	CSC(3,16)=DWS(3)/BS*(1+BR)
	CSC(3,5) =DWS(4)/BS*(1+BR)
	CSC(3,11)=DWS(5)/BS*(1+BR)
	CSC(3,17)=DWS(6)/BS*(1+BR)
	CSC(3,3) =DWS(7)/BS*(1+BR)
	CSC(3,9) =DWS(8)/BS*(1+BR)
	CSC(3,15)=DWS(9)/BS*(1+BR)

	CSC(4,4) =-ANT(1)*(1+BR)
	CSC(4,10)=-ANT(1)*(1+BR)
	CSC(4,16)=-ANT(1)*(1+BR)
	CSC(4,5) =0.0
	CSC(4,11)=0.0
	CSC(4,17)=0.0
	CSC(4,3) =0.0
	CSC(4,9) =0.0
	CSC(4,15)=0.0

	RETURN
	END
C=====================================================================
      SUBROUTINE TABCS(AREA,CSC,ACS,ACSC)

	IMPLICIT REAL*8 (A-H,O-Z)
C	----------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR INVERSE(Ab)*(Cb)*TRANSPOSE(T)=ACB
C
C	INPUT VARIABLES
C	AREA				= AREA OF TRIANGULAR ELEMENT
C	CS(2,18)			= TRANSVERSE SHEAR STRAIN PARAMETER
C	CSC(4,18)			= COMPONENTS OF CS
C
C	LOCAL VARIALBLE
C	ASI(2,2)			= INVERSE OF As
C
C	OUTPUT VARIAVBLES
C	ACS(2,18)			= INVERSE(As)*(Cs)*TRANSPOSE(T)
C	ACSC(4,18)			= COMPONENTS OF ACSC
C	----------------------------------------------------------
	DIMENSION CS(2,18),CSC(4,18)
	DIMENSION ASI(2,2)
	DIMENSION ACS(2,18),ACSC(4,18)

C	---------------------------
C	SOLVE FOR INVERSE(As) - ASI
C	---------------------------
	ASI(1,1)=1.0/AREA
	ASI(2,2)=ASI(1,1)

C	-------------------------------
C	SOLVE FOR ACS COMPONENTS - ACSC
C	-------------------------------
      ACSC(1,3) = ASI(1,1)*CSC(1,3)
      ACSC(1,4) = ASI(1,1)*CSC(1,4)
      ACSC(1,5) = ASI(1,1)*CSC(1,5)
      ACSC(1,9) = ASI(1,1)*CSC(1,9)
      ACSC(1,10) = ASI(1,1)*CSC(1,10)
      ACSC(1,11) = ASI(1,1)*CSC(1,11)
      ACSC(1,15) = ASI(1,1)*CSC(1,15)
      ACSC(1,16) = ASI(1,1)*CSC(1,16)
      ACSC(1,17) = ASI(1,1)*CSC(1,17)

      ACSC(2,3) = ASI(1,1)*CSC(2,3)
      ACSC(2,4) = ASI(1,1)*CSC(2,4)
      ACSC(2,5) = ASI(1,1)*CSC(2,5)
      ACSC(2,9) = ASI(1,1)*CSC(2,9)
      ACSC(2,10) = ASI(1,1)*CSC(2,10)
      ACSC(2,11) = ASI(1,1)*CSC(2,11)
      ACSC(2,15) = ASI(1,1)*CSC(2,15)
      ACSC(2,16) = ASI(1,1)*CSC(2,16)
      ACSC(2,17) = ASI(1,1)*CSC(2,17)

      ACSC(3,3) = ASI(2,2)*CSC(3,3)
      ACSC(3,4) = ASI(2,2)*CSC(3,4)
      ACSC(3,5) = ASI(2,2)*CSC(3,5)
      ACSC(3,9) = ASI(2,2)*CSC(3,9)
      ACSC(3,10) = ASI(2,2)*CSC(3,10)
      ACSC(3,11) = ASI(2,2)*CSC(3,11)
      ACSC(3,15) = ASI(2,2)*CSC(3,15)
      ACSC(3,16) = ASI(2,2)*CSC(3,16)
      ACSC(3,17) = ASI(2,2)*CSC(3,17)

      ACSC(4,3) = ASI(2,2)*CSC(4,3)
      ACSC(4,4) = ASI(2,2)*CSC(4,4)
      ACSC(4,5) = ASI(2,2)*CSC(4,5)
      ACSC(4,9) = ASI(2,2)*CSC(4,9)
      ACSC(4,10) = ASI(2,2)*CSC(4,10)
      ACSC(4,11) = ASI(2,2)*CSC(4,11)
      ACSC(4,15) = ASI(2,2)*CSC(4,15)
      ACSC(4,16) = ASI(2,2)*CSC(4,16)
      ACSC(4,17) = ASI(2,2)*CSC(4,17)

C	--------------------------------
C	SOLVE FOR INVERSE(As)*(Cs) - ACS
C	--------------------------------
      ACS(1,3) = ACSC(1,3)+ACSC(2,3)
      ACS(1,4) = ACSC(1,4)+ACSC(2,4)
      ACS(1,5) = ACSC(1,5)+ACSC(2,5)
      ACS(1,9) = ACSC(1,9)+ACSC(2,9)
      ACS(1,10) = ACSC(1,10)+ACSC(2,10)
      ACS(1,11) = ACSC(1,11)+ACSC(2,11)
      ACS(1,15) = ACSC(1,15)+ACSC(2,15)
      ACS(1,16) = ACSC(1,16)+ACSC(2,16)
      ACS(1,17) = ACSC(1,17)+ACSC(2,17)

      ACS(2,3) = ACSC(3,3)+ACSC(4,3)
      ACS(2,4) = ACSC(3,4)+ACSC(4,4)
      ACS(2,5) = ACSC(3,5)+ACSC(4,5)
      ACS(2,9) = ACSC(3,9)+ACSC(4,9)
      ACS(2,10) = ACSC(3,10)+ACSC(4,10)
      ACS(2,11) = ACSC(3,11)+ACSC(4,11)
      ACS(2,15) = ACSC(3,15)+ACSC(4,15)
      ACS(2,16) = ACSC(3,16)+ACSC(4,16)
      ACS(2,17) = ACSC(3,17)+ACSC(4,17)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE APHAS18(ACSC,TEDIS,APS)

	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR PARAMETER ALPHA
C
C	INPUT VARIABLES
C	ACSC(4,24)	= COMPONENTS OF ACS
C	TEDIS(48)	= TRANSPOSE(T)*REDIS
C
C	OUTPUT VARIAVBLES
C	APS(4)		= MEMBRANE STRAIN PARAMETER
C				-> APS(1) - FOR (DW/Dr)
C				-> APS(2) - FOR D(PHI)s
C				-> APS(3) - FOR DW/Ds
C				-> APS(4) - FOR -D(PHI)r
C	--------------------------------------------------------
      DIMENSION ACSC(4,18),TEDIS(18)
	DIMENSION APS(4)

      APS(1) = ACSC(1,3)*TEDIS(3)+ACSC(1,4)*TEDIS(4)+ACSC(1,5)*TEDIS(5)+
     #ACSC(1,9)*TEDIS(9)+ACSC(1,10)*TEDIS(10)+ACSC(1,11)*TEDIS(11)+ACSC(
     #1,15)*TEDIS(15)+ACSC(1,16)*TEDIS(16)+ACSC(1,17)*TEDIS(17)

      APS(2) = ACSC(2,3)*TEDIS(3)+ACSC(2,4)*TEDIS(4)+ACSC(2,5)*TEDIS(5)+
     #ACSC(2,9)*TEDIS(9)+ACSC(2,10)*TEDIS(10)+ACSC(2,11)*TEDIS(11)+ACSC(
     #2,15)*TEDIS(15)+ACSC(2,16)*TEDIS(16)+ACSC(2,17)*TEDIS(17)

      APS(3) = ACSC(3,3)*TEDIS(3)+ACSC(3,4)*TEDIS(4)+ACSC(3,5)*TEDIS(5)+
     #ACSC(3,9)*TEDIS(9)+ACSC(3,10)*TEDIS(10)+ACSC(3,11)*TEDIS(11)+ACSC(
     #3,15)*TEDIS(15)+ACSC(3,16)*TEDIS(16)+ACSC(3,17)*TEDIS(17)

      APS(4) = ACSC(4,3)*TEDIS(3)+ACSC(4,4)*TEDIS(4)+ACSC(4,5)*TEDIS(5)+
     #ACSC(4,9)*TEDIS(9)+ACSC(4,10)*TEDIS(10)+ACSC(4,11)*TEDIS(11)+ACSC(
     #4,15)*TEDIS(15)+ACSC(4,16)*TEDIS(16)+ACSC(4,17)*TEDIS(17)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE KSHR18(ACS,PTP,ST)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-----------------------------------------------------------
C	PURPOSE:	TO ADD THE TRANSVERSE SHEAR CONTRIBUTION TO THE
C					ELEMENT STIFFNESS MATRIX
C
C	INPUT VARIABLES
C	ACS(2,18)	= INVERSE(As)*(Cs)
C	PTP(2,2)	= Ns*T*Ps -> PTP(1,1)=PTP(2,2)=AB1(1,1)*DR(55)
C				= A DIAGONAL SYMMETRIC MATRIX
C
C	OUTPUT VARIAVBLES
C	ST(171)		= STIFFNESS MATRIX IN ROW FORMAT
C	-----------------------------------------------------------
	DIMENSION PTP(2,2),ACS(2,18)
	DIMENSION ST(171)

      t1 = ACS(1,3)**2
      t3 = ACS(2,3)**2
      t6 = ACS(1,3)*PTP(1,1)
      t8 = ACS(2,3)*PTP(2,2)
      t32 = ACS(1,4)**2
      t34 = ACS(2,4)**2
      t37 = ACS(1,4)*PTP(1,1)
      t39 = ACS(2,4)*PTP(2,2)
      t60 = ACS(1,5)**2
      t62 = ACS(2,5)**2
      t65 = ACS(1,5)*PTP(1,1)
      t67 = ACS(2,5)*PTP(2,2)
      t85 = ACS(1,9)**2
      t87 = ACS(2,9)**2
      t90 = ACS(1,9)*PTP(1,1)
      t92 = ACS(2,9)*PTP(2,2)
      t107 = ACS(1,10)**2
      t109 = ACS(2,10)**2
      t112 = ACS(1,10)*PTP(1,1)
      t114 = ACS(2,10)*PTP(2,2)
      t126 = ACS(1,11)**2
      t128 = ACS(2,11)**2
      t131 = ACS(1,11)*PTP(1,1)
      t133 = ACS(2,11)*PTP(2,2)
      t142 = ACS(1,15)**2
      t144 = ACS(2,15)**2
      t147 = ACS(1,15)*PTP(1,1)
      t149 = ACS(2,15)*PTP(2,2)
      t155 = ACS(1,16)**2
      t157 = ACS(2,16)**2
      t165 = ACS(1,17)**2
      t167 = ACS(2,17)**2
      ST(1) = ST(1)
      ST(2) = ST(2)
      ST(3) = ST(3)
      ST(4) = ST(4)
      ST(5) = ST(5)
      ST(6) = ST(6)
      ST(7) = ST(7)
      ST(8) = ST(8)
      ST(9) = ST(9)
      ST(10) = ST(10)
      ST(11) = ST(11)
      ST(12) = ST(12)
      ST(13) = ST(13)
      ST(14) = ST(14)
      ST(15) = ST(15)
      ST(16) = ST(16)
      ST(17) = ST(17)
      ST(18) = ST(18)
      ST(19) = ST(19)
      ST(20) = ST(20)
      ST(21) = ST(21)
      ST(22) = ST(22)
      ST(23) = ST(23)
      ST(24) = ST(24)
      ST(25) = ST(25)
      ST(26) = ST(26)
      ST(27) = ST(27)
      ST(28) = ST(28)
      ST(29) = ST(29)
      ST(30) = ST(30)
      ST(31) = ST(31)
      ST(32) = ST(32)
      ST(33) = ST(33)
      ST(34) = ST(34)
      ST(35) = ST(35)
      ST(36) = ST(36)+t1*PTP(1,1)+t3*PTP(2,2)
      ST(37) = ST(37)+t6*ACS(1,4)+t8*ACS(2,4)
      ST(38) = ST(38)+t6*ACS(1,5)+t8*ACS(2,5)
      ST(39) = ST(39)
      ST(40) = ST(40)
      ST(41) = ST(41)
      ST(42) = ST(42)+t6*ACS(1,9)+t8*ACS(2,9)
      ST(43) = ST(43)+t6*ACS(1,10)+t8*ACS(2,10)
      ST(44) = ST(44)+t6*ACS(1,11)+t8*ACS(2,11)
      ST(45) = ST(45)
      ST(46) = ST(46)
      ST(47) = ST(47)
      ST(48) = ST(48)+t6*ACS(1,15)+t8*ACS(2,15)
      ST(49) = ST(49)+t6*ACS(1,16)+t8*ACS(2,16)
      ST(50) = ST(50)+t6*ACS(1,17)+t8*ACS(2,17)
      ST(51) = ST(51)
      ST(52) = ST(52)+t32*PTP(1,1)+t34*PTP(2,2)
      ST(53) = ST(53)+t37*ACS(1,5)+t39*ACS(2,5)
      ST(54) = ST(54)
      ST(55) = ST(55)
      ST(56) = ST(56)
      ST(57) = ST(57)+t37*ACS(1,9)+t39*ACS(2,9)
      ST(58) = ST(58)+t37*ACS(1,10)+t39*ACS(2,10)
      ST(59) = ST(59)+t37*ACS(1,11)+t39*ACS(2,11)
      ST(60) = ST(60)
      ST(61) = ST(61)
      ST(62) = ST(62)
      ST(63) = ST(63)+t37*ACS(1,15)+t39*ACS(2,15)
      ST(64) = ST(64)+t37*ACS(1,16)+t39*ACS(2,16)
      ST(65) = ST(65)+t37*ACS(1,17)+t39*ACS(2,17)
      ST(66) = ST(66)
      ST(67) = ST(67)+t60*PTP(1,1)+t62*PTP(2,2)
      ST(68) = ST(68)
      ST(69) = ST(69)
      ST(70) = ST(70)
      ST(71) = ST(71)+t65*ACS(1,9)+t67*ACS(2,9)
      ST(72) = ST(72)+t65*ACS(1,10)+t67*ACS(2,10)
      ST(73) = ST(73)+t65*ACS(1,11)+t67*ACS(2,11)
      ST(74) = ST(74)
      ST(75) = ST(75)
      ST(76) = ST(76)
      ST(77) = ST(77)+t65*ACS(1,15)+t67*ACS(2,15)
      ST(78) = ST(78)+t65*ACS(1,16)+t67*ACS(2,16)
      ST(79) = ST(79)+t65*ACS(1,17)+t67*ACS(2,17)
      ST(80) = ST(80)
      ST(81) = ST(81)
      ST(82) = ST(82)
      ST(83) = ST(83)
      ST(84) = ST(84)
      ST(85) = ST(85)
      ST(86) = ST(86)
      ST(87) = ST(87)
      ST(88) = ST(88)
      ST(89) = ST(89)
      ST(90) = ST(90)
      ST(91) = ST(91)
      ST(92) = ST(92)
      ST(93) = ST(93)
      ST(94) = ST(94)
      ST(95) = ST(95)
      ST(96) = ST(96)
      ST(97) = ST(97)
      ST(98) = ST(98)
      ST(99) = ST(99)
      ST(100) = ST(100)
      ST(101) = ST(101)
      ST(102) = ST(102)
      ST(103) = ST(103)
      ST(104) = ST(104)
      ST(105) = ST(105)
      ST(106) = ST(106)
      ST(107) = ST(107)
      ST(108) = ST(108)
      ST(109) = ST(109)
      ST(110) = ST(110)
      ST(111) = ST(111)
      ST(112) = ST(112)
      ST(113) = ST(113)
      ST(114) = ST(114)
      ST(115) = ST(115)
      ST(116) = ST(116)
      ST(117) = ST(117)+t85*PTP(1,1)+t87*PTP(2,2)
      ST(118) = ST(118)+t90*ACS(1,10)+t92*ACS(2,10)
      ST(119) = ST(119)+t90*ACS(1,11)+t92*ACS(2,11)
      ST(120) = ST(120)
      ST(121) = ST(121)
      ST(122) = ST(122)
      ST(123) = ST(123)+t90*ACS(1,15)+t92*ACS(2,15)
      ST(124) = ST(124)+t90*ACS(1,16)+t92*ACS(2,16)
      ST(125) = ST(125)+t90*ACS(1,17)+t92*ACS(2,17)
      ST(126) = ST(126)
      ST(127) = ST(127)+t107*PTP(1,1)+t109*PTP(2,2)
      ST(128) = ST(128)+t112*ACS(1,11)+t114*ACS(2,11)
      ST(129) = ST(129)
      ST(130) = ST(130)
      ST(131) = ST(131)
      ST(132) = ST(132)+t112*ACS(1,15)+t114*ACS(2,15)
      ST(133) = ST(133)+t112*ACS(1,16)+t114*ACS(2,16)
      ST(134) = ST(134)+t112*ACS(1,17)+t114*ACS(2,17)
      ST(135) = ST(135)
      ST(136) = ST(136)+t126*PTP(1,1)+t128*PTP(2,2)
      ST(137) = ST(137)
      ST(138) = ST(138)
      ST(139) = ST(139)
      ST(140) = ST(140)+t131*ACS(1,15)+t133*ACS(2,15)
      ST(141) = ST(141)+t131*ACS(1,16)+t133*ACS(2,16)
      ST(142) = ST(142)+t131*ACS(1,17)+t133*ACS(2,17)
      ST(143) = ST(143)
      ST(144) = ST(144)
      ST(145) = ST(145)
      ST(146) = ST(146)
      ST(147) = ST(147)
      ST(148) = ST(148)
      ST(149) = ST(149)
      ST(150) = ST(150)
      ST(151) = ST(151)
      ST(152) = ST(152)
      ST(153) = ST(153)
      ST(154) = ST(154)
      ST(155) = ST(155)
      ST(156) = ST(156)
      ST(157) = ST(157)
      ST(158) = ST(158)
      ST(159) = ST(159)
      ST(160) = ST(160)
      ST(161) = ST(161)
      ST(162) = ST(162)+t142*PTP(1,1)+t144*PTP(2,2)
      ST(163) = ST(163)+t147*ACS(1,16)+t149*ACS(2,16)
      ST(164) = ST(164)+t147*ACS(1,17)+t149*ACS(2,17)
      ST(165) = ST(165)
      ST(166) = ST(166)+t155*PTP(1,1)+t157*PTP(2,2)
      ST(167) = ST(167)+ACS(1,16)*PTP(1,1)*ACS(1,17)+ACS(2,16)*
     #PTP(2,2)*ACS(2,17)
      ST(168) = ST(168)
      ST(169) = ST(169)+t165*PTP(1,1)+t167*PTP(2,2)
      ST(170) = ST(170)
      ST(171) = ST(171)

	RETURN
	END
C
C=====================================================================





C======================================================================
C	START 3 NODE SHELL ELEMENT - DRILLING SUBROUTINES
C======================================================================
C=====================================================================
      SUBROUTINE TRICD(ANT,CM,CD)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cd FOR DRILLING STRAINS
C
C	INPUT VARIABLES
C	ANT(4)		= AREA INTEGRATION FACTORS
C	CM(5,24)	= MEMBRANE STIFFNESS PARAMENTERS HAVING COMMON TERMS
C					WITH THE DRILLING STRAIN
C
C	OUTPUT VARIABLES
C	CD(24)		= DRILLING STRAIN PARAMETER
C	--------------------------------------------------------------
	DIMENSION ANT(3),CM(3,18)
	DIMENSION CD(18)

C	---------------------------------
C	SOLVE FOR Cd FOR DRILLING STRAINS
C	---------------------------------
C	ELEMENTS FOR DISPLACEMENT U
	CD(1) = 0.5*CM(3,1)
	CD(7) = 0.5*CM(3,7)
	CD(13) = 0.5*CM(3,13)
C	ELEMENTS FOR DISPLACEMENT V
	CD(2) = -0.5*CM(3,2)
	CD(8) = -0.5*CM(3,8)
	CD(14) = -0.5*CM(3,14)
C	ELEMENTS FOR DISPLACEMENT THETA-t
	CD(6) = ANT(1)
	CD(12) = ANT(2)
	CD(18) = ANT(3)

      RETURN
      END
C
C=====================================================================
      SUBROUTINE KDRL18(ACD,DFAC,ST)

	IMPLICIT REAL*8 (A-H,O-Z)
C	-----------------------------------------------------------
C	PURPOSE:	TO ADD THE DRILLING CONTRIBUTION TO THE
C					ELEMENT STIFFNESS MATRIX
C
C	INPUT VARIABLES
C	ACD(18)	= INVERSE(Ad)*(Cd)
C	DFAC	= DRILLING RIGIDITY
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIAVBLES
C	ST(171)		= STIFFNESS MATRIX IN ROW FORMAT
C	-----------------------------------------------------------
      DIMENSION ACD(18)
	DIMENSION ST(171)

      t1 = ACD(1)**2
      t4 = DFAC*ACD(1)
      t21 = ACD(2)**2
      t24 = DFAC*ACD(2)
      t39 = ACD(6)**2
      t42 = DFAC*ACD(6)
      t55 = ACD(7)**2
      t58 = DFAC*ACD(7)
      t69 = ACD(8)**2
      t72 = DFAC*ACD(8)
      t81 = ACD(12)**2
      t84 = DFAC*ACD(12)
      t91 = ACD(13)**2
      t94 = DFAC*ACD(13)
      t99 = ACD(14)**2
      t105 = ACD(18)**2

      ST(1) = ST(1)+DFAC*t1
      ST(2) = ST(2)+t4*ACD(2)

      ST(6) = ST(6)+t4*ACD(6)
      ST(7) = ST(7)+t4*ACD(7)
      ST(8) = ST(8)+t4*ACD(8)

      ST(12) = ST(12)+t4*ACD(12)
      ST(13) = ST(13)+t4*ACD(13)
      ST(14) = ST(14)+t4*ACD(14)

      ST(18) = ST(18)+t4*ACD(18)
      ST(19) = ST(19)+DFAC*t21

      ST(23) = ST(23)+t24*ACD(6)
      ST(24) = ST(24)+t24*ACD(7)
      ST(25) = ST(25)+t24*ACD(8)

      ST(29) = ST(29)+t24*ACD(12)
      ST(30) = ST(30)+t24*ACD(13)
      ST(31) = ST(31)+t24*ACD(14)

      ST(35) = ST(35)+t24*ACD(18)

      ST(81) = ST(81)+DFAC*t39
      ST(82) = ST(82)+t42*ACD(7)
      ST(83) = ST(83)+t42*ACD(8)

      ST(87) = ST(87)+t42*ACD(12)
      ST(88) = ST(88)+t42*ACD(13)
      ST(89) = ST(89)+t42*ACD(14)

      ST(93) = ST(93)+t42*ACD(18)
      ST(94) = ST(94)+DFAC*t55
      ST(95) = ST(95)+t58*ACD(8)

      ST(99) = ST(99)+t58*ACD(12)
      ST(100) = ST(100)+t58*ACD(13)
      ST(101) = ST(101)+t58*ACD(14)

      ST(105) = ST(105)+t58*ACD(18)
      ST(106) = ST(106)+DFAC*t69

      ST(110) = ST(110)+t72*ACD(12)
      ST(111) = ST(111)+t72*ACD(13)
      ST(112) = ST(112)+t72*ACD(14)

      ST(116) = ST(116)+t72*ACD(18)

      ST(144) = ST(144)+DFAC*t81
      ST(145) = ST(145)+t84*ACD(13)
      ST(146) = ST(146)+t84*ACD(14)

      ST(150) = ST(150)+t84*ACD(18)
      ST(151) = ST(151)+DFAC*t91
      ST(152) = ST(152)+t94*ACD(14)

      ST(156) = ST(156)+t94*ACD(18)
      ST(157) = ST(157)+DFAC*t99

      ST(161) = ST(161)+DFAC*ACD(14)*ACD(18)

      ST(171) = ST(171)+DFAC*t105

      RETURN
      END
C
C=====================================================================
      SUBROUTINE PDRPT(DR,R,S,PMPM,PMPB,PBPM,PBPB,PTP,MTMOD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTES FOR INTEGRAL(TRANS(P)*DR*P)dA
C
C	INPUT VARIABLES
C	DR(64)  = NODAL STRESS RESULTANTS
C	A       = FACTOR FOR INT(DR)dA
C     AR      = FACTOR FOR INT(R*DR)dA
C     AS      = FACTOR FOR INT(S*DR)dA
C     ARS     = FACTOR FOR INT(R*S*DR)dA
C     ARR     = FACTOR FOR INT(R*R*DR)dA
C     ASS     = FACTOR FOR INT(S*S*DR)dA
C     ARRS    = FACTOR FOR INT(R*R*S*DR)dA
C     ASSR    = FACTOR FOR INT(S*S*R*DR)dA
C     ARRSS   = FACTOR FOR INT(R*R*S*S*DR)dA
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C	PMPM(5,5)   = INTEGRAL(TRANS(Pm)*DR*Pm) - MEMBRANE PART
C	PMPB(5,11)  = INTEGRAL(TRANS(Pm)*DR*Pb) - MEMBRANE-BENDING PART
C	PBPB(11,11) = INTEGRAL(TRANS(Pb)*DR*Pb) - BENDING PART
C	PTP(2,2)	= INTEGRAL(TRANS(Ps)*DR*Ps) - TRANSVERSE SHEAR PART
C     ----------------------------------------------------------------
	DIMENSION DR(64)
	DIMENSION PMPM(7,7),PMPB(7,9),PBPM(9,7),PBPB(9,9), PTP(2,2)

C	----------------------------
C	COMPUTE PMPM - MEMBRANE PART
C	----------------------------
	PMPM(1,1) = PMPM(1,1)+DR(1)
      PMPM(1,2) = PMPM(1,2)+R*DR(1)-S*DR(3)
      PMPM(1,3) = PMPM(1,3)+S*DR(1)
      PMPM(1,4) = PMPM(1,4)+DR(2)
      PMPM(1,5) = PMPM(1,5)+R*DR(2)
      PMPM(1,6) = PMPM(1,6)+S*DR(2)-R*DR(3)
      PMPM(1,7) = PMPM(1,7)+DR(3)
      PMPM(2,1) = PMPM(2,1)+R*DR(1)-S*DR(17)
      PMPM(2,2) = PMPM(2,2)+R**2*DR(1)-R*S*DR(17)-S*DR(3)*R+S**2*DR(19)
      PMPM(2,3) = PMPM(2,3)+S*DR(1)*R-S**2*DR(17)
      PMPM(2,4) = PMPM(2,4)+R*DR(2)-S*DR(18)
      PMPM(2,5) = PMPM(2,5)+R**2*DR(2)-R*S*DR(18)
      PMPM(2,6) = PMPM(2,6)+S*DR(2)*R-S**2*DR(18)-R**2*DR(3)+R*S*DR(19)
      PMPM(2,7) = PMPM(2,7)+R*DR(3)-S*DR(19)
      PMPM(3,1) = PMPM(3,1)+S*DR(1)
      PMPM(3,2) = PMPM(3,2)+S*DR(1)*R-S**2*DR(3)
      PMPM(3,3) = PMPM(3,3)+S**2*DR(1)
      PMPM(3,4) = PMPM(3,4)+S*DR(2)
      PMPM(3,5) = PMPM(3,5)+S*DR(2)*R
      PMPM(3,6) = PMPM(3,6)+S**2*DR(2)-S*DR(3)*R
      PMPM(3,7) = PMPM(3,7)+S*DR(3)
      PMPM(4,1) = PMPM(4,1)+DR(9)
      PMPM(4,2) = PMPM(4,2)+R*DR(9)-S*DR(11)
      PMPM(4,3) = PMPM(4,3)+S*DR(9)
      PMPM(4,4) = PMPM(4,4)+DR(10)
      PMPM(4,5) = PMPM(4,5)+R*DR(10)
      PMPM(4,6) = PMPM(4,6)+S*DR(10)-R*DR(11)
      PMPM(4,7) = PMPM(4,7)+DR(11)
      PMPM(5,1) = PMPM(5,1)+R*DR(9)
      PMPM(5,2) = PMPM(5,2)+R**2*DR(9)-R*DR(11)*S
      PMPM(5,3) = PMPM(5,3)+R*DR(9)*S
      PMPM(5,4) = PMPM(5,4)+R*DR(10)
      PMPM(5,5) = PMPM(5,5)+R**2*DR(10)
      PMPM(5,6) = PMPM(5,6)+R*DR(10)*S-R**2*DR(11)
      PMPM(5,7) = PMPM(5,7)+R*DR(11)
      PMPM(6,1) = PMPM(6,1)+S*DR(9)-R*DR(17)
      PMPM(6,2) = PMPM(6,2)+R*DR(9)*S-R**2*DR(17)-S**2*DR(11)+R*S*DR(19)
      PMPM(6,3) = PMPM(6,3)+S**2*DR(9)-R*S*DR(17)
      PMPM(6,4) = PMPM(6,4)+S*DR(10)-R*DR(18)
      PMPM(6,5) = PMPM(6,5)+R*DR(10)*S-R**2*DR(18)
      PMPM(6,6) = PMPM(6,6)+S**2*DR(10)-R*S*DR(18)-R*DR(11)*S+
     #                      R**2*DR(19)
      PMPM(6,7) = PMPM(6,7)+S*DR(11)-R*DR(19)
      PMPM(7,1) = PMPM(7,1)+DR(17)
      PMPM(7,2) = PMPM(7,2)+R*DR(17)-S*DR(19)
      PMPM(7,3) = PMPM(7,3)+S*DR(17)
      PMPM(7,4) = PMPM(7,4)+DR(18)
      PMPM(7,5) = PMPM(7,5)+R*DR(18)
      PMPM(7,6) = PMPM(7,6)+S*DR(18)-R*DR(19)
      PMPM(7,7) = PMPM(7,7)+DR(19)

C	---------------------------------------------
C	COMPUTE PMPB - MEMBRANE-BENDING COUPLING PART
C	---------------------------------------------
	IF (MTMOD.GT.1) THEN
			PMPB(1,1) = PMPB(1,1)+DR(4)
      		PMPB(1,2) = PMPB(1,2)+R*DR(4)
      		PMPB(1,3) = PMPB(1,3)+S*DR(4)
      		PMPB(1,4) = PMPB(1,4)+DR(5)
      		PMPB(1,5) = PMPB(1,5)+R*DR(5)
      		PMPB(1,6) = PMPB(1,6)+S*DR(5)
      		PMPB(1,7) = PMPB(1,7)+DR(6)
      		PMPB(1,8) = PMPB(1,8)+R*DR(6)
      		PMPB(1,9) = PMPB(1,9)+S*DR(6)
      		PMPB(2,1) = PMPB(2,1)+R*DR(4)-S*DR(20)
      		PMPB(2,2) = PMPB(2,2)+R**2*DR(4)-R*S*DR(20)
      		PMPB(2,3) = PMPB(2,3)+S*DR(4)*R-S**2*DR(20)
      		PMPB(2,4) = PMPB(2,4)+R*DR(5)-S*DR(21)
      		PMPB(2,5) = PMPB(2,5)+R**2*DR(5)-R*S*DR(21)
      		PMPB(2,6) = PMPB(2,6)+S*DR(5)*R-S**2*DR(21)
      		PMPB(2,7) = PMPB(2,7)+R*DR(6)-S*DR(22)
      		PMPB(2,8) = PMPB(2,8)+R**2*DR(6)-R*S*DR(22)
      		PMPB(2,9) = PMPB(2,9)+S*DR(6)*R-S**2*DR(22)
      		PMPB(3,1) = PMPB(3,1)+S*DR(4)
      		PMPB(3,2) = PMPB(3,2)+S*DR(4)*R
      		PMPB(3,3) = PMPB(3,3)+S**2*DR(4)
      		PMPB(3,4) = PMPB(3,4)+S*DR(5)
      		PMPB(3,5) = PMPB(3,5)+S*DR(5)*R
      		PMPB(3,6) = PMPB(3,6)+S**2*DR(5)
      		PMPB(3,7) = PMPB(3,7)+S*DR(6)
      		PMPB(3,8) = PMPB(3,8)+S*DR(6)*R
      		PMPB(3,9) = PMPB(3,9)+S**2*DR(6)
      		PMPB(4,1) = PMPB(4,1)+DR(12)
      		PMPB(4,2) = PMPB(4,2)+R*DR(12)
      		PMPB(4,3) = PMPB(4,3)+S*DR(12)
      		PMPB(4,4) = PMPB(4,4)+DR(13)
      		PMPB(4,5) = PMPB(4,5)+R*DR(13)
      		PMPB(4,6) = PMPB(4,6)+S*DR(13)
      		PMPB(4,7) = PMPB(4,7)+DR(14)
      		PMPB(4,8) = PMPB(4,8)+R*DR(14)
      		PMPB(4,9) = PMPB(4,9)+S*DR(14)
      		PMPB(5,1) = PMPB(5,1)+R*DR(12)
      		PMPB(5,2) = PMPB(5,2)+R**2*DR(12)
      		PMPB(5,3) = PMPB(5,3)+R*DR(12)*S
      		PMPB(5,4) = PMPB(5,4)+R*DR(13)
      		PMPB(5,5) = PMPB(5,5)+R**2*DR(13)
      		PMPB(5,6) = PMPB(5,6)+R*DR(13)*S
      		PMPB(5,7) = PMPB(5,7)+R*DR(14)
      		PMPB(5,8) = PMPB(5,8)+R**2*DR(14)
      		PMPB(5,9) = PMPB(5,9)+R*DR(14)*S
      		PMPB(6,1) = PMPB(6,1)+S*DR(12)-R*DR(20)
      		PMPB(6,2) = PMPB(6,2)+R*DR(12)*S-R**2*DR(20)
      		PMPB(6,3) = PMPB(6,3)+S**2*DR(12)-R*S*DR(20)
      		PMPB(6,4) = PMPB(6,4)+S*DR(13)-R*DR(21)
      		PMPB(6,5) = PMPB(6,5)+R*DR(13)*S-R**2*DR(21)
      		PMPB(6,6) = PMPB(6,6)+S**2*DR(13)-R*S*DR(21)
      		PMPB(6,7) = PMPB(6,7)+S*DR(14)-R*DR(22)
      		PMPB(6,8) = PMPB(6,8)+R*DR(14)*S-R**2*DR(22)
      		PMPB(6,9) = PMPB(6,9)+S**2*DR(14)-R*S*DR(22)
      		PMPB(7,1) = PMPB(7,1)+DR(20)
      		PMPB(7,2) = PMPB(7,2)+R*DR(20)
      		PMPB(7,3) = PMPB(7,3)+S*DR(20)
      		PMPB(7,4) = PMPB(7,4)+DR(21)
      		PMPB(7,5) = PMPB(7,5)+R*DR(21)
      		PMPB(7,6) = PMPB(7,6)+S*DR(21)
      		PMPB(7,7) = PMPB(7,7)+DR(22)
      		PMPB(7,8) = PMPB(7,8)+R*DR(22)
      		PMPB(7,9) = PMPB(7,9)+S*DR(22)
	ENDIF

C	---------------------------------------------
C	COMPUTE PBPM - MEMBRANE-BENDING COUPLING PART
C	---------------------------------------------
	IF (MTMOD.GT.1) THEN
			PBPM(1,1) = PBPM(1,1)+DR(25)
      		PBPM(1,2) = PBPM(1,2)+R*DR(25)-S*DR(27)
      		PBPM(1,3) = PBPM(1,3)+S*DR(25)
      		PBPM(1,4) = PBPM(1,4)+DR(26)
      		PBPM(1,5) = PBPM(1,5)+R*DR(26)
      		PBPM(1,6) = PBPM(1,6)+S*DR(26)-R*DR(27)
      		PBPM(1,7) = PBPM(1,7)+DR(27)
      		PBPM(2,1) = PBPM(2,1)+R*DR(25)
      		PBPM(2,2) = PBPM(2,2)+R**2*DR(25)-R*DR(27)*S
      		PBPM(2,3) = PBPM(2,3)+R*DR(25)*S
      		PBPM(2,4) = PBPM(2,4)+R*DR(26)
      		PBPM(2,5) = PBPM(2,5)+R**2*DR(26)
      		PBPM(2,6) = PBPM(2,6)+R*DR(26)*S-R**2*DR(27)
      		PBPM(2,7) = PBPM(2,7)+R*DR(27)
      		PBPM(3,1) = PBPM(3,1)+S*DR(25)
      		PBPM(3,2) = PBPM(3,2)+R*DR(25)*S-S**2*DR(27)
      		PBPM(3,3) = PBPM(3,3)+S**2*DR(25)
      		PBPM(3,4) = PBPM(3,4)+S*DR(26)
      		PBPM(3,5) = PBPM(3,5)+R*DR(26)*S
      		PBPM(3,6) = PBPM(3,6)+S**2*DR(26)-R*DR(27)*S
      		PBPM(3,7) = PBPM(3,7)+S*DR(27)
      		PBPM(4,1) = PBPM(4,1)+DR(33)
      		PBPM(4,2) = PBPM(4,2)+R*DR(33)-S*DR(35)
      		PBPM(4,3) = PBPM(4,3)+S*DR(33)
      		PBPM(4,4) = PBPM(4,4)+DR(34)
      		PBPM(4,5) = PBPM(4,5)+R*DR(34)
      		PBPM(4,6) = PBPM(4,6)+S*DR(34)-R*DR(35)
      		PBPM(4,7) = PBPM(4,7)+DR(35)
      		PBPM(5,1) = PBPM(5,1)+R*DR(33)
      		PBPM(5,2) = PBPM(5,2)+R**2*DR(33)-R*DR(35)*S
      		PBPM(5,3) = PBPM(5,3)+R*DR(33)*S
      		PBPM(5,4) = PBPM(5,4)+R*DR(34)
      		PBPM(5,5) = PBPM(5,5)+R**2*DR(34)
      		PBPM(5,6) = PBPM(5,6)+R*DR(34)*S-R**2*DR(35)
      		PBPM(5,7) = PBPM(5,7)+R*DR(35)
      		PBPM(6,1) = PBPM(6,1)+S*DR(33)
      		PBPM(6,2) = PBPM(6,2)+R*DR(33)*S-S**2*DR(35)
      		PBPM(6,3) = PBPM(6,3)+S**2*DR(33)
      		PBPM(6,4) = PBPM(6,4)+S*DR(34)
      		PBPM(6,5) = PBPM(6,5)+R*DR(34)*S
      		PBPM(6,6) = PBPM(6,6)+S**2*DR(34)-R*DR(35)*S
      		PBPM(6,7) = PBPM(6,7)+S*DR(35)
      		PBPM(7,1) = PBPM(7,1)+DR(41)
      		PBPM(7,2) = PBPM(7,2)+R*DR(41)-S*DR(43)
      		PBPM(7,3) = PBPM(7,3)+S*DR(41)
      		PBPM(7,4) = PBPM(7,4)+DR(42)
      		PBPM(7,5) = PBPM(7,5)+R*DR(42)
      		PBPM(7,6) = PBPM(7,6)+S*DR(42)-R*DR(43)
      		PBPM(7,7) = PBPM(7,7)+DR(43)
      		PBPM(8,1) = PBPM(8,1)+R*DR(41)
      		PBPM(8,2) = PBPM(8,2)+R**2*DR(41)-R*DR(43)*S
      		PBPM(8,3) = PBPM(8,3)+R*DR(41)*S
      		PBPM(8,4) = PBPM(8,4)+R*DR(42)
      		PBPM(8,5) = PBPM(8,5)+R**2*DR(42)
      		PBPM(8,6) = PBPM(8,6)+R*DR(42)*S-R**2*DR(43)
      		PBPM(8,7) = PBPM(8,7)+R*DR(43)
      		PBPM(9,1) = PBPM(9,1)+S*DR(41)
      		PBPM(9,2) = PBPM(9,2)+R*DR(41)*S-S**2*DR(43)
      		PBPM(9,3) = PBPM(9,3)+S**2*DR(41)
      		PBPM(9,4) = PBPM(9,4)+S*DR(42)
      		PBPM(9,5) = PBPM(9,5)+R*DR(42)*S
      		PBPM(9,6) = PBPM(9,6)+S**2*DR(42)-R*DR(43)*S
      		PBPM(9,7) = PBPM(9,7)+S*DR(43)
	ENDIF

C	---------------------------
C	COMPUTE PBPB - BENDING PART
C	---------------------------
	PBPB(1,1) = PBPB(1,1)+DR(28)
      PBPB(1,2) = PBPB(1,2)+R*DR(28)
      PBPB(1,3) = PBPB(1,3)+S*DR(28)
      PBPB(1,4) = PBPB(1,4)+DR(29)
      PBPB(1,5) = PBPB(1,5)+R*DR(29)
      PBPB(1,6) = PBPB(1,6)+S*DR(29)
      PBPB(1,7) = PBPB(1,7)+DR(30)
      PBPB(1,8) = PBPB(1,8)+R*DR(30)
      PBPB(1,9) = PBPB(1,9)+S*DR(30)
      PBPB(2,1) = PBPB(2,1)+R*DR(28)
      PBPB(2,2) = PBPB(2,2)+R**2*DR(28)
      PBPB(2,3) = PBPB(2,3)+R*DR(28)*S
      PBPB(2,4) = PBPB(2,4)+R*DR(29)
      PBPB(2,5) = PBPB(2,5)+R**2*DR(29)
      PBPB(2,6) = PBPB(2,6)+R*DR(29)*S
      PBPB(2,7) = PBPB(2,7)+R*DR(30)
      PBPB(2,8) = PBPB(2,8)+R**2*DR(30)
      PBPB(2,9) = PBPB(2,9)+R*DR(30)*S
      PBPB(3,1) = PBPB(3,1)+S*DR(28)
      PBPB(3,2) = PBPB(3,2)+R*DR(28)*S
      PBPB(3,3) = PBPB(3,3)+S**2*DR(28)
      PBPB(3,4) = PBPB(3,4)+S*DR(29)
      PBPB(3,5) = PBPB(3,5)+R*DR(29)*S
      PBPB(3,6) = PBPB(3,6)+S**2*DR(29)
      PBPB(3,7) = PBPB(3,7)+S*DR(30)
      PBPB(3,8) = PBPB(3,8)+R*DR(30)*S
      PBPB(3,9) = PBPB(3,9)+S**2*DR(30)
      PBPB(4,1) = PBPB(4,1)+DR(36)
      PBPB(4,2) = PBPB(4,2)+R*DR(36)
      PBPB(4,3) = PBPB(4,3)+S*DR(36)
      PBPB(4,4) = PBPB(4,4)+DR(37)
      PBPB(4,5) = PBPB(4,5)+R*DR(37)
      PBPB(4,6) = PBPB(4,6)+S*DR(37)
      PBPB(4,7) = PBPB(4,7)+DR(38)
      PBPB(4,8) = PBPB(4,8)+R*DR(38)
      PBPB(4,9) = PBPB(4,9)+S*DR(38)
      PBPB(5,1) = PBPB(5,1)+R*DR(36)
      PBPB(5,2) = PBPB(5,2)+R**2*DR(36)
      PBPB(5,3) = PBPB(5,3)+R*DR(36)*S
      PBPB(5,4) = PBPB(5,4)+R*DR(37)
      PBPB(5,5) = PBPB(5,5)+R**2*DR(37)
      PBPB(5,6) = PBPB(5,6)+R*DR(37)*S
      PBPB(5,7) = PBPB(5,7)+R*DR(38)
      PBPB(5,8) = PBPB(5,8)+R**2*DR(38)
      PBPB(5,9) = PBPB(5,9)+R*DR(38)*S
      PBPB(6,1) = PBPB(6,1)+S*DR(36)
      PBPB(6,2) = PBPB(6,2)+R*DR(36)*S
      PBPB(6,3) = PBPB(6,3)+S**2*DR(36)
      PBPB(6,4) = PBPB(6,4)+S*DR(37)
      PBPB(6,5) = PBPB(6,5)+R*DR(37)*S
      PBPB(6,6) = PBPB(6,6)+S**2*DR(37)
      PBPB(6,7) = PBPB(6,7)+S*DR(38)
      PBPB(6,8) = PBPB(6,8)+R*DR(38)*S
      PBPB(6,9) = PBPB(6,9)+S**2*DR(38)
      PBPB(7,1) = PBPB(7,1)+DR(44)
      PBPB(7,2) = PBPB(7,2)+R*DR(44)
      PBPB(7,3) = PBPB(7,3)+S*DR(44)
      PBPB(7,4) = PBPB(7,4)+DR(45)
      PBPB(7,5) = PBPB(7,5)+R*DR(45)
      PBPB(7,6) = PBPB(7,6)+S*DR(45)
      PBPB(7,7) = PBPB(7,7)+DR(46)
      PBPB(7,8) = PBPB(7,8)+R*DR(46)
      PBPB(7,9) = PBPB(7,9)+S*DR(46)
      PBPB(8,1) = PBPB(8,1)+R*DR(44)
      PBPB(8,2) = PBPB(8,2)+R**2*DR(44)
      PBPB(8,3) = PBPB(8,3)+R*DR(44)*S
      PBPB(8,4) = PBPB(8,4)+R*DR(45)
      PBPB(8,5) = PBPB(8,5)+R**2*DR(45)
      PBPB(8,6) = PBPB(8,6)+R*DR(45)*S
      PBPB(8,7) = PBPB(8,7)+R*DR(46)
      PBPB(8,8) = PBPB(8,8)+R**2*DR(46)
      PBPB(8,9) = PBPB(8,9)+R*DR(46)*S
      PBPB(9,1) = PBPB(9,1)+S*DR(44)
      PBPB(9,2) = PBPB(9,2)+R*DR(44)*S
      PBPB(9,3) = PBPB(9,3)+S**2*DR(44)
      PBPB(9,4) = PBPB(9,4)+S*DR(45)
      PBPB(9,5) = PBPB(9,5)+R*DR(45)*S
      PBPB(9,6) = PBPB(9,6)+S**2*DR(45)
      PBPB(9,7) = PBPB(9,7)+S*DR(46)
      PBPB(9,8) = PBPB(9,8)+R*DR(46)*S
	PBPB(9,9) = PBPB(9,9)+S**2*DR(46)

C	-----------------------------------
C	COMPUTE PTP - TRANSVERSE SHEAR PART
C	-----------------------------------
      PTP(1,1) = PTP(1,1)+DR(55)
      PTP(1,2) = PTP(1,2)+DR(56)
      PTP(2,1) = PTP(2,1)+DR(63)
      PTP(2,2) = PTP(2,2)+DR(64)

      RETURN
      END
C
C=====================================================================
      SUBROUTINE PSIGRT(SIGR,R,S,PAPM,PAPB,PAPSH)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTES FOR INTEGRAL(TRANS(P)*SIGR)dA
C
C	INPUT VARIABLES
C	SIGR(8)  = NODAL STRESS RESULTANTS
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C	PAPM(5)  = INTEGRAL(TRANS(Pm)*N)
C	PAPB(11) = INTEGRAL(TRANS(Pb)*M)
C	PAPSH(2) = INTEGRAL(TRANS(Ps)*Q)
C     ----------------------------------------------------------------
	DIMENSION SIGR(8)
	DIMENSION PAPM(7),PAPB(9),PAPSH(2)

C	----------------------------
C	COMPUTE PAPM - MEMBRANE PART
C	----------------------------
      PAPM(1) = PAPM(1)+SIGR(1)
      PAPM(2) = PAPM(2)+R*SIGR(1)-S*SIGR(3)
      PAPM(3) = PAPM(3)+S*SIGR(1)
      PAPM(4) = PAPM(4)+SIGR(2)
      PAPM(5) = PAPM(5)+R*SIGR(2)
      PAPM(6) = PAPM(6)+S*SIGR(2)-R*SIGR(3)
      PAPM(7) = PAPM(7)+SIGR(3)

C	---------------------------
C	COMPUTE PAPB - BENDING PART
C	---------------------------
      PAPB(1) = PAPB(1)+SIGR(4)
      PAPB(2) = PAPB(2)+R*SIGR(4)
      PAPB(3) = PAPB(3)+S*SIGR(4)
      PAPB(4) = PAPB(4)+SIGR(5)
      PAPB(5) = PAPB(5)+R*SIGR(5)
      PAPB(6) = PAPB(6)+S*SIGR(5)
      PAPB(7) = PAPB(7)+SIGR(6)
      PAPB(8) = PAPB(8)+R*SIGR(6)
      PAPB(9) = PAPB(9)+S*SIGR(6)

C	-------------------------------------
C	COMPUTE PAPSH - TRANSVERSE SHEAR PART
C	-------------------------------------
      PAPSH(1) = PAPSH(1)+SIGR(7)
      PAPSH(2) = PAPSH(2)+SIGR(8)

      RETURN
      END
C
C=====================================================================
      SUBROUTINE PFPT(SIGR,R,S,PFP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTES FOR INTEGRAL(TRANS(P)*SIGR)dA
C
C	INPUT VARIABLES
C	SIGR(8)  = NODAL STRESS RESULTANTS
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C	PAPM(5)  = INTEGRAL(TRANS(Pm)*N)
C	PAPB(11) = INTEGRAL(TRANS(Pb)*M)
C	PAPSH(2) = INTEGRAL(TRANS(Ps)*Q)
C     ----------------------------------------------------------------
	DIMENSION SIGR(8),PFP(15,15)

      PFP(1,1) = PFP(1,1)+SIGR(1)
      PFP(1,2) = PFP(1,2)+SIGR(3)
      PFP(1,9) = PFP(1,9)+SIGR(4)
      PFP(1,10) = PFP(1,10)+SIGR(6)
      PFP(1,14) = PFP(1,14)+SIGR(7)

      PFP(2,1) = PFP(2,1)+SIGR(3)
      PFP(2,2) = PFP(2,2)+SIGR(2)
      PFP(2,9) = PFP(2,9)+SIGR(6)
      PFP(2,10) = PFP(2,10)+SIGR(5)
      PFP(2,14) = PFP(2,14)+SIGR(8)

      PFP(3,3) = PFP(3,3)+SIGR(1)
      PFP(3,4) = PFP(3,4)+SIGR(3)
      PFP(3,7) = PFP(3,7)-SIGR(4)
      PFP(3,8) = PFP(3,8)-SIGR(6)
      PFP(3,13) = PFP(3,13)-SIGR(7)

      PFP(4,3) = PFP(4,3)+SIGR(3)
      PFP(4,4) = PFP(4,4)+SIGR(2)
      PFP(4,7) = PFP(4,7)-SIGR(6)
      PFP(4,8) = PFP(4,8)-SIGR(5)
      PFP(4,13) = PFP(4,13)-SIGR(8)

      PFP(5,5) = PFP(5,5)+SIGR(1)
      PFP(5,6) = PFP(5,6)+SIGR(3)

      PFP(6,5) = PFP(6,5)+SIGR(3)
      PFP(6,6) = PFP(6,6)+SIGR(2)

      PFP(7,3) = PFP(7,3)-SIGR(4)
      PFP(7,4) = PFP(7,4)-SIGR(6)
      PFP(7,15) = PFP(7,15)+SIGR(4)/2.0D0

      PFP(8,3) = PFP(8,3)-SIGR(6)
      PFP(8,4) = PFP(8,4)-SIGR(5)
      PFP(8,15) = PFP(8,15)+SIGR(6)/2.0D0

      PFP(9,1) = PFP(9,1)+SIGR(4)
      PFP(9,2) = PFP(9,2)+SIGR(6)
      PFP(9,15) = PFP(9,15)+SIGR(6)/2.0D0

      PFP(10,1) = PFP(10,1)+SIGR(6)
      PFP(10,2) = PFP(10,2)+SIGR(5)
      PFP(10,15) = PFP(10,15)+SIGR(5)/2.0D0

      PFP(11,13) = PFP(11,13)+SIGR(4)/2.0D0
      PFP(11,14) = PFP(11,14)+SIGR(6)/2.0D0

      PFP(12,13) = PFP(12,13)+SIGR(6)/2.0D0
      PFP(12,14) = PFP(12,14)+SIGR(5)/2.0D0

      PFP(13,3) = PFP(13,3)-SIGR(7)
      PFP(13,4) = PFP(13,4)-SIGR(8)
      PFP(13,11) = PFP(13,11)+SIGR(4)/2.0D0
      PFP(13,12) = PFP(13,12)+SIGR(6)/2.0D0
      PFP(13,15) = PFP(13,15)+SIGR(7)/2.0D0

      PFP(14,1) = PFP(14,1)+SIGR(7)
      PFP(14,2) = PFP(14,2)+SIGR(8)
      PFP(14,11) = PFP(14,11)+SIGR(6)/2.0D0
      PFP(14,12) = PFP(14,12)+SIGR(5)/2.0D0
      PFP(14,15) = PFP(14,15)+SIGR(8)/2.0D0

      PFP(15,7) = PFP(15,7)+SIGR(4)/2.0D0
      PFP(15,8) = PFP(15,8)+SIGR(6)/2.0D0
      PFP(15,9) = PFP(15,9)+SIGR(6)/2.0D0
      PFP(15,10) = PFP(15,10)+SIGR(5)/2.0D0
      PFP(15,13) = PFP(15,13)+SIGR(7)/2.0D0
      PFP(15,14) = PFP(15,14)+SIGR(8)/2.0D0

      RETURN
      END
C
C=====================================================================
      SUBROUTINE CGEOT(AB1,ANT,RN,SN,CB,CM,CS,CSC,RS,SLN,BL,
	1                DWR,DWS,FLR,FLS,ACG)

	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cs FOR TRANSVERSE SHEAR STRAINS
C
C	INPUT VARIABLES
C	ANT(4)		= AREA INTEGRATION FACTORS
C	CM(9,24)	= MEMBRANE AND DRILLING STRAIN PARAMETER
C	CB(11,24)	= BENDING STRAIN PARAMETER
C	CS(2,24)	= TRANSVERSE SHEAR STRAIN PARAMETER
C	CSC(4,24)	= COMPONENTS OF CS
C				= CSC(1,*) -> LINE INTEGRAL (DW/Dr)
C				= CSC(2,*) -> AREA INTEGRAL D(PHI)s
C				= CSC(3,*) -> INT(DW/Ds)
C				= CSC(4,*) -> AREA INTEGRAL -D(PHI)r
C	RS(2,4)		= ELEMENT COORDINATES IN LOCAL SPACE
C	RN(4),SN(4)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	SLN(4)		= BOUNDARY LENGTHS
C	BL,BLR,BLS	= SHEAR DEFORMATION PARAMETERS ALONG THE BOUNDARY,
C                     R, S, RESP.
C
C	LOACAL VARIABLES
C	ANR(4)		= AREA INTEGRATION WITH FACTOR R
C	ANS(4)		= AREA INTEGRATION WITH FACTOR S
C
C	OUTPUT VARIAVBLES
C	ACG(51,24)	= inverse(Ag)*(Cg)
C	--------------------------------------------------------------
	DIMENSION AB1(3,3),ANT(3),RN(3),SN(3)
	DIMENSION CB(9,18),CM(7,18),CS(2,18),CSC(4,18),RS(2,3)
      DIMENSION SLN(3),BL(3)
	DIMENSION DWR(9),DWS(9),FLR(3),FLS(3)
	DIMENSION ACG(15,18)

C	--------------------
C	COMPONENTS FOR DU/DR
C	--------------------
	ACG(1,1) =CM(1,1)/AB1(1,1)
	ACG(1,7) =CM(1,7)/AB1(1,1)
	ACG(1,13)=CM(1,13)/AB1(1,1)

	ACG(1,6) =CM(1,6)/AB1(1,1)
	ACG(1,12)=CM(1,12)/AB1(1,1)
	ACG(1,18)=CM(1,18)/AB1(1,1)

C	--------------------
C	COMPONENTS FOR DU/DS
C	--------------------
	ACG(2,1) =(SLN(1)*SN(1)/2+SLN(3)*SN(3)/2)/AB1(1,1) !*du(1)
	ACG(2,7) =(SLN(2)*SN(2)/2+SLN(1)*SN(1)/2)/AB1(1,1) !*du(2)
	ACG(2,13)=(SLN(2)*SN(2)/2+SLN(3)*SN(3)/2)/AB1(1,1) !*du(3)

	ACG(2,6) =((-RS(2,2)/12+RS(2,1)/12)*SN(1)*SLN(1)+
     #           (RS(2,1)/12-RS(2,3)/12)*SN(3)*SLN(3))/AB1(1,1) !*dpt(1)
	ACG(2,12)=((-RS(2,1)/12+RS(2,2)/12)*SN(1)*SLN(1)+
     #           (RS(2,2)/12-RS(2,3)/12)*SN(2)*SLN(2))/AB1(1,1) !*dpt(2)
	ACG(2,18)=((-RS(2,2)/12+RS(2,3)/12)*SN(2)*SLN(2)+
     #           (RS(2,3)/12-RS(2,1)/12)*SN(3)*SLN(3))/AB1(1,1) !*dpt(3)

C	--------------------
C	COMPONENTS FOR DV/DR
C	--------------------
	ACG(3,2) =(SLN(1)*RN(1)/2+SLN(3)*RN(3)/2)/AB1(1,1) !*dv(1)
	ACG(3,8) =(SLN(2)*RN(2)/2+SLN(1)*RN(1)/2)/AB1(1,1) !*dv(2)
	ACG(3,14)=(SLN(2)*RN(2)/2+SLN(3)*RN(3)/2)/AB1(1,1) !*dv(3)

	ACG(3,6) =((RS(1,2)/12-RS(1,1)/12)*RN(1)*SLN(1)+
     #           (-RS(1,1)/12+RS(1,3)/12)*RN(3)*SLN(3))/AB1(1,1) !*dpt(1)
	ACG(3,12)=((RS(1,1)/12-RS(1,2)/12)*RN(1)*SLN(1)+
     #           (-RS(1,2)/12+RS(1,3)/12)*RN(2)*SLN(2))/AB1(1,1) !*dpt(2)
	ACG(3,18)=((RS(1,2)/12-RS(1,3)/12)*RN(2)*SLN(2)+
     #           (-RS(1,3)/12+RS(1,1)/12)*RN(3)*SLN(3))/AB1(1,1) !*dpt(3)

C	--------------------
C	COMPONENTS FOR DV/DS
C	--------------------
	ACG(4,2) =CM(4,2)/AB1(1,1)
	ACG(4,8) =CM(4,8)/AB1(1,1)
	ACG(4,14)=CM(4,14)/AB1(1,1)

	ACG(4,6) =CM(4,6)/AB1(1,1)
	ACG(4,12)=CM(4,12)/AB1(1,1)
	ACG(4,18)=CM(4,18)/AB1(1,1)

C	--------------------
C	COMPONENTS FOR DW/DR
C	--------------------
	ACG(5,4) =CSC(1,4)/AB1(1,1)
	ACG(5,10)=CSC(1,10)/AB1(1,1)
	ACG(5,16)=CSC(1,16)/AB1(1,1)

	ACG(5,5) =CSC(1,5)/AB1(1,1)
	ACG(5,11)=CSC(1,11)/AB1(1,1)
	ACG(5,17)=CSC(1,17)/AB1(1,1)

	ACG(5,3) =CSC(1,3)/AB1(1,1)
	ACG(5,9) =CSC(1,9)/AB1(1,1)
	ACG(5,15)=CSC(1,15)/AB1(1,1)

C	--------------------
C	COMPONENTS FOR DW/DS
C	--------------------
	ACG(6,4) =CSC(3,4)/AB1(1,1)
	ACG(6,10)=CSC(3,10)/AB1(1,1)
	ACG(6,16)=CSC(3,16)/AB1(1,1)

	ACG(6,5) =CSC(3,5)/AB1(1,1)
	ACG(6,11)=CSC(3,11)/AB1(1,1)
	ACG(6,17)=CSC(3,17)/AB1(1,1)

	ACG(6,3) =CSC(3,3)/AB1(1,1)
	ACG(6,9) =CSC(3,9)/AB1(1,1)
	ACG(6,15)=CSC(3,15)/AB1(1,1)

C	---------------------
C	COMPONENTS FOR DJr/DR
C	---------------------
	ACG(7,4) =
     #((-BL(1)/2+1.E0/2.E0)*SLN(1)*RN(1)**3+SLN(1)*RN(1)*SN(1)**2/2
     #+(1.E0/2.E0-BL(3)/2)*SLN(3)*RN(3)**3+SLN(3)*RN(3)*SN(3)**2/2) !*dpr(1)
     #/AB1(1,1)
	ACG(7,10)=
     #((-BL(1)/2+1.E0/2.E0)*SLN(1)*RN(1)**3+SLN(1)*RN(1)*SN(1)**2/2
     #+(1.E0/2.E0-BL(2)/2)*SLN(2)*RN(2)**3+SLN(2)*RN(2)*SN(2)**2/2) !*dpr(2)
     #/AB1(1,1)
	ACG(7,16)= 
     #((1.E0/2.E0-BL(2)/2)*SLN(2)*RN(2)**3+SLN(2)*RN(2)*SN(2)**2/2
     #+(1.E0/2.E0-BL(3)/2)*SLN(3)*RN(3)**3+SLN(3)*RN(3)*SN(3)**2/2) !*dpr(3)
     #/AB1(1,1)

	ACG(7,5) =
     #(-SLN(3)*RN(3)**2*BL(3)*SN(3)/2-SLN(1)*RN(1)**2*BL(1)*SN(1)/2) !*dps(1)
     #/AB1(1,1)
	ACG(7,11)=
     #(-SLN(1)*RN(1)**2*BL(1)*SN(1)/2-SLN(2)*RN(2)**2*BL(2)*SN(2)/2) !*dps(2)
     #/AB1(1,1)
	ACG(7,17)=
     #(-SLN(3)*RN(3)**2*BL(3)*SN(3)/2-SLN(2)*RN(2)**2*BL(2)*SN(2)/2) !*dps(3)
     #/AB1(1,1)

	ACG(7,3) =(-RN(1)**2*BL(1)+RN(3)**2*BL(3))/AB1(1,1) !*dw(1)
	ACG(7,9) =(-RN(2)**2*BL(2)+RN(1)**2*BL(1))/AB1(1,1) !*dw(2)
	ACG(7,15)=( RN(2)**2*BL(2)-RN(3)**2*BL(3))/AB1(1,1) !*dw(3)

C	---------------------
C	COMPONENTS FOR DJr/DS
C	---------------------
	ACG(8,4) =-CB(4,4)/AB1(1,1)
	ACG(8,10)=-CB(4,10)/AB1(1,1)
	ACG(8,16)=-CB(4,16)/AB1(1,1)

	ACG(8,5) =-CB(4,5)/AB1(1,1)
	ACG(8,11)=-CB(4,11)/AB1(1,1)
	ACG(8,17)=-CB(4,17)/AB1(1,1)

	ACG(8,3) =-CB(4,3)/AB1(1,1)
	ACG(8,9) =-CB(4,9)/AB1(1,1)
	ACG(8,15)=-CB(4,15)/AB1(1,1)

C	---------------------
C	COMPONENTS FOR DJs/DR
C	---------------------
	ACG(9,4) =CB(1,4)/AB1(1,1)
	ACG(9,10)=CB(1,10)/AB1(1,1)
	ACG(9,16)=CB(1,16)/AB1(1,1)

	ACG(9,5) =CB(1,5)/AB1(1,1)
	ACG(9,11)=CB(1,11)/AB1(1,1)
	ACG(9,17)=CB(1,17)/AB1(1,1)

	ACG(9,3) =CB(1,3)/AB1(1,1)
	ACG(9,9) =CB(1,9)/AB1(1,1)
	ACG(9,15)=CB(1,15)/AB1(1,1)

C	---------------------
C	COMPONENTS FOR DJs/DS
C	---------------------
	ACG(10,4) = 
     #(-SLN(3)*SN(3)**2*BL(3)*RN(3)/2-SLN(1)*SN(1)**2*BL(1)*RN(1)/2) !*dpr(1)
     #/AB1(1,1)
	ACG(10,10)= 
     #(-SLN(1)*SN(1)**2*BL(1)*RN(1)/2-SLN(2)*SN(2)**2*BL(2)*RN(2)/2) !*dpr(2)
     #/AB1(1,1)
	ACG(10,16)= 
     #(-SLN(2)*SN(2)**2*BL(2)*RN(2)/2-SLN(3)*SN(3)**2*BL(3)*RN(3)/2) !*dpr(3)
     #/AB1(1,1)


	ACG(10,5) = 
     #(SLN(1)*RN(1)**2*SN(1)/2+SLN(3)*RN(3)**2*SN(3)/2+
     #(-BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3+
     #(1.E0/2.E0-BL(3)/2)*SLN(3)*SN(3)**3) !*dps(1)
     #/AB1(1,1)
	ACG(10,11)= 
     #(SLN(1)*RN(1)**2*SN(1)/2+SLN(2)*RN(2)**2*SN(2)/2+
     #(-BL(1)/2+1.E0/2.E0)*SLN(1)*SN(1)**3+
     #(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3) !*dps(2)
     #/AB1(1,1)
	ACG(10,17)= 
     #(SLN(2)*RN(2)**2*SN(2)/2+SLN(3)*RN(3)**2*SN(3)/2+
     #(1.E0/2.E0-BL(2)/2)*SLN(2)*SN(2)**3+
     #(1.E0/2.E0-BL(3)/2)*SLN(3)*SN(3)**3) !*dps(3)
     #/AB1(1,1)


	ACG(10,3) = (-SN(1)**2*BL(1)+SN(3)**2*BL(3))/AB1(1,1) !dw(1)
	ACG(10,9) = (-SN(2)**2*BL(2)+SN(1)**2*BL(1))/AB1(1,1) !*dw(2)
	ACG(10,15)= ( SN(2)**2*BL(2)-SN(3)**2*BL(3))/AB1(1,1) !*dw(3)

C	---------------------
C	COMPONENTS FOR DJt/DR
C	---------------------
	ACG(11,6) =CM(1,1)/AB1(1,1)
	ACG(11,12)=CM(1,7)/AB1(1,1)
	ACG(11,18)=CM(1,13)/AB1(1,1)

C	---------------------
C	COMPONENTS FOR DJt/DS
C	---------------------
	ACG(12,6) =CM(4,2)/AB1(1,1)
	ACG(12,12)=CM(4,8)/AB1(1,1)
	ACG(12,18)=CM(4,14)/AB1(1,1)

C	-----------------
C	COMPONENTS FOR Jr
C	-----------------
	ACG(13,4) =(-DWS(1)-FLR(1))/AB1(1,1)
	ACG(13,10)=(-DWS(2)-FLR(2))/AB1(1,1)
	ACG(13,16)=(-DWS(3)-FLR(3))/AB1(1,1)

	ACG(13,5) = -DWS(4)/AB1(1,1)
	ACG(13,11)= -DWS(5)/AB1(1,1)
	ACG(13,17)= -DWS(6)/AB1(1,1)

	ACG(13,3) = -DWS(7)/AB1(1,1)
	ACG(13,9) = -DWS(8)/AB1(1,1)
	ACG(13,15)= -DWS(9)/AB1(1,1)

C	-----------------
C	COMPONENTS FOR Js
C	-----------------
	ACG(14,4) =DWR(1)/AB1(1,1)
	ACG(14,10)=DWR(2)/AB1(1,1)
	ACG(14,16)=DWR(3)/AB1(1,1)

	ACG(14,5) =(DWR(4)+FLS(1))/AB1(1,1)
	ACG(14,11)=(DWR(5)+FLS(2))/AB1(1,1)
	ACG(14,17)=(DWR(6)+FLS(3))/AB1(1,1)

	ACG(14,3) =DWR(7)/AB1(1,1)
	ACG(14,9) =DWR(8)/AB1(1,1)
	ACG(14,15)=DWR(9)/AB1(1,1)

C	-----------------
C	COMPONENTS FOR Jt
C	-----------------
	ACG(15,6) =ANT(1)/AB1(1,1)
	ACG(15,12)=ANT(2)/AB1(1,1)
	ACG(15,18)=ANT(3)/AB1(1,1)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE KGEO18(ACG,PFP,ST)

	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR THE GEOMETRIC STIFFNESS - S
C
C	INPUT VARIABLES
C	ACG(51,24)	= INVERSE(Ag)*(Cg)*TRANSPOSE(T)
C	PFP(51,51)	= Ng*F*Pg
C
C	LOACAL VARIABLES
C
C	OUTPUT VARIAVBLES
C	ST(*)		= GEOMETRIC STIFFNESS - UPPER TRIANG - VECTOR FORM
C	--------------------------------------------------------------
      DIMENSION ACG(15,18),PFP(15,15)
	DIMENSION ST(171)

      ST(1) = ST(1)+(ACG(1,1)*PFP(1,1)+ACG(2,1)*PFP(1,2))*ACG(1,1)+(ACG(
     #1,1)*PFP(1,2)+ACG(2,1)*PFP(2,2))*ACG(2,1)
      ST(2) = ST(2)
      ST(3) = ST(3)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,3)+(ACG(
     #1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,3)+(ACG(1,1)*PFP(1,14)+A
     #CG(2,1)*PFP(2,14))*ACG(14,3)
      ST(4) = ST(4)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,4)+(ACG(
     #1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,4)+(ACG(1,1)*PFP(1,14)+A
     #CG(2,1)*PFP(2,14))*ACG(14,4)
      ST(5) = ST(5)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,5)+(ACG(
     #1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,5)+(ACG(1,1)*PFP(1,14)+A
     #CG(2,1)*PFP(2,14))*ACG(14,5)
      ST(6) = ST(6)+(ACG(1,1)*PFP(1,1)+ACG(2,1)*PFP(1,2))*ACG(1,6)+(ACG(
     #1,1)*PFP(1,2)+ACG(2,1)*PFP(2,2))*ACG(2,6)
      ST(7) = ST(7)+(ACG(1,1)*PFP(1,1)+ACG(2,1)*PFP(1,2))*ACG(1,7)+(ACG(
     #1,1)*PFP(1,2)+ACG(2,1)*PFP(2,2))*ACG(2,7)
      ST(8) = ST(8)
      ST(9) = ST(9)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,9)+(ACG(
     #1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,9)+(ACG(1,1)*PFP(1,14)+A
     #CG(2,1)*PFP(2,14))*ACG(14,9)
      ST(10) = ST(10)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,10)+(A
     #CG(1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,10)+(ACG(1,1)*PFP(1,1
     #4)+ACG(2,1)*PFP(2,14))*ACG(14,10)
      ST(11) = ST(11)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,11)+(A
     #CG(1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,11)+(ACG(1,1)*PFP(1,1
     #4)+ACG(2,1)*PFP(2,14))*ACG(14,11)
      ST(12) = ST(12)+(ACG(1,1)*PFP(1,1)+ACG(2,1)*PFP(1,2))*ACG(1,12)+(A
     #CG(1,1)*PFP(1,2)+ACG(2,1)*PFP(2,2))*ACG(2,12)
      ST(13) = ST(13)+(ACG(1,1)*PFP(1,1)+ACG(2,1)*PFP(1,2))*ACG(1,13)+(A
     #CG(1,1)*PFP(1,2)+ACG(2,1)*PFP(2,2))*ACG(2,13)
      ST(14) = ST(14)
      ST(15) = ST(15)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,15)+(A
     #CG(1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,15)+(ACG(1,1)*PFP(1,1
     #4)+ACG(2,1)*PFP(2,14))*ACG(14,15)
      ST(16) = ST(16)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,16)+(A
     #CG(1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,16)+(ACG(1,1)*PFP(1,1
     #4)+ACG(2,1)*PFP(2,14))*ACG(14,16)
      ST(17) = ST(17)+(ACG(1,1)*PFP(1,9)+ACG(2,1)*PFP(2,9))*ACG(9,17)+(A
     #CG(1,1)*PFP(1,10)+ACG(2,1)*PFP(2,10))*ACG(10,17)+(ACG(1,1)*PFP(1,1
     #4)+ACG(2,1)*PFP(2,14))*ACG(14,17)
      ST(18) = ST(18)+(ACG(1,1)*PFP(1,1)+ACG(2,1)*PFP(1,2))*ACG(1,18)+(A
     #CG(1,1)*PFP(1,2)+ACG(2,1)*PFP(2,2))*ACG(2,18)
      ST(19) = ST(19)+(ACG(3,2)*PFP(3,3)+ACG(4,2)*PFP(3,4))*ACG(3,2)+(AC
     #G(3,2)*PFP(3,4)+ACG(4,2)*PFP(4,4))*ACG(4,2)
      ST(20) = ST(20)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,3)+(AC
     #G(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,3)+(ACG(3,2)*PFP(3,13)+AC
     #G(4,2)*PFP(4,13))*ACG(13,3)
      ST(21) = ST(21)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,4)+(AC
     #G(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,4)+(ACG(3,2)*PFP(3,13)+AC
     #G(4,2)*PFP(4,13))*ACG(13,4)
      ST(22) = ST(22)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,5)+(AC
     #G(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,5)+(ACG(3,2)*PFP(3,13)+AC
     #G(4,2)*PFP(4,13))*ACG(13,5)
      ST(23) = ST(23)+(ACG(3,2)*PFP(3,3)+ACG(4,2)*PFP(3,4))*ACG(3,6)+(AC
     #G(3,2)*PFP(3,4)+ACG(4,2)*PFP(4,4))*ACG(4,6)
      ST(24) = ST(24)
      ST(25) = ST(25)+(ACG(3,2)*PFP(3,3)+ACG(4,2)*PFP(3,4))*ACG(3,8)+(AC
     #G(3,2)*PFP(3,4)+ACG(4,2)*PFP(4,4))*ACG(4,8)
      ST(26) = ST(26)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,9)+(AC
     #G(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,9)+(ACG(3,2)*PFP(3,13)+AC
     #G(4,2)*PFP(4,13))*ACG(13,9)
      ST(27) = ST(27)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,10)+(A
     #CG(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,10)+(ACG(3,2)*PFP(3,13)+
     #ACG(4,2)*PFP(4,13))*ACG(13,10)
      ST(28) = ST(28)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,11)+(A
     #CG(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,11)+(ACG(3,2)*PFP(3,13)+
     #ACG(4,2)*PFP(4,13))*ACG(13,11)
      ST(29) = ST(29)+(ACG(3,2)*PFP(3,3)+ACG(4,2)*PFP(3,4))*ACG(3,12)+(A
     #CG(3,2)*PFP(3,4)+ACG(4,2)*PFP(4,4))*ACG(4,12)
      ST(30) = ST(30)
      ST(31) = ST(31)+(ACG(3,2)*PFP(3,3)+ACG(4,2)*PFP(3,4))*ACG(3,14)+(A
     #CG(3,2)*PFP(3,4)+ACG(4,2)*PFP(4,4))*ACG(4,14)
      ST(32) = ST(32)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,15)+(A
     #CG(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,15)+(ACG(3,2)*PFP(3,13)+
     #ACG(4,2)*PFP(4,13))*ACG(13,15)
      ST(33) = ST(33)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,16)+(A
     #CG(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,16)+(ACG(3,2)*PFP(3,13)+
     #ACG(4,2)*PFP(4,13))*ACG(13,16)
      ST(34) = ST(34)+(ACG(3,2)*PFP(3,7)+ACG(4,2)*PFP(4,7))*ACG(7,17)+(A
     #CG(3,2)*PFP(3,8)+ACG(4,2)*PFP(4,8))*ACG(8,17)+(ACG(3,2)*PFP(3,13)+
     #ACG(4,2)*PFP(4,13))*ACG(13,17)
      ST(35) = ST(35)+(ACG(3,2)*PFP(3,3)+ACG(4,2)*PFP(3,4))*ACG(3,18)+(A
     #CG(3,2)*PFP(3,4)+ACG(4,2)*PFP(4,4))*ACG(4,18)
      ST(36) = ST(36)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,3)+(AC
     #G(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,3)
      ST(37) = ST(37)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,4)+(AC
     #G(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,4)
      ST(38) = ST(38)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,5)+(AC
     #G(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,5)
      ST(39) = ST(39)+(ACG(9,3)*PFP(1,9)+ACG(10,3)*PFP(1,10)+ACG(14,3)*P
     #FP(1,14))*ACG(1,6)+(ACG(9,3)*PFP(2,9)+ACG(10,3)*PFP(2,10)+ACG(14,3
     #)*PFP(2,14))*ACG(2,6)+(ACG(7,3)*PFP(3,7)+ACG(8,3)*PFP(3,8)+ACG(13,
     #3)*PFP(3,13))*ACG(3,6)+(ACG(7,3)*PFP(4,7)+ACG(8,3)*PFP(4,8)+ACG(13
     #,3)*PFP(4,13))*ACG(4,6)+(ACG(13,3)*PFP(11,13)+ACG(14,3)*PFP(11,14)
     #)*ACG(11,6)+(ACG(13,3)*PFP(12,13)+ACG(14,3)*PFP(12,14))*ACG(12,6)+
     #(ACG(7,3)*PFP(7,15)+ACG(8,3)*PFP(8,15)+ACG(9,3)*PFP(9,15)+ACG(10,3
     #)*PFP(10,15)+ACG(13,3)*PFP(13,15)+ACG(14,3)*PFP(14,15))*ACG(15,6)
      ST(40) = ST(40)+(ACG(9,3)*PFP(1,9)+ACG(10,3)*PFP(1,10)+ACG(14,3)*P
     #FP(1,14))*ACG(1,7)+(ACG(9,3)*PFP(2,9)+ACG(10,3)*PFP(2,10)+ACG(14,3
     #)*PFP(2,14))*ACG(2,7)
      ST(41) = ST(41)+(ACG(7,3)*PFP(3,7)+ACG(8,3)*PFP(3,8)+ACG(13,3)*PFP
     #(3,13))*ACG(3,8)+(ACG(7,3)*PFP(4,7)+ACG(8,3)*PFP(4,8)+ACG(13,3)*PF
     #P(4,13))*ACG(4,8)
      ST(42) = ST(42)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,9)+(AC
     #G(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,9)
      ST(43) = ST(43)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,10)+(A
     #CG(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,10)
      ST(44) = ST(44)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,11)+(A
     #CG(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,11)
      ST(45) = ST(45)+(ACG(9,3)*PFP(1,9)+ACG(10,3)*PFP(1,10)+ACG(14,3)*P
     #FP(1,14))*ACG(1,12)+(ACG(9,3)*PFP(2,9)+ACG(10,3)*PFP(2,10)+ACG(14,
     #3)*PFP(2,14))*ACG(2,12)+(ACG(7,3)*PFP(3,7)+ACG(8,3)*PFP(3,8)+ACG(1
     #3,3)*PFP(3,13))*ACG(3,12)+(ACG(7,3)*PFP(4,7)+ACG(8,3)*PFP(4,8)+ACG
     #(13,3)*PFP(4,13))*ACG(4,12)+(ACG(13,3)*PFP(11,13)+ACG(14,3)*PFP(11
     #,14))*ACG(11,12)+(ACG(13,3)*PFP(12,13)+ACG(14,3)*PFP(12,14))*ACG(1
     #2,12)+(ACG(7,3)*PFP(7,15)+ACG(8,3)*PFP(8,15)+ACG(9,3)*PFP(9,15)+AC
     #G(10,3)*PFP(10,15)+ACG(13,3)*PFP(13,15)+ACG(14,3)*PFP(14,15))*ACG(
     #15,12)
      ST(46) = ST(46)+(ACG(9,3)*PFP(1,9)+ACG(10,3)*PFP(1,10)+ACG(14,3)*P
     #FP(1,14))*ACG(1,13)+(ACG(9,3)*PFP(2,9)+ACG(10,3)*PFP(2,10)+ACG(14,
     #3)*PFP(2,14))*ACG(2,13)
      ST(47) = ST(47)+(ACG(7,3)*PFP(3,7)+ACG(8,3)*PFP(3,8)+ACG(13,3)*PFP
     #(3,13))*ACG(3,14)+(ACG(7,3)*PFP(4,7)+ACG(8,3)*PFP(4,8)+ACG(13,3)*P
     #FP(4,13))*ACG(4,14)
      ST(48) = ST(48)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,15)+(A
     #CG(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,15)
      ST(49) = ST(49)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,16)+(A
     #CG(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,16)
      ST(50) = ST(50)+(ACG(5,3)*PFP(5,5)+ACG(6,3)*PFP(5,6))*ACG(5,17)+(A
     #CG(5,3)*PFP(5,6)+ACG(6,3)*PFP(6,6))*ACG(6,17)
      ST(51) = ST(51)+(ACG(9,3)*PFP(1,9)+ACG(10,3)*PFP(1,10)+ACG(14,3)*P
     #FP(1,14))*ACG(1,18)+(ACG(9,3)*PFP(2,9)+ACG(10,3)*PFP(2,10)+ACG(14,
     #3)*PFP(2,14))*ACG(2,18)+(ACG(7,3)*PFP(3,7)+ACG(8,3)*PFP(3,8)+ACG(1
     #3,3)*PFP(3,13))*ACG(3,18)+(ACG(7,3)*PFP(4,7)+ACG(8,3)*PFP(4,8)+ACG
     #(13,3)*PFP(4,13))*ACG(4,18)+(ACG(13,3)*PFP(11,13)+ACG(14,3)*PFP(11
     #,14))*ACG(11,18)+(ACG(13,3)*PFP(12,13)+ACG(14,3)*PFP(12,14))*ACG(1
     #2,18)+(ACG(7,3)*PFP(7,15)+ACG(8,3)*PFP(8,15)+ACG(9,3)*PFP(9,15)+AC
     #G(10,3)*PFP(10,15)+ACG(13,3)*PFP(13,15)+ACG(14,3)*PFP(14,15))*ACG(
     #15,18)
      ST(52) = ST(52)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,4)+(AC
     #G(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,4)
      ST(53) = ST(53)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,5)+(AC
     #G(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,5)
      ST(54) = ST(54)+(ACG(9,4)*PFP(1,9)+ACG(10,4)*PFP(1,10)+ACG(14,4)*P
     #FP(1,14))*ACG(1,6)+(ACG(9,4)*PFP(2,9)+ACG(10,4)*PFP(2,10)+ACG(14,4
     #)*PFP(2,14))*ACG(2,6)+(ACG(7,4)*PFP(3,7)+ACG(8,4)*PFP(3,8)+ACG(13,
     #4)*PFP(3,13))*ACG(3,6)+(ACG(7,4)*PFP(4,7)+ACG(8,4)*PFP(4,8)+ACG(13
     #,4)*PFP(4,13))*ACG(4,6)+(ACG(13,4)*PFP(11,13)+ACG(14,4)*PFP(11,14)
     #)*ACG(11,6)+(ACG(13,4)*PFP(12,13)+ACG(14,4)*PFP(12,14))*ACG(12,6)+
     #(ACG(7,4)*PFP(7,15)+ACG(8,4)*PFP(8,15)+ACG(9,4)*PFP(9,15)+ACG(10,4
     #)*PFP(10,15)+ACG(13,4)*PFP(13,15)+ACG(14,4)*PFP(14,15))*ACG(15,6)
      ST(55) = ST(55)+(ACG(9,4)*PFP(1,9)+ACG(10,4)*PFP(1,10)+ACG(14,4)*P
     #FP(1,14))*ACG(1,7)+(ACG(9,4)*PFP(2,9)+ACG(10,4)*PFP(2,10)+ACG(14,4
     #)*PFP(2,14))*ACG(2,7)
      ST(56) = ST(56)+(ACG(7,4)*PFP(3,7)+ACG(8,4)*PFP(3,8)+ACG(13,4)*PFP
     #(3,13))*ACG(3,8)+(ACG(7,4)*PFP(4,7)+ACG(8,4)*PFP(4,8)+ACG(13,4)*PF
     #P(4,13))*ACG(4,8)
      ST(57) = ST(57)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,9)+(AC
     #G(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,9)
      ST(58) = ST(58)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,10)+(A
     #CG(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,10)
      ST(59) = ST(59)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,11)+(A
     #CG(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,11)
      ST(60) = ST(60)+(ACG(9,4)*PFP(1,9)+ACG(10,4)*PFP(1,10)+ACG(14,4)*P
     #FP(1,14))*ACG(1,12)+(ACG(9,4)*PFP(2,9)+ACG(10,4)*PFP(2,10)+ACG(14,
     #4)*PFP(2,14))*ACG(2,12)+(ACG(7,4)*PFP(3,7)+ACG(8,4)*PFP(3,8)+ACG(1
     #3,4)*PFP(3,13))*ACG(3,12)+(ACG(7,4)*PFP(4,7)+ACG(8,4)*PFP(4,8)+ACG
     #(13,4)*PFP(4,13))*ACG(4,12)+(ACG(13,4)*PFP(11,13)+ACG(14,4)*PFP(11
     #,14))*ACG(11,12)+(ACG(13,4)*PFP(12,13)+ACG(14,4)*PFP(12,14))*ACG(1
     #2,12)+(ACG(7,4)*PFP(7,15)+ACG(8,4)*PFP(8,15)+ACG(9,4)*PFP(9,15)+AC
     #G(10,4)*PFP(10,15)+ACG(13,4)*PFP(13,15)+ACG(14,4)*PFP(14,15))*ACG(
     #15,12)
      ST(61) = ST(61)+(ACG(9,4)*PFP(1,9)+ACG(10,4)*PFP(1,10)+ACG(14,4)*P
     #FP(1,14))*ACG(1,13)+(ACG(9,4)*PFP(2,9)+ACG(10,4)*PFP(2,10)+ACG(14,
     #4)*PFP(2,14))*ACG(2,13)
      ST(62) = ST(62)+(ACG(7,4)*PFP(3,7)+ACG(8,4)*PFP(3,8)+ACG(13,4)*PFP
     #(3,13))*ACG(3,14)+(ACG(7,4)*PFP(4,7)+ACG(8,4)*PFP(4,8)+ACG(13,4)*P
     #FP(4,13))*ACG(4,14)
      ST(63) = ST(63)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,15)+(A
     #CG(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,15)
      ST(64) = ST(64)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,16)+(A
     #CG(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,16)
      ST(65) = ST(65)+(ACG(5,4)*PFP(5,5)+ACG(6,4)*PFP(5,6))*ACG(5,17)+(A
     #CG(5,4)*PFP(5,6)+ACG(6,4)*PFP(6,6))*ACG(6,17)
      ST(66) = ST(66)+(ACG(9,4)*PFP(1,9)+ACG(10,4)*PFP(1,10)+ACG(14,4)*P
     #FP(1,14))*ACG(1,18)+(ACG(9,4)*PFP(2,9)+ACG(10,4)*PFP(2,10)+ACG(14,
     #4)*PFP(2,14))*ACG(2,18)+(ACG(7,4)*PFP(3,7)+ACG(8,4)*PFP(3,8)+ACG(1
     #3,4)*PFP(3,13))*ACG(3,18)+(ACG(7,4)*PFP(4,7)+ACG(8,4)*PFP(4,8)+ACG
     #(13,4)*PFP(4,13))*ACG(4,18)+(ACG(13,4)*PFP(11,13)+ACG(14,4)*PFP(11
     #,14))*ACG(11,18)+(ACG(13,4)*PFP(12,13)+ACG(14,4)*PFP(12,14))*ACG(1
     #2,18)+(ACG(7,4)*PFP(7,15)+ACG(8,4)*PFP(8,15)+ACG(9,4)*PFP(9,15)+AC
     #G(10,4)*PFP(10,15)+ACG(13,4)*PFP(13,15)+ACG(14,4)*PFP(14,15))*ACG(
     #15,18)
      ST(67) = ST(67)+(ACG(5,5)*PFP(5,5)+ACG(6,5)*PFP(5,6))*ACG(5,5)+(AC
     #G(5,5)*PFP(5,6)+ACG(6,5)*PFP(6,6))*ACG(6,5)
      ST(68) = ST(68)+(ACG(9,5)*PFP(1,9)+ACG(10,5)*PFP(1,10)+ACG(14,5)*P
     #FP(1,14))*ACG(1,6)+(ACG(9,5)*PFP(2,9)+ACG(10,5)*PFP(2,10)+ACG(14,5
     #)*PFP(2,14))*ACG(2,6)+(ACG(7,5)*PFP(3,7)+ACG(8,5)*PFP(3,8)+ACG(13,
     #5)*PFP(3,13))*ACG(3,6)+(ACG(7,5)*PFP(4,7)+ACG(8,5)*PFP(4,8)+ACG(13
     #,5)*PFP(4,13))*ACG(4,6)+(ACG(13,5)*PFP(11,13)+ACG(14,5)*PFP(11,14)
     #)*ACG(11,6)+(ACG(13,5)*PFP(12,13)+ACG(14,5)*PFP(12,14))*ACG(12,6)+
     #(ACG(7,5)*PFP(7,15)+ACG(8,5)*PFP(8,15)+ACG(9,5)*PFP(9,15)+ACG(10,5
     #)*PFP(10,15)+ACG(13,5)*PFP(13,15)+ACG(14,5)*PFP(14,15))*ACG(15,6)
      ST(69) = ST(69)+(ACG(9,5)*PFP(1,9)+ACG(10,5)*PFP(1,10)+ACG(14,5)*P
     #FP(1,14))*ACG(1,7)+(ACG(9,5)*PFP(2,9)+ACG(10,5)*PFP(2,10)+ACG(14,5
     #)*PFP(2,14))*ACG(2,7)
      ST(70) = ST(70)+(ACG(7,5)*PFP(3,7)+ACG(8,5)*PFP(3,8)+ACG(13,5)*PFP
     #(3,13))*ACG(3,8)+(ACG(7,5)*PFP(4,7)+ACG(8,5)*PFP(4,8)+ACG(13,5)*PF
     #P(4,13))*ACG(4,8)
      ST(71) = ST(71)+(ACG(5,5)*PFP(5,5)+ACG(6,5)*PFP(5,6))*ACG(5,9)+(AC
     #G(5,5)*PFP(5,6)+ACG(6,5)*PFP(6,6))*ACG(6,9)
      ST(72) = ST(72)+(ACG(5,5)*PFP(5,5)+ACG(6,5)*PFP(5,6))*ACG(5,10)+(A
     #CG(5,5)*PFP(5,6)+ACG(6,5)*PFP(6,6))*ACG(6,10)
      ST(73) = ST(73)+(ACG(5,5)*PFP(5,5)+ACG(6,5)*PFP(5,6))*ACG(5,11)+(A
     #CG(5,5)*PFP(5,6)+ACG(6,5)*PFP(6,6))*ACG(6,11)
      ST(74) = ST(74)+(ACG(9,5)*PFP(1,9)+ACG(10,5)*PFP(1,10)+ACG(14,5)*P
     #FP(1,14))*ACG(1,12)+(ACG(9,5)*PFP(2,9)+ACG(10,5)*PFP(2,10)+ACG(14,
     #5)*PFP(2,14))*ACG(2,12)+(ACG(7,5)*PFP(3,7)+ACG(8,5)*PFP(3,8)+ACG(1
     #3,5)*PFP(3,13))*ACG(3,12)+(ACG(7,5)*PFP(4,7)+ACG(8,5)*PFP(4,8)+ACG
     #(13,5)*PFP(4,13))*ACG(4,12)+(ACG(13,5)*PFP(11,13)+ACG(14,5)*PFP(11
     #,14))*ACG(11,12)+(ACG(13,5)*PFP(12,13)+ACG(14,5)*PFP(12,14))*ACG(1
     #2,12)+(ACG(7,5)*PFP(7,15)+ACG(8,5)*PFP(8,15)+ACG(9,5)*PFP(9,15)+AC
     #G(10,5)*PFP(10,15)+ACG(13,5)*PFP(13,15)+ACG(14,5)*PFP(14,15))*ACG(
     #15,12)
      ST(75) = ST(75)+(ACG(9,5)*PFP(1,9)+ACG(10,5)*PFP(1,10)+ACG(14,5)*P
     #FP(1,14))*ACG(1,13)+(ACG(9,5)*PFP(2,9)+ACG(10,5)*PFP(2,10)+ACG(14,
     #5)*PFP(2,14))*ACG(2,13)
      ST(76) = ST(76)+(ACG(7,5)*PFP(3,7)+ACG(8,5)*PFP(3,8)+ACG(13,5)*PFP
     #(3,13))*ACG(3,14)+(ACG(7,5)*PFP(4,7)+ACG(8,5)*PFP(4,8)+ACG(13,5)*P
     #FP(4,13))*ACG(4,14)
      ST(77) = ST(77)+(ACG(5,5)*PFP(5,5)+ACG(6,5)*PFP(5,6))*ACG(5,15)+(A
     #CG(5,5)*PFP(5,6)+ACG(6,5)*PFP(6,6))*ACG(6,15)
      ST(78) = ST(78)+(ACG(5,5)*PFP(5,5)+ACG(6,5)*PFP(5,6))*ACG(5,16)+(A
     #CG(5,5)*PFP(5,6)+ACG(6,5)*PFP(6,6))*ACG(6,16)
      ST(79) = ST(79)+(ACG(5,5)*PFP(5,5)+ACG(6,5)*PFP(5,6))*ACG(5,17)+(A
     #CG(5,5)*PFP(5,6)+ACG(6,5)*PFP(6,6))*ACG(6,17)
      ST(80) = ST(80)+(ACG(9,5)*PFP(1,9)+ACG(10,5)*PFP(1,10)+ACG(14,5)*P
     #FP(1,14))*ACG(1,18)+(ACG(9,5)*PFP(2,9)+ACG(10,5)*PFP(2,10)+ACG(14,
     #5)*PFP(2,14))*ACG(2,18)+(ACG(7,5)*PFP(3,7)+ACG(8,5)*PFP(3,8)+ACG(1
     #3,5)*PFP(3,13))*ACG(3,18)+(ACG(7,5)*PFP(4,7)+ACG(8,5)*PFP(4,8)+ACG
     #(13,5)*PFP(4,13))*ACG(4,18)+(ACG(13,5)*PFP(11,13)+ACG(14,5)*PFP(11
     #,14))*ACG(11,18)+(ACG(13,5)*PFP(12,13)+ACG(14,5)*PFP(12,14))*ACG(1
     #2,18)+(ACG(7,5)*PFP(7,15)+ACG(8,5)*PFP(8,15)+ACG(9,5)*PFP(9,15)+AC
     #G(10,5)*PFP(10,15)+ACG(13,5)*PFP(13,15)+ACG(14,5)*PFP(14,15))*ACG(
     #15,18)
      ST(81) = ST(81)+(ACG(1,6)*PFP(1,1)+ACG(2,6)*PFP(1,2))*ACG(1,6)+(AC
     #G(1,6)*PFP(1,2)+ACG(2,6)*PFP(2,2))*ACG(2,6)+(ACG(3,6)*PFP(3,3)+ACG
     #(4,6)*PFP(3,4))*ACG(3,6)+(ACG(3,6)*PFP(3,4)+ACG(4,6)*PFP(4,4))*ACG
     #(4,6)
      ST(82) = ST(82)+(ACG(1,6)*PFP(1,1)+ACG(2,6)*PFP(1,2))*ACG(1,7)+(AC
     #G(1,6)*PFP(1,2)+ACG(2,6)*PFP(2,2))*ACG(2,7)
      ST(83) = ST(83)+(ACG(3,6)*PFP(3,3)+ACG(4,6)*PFP(3,4))*ACG(3,8)+(AC
     #G(3,6)*PFP(3,4)+ACG(4,6)*PFP(4,4))*ACG(4,8)
      ST(84) = ST(84)+(ACG(3,6)*PFP(3,7)+ACG(4,6)*PFP(4,7)+ACG(15,6)*PFP
     #(7,15))*ACG(7,9)+(ACG(3,6)*PFP(3,8)+ACG(4,6)*PFP(4,8)+ACG(15,6)*PF
     #P(8,15))*ACG(8,9)+(ACG(1,6)*PFP(1,9)+ACG(2,6)*PFP(2,9)+ACG(15,6)*P
     #FP(9,15))*ACG(9,9)+(ACG(1,6)*PFP(1,10)+ACG(2,6)*PFP(2,10)+ACG(15,6
     #)*PFP(10,15))*ACG(10,9)+(ACG(3,6)*PFP(3,13)+ACG(4,6)*PFP(4,13)+ACG
     #(11,6)*PFP(11,13)+ACG(12,6)*PFP(12,13)+ACG(15,6)*PFP(13,15))*ACG(1
     #3,9)+(ACG(1,6)*PFP(1,14)+ACG(2,6)*PFP(2,14)+ACG(11,6)*PFP(11,14)+A
     #CG(12,6)*PFP(12,14)+ACG(15,6)*PFP(14,15))*ACG(14,9)
      ST(85) = ST(85)+(ACG(3,6)*PFP(3,7)+ACG(4,6)*PFP(4,7)+ACG(15,6)*PFP
     #(7,15))*ACG(7,10)+(ACG(3,6)*PFP(3,8)+ACG(4,6)*PFP(4,8)+ACG(15,6)*P
     #FP(8,15))*ACG(8,10)+(ACG(1,6)*PFP(1,9)+ACG(2,6)*PFP(2,9)+ACG(15,6)
     #*PFP(9,15))*ACG(9,10)+(ACG(1,6)*PFP(1,10)+ACG(2,6)*PFP(2,10)+ACG(1
     #5,6)*PFP(10,15))*ACG(10,10)+(ACG(3,6)*PFP(3,13)+ACG(4,6)*PFP(4,13)
     #+ACG(11,6)*PFP(11,13)+ACG(12,6)*PFP(12,13)+ACG(15,6)*PFP(13,15))*A
     #CG(13,10)+(ACG(1,6)*PFP(1,14)+ACG(2,6)*PFP(2,14)+ACG(11,6)*PFP(11,
     #14)+ACG(12,6)*PFP(12,14)+ACG(15,6)*PFP(14,15))*ACG(14,10)
      ST(86) = ST(86)+(ACG(3,6)*PFP(3,7)+ACG(4,6)*PFP(4,7)+ACG(15,6)*PFP
     #(7,15))*ACG(7,11)+(ACG(3,6)*PFP(3,8)+ACG(4,6)*PFP(4,8)+ACG(15,6)*P
     #FP(8,15))*ACG(8,11)+(ACG(1,6)*PFP(1,9)+ACG(2,6)*PFP(2,9)+ACG(15,6)
     #*PFP(9,15))*ACG(9,11)+(ACG(1,6)*PFP(1,10)+ACG(2,6)*PFP(2,10)+ACG(1
     #5,6)*PFP(10,15))*ACG(10,11)+(ACG(3,6)*PFP(3,13)+ACG(4,6)*PFP(4,13)
     #+ACG(11,6)*PFP(11,13)+ACG(12,6)*PFP(12,13)+ACG(15,6)*PFP(13,15))*A
     #CG(13,11)+(ACG(1,6)*PFP(1,14)+ACG(2,6)*PFP(2,14)+ACG(11,6)*PFP(11,
     #14)+ACG(12,6)*PFP(12,14)+ACG(15,6)*PFP(14,15))*ACG(14,11)
      ST(87) = ST(87)+(ACG(1,6)*PFP(1,1)+ACG(2,6)*PFP(1,2))*ACG(1,12)+(A
     #CG(1,6)*PFP(1,2)+ACG(2,6)*PFP(2,2))*ACG(2,12)+(ACG(3,6)*PFP(3,3)+A
     #CG(4,6)*PFP(3,4))*ACG(3,12)+(ACG(3,6)*PFP(3,4)+ACG(4,6)*PFP(4,4))*
     #ACG(4,12)
      ST(88) = ST(88)+(ACG(1,6)*PFP(1,1)+ACG(2,6)*PFP(1,2))*ACG(1,13)+(A
     #CG(1,6)*PFP(1,2)+ACG(2,6)*PFP(2,2))*ACG(2,13)
      ST(89) = ST(89)+(ACG(3,6)*PFP(3,3)+ACG(4,6)*PFP(3,4))*ACG(3,14)+(A
     #CG(3,6)*PFP(3,4)+ACG(4,6)*PFP(4,4))*ACG(4,14)
      ST(90) = ST(90)+(ACG(3,6)*PFP(3,7)+ACG(4,6)*PFP(4,7)+ACG(15,6)*PFP
     #(7,15))*ACG(7,15)+(ACG(3,6)*PFP(3,8)+ACG(4,6)*PFP(4,8)+ACG(15,6)*P
     #FP(8,15))*ACG(8,15)+(ACG(1,6)*PFP(1,9)+ACG(2,6)*PFP(2,9)+ACG(15,6)
     #*PFP(9,15))*ACG(9,15)+(ACG(1,6)*PFP(1,10)+ACG(2,6)*PFP(2,10)+ACG(1
     #5,6)*PFP(10,15))*ACG(10,15)+(ACG(3,6)*PFP(3,13)+ACG(4,6)*PFP(4,13)
     #+ACG(11,6)*PFP(11,13)+ACG(12,6)*PFP(12,13)+ACG(15,6)*PFP(13,15))*A
     #CG(13,15)+(ACG(1,6)*PFP(1,14)+ACG(2,6)*PFP(2,14)+ACG(11,6)*PFP(11,
     #14)+ACG(12,6)*PFP(12,14)+ACG(15,6)*PFP(14,15))*ACG(14,15)
      ST(91) = ST(91)+(ACG(3,6)*PFP(3,7)+ACG(4,6)*PFP(4,7)+ACG(15,6)*PFP
     #(7,15))*ACG(7,16)+(ACG(3,6)*PFP(3,8)+ACG(4,6)*PFP(4,8)+ACG(15,6)*P
     #FP(8,15))*ACG(8,16)+(ACG(1,6)*PFP(1,9)+ACG(2,6)*PFP(2,9)+ACG(15,6)
     #*PFP(9,15))*ACG(9,16)+(ACG(1,6)*PFP(1,10)+ACG(2,6)*PFP(2,10)+ACG(1
     #5,6)*PFP(10,15))*ACG(10,16)+(ACG(3,6)*PFP(3,13)+ACG(4,6)*PFP(4,13)
     #+ACG(11,6)*PFP(11,13)+ACG(12,6)*PFP(12,13)+ACG(15,6)*PFP(13,15))*A
     #CG(13,16)+(ACG(1,6)*PFP(1,14)+ACG(2,6)*PFP(2,14)+ACG(11,6)*PFP(11,
     #14)+ACG(12,6)*PFP(12,14)+ACG(15,6)*PFP(14,15))*ACG(14,16)
      ST(92) = ST(92)+(ACG(3,6)*PFP(3,7)+ACG(4,6)*PFP(4,7)+ACG(15,6)*PFP
     #(7,15))*ACG(7,17)+(ACG(3,6)*PFP(3,8)+ACG(4,6)*PFP(4,8)+ACG(15,6)*P
     #FP(8,15))*ACG(8,17)+(ACG(1,6)*PFP(1,9)+ACG(2,6)*PFP(2,9)+ACG(15,6)
     #*PFP(9,15))*ACG(9,17)+(ACG(1,6)*PFP(1,10)+ACG(2,6)*PFP(2,10)+ACG(1
     #5,6)*PFP(10,15))*ACG(10,17)+(ACG(3,6)*PFP(3,13)+ACG(4,6)*PFP(4,13)
     #+ACG(11,6)*PFP(11,13)+ACG(12,6)*PFP(12,13)+ACG(15,6)*PFP(13,15))*A
     #CG(13,17)+(ACG(1,6)*PFP(1,14)+ACG(2,6)*PFP(2,14)+ACG(11,6)*PFP(11,
     #14)+ACG(12,6)*PFP(12,14)+ACG(15,6)*PFP(14,15))*ACG(14,17)
      ST(93) = ST(93)+(ACG(1,6)*PFP(1,1)+ACG(2,6)*PFP(1,2))*ACG(1,18)+(A
     #CG(1,6)*PFP(1,2)+ACG(2,6)*PFP(2,2))*ACG(2,18)+(ACG(3,6)*PFP(3,3)+A
     #CG(4,6)*PFP(3,4))*ACG(3,18)+(ACG(3,6)*PFP(3,4)+ACG(4,6)*PFP(4,4))*
     #ACG(4,18)
      ST(94) = ST(94)+(ACG(1,7)*PFP(1,1)+ACG(2,7)*PFP(1,2))*ACG(1,7)+(AC
     #G(1,7)*PFP(1,2)+ACG(2,7)*PFP(2,2))*ACG(2,7)
      ST(95) = ST(95)
      ST(96) = ST(96)+(ACG(1,7)*PFP(1,9)+ACG(2,7)*PFP(2,9))*ACG(9,9)+(AC
     #G(1,7)*PFP(1,10)+ACG(2,7)*PFP(2,10))*ACG(10,9)+(ACG(1,7)*PFP(1,14)
     #+ACG(2,7)*PFP(2,14))*ACG(14,9)
      ST(97) = ST(97)+(ACG(1,7)*PFP(1,9)+ACG(2,7)*PFP(2,9))*ACG(9,10)+(A
     #CG(1,7)*PFP(1,10)+ACG(2,7)*PFP(2,10))*ACG(10,10)+(ACG(1,7)*PFP(1,1
     #4)+ACG(2,7)*PFP(2,14))*ACG(14,10)
      ST(98) = ST(98)+(ACG(1,7)*PFP(1,9)+ACG(2,7)*PFP(2,9))*ACG(9,11)+(A
     #CG(1,7)*PFP(1,10)+ACG(2,7)*PFP(2,10))*ACG(10,11)+(ACG(1,7)*PFP(1,1
     #4)+ACG(2,7)*PFP(2,14))*ACG(14,11)
      ST(99) = ST(99)+(ACG(1,7)*PFP(1,1)+ACG(2,7)*PFP(1,2))*ACG(1,12)+(A
     #CG(1,7)*PFP(1,2)+ACG(2,7)*PFP(2,2))*ACG(2,12)
      ST(100) = ST(100)+(ACG(1,7)*PFP(1,1)+ACG(2,7)*PFP(1,2))*ACG(1,13)+
     #(ACG(1,7)*PFP(1,2)+ACG(2,7)*PFP(2,2))*ACG(2,13)
      ST(101) = ST(101)
      ST(102) = ST(102)+(ACG(1,7)*PFP(1,9)+ACG(2,7)*PFP(2,9))*ACG(9,15)+
     #(ACG(1,7)*PFP(1,10)+ACG(2,7)*PFP(2,10))*ACG(10,15)+(ACG(1,7)*PFP(1
     #,14)+ACG(2,7)*PFP(2,14))*ACG(14,15)
      ST(103) = ST(103)+(ACG(1,7)*PFP(1,9)+ACG(2,7)*PFP(2,9))*ACG(9,16)+
     #(ACG(1,7)*PFP(1,10)+ACG(2,7)*PFP(2,10))*ACG(10,16)+(ACG(1,7)*PFP(1
     #,14)+ACG(2,7)*PFP(2,14))*ACG(14,16)
      ST(104) = ST(104)+(ACG(1,7)*PFP(1,9)+ACG(2,7)*PFP(2,9))*ACG(9,17)+
     #(ACG(1,7)*PFP(1,10)+ACG(2,7)*PFP(2,10))*ACG(10,17)+(ACG(1,7)*PFP(1
     #,14)+ACG(2,7)*PFP(2,14))*ACG(14,17)
      ST(105) = ST(105)+(ACG(1,7)*PFP(1,1)+ACG(2,7)*PFP(1,2))*ACG(1,18)+
     #(ACG(1,7)*PFP(1,2)+ACG(2,7)*PFP(2,2))*ACG(2,18)
      ST(106) = ST(106)+(ACG(3,8)*PFP(3,3)+ACG(4,8)*PFP(3,4))*ACG(3,8)+(
     #ACG(3,8)*PFP(3,4)+ACG(4,8)*PFP(4,4))*ACG(4,8)
      ST(107) = ST(107)+(ACG(3,8)*PFP(3,7)+ACG(4,8)*PFP(4,7))*ACG(7,9)+(
     #ACG(3,8)*PFP(3,8)+ACG(4,8)*PFP(4,8))*ACG(8,9)+(ACG(3,8)*PFP(3,13)+
     #ACG(4,8)*PFP(4,13))*ACG(13,9)
      ST(108) = ST(108)+(ACG(3,8)*PFP(3,7)+ACG(4,8)*PFP(4,7))*ACG(7,10)+
     #(ACG(3,8)*PFP(3,8)+ACG(4,8)*PFP(4,8))*ACG(8,10)+(ACG(3,8)*PFP(3,13
     #)+ACG(4,8)*PFP(4,13))*ACG(13,10)
      ST(109) = ST(109)+(ACG(3,8)*PFP(3,7)+ACG(4,8)*PFP(4,7))*ACG(7,11)+
     #(ACG(3,8)*PFP(3,8)+ACG(4,8)*PFP(4,8))*ACG(8,11)+(ACG(3,8)*PFP(3,13
     #)+ACG(4,8)*PFP(4,13))*ACG(13,11)
      ST(110) = ST(110)+(ACG(3,8)*PFP(3,3)+ACG(4,8)*PFP(3,4))*ACG(3,12)+
     #(ACG(3,8)*PFP(3,4)+ACG(4,8)*PFP(4,4))*ACG(4,12)
      ST(111) = ST(111)
      ST(112) = ST(112)+(ACG(3,8)*PFP(3,3)+ACG(4,8)*PFP(3,4))*ACG(3,14)+
     #(ACG(3,8)*PFP(3,4)+ACG(4,8)*PFP(4,4))*ACG(4,14)
      ST(113) = ST(113)+(ACG(3,8)*PFP(3,7)+ACG(4,8)*PFP(4,7))*ACG(7,15)+
     #(ACG(3,8)*PFP(3,8)+ACG(4,8)*PFP(4,8))*ACG(8,15)+(ACG(3,8)*PFP(3,13
     #)+ACG(4,8)*PFP(4,13))*ACG(13,15)
      ST(114) = ST(114)+(ACG(3,8)*PFP(3,7)+ACG(4,8)*PFP(4,7))*ACG(7,16)+
     #(ACG(3,8)*PFP(3,8)+ACG(4,8)*PFP(4,8))*ACG(8,16)+(ACG(3,8)*PFP(3,13
     #)+ACG(4,8)*PFP(4,13))*ACG(13,16)
      ST(115) = ST(115)+(ACG(3,8)*PFP(3,7)+ACG(4,8)*PFP(4,7))*ACG(7,17)+
     #(ACG(3,8)*PFP(3,8)+ACG(4,8)*PFP(4,8))*ACG(8,17)+(ACG(3,8)*PFP(3,13
     #)+ACG(4,8)*PFP(4,13))*ACG(13,17)
      ST(116) = ST(116)+(ACG(3,8)*PFP(3,3)+ACG(4,8)*PFP(3,4))*ACG(3,18)+
     #(ACG(3,8)*PFP(3,4)+ACG(4,8)*PFP(4,4))*ACG(4,18)
      ST(117) = ST(117)+(ACG(5,9)*PFP(5,5)+ACG(6,9)*PFP(5,6))*ACG(5,9)+(
     #ACG(5,9)*PFP(5,6)+ACG(6,9)*PFP(6,6))*ACG(6,9)
      ST(118) = ST(118)+(ACG(5,9)*PFP(5,5)+ACG(6,9)*PFP(5,6))*ACG(5,10)+
     #(ACG(5,9)*PFP(5,6)+ACG(6,9)*PFP(6,6))*ACG(6,10)
      ST(119) = ST(119)+(ACG(5,9)*PFP(5,5)+ACG(6,9)*PFP(5,6))*ACG(5,11)+
     #(ACG(5,9)*PFP(5,6)+ACG(6,9)*PFP(6,6))*ACG(6,11)
      ST(120) = ST(120)+(ACG(9,9)*PFP(1,9)+ACG(10,9)*PFP(1,10)+ACG(14,9)
     #*PFP(1,14))*ACG(1,12)+(ACG(9,9)*PFP(2,9)+ACG(10,9)*PFP(2,10)+ACG(1
     #4,9)*PFP(2,14))*ACG(2,12)+(ACG(7,9)*PFP(3,7)+ACG(8,9)*PFP(3,8)+ACG
     #(13,9)*PFP(3,13))*ACG(3,12)+(ACG(7,9)*PFP(4,7)+ACG(8,9)*PFP(4,8)+A
     #CG(13,9)*PFP(4,13))*ACG(4,12)+(ACG(13,9)*PFP(11,13)+ACG(14,9)*PFP(
     #11,14))*ACG(11,12)+(ACG(13,9)*PFP(12,13)+ACG(14,9)*PFP(12,14))*ACG
     #(12,12)+(ACG(7,9)*PFP(7,15)+ACG(8,9)*PFP(8,15)+ACG(9,9)*PFP(9,15)+
     #ACG(10,9)*PFP(10,15)+ACG(13,9)*PFP(13,15)+ACG(14,9)*PFP(14,15))*AC
     #G(15,12)
      ST(121) = ST(121)+(ACG(9,9)*PFP(1,9)+ACG(10,9)*PFP(1,10)+ACG(14,9)
     #*PFP(1,14))*ACG(1,13)+(ACG(9,9)*PFP(2,9)+ACG(10,9)*PFP(2,10)+ACG(1
     #4,9)*PFP(2,14))*ACG(2,13)
      ST(122) = ST(122)+(ACG(7,9)*PFP(3,7)+ACG(8,9)*PFP(3,8)+ACG(13,9)*P
     #FP(3,13))*ACG(3,14)+(ACG(7,9)*PFP(4,7)+ACG(8,9)*PFP(4,8)+ACG(13,9)
     #*PFP(4,13))*ACG(4,14)
      ST(123) = ST(123)+(ACG(5,9)*PFP(5,5)+ACG(6,9)*PFP(5,6))*ACG(5,15)+
     #(ACG(5,9)*PFP(5,6)+ACG(6,9)*PFP(6,6))*ACG(6,15)
      ST(124) = ST(124)+(ACG(5,9)*PFP(5,5)+ACG(6,9)*PFP(5,6))*ACG(5,16)+
     #(ACG(5,9)*PFP(5,6)+ACG(6,9)*PFP(6,6))*ACG(6,16)
      ST(125) = ST(125)+(ACG(5,9)*PFP(5,5)+ACG(6,9)*PFP(5,6))*ACG(5,17)+
     #(ACG(5,9)*PFP(5,6)+ACG(6,9)*PFP(6,6))*ACG(6,17)
      ST(126) = ST(126)+(ACG(9,9)*PFP(1,9)+ACG(10,9)*PFP(1,10)+ACG(14,9)
     #*PFP(1,14))*ACG(1,18)+(ACG(9,9)*PFP(2,9)+ACG(10,9)*PFP(2,10)+ACG(1
     #4,9)*PFP(2,14))*ACG(2,18)+(ACG(7,9)*PFP(3,7)+ACG(8,9)*PFP(3,8)+ACG
     #(13,9)*PFP(3,13))*ACG(3,18)+(ACG(7,9)*PFP(4,7)+ACG(8,9)*PFP(4,8)+A
     #CG(13,9)*PFP(4,13))*ACG(4,18)+(ACG(13,9)*PFP(11,13)+ACG(14,9)*PFP(
     #11,14))*ACG(11,18)+(ACG(13,9)*PFP(12,13)+ACG(14,9)*PFP(12,14))*ACG
     #(12,18)+(ACG(7,9)*PFP(7,15)+ACG(8,9)*PFP(8,15)+ACG(9,9)*PFP(9,15)+
     #ACG(10,9)*PFP(10,15)+ACG(13,9)*PFP(13,15)+ACG(14,9)*PFP(14,15))*AC
     #G(15,18)
      ST(127) = ST(127)+(ACG(5,10)*PFP(5,5)+ACG(6,10)*PFP(5,6))*ACG(5,10
     #)+(ACG(5,10)*PFP(5,6)+ACG(6,10)*PFP(6,6))*ACG(6,10)
      ST(128) = ST(128)+(ACG(5,10)*PFP(5,5)+ACG(6,10)*PFP(5,6))*ACG(5,11
     #)+(ACG(5,10)*PFP(5,6)+ACG(6,10)*PFP(6,6))*ACG(6,11)
      ST(129) = ST(129)+(ACG(9,10)*PFP(1,9)+ACG(10,10)*PFP(1,10)+ACG(14,
     #10)*PFP(1,14))*ACG(1,12)+(ACG(9,10)*PFP(2,9)+ACG(10,10)*PFP(2,10)+
     #ACG(14,10)*PFP(2,14))*ACG(2,12)+(ACG(7,10)*PFP(3,7)+ACG(8,10)*PFP(
     #3,8)+ACG(13,10)*PFP(3,13))*ACG(3,12)+(ACG(7,10)*PFP(4,7)+ACG(8,10)
     #*PFP(4,8)+ACG(13,10)*PFP(4,13))*ACG(4,12)+(ACG(13,10)*PFP(11,13)+A
     #CG(14,10)*PFP(11,14))*ACG(11,12)+(ACG(13,10)*PFP(12,13)+ACG(14,10)
     #*PFP(12,14))*ACG(12,12)+(ACG(7,10)*PFP(7,15)+ACG(8,10)*PFP(8,15)+A
     #CG(9,10)*PFP(9,15)+ACG(10,10)*PFP(10,15)+ACG(13,10)*PFP(13,15)+ACG
     #(14,10)*PFP(14,15))*ACG(15,12)
      ST(130) = ST(130)+(ACG(9,10)*PFP(1,9)+ACG(10,10)*PFP(1,10)+ACG(14,
     #10)*PFP(1,14))*ACG(1,13)+(ACG(9,10)*PFP(2,9)+ACG(10,10)*PFP(2,10)+
     #ACG(14,10)*PFP(2,14))*ACG(2,13)
      ST(131) = ST(131)+(ACG(7,10)*PFP(3,7)+ACG(8,10)*PFP(3,8)+ACG(13,10
     #)*PFP(3,13))*ACG(3,14)+(ACG(7,10)*PFP(4,7)+ACG(8,10)*PFP(4,8)+ACG(
     #13,10)*PFP(4,13))*ACG(4,14)
      ST(132) = ST(132)+(ACG(5,10)*PFP(5,5)+ACG(6,10)*PFP(5,6))*ACG(5,15
     #)+(ACG(5,10)*PFP(5,6)+ACG(6,10)*PFP(6,6))*ACG(6,15)
      ST(133) = ST(133)+(ACG(5,10)*PFP(5,5)+ACG(6,10)*PFP(5,6))*ACG(5,16
     #)+(ACG(5,10)*PFP(5,6)+ACG(6,10)*PFP(6,6))*ACG(6,16)
      ST(134) = ST(134)+(ACG(5,10)*PFP(5,5)+ACG(6,10)*PFP(5,6))*ACG(5,17
     #)+(ACG(5,10)*PFP(5,6)+ACG(6,10)*PFP(6,6))*ACG(6,17)
      ST(135) = ST(135)+(ACG(9,10)*PFP(1,9)+ACG(10,10)*PFP(1,10)+ACG(14,
     #10)*PFP(1,14))*ACG(1,18)+(ACG(9,10)*PFP(2,9)+ACG(10,10)*PFP(2,10)+
     #ACG(14,10)*PFP(2,14))*ACG(2,18)+(ACG(7,10)*PFP(3,7)+ACG(8,10)*PFP(
     #3,8)+ACG(13,10)*PFP(3,13))*ACG(3,18)+(ACG(7,10)*PFP(4,7)+ACG(8,10)
     #*PFP(4,8)+ACG(13,10)*PFP(4,13))*ACG(4,18)+(ACG(13,10)*PFP(11,13)+A
     #CG(14,10)*PFP(11,14))*ACG(11,18)+(ACG(13,10)*PFP(12,13)+ACG(14,10)
     #*PFP(12,14))*ACG(12,18)+(ACG(7,10)*PFP(7,15)+ACG(8,10)*PFP(8,15)+A
     #CG(9,10)*PFP(9,15)+ACG(10,10)*PFP(10,15)+ACG(13,10)*PFP(13,15)+ACG
     #(14,10)*PFP(14,15))*ACG(15,18)
      ST(136) = ST(136)+(ACG(5,11)*PFP(5,5)+ACG(6,11)*PFP(5,6))*ACG(5,11
     #)+(ACG(5,11)*PFP(5,6)+ACG(6,11)*PFP(6,6))*ACG(6,11)
      ST(137) = ST(137)+(ACG(9,11)*PFP(1,9)+ACG(10,11)*PFP(1,10)+ACG(14,
     #11)*PFP(1,14))*ACG(1,12)+(ACG(9,11)*PFP(2,9)+ACG(10,11)*PFP(2,10)+
     #ACG(14,11)*PFP(2,14))*ACG(2,12)+(ACG(7,11)*PFP(3,7)+ACG(8,11)*PFP(
     #3,8)+ACG(13,11)*PFP(3,13))*ACG(3,12)+(ACG(7,11)*PFP(4,7)+ACG(8,11)
     #*PFP(4,8)+ACG(13,11)*PFP(4,13))*ACG(4,12)+(ACG(13,11)*PFP(11,13)+A
     #CG(14,11)*PFP(11,14))*ACG(11,12)+(ACG(13,11)*PFP(12,13)+ACG(14,11)
     #*PFP(12,14))*ACG(12,12)+(ACG(7,11)*PFP(7,15)+ACG(8,11)*PFP(8,15)+A
     #CG(9,11)*PFP(9,15)+ACG(10,11)*PFP(10,15)+ACG(13,11)*PFP(13,15)+ACG
     #(14,11)*PFP(14,15))*ACG(15,12)
      ST(138) = ST(138)+(ACG(9,11)*PFP(1,9)+ACG(10,11)*PFP(1,10)+ACG(14,
     #11)*PFP(1,14))*ACG(1,13)+(ACG(9,11)*PFP(2,9)+ACG(10,11)*PFP(2,10)+
     #ACG(14,11)*PFP(2,14))*ACG(2,13)
      ST(139) = ST(139)+(ACG(7,11)*PFP(3,7)+ACG(8,11)*PFP(3,8)+ACG(13,11
     #)*PFP(3,13))*ACG(3,14)+(ACG(7,11)*PFP(4,7)+ACG(8,11)*PFP(4,8)+ACG(
     #13,11)*PFP(4,13))*ACG(4,14)
      ST(140) = ST(140)+(ACG(5,11)*PFP(5,5)+ACG(6,11)*PFP(5,6))*ACG(5,15
     #)+(ACG(5,11)*PFP(5,6)+ACG(6,11)*PFP(6,6))*ACG(6,15)
      ST(141) = ST(141)+(ACG(5,11)*PFP(5,5)+ACG(6,11)*PFP(5,6))*ACG(5,16
     #)+(ACG(5,11)*PFP(5,6)+ACG(6,11)*PFP(6,6))*ACG(6,16)
      ST(142) = ST(142)+(ACG(5,11)*PFP(5,5)+ACG(6,11)*PFP(5,6))*ACG(5,17
     #)+(ACG(5,11)*PFP(5,6)+ACG(6,11)*PFP(6,6))*ACG(6,17)
      ST(143) = ST(143)+(ACG(9,11)*PFP(1,9)+ACG(10,11)*PFP(1,10)+ACG(14,
     #11)*PFP(1,14))*ACG(1,18)+(ACG(9,11)*PFP(2,9)+ACG(10,11)*PFP(2,10)+
     #ACG(14,11)*PFP(2,14))*ACG(2,18)+(ACG(7,11)*PFP(3,7)+ACG(8,11)*PFP(
     #3,8)+ACG(13,11)*PFP(3,13))*ACG(3,18)+(ACG(7,11)*PFP(4,7)+ACG(8,11)
     #*PFP(4,8)+ACG(13,11)*PFP(4,13))*ACG(4,18)+(ACG(13,11)*PFP(11,13)+A
     #CG(14,11)*PFP(11,14))*ACG(11,18)+(ACG(13,11)*PFP(12,13)+ACG(14,11)
     #*PFP(12,14))*ACG(12,18)+(ACG(7,11)*PFP(7,15)+ACG(8,11)*PFP(8,15)+A
     #CG(9,11)*PFP(9,15)+ACG(10,11)*PFP(10,15)+ACG(13,11)*PFP(13,15)+ACG
     #(14,11)*PFP(14,15))*ACG(15,18)
      ST(144) = ST(144)+(ACG(1,12)*PFP(1,1)+ACG(2,12)*PFP(1,2))*ACG(1,12
     #)+(ACG(1,12)*PFP(1,2)+ACG(2,12)*PFP(2,2))*ACG(2,12)+(ACG(3,12)*PFP
     #(3,3)+ACG(4,12)*PFP(3,4))*ACG(3,12)+(ACG(3,12)*PFP(3,4)+ACG(4,12)*
     #PFP(4,4))*ACG(4,12)
      ST(145) = ST(145)+(ACG(1,12)*PFP(1,1)+ACG(2,12)*PFP(1,2))*ACG(1,13
     #)+(ACG(1,12)*PFP(1,2)+ACG(2,12)*PFP(2,2))*ACG(2,13)
      ST(146) = ST(146)+(ACG(3,12)*PFP(3,3)+ACG(4,12)*PFP(3,4))*ACG(3,14
     #)+(ACG(3,12)*PFP(3,4)+ACG(4,12)*PFP(4,4))*ACG(4,14)
      ST(147) = ST(147)+(ACG(3,12)*PFP(3,7)+ACG(4,12)*PFP(4,7)+ACG(15,12
     #)*PFP(7,15))*ACG(7,15)+(ACG(3,12)*PFP(3,8)+ACG(4,12)*PFP(4,8)+ACG(
     #15,12)*PFP(8,15))*ACG(8,15)+(ACG(1,12)*PFP(1,9)+ACG(2,12)*PFP(2,9)
     #+ACG(15,12)*PFP(9,15))*ACG(9,15)+(ACG(1,12)*PFP(1,10)+ACG(2,12)*PF
     #P(2,10)+ACG(15,12)*PFP(10,15))*ACG(10,15)+(ACG(3,12)*PFP(3,13)+ACG
     #(4,12)*PFP(4,13)+ACG(11,12)*PFP(11,13)+ACG(12,12)*PFP(12,13)+ACG(1
     #5,12)*PFP(13,15))*ACG(13,15)+(ACG(1,12)*PFP(1,14)+ACG(2,12)*PFP(2,
     #14)+ACG(11,12)*PFP(11,14)+ACG(12,12)*PFP(12,14)+ACG(15,12)*PFP(14,
     #15))*ACG(14,15)
      ST(148) = ST(148)+(ACG(3,12)*PFP(3,7)+ACG(4,12)*PFP(4,7)+ACG(15,12
     #)*PFP(7,15))*ACG(7,16)+(ACG(3,12)*PFP(3,8)+ACG(4,12)*PFP(4,8)+ACG(
     #15,12)*PFP(8,15))*ACG(8,16)+(ACG(1,12)*PFP(1,9)+ACG(2,12)*PFP(2,9)
     #+ACG(15,12)*PFP(9,15))*ACG(9,16)+(ACG(1,12)*PFP(1,10)+ACG(2,12)*PF
     #P(2,10)+ACG(15,12)*PFP(10,15))*ACG(10,16)+(ACG(3,12)*PFP(3,13)+ACG
     #(4,12)*PFP(4,13)+ACG(11,12)*PFP(11,13)+ACG(12,12)*PFP(12,13)+ACG(1
     #5,12)*PFP(13,15))*ACG(13,16)+(ACG(1,12)*PFP(1,14)+ACG(2,12)*PFP(2,
     #14)+ACG(11,12)*PFP(11,14)+ACG(12,12)*PFP(12,14)+ACG(15,12)*PFP(14,
     #15))*ACG(14,16)
      ST(149) = ST(149)+(ACG(3,12)*PFP(3,7)+ACG(4,12)*PFP(4,7)+ACG(15,12
     #)*PFP(7,15))*ACG(7,17)+(ACG(3,12)*PFP(3,8)+ACG(4,12)*PFP(4,8)+ACG(
     #15,12)*PFP(8,15))*ACG(8,17)+(ACG(1,12)*PFP(1,9)+ACG(2,12)*PFP(2,9)
     #+ACG(15,12)*PFP(9,15))*ACG(9,17)+(ACG(1,12)*PFP(1,10)+ACG(2,12)*PF
     #P(2,10)+ACG(15,12)*PFP(10,15))*ACG(10,17)+(ACG(3,12)*PFP(3,13)+ACG
     #(4,12)*PFP(4,13)+ACG(11,12)*PFP(11,13)+ACG(12,12)*PFP(12,13)+ACG(1
     #5,12)*PFP(13,15))*ACG(13,17)+(ACG(1,12)*PFP(1,14)+ACG(2,12)*PFP(2,
     #14)+ACG(11,12)*PFP(11,14)+ACG(12,12)*PFP(12,14)+ACG(15,12)*PFP(14,
     #15))*ACG(14,17)
      ST(150) = ST(150)+(ACG(1,12)*PFP(1,1)+ACG(2,12)*PFP(1,2))*ACG(1,18
     #)+(ACG(1,12)*PFP(1,2)+ACG(2,12)*PFP(2,2))*ACG(2,18)+(ACG(3,12)*PFP
     #(3,3)+ACG(4,12)*PFP(3,4))*ACG(3,18)+(ACG(3,12)*PFP(3,4)+ACG(4,12)*
     #PFP(4,4))*ACG(4,18)
      ST(151) = ST(151)+(ACG(1,13)*PFP(1,1)+ACG(2,13)*PFP(1,2))*ACG(1,13
     #)+(ACG(1,13)*PFP(1,2)+ACG(2,13)*PFP(2,2))*ACG(2,13)
      ST(152) = ST(152)
      ST(153) = ST(153)+(ACG(1,13)*PFP(1,9)+ACG(2,13)*PFP(2,9))*ACG(9,15
     #)+(ACG(1,13)*PFP(1,10)+ACG(2,13)*PFP(2,10))*ACG(10,15)+(ACG(1,13)*
     #PFP(1,14)+ACG(2,13)*PFP(2,14))*ACG(14,15)
      ST(154) = ST(154)+(ACG(1,13)*PFP(1,9)+ACG(2,13)*PFP(2,9))*ACG(9,16
     #)+(ACG(1,13)*PFP(1,10)+ACG(2,13)*PFP(2,10))*ACG(10,16)+(ACG(1,13)*
     #PFP(1,14)+ACG(2,13)*PFP(2,14))*ACG(14,16)
      ST(155) = ST(155)+(ACG(1,13)*PFP(1,9)+ACG(2,13)*PFP(2,9))*ACG(9,17
     #)+(ACG(1,13)*PFP(1,10)+ACG(2,13)*PFP(2,10))*ACG(10,17)+(ACG(1,13)*
     #PFP(1,14)+ACG(2,13)*PFP(2,14))*ACG(14,17)
      ST(156) = ST(156)+(ACG(1,13)*PFP(1,1)+ACG(2,13)*PFP(1,2))*ACG(1,18
     #)+(ACG(1,13)*PFP(1,2)+ACG(2,13)*PFP(2,2))*ACG(2,18)
      ST(157) = ST(157)+(ACG(3,14)*PFP(3,3)+ACG(4,14)*PFP(3,4))*ACG(3,14
     #)+(ACG(3,14)*PFP(3,4)+ACG(4,14)*PFP(4,4))*ACG(4,14)
      ST(158) = ST(158)+(ACG(3,14)*PFP(3,7)+ACG(4,14)*PFP(4,7))*ACG(7,15
     #)+(ACG(3,14)*PFP(3,8)+ACG(4,14)*PFP(4,8))*ACG(8,15)+(ACG(3,14)*PFP
     #(3,13)+ACG(4,14)*PFP(4,13))*ACG(13,15)
      ST(159) = ST(159)+(ACG(3,14)*PFP(3,7)+ACG(4,14)*PFP(4,7))*ACG(7,16
     #)+(ACG(3,14)*PFP(3,8)+ACG(4,14)*PFP(4,8))*ACG(8,16)+(ACG(3,14)*PFP
     #(3,13)+ACG(4,14)*PFP(4,13))*ACG(13,16)
      ST(160) = ST(160)+(ACG(3,14)*PFP(3,7)+ACG(4,14)*PFP(4,7))*ACG(7,17
     #)+(ACG(3,14)*PFP(3,8)+ACG(4,14)*PFP(4,8))*ACG(8,17)+(ACG(3,14)*PFP
     #(3,13)+ACG(4,14)*PFP(4,13))*ACG(13,17)
      ST(161) = ST(161)+(ACG(3,14)*PFP(3,3)+ACG(4,14)*PFP(3,4))*ACG(3,18
     #)+(ACG(3,14)*PFP(3,4)+ACG(4,14)*PFP(4,4))*ACG(4,18)
      ST(162) = ST(162)+(ACG(5,15)*PFP(5,5)+ACG(6,15)*PFP(5,6))*ACG(5,15
     #)+(ACG(5,15)*PFP(5,6)+ACG(6,15)*PFP(6,6))*ACG(6,15)
      ST(163) = ST(163)+(ACG(5,15)*PFP(5,5)+ACG(6,15)*PFP(5,6))*ACG(5,16
     #)+(ACG(5,15)*PFP(5,6)+ACG(6,15)*PFP(6,6))*ACG(6,16)
      ST(164) = ST(164)+(ACG(5,15)*PFP(5,5)+ACG(6,15)*PFP(5,6))*ACG(5,17
     #)+(ACG(5,15)*PFP(5,6)+ACG(6,15)*PFP(6,6))*ACG(6,17)
      ST(165) = ST(165)+(ACG(9,15)*PFP(1,9)+ACG(10,15)*PFP(1,10)+ACG(14,
     #15)*PFP(1,14))*ACG(1,18)+(ACG(9,15)*PFP(2,9)+ACG(10,15)*PFP(2,10)+
     #ACG(14,15)*PFP(2,14))*ACG(2,18)+(ACG(7,15)*PFP(3,7)+ACG(8,15)*PFP(
     #3,8)+ACG(13,15)*PFP(3,13))*ACG(3,18)+(ACG(7,15)*PFP(4,7)+ACG(8,15)
     #*PFP(4,8)+ACG(13,15)*PFP(4,13))*ACG(4,18)+(ACG(13,15)*PFP(11,13)+A
     #CG(14,15)*PFP(11,14))*ACG(11,18)+(ACG(13,15)*PFP(12,13)+ACG(14,15)
     #*PFP(12,14))*ACG(12,18)+(ACG(7,15)*PFP(7,15)+ACG(8,15)*PFP(8,15)+A
     #CG(9,15)*PFP(9,15)+ACG(10,15)*PFP(10,15)+ACG(13,15)*PFP(13,15)+ACG
     #(14,15)*PFP(14,15))*ACG(15,18)
      ST(166) = ST(166)+(ACG(5,16)*PFP(5,5)+ACG(6,16)*PFP(5,6))*ACG(5,16
     #)+(ACG(5,16)*PFP(5,6)+ACG(6,16)*PFP(6,6))*ACG(6,16)
      ST(167) = ST(167)+(ACG(5,16)*PFP(5,5)+ACG(6,16)*PFP(5,6))*ACG(5,17
     #)+(ACG(5,16)*PFP(5,6)+ACG(6,16)*PFP(6,6))*ACG(6,17)
      ST(168) = ST(168)+(ACG(9,16)*PFP(1,9)+ACG(10,16)*PFP(1,10)+ACG(14,
     #16)*PFP(1,14))*ACG(1,18)+(ACG(9,16)*PFP(2,9)+ACG(10,16)*PFP(2,10)+
     #ACG(14,16)*PFP(2,14))*ACG(2,18)+(ACG(7,16)*PFP(3,7)+ACG(8,16)*PFP(
     #3,8)+ACG(13,16)*PFP(3,13))*ACG(3,18)+(ACG(7,16)*PFP(4,7)+ACG(8,16)
     #*PFP(4,8)+ACG(13,16)*PFP(4,13))*ACG(4,18)+(ACG(13,16)*PFP(11,13)+A
     #CG(14,16)*PFP(11,14))*ACG(11,18)+(ACG(13,16)*PFP(12,13)+ACG(14,16)
     #*PFP(12,14))*ACG(12,18)+(ACG(7,16)*PFP(7,15)+ACG(8,16)*PFP(8,15)+A
     #CG(9,16)*PFP(9,15)+ACG(10,16)*PFP(10,15)+ACG(13,16)*PFP(13,15)+ACG
     #(14,16)*PFP(14,15))*ACG(15,18)
      ST(169) = ST(169)+(ACG(5,17)*PFP(5,5)+ACG(6,17)*PFP(5,6))*ACG(5,17
     #)+(ACG(5,17)*PFP(5,6)+ACG(6,17)*PFP(6,6))*ACG(6,17)
      ST(170) = ST(170)+(ACG(9,17)*PFP(1,9)+ACG(10,17)*PFP(1,10)+ACG(14,
     #17)*PFP(1,14))*ACG(1,18)+(ACG(9,17)*PFP(2,9)+ACG(10,17)*PFP(2,10)+
     #ACG(14,17)*PFP(2,14))*ACG(2,18)+(ACG(7,17)*PFP(3,7)+ACG(8,17)*PFP(
     #3,8)+ACG(13,17)*PFP(3,13))*ACG(3,18)+(ACG(7,17)*PFP(4,7)+ACG(8,17)
     #*PFP(4,8)+ACG(13,17)*PFP(4,13))*ACG(4,18)+(ACG(13,17)*PFP(11,13)+A
     #CG(14,17)*PFP(11,14))*ACG(11,18)+(ACG(13,17)*PFP(12,13)+ACG(14,17)
     #*PFP(12,14))*ACG(12,18)+(ACG(7,17)*PFP(7,15)+ACG(8,17)*PFP(8,15)+A
     #CG(9,17)*PFP(9,15)+ACG(10,17)*PFP(10,15)+ACG(13,17)*PFP(13,15)+ACG
     #(14,17)*PFP(14,15))*ACG(15,18)
      ST(171) = ST(171)+(ACG(1,18)*PFP(1,1)+ACG(2,18)*PFP(1,2))*ACG(1,18
     #)+(ACG(1,18)*PFP(1,2)+ACG(2,18)*PFP(2,2))*ACG(2,18)+(ACG(3,18)*PFP
     #(3,3)+ACG(4,18)*PFP(3,4))*ACG(3,18)+(ACG(3,18)*PFP(3,4)+ACG(4,18)*
     #PFP(4,4))*ACG(4,18)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE SHMDSP3(COORD,COORDI,EDIS,VR,VS,VT,REDIS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     PURPOSE:	MODIFIES TOTAL DISPLACEMENT VECTOR BY DEDUCTING
C				RIGID BODY TRANSLATIONS AND ROTATIONS
C
C	INPUT VARIABLES
C     COORD(3,NNO)      = CURRENT NODAL COORDINATES
C     COORDI(3,NNO)     = INITIAL NODAL COORDINATES
C     EDIS(NEF)         = CURRENT NODAL DISPLACEMENTS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     VR(3),VS(3),VT(3) = CURRENT DIRECTION COSINE VECTORS
C
C	LOCAL VARIABLES
C	ARC				  = 2*PI
C	DUM(2,4)		  = DUMMY VARIABLE FOR INITIAL LOCAL COORD.
C
C	OUTPUT VARIABLES
C     REDIS(48)         = COROTATIONAL FORM OF EDIS
C     ----------------------------------------------------
      DIMENSION COORD(3,3),COORDI(3,3),EDIS(18),VR(3),VS(3),VT(3)
      DIMENSION TM(3,3),XJI(3),CDO(3,3),DUM(2,3),XYZ(3),DUM2(3)
      DIMENSION VV(3),VRO(3),VSO(3),VTO(3),ROT(3),ROV(3)
	DIMENSION REDIS(18)
C
      ARC=6.2831853071796

C     ---------------------------------------------------------------
C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
C     ---------------------------------------------------------------
C	DETERMINE ROTATIONAL COMPONENTS A THE ELEMENT CENTER/REF POINT
      CALL CLEARA (ROT,3)
C      ROT(1)=(EDIS(4)+EDIS(10)+EDIS(16))/3
C      ROT(2)=(EDIS(5)+EDIS(11)+EDIS(17))/3
C      ROT(3)=(EDIS(6)+EDIS(12)+EDIS(18))/3
      ROT(1)=EDIS(4)
      ROT(2)=EDIS(5)
      ROT(3)=EDIS(6)

C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
      CALL SCALEN (ROT,ROV,RO,3)
CCC      IF (RO.EQ.0.) RETURN

C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
C	ELEMENT INITIAL CENTER COORD/INITIAL REF PT COORD - XYZ
      CALL CLEARA (XYZ,3)
	XYZ(1)=COORDI(1,1)
	XYZ(2)=COORDI(2,1)
	XYZ(3)=COORDI(3,1)
C	XYZ(1)=(COORDI(1,1)+COORDI(1,2)+COORDI(1,3))/3
C	XYZ(2)=(COORDI(2,1)+COORDI(2,2)+COORDI(2,3))/3
C	XYZ(3)=(COORDI(3,1)+COORDI(3,2)+COORDI(3,3))/3

C	INITIAL NODAL DISTANCE RELATIVE TO INITIAL REF POINT - CDO
      DO 150 I=1,NNO
      DO 150 J=1,3
150	CDO(J,I)=COORDI(J,I)-XYZ(J)

C	FIND INITIAL COORDINATE DIRECTION VECTORS - VRO,VSO,VTO
	CALL TRIVRS(COORDI,VRO,VSO,VTO,DUM,NNO)

C	DETERMINE ORTHOGONAL ROTATION MATRIX - TM
      DO 160 I=1,3
      VV(I)=ROV(I)
      DO 160 J=1,3
160	TM(I,J)=VR(I)*VRO(J)+VS(I)*VSO(J)+VT(I)*VTO(J)

C	DETERMINE ANGLE OF ROTATION - RR
      ABCD=(.5*(TM(1,1)+TM(2,2)+TM(3,3)-1.))
	IF (ABCD .GT.  1.00000000000000000 ) ABCD= 1.000000000000000000
	IF (ABCD .LT. -1.00000000000000000 ) ABCD=-1.000000000000000000
      RR=ACOS(ABCD)
      SN=SIN(RR)
      IF (SN.EQ.0.) GOTO 190

C	DETERMINE ROTATION AXIS VECTOR - VV
      F=0.5/SN
      VV(1)=F*(TM(3,2)-TM(2,3))
      VV(2)=F*(TM(1,3)-TM(3,1))
      VV(3)=F*(TM(2,1)-TM(1,2))
      CALL SCAPRD (ROV,VV,CS,3)
      IF (CS.GE.0.) GOTO 190
      RR=-RR
      DO 180 I=1,3
180	VV(I)=-VV(I)

C	???
190	RD=RO-RR
      N=(RD+SIGN(ARC/2.01,RD))/ARC
      RR=RR+N*ARC

C	DETERMINE CO-ROTATIONAL DISPLACEMENTS
      K=1
      DO 220 I=1,NNO
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)
c      REDIS(K)=COORD(J,I)-TCDO
      REDIS(K)= (COORD(J,I)-COORD(J,1))-TCDO
      REDIS(K+3)=EDIS(K+3)-RR*VV(J)
210	K=K+1
220	K=K+3

      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHMDSP3B(COORD,COORDI,EDIS,VR,VS,VT,REDIS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     PURPOSE:	MODIFIES TOTAL DISPLACEMENT VECTOR BY DEDUCTING
C				RIGID BODY TRANSLATIONS AND ROTATIONS
C
C	INPUT VARIABLES
C     COORD(3,NNO)      = CURRENT NODAL COORDINATES
C     COORDI(3,NNO)     = INITIAL NODAL COORDINATES
C     EDIS(NEF)         = CURRENT NODAL DISPLACEMENTS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     VR(3),VS(3),VT(3) = CURRENT DIRECTION COSINE VECTORS
C
C	LOCAL VARIABLES
C	ARC				  = 2*PI
C	DUM(2,4)		  = DUMMY VARIABLE FOR INITIAL LOCAL COORD.
C
C	OUTPUT VARIABLES
C     REDIS(48)         = COROTATIONAL FORM OF EDIS
C     ----------------------------------------------------
      DIMENSION COORD(3,3),COORDI(3,3),EDIS(18),VR(3),VS(3),VT(3)
      DIMENSION TM(3,3),XJI(3),CDO(3,3),DUM(2,3),XYZ(3),DUM2(3)
      DIMENSION VV(3),VRO(3),VSO(3),VTO(3),ROT(3),ROV(3)
	DIMENSION REDIS(18)
	DIMENSION V15(3),XYZI(3)
C
      ARC=6.2831853071796

C     ---------------------------------------------------------------
C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
C     ---------------------------------------------------------------
C	DETERMINE ROTATIONAL COMPONENTS A THE ELEMENT CENTER/REF POINT
      CALL CLEARA (ROT,3)
      ROT(1)=(EDIS(4)+EDIS(10)+EDIS(16))/3
      ROT(2)=(EDIS(5)+EDIS(11)+EDIS(17))/3
      ROT(3)=(EDIS(6)+EDIS(12)+EDIS(18))/3

C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
      CALL SCALEN (ROT,ROV,RO,3)
      IF (RO.EQ.0.) RETURN


C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
C
C	ELEMENT INITIAL CENTER COORD/INITIAL REF PT COORD - XYZ
C     -------------------------------------------------------
C	COMPUTE FOR VECTOR PASSING THROUGH NODE 1 AND MIDPOINT OF SIDE 2-3
C	------------------------------------------------------------------
	V15(1) = (COORDI(1,2)+COORDI(1,3))/2.0d0 - COORDI(1,1)
	V15(2) = (COORDI(2,2)+COORDI(2,3))/2.0d0 - COORDI(2,1)
	V15(3) = (COORDI(3,2)+COORDI(3,3))/2.0d0 - COORDI(3,1)
	CALL SCALEN(V15,V15,DUMMY,3)

C	COMPUTE FOR ELEMENT CENTER COORDINATES - EC
C	-------------------------------------------
	XYZI(1) = COORDI(1,1) + V15(1)*DUMMY*(2.0D0/3.0D0)
	XYZI(2) = COORDI(2,1) + V15(2)*DUMMY*(2.0D0/3.0D0)
	XYZI(3) = COORDI(3,1) + V15(3)*DUMMY*(2.0D0/3.0D0)

C	INITIAL NODAL DISTANCE RELATIVE TO INITIAL REF POINT - CDO
      DO 150 I=1,NNO
      DO 150 J=1,3
150	CDO(J,I)=COORDI(J,I)-XYZI(J)

C	FIND INITIAL COORDINATE DIRECTION VECTORS - VRO,VSO,VTO
	CALL TRIVRS2(COORDI,VRO,VSO,VTO,DUM,NNO)

C	DETERMINE ORTHOGONAL ROTATION MATRIX - TM
      DO 160 I=1,3
      VV(I)=ROV(I)
      DO 160 J=1,3
160	TM(I,J)=VR(I)*VRO(J)+VS(I)*VSO(J)+VT(I)*VTO(J)

C	DETERMINE ANGLE OF ROTATION - RR
      ABCD=(.5*(TM(1,1)+TM(2,2)+TM(3,3)-1.))
	IF (ABCD .GT.  1.00000000000000000 ) ABCD= 1.000000000000000000
	IF (ABCD .LT. -1.00000000000000000 ) ABCD=-1.000000000000000000
      RR=ACOS(ABCD)
      SN=SIN(RR)
      IF (SN.EQ.0.) GOTO 190

C	DETERMINE ROTATION AXIS VECTOR - VV
      F=0.5/SN
      VV(1)=F*(TM(3,2)-TM(2,3))
      VV(2)=F*(TM(1,3)-TM(3,1))
      VV(3)=F*(TM(2,1)-TM(1,2))
      CALL SCAPRD (ROV,VV,CS,3)
      IF (CS.GE.0.) GOTO 190
      RR=-RR
      DO 180 I=1,3
180	VV(I)=-VV(I)

C	???
190	RD=RO-RR
      N=(RD+SIGN(ARC/2.01,RD))/ARC
      RR=RR+N*ARC

C	ELEMENT NEW CENTER COORD/NEW REF PT COORD - XYZ
C     -------------------------------------------------------
C	COMPUTE FOR VECTOR PASSING THROUGH NODE 1 AND MIDPOINT OF SIDE 2-3
C	------------------------------------------------------------------
	V15(1) = (COORD(1,2)+COORD(1,3))/2.0d0 - COORD(1,1)
	V15(2) = (COORD(2,2)+COORD(2,3))/2.0d0 - COORD(2,1)
	V15(3) = (COORD(3,2)+COORD(3,3))/2.0d0 - COORD(3,1)
	CALL SCALEN(V15,V15,DUMMY,3)

C	COMPUTE FOR ELEMENT CENTER COORDINATES - EC
C	-------------------------------------------
	XYZ(1) = COORD(1,1) + V15(1)*DUMMY*(2.0D0/3.0D0)
	XYZ(2) = COORD(2,1) + V15(2)*DUMMY*(2.0D0/3.0D0)
	XYZ(3) = COORD(3,1) + V15(3)*DUMMY*(2.0D0/3.0D0)

C	DETERMINE CO-ROTATIONAL DISPLACEMENTS
      K=1
      DO 220 I=1,NNO
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)
c      REDIS(K)=COORD(J,I)-TCDO
      REDIS(K)= (COORD(J,I)-XYZ(J))-TCDO
      REDIS(K+3)=EDIS(K+3)-RR*VV(J)
210	K=K+1
220	K=K+3

      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHMDSP3O(COORD,COORDI,EDIS,VR,VS,VT,REDIS,NNO) !original
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     PURPOSE:	MODIFIES TOTAL DISPLACEMENT VECTOR BY DEDUCTING
C				RIGID BODY TRANSLATIONS AND ROTATIONS
C
C	INPUT VARIABLES
C     COORD(3,NNO)      = CURRENT NODAL COORDINATES
C     COORDI(3,NNO)     = INITIAL NODAL COORDINATES
C     EDIS(NEF)         = CURRENT NODAL DISPLACEMENTS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     VR(3),VS(3),VT(3) = CURRENT DIRECTION COSINE VECTORS
C
C	LOCAL VARIABLES
C	ARC				  = 2*PI
C	DUM(2,4)		  = DUMMY VARIABLE FOR INITIAL LOCAL COORD.
C
C	OUTPUT VARIABLES
C     REDIS(48)         = COROTATIONAL FORM OF EDIS
C     ----------------------------------------------------
      DIMENSION COORD(3,3),COORDI(3,3),EDIS(18),VR(3),VS(3),VT(3)
      DIMENSION TM(3,3),XJI(3),CDO(3,3),DUM(2,3),XYZ(3),DUM2(3)
      DIMENSION VV(3),VRO(3),VSO(3),VTO(3),ROT(3),ROV(3)
	DIMENSION REDIS(18)
C
      ARC=6.2831853071796

C     ---------------------------------------------------------------
C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
C     ---------------------------------------------------------------
C	DETERMINE ROTATIONAL COMPONENTS A THE ELEMENT CENTER/REF POINT
      CALL CLEARA (ROT,3)
      ROT(1)=(EDIS(4)+EDIS(10)+EDIS(16))/3
      ROT(2)=(EDIS(5)+EDIS(11)+EDIS(17))/3
      ROT(3)=(EDIS(6)+EDIS(12)+EDIS(18))/3

C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
      CALL SCALEN (ROT,ROV,RO,3)
      IF (RO.EQ.0.) RETURN

C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
C	ELEMENT INITIAL CENTER COORD/INITIAL REF PT COORD - XYZ
      CALL CLEARA (XYZ,3)
	XYZ(1)=COORDI(1,1)
	XYZ(2)=COORDI(2,1)
	XYZ(3)=COORDI(3,1)
C	XYZ(1)=(COORDI(1,1)+COORDI(1,2)+COORDI(1,3))/3
C	XYZ(2)=(COORDI(2,1)+COORDI(2,2)+COORDI(2,3))/3
C	XYZ(3)=(COORDI(3,1)+COORDI(3,2)+COORDI(3,3))/3

C	INITIAL NODAL DISTANCE RELATIVE TO INITIAL REF POINT - CDO
      DO 150 I=1,NNO
      DO 150 J=1,3
150	CDO(J,I)=COORDI(J,I)-XYZ(J)

C	FIND INITIAL COORDINATE DIRECTION VECTORS - VRO,VSO,VTO
C	NEXT LINE CHANGED BY GILSON - AUG2002
C	CALL VRST(COORDI,VRO,VSO,VTO,DUM,NNO)
	CALL VRST(COORDI,VRO,VSO,VTO,DUM,DUM2,NNO)

C	DETERMINE ORTHOGONAL ROTATION MATRIX - TM
      DO 160 I=1,3
      VV(I)=ROV(I)
      DO 160 J=1,3
160	TM(I,J)=VR(I)*VRO(J)+VS(I)*VSO(J)+VT(I)*VTO(J)

C	DETERMINE ANGLE OF ROTATION - RR
      ABCD=(.5*(TM(1,1)+TM(2,2)+TM(3,3)-1.))
	IF (ABCD .GT.  1.00000000000000000 ) ABCD= 1.000000000000000000
	IF (ABCD .LT. -1.00000000000000000 ) ABCD=-1.000000000000000000
      RR=ACOS(ABCD)
      SN=SIN(RR)
      IF (SN.EQ.0.) GOTO 190

C	DETERMINE ROTATION AXIS VECTOR - VV
      F=0.5/SN
      VV(1)=F*(TM(3,2)-TM(2,3))
      VV(2)=F*(TM(1,3)-TM(3,1))
      VV(3)=F*(TM(2,1)-TM(1,2))
      CALL SCAPRD (ROV,VV,CS,3)
      IF (CS.GE.0.) GOTO 190
      RR=-RR
      DO 180 I=1,3
180	VV(I)=-VV(I)

C	???
190	RD=RO-RR
      N=(RD+SIGN(ARC/2.01,RD))/ARC
      RR=RR+N*ARC

C	DETERMINE CO-ROTATIONAL DISPLACEMENTS
      K=1
      DO 220 I=1,NNO
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)
      REDIS(K)=COORD(J,I)-TCDO
      REDIS(K+3)=EDIS(K+3)-RR*VV(J)
210	K=K+1
220	K=K+3

      RETURN
      END
C
C=====================================================================
	SUBROUTINE NGFPG18(AB1,APM,APB,APS,DR,PFP)

	IMPLICIT REAL*8 (A-H,O-Z)

C	--------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Ng*F*Pg - PFP
C
C	INPUT VARIABLES
C	AB1(4,4)	= SUB-MATRIX OF INTEGRAL{TRANS(Pb)*Pb}
C	APM(9)		= MEMBRANE STRAIN PARAMETER
C	APB(11)		= BENDING STRAIN PARAMETER
C	APS(4)		= TRANSVERSE SHEAR STRAIN PARAMETER
C	DR(64)		= MODULI MATRIX
C
C
C	OUTPUT VARIAVBLE
C	PFP(51,51)	= Ng*F*Pg
C	--------------------------------------------------

      DIMENSION AB1(3,3),APM(7),APB(9),APS(4),DR(64)
	DIMENSION PFP(15,15)

	dvol = 1.0d0
	SLR = 1.0d0
	SLS = 1.0d0
	CALL SHDELA (DR,DVOL,SLR,SLS)

C	-----------------------
C	SOLVE FOR Ng*F*Pg - PFP
C	-----------------------
      t3 = AB1(1,3)*APM(3)+AB1(1,1)*APM(1)+AB1(1,2)*APM(2)
      t7 = AB1(1,2)*APM(5)+AB1(1,3)*APM(6)+AB1(1,1)*APM(4)
      t9 = t3*DR(1)+t7*DR(2)
      t13 = (-AB1(1,3)*APM(2)-AB1(1,2)*APM(6)+AB1(1,1)*APM(7))*DR(19)
      t16 = AB1(1,3)*APB(3)+AB1(1,1)*APB(1)+AB1(1,2)*APB(2)
      t20 = AB1(1,2)*APB(5)+AB1(1,3)*APB(6)+AB1(1,1)*APB(4)
      t22 = t16*DR(28)+t20*DR(29)
      t25 = AB1(1,1)*APB(7)+AB1(1,2)*APB(8)+AB1(1,3)*APB(9)
      t26 = t25*DR(46)
      t27 = AB1(1,1)*(APS(1)+APS(2))
      t28 = DR(55)*t27
      t31 = t7*DR(1)+t3*DR(2)
      t34 = t20*DR(28)+t16*DR(29)
      t35 = AB1(1,1)*(APS(3)+APS(4))
      t36 = DR(55)*t35
      t39 = -t16*DR(28)-t20*DR(29)
      t40 = -t25*DR(46)
      t41 = -DR(55)*t27
      t44 = -t20*DR(28)-t16*DR(29)
      t45 = -DR(55)*t35
      t48 = t16*DR(28)/2+t20*DR(29)/2
      t49 = t25*DR(46)/2
      t52 = t20*DR(28)/2+t16*DR(29)/2
      t53 = DR(55)*t27/2
      t54 = DR(55)*t35/2
      PFP(1,1) = t9
      PFP(1,2) = t13
      PFP(1,9) = t22
      PFP(1,10) = t26
      PFP(1,14) = t28
      PFP(2,1) = t13
      PFP(2,2) = t31
      PFP(2,9) = t26
      PFP(2,10) = t34
      PFP(2,14) = t36
      PFP(3,3) = t9
      PFP(3,4) = t13
      PFP(3,7) = t39
      PFP(3,8) = t40
      PFP(3,13) = t41
      PFP(4,3) = t13
      PFP(4,4) = t31
      PFP(4,7) = t40
      PFP(4,8) = t44
      PFP(4,13) = t45
      PFP(5,5) = t9
      PFP(5,6) = t13
      PFP(6,5) = t13
      PFP(6,6) = t31
      PFP(7,3) = t39
      PFP(7,4) = t40
      PFP(7,15) = t48
      PFP(8,3) = t40
      PFP(8,4) = t44
      PFP(8,15) = t49
      PFP(9,1) = t22
      PFP(9,2) = t26
      PFP(9,15) = t49
      PFP(10,1) = t26
      PFP(10,2) = t34
      PFP(10,15) = t52
      PFP(11,13) = t48
      PFP(11,14) = t49
      PFP(12,13) = t49
      PFP(12,14) = t52
      PFP(13,3) = t41
      PFP(13,4) = t45
      PFP(13,11) = t48
      PFP(13,12) = t49
      PFP(13,15) = t53
      PFP(14,1) = t28
      PFP(14,2) = t36
      PFP(14,11) = t49
      PFP(14,12) = t52
      PFP(14,15) = t54
      PFP(15,7) = t48
      PFP(15,8) = t49
      PFP(15,9) = t49
      PFP(15,10) = t52
      PFP(15,13) = t53
      PFP(15,14) = t54

	RETURN
	END
C
C=====================================================================
	SUBROUTINE TRICF_B(RS,BL,BLR2,BLS2,SLN,RN,SN,FLR,FLS,FLSS,FLRR,
	1                   DWR,DWRS,DWS,DWSR,ANT)
	IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION RS(2,3),BL(3),FLR(3),FLS(3),FLSS(3),FLRR(3)
	DIMENSION SLN(3),SN(3),RN(3)
	DIMENSION DWR(9),DWRS(9),DWS(9),DWSR(9),ANT(3)


	BLR = BLR2
	BLS = BLS2

	F1 = 1.0D0-BLR
	F2 = 1.0D0-BLS

C	-------------------
C	AREA OF THE ELEMENT
C	-------------------
	AREA = 0.50*( (RS(1,1)*RS(2,2)+RS(1,2)*RS(2,3)+RS(1,3)*RS(2,1))-
     #              (RS(2,1)*RS(1,2)+RS(2,2)*RS(1,3)+RS(2,3)*RS(1,1))  )

C	------------------------------------------------
C	COMPUTE FOR COMMON FACTORS FOR AREA INTEGRATIONS
C	------------------------------------------------

C	AREA INTEGRATION FACTORS - ANT
C	------------------------------
	ANT(1) = 
     #((RS(2,2)/6-RS(2,3)/6)*RS(1,1)+(RS(2,3)/6-RS(2,1)/6)*RS(1,2)+
     #(-RS(2,2)/6+RS(2,1)/6)*RS(1,3)) !*dp(1)
	ANT(2) = 
     #((RS(2,2)/6-RS(2,3)/6)*RS(1,1)+(RS(2,3)/6-RS(2,1)/6)*RS(1,2)+
     #(-RS(2,2)/6+RS(2,1)/6)*RS(1,3)) !*dp(2)
	ANT(3) = 
     #((RS(2,2)/6-RS(2,3)/6)*RS(1,1)+(RS(2,3)/6-RS(2,1)/6)*RS(1,2)+
     #(-RS(2,2)/6+RS(2,1)/6)*RS(1,3)) !*dp(3)

C	FACTORS FOR ROTATION ABOUT S
C	----------------------------
	FLS(1) = F1*ANT(1)
	FLS(2) = F1*ANT(2)
	FLS(3) = F1*ANT(3)

C	FACTORS FOR ROTATION ABOUT R
C	----------------------------
	FLR(1) = -F2*ANT(1)
	FLR(2) = -F2*ANT(2)
	FLR(3) = -F2*ANT(3)

C	FACTORS FOR ROTATION ABOUT S MULTIPLIED BY S
C	--------------------------------------------
	FLSS(1)= F1*(
     #(((RS(2,2)/12-RS(2,3)/12)*RS(2,1)+RS(2,2)**2/24-RS(2,3)**2/24
     #)*RS(1,1)+(-RS(2,1)**2/12+(RS(2,3)/24-RS(2,2)/24)*RS(2,1)+RS(2,3)*
     #*2/24+RS(2,2)*RS(2,3)/24)*RS(1,2)+(RS(2,1)**2/12+(RS(2,3)/24-RS(2,
     #2)/24)*RS(2,1)-RS(2,2)**2/24-RS(2,2)*RS(2,3)/24)*RS(1,3)) ) !*dp(1)
	FLSS(2)= F1*(
     #(((-RS(2,3)/24+RS(2,2)/24)*RS(2,1)-RS(2,3)**2/24+RS(2,2)**2/1
     #2-RS(2,2)*RS(2,3)/24)*RS(1,1)+(RS(2,3)**2/24-RS(2,1)**2/24+RS(2,2)
     #*RS(2,3)/12-RS(2,1)*RS(2,2)/12)*RS(1,2)+(RS(2,1)**2/24+(RS(2,3)/24
     #+RS(2,2)/24)*RS(2,1)-RS(2,2)*RS(2,3)/24-RS(2,2)**2/12)*RS(1,3)) ) !*dp(2)
	FLSS(3)= F1*(
     #(((-RS(2,3)/24+RS(2,2)/24)*RS(2,1)-RS(2,3)**2/12+RS(2,2)*RS(2,
     #3)/24+RS(2,2)**2/24)*RS(1,1)+(-RS(2,1)**2/24+(-RS(2,3)/24-RS(2,2)/
     #24)*RS(2,1)+RS(2,2)*RS(2,3)/24+RS(2,3)**2/12)*RS(1,2)+(RS(2,1)**2/
     #24-RS(2,2)**2/24-RS(2,2)*RS(2,3)/12+RS(2,1)*RS(2,3)/12)*RS(1,3)) ) !*dp(3)

C	FACTORS FOR ROTATION ABOUT R MULTIPLIED BY R
C	--------------------------------------------
	FLRR(1)= -F2*(
     #((RS(2,2)/12-RS(2,3)/12)*RS(1,1)**2+((-RS(2,1)/12+RS(2,2)/24+
     #RS(2,3)/24)*RS(1,2)+(-RS(2,3)/24+RS(2,1)/12-RS(2,2)/24)*RS(1,3))*R
     #S(1,1)+(-RS(2,1)/24+RS(2,3)/24)*RS(1,2)**2+(RS(2,3)/24-RS(2,2)/24)
     #*RS(1,3)*RS(1,2)+(-RS(2,2)/24+RS(2,1)/24)*RS(1,3)**2) ) !*dp(1)
	FLRR(2)= -F2*(
     #((-RS(2,3)/24+RS(2,2)/24)*RS(1,1)**2+((-RS(2,1)/24-RS(2,3)/24
     #+RS(2,2)/12)*RS(1,2)+(-RS(2,3)/24+RS(2,1)/24)*RS(1,3))*RS(1,1)+(RS
     #(2,3)/12-RS(2,1)/12)*RS(1,2)**2+(-RS(2,2)/12+RS(2,3)/24+RS(2,1)/24
     #)*RS(1,3)*RS(1,2)+(-RS(2,2)/24+RS(2,1)/24)*RS(1,3)**2) ) !*dp(2)
	FLRR(3)= -F2*(
     #((-RS
     #(2,3)/24+RS(2,2)/24)*RS(1,1)**2+((RS(2,2)/24-RS(2,1)/24)*RS(1,2)+(
     #RS(2,1)/24+RS(2,2)/24-RS(2,3)/12)*RS(1,3))*RS(1,1)+(-RS(2,1)/24+RS
     #(2,3)/24)*RS(1,2)**2+(-RS(2,2)/24+RS(2,3)/12-RS(2,1)/24)*RS(1,3)*R
     #S(1,2)+(RS(2,1)/12-RS(2,2)/12)*RS(1,3)**2) ) !*dp(3)

C	FACTORS FOR W WITH RN
C	---------------------
	DWR(1) =(-RN(3)**2*SLN(3)**2/12+RN(1)**2*SLN(1)**2/12)*(-BLR) !*dpr(1)
	DWR(2) =(RN(2)**2*SLN(2)**2/12-RN(1)**2*SLN(1)**2/12)*(-BLR) !*dpr(2)
	DWR(3) =(-RN(2)**2*SLN(2)**2/12+RN(3)**2*SLN(3)**2/12)*(-BLR) !*dpr(3)

	DWR(4) =(RN(1)*SLN(1)**2*SN(1)/12-RN(3)*SLN(3)**2*SN(3)/12)*(-BLR) !*dps(1)
	DWR(5) =(RN(2)*SLN(2)**2*SN(2)/12-RN(1)*SLN(1)**2*SN(1)/12)*(-BLR) !*dps(2)
	DWR(6) =(RN(3)*SLN(3)**2*SN(3)/12-RN(2)*SLN(2)**2*SN(2)/12)*(-BLR) !*dps(3)

	DWR(7) =(SLN(3)*RN(3)/2+SLN(1)*RN(1)/2)*(-BLR) !*dw(1)
	DWR(8) =(SLN(2)*RN(2)/2+SLN(1)*RN(1)/2)*(-BLR) !*dw(2)
	DWR(9) =(SLN(3)*RN(3)/2+SLN(2)*RN(2)/2)*(-BLR) !*dw(3)

C	FACTORS FOR W WITH RN MULTIPLIED BY S
C	-------------------------------------
	DWRS(1)= -BLR*(
     #((RS(2,2)/30+RS(2,1)/20)*SLN(1)**2*RN(1)**2+(-RS(2,1)/20-RS(2
     #,3)/30)*SLN(3)**2*RN(3)**2) ) !*dpr(1)
	DWRS(2)= -BLR*(
     #((-RS(2,2)/20-RS(2,1)/30)*SLN(1
     #)**2*RN(1)**2+(RS(2,3)/30+RS(2,2)/20)*SLN(2)**2*RN(2)**2) ) !*dpr(2)
	DWRS(3)= -BLR*(
     #(
     #(-RS(2,3)/20-RS(2,2)/30)*SLN(2)**2*RN(2)**2+(RS(2,1)/30+RS(2,3)/20
     #)*SLN(3)**2*RN(3)**2) ) !*dpr(3)

	DWRS(4)= -BLR*(
     #((RS(2,2)/30+RS(2,1)/20)*SLN(1)**2*SN
     #(1)*RN(1)+(-RS(2,1)/20-RS(2,3)/30)*SLN(3)**2*SN(3)*RN(3)) ) !*dps(1)
	DWRS(5)= -BLR*(
     #((-RS(2,2)/20-RS(2,1)/30)*SLN(1)**2*SN(1)*RN(1)+(RS(2,3)/3
     #0+RS(2,2)/20)*SLN(2)**2*SN(2)*RN(2)) ) !*dps(2)
	DWRS(6)= -BLR*(
     #((-RS(2,3)/20-RS(2,2)/
     #30)*SLN(2)**2*SN(2)*RN(2)+(RS(2,1)/30+RS(2,3)/20)*SLN(3)**2*SN(3)*
     #RN(3)) ) !*dps(3)

	DWRS(7)= -BLR*(
     #((7.E0/20.E0*RS(2,1)+3.E0/20.E0*RS(2,2))*SLN(1)*RN(1
     #)+(3.E0/20.E0*RS(2,3)+7.E0/20.E0*RS(2,1))*SLN(3)*RN(3)) ) !*dw(1)
	DWRS(8)= -BLR*(
     #((3.
     #E0/20.E0*RS(2,1)+7.E0/20.E0*RS(2,2))*SLN(1)*RN(1)+(3.E0/20.E0*RS(2
     #,3)+7.E0/20.E0*RS(2,2))*SLN(2)*RN(2)) ) !*dw(2)
	DWRS(9)= -BLR*(
     #((7.E0/20.E0*RS(2,3)+3
     #.E0/20.E0*RS(2,2))*SLN(2)*RN(2)+(3.E0/20.E0*RS(2,1)+7.E0/20.E0*RS(
     #2,3))*SLN(3)*RN(3)) ) !*dw(3)

C	FACTORS FOR W WITH NS
C	---------------------
	DWS(1) =-BLS*((RN(1)*SLN(1)**2*SN(1)/12-RN(3)*SLN(3)**2*SN(3)/12))!*dpr(1)
	DWS(2) =-BLS*((RN(2)*SLN(2)**2*SN(2)/12-RN(1)*SLN(1)**2*SN(1)/12))!*dpr(2)
	DWS(3) =-BLS*((RN(3)*SLN(3)**2*SN(3)/12-RN(2)*SLN(2)**2*SN(2)/12))!*dpr(3)

	DWS(4) =-BLS*((SN(1)**2*SLN(1)**2/12-SN(3)**2*SLN(3)**2/12))!*dps(1)
	DWS(5) =-BLS*((SN(2)**2*SLN(2)**2/12-SN(1)**2*SLN(1)**2/12))!*dps(2)
	DWS(6) =-BLS*((SN(3)**2*SLN(3)**2/12-SN(2)**2*SLN(2)**2/12))!*dps(3)

	DWS(7) =-BLS*((SLN(3)*SN(3)/2+SLN(1)*SN(1)/2))!*dw(1)
	DWS(8) =-BLS*((SLN(2)*SN(2)/2+SLN(1)*SN(1)/2))!*dw(2)
	DWS(9) =-BLS*((SLN(3)*SN(3)/2+SLN(2)*SN(2)/2))!*dw(3)

C	FACTORS FOR W WITH NS MULTIPLIED BY R
C	-------------------------------------
	DWSR(1)= -BLS*(
     #((RS(1,2)/30+RS(1,1)/20)*SLN(1)**2*SN(1)*RN(1)+(-RS(1,3)/30-R
     #S(1,1)/20)*SLN(3)**2*SN(3)*RN(3)) ) !*dpr(1)
	DWSR(2)= -BLS*(
     #((-RS(1,1)/30-RS(1,2)/20)
     #*SLN(1)**2*SN(1)*RN(1)+(RS(1,2)/20+RS(1,3)/30)*SLN(2)**2*SN(2)*RN(
     #2)) ) !*dpr(2)
	DWSR(3)= -BLS*(
     #((-RS(1,3)/20-RS(1,2)/30)*SLN(2)**2*SN(2)*RN(2)+(RS(1,1
     #)/30+RS(1,3)/20)*SLN(3)**2*SN(3)*RN(3)) ) !*dpr(3)

	DWSR(4)= -BLS*(
     #((RS(1,2)/30+RS(1,1
     #)/20)*SLN(1)**2*SN(1)**2+(-RS(1,3)/30-RS(1,1)/20)*SLN(3)**2*SN(3)*
     #*2) ) !*dps(1)
	DWSR(5)= -BLS*(
     #((-RS(1,1)/30-RS(1,2)/20)*SLN(1)**2*SN(1)**2+(RS(1,2)/20+R
     #S(1,3)/30)*SLN(2)**2*SN(2)**2) ) !*dps(2)
	DWSR(6)= -BLS*(
     #((-RS(1,3)/20-RS(1,2)/30)*SL
     #N(2)**2*SN(2)**2+(RS(1,1)/30+RS(1,3)/20)*SLN(3)**2*SN(3)**2) ) !*dps(3)

	DWSR(7)= -BLS*(
     #((3.E0/20.E0*RS(1,2)+7.E0/20.E0*RS(1,1))*SLN(1)*SN(1)+(7.E0/20.E
     #0*RS(1,1)+3.E0/20.E0*RS(1,3))*SLN(3)*SN(3)) ) !*dw(1)
	DWSR(8)= -BLS*(
     #((3.E0/20.E0*RS(
     #1,1)+7.E0/20.E0*RS(1,2))*SLN(1)*SN(1)+(3.E0/20.E0*RS(1,3)+7.E0/20.
     #E0*RS(1,2))*SLN(2)*SN(2)) ) !*dw(2)
	DWSR(9)= -BLS*(
     #((3.E0/20.E0*RS(1,2)+7.E0/20.E0*RS
     #(1,3))*SLN(2)*SN(2)+(7.E0/20.E0*RS(1,3)+3.E0/20.E0*RS(1,1))*SLN(3)
     #*SN(3)) ) !*dw(3)

	RETURN
	END
C
C=====================================================================
	SUBROUTINE SHMDSP3D(COORD,COORDI,EDIS,TEDIS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     PURPOSE:	MODIFIES TOTAL DISPLACEMENT VECTOR BY DEDUCTING
C				RIGID BODY TRANSLATIONS AND ROTATIONS
C
C	INPUT VARIABLES
C     COORD(3,NNO)      = CURRENT NODAL COORDINATES
C     COORDI(3,NNO)     = INITIAL NODAL COORDINATES
C     EDIS(NEF)         = CURRENT NODAL DISPLACEMENTS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C
C	LOCAL VARIABLES
C	A(3,3),AI(3,3)			= SIMULTANEOUS EQUATION COEFFICIENT MATRIX
C	APL(3)					= SINE(ALPHA) ROTATION AXIS
C	ARC						= 2*PI
C	B(3)					= SIMULTANEOUS EQUATION RHS VECTOR
C	Q1(3)					= NODE ROTATION VECTOR
C	QG(9)					= GLOBAL NODE DEFORMATION
C	RSI(2,3),RS(2,3)		= INITIAL AND CURRENT LOCAL COORDINATES
C	UL(6)					= LOCAL DEFORMATION, r AND s
C     VR(3),VS(3),VT(3)		= CURRENT DIRECTION COSINE VECTORS
C     VRO(3),VSO(3),VTO(3)	= INITIAL DIRECTION COSINE VECTORS
C	VN(NODE,VECTOR)			= CURRENT NODE VECTOR DIRECTION
C
C	OUTPUT VARIABLES
C     TEDIS(18)         = COROTATIONAL FORM OF EDIS IN LOCAL COORDINATES
C
C     ----------------------------------------------------
	DIMENSION COORDI(3,3),VRO(3),VSO(3),VTO(3),RSI(2,3)
	DIMENSION COORD(3,3),VR(3),VS(3),VT(3),RS(2,3),UL(6)

	DIMENSION A(3,3),AI(3,3),VN(3,3),B(3),Q1(3),QG(9),APL(3)
	DIMENSION EDIS(18),TEDIS(18)

	DIMENSION VRQ(3),VSQ(3),VTQ(3),VNL(3)
	DIMENSION XB(3),VNT(3)
C
      ARC=6.2831853071796D0

C	ELEMENT INITIAL LOCAL COORDINATES
C	---------------------------------
	CALL TRIVRS2(COORDI,VRO,VSO,VTO,RSI,NNO)

C	ELEMENT CURRENT LOCAL COORDINATES
C	---------------------------------
	CALL TRIVRS2(COORD,VR,VS,VT,RS,NNO)

C	LOCAL DISPLACEMENTS
C	-------------------
	UL(1) = RS(1,1) - RSI(1,1) !UL1
	UL(2) = RS(2,1) - RSI(2,1) !VL1
	UL(3) = RS(1,2) - RSI(1,2)
	UL(4) = RS(2,2) - RSI(2,2)
	UL(5) = RS(1,3) - RSI(1,3)
	UL(6) = RS(2,3) - RSI(2,3)


C	DETERMINE CURRENT NODE VECTOR
C	-----------------------------
	DO I = 1,3 !NODE LOOP
C		NODE ROTATION VECTOR
		Q1(1) = EDIS((I-1)*6+4)
		Q1(2) = EDIS((I-1)*6+5)
		Q1(3) = EDIS((I-1)*6+6)
C		NODE ROTATION ANGLE
		ALPHA = DSQRT(Q1(1)**2+Q1(2)**2+Q1(3)**2)
C		COMPUTE FOR CURRENT NODE VECTOR
		IF (ALPHA.LE.0.0D0) THEN
C			CURRENT NODE VECTOR
			VN(I,1) = VTO(1)
			VN(I,2) = VTO(2)
			VN(I,3) = VTO(3)
		ELSE
C			Solve for Axis VS
			VSQ(1) = Q1(1)/ALPHA
			VSQ(2) = Q1(2)/ALPHA
			VSQ(3) = Q1(3)/ALPHA
C			Solve for VT
			VTQ(1) = VTO(1)
			VTQ(2) = VTO(2)
			VTQ(3) = VTO(3)
C			Solve for VR
			CALL VECPRD(VSQ,VTQ,VRQ)
			CALL SCALEN(VRQ,VRQ,DUM,3)
C			Solve for Orthogonal VT
			CALL VECPRD(VRQ,VSQ,VTQ)
			CALL SCALEN(VTQ,VTQ,DUM,3)
C			Solve for Angle bet Initial VT and Rotation Axis
			BETA = DACOS( VTO(1)*VSQ(1)+VTO(2)*VSQ(2)+VTO(3)*VSQ(3) )
C			Solve for Rotation Radius
			RAD = DSIN(BETA)
C			Define Node New Node Normal in Local
			VNL(1) = RAD*DSIN(ALPHA)
			VNL(2) = VTO(1)*VSQ(1)+VTO(2)*VSQ(2)+VTO(3)*VSQ(3)
			VNL(3) = RAD*DCOS(ALPHA)
C			Transform Point BQ to Global Coordinates
			VNT(1) = ( VRQ(1)*VNL(1)+VSQ(1)*VNL(2)+VTQ(1)*VNL(3) ) 
			VNT(2) = ( VRQ(2)*VNL(1)+VSQ(2)*VNL(2)+VTQ(2)*VNL(3) ) 
			VNT(3) = ( VRQ(3)*VNL(1)+VSQ(3)*VNL(2)+VTQ(3)*VNL(3) ) 
			CALL SCALEN(VNT,VNT,DUM,3)

			VN(I,1) = VNT(1)
			VN(I,2) = VNT(2)
			VN(I,3) = VNT(3)
		ENDIF
	ENDDO


C	CURRENT NODE DEFORMATION
	DO I = 1,3
C		ANGLE OF ROTATION
		ALPHA = VT(1)*VN(I,1)+VT(2)*VN(I,2)+VT(3)*VN(I,3)
		IF (ALPHA.GE. 0.999999999999999D0) ALPHA =  1.0D0
		IF (ALPHA.LE.-0.999999999999999D0) ALPHA = -1.0D0
		ALPHA = DACOS(ALPHA)
C		COMPUTE FOR GLOBAL NODE DEFORMATION
C		SINE(ALPHA) ROTATION AXIS
		APL(1) = ( VT(2)*VN(I,3)-VT(3)*VN(I,2) )
		APL(2) = ( VT(3)*VN(I,1)-VT(1)*VN(I,3) )
		APL(3) = ( VT(1)*VN(I,2)-VT(2)*VN(I,1) ) 
		CALL SCALEN(APL,APL,DUM,3)
C		GLOBAL NODE DEFORMATION
		QG((I-1)*3+1) = ALPHA*APL(1)		
		QG((I-1)*3+2) = ALPHA*APL(2)		
		QG((I-1)*3+3) = ALPHA*APL(3)		
	ENDDO


C	ASSEMBLE LOCAL DEFORMATION VACTOR
	TEDIS(1)  = UL(1)
	TEDIS(2)  = UL(2)
	TEDIS(3)  = 0.0D0
	TEDIS(4)  = VR(1)*QG(1)+VR(2)*QG(2)+VR(3)*QG(3)
	TEDIS(5)  = VS(1)*QG(1)+VS(2)*QG(2)+VS(3)*QG(3)
	TEDIS(6)  = VT(1)*QG(1)+VT(2)*QG(2)+VT(3)*QG(3)
	TEDIS(7)  = UL(3)
	TEDIS(8)  = UL(4)
	TEDIS(9)  = 0.0D0
	TEDIS(10) = VR(1)*QG(4)+VR(2)*QG(5)+VR(3)*QG(6)
	TEDIS(11) = VS(1)*QG(4)+VS(2)*QG(5)+VS(3)*QG(6)
	TEDIS(12) = VT(1)*QG(4)+VT(2)*QG(5)+VT(3)*QG(6)
	TEDIS(13) = UL(5)
	TEDIS(14) = UL(6)
	TEDIS(15) = 0.0D0
	TEDIS(16) = VR(1)*QG(7)+VR(2)*QG(8)+VR(3)*QG(9)
	TEDIS(17) = VS(1)*QG(7)+VS(2)*QG(8)+VS(3)*QG(9)
	TEDIS(18) = VT(1)*QG(7)+VT(2)*QG(8)+VT(3)*QG(9)


1000	FORMAT (9E30.20)

      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE SHAP2D3(R,S,H,P,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION  H(NNO),P(2,NNO)
C
      H = 0.0D0
	P = 0.0D0
	
	X1 = 1.0-R-S
	X2 = R
	X3 = S
	P1R = -1.0D0
	P1S = -1.0D0
	P2R =  1.0D0
	P2S =  0.0D0
	P3R =  0.0D0
	P3S =  1.0D0

	SELECT CASE (NNO)

	CASE(3)
	H(1) = X1
	H(2) = X2
	H(3) = X3
	P(1,1) =  P1R
	P(2,1) =  P1S
	P(1,2) =  P2R
	P(2,2) =  P2S
	P(1,3) =  P3R
	P(2,3) =  P3S

	CASE(6)
	H(1) = X1*(2.0*X1-1.0)
	H(2) = X2*(2.0*X2-1.0)
	H(3) = X3*(2.0*X3-1.0)
	H(4) = 4.0*X1*X2
	H(5) = 4.0*X2*X3
	H(6) = 4.0*X1*X3

	P(1,1) = (4.0*X1-1.0)*P1R
	P(1,2) = (4.0*X2-1.0)*P2R
	P(1,3) = (4.0*X3-1.0)*P3R
	P(1,4) = 4.0*X2*P1R + 4.0*X1*P2R
	P(1,5) = 4.0*X3*P3R + 4.0*X2*P2R
	P(1,6) = 4.0*X3*P3R + 4.0*X1*P1R

	P(2,1) = (4.0*X1-1.0)*P1S
	P(2,2) = (4.0*X2-1.0)*P2S
	P(2,3) = (4.0*X3-1.0)*P3S
	P(2,4) = 4.0*X2*P1S + 4.0*X1*P2S
	P(2,5) = 4.0*X3*P3S + 4.0*X2*P2S
	P(2,6) = 4.0*X3*P3S + 4.0*X1*P1S

	END SELECT

      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE JACO2D3(XY,P,VR,VS,VT,XJ,XJI,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(3,NNO),P(2,NNO),XJ(2,2),XJI(4)
	DIMENSION RS(2,NNO),VR(3),VS(3),VT(3)

	VR(1) = XY(1,2) - XY(1,1)
	VR(2) = XY(2,2) - XY(2,1)
	VR(3) = XY(3,2) - XY(3,1)
	CALL SCALEN(VR,VR,DUM,3)
	VS(1) = XY(1,3) - XY(1,1)
	VS(2) = XY(2,3) - XY(2,1)
	VS(3) = XY(3,3) - XY(3,1)
	CALL SCALEN(VS,VS,DUM,3)
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
	CALL VECPRD(VT,VR,VS)

	DO I = 1,NNO
	RS(1,I) = XY(1,I)*VR(1) + XY(2,I)*VR(2) + XY(3,I)*VR(3) 
	RS(2,I) = XY(1,I)*VS(1) + XY(2,I)*VS(2) + XY(3,I)*VS(3) 
	ENDDO
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  I=1,2
      DO 100  J=1,2
      DUM = 0.0
      DO 90   K=1,NNO
 90   DUM = DUM + P(I,K)*RS(J,K)
 100  XJ(I,J) = DUM
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = XJ(1,1)*XJ(2,2) - XJ(2,1)*XJ(1,2)
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')
C     -----------------------------
C     INVERSE (XJI) OF THE JACOBIAN
C     -----------------------------
      DUM = 1.0/DET
      XJI(1) =  XJ(2,2)*DUM
      XJI(2) = -XJ(2,1)*DUM
      XJI(3) = -XJ(1,2)*DUM
      XJI(4) =  XJ(1,1)*DUM

C	MODIFIED FOR TRIANGULAR SHAPE
	DET = 0.5*DET      
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE SHAP1D3(R,H,P,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------------------------------------------------
      DIMENSION H(NNO),P(NNO)

	H = 0.0
	P = 0.0

	SELECT CASE (NNO)
	
	CASE(2)
	H(1) = 0.5*(1.0 + R)
	H(2) = 0.5*(1.0 - R)
	P(1) = 0.5
	P(2) =-0.5

	CASE(3)
	H(1) =  0.5*R + 0.5*R*R
	H(2) =  1.0   -     R*R
	H(3) = -0.5*R + 0.5*R*R
	P(1) =  0.5   + 1.0*R
	P(2) =        - 2.0*R
	P(3) = -0.5   + 1.0*R


	END SELECT
	

	RETURN

	END



C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE JACO1D3(IFACE,XY,P,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(3,NNO),P(NNO),XJ(3),NN(3)
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R + S = 1.0
C	FACE 2; R= 0.0
C	FACE 3: S= 0.0
C	---------------------

	SELECT CASE(IFACE)
	CASE(2)
	IF(NNO.EQ.2) NN(1:2) = [3,2]
	CASE(3)
	IF(NNO.EQ.2) NN(1:2) = [3,1]
	CASE(1)
	IF(NNO.EQ.2) NN(1:2) = [2,1] 
	END SELECT
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  J=1,3
      DUM = 0.0
      DO 90   K=1,NNO
	KK  = NN(K)
 90   DUM = DUM + P(K)*XY(J,KK)
 100  XJ(J) = DUM
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = SQRT(XJ(1)*XJ(1) + XJ(2)*XJ(2) + XJ(3)*XJ(3)) 
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')   
C
      RETURN
      END


C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE RESHPT3(H,P,IFACE,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------------------------------------------------
      DIMENSION H(NNO),HH(NNO),P(NNO),PP(NNO)


	HH = 0.0
	PP = 0.0

	SELECT CASE(IFACE)

	CASE(2)
	IF(NNO.EQ.3) THEN
	HH(3) = H(1)
	HH(2) = H(2)
	PP(3) = P(1)
	PP(2) = P(2)
	ENDIF

	CASE(3)
	IF(NNO.EQ.3) THEN
	HH(3) = H(1)
	HH(1) = H(2)
	PP(3) = P(1)
	PP(1) = P(2)
	ENDIF

	CASE(1)
	IF(NNO.EQ.3) THEN
	HH(2) = H(1)
	HH(1) = H(2)
	PP(2) = P(1)
	PP(1) = P(2)
	ENDIF

	END SELECT

	H(1:NNO) = HH(1:NNO)
	P(1:NNO) = PP(1:NNO)	

	RETURN

	END



C	=====================================================================
C	=====================================================================	
C	=====================================================================
