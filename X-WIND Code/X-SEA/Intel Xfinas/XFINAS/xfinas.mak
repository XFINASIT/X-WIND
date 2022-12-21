# Microsoft Developer Studio Generated NMAKE File, Based on xfinas.dsp
!IF "$(CFG)" == ""
CFG=xfinas - Win32 Debug
!MESSAGE No configuration specified. Defaulting to xfinas - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "xfinas - Win32 Release" && "$(CFG)" != "xfinas - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "xfinas.mak" CFG="xfinas - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "xfinas - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "xfinas - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "xfinas - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\xfinas.exe" "$(OUTDIR)\xfinas.bsc"

!ELSE 

ALL : "$(OUTDIR)\xfinas.exe" "$(OUTDIR)\xfinas.bsc"

!ENDIF 

CLEAN :
	-@erase "$(INTDIR)\2D Heat main.obj"
	-@erase "$(INTDIR)\2D Heat main.sbr"
	-@erase "$(INTDIR)\3D Heat main.obj"
	-@erase "$(INTDIR)\3D Heat main.sbr"
	-@erase "$(INTDIR)\ADS1.OBJ"
	-@erase "$(INTDIR)\ADS1.SBR"
	-@erase "$(INTDIR)\ADS2.OBJ"
	-@erase "$(INTDIR)\ADS2.SBR"
	-@erase "$(INTDIR)\ADS3.OBJ"
	-@erase "$(INTDIR)\ADS3.SBR"
	-@erase "$(INTDIR)\baseiso.obj"
	-@erase "$(INTDIR)\baseiso.sbr"
	-@erase "$(INTDIR)\beam close section.obj"
	-@erase "$(INTDIR)\beam close section.sbr"
	-@erase "$(INTDIR)\beam21f.obj"
	-@erase "$(INTDIR)\beam21f.sbr"
	-@erase "$(INTDIR)\BOND3D.OBJ"
	-@erase "$(INTDIR)\BOND3D.SBR"
	-@erase "$(INTDIR)\BondLink.obj"
	-@erase "$(INTDIR)\BondLink.sbr"
	-@erase "$(INTDIR)\CombSort.obj"
	-@erase "$(INTDIR)\CombSort.sbr"
	-@erase "$(INTDIR)\comm.obj"
	-@erase "$(INTDIR)\comm.sbr"
	-@erase "$(INTDIR)\comrgd.obj"
	-@erase "$(INTDIR)\comrgd.sbr"
	-@erase "$(INTDIR)\Consoli.obj"
	-@erase "$(INTDIR)\Consoli.sbr"
	-@erase "$(INTDIR)\Damage 3D.obj"
	-@erase "$(INTDIR)\Damage 3D.sbr"
	-@erase "$(INTDIR)\DM1DForce.obj"
	-@erase "$(INTDIR)\DM1DForce.sbr"
	-@erase "$(INTDIR)\DM3D.OBJ"
	-@erase "$(INTDIR)\DM3D.SBR"
	-@erase "$(INTDIR)\DMINOUT.obj"
	-@erase "$(INTDIR)\DMINOUT.sbr"
	-@erase "$(INTDIR)\DMLOAD1.obj"
	-@erase "$(INTDIR)\DMLOAD1.sbr"
	-@erase "$(INTDIR)\Drucker3D.obj"
	-@erase "$(INTDIR)\Drucker3D.sbr"
	-@erase "$(INTDIR)\dyna.obj"
	-@erase "$(INTDIR)\dyna.sbr"
	-@erase "$(INTDIR)\elem.obj"
	-@erase "$(INTDIR)\elem.sbr"
	-@erase "$(INTDIR)\EXTRUDE.OBJ"
	-@erase "$(INTDIR)\EXTRUDE.SBR"
	-@erase "$(INTDIR)\FMPLASTIC.OBJ"
	-@erase "$(INTDIR)\FMPLASTIC.SBR"
	-@erase "$(INTDIR)\Frame1E.obj"
	-@erase "$(INTDIR)\Frame1E.sbr"
	-@erase "$(INTDIR)\FrameNL.obj"
	-@erase "$(INTDIR)\FrameNL.sbr"
	-@erase "$(INTDIR)\Genfiber.obj"
	-@erase "$(INTDIR)\Genfiber.sbr"
	-@erase "$(INTDIR)\infine.obj"
	-@erase "$(INTDIR)\infine.sbr"
	-@erase "$(INTDIR)\input.obj"
	-@erase "$(INTDIR)\input.sbr"
	-@erase "$(INTDIR)\Interface Element.obj"
	-@erase "$(INTDIR)\Interface Element.sbr"
	-@erase "$(INTDIR)\Jconstraint.obj"
	-@erase "$(INTDIR)\Jconstraint.sbr"
	-@erase "$(INTDIR)\lanc.obj"
	-@erase "$(INTDIR)\lanc.sbr"
	-@erase "$(INTDIR)\LNDSUP.obj"
	-@erase "$(INTDIR)\LNDSUP.sbr"
	-@erase "$(INTDIR)\membrane.obj"
	-@erase "$(INTDIR)\membrane.sbr"
	-@erase "$(INTDIR)\Membrane_Bo.obj"
	-@erase "$(INTDIR)\Membrane_Bo.sbr"
	-@erase "$(INTDIR)\MPSmooth.obj"
	-@erase "$(INTDIR)\MPSmooth.sbr"
	-@erase "$(INTDIR)\NewCable.obj"
	-@erase "$(INTDIR)\NewCable.sbr"
	-@erase "$(INTDIR)\newsoli.obj"
	-@erase "$(INTDIR)\newsoli.sbr"
	-@erase "$(INTDIR)\NewSolver.obj"
	-@erase "$(INTDIR)\NewSolver.sbr"
	-@erase "$(INTDIR)\NLDSUP.obj"
	-@erase "$(INTDIR)\NLDSUP.sbr"
	-@erase "$(INTDIR)\nldyna.obj"
	-@erase "$(INTDIR)\nldyna.sbr"
	-@erase "$(INTDIR)\Nonconso.obj"
	-@erase "$(INTDIR)\Nonconso.sbr"
	-@erase "$(INTDIR)\Nonlinearlink.obj"
	-@erase "$(INTDIR)\Nonlinearlink.sbr"
	-@erase "$(INTDIR)\Password.obj"
	-@erase "$(INTDIR)\Password.sbr"
	-@erase "$(INTDIR)\Plasticity Hard.obj"
	-@erase "$(INTDIR)\Plasticity Hard.sbr"
	-@erase "$(INTDIR)\PLVISCO.OBJ"
	-@erase "$(INTDIR)\PLVISCO.SBR"
	-@erase "$(INTDIR)\RCShell.obj"
	-@erase "$(INTDIR)\RCShell.sbr"
	-@erase "$(INTDIR)\REACTION.OBJ"
	-@erase "$(INTDIR)\REACTION.SBR"
	-@erase "$(INTDIR)\Section prop.obj"
	-@erase "$(INTDIR)\Section prop.sbr"
	-@erase "$(INTDIR)\Seepage.obj"
	-@erase "$(INTDIR)\Seepage.sbr"
	-@erase "$(INTDIR)\SEISMA.obj"
	-@erase "$(INTDIR)\SEISMA.sbr"
	-@erase "$(INTDIR)\SelfweightLoad.obj"
	-@erase "$(INTDIR)\SelfweightLoad.sbr"
	-@erase "$(INTDIR)\Settlement.obj"
	-@erase "$(INTDIR)\Settlement.sbr"
	-@erase "$(INTDIR)\sh3nod.obj"
	-@erase "$(INTDIR)\sh3nod.sbr"
	-@erase "$(INTDIR)\shdrill.obj"
	-@erase "$(INTDIR)\shdrill.sbr"
	-@erase "$(INTDIR)\Shell4 Jacob.obj"
	-@erase "$(INTDIR)\Shell4 Jacob.sbr"
	-@erase "$(INTDIR)\Shell6.obj"
	-@erase "$(INTDIR)\Shell6.sbr"
	-@erase "$(INTDIR)\SOLID.OBJ"
	-@erase "$(INTDIR)\SOLID.SBR"
	-@erase "$(INTDIR)\Solid_Bo.obj"
	-@erase "$(INTDIR)\Solid_Bo.sbr"
	-@erase "$(INTDIR)\Solidaverstr.obj"
	-@erase "$(INTDIR)\Solidaverstr.sbr"
	-@erase "$(INTDIR)\SOLIDECONSO.OBJ"
	-@erase "$(INTDIR)\SOLIDECONSO.SBR"
	-@erase "$(INTDIR)\SOLIDHBD3.OBJ"
	-@erase "$(INTDIR)\SOLIDHBD3.SBR"
	-@erase "$(INTDIR)\SPConcrete.obj"
	-@erase "$(INTDIR)\SPConcrete.sbr"
	-@erase "$(INTDIR)\Spring1.obj"
	-@erase "$(INTDIR)\Spring1.sbr"
	-@erase "$(INTDIR)\TaperSec.obj"
	-@erase "$(INTDIR)\TaperSec.sbr"
	-@erase "$(INTDIR)\TendonForce.obj"
	-@erase "$(INTDIR)\TendonForce.sbr"
	-@erase "$(INTDIR)\Thin3D.obj"
	-@erase "$(INTDIR)\Thin3D.sbr"
	-@erase "$(INTDIR)\ThinSolid.obj"
	-@erase "$(INTDIR)\ThinSolid.sbr"
	-@erase "$(INTDIR)\truss.obj"
	-@erase "$(INTDIR)\truss.sbr"
	-@erase "$(INTDIR)\Warray.obj"
	-@erase "$(INTDIR)\Warray.sbr"
	-@erase "$(INTDIR)\xfinas.obj"
	-@erase "$(INTDIR)\xfinas.sbr"
	-@erase "$(INTDIR)\XPCCHANGE.OBJ"
	-@erase "$(INTDIR)\XPCCHANGE.SBR"
	-@erase "$(INTDIR)\XPCCONCRETE.OBJ"
	-@erase "$(INTDIR)\XPCCONCRETE.SBR"
	-@erase "$(INTDIR)\XPCCSTATE.OBJ"
	-@erase "$(INTDIR)\XPCCSTATE.SBR"
	-@erase "$(INTDIR)\XPCINPUT.OBJ"
	-@erase "$(INTDIR)\XPCINPUT.SBR"
	-@erase "$(INTDIR)\XPCLOAD.OBJ"
	-@erase "$(INTDIR)\XPCLOAD.SBR"
	-@erase "$(INTDIR)\XPCNASCAB.OBJ"
	-@erase "$(INTDIR)\XPCNASCAB.SBR"
	-@erase "$(INTDIR)\XPCOUTPUT.OBJ"
	-@erase "$(INTDIR)\XPCOUTPUT.SBR"
	-@erase "$(INTDIR)\XPCSOLVE.OBJ"
	-@erase "$(INTDIR)\XPCSOLVE.SBR"
	-@erase "$(INTDIR)\XPCSTATE.OBJ"
	-@erase "$(INTDIR)\XPCSTATE.SBR"
	-@erase "$(INTDIR)\XPCSTIFF.OBJ"
	-@erase "$(INTDIR)\XPCSTIFF.SBR"
	-@erase "$(INTDIR)\XPCSTYILD.OBJ"
	-@erase "$(INTDIR)\XPCSTYILD.SBR"
	-@erase "$(INTDIR)\XPCSTYINP.OBJ"
	-@erase "$(INTDIR)\XPCSTYINP.SBR"
	-@erase "$(INTDIR)\XPCUTILITY.OBJ"
	-@erase "$(INTDIR)\XPCUTILITY.SBR"
	-@erase "$(OUTDIR)\xfinas.bsc"
	-@erase "$(OUTDIR)\xfinas.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=df.exe
F90_PROJ=/browser:"Release/" /include:"$(INTDIR)\\" /compile_only /nologo\
 /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\xfinas.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\2D Heat main.sbr" \
	"$(INTDIR)\3D Heat main.sbr" \
	"$(INTDIR)\ADS1.SBR" \
	"$(INTDIR)\ADS2.SBR" \
	"$(INTDIR)\ADS3.SBR" \
	"$(INTDIR)\baseiso.sbr" \
	"$(INTDIR)\beam close section.sbr" \
	"$(INTDIR)\beam21f.sbr" \
	"$(INTDIR)\BOND3D.SBR" \
	"$(INTDIR)\BondLink.sbr" \
	"$(INTDIR)\CombSort.sbr" \
	"$(INTDIR)\comm.sbr" \
	"$(INTDIR)\comrgd.sbr" \
	"$(INTDIR)\Consoli.sbr" \
	"$(INTDIR)\Damage 3D.sbr" \
	"$(INTDIR)\DM1DForce.sbr" \
	"$(INTDIR)\DM3D.SBR" \
	"$(INTDIR)\DMINOUT.sbr" \
	"$(INTDIR)\DMLOAD1.sbr" \
	"$(INTDIR)\Drucker3D.sbr" \
	"$(INTDIR)\dyna.sbr" \
	"$(INTDIR)\elem.sbr" \
	"$(INTDIR)\EXTRUDE.SBR" \
	"$(INTDIR)\FMPLASTIC.SBR" \
	"$(INTDIR)\Frame1E.sbr" \
	"$(INTDIR)\FrameNL.sbr" \
	"$(INTDIR)\Genfiber.sbr" \
	"$(INTDIR)\infine.sbr" \
	"$(INTDIR)\input.sbr" \
	"$(INTDIR)\Interface Element.sbr" \
	"$(INTDIR)\Jconstraint.sbr" \
	"$(INTDIR)\lanc.sbr" \
	"$(INTDIR)\LNDSUP.sbr" \
	"$(INTDIR)\membrane.sbr" \
	"$(INTDIR)\Membrane_Bo.sbr" \
	"$(INTDIR)\MPSmooth.sbr" \
	"$(INTDIR)\NewCable.sbr" \
	"$(INTDIR)\newsoli.sbr" \
	"$(INTDIR)\NewSolver.sbr" \
	"$(INTDIR)\NLDSUP.sbr" \
	"$(INTDIR)\nldyna.sbr" \
	"$(INTDIR)\Nonconso.sbr" \
	"$(INTDIR)\Nonlinearlink.sbr" \
	"$(INTDIR)\Password.sbr" \
	"$(INTDIR)\Plasticity Hard.sbr" \
	"$(INTDIR)\PLVISCO.SBR" \
	"$(INTDIR)\RCShell.sbr" \
	"$(INTDIR)\REACTION.SBR" \
	"$(INTDIR)\Section prop.sbr" \
	"$(INTDIR)\Seepage.sbr" \
	"$(INTDIR)\SEISMA.sbr" \
	"$(INTDIR)\SelfweightLoad.sbr" \
	"$(INTDIR)\Settlement.sbr" \
	"$(INTDIR)\sh3nod.sbr" \
	"$(INTDIR)\shdrill.sbr" \
	"$(INTDIR)\Shell4 Jacob.sbr" \
	"$(INTDIR)\Shell6.sbr" \
	"$(INTDIR)\SOLID.SBR" \
	"$(INTDIR)\Solid_Bo.sbr" \
	"$(INTDIR)\Solidaverstr.sbr" \
	"$(INTDIR)\SOLIDECONSO.SBR" \
	"$(INTDIR)\SOLIDHBD3.SBR" \
	"$(INTDIR)\SPConcrete.sbr" \
	"$(INTDIR)\Spring1.sbr" \
	"$(INTDIR)\TaperSec.sbr" \
	"$(INTDIR)\TendonForce.sbr" \
	"$(INTDIR)\Thin3D.sbr" \
	"$(INTDIR)\ThinSolid.sbr" \
	"$(INTDIR)\truss.sbr" \
	"$(INTDIR)\Warray.sbr" \
	"$(INTDIR)\xfinas.sbr" \
	"$(INTDIR)\XPCCHANGE.SBR" \
	"$(INTDIR)\XPCCONCRETE.SBR" \
	"$(INTDIR)\XPCCSTATE.SBR" \
	"$(INTDIR)\XPCINPUT.SBR" \
	"$(INTDIR)\XPCLOAD.SBR" \
	"$(INTDIR)\XPCNASCAB.SBR" \
	"$(INTDIR)\XPCOUTPUT.SBR" \
	"$(INTDIR)\XPCSOLVE.SBR" \
	"$(INTDIR)\XPCSTATE.SBR" \
	"$(INTDIR)\XPCSTIFF.SBR" \
	"$(INTDIR)\XPCSTYILD.SBR" \
	"$(INTDIR)\XPCSTYINP.SBR" \
	"$(INTDIR)\XPCUTILITY.SBR"

"$(OUTDIR)\xfinas.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /stack:0x35a4e900,0x35a4e900\
 /subsystem:console /incremental:no /pdb:"$(OUTDIR)\xfinas.pdb" /machine:I386\
 /out:"$(OUTDIR)\xfinas.exe" 
LINK32_OBJS= \
	"$(INTDIR)\2D Heat main.obj" \
	"$(INTDIR)\3D Heat main.obj" \
	"$(INTDIR)\ADS1.OBJ" \
	"$(INTDIR)\ADS2.OBJ" \
	"$(INTDIR)\ADS3.OBJ" \
	"$(INTDIR)\baseiso.obj" \
	"$(INTDIR)\beam close section.obj" \
	"$(INTDIR)\beam21f.obj" \
	"$(INTDIR)\BOND3D.OBJ" \
	"$(INTDIR)\BondLink.obj" \
	"$(INTDIR)\CombSort.obj" \
	"$(INTDIR)\comm.obj" \
	"$(INTDIR)\comrgd.obj" \
	"$(INTDIR)\Consoli.obj" \
	"$(INTDIR)\Damage 3D.obj" \
	"$(INTDIR)\DM1DForce.obj" \
	"$(INTDIR)\DM3D.OBJ" \
	"$(INTDIR)\DMINOUT.obj" \
	"$(INTDIR)\DMLOAD1.obj" \
	"$(INTDIR)\Drucker3D.obj" \
	"$(INTDIR)\dyna.obj" \
	"$(INTDIR)\elem.obj" \
	"$(INTDIR)\EXTRUDE.OBJ" \
	"$(INTDIR)\FMPLASTIC.OBJ" \
	"$(INTDIR)\Frame1E.obj" \
	"$(INTDIR)\FrameNL.obj" \
	"$(INTDIR)\Genfiber.obj" \
	"$(INTDIR)\infine.obj" \
	"$(INTDIR)\input.obj" \
	"$(INTDIR)\Interface Element.obj" \
	"$(INTDIR)\Jconstraint.obj" \
	"$(INTDIR)\lanc.obj" \
	"$(INTDIR)\LNDSUP.obj" \
	"$(INTDIR)\membrane.obj" \
	"$(INTDIR)\Membrane_Bo.obj" \
	"$(INTDIR)\MPSmooth.obj" \
	"$(INTDIR)\NewCable.obj" \
	"$(INTDIR)\newsoli.obj" \
	"$(INTDIR)\NewSolver.obj" \
	"$(INTDIR)\NLDSUP.obj" \
	"$(INTDIR)\nldyna.obj" \
	"$(INTDIR)\Nonconso.obj" \
	"$(INTDIR)\Nonlinearlink.obj" \
	"$(INTDIR)\Password.obj" \
	"$(INTDIR)\Plasticity Hard.obj" \
	"$(INTDIR)\PLVISCO.OBJ" \
	"$(INTDIR)\RCShell.obj" \
	"$(INTDIR)\REACTION.OBJ" \
	"$(INTDIR)\Section prop.obj" \
	"$(INTDIR)\Seepage.obj" \
	"$(INTDIR)\SEISMA.obj" \
	"$(INTDIR)\SelfweightLoad.obj" \
	"$(INTDIR)\Settlement.obj" \
	"$(INTDIR)\sh3nod.obj" \
	"$(INTDIR)\shdrill.obj" \
	"$(INTDIR)\Shell4 Jacob.obj" \
	"$(INTDIR)\Shell6.obj" \
	"$(INTDIR)\SOLID.OBJ" \
	"$(INTDIR)\Solid_Bo.obj" \
	"$(INTDIR)\Solidaverstr.obj" \
	"$(INTDIR)\SOLIDECONSO.OBJ" \
	"$(INTDIR)\SOLIDHBD3.OBJ" \
	"$(INTDIR)\SPConcrete.obj" \
	"$(INTDIR)\Spring1.obj" \
	"$(INTDIR)\TaperSec.obj" \
	"$(INTDIR)\TendonForce.obj" \
	"$(INTDIR)\Thin3D.obj" \
	"$(INTDIR)\ThinSolid.obj" \
	"$(INTDIR)\truss.obj" \
	"$(INTDIR)\Warray.obj" \
	"$(INTDIR)\xfinas.obj" \
	"$(INTDIR)\XPCCHANGE.OBJ" \
	"$(INTDIR)\XPCCONCRETE.OBJ" \
	"$(INTDIR)\XPCCSTATE.OBJ" \
	"$(INTDIR)\XPCINPUT.OBJ" \
	"$(INTDIR)\XPCLOAD.OBJ" \
	"$(INTDIR)\XPCNASCAB.OBJ" \
	"$(INTDIR)\XPCOUTPUT.OBJ" \
	"$(INTDIR)\XPCSOLVE.OBJ" \
	"$(INTDIR)\XPCSTATE.OBJ" \
	"$(INTDIR)\XPCSTIFF.OBJ" \
	"$(INTDIR)\XPCSTYILD.OBJ" \
	"$(INTDIR)\XPCSTYINP.OBJ" \
	"$(INTDIR)\XPCUTILITY.OBJ"

"$(OUTDIR)\xfinas.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "xfinas - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\xfinas.exe" "$(OUTDIR)\DF50.PDB" "$(OUTDIR)\xfinas.bsc"

!ELSE 

ALL : "$(OUTDIR)\xfinas.exe" "$(OUTDIR)\DF50.PDB" "$(OUTDIR)\xfinas.bsc"

!ENDIF 

CLEAN :
	-@erase "$(INTDIR)\2D Heat main.obj"
	-@erase "$(INTDIR)\2D Heat main.sbr"
	-@erase "$(INTDIR)\3D Heat main.obj"
	-@erase "$(INTDIR)\3D Heat main.sbr"
	-@erase "$(INTDIR)\ADS1.OBJ"
	-@erase "$(INTDIR)\ADS1.SBR"
	-@erase "$(INTDIR)\ADS2.OBJ"
	-@erase "$(INTDIR)\ADS2.SBR"
	-@erase "$(INTDIR)\ADS3.OBJ"
	-@erase "$(INTDIR)\ADS3.SBR"
	-@erase "$(INTDIR)\baseiso.obj"
	-@erase "$(INTDIR)\baseiso.sbr"
	-@erase "$(INTDIR)\beam close section.obj"
	-@erase "$(INTDIR)\beam close section.sbr"
	-@erase "$(INTDIR)\beam21f.obj"
	-@erase "$(INTDIR)\beam21f.sbr"
	-@erase "$(INTDIR)\BOND3D.OBJ"
	-@erase "$(INTDIR)\BOND3D.SBR"
	-@erase "$(INTDIR)\BondLink.obj"
	-@erase "$(INTDIR)\BondLink.sbr"
	-@erase "$(INTDIR)\CombSort.obj"
	-@erase "$(INTDIR)\CombSort.sbr"
	-@erase "$(INTDIR)\comm.obj"
	-@erase "$(INTDIR)\comm.sbr"
	-@erase "$(INTDIR)\comrgd.obj"
	-@erase "$(INTDIR)\comrgd.sbr"
	-@erase "$(INTDIR)\Consoli.obj"
	-@erase "$(INTDIR)\Consoli.sbr"
	-@erase "$(INTDIR)\Damage 3D.obj"
	-@erase "$(INTDIR)\Damage 3D.sbr"
	-@erase "$(INTDIR)\DF50.PDB"
	-@erase "$(INTDIR)\DM1DForce.obj"
	-@erase "$(INTDIR)\DM1DForce.sbr"
	-@erase "$(INTDIR)\DM3D.OBJ"
	-@erase "$(INTDIR)\DM3D.SBR"
	-@erase "$(INTDIR)\DMINOUT.obj"
	-@erase "$(INTDIR)\DMINOUT.sbr"
	-@erase "$(INTDIR)\DMLOAD1.obj"
	-@erase "$(INTDIR)\DMLOAD1.sbr"
	-@erase "$(INTDIR)\Drucker3D.obj"
	-@erase "$(INTDIR)\Drucker3D.sbr"
	-@erase "$(INTDIR)\dyna.obj"
	-@erase "$(INTDIR)\dyna.sbr"
	-@erase "$(INTDIR)\elem.obj"
	-@erase "$(INTDIR)\elem.sbr"
	-@erase "$(INTDIR)\EXTRUDE.OBJ"
	-@erase "$(INTDIR)\EXTRUDE.SBR"
	-@erase "$(INTDIR)\FMPLASTIC.OBJ"
	-@erase "$(INTDIR)\FMPLASTIC.SBR"
	-@erase "$(INTDIR)\Frame1E.obj"
	-@erase "$(INTDIR)\Frame1E.sbr"
	-@erase "$(INTDIR)\FrameNL.obj"
	-@erase "$(INTDIR)\FrameNL.sbr"
	-@erase "$(INTDIR)\Genfiber.obj"
	-@erase "$(INTDIR)\Genfiber.sbr"
	-@erase "$(INTDIR)\infine.obj"
	-@erase "$(INTDIR)\infine.sbr"
	-@erase "$(INTDIR)\input.obj"
	-@erase "$(INTDIR)\input.sbr"
	-@erase "$(INTDIR)\Interface Element.obj"
	-@erase "$(INTDIR)\Interface Element.sbr"
	-@erase "$(INTDIR)\Jconstraint.obj"
	-@erase "$(INTDIR)\Jconstraint.sbr"
	-@erase "$(INTDIR)\lanc.obj"
	-@erase "$(INTDIR)\lanc.sbr"
	-@erase "$(INTDIR)\LNDSUP.obj"
	-@erase "$(INTDIR)\LNDSUP.sbr"
	-@erase "$(INTDIR)\membrane.obj"
	-@erase "$(INTDIR)\membrane.sbr"
	-@erase "$(INTDIR)\Membrane_Bo.obj"
	-@erase "$(INTDIR)\Membrane_Bo.sbr"
	-@erase "$(INTDIR)\MPSmooth.obj"
	-@erase "$(INTDIR)\MPSmooth.sbr"
	-@erase "$(INTDIR)\NewCable.obj"
	-@erase "$(INTDIR)\NewCable.sbr"
	-@erase "$(INTDIR)\newsoli.obj"
	-@erase "$(INTDIR)\newsoli.sbr"
	-@erase "$(INTDIR)\NewSolver.obj"
	-@erase "$(INTDIR)\NewSolver.sbr"
	-@erase "$(INTDIR)\NLDSUP.obj"
	-@erase "$(INTDIR)\NLDSUP.sbr"
	-@erase "$(INTDIR)\nldyna.obj"
	-@erase "$(INTDIR)\nldyna.sbr"
	-@erase "$(INTDIR)\Nonconso.obj"
	-@erase "$(INTDIR)\Nonconso.sbr"
	-@erase "$(INTDIR)\Nonlinearlink.obj"
	-@erase "$(INTDIR)\Nonlinearlink.sbr"
	-@erase "$(INTDIR)\Password.obj"
	-@erase "$(INTDIR)\Password.sbr"
	-@erase "$(INTDIR)\Plasticity Hard.obj"
	-@erase "$(INTDIR)\Plasticity Hard.sbr"
	-@erase "$(INTDIR)\PLVISCO.OBJ"
	-@erase "$(INTDIR)\PLVISCO.SBR"
	-@erase "$(INTDIR)\RCShell.obj"
	-@erase "$(INTDIR)\RCShell.sbr"
	-@erase "$(INTDIR)\REACTION.OBJ"
	-@erase "$(INTDIR)\REACTION.SBR"
	-@erase "$(INTDIR)\Section prop.obj"
	-@erase "$(INTDIR)\Section prop.sbr"
	-@erase "$(INTDIR)\Seepage.obj"
	-@erase "$(INTDIR)\Seepage.sbr"
	-@erase "$(INTDIR)\SEISMA.obj"
	-@erase "$(INTDIR)\SEISMA.sbr"
	-@erase "$(INTDIR)\SelfweightLoad.obj"
	-@erase "$(INTDIR)\SelfweightLoad.sbr"
	-@erase "$(INTDIR)\Settlement.obj"
	-@erase "$(INTDIR)\Settlement.sbr"
	-@erase "$(INTDIR)\sh3nod.obj"
	-@erase "$(INTDIR)\sh3nod.sbr"
	-@erase "$(INTDIR)\shdrill.obj"
	-@erase "$(INTDIR)\shdrill.sbr"
	-@erase "$(INTDIR)\Shell4 Jacob.obj"
	-@erase "$(INTDIR)\Shell4 Jacob.sbr"
	-@erase "$(INTDIR)\Shell6.obj"
	-@erase "$(INTDIR)\Shell6.sbr"
	-@erase "$(INTDIR)\SOLID.OBJ"
	-@erase "$(INTDIR)\SOLID.SBR"
	-@erase "$(INTDIR)\Solid_Bo.obj"
	-@erase "$(INTDIR)\Solid_Bo.sbr"
	-@erase "$(INTDIR)\Solidaverstr.obj"
	-@erase "$(INTDIR)\Solidaverstr.sbr"
	-@erase "$(INTDIR)\SOLIDECONSO.OBJ"
	-@erase "$(INTDIR)\SOLIDECONSO.SBR"
	-@erase "$(INTDIR)\SOLIDHBD3.OBJ"
	-@erase "$(INTDIR)\SOLIDHBD3.SBR"
	-@erase "$(INTDIR)\SPConcrete.obj"
	-@erase "$(INTDIR)\SPConcrete.sbr"
	-@erase "$(INTDIR)\Spring1.obj"
	-@erase "$(INTDIR)\Spring1.sbr"
	-@erase "$(INTDIR)\TaperSec.obj"
	-@erase "$(INTDIR)\TaperSec.sbr"
	-@erase "$(INTDIR)\TendonForce.obj"
	-@erase "$(INTDIR)\TendonForce.sbr"
	-@erase "$(INTDIR)\Thin3D.obj"
	-@erase "$(INTDIR)\Thin3D.sbr"
	-@erase "$(INTDIR)\ThinSolid.obj"
	-@erase "$(INTDIR)\ThinSolid.sbr"
	-@erase "$(INTDIR)\truss.obj"
	-@erase "$(INTDIR)\truss.sbr"
	-@erase "$(INTDIR)\Warray.obj"
	-@erase "$(INTDIR)\Warray.sbr"
	-@erase "$(INTDIR)\xfinas.obj"
	-@erase "$(INTDIR)\xfinas.sbr"
	-@erase "$(INTDIR)\XPCCHANGE.OBJ"
	-@erase "$(INTDIR)\XPCCHANGE.SBR"
	-@erase "$(INTDIR)\XPCCONCRETE.OBJ"
	-@erase "$(INTDIR)\XPCCONCRETE.SBR"
	-@erase "$(INTDIR)\XPCCSTATE.OBJ"
	-@erase "$(INTDIR)\XPCCSTATE.SBR"
	-@erase "$(INTDIR)\XPCINPUT.OBJ"
	-@erase "$(INTDIR)\XPCINPUT.SBR"
	-@erase "$(INTDIR)\XPCLOAD.OBJ"
	-@erase "$(INTDIR)\XPCLOAD.SBR"
	-@erase "$(INTDIR)\XPCNASCAB.OBJ"
	-@erase "$(INTDIR)\XPCNASCAB.SBR"
	-@erase "$(INTDIR)\XPCOUTPUT.OBJ"
	-@erase "$(INTDIR)\XPCOUTPUT.SBR"
	-@erase "$(INTDIR)\XPCSOLVE.OBJ"
	-@erase "$(INTDIR)\XPCSOLVE.SBR"
	-@erase "$(INTDIR)\XPCSTATE.OBJ"
	-@erase "$(INTDIR)\XPCSTATE.SBR"
	-@erase "$(INTDIR)\XPCSTIFF.OBJ"
	-@erase "$(INTDIR)\XPCSTIFF.SBR"
	-@erase "$(INTDIR)\XPCSTYILD.OBJ"
	-@erase "$(INTDIR)\XPCSTYILD.SBR"
	-@erase "$(INTDIR)\XPCSTYINP.OBJ"
	-@erase "$(INTDIR)\XPCSTYINP.SBR"
	-@erase "$(INTDIR)\XPCUTILITY.OBJ"
	-@erase "$(INTDIR)\XPCUTILITY.SBR"
	-@erase "$(OUTDIR)\xfinas.bsc"
	-@erase "$(OUTDIR)\xfinas.exe"
	-@erase "$(OUTDIR)\xfinas.ilk"
	-@erase "$(OUTDIR)\xfinas.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=df.exe
F90_PROJ=/browser:"Debug/" /include:"$(INTDIR)\\" /compile_only /nologo\
 /debug:full /warn:nofileopt /module:"Debug/" /object:"Debug/"\
 /pdbfile:"Debug/DF50.PDB" 
F90_OBJS=.\Debug/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\xfinas.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\2D Heat main.sbr" \
	"$(INTDIR)\3D Heat main.sbr" \
	"$(INTDIR)\ADS1.SBR" \
	"$(INTDIR)\ADS2.SBR" \
	"$(INTDIR)\ADS3.SBR" \
	"$(INTDIR)\baseiso.sbr" \
	"$(INTDIR)\beam close section.sbr" \
	"$(INTDIR)\beam21f.sbr" \
	"$(INTDIR)\BOND3D.SBR" \
	"$(INTDIR)\BondLink.sbr" \
	"$(INTDIR)\CombSort.sbr" \
	"$(INTDIR)\comm.sbr" \
	"$(INTDIR)\comrgd.sbr" \
	"$(INTDIR)\Consoli.sbr" \
	"$(INTDIR)\Damage 3D.sbr" \
	"$(INTDIR)\DM1DForce.sbr" \
	"$(INTDIR)\DM3D.SBR" \
	"$(INTDIR)\DMINOUT.sbr" \
	"$(INTDIR)\DMLOAD1.sbr" \
	"$(INTDIR)\Drucker3D.sbr" \
	"$(INTDIR)\dyna.sbr" \
	"$(INTDIR)\elem.sbr" \
	"$(INTDIR)\EXTRUDE.SBR" \
	"$(INTDIR)\FMPLASTIC.SBR" \
	"$(INTDIR)\Frame1E.sbr" \
	"$(INTDIR)\FrameNL.sbr" \
	"$(INTDIR)\Genfiber.sbr" \
	"$(INTDIR)\infine.sbr" \
	"$(INTDIR)\input.sbr" \
	"$(INTDIR)\Interface Element.sbr" \
	"$(INTDIR)\Jconstraint.sbr" \
	"$(INTDIR)\lanc.sbr" \
	"$(INTDIR)\LNDSUP.sbr" \
	"$(INTDIR)\membrane.sbr" \
	"$(INTDIR)\Membrane_Bo.sbr" \
	"$(INTDIR)\MPSmooth.sbr" \
	"$(INTDIR)\NewCable.sbr" \
	"$(INTDIR)\newsoli.sbr" \
	"$(INTDIR)\NewSolver.sbr" \
	"$(INTDIR)\NLDSUP.sbr" \
	"$(INTDIR)\nldyna.sbr" \
	"$(INTDIR)\Nonconso.sbr" \
	"$(INTDIR)\Nonlinearlink.sbr" \
	"$(INTDIR)\Password.sbr" \
	"$(INTDIR)\Plasticity Hard.sbr" \
	"$(INTDIR)\PLVISCO.SBR" \
	"$(INTDIR)\RCShell.sbr" \
	"$(INTDIR)\REACTION.SBR" \
	"$(INTDIR)\Section prop.sbr" \
	"$(INTDIR)\Seepage.sbr" \
	"$(INTDIR)\SEISMA.sbr" \
	"$(INTDIR)\SelfweightLoad.sbr" \
	"$(INTDIR)\Settlement.sbr" \
	"$(INTDIR)\sh3nod.sbr" \
	"$(INTDIR)\shdrill.sbr" \
	"$(INTDIR)\Shell4 Jacob.sbr" \
	"$(INTDIR)\Shell6.sbr" \
	"$(INTDIR)\SOLID.SBR" \
	"$(INTDIR)\Solid_Bo.sbr" \
	"$(INTDIR)\Solidaverstr.sbr" \
	"$(INTDIR)\SOLIDECONSO.SBR" \
	"$(INTDIR)\SOLIDHBD3.SBR" \
	"$(INTDIR)\SPConcrete.sbr" \
	"$(INTDIR)\Spring1.sbr" \
	"$(INTDIR)\TaperSec.sbr" \
	"$(INTDIR)\TendonForce.sbr" \
	"$(INTDIR)\Thin3D.sbr" \
	"$(INTDIR)\ThinSolid.sbr" \
	"$(INTDIR)\truss.sbr" \
	"$(INTDIR)\Warray.sbr" \
	"$(INTDIR)\xfinas.sbr" \
	"$(INTDIR)\XPCCHANGE.SBR" \
	"$(INTDIR)\XPCCONCRETE.SBR" \
	"$(INTDIR)\XPCCSTATE.SBR" \
	"$(INTDIR)\XPCINPUT.SBR" \
	"$(INTDIR)\XPCLOAD.SBR" \
	"$(INTDIR)\XPCNASCAB.SBR" \
	"$(INTDIR)\XPCOUTPUT.SBR" \
	"$(INTDIR)\XPCSOLVE.SBR" \
	"$(INTDIR)\XPCSTATE.SBR" \
	"$(INTDIR)\XPCSTIFF.SBR" \
	"$(INTDIR)\XPCSTYILD.SBR" \
	"$(INTDIR)\XPCSTYINP.SBR" \
	"$(INTDIR)\XPCUTILITY.SBR"

"$(OUTDIR)\xfinas.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /stack:0x5f5e100,0x5f5e100 /subsystem:console\
 /incremental:yes /pdb:"$(OUTDIR)\xfinas.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)\xfinas.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\2D Heat main.obj" \
	"$(INTDIR)\3D Heat main.obj" \
	"$(INTDIR)\ADS1.OBJ" \
	"$(INTDIR)\ADS2.OBJ" \
	"$(INTDIR)\ADS3.OBJ" \
	"$(INTDIR)\baseiso.obj" \
	"$(INTDIR)\beam close section.obj" \
	"$(INTDIR)\beam21f.obj" \
	"$(INTDIR)\BOND3D.OBJ" \
	"$(INTDIR)\BondLink.obj" \
	"$(INTDIR)\CombSort.obj" \
	"$(INTDIR)\comm.obj" \
	"$(INTDIR)\comrgd.obj" \
	"$(INTDIR)\Consoli.obj" \
	"$(INTDIR)\Damage 3D.obj" \
	"$(INTDIR)\DM1DForce.obj" \
	"$(INTDIR)\DM3D.OBJ" \
	"$(INTDIR)\DMINOUT.obj" \
	"$(INTDIR)\DMLOAD1.obj" \
	"$(INTDIR)\Drucker3D.obj" \
	"$(INTDIR)\dyna.obj" \
	"$(INTDIR)\elem.obj" \
	"$(INTDIR)\EXTRUDE.OBJ" \
	"$(INTDIR)\FMPLASTIC.OBJ" \
	"$(INTDIR)\Frame1E.obj" \
	"$(INTDIR)\FrameNL.obj" \
	"$(INTDIR)\Genfiber.obj" \
	"$(INTDIR)\infine.obj" \
	"$(INTDIR)\input.obj" \
	"$(INTDIR)\Interface Element.obj" \
	"$(INTDIR)\Jconstraint.obj" \
	"$(INTDIR)\lanc.obj" \
	"$(INTDIR)\LNDSUP.obj" \
	"$(INTDIR)\membrane.obj" \
	"$(INTDIR)\Membrane_Bo.obj" \
	"$(INTDIR)\MPSmooth.obj" \
	"$(INTDIR)\NewCable.obj" \
	"$(INTDIR)\newsoli.obj" \
	"$(INTDIR)\NewSolver.obj" \
	"$(INTDIR)\NLDSUP.obj" \
	"$(INTDIR)\nldyna.obj" \
	"$(INTDIR)\Nonconso.obj" \
	"$(INTDIR)\Nonlinearlink.obj" \
	"$(INTDIR)\Password.obj" \
	"$(INTDIR)\Plasticity Hard.obj" \
	"$(INTDIR)\PLVISCO.OBJ" \
	"$(INTDIR)\RCShell.obj" \
	"$(INTDIR)\REACTION.OBJ" \
	"$(INTDIR)\Section prop.obj" \
	"$(INTDIR)\Seepage.obj" \
	"$(INTDIR)\SEISMA.obj" \
	"$(INTDIR)\SelfweightLoad.obj" \
	"$(INTDIR)\Settlement.obj" \
	"$(INTDIR)\sh3nod.obj" \
	"$(INTDIR)\shdrill.obj" \
	"$(INTDIR)\Shell4 Jacob.obj" \
	"$(INTDIR)\Shell6.obj" \
	"$(INTDIR)\SOLID.OBJ" \
	"$(INTDIR)\Solid_Bo.obj" \
	"$(INTDIR)\Solidaverstr.obj" \
	"$(INTDIR)\SOLIDECONSO.OBJ" \
	"$(INTDIR)\SOLIDHBD3.OBJ" \
	"$(INTDIR)\SPConcrete.obj" \
	"$(INTDIR)\Spring1.obj" \
	"$(INTDIR)\TaperSec.obj" \
	"$(INTDIR)\TendonForce.obj" \
	"$(INTDIR)\Thin3D.obj" \
	"$(INTDIR)\ThinSolid.obj" \
	"$(INTDIR)\truss.obj" \
	"$(INTDIR)\Warray.obj" \
	"$(INTDIR)\xfinas.obj" \
	"$(INTDIR)\XPCCHANGE.OBJ" \
	"$(INTDIR)\XPCCONCRETE.OBJ" \
	"$(INTDIR)\XPCCSTATE.OBJ" \
	"$(INTDIR)\XPCINPUT.OBJ" \
	"$(INTDIR)\XPCLOAD.OBJ" \
	"$(INTDIR)\XPCNASCAB.OBJ" \
	"$(INTDIR)\XPCOUTPUT.OBJ" \
	"$(INTDIR)\XPCSOLVE.OBJ" \
	"$(INTDIR)\XPCSTATE.OBJ" \
	"$(INTDIR)\XPCSTIFF.OBJ" \
	"$(INTDIR)\XPCSTYILD.OBJ" \
	"$(INTDIR)\XPCSTYINP.OBJ" \
	"$(INTDIR)\XPCUTILITY.OBJ"

"$(OUTDIR)\xfinas.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(CFG)" == "xfinas - Win32 Release" || "$(CFG)" == "xfinas - Win32 Debug"
SOURCE=".\2D Heat main.for"

"$(INTDIR)\2D Heat main.obj"	"$(INTDIR)\2D Heat main.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=".\3D Heat main.for"

"$(INTDIR)\3D Heat main.obj"	"$(INTDIR)\3D Heat main.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\ADS1.FOR

"$(INTDIR)\ADS1.OBJ"	"$(INTDIR)\ADS1.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\ADS2.FOR

"$(INTDIR)\ADS2.OBJ"	"$(INTDIR)\ADS2.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\ADS3.FOR

"$(INTDIR)\ADS3.OBJ"	"$(INTDIR)\ADS3.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\baseiso.for

"$(INTDIR)\baseiso.obj"	"$(INTDIR)\baseiso.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=".\beam close section.for"

"$(INTDIR)\beam close section.obj"	"$(INTDIR)\beam close section.sbr" : \
$(SOURCE) "$(INTDIR)"


SOURCE=.\beam21f.for

"$(INTDIR)\beam21f.obj"	"$(INTDIR)\beam21f.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\BOND3D.FOR

"$(INTDIR)\BOND3D.OBJ"	"$(INTDIR)\BOND3D.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\BondLink.for

"$(INTDIR)\BondLink.obj"	"$(INTDIR)\BondLink.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\CombSort.for

"$(INTDIR)\CombSort.obj"	"$(INTDIR)\CombSort.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\comm.for

"$(INTDIR)\comm.obj"	"$(INTDIR)\comm.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\comrgd.for

"$(INTDIR)\comrgd.obj"	"$(INTDIR)\comrgd.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Consoli.for

"$(INTDIR)\Consoli.obj"	"$(INTDIR)\Consoli.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=".\Damage 3D.for"

"$(INTDIR)\Damage 3D.obj"	"$(INTDIR)\Damage 3D.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\DM1DForce.FOR

"$(INTDIR)\DM1DForce.obj"	"$(INTDIR)\DM1DForce.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\DM3D.FOR

"$(INTDIR)\DM3D.OBJ"	"$(INTDIR)\DM3D.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\DMINOUT.for

"$(INTDIR)\DMINOUT.obj"	"$(INTDIR)\DMINOUT.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\DMLOAD1.for

"$(INTDIR)\DMLOAD1.obj"	"$(INTDIR)\DMLOAD1.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Drucker3D.for

"$(INTDIR)\Drucker3D.obj"	"$(INTDIR)\Drucker3D.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\dyna.for

"$(INTDIR)\dyna.obj"	"$(INTDIR)\dyna.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\elem.for

"$(INTDIR)\elem.obj"	"$(INTDIR)\elem.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\EXTRUDE.FOR

"$(INTDIR)\EXTRUDE.OBJ"	"$(INTDIR)\EXTRUDE.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\FMPLASTIC.FOR

"$(INTDIR)\FMPLASTIC.OBJ"	"$(INTDIR)\FMPLASTIC.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Frame1E.for

"$(INTDIR)\Frame1E.obj"	"$(INTDIR)\Frame1E.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\FrameNL.for

"$(INTDIR)\FrameNL.obj"	"$(INTDIR)\FrameNL.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Genfiber.for

"$(INTDIR)\Genfiber.obj"	"$(INTDIR)\Genfiber.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\infine.for

"$(INTDIR)\infine.obj"	"$(INTDIR)\infine.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\input.for

"$(INTDIR)\input.obj"	"$(INTDIR)\input.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=".\Interface Element.for"

"$(INTDIR)\Interface Element.obj"	"$(INTDIR)\Interface Element.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\Jconstraint.for

"$(INTDIR)\Jconstraint.obj"	"$(INTDIR)\Jconstraint.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\lanc.for
DEP_F90_LANC_=\
	"..\..\..\..\..\..\Program Files\DevStudio\DF\IMSL\Include\MSIMSL.mod"\
	

"$(INTDIR)\lanc.obj"	"$(INTDIR)\lanc.sbr" : $(SOURCE) $(DEP_F90_LANC_)\
 "$(INTDIR)"


SOURCE=.\LNDSUP.for

"$(INTDIR)\LNDSUP.obj"	"$(INTDIR)\LNDSUP.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\membrane.for

"$(INTDIR)\membrane.obj"	"$(INTDIR)\membrane.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Membrane_Bo.for

"$(INTDIR)\Membrane_Bo.obj"	"$(INTDIR)\Membrane_Bo.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\MPSmooth.for

"$(INTDIR)\MPSmooth.obj"	"$(INTDIR)\MPSmooth.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\NewCable.for

"$(INTDIR)\NewCable.obj"	"$(INTDIR)\NewCable.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\newsoli.for

"$(INTDIR)\newsoli.obj"	"$(INTDIR)\newsoli.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\NewSolver.for

"$(INTDIR)\NewSolver.obj"	"$(INTDIR)\NewSolver.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\NLDSUP.for

"$(INTDIR)\NLDSUP.obj"	"$(INTDIR)\NLDSUP.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\nldyna.for

"$(INTDIR)\nldyna.obj"	"$(INTDIR)\nldyna.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Nonconso.for

"$(INTDIR)\Nonconso.obj"	"$(INTDIR)\Nonconso.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Nonlinearlink.for

"$(INTDIR)\Nonlinearlink.obj"	"$(INTDIR)\Nonlinearlink.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\Password.for

"$(INTDIR)\Password.obj"	"$(INTDIR)\Password.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=".\Plasticity Hard.for"

"$(INTDIR)\Plasticity Hard.obj"	"$(INTDIR)\Plasticity Hard.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\PLVISCO.FOR

"$(INTDIR)\PLVISCO.OBJ"	"$(INTDIR)\PLVISCO.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\RCShell.for

"$(INTDIR)\RCShell.obj"	"$(INTDIR)\RCShell.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\REACTION.FOR

"$(INTDIR)\REACTION.OBJ"	"$(INTDIR)\REACTION.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=".\Section prop.for"

"$(INTDIR)\Section prop.obj"	"$(INTDIR)\Section prop.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\Seepage.for

"$(INTDIR)\Seepage.obj"	"$(INTDIR)\Seepage.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\SEISMA.for

"$(INTDIR)\SEISMA.obj"	"$(INTDIR)\SEISMA.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\SelfweightLoad.for

"$(INTDIR)\SelfweightLoad.obj"	"$(INTDIR)\SelfweightLoad.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\Settlement.for

"$(INTDIR)\Settlement.obj"	"$(INTDIR)\Settlement.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\sh3nod.for

"$(INTDIR)\sh3nod.obj"	"$(INTDIR)\sh3nod.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\shdrill.for

"$(INTDIR)\shdrill.obj"	"$(INTDIR)\shdrill.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=".\Shell4 Jacob.for"

"$(INTDIR)\Shell4 Jacob.obj"	"$(INTDIR)\Shell4 Jacob.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\Shell6.for

"$(INTDIR)\Shell6.obj"	"$(INTDIR)\Shell6.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\SOLID.FOR

"$(INTDIR)\SOLID.OBJ"	"$(INTDIR)\SOLID.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Solid_Bo.for

"$(INTDIR)\Solid_Bo.obj"	"$(INTDIR)\Solid_Bo.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Solidaverstr.for

"$(INTDIR)\Solidaverstr.obj"	"$(INTDIR)\Solidaverstr.sbr" : $(SOURCE)\
 "$(INTDIR)"


SOURCE=.\SOLIDECONSO.FOR

"$(INTDIR)\SOLIDECONSO.OBJ"	"$(INTDIR)\SOLIDECONSO.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\SOLIDHBD3.FOR

"$(INTDIR)\SOLIDHBD3.OBJ"	"$(INTDIR)\SOLIDHBD3.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\SPConcrete.for

"$(INTDIR)\SPConcrete.obj"	"$(INTDIR)\SPConcrete.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Spring1.for

"$(INTDIR)\Spring1.obj"	"$(INTDIR)\Spring1.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\TaperSec.for

"$(INTDIR)\TaperSec.obj"	"$(INTDIR)\TaperSec.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\TendonForce.for

"$(INTDIR)\TendonForce.obj"	"$(INTDIR)\TendonForce.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Thin3D.for

"$(INTDIR)\Thin3D.obj"	"$(INTDIR)\Thin3D.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\ThinSolid.for

"$(INTDIR)\ThinSolid.obj"	"$(INTDIR)\ThinSolid.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\truss.for

"$(INTDIR)\truss.obj"	"$(INTDIR)\truss.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\Warray.for

"$(INTDIR)\Warray.obj"	"$(INTDIR)\Warray.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\xfinas.for

"$(INTDIR)\xfinas.obj"	"$(INTDIR)\xfinas.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCCHANGE.FOR

"$(INTDIR)\XPCCHANGE.OBJ"	"$(INTDIR)\XPCCHANGE.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCCONCRETE.FOR

"$(INTDIR)\XPCCONCRETE.OBJ"	"$(INTDIR)\XPCCONCRETE.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCCSTATE.FOR

"$(INTDIR)\XPCCSTATE.OBJ"	"$(INTDIR)\XPCCSTATE.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCINPUT.FOR

"$(INTDIR)\XPCINPUT.OBJ"	"$(INTDIR)\XPCINPUT.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCLOAD.FOR

"$(INTDIR)\XPCLOAD.OBJ"	"$(INTDIR)\XPCLOAD.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCNASCAB.FOR

"$(INTDIR)\XPCNASCAB.OBJ"	"$(INTDIR)\XPCNASCAB.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCOUTPUT.FOR

"$(INTDIR)\XPCOUTPUT.OBJ"	"$(INTDIR)\XPCOUTPUT.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCSOLVE.FOR

"$(INTDIR)\XPCSOLVE.OBJ"	"$(INTDIR)\XPCSOLVE.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCSTATE.FOR

"$(INTDIR)\XPCSTATE.OBJ"	"$(INTDIR)\XPCSTATE.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCSTIFF.FOR

"$(INTDIR)\XPCSTIFF.OBJ"	"$(INTDIR)\XPCSTIFF.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCSTYILD.FOR

"$(INTDIR)\XPCSTYILD.OBJ"	"$(INTDIR)\XPCSTYILD.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCSTYINP.FOR

"$(INTDIR)\XPCSTYINP.OBJ"	"$(INTDIR)\XPCSTYINP.SBR" : $(SOURCE) "$(INTDIR)"


SOURCE=.\XPCUTILITY.FOR

"$(INTDIR)\XPCUTILITY.OBJ"	"$(INTDIR)\XPCUTILITY.SBR" : $(SOURCE) "$(INTDIR)"



!ENDIF 

