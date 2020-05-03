# Microsoft Developer Studio Project File - Name="Diffusion3D" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Diffusion3D - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Diffusion3D.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Diffusion3D.mak" CFG="Diffusion3D - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Diffusion3D - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Diffusion3D - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "Diffusion3D - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /check:bounds /compile_only /define:"MPI" /include:"Release/" /nologo /real_size:64 /warn:nofileopt
# SUBTRACT F90 /threads
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /GX /YX /FD /c
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /stack:0x3b9aca00 /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Diffusion3D - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /real_size:64 /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /stack:0x3b9aca00 /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Diffusion3D - Win32 Release"
# Name "Diffusion3D - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\DoInR.f90
DEP_F90_DOINR=\
	".\Release\isDiagPredomFunction.mod"\
	{$(INCLUDE)}"linear_operators.mod"\
	{$(INCLUDE)}"mpif.h"\
	
# End Source File
# Begin Source File

SOURCE=.\DoInX.f90
DEP_F90_DOINX=\
	".\Release\isDiagPredomFunction.mod"\
	{$(INCLUDE)}"linear_operators.mod"\
	{$(INCLUDE)}"mpif.h"\
	
# End Source File
# Begin Source File

SOURCE=.\DoInY.f90
DEP_F90_DOINY=\
	".\Release\isDiagPredomFunction.mod"\
	{$(INCLUDE)}"linear_operators.mod"\
	{$(INCLUDE)}"mpif.h"\
	
# End Source File
# Begin Source File

SOURCE=.\fRHS.f90
# End Source File
# Begin Source File

SOURCE=.\GetMass.f90
# End Source File
# Begin Source File

SOURCE=.\isDiagPredom.f90
# End Source File
# Begin Source File

SOURCE=.\kron.f90
# End Source File
# Begin Source File

SOURCE=.\SaveToDisk.f90
# End Source File
# Begin Source File

SOURCE=.\SwapXY.f90
# End Source File
# Begin Source File

SOURCE=.\SwapYX.f90
# End Source File
# Begin Source File

SOURCE=.\Test.f90
DEP_F90_TEST_=\
	{$(INCLUDE)}"linear_operators.mod"\
	{$(INCLUDE)}"mpif.h"\
	
NODEP_F90_TEST_=\
	".\Release\fRHSFunction.mod"\
	".\Release\GetMassSubroutine.mod"\
	".\Release\kronFunction.mod"\
	
# End Source File
# End Group
# Begin Group "Lib"

# PROP Default_Filter ""
# Begin Source File

SOURCE="C:\Program Files\DevStudio\DF\IMSL\LIB\SMATHD.LIB"
# End Source File
# Begin Source File

SOURCE="C:\Program Files\DevStudio\DF\IMSL\LIB\SF90MP.LIB"
# End Source File
# Begin Source File

SOURCE="C:\Program Files\Microsoft HPC Pack 2008 R2\Lib\i386\msmpi.lib"
# End Source File
# Begin Source File

SOURCE="C:\Program Files\Microsoft HPC Pack 2008 R2\Lib\i386\msmpifms.lib"
# End Source File
# Begin Source File

SOURCE="C:\Program Files\Microsoft HPC Pack 2008 R2\Lib\i386\msmpifec.lib"
# End Source File
# Begin Source File

SOURCE="C:\Program Files\Microsoft HPC Pack 2008 R2\Lib\i386\msmpifes.lib"
# End Source File
# Begin Source File

SOURCE="C:\Program Files\Microsoft HPC Pack 2008 R2\Lib\i386\msmpifmc.lib"
# End Source File
# End Group
# End Target
# End Project
