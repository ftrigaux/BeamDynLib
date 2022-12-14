#================================================================================#
# This makefile copied from B. Jonkman by J. Michalakes on 29-Jan-2013,          #
# adapted from Crunch (M. Buhl on 25-Jan-2013).                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
# It was designed to be used with:                                               #
#     Module1        (v1.00.04,      ???????????)                       #
#     NWTC Subroutine Library (CSC version    29-Jan-2013)                       #
#                                                                                #
# Older versions of ModuleName and NWTC Library may not work with this makefile. #
#================================================================================#

	# 32-bit or 64-bit?
#BITS = 32
BITS = 64
OS = Linux


	# Location of source files.  You will probably need to change these for your system.

ifeq ($(OS),Windows_NT)
	FAST_dir = C:/Users/bjonkman/Documents/DATA/DesignCodes/simulators/FAST/SVNdirectory/branches/BJonkman/

	REGISTRY = $(FAST_dir)/bin/Registry_Win32.exe
	LAPACK_LINK  = -llapack -lblas -LC:/LAPACK/win32
else
	FAST_dir=

	REGISTRY   =
	LAPACK_LINK  = -llapack -lblas
endif

BD_DIR = src
LIB_DIR = src/nwtc-library/src
NETLIB_DIR = src/nwtc-library/src/NetLib/lapack/
VERSION_DIR = src/version/src/

	# Name of compiler to use and flags to use.

FC     = gfortran
OPT    =
#FFLAGS = -O2 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input  -fcheck=bounds -C
#FFLAGS = -O3 -fbacktrace -ffree-line-length-none -x f95-cpp-input
#LDFLAGS = -O2 -m$(BITS)  -fbacktrace $(BLAS_LAPACK_LIBS)
FFLAGS  = $(OPT) -g -m$(BITS)  -fdefault-real-8 -fbacktrace -ffree-line-length-none -x f95-cpp-input -Wsurprising -DDOUBLE_PRECISION -fpic #-fcheck=all
LDFLAGS = $(OPT) -g -m$(BITS) -fbacktrace


	# Precision.

	#==========================================================#
	# You should not need to change anything beyond this point #
	#==========================================================#

	# System-specific settings.

ifeq ($(OS),Windows_NT)
		# Windows
	DEL_CMD   = del
	EXE_EXT   = _gwin$(BITS).exe
	INTER_DIR = Obj_win$(BITS)
	MD_CMD    = @mkdir
	OBJ_EXT   = .obj
	PATH_SEP  = \\
	SYS_FILE  = SysGnuWin
else
		# Linux
		# Linux
	DEL_CMD   = rm -f
	EXE_EXT   = 
	INTER_DIR = obj
	MD_CMD    = @mkdir -p
	OBJ_EXT   = .o
	PATH_SEP  = /
	SYS_FILE  = SysGnuLinux
endif

	# Destination and RootName for executable

OUTPUT_NAME = BeamDyn
DEST_DIR    = .

	# Library files.

LIB_SOURCES =           \
	SingPrec.f90         \
	NWTC_Base.f90        \
	$(SYS_FILE).f90      \
	 NWTC_Library_Types.f90   \
	NWTC_IO.f90          \
	NWTC_Num.f90         \
	ModMesh_Types.f90    \
	ModMesh.f90          \
	ModMesh_Mapping.f90  \
	NWTC_Library.f90     \

NETLIB_SOURCES=             \
		  NWTC_LAPACK.f90

VERSION_SOURCES =	\
	VersionInfo.f90		

BD_SOURCES   =            \
	BeamDyn.f90           \
	BeamDyn_BldNdOuts_IO.f90 \
	BeamDyn_IO.f90        \
	BeamDyn_Subs.f90      \
	BeamDyn_Types.f90     \
	BeamDynLib_Types.f90  \
	BeamDynLib.f90        \
	BeamDynLib_CBind.f90 

vpath %.f90 $(LIB_DIR) $(NETLIB_DIR) $(VERSION_DIR) $(BD_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)
vpath %.txt $(BD_DIR)

ALL_SOURCES = $(LIB_SOURCES) $(NETLIB_SOURCES) $(VERSION_SOURCES) $(BD_SOURCES)
ALL_OBJS    = $(ALL_SOURCES:.f90=.obj)
ALL_OBJS   := $(ALL_OBJS:.F90=.obj)       #note the upper case here (from IceFloe)
ALL_OBJS   := $(ALL_OBJS:.f=.obj)



	# Rule to do everything.
all:     default
default: $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)

	# General rule for making the files.

# -B is needed for MinGW version of Gfortran
%.obj: %.f90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -g -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)
#-B $(INTER_DIR)

	#  Dependency rules.

#NWTC Library dependency rules:
NWTC_Base.obj:              SingPrec.obj
$(SYS_FILE).obj:            NWTC_Base.obj
NWTC_Library_Types.obj:     $(SYS_FILE).obj
NWTC_IO.obj:                NWTC_Library_Types.obj VersionInfo.obj
NWTC_Num.obj:               NWTC_IO.obj
ModMesh_Types.obj:          NWTC_Num.obj
ModMesh.obj:                ModMesh_Types.obj
ModMesh_Mapping.obj:        ModMesh.obj NWTC_LAPACK.obj
NWTC_Library.obj:           ModMesh.obj  ModMesh_Mapping.obj

NWTC_LAPACK.obj:            NWTC_Base.obj

BeamDyn_Types.obj:       NWTC_Library.obj  $(BD_DIR)/BeamDyn_Types.f90
BeamDyn_Subs.obj:        BeamDyn_Types.obj
BeamDyn_IO.obj:          BeamDyn_Types.obj  BeamDyn_Subs.obj BeamDyn_BldNdOuts_IO.obj
BeamDyn.obj:             BeamDyn_IO.obj  BeamDyn_Subs.obj

Driver_Beam_Subs.obj:   BeamDyn.obj NWTC_Library.obj
Driver_Beam.obj:        Driver_Beam_Subs.obj

	# Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


	# Run the registry if the input file changes.

#$(BD_DIR)/BeamDyn_Types.f90: $(BD_DIR)/Registry_BeamDyn.txt
#	$(REGISTRY) $< -I $(LIB_DIR) -O $(BD_DIR)


	# For compiling the driver/glue code.

OBJC  = CBeamDyn.o
OBJDIR=obj
DIRC=src
OBJC := $(OBJC:%.o=$(OBJDIR)/%.o)
CC=gcc
OPTIFLAG=-g
CFLAGS=$(OPTIFLAG) -fPIC
$(OBJC): $(OBJDIR)/%.o: $(DIRC)/%.c
	$(CC) $(CFLAGS) -c  $< -o $@ 
INSTALL_DIR=/usr/local

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(ALL_OBJS) $(OBJC) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src))) obj/CBeamDyn.o $(LAPACK_LINK)

	# Cleanup afterwards.

lib: $(ALL_OBJS) $(OBJC) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -shared -o $(DEST_DIR)/lib$(OUTPUT_NAME).so \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src))) obj/CBeamDyn.o $(LAPACK_LINK)

fprog: $(ALL_OBJS) BeamDynLib_Program.obj $(OBJC) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src))) $(INTER_DIR)/BeamDynLib_Program.obj $(LAPACK_LINK)

install:
	ln -sf `pwd`/$(DEST_DIR)/lib$(OUTPUT_NAME).so $(INSTALL_DIR)/lib/
	ln -sf `pwd`/$(DIRC)/CBeamDyn.h $(INSTALL_DIR)/include/

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj $(INTER_DIR)$(PATH_SEP)*.o $(OUTPUT_NAME)$(EXE_EXT) \
	*.dat *.out $(DEST_DIR)/lib$(OUTPUT_NAME).so
