# *********************************************************************
#
#  Makefile for ADERDG3D 
#  
# *********************************************************************
#
.SUFFIXES : .o .f90 .f .c .F90
.KEEP_STATE:
#
EXE  = ADERDG  
# *********************************************************************

F90 = gfortran
LD = g++

#------------------------------------------------------------------------------
# LIBRARIES
#------------------------------------------------------------------------------

#-LIBDIR      = ./lib
#-LTEC   = -ltecio

# download the tecio at http://www.tecplot.com/my/tecio-library/
# it is not open source but freeware. In the ExaHyPE repo we have
# the binary for Linux 64.26 included, it's 1/2MB of stuff.
LTEC   = ./libtecio.a

LFLAGS = $(LTEC) -lgfortran -pthread
# when linking with gfortran: -lstdc++
# with intel:  -cxxlib

# COMPILATION OPTIMIZATION FLAGS for ifort
F90FLAGS = -c -cpp  -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fcray-pointer
# with intel: -O2 -r8 -align -fpp  

# PREPROCESSOR DEFINITIONS 
FCFLAGS = -DEULER
# add -DLIMITER to activate limiter stuff. we don't compile it in here
# to be able to compute without limiter
# same with point sources

#
.SUFFIXES: .f90 .o .f77 .F90

.f90.o:
	@echo "   "
	@echo " Compiling Source File: " $<
	@echo " ---------------------- ===================="
	$(F90) $(F90FLAGS) $(FCFLAGS) $(INCLUDES) $<

.F90.o:
	@echo "   "
	@echo " Compiling Source File: " $<
	@echo " ---------------------- ===================="
	$(F90) $(F90FLAGS) $(FCFLAGS) $(INCLUDES) $<

.f.o:
	@echo "   "
	@echo " Compiling Source File: " $<
	@echo " ---------------------- ===================="
	$(F90) $(F90FLAGS) $(FCFLAGS) $(INCLUDES) $<

# Fortran files with no optimization
NO_OPT_OBJ = 

# Compiled object files

F90_OBJ = gauleg.o                  \
          GAUSS.o                   \
          typesDef.o                \
          DGBasis.o                 \
          PDE.o                     \
          CalcTimeStep.o            \
          IO.o                      \
          Init.o                    \
          SpaceTimePredictor.o      \
          VolumeIntegral.o          \
          BoundaryConditions.o      \
          RiemannSolver.o           \
          SurfaceIntegral.o         \
          ElementUpdate.o           \
          PointSources.o            \
          AnalyseDG.o               \
          ADERDG3D.o                

OBJS       = $(F90_OBJ)
FOR_FILES  = $(OBJS:.o=.f90)

Executable: $(EXE)

$(EXE): Makefile $(OBJS)
	@echo "   "
	@echo " Now Linking..."
	@echo " --------------"
	$(LD) $(OBJS) $(LFLAGS)  -o $(EXE)
	@echo "   "
	@echo " Compiled ! "
	@echo " ---------- "

#
clean:
	rm -f core* *.M *.o $(EXE) *.kmo *.mod 


