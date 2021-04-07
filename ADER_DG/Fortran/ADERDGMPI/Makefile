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

F90 = mpiifort

#------------------------------------------------------------------------------
# LIBRARIES
#------------------------------------------------------------------------------

LIBDIR      = ./lib
LTEC   = -ltecio

LFLAGS = $(LTEC) -cxxlib

# COMPILATION OPTIMIZATION FLAGS for ifort
F90FLAGS = -c -O2 -r8 -align -fpp  

# PREPROCESSOR DEFINITIONS 
FCFLAGS = -DELASTICITY -DNOLIMITER -DPARALLEL             

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
          f2kcli.o                  \
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
          Limiter.o                 \
          AnalyseDG.o               \
          MPITools.o                \
          PointSources.o            \
          ADERDG3D.o                

OBJS       = $(F90_OBJ)
FOR_FILES  = $(OBJS:.o=.f90)

Executable: $(EXE)

$(EXE): Makefile $(OBJS)
	@echo "   "
	@echo " Now Linking..."
	@echo " --------------"
	$(F90) $(OBJS) $(LFLAGS)  -o $(EXE)
	@echo "   "
	@echo " Compiled ! "
	@echo " ---------- "

#
clean:
	rm -f core* *.M *.o $(EXE) *.kmo *.mod 


