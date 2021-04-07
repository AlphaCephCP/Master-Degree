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

F90 = mpifort 

# OMP = -fopenmp 

#------------------------------------------------------------------------------
# LIBRARIES
#------------------------------------------------------------------------------

LIBDIR      = ./lib
#LTEC   = -ltecio
LTEC = 

LFLAGS = $(LTEC) 

# COMPILATION OPTIMIZATION FLAGS for ifort
# F90FLAGS = -c -O2 -r8 -align -fpp  

# COMPILATION OPTIMIZATION FLAGS for gfortran
F90FLAGS = -c -O2 -cpp -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fcray-pointer


# PREPROCESSOR DEFINITIONS 
#FCFLAGS = -DGRMHD -DFLUX -DNOLIMITER -DSOURCE -DNONCONSERVATIVE -DNOIO -DPARALLEL -DRKDG               
FCFLAGS = -DCCZ4EINSTEIN -DNOFLUX -DNOLIMITER -DSOURCE -DNONCONSERVATIVE -DNOIO -DPARALLEL

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
          Bessel.o                  \
          f2kcli.o                  \
          typesDef.o                \
          DGBasis.o                 \
          recipies.o                \
          Metric.o                  \
          GRMHD.o                   \
          PDE.o                     \
          PDEPrim.o                 \
          EinsteinConstraints.o     \
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
          RKDG3D.o                  \
          ADERDG3D.o                

OBJS       = $(F90_OBJ)
FOR_FILES  = $(OBJS:.o=.f90)

Executable: $(EXE)

$(EXE): Makefile $(OBJS)
	@echo "   "
	@echo " Now Linking..."
	@echo " --------------"
	$(F90) $(OBJS) -L$(LIBDIR) $(LFLAGS) ${OMP} -o $(EXE)
	@echo "   "
	@echo " Compiled ! "
	@echo " ---------- "

#
clean:
	rm -f core* *.M *.o $(EXE) *.kmo *.mod 


