# Makefile for SPECFEM3D_GEOTECH
# REVISION
#   HNG, Jul 12,2011; HNG, Apr 09,2010
# binary directory
BIN = ../bin

# serial fortran compiler
#FC = gfortran
#FFLAGS = -g -std=f2003 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wline-truncation -fno-trapping-math -fno-backslash

# parallel fortran compiler
#FC = mpif90 # OR path to appropriate MPI forrtan compiler
FC = mpiifort # OR path to appropriate MPI forrtan compiler
#FFLAGS = -g -std=f2003 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wline-truncation -fno-trapping-math -fno-backslash -Wunused
FFLAGS = -std08 -O3 -g -axAVX -xHost -mcmodel medium -shared-intel -check all -check bounds -fp-stack-check -traceback -ftz -funroll-loops -unroll5 -align sequence -qopt-report0 -qopt-report-phase=vec -implicitnone -warn truncated_source -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds
FC_COMPILE = $(FC) $(FFLAGS)

# C compiler
#CC = gcc
CC = icc
CFLAGS = -g -O3 -Werror -Wall -pedantic-errors
CC_COMPILE = $(CC) $(CFLAGS)

LIB =
#LIBMPI = -L/~/openmpi-1.4.3/lib/lmpi # OR path to appropriate MPI library
#LIBSCOTCH = -L/home/hgharti/lsoft/scotch_6.0.4-gcc-5.3.1/lib -lscotch -lscotcherr
#INCSCOTCH = -I$(HOME)/lsoft/scotch_6.0.4-gcc-5.3.1/include
LIBSCOTCH = -L/home/hgharti/lsoft/scotch_6.0.4-intel17/lib -lscotch -lscotcherr
INCSCOTCH = -I$(HOME)/lsoft/scotch_6.0.4-intel17/include

# object tree
OBJ = $(BIN)/global.o           \
      $(BIN)/string_library.o   \
      $(BIN)/math_library.o     \
      $(BIN)/shape_library.o    \
      $(BIN)/element.o          \
      $(BIN)/elastic.o          \
      $(BIN)/read_input.o       \
      $(BIN)/gll_library.o      \
      $(BIN)/mesh_spec.o        \
      $(BIN)/integration.o      \
      $(BIN)/apply_bc.o         \
      $(BIN)/preprocess.o       \
      $(BIN)/apply_traction.o   \
      $(BIN)/compute_pressure.o \
      $(BIN)/plastic_library.o  \
      $(BIN)/visual.o           \
      $(BIN)/postprocess.o

# object tree for EXCAVATION
OBJ_EXCAV = $(BIN)/excavation_library.o

# object tree for SEM
OBJ_SEM = $(BIN)/serial_library.o      \
          $(BIN)/math_library_serial.o \
          $(BIN)/solver.o

# object tree for SEMMPI parallel
OBJ_SEMMPI = $(BIN)/set_precision_mpi.o \
             $(BIN)/mpi_library.o       \
             $(BIN)/math_library_mpi.o  \
             $(BIN)/ghost_library_mpi.o \
             $(BIN)/solver_mpi.o

# object tree for SEM
OBJ_SEM = $(BIN)/serial_library.o      \
          $(BIN)/math_library_serial.o \
          $(BIN)/solver.o

# object tree for PARTMESH
OBJ_PARTMESH = $(BIN)/global.o           \
               $(BIN)/string_library.o   \
               $(BIN)/serial_library.o   \
               $(BIN)/math_library.o     \
               $(BIN)/partmesh_library.o \
               $(BIN)/read_input.o       \
	             $(BIN)/partmesh_scotch.o  \
               $(BIN)/partmesh.o

OBJ_semgeotech = $(OBJ) $(OBJ_EXCAV) $(OBJ_SEM) $(BIN)/slope.o $(BIN)/excavation.o $(BIN)/semgeotech.o
OBJ_psemgeotech = $(OBJ) $(OBJ_SEMMPI) $(BIN)/pslope.o $(OBJ_EXCAV) $(BIN)/pexcavation.o $(BIN)/psemgeotech.o

# main executable for spectral-element method
semgeotech = $(BIN)/semgeotech
psemgeotech = $(BIN)/psemgeotech
PARTMESH = $(BIN)/partmesh

default: semgeotech psemgeotech

clean:
	rm -f *.mod $(BIN)/*.o core

all: partmesh semgeotech psemgeotech

semgeotech: $(OBJ_semgeotech)
	$(FC_COMPILE) -o $(semgeotech) $(OBJ_semgeotech) $(LIB)

psemgeotech: $(OBJ_psemgeotech)
	$(FC_COMPILE) -o $(psemgeotech) $(OBJ_psemgeotech) $(LIBMPI)

partmesh: $(OBJ_PARTMESH)
	$(FC_COMPILE) -o ${PARTMESH} $(OBJ_PARTMESH) $(LIBSCOTCH)

# compile instruction for each source
$(BIN)/set_precision_mpi.o: set_precision_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/set_precision_mpi.o set_precision_mpi.f90

$(BIN)/global.o: global.f90
	$(FC_COMPILE) -c -o $(BIN)/global.o global.f90

$(BIN)/string_library.o: string_library.f90
	$(FC_COMPILE) -c -o $(BIN)/string_library.o string_library.f90

$(BIN)/math_library.o: math_library.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library.o math_library.f90

$(BIN)/mpi_library.o: mpi_library.f90
	$(FC_COMPILE) -c -o $(BIN)/mpi_library.o mpi_library.f90

$(BIN)/serial_library.o: serial_library.f90
	$(FC_COMPILE) -c -o $(BIN)/serial_library.o serial_library.f90

$(BIN)/math_library_serial.o: math_library_serial.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library_serial.o math_library_serial.f90

$(BIN)/math_library_mpi.o: math_library_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library_mpi.o math_library_mpi.f90

$(BIN)/read_input.o: read_input.f90
	$(FC_COMPILE) -c -o $(BIN)/read_input.o read_input.f90

$(BIN)/gll_library.o: gll_library.f90
	$(FC_COMPILE) -c -o $(BIN)/gll_library.o gll_library.f90

$(BIN)/mesh_spec.o: mesh_spec.f90
	$(FC_COMPILE) -c -o $(BIN)/mesh_spec.o mesh_spec.f90

$(BIN)/shape_library.o: shape_library.f90
	$(FC_COMPILE) -c -o $(BIN)/shape_library.o shape_library.f90

$(BIN)/element.o: element.f90
	$(FC_COMPILE) -c -o $(BIN)/element.o element.f90

$(BIN)/integration.o: integration.f90
	$(FC_COMPILE) -c -o $(BIN)/integration.o integration.f90

$(BIN)/elastic.o: elastic.f90
	$(FC_COMPILE) -c -o $(BIN)/elastic.o elastic.f90

$(BIN)/gauss_library.o: gauss_library.f90
	$(FC_COMPILE) -c -o $(BIN)/gauss_library.o gauss_library.f90

$(BIN)/apply_bc.o: apply_bc.f90
	$(FC_COMPILE) -c -o $(BIN)/apply_bc.o apply_bc.f90

$(BIN)/preprocess.o: preprocess.f90
	$(FC_COMPILE) -c -o $(BIN)/preprocess.o preprocess.f90

$(BIN)/apply_traction.o: apply_traction.f90
	$(FC_COMPILE) -c -o $(BIN)/apply_traction.o apply_traction.f90

$(BIN)/compute_pressure.o: compute_pressure.f90
	$(FC_COMPILE) -c -o $(BIN)/compute_pressure.o compute_pressure.f90

$(BIN)/ghost_library_mpi.o: ghost_library_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/ghost_library_mpi.o ghost_library_mpi.f90

$(BIN)/plastic_library.o: plastic_library.f90
	$(FC_COMPILE) -c -o $(BIN)/plastic_library.o plastic_library.f90

$(BIN)/visual.o: visual.f90
	$(FC_COMPILE) -c -o $(BIN)/visual.o visual.f90

$(BIN)/postprocess.o: postprocess.f90
	$(FC_COMPILE) -c -o $(BIN)/postprocess.o postprocess.f90

$(BIN)/solver.o: solver.f90
	$(FC_COMPILE) -c -o $(BIN)/solver.o solver.f90

$(BIN)/solver_mpi.o: solver_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/solver_mpi.o solver_mpi.f90

# following objective files are compiled from the same file but with
# different compiler preprocessing directives options.
# -DUSE_MPI provides the flag to use the MPI library for the compilation.
# the flag -DUSE_MPI is equivalent to #define USE_MPI 1.
$(BIN)/slope.o: slope.F90
	$(FC_COMPILE) -c -o $(BIN)/slope.o slope.F90

$(BIN)/pslope.o: slope.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/pslope.o slope.F90

$(BIN)/excavation_library.o: excavation_library.f90
	$(FC_COMPILE) -c -o $(BIN)/excavation_library.o excavation_library.f90

# following objective files are compiled from the same file but with
# different compiler preprocessing directives options.
# -DUSE_MPI provides the flag to use the MPI library for the compilation.
# the flag -DUSE_MPI is equivalent to #define USE_MPI 1.
$(BIN)/excavation.o: excavation.F90
	$(FC_COMPILE) -c -o $(BIN)/excavation.o excavation.F90

$(BIN)/pexcavation.o: excavation.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/pexcavation.o excavation.F90

# following objective files are compiled from the same file but with
# different compiler preprocessing directives options.
$(BIN)/semgeotech.o: semgeotech.F90
	$(FC_COMPILE) -c -o $(BIN)/semgeotech.o semgeotech.F90

$(BIN)/psemgeotech.o: semgeotech.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/psemgeotech.o semgeotech.F90

# mesh partition
$(BIN)/partmesh_library.o: partmesh_library.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh_library.o partmesh_library.f90

$(BIN)/partmesh_scotch.o: partmesh_scotch.f90 partmesh_library.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh_scotch.o partmesh_scotch.f90 $(INCSCOTCH)

$(BIN)/partmesh.o: partmesh.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh.o partmesh.f90


