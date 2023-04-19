#============================================================
# Makefile variable for Thomas (MMM Hub) using gfortran 4.9.2
# and OpenMPI
#
# Unload any MPI or compiler modules which may conflict with
# the GCC 4.9.2 and OpenMPI modules.
#
# Load the GCC 4.9.2 and corresponding OpenMPI modules:
# 	module load compilers/gnu/4.9.2
# 	module load mpi/openmpi/1.10.1/gnu-4.9.2
#
# James C. Womack, 17/10/2017
#
# Based on parallel_laptop.inc created by
# Lucian Anton, 20/04/12
#============================================================

FC := mpif90
F77 := mpif77

COMP := gnu
BUILD := opt
#BUILD := profile

USE_OPENMP = yes

MPIFLAGS          := -DMPI

ifdef USE_INCLUDE_MPI
  MPIFLAGS += -DUSE_INCLUDE_MPIF
endif


OMPFLAGS_gnu_yes  := -fopenmp
OMPFLAGS_gnu_no   :=
FFLAGS_gnu_opt    := -O3

FFLAGS_gnu_profile  := -DUSE_TIMER $(FFLAGS_gnu_opt)

# enable lots of prints in V cycle
ifdef RELAX_PRINT
  FRELAX_PRINT := -DDLMG_RELAX_PRINT
endif

FFLAGS_gnu_debug  :=  $(FRELAX_PRINT) -DDL_MG_TEST_ERR_CODE -O0 -g -fbounds-check -Wall -Wuninitialized -Wno-unused -finit-real=nan -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -fbacktrace

# -DDUMP_DATA -DUSE_DUMP
#FFLAGS_gnu_debug := -g -fcheck=all -Wall -finit-real=snan -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core
#FFLAGS_gnu_debug := -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fsignaling-nans

FFLAGS       :=  -J$(OBJDIR) -I$(LIBDIR) -I$(OBJDIR) $(MPIFLAGS) $(FFLAGS_$(COMP)_$(BUILD)) $(OMPFLAGS_$(COMP)_$(USE_OPENMP))
