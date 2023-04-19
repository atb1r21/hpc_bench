#============================================================
# Makefile variable for Redhat 6.5 (and later) workstations
# in Southampton (gfortran + Intel MPI)
#
# Load the gcc 4.9.1 and Intel Parallel Studio XE 2016
# modules:
#
# module load gcc/4.9.1
# module load intel/parallel_studio_xe_2016
#
# James C. Womack, 23/09/16
#
# Based on parallel_laptop.inc created by
# Lucian Anton, 20/04/12
#============================================================

FC := mpif90
F77 := mpif77

COMP := gnu
BUILD := opt

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


