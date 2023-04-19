#============================================================
# Makefile variable for Thomas (MMM Hub) using ifort 17 and
# and Intel MPI
#
# Unload any MPI or compiler modules which may conflict with
# the Intel compiler 17 and Intel MPI modules
#
# Load the Intel compiler 17 and corresponding Intel MPI
# modules:
# 	module load compilers/intel/2017/update1
# 	module load mpi/intel/2017/update1/intel
#
# These are in the set of default modules loaded by thomas
# when /shared/ucl/apps/bin/defmods is sourced in .bashrc.
#
# James C. Womack, 17/10/2017
#
# Based on parallel_laptop.inc created by
# Lucian Anton, 20/04/12
#============================================================

FC := mpiifort

COMP  := intel
BUILD := opt
#BUILD := profile

USE_OPENMP = yes

MPIFLAGS          := -DMPI

ifdef USE_INCLUDE_MPI
  MPIFLAGS += -DUSE_INCLUDE_MPIF
endif


OMPFLAGS_intel_yes  := -qopenmp
OMPFLAGS_intel_no   :=
FFLAGS_intel_opt    := -O3

FFLAGS_intel_profile  := -DUSE_TIMER $(FFLAGS_intel_opt)

# enable lots of prints in V cycle
ifdef RELAX_PRINT
  FRELAX_PRINT := -DDLMG_RELAX_PRINT
endif

FFLAGS_intel_debug  :=  $(FRELAX_PRINT) -DDL_MG_TEST_ERR_CODE -O0 -g -traceback -check all,noarg_temp_created -ftrapuv

FFLAGS       := $(OMPFLAGS_$(COMP)_$(USE_OPENMP)) $(FFLAGS_$(COMP)_$(BUILD)) -module $(OBJDIR) -I$(LIBDIR) -I$(OBJDIR) $(MPIFLAGS)