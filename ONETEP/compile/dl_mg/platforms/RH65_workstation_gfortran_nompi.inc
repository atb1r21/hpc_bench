#============================================================
# Makefile variable for Redhat 6.5 (and later) workstations
# in Southampton (gfortran, no MPI)
#
# Load the gcc 4.9.1
#
# module load gcc/4.9.1
#
# James C. Womack, 04/10/17
#
# Based on parallel_laptop.inc created by
# Lucian Anton, 20/04/12
#============================================================

FC := gfortran

COMP := gnu
BUILD := debug

USE_OPENMP = yes

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

FFLAGS       :=  $(OMPFLAGS_$(COMP)_$(USE_OPENMP)) $(FFLAGS_$(COMP)_$(BUILD)) -J$(OBJDIR) -I$(LIBDIR) -I$(OBJDIR)
