#$(error deprecated platform!!!)

#============================================================
# Makefile variable for laptop (gfortran + MPICH2)
#
# Lucian Anton
# created: 20/04/12
#============================================================

FC := mpif90
F77 := mpif77

COMP := gnu

USE_OPENMP = yes

MPIFLAGS          := -DMPI

ifdef USE_INCLUDE_MPI
  MPIFLAGS += -DUSE_INCLUDE_MPIF
endif

ifdef TEST_ERR_CODES
  TEST_ERR := -DDL_MG_TEST_ERR_CODES
endif


OMPFLAGS_gnu_yes  := -fopenmp
OMPFLAGS_gnu_no   :=
FFLAGS_gnu_opt    := -g -O3

FFLAGS_gnu_profile  := -DUSE_TIMER $(FFLAGS_gnu_opt)

# enable lots of prints in V cycle
ifdef RELAX_PRINT
  FRELAX_PRINT := -DDLMG_RELAX_PRINT
endif

MODS_DESTINATION_gnu = -J$(OBJDIR)

# -DDUMP_DATA -DUSE_DUMP
FFLAGS_gnu_debug := -g -fcheck=all -Wall -finit-real=snan -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -DUSE_BACKTRACE
#FFLAGS_gnu_debug := -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fsignaling-nans
#FFLAGS_gnu_debug  :=  $(FRELAX_PRINT) -O0 -g -fbounds-check -Wall -Wuninitialized -Wno-unused -finit-real=nan -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -fbacktrace

# intel flags
OMPFLAGS_intel_yes  := -qopenmp
OMPFLAGS_intel_no   :=
FFLAGS_intel_opt    := -g -O3

FFLAGS_intel_profile  := -DUSE_TIMER $(FFLAGS_intel_opt)

FFLAGS_intel_debug := -g -check noarg_temp_created -traceback -DUSE_BACKTARCE

MODS_DESTINATION_intel := -module $(OBJDIR)


FFLAGS  :=  $(MPIFLAGS) $(FFLAGS_$(COMP)_$(BUILD)) $(OMPFLAGS_$(COMP)_$(USE_OPENMP)) -DHAVE_CONTIGUOUS $(TEST_ERR)

INCPATHS := $(MODS_DESTINATION_$(COMP)) -I$(OBJDIR)
