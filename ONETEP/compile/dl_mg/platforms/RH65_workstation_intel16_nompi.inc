#============================================================
# Makefile variable for Redhat 6.5 (and later) workstations
# in Southampton (Intel, no MPI)
#
# Load the Intel Parallel Studio XE 2016 module:
#
# module load intel/parallel_studio_xe_2016
#
# James C. Womack, 02/10/17
#
# Based on parallel_laptop.inc created by
# Lucian Anton, 20/04/12
#============================================================

FC := ifort

COMP  := intel
BUILD := opt

USE_OPENMP = yes

OMPFLAGS_intel_yes  := -openmp
OMPFLAGS_intel_no   :=
FFLAGS_intel_opt    := -O3

FFLAGS_intel_profile  := -DUSE_TIMER $(FFLAGS_intel_opt)

# enable lots of prints in V cycle
ifdef RELAX_PRINT
  FRELAX_PRINT := -DDLMG_RELAX_PRINT
endif

FFLAGS_intel_debug  :=  $(FRELAX_PRINT) -DDL_MG_TEST_ERR_CODE -O0 -g -traceback -check all,noarg_temp_created -ftrapuv

FFLAGS       := $(OMPFLAGS_$(COMP)_$(USE_OPENMP)) $(FFLAGS_$(COMP)_$(BUILD)) -module $(OBJDIR) -I$(LIBDIR) -I$(OBJDIR)
