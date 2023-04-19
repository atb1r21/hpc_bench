$(error deprecated platform !!!)

#============================================================
# Makefile variable for  intel mpi
#
# Lucian Anton
# created: 5/06/12
#============================================================

FC := source /opt/intel/composerxe/bin/compilervars.sh intel64 && \
      source /opt/intel/impi/4.1.0/intel64/bin/mpivars.sh && mpiifort
#F77 := $(FC)

#FC = module add intel/comp intel/mpi && mpiifort
FFLAGS := -DMPI

COMP := intel

USE_OPENMP = yes

MPIFLAGS         :=
OMPFLAGS_intel_yes := -openmp
OMPFLAGS_intel_no  :=
FFLAGS_intel_opt   := -O3
FFLAGS_intel_debug  := -g -C -fpe0 -check noarg_temp_created -traceback
#-DDLMG_RELAX_PRINT
FFLAGS  += $(FFLAGS_$(COMP)_$(BUILD)) $(OMPFLAGS_$(COMP)_$(USE_OPENMP))
