# Wrapper for CASTEP config files
#
# James C. Womack, 10/2016
# Updated by James C. Womack, 02/2017

include $(CASTEP_ROOT)/obj/$(CASTEP_ARCH).mk

ifeq ($(COMMS_ARCH),mpi)
MPIFLAGS = -DMPI
else ifeq ($(COMMS_ARCH),serial)
MPIFLAGS =
else
$(error COMMS_ARCH should be one of "mpi","serial" (not $(COMMS_ARCH)))
endif

# To compile DL_MG object files and module we need to add the $(OBJDIR)
# and $(LIBDIR) (should be passed to the call to make from ONETEP's
# Makefile) to include path in $(FFLAGS).
FFLAGS +=  -I$(OBJDIR) -I$(LIBDIR)

# We also need to inform the compiler that module files should be placed
# in $(OBJDIR). This has a compiler-specific flag.
ifeq ($(COMPILER),INTEL-ifort-on-LINUX)
# Use Intel's syntax (different from GNU)
  MODFLAGS += -module $(OBJDIR)
# Use 'include "mpif.h"' rather than 'use mpi' because CASTEP's default
# platform files sometimes use the mpif90 wrapper and use the I_MPI_F90
# environment variable to set the compiled to ifort. In Intel Parallel Studio
# 2017, the mpif90 wrapper always includes mpi.mod for gfortran, even when
# I_MPI_F90=ifort. Using 'include "mpif.h"' avoids this issue.
ifeq ($(COMMS_ARCH),mpi)
  MPIFLAGS += -DUSE_INCLUDE_MPIF
endif
else
# Default to the the GNU syntax (may not work for all compilers)
  MODFLAGS += -J$(OBJDIR)
endif

FFLAGS += $(MPIFLAGS) $(MODFLAGS)

FC := $(F90)
