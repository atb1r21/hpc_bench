#============================================================
# Makefile variables for gfortran (serial)
#
# Lucian Anton
# created: 7/02/14
#============================================================

FC := gfortran

COMP := gnu

USE_OPENMP = yes

OMPFLAGS_gnu_yes  := -fopenmp
OMPFLAGS_gnu_no   :=
FFLAGS_gnu_opt    := -O3

FFLAGS_gnu_profile  := -DUSE_TIMER $(FFLAGS_gnu_opt)

# enable lots of prints in V cycle
ifdef RELAX_PRINT
  FRELAX_PRINT := -DDLMG_RELAX_PRINT
endif

# differit gfortran flags for differit compiler version
# is there a environment version variable?
FFLAGS_gnu_debug  := $(FRELAX_PRINT) -O0 -g -fbounds-check -Wall -Wuninitialized -Wno-unused -finit-real=nan -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -fbacktrace


#FFLAGS_gnu_debug :=  -g -fcheck=all -Wall -Wno-unused -finit-real=snan -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -fbacktrace
#FFLAGS_gnu_debug := -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fsignaling-nans


FFLAGS       :=  -J$(OBJDIR) -I$(OBJDIR) -I$(LIBDIR) $(FFLAGS_$(COMP)_$(BUILD)) $(OMPFLAGS_$(COMP)_$(USE_OPENMP))
