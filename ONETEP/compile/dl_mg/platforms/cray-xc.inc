FC := ftn
CC := cc

BUILD=opt

USE_OPENMP := yes

# select compiler flags according the PE_ENV
ifeq ($(strip $(PE_ENV)), CRAY)
  COMP = cray
else ifeq ($(strip $(PE_ENV)), GNU)
  COMP = gcc
else ifeq ($(strip $(PE_ENV)), INTEL)
  COMP = intel
else
  $(error unknown PE_ENV value : $(PE_ENV)
endif

# Cray
FLAGS_cray_opt :=  -ram
FLAGS_cray_debug :=  -e D
FLAGS_cray_profile := -DUSE_TIMER $(FLAGS_cray_opt)
OMPFLAGS_cray_yes :=
OMPFLAGS_cray_no  := -O noomp
FLAGS_cray := -DMPI -J $(OBJDIR) -I$(OBJDIR) -I$(LIBDIR) $(FLAGS_cray_$(BUILD))

# INTEL
FLAGS_intel_opt := -O3 -diag-disable 8290,8291
FLAGS_intel_debug := -g -C -check all -traceback -diag-enable openmp  -diag-disable 8290,8291
FLAGS_intel_profile := -DUSE_TIMER $(FLAGS_intel_opt)
OMPFLAGS_intel_yes := -openmp
OMPFLAGS_intel_no :=
FLAGS_intel := -DMPI -module $(OBJDIR) -I$(OBJDIR) -I$(LIBDIR) $(FLAGS_intel_$(BUILD))

# GCC
FLAGS_gcc_opt := -Ofast
FLAGS_gcc_debug := -g -fcheck=all -finit-real=snan -Wall -Wtabs -Wno-unused
#-DDLMG_RELAX_PRINT
#-fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core
FLAGS_gcc_profile := -DUSE_TIMER $(FLAGS_gcc_opt)
OMPFLAGS_gcc_yes  := -fopenmp
OMPFLAGS_gcc_no   :=
FLAGS_gcc = -DMPI -J$(OBJDIR) -I$(OBJDIR) -I$(LIBDIR) $(FLAGS_gcc_$(BUILD))

# the flags are:
FFLAGS := $(FLAGS_$(COMP)) $(OMPFLAGS_$(COMP)_$(USE_OPENMP))

