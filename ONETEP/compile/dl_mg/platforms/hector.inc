$(error deprecated platform !!!)

FC := ftn
CC := cc

COMP=gcc
BUILD=opt

USE_OPENMP := yes

# Compiler flags

# Cray
FLAGS_cray_opt := -DMPI -O ipa1 
FLAGS_cray_debug := -DMPI -e D 
FLAGS_cray_profile := -DUSE_TIMER $(FLAGS_cray_opt)
OMPFLAGS_cray_yes :=
OMPFLAGS_cray_no  := -O noomp

# PGI 
FLAGS_pgi_opt := -DMPI -fast 
FLAGS_pgi_debug := -DMPI -g -Mbounds -traceback 
FLAGS_pgi_profile := -DUSE_TIMER $(FLAGS_pgi_opt)

# GCC
FLAGS_gcc_opt := -DMPI -Ofast
FLAGS_gcc_debug := -DMPI -g -fcheck=all -finit-real=snan -Wall -Wtabs -Wno-unused
#-DDLMG_RELAX_PRINT
-fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core 
FLAGS_gcc_profile := -DUSE_TIMER $(FLAGS_gcc_opt)
OMPFLAGS_gcc_yes  := -fopenmp
OMPFLAGS_gcc_no   :=

# NAG
FLAGS_nag_opt   := -O3 
FLAGS_nag_debug := -g90 -gline -C=all 
FLAGS_nag_profile := -DUSE_TIMER $(FLAGS_nag_opt)

# the flags are:
FFLAGS = $(FLAGS_$(COMP)_$(BUILD)) $(OMPFLAGS_$(COMP)_$(USE_OPENMP))
CFLAGS = $(CFLAGS_$(COMP)_$(BUILD))

# change environment acorrding to COMP value
ifeq ($(strip $(COMP)),gcc)
  MODULE_PRG := PrgEnv-gnu
else ifeq ($(strip $(COMP)),pgi)
  MODULE_PRG := PrgEnv-pgi
else ifeq ($(strip $(COMP)),cray)
  MODULE_PRG := PrgEnv-cray
else
  $(error unknown value of COMP = $(COMP))
endif

LOAD_MODULES := . /etc/bash.bashrc.local && module unload PrgEnv-gnu && module unload PrgEnv-pgi &&  module unload PrgEnv-cray && module load $(MODULE_PRG) && 

# add perftools module for profiling build
ifeq ($(BUILD), craypat)
  FFLAGS = $(FLAGS_$(COMP)_opt)
  CFLAGS = $(CFLAGS_$(COMP)_opt)
  LOAD_MODULES +=  module load perftools &&
endif
