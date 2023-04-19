$(error deprecated platform !!!)

FC  := mpif90
F77 := mpif77
CC  := mpicc

COMP := ibm
BUILD := opt
USE_OPENMP := yes

FLAGS_ibm_debug := -WF,-DMPI -g -C  
FLAGS_ibm_opt   := -WF,-DMPI  -O3  -qstrict
FLAGS_ibm_profile := -WF,-DUSE_TIMER $(FLAGS_ibm_opt)
OMPFLAGS_ibm_yes := -qsmp=omp
OMPFLAGS_ibm_no :=

FLAGS_gcc_opt := -DMPI -O3
FLAGS_gcc_debug := -DMPI -g -fcheck=all -finit-real=snan -Wall -Wtabs -Wno-unused
-fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core 
FLAGS_gcc_profile := -DUSE_TIMER $(FLAGS_gcc_opt)
OMPFLAGS_gcc_yes  := -fopenmp
OMPFLAGS_gcc_no   :=

FFLAGS := $(FLAGS_$(COMP)_$(BUILD))  $(OMPFLAGS_$(COMP)_$(USE_OPENMP))
#CFLAGS := -DNUM_UNDERSCORE=0

LD_FLAGS := 

#LD_FLAGS := -L/gpfs/packages/ibm/fftw/2.1.5/lib -ldfftw_mpi -ldrfftw_mpi -ldfftw -ldrf#ftw \
#            -L/gpfs/home/SCD/jpf02/lxa24-jpf02/MUMPS_4.10.0/lib -ldmumps -lmumps_commo#n  \
#            -L/gpfs/home/SCD/jpf02/lxa24-jpf02/ParMetis-3.2.0 -lparmetis -lmetis \
#            -L/gpfs/packages/ibm/scalapack/2.0.2/lib -lscalapack \
#            -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64 -lesslbg \
#            -L/gpfs/packages/ibm/lapack/3.4.2/lib  -llapack

# alternative  fftw 
# -L/gpfs/packages/ibm/fftw/2.1.5/lib -ldfftw_mpi -ldrfftw_mpi -ldfftw -ldrfftw \
#-L/gpfs/home/SCD/jpf02/lxa24-jpf02/fftw-2.1.5/mpi/.libs -ldfftw_mpi -ldrfftw_mpi \
#            -L/gpfs/home/SCD/jpf02/lxa24-jpf02/fftw-2.1.5/fftw/.libs -ldfftw \
#            -L/gpfs/home/SCD/jpf02/lxa24-jpf02/fftw-2.1.5/rfftw/.libs -ldrfftw  \
#            -L/gpfs/home/SCD/jpf02/lxa24-jpf02/MUMPS_4.10.0/lib/ -ldmumps -lmumps_common -lpord  \

ifeq ($(strip $(COMP)),gcc)
  MODULE_PRG := gccmpi
else ifeq ($(strip $(COMP)),ibm)
  MODULE_PRG := ibmmpi
endif

LOAD_MODULES := module unload ibmmpi && module unload gccmpi && module load $(MODULE_PRG) &&

