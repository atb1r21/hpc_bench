# Environment for compilation on CSC Desktops, University of Warwick
# Tested with module load intel impi imkl on 05/09/17

F90 = mpiifort -qopenmp -fPIC -assume buffered_io
FFLAGS = -DMKL_FFTW3 -DMPI -DSCALAPACK -I${MKLROOT}/include/intel64/lp64/
LIBS = ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
OPTFLAGS = -O2
DEBUGFLAGS = -g -C -check all -warn all -traceback -DDEBUG_ARRAYS
COMPILER = INTEL-ifort-on-LINUX
