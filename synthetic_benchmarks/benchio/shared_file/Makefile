MF=     Makefile

# Compilers
FC=     mpiifort
CC=     mpiicc
#
#
FFLAGS= -g

# Comment out a line below to disable the named tests
FFLAGS+= -DWITH_SERIAL
FFLAGS+= -DWITH_MPIIO
# FFLAGS+= -DWITH_HDF5
# FFLAGS+= -DWITH_NETCDF

CFLAGS= $(FFLAGS)
LFLAGS=
#
EXE=	benchio.x
#
FSRC= \
	mpiio.F90      \
	netcdf.F90      \
	hdf5.F90      \
	benchclock.f90 \
	benchio.F90

CSRC=

# .SUFFIXES:
.SUFFIXES: .F90 .f90 .c .o

FOBJ=	$(FSRC:.F90=.o)
OBJ=    $(FOBJ:.f90=.o) $(CSRC:.c=.o)

FMOD=	$(FSRC:.F90=.mod)
MOD=    $(FMOD:.f90=.mod)

.F90.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ): $(MF)

clean:
	rm -f $(OBJ) $(EXE) $(MOD) core
