# write info
#
# Written by Jolyon Aarons, Southampton, 02/2014
# Modified by Nicholas Hine to skip compilation when MPI and SCALAPACK are
# unavailable, 03/2014,


include $(ONETEP_ROOT)/config/conf.$(ONETEP_ARCH)

ifeq (,$(findstring -DESSL,$(FFLAGS)))
SOURCES = elpa1.F90 elpa2.F90 ./elpa2_kernels/elpa2_kernels_asm_x86_64.s
else
SOURCES = elpa1.F90 elpa2.F90 ./elpa2_kernels/elpa2_kernels_real_bgq.f90 ./elpa2_kernels/elpa2_kernels_complex.f90
endif

tmp = $(SOURCES:.F90=.o)
tmp2 = $(tmp:.f90=.o)
tmp3 = $(tmp2:.s=.o)
OBJECTS = $(tmp3:./elpa2_kernels/%=%)

LOC_LIBS = $(LIBS:-lelpa=)

BUILD_CMD = $(F90) $(FFLAGS) $< -c # $(INCPATH) -fPIC 

elpa1.o: elpa1.F90
	$(BUILD_CMD)
elpa2.o: elpa2.F90 elpa1.o
	$(BUILD_CMD)
elpa2_kernels_asm_x86_64.o: ./elpa2_kernels/elpa2_kernels_asm_x86_64.s
	$(BUILD_CMD)
elpa2_kernels_complex.o: ./elpa2_kernels/elpa2_kernels_complex.f90
	$(BUILD_CMD)
elpa2_kernels_complex_simple.o: ./elpa2_kernels/elpa2_kernels_complex_simple.f90
	$(BUILD_CMD)
elpa2_kernels_real.o: ./elpa2_kernels/elpa2_kernels_real.f90
	$(BUILD_CMD)
elpa2_kernels_real_bgq.o: ./elpa2_kernels/elpa2_kernels_real_bgq.f90
	$(BUILD_CMD)

all: elpa1.o elpa2.o elpa2_kernels.o
	$(F90) $^ -shared -Wl,-soname,libelpa.so $(LOC_LIBS) -o libelpa.so.2 
	ln -sf libelpa.so.2 libelpa.so

ifeq (,$(findstring -DSCALAPACK,$(FFLAGS)))
elpa2:
	@echo NOT BUILDING ELPA AS -DSCALAPACK NOT INCLUDED IN FFLAGS
	BUILD_ELPA=no
else
elpa2: $(OBJECTS)
	$(BUILD_CMD)
endif

clean:
	rm -f elpa*.o elpa*.mod
