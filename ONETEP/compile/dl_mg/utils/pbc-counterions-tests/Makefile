FILES := main.F90 modules.F90 ../mg_utils.F90
OBJ := $(subst ../,, $(FILES:.F90=.o))

$(info fc is $(FC) libdir is $(LIBDIR) obj is $(OBJ))

cntion_test.x : $(OBJ) $(MG_ROOTDIR)/$(LIBDIR)/libdlmg.a
	$(FC) -o $@ -I$(MG_ROOTDIR)/$(LIBDIR) $(FFLAGS) $^ $(MG_ROOTDIR)/$(LIBDIR)/libdlmg.a

.PHONY: clean
clean :
	rm -fr $(OBJ)

%.o : %.F90
	${FC} -o $@ -c -I/$(MG_ROOTDIR)/$(LIBDIR) ${FFLAGS} $<

%.o : ../%.F90
	${FC} -o $@ -c -I/$(MG_ROOTDIR)/$(LIBDIR) ${FFLAGS} $<

main.o : main.F90 modules.o mg_utils.o
modules.o : modules.F90
mg_utils.o : ../mg_utils.F90
