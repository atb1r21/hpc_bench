FILES := main.F90 ../utils/grid_data.F90
OBJ := $(notdir $(FILES:.F90=.o))

$(info fc is $(FC) dl_mg obj are  $(MG_SOLVER_OBJ) obj is $(OBJ))

unit_test.x : $(OBJ) $(MG_SOLVER_OBJ)
	$(FC) -o $@ -I$(OBJDIR) $(FFLAGS) $^

.PHONY: clean
clean :
	rm -fr $(OBJ)

%.o : %.F90
	${FC} -o $@ -c -I/$(OBJDIR) ${FFLAGS} $<

%.o : ../utils/%.F90
	${FC} -o $@ -c -I/$(OBJDIR) ${FFLAGS} $<

main.o : main.F90 grid_data.o
grid_data.o : ../utils/grid_data.F90
