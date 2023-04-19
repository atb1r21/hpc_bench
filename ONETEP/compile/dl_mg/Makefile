THIS_MAKEFILE := $(lastword $(MAKEFILE_LIST))

MG_ROOTDIR = $(CURDIR)

ifndef EXE
  EXE := dl_mg_test.exe
endif

OBJDIR ?= $(MG_ROOTDIR)/o
LIBDIR ?= lib

vpath %.F90 src:utils
vpath %.inc src:utils

ifdef PLATFORM
  ifeq (,$(wildcard $(MG_ROOTDIR)/platforms/$(PLATFORM).inc))
    $(error platforms/$(PLATFORM).inc does not exist !!!)
  endif
  include $(MG_ROOTDIR)/platforms/$(PLATFORM).inc
endif

# CP_MOD is run for every compilation of an object file and allows module files
# to be copied to alternative locations / with different filenames, if required.
# This is required in ONETEP, which relies upon lower case module names. CP_MOD
# should be set on a per-platform basis, and if CP_MOD is not set in the
# platform file, set it to empty.
CP_MOD ?=
# Note: CP_MOD is executed only if it is non-empty and command output is suppressed
# 			so will not appear in the output of make.

MG_SOLVER_F90 := dl_mg_params.F90 dl_mg_types.F90 dl_mg_common_data.F90 dl_mg_mpi_header.F90 dl_mg_timer.F90 dl_mg_grids.F90 dl_mg_mpi.F90  dl_mg_info.F90 dl_mg_stencil_coefficients.F90 dl_mg_kernels.F90 dl_mg_multigrid_method.F90 dl_mg.F90 dl_mg_nonlinear_model.F90 dl_mg_errors.F90 dl_mg_defco_mod.F90 dl_mg_defco_fd_mod.F90 dl_mg_defco_utils_mod.F90 dl_mg_convergence_params.F90 dl_mg_utils.F90 dl_mg_neutralisation_with_ions.F90 dl_mg_newton_method.F90 dl_mg_second_order_solvers.F90 dl_mg_conjugate_gradient_method.F90 dl_mg_alloc.F90
MG_SOLVER_INC := dl_mg_common_params.inc

TEST_F90 := mg_utils.F90 mg_tests.F90 main_mg.F90


MG_SOLVER_OBJ := $(addprefix $(OBJDIR)/, $(MG_SOLVER_F90:.F90=.o))
TEST_OBJ      := $(addprefix $(OBJDIR)/, $(TEST_F90:.F90=.o))

MG_SOLVER_INC := $(addprefix $(OBJDIR)/, $(MG_SOLVER_INC))

default : $(LIBDIR)/libdlmg.a

$(LIBDIR)/libdlmg.a : checkplatform $(MG_SOLVER_OBJ) $(MG_SOLVER_INC) | $(LIBDIR)
	$(PLATFORM_ERR) ar -cr $(LIBDIR)/libdlmg.a $(MG_SOLVER_OBJ) && \
        if [ -f $(OBJDIR)/dl_mg.mod ] ; then mv $(OBJDIR)/dl_mg.mod $(LIBDIR) ; fi

# build tests program
test :  checkplatform $(TEST_OBJ) $(MG_SOLVER_OBJ) $(MG_SOLVER_INC)
	$(LOAD_MODULES) $(FC) -o $(EXE) $(INCPATHS) $(FFLAGS) $(TEST_OBJ) $(MG_SOLVER_OBJ)

testlib : checkplatform $(LIBDIR)/libdlmg.a $(TEST_OBJ)
	$(LOAD_MODULES) $(FC) -o $(EXE) $(INCPATHS) $(FFLAGS) $(TEST_OBJ) -I$(LIBDIR) $(LIBDIR)/libdlmg.a

# test for onetep input data
test_onetep_ls : checkplatform $(OBJDIR)/test_onetep_ls.o $(MG_SOLVER_OBJ) $(MG_SOLVER_INC)
	$(LOAD_MODULES) $(FC) -o $@.exe $(INCPATHS) $(FFLAGS) $(OBJDIR)/test_onetep_ls.o $(MG_SOLVER_OBJ)

# test for onetep input data (numerical, defect correction)
test_onetep_num : checkplatform $(OBJDIR)/test_onetep_num.o $(MG_SOLVER_OBJ) $(MG_SOLVER_INC)
	$(LOAD_MODULES) $(FC) -o $@.exe $(INCPATHS) $(FFLAGS) $(OBJDIR)/test_onetep_num.o $(MG_SOLVER_OBJ)

obj_for_onetep : checkplatform $(MG_SOLVER_OBJ) $(MG_SOLVER_INC)

.PHONY : clean vclean clean_onetep obj_for_onetep doc docdev checkplatform release doc docdev

DOCTYPE := _user
VERSION_OPT := -no-hash
docdev : DOCTYPE = _developer
docdev : VERSION_OPT=-with-date
docdev : doc
doc :  VERSION_OPT=-with-date
doc:
	$$(sh utils/get_version_string.sh $(VERSION_OPT) > VERSION)
	doxygen $(MG_ROOTDIR)/docs/doxygen_doc/Doxyfile$(DOCTYPE)

release: doc
	if [ $$(git rev-parse --abbrev-ref HEAD) != master ] ; then \
          echo "release task should be used only on master branch. Quitting ..."; exit 3 ;\
        else cd docs/doxygen_doc/latex && $(MAKE) clean && $(MAKE) && \
          mv refman.pdf $(MG_ROOTDIR)/docs/dl_mg.pdf && cd $(MG_ROOTDIR) && \
          vs=$$(sh utils/get_version_string.sh) && \
          git archive -o dl_mg_$${vs}.tar --prefix=dl_mg_$${vs}/ HEAD && \
          tar --xform='s,^,dl_mg_'$${vs}'/,' -r -f dl_mg_$${vs}.tar docs/dl_mg.pdf ; \
        fi


export FC FFLAGS MG_ROOTDIR LIBDIR
cntion_test: checkplatform $(OBJDIR)/dl_mg_utils.o
	cd utils/pbc-counterions-tests && $(MAKE)

export MG_SOLVER_OBJ OBJDIR
unit_test: checkplatform $(MG_SOLVER_OBJ)
	cd src/tests/unit_tests && $(MAKE)

clean_onetep:
	@rm -f dl_mg*.o dl_mg*.mod

clean :
	@rm -fr $(OBJDIR)
vclean : clean
	@rm -fr *.exe $(LIBDIR)

checkplatform :
	@if [ -z "$(PLATFORM)" ] ; then echo "PLATFORM needs to be defined! See ./platforms folder and the build documetation"; exit 2 ; fi

$(OBJDIR)/%.o: %.F90 | $(MG_SOLVER_INC)
	$(LOAD_MODULES) ${FC} -o $@ -c $(INCPATHS) ${FFLAGS} $<
ifneq ($(strip $(CP_MOD)),)
	@($(CP_MOD))
endif

$(OBJDIR)/%.inc: %.inc
	cp $< $@

$(MG_SOLVER_OBJ) $(TEST_OBJ) $(MG_SOLVER_INC) : | $(OBJDIR)

$(OBJDIR) :
	mkdir -p $(OBJDIR)

# Avoid warnings about overriding existing recipes
ifneq ($(strip $(OBJDIR)),$(strip $(LIBDIR)))
$(LIBDIR) :
	mkdir -p $(LIBDIR)
else
$(LIBDIR) :
endif

$(TEST_OBJ) : INCPATHS += -I$(LIBDIR)

$(OBJDIR)/dl_mg.o: dl_mg.F90 dl_mg_common_params.inc \
            $(addprefix $(OBJDIR)/, \
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_grids.o \
            dl_mg_info.o \
            dl_mg_mpi.o \
            dl_mg_mpi_header.o \
            dl_mg_nonlinear_model.o \
            dl_mg_params.o \
            dl_mg_stencil_coefficients.o \
            dl_mg_timer.o \
            dl_mg_types.o \
            dl_mg_defco_mod.o \
            dl_mg_defco_fd_mod.o \
            dl_mg_convergence_params.o\
            dl_mg_utils.o\
            dl_mg_conjugate_gradient_method.o)

$(OBJDIR)/dl_mg_params.o: dl_mg_params.F90 \
            dl_mg_common_params.inc \
            dl_mg_neutralisation_with_ions_params.inc

$(OBJDIR)/dl_mg_types.o: dl_mg_types.F90 $(OBJDIR)/dl_mg_params.o

$(OBJDIR)/dl_mg_common_data.o: dl_mg_common_data.F90 dl_mg_common_params.inc \
                              $(addprefix $(OBJDIR)/, dl_mg_params.o dl_mg_types.o )

$(OBJDIR)/dl_mg_mpi.o: dl_mg_mpi.F90 dl_mg_common_params.inc \
            $(addprefix $(OBJDIR)/, \
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_mpi_header.o \
            dl_mg_params.o \
            dl_mg_types.o )
ifdef DUMP_ARRAYS
$(OBJDIR)/dl_mg_mpi.o: $(OBJDIR)/mg_utils.o
endif

$(OBJDIR)/dl_mg_grids.o: dl_mg_grids.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_info.o \
            dl_mg_kernels.o \
            dl_mg_mpi.o \
            dl_mg_mpi_header.o \
            dl_mg_params.o \
            dl_mg_types.o \
            dl_mg_defco_utils_mod.o)

$(OBJDIR)/dl_mg_info.o: dl_mg_info.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_mpi_header.o \
            dl_mg_nonlinear_model.o \
            dl_mg_params.o \
            dl_mg_mpi.o \
            dl_mg_types.o \
            dl_mg_neutralisation_with_ions.o \
            dl_mg_convergence_params.o)

$(OBJDIR)/dl_mg_stencil_coefficients.o: dl_mg_stencil_coefficients.F90 \
            $(addprefix $(OBJDIR)/, \
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_mpi.o \
            dl_mg_params.o \
            dl_mg_types.o )

$(OBJDIR)/dl_mg_multigrid_method.o: dl_mg_multigrid_method.F90 \
          $(addprefix $(OBJDIR)/, \
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_info.o \
            dl_mg_kernels.o \
            dl_mg_params.o \
            dl_mg_mpi.o \
            dl_mg_mpi_header.o \
            dl_mg_stencil_coefficients.o \
            dl_mg_timer.o \
            dl_mg_types.o \
            dl_mg_utils.o \
            dl_mg_convergence_params.o \
            dl_mg_nonlinear_model.o \
            dl_mg_utils.o \
            dl_mg_neutralisation_with_ions.o)

$(OBJDIR)/dl_mg_timer.o : dl_mg_timer.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_common_data.o \
            dl_mg_grids.o \
            dl_mg_mpi_header.o )

$(OBJDIR)/dl_mg_kernels.o : dl_mg_kernels.F90 \
            implementations/blocked_dotproduct_2.inc \
            implementations/blocked_dotproduct_2.inc \
            $(addprefix $(OBJDIR)/,\
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_mpi_header.o \
            dl_mg_nonlinear_model.o \
            dl_mg_params.o \
            dl_mg_types.o \
            dl_mg_alloc.o)

$(OBJDIR)/dl_mg_nonlinear_model.o : dl_mg_nonlinear_model.F90 \
            $(addprefix $(OBJDIR)/, \
            dl_mg_params.o dl_mg_errors.o )

$(OBJDIR)/dl_mg_errors.o : dl_mg_errors.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_common_data.o \
            dl_mg_mpi_header.o)

$(OBJDIR)/dl_mg_defco_mod.o : dl_mg_defco_mod.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_common_data.o \
            dl_mg_errors.o \
            dl_mg_grids.o \
            dl_mg_info.o \
            dl_mg_mpi.o \
            dl_mg_mpi_header.o \
            dl_mg_kernels.o \
            dl_mg_nonlinear_model.o \
            dl_mg_params.o \
            dl_mg_types.o \
            dl_mg_timer.o \
            dl_mg_defco_fd_mod.o \
            dl_mg_utils.o \
            dl_mg_defco_utils_mod.o \
            dl_mg_convergence_params.o \
            dl_mg_neutralisation_with_ions.o \
            dl_mg_second_order_solvers.o)

$(OBJDIR)/dl_mg_defco_fd_mod.o : dl_mg_defco_mod.F90 $(addprefix $(OBJDIR)/, \
            dl_mg_common_data.o \
            dl_mg_grids.o \
            dl_mg_mpi.o \
            dl_mg_mpi_header.o \
            dl_mg_params.o \
	    dl_mg_timer.o \
            dl_mg_defco_utils_mod.o )

$(OBJDIR)/dl_mg_defco_utils_mod.o : dl_mg_defco_utils_mod.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_errors.o \
            dl_mg_mpi_header.o \
            dl_mg_params.o)

$(OBJDIR)/dl_mg_convergence_params.o : dl_mg_convergence_params.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_params.o \
            dl_mg_errors.o \
            dl_mg_common_data.o \
            dl_mg_mpi.o)

$(OBJDIR)/dl_mg_utils.o : dl_mg_utils.F90 $(addprefix $(OBJDIR)/,\
            dl_mg_info.o \
            dl_mg_params.o \
            dl_mg_mpi.o \
            dl_mg_grids.o \
            dl_mg_common_data.o \
            dl_mg_timer.o \
            dl_mg_types.o \
            dl_mg_nonlinear_model.o \
            dl_mg_info.o \
            dl_mg_convergence_params.o \
            dl_mg_neutralisation_with_ions.o )

$(OBJDIR)/dl_mg_neutralisation_with_ions.o : dl_mg_neutralisation_with_ions.F90 \
            $(addprefix $(OBJDIR)/,\
            dl_mg_params.o \
            dl_mg_common_data.o \
            dl_mg_nonlinear_model.o \
            dl_mg_errors.o)

$(OBJDIR)/dl_mg_second_order_solvers.o: dl_mg_second_order_solvers.F90 \
            $(addprefix $(OBJDIR)/, \
              dl_mg_mpi_header.o \
              dl_mg_types.o \
              dl_mg_common_data.o \
              dl_mg_grids.o \
              dl_mg_timer.o \
              dl_mg_errors.o \
              dl_mg_params.o \
              dl_mg_stencil_coefficients.o \
              dl_mg_multigrid_method.o \
              dl_mg_newton_method.o \
              dl_mg_conjugate_gradient_method.o)


$(OBJDIR)/dl_mg_newton_method.o: dl_mg_newton_method.F90 \
            $(addprefix $(OBJDIR)/,\
              dl_mg_mpi.o \
              dl_mg_common_data.o \
              dl_mg_types.o \
              dl_mg_kernels.o \
              dl_mg_errors.o \
              dl_mg_mpi_header.o \
              dl_mg_convergence_params.o \
              dl_mg_utils.o \
              dl_mg_nonlinear_model.o \
              dl_mg_utils.o \
              dl_mg_mpi.o \
              dl_mg_conjugate_gradient_method.o)

$(OBJDIR)/dl_mg_conjugate_gradient_method.o: dl_mg_conjugate_gradient_method.F90 \
            $(addprefix $(OBJDIR)/,\
              dl_mg_params.o \
              dl_mg_types.o \
              dl_mg_utils.o \
              dl_mg_kernels.o \
              dl_mg_multigrid_method.o \
              dl_mg_mpi.o\
              dl_mg_errors.o)

$(OBJDIR)/dl_mg_alloc.o: dl_mg_alloc.F90 \
            $(addprefix $(OBJDIR)/,\
              dl_mg_mpi.o \
              dl_mg_params.o \
              dl_mg_types.o \
              dl_mg_errors.o)

#dependencies for the internal tests
$(OBJDIR)/main_mg.o : main_mg.F90 $(addprefix $(OBJDIR)/, \
            dl_mg.o \
            dl_mg_mpi.o \
            dl_mg_params.o \
            mg_tests.o \
            mg_utils.o )
$(OBJDIR)/mg_tests.o : mg_tests.F90 $(addprefix $(OBJDIR)/, \
            dl_mg.o \
            dl_mg_params.o \
            dl_mg_mpi_header.o \
            mg_utils.o \
            dl_mg_nonlinear_model.o )

$(OBJDIR)/mg_utils.o : mg_utils.F90 $(addprefix $(OBJDIR)/, dl_mg_mpi_header.o)
$(OBJDIR)/test_onetep_ls.o : test_onetep_ls.F90 $(addprefix $(OBJDIR)/, dl_mg.o)
$(OBJDIR)/test_onetep_num.o : test_onetep_num.F90 $(addprefix $(OBJDIR)/, dl_mg.o dl_mg_params.o)


