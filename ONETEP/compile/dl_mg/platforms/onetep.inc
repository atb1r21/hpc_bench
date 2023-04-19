# wrapper for ONETEP config files
#
#Lucian Anton, May 2013
#Updated by James Womack, 2017

# JCW: Set ROOTDIR to ONETEP_ROOT, so conf file can
# JCW: use this as it would when included in standard
# JCW: ONETEP Makefile
ROOTDIR = $(ONETEP_ROOT)

include $(ONETEP_ROOT)/config/conf.$(ONETEP_ARCH)
# JCW: ONETEP conf file brings $(COMPILER), $(DEBUGFLAGS), $(FFLAGS)
# JCW: $(OPTFLAGS), $(LIBS).


FC := $(F90)

# JCW: To compile DL_MG object files and module we need to add the $(OBJDIR)
# JCW: and $(LIBDIR) (should be passed to the call to make from ONETEP's
# JCW: Makefile) to include path in $(FFLAGS).
FFLAGS +=  -I$(OBJDIR) -I$(LIBDIR)

# JCW: We also need to inform the compiler that module files should be placed
# JCW: in $(OBJDIR). This has a compiler-specific flag.
ifeq ($(COMPILER),INTEL-ifort-on-LINUX)
  FFLAGS += -module $(OBJDIR)
else ifeq ($(COMPILER),NAG-f95-on-LINUX)
  FFLAGS += -mdir $(OBJDIR)
else
# Default to the GNU (and CRAY) syntax (may not work for all compilers)
  FFLAGS += -J$(OBJDIR)
endif

# JCW: Add $(DEBUGFLAGS) or $(OPTFLAGS) to $(FFLAGS) as appropriate based on
# JCW: $(ONETEP_TARGET), which is passed to DL_MG's invocation of make by ONETEP's
# JCW: Makefile multigrid target.
ifeq ($(ONETEP_TARGET),debug)
  FFLAGS += $(DEBUGFLAGS)
else
  FFLAGS += $(OPTFLAGS) $(OPENMPFLAGS)
endif

# If compiler produces uppercase modules, copy to lowercase in order that ONETEP's
# dependencies (on dl_mg.mod and dl_mg_nonlin_model.mod) are correctly recognized
ifeq ($(UPPERCASE_MODULE),yes)
  CP_MOD = UPPERMOD=$(dir $@)$$(echo $(basename $(notdir $@)) | tr '[:lower:]' '[:upper:]'); \
  for EXT in "mod" "MOD"; do \
    if [[ -e $${UPPERMOD}.$${EXT} ]]; then \
      UPPERMOD=$${UPPERMOD}.$${EXT}; cp $${UPPERMOD} $(basename $@).mod; \
    fi; \
  done
else
  CP_MOD =
endif


