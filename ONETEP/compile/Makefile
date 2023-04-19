#
# Top-level makefile for ONETEP
#
# This version for multiple architectures by Peter Haynes
#
# System-specific changes should be made to the appropriate
#  conf.arch files
#

ROOTDIR = $(PWD)

ARCH = `$(ROOTDIR)/bin/arch`

MULTIGRID = dl_mg

BUILD_ELPA = no

TARGET ?= onetep

.PHONY: default onetep debug timing profile warnings utils clean cleanall \
        cleanutils dist help checkdeps checkesdf checkuses objdir libdir \
        utildir multigrid cleanmultigrid elpa2 cleanelpa

default: help

onetep: objdir libdir multigrid elpa2
	@( cd obj/$(ARCH) ; \
	$(MAKE) -f $(ROOTDIR)/src/Makefile $(TARGET) VPATH=$(ROOTDIR)/src:$(VPATH) \
		ARCH=$(ARCH) ROOTDIR=$(ROOTDIR) MULTIGRID=$(MULTIGRID) \
		BUILD_ELPA=$(BUILD_ELPA))
	@( echo; echo ONETEP compiled and linked successfully.; echo If this \
	is the first time you have built this version on this system,; echo make \
	sure to run the quality-check \(QC\) tests -- see the tests directory.)

libonetep: onetep
	@( cd obj/$(ARCH) ; \
        ar src libonetep.a *.o)

python: libonetep
	@( cd obj/$(ARCH) ; \
	$(MAKE) -f $(ROOTDIR)/src/Makefile python VPATH=$(ROOTDIR)/src:$(VPATH) \
		ARCH=$(ARCH) ROOTDIR=$(ROOTDIR) MULTIGRID=$(MULTIGRID) \
		BUILD_ELPA=$(BUILD_ELPA))

debug:
	$(MAKE) onetep TARGET=debug

timing:
	$(MAKE) onetep TARGET=timing

profile:
	$(MAKE) onetep TARGET=profile

warnings:
	$(MAKE) onetep TARGET=warnings

utils: utildir
	@( cd utils ; \
	$(MAKE) -f Makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR) EXTFLAGS=$(EXTFLAGS) )

clean:
	@$(RM) -fr obj/$(ARCH)

cleanall: cleanutils cleanmultigrid
	@$(RM) -fr obj/* bin/onetep.* *~ src/*~

cleanutils:
	@(cd utils ; \
	$(MAKE) -f Makefile clean ARCH=$(ARCH) ROOTDIR=$(ROOTDIR) )

dist:
	@( date > Archived.txt ; \
	tar -cf onetep.tar Archived.txt Makefile bin config include lib \
        obj src utils INSTALL )

help:
	@echo
	@echo '        ####### #     # ####### ####### ####### ######   '
	@echo '        #     # ##    # #          #    #       #     #  '
	@echo '        #     # # #   # #          #    #       #     #  '
	@echo '        #     # #  #  # #####      #    #####   ######   '
	@echo '        #     # #   # # #          #    #       #        '
	@echo '        #     # #    ## #          #    #       #        '
	@echo '        ####### #     # #######    #    ####### #        '
	@echo
	@echo 'To install type:'
	@echo '> gmake onetep'
	@echo
	@echo 'Other options:'
	@echo 'To compile a version for debugging:'
	@echo '> gmake clean ; gmake debug'
	@echo 'To compile a version for profiling:'
	@echo '> gmake clean ; gmake profile'
	@echo 'To remove object files for this architecture:'
	@echo '> gmake clean'
	@echo 'To remove object files for all architectures:'
	@echo '> gmake cleanall'
	@echo 'To create an archive for distribution:'
	@echo '> gmake dist'
	@echo 'To make utilities associated with ONETEP:'
	@echo '> gmake utils'
	@echo

checkdeps:
	utils/check_dependencies

checkesdf:
	utils/check_esdf

checkuses:
	utils/check_use_clauses

objdir:
	@( if [ ! -d obj ] ; \
	        then mkdir obj ; \
	fi )
	@( if [ ! -d obj/$(ARCH) ] ; \
		then mkdir obj/$(ARCH) ; \
	fi )

libdir:
	@( if [ ! -d lib ] ; \
	        then mkdir lib ; \
	fi )
	@( if [ ! -d lib/$(ARCH) ] ; \
	        then mkdir lib/$(ARCH) ; \
	fi )


utildir:
	@( if [ ! -d utils ] ; \
		then echo 'utils directory not found!' ; exit 1 ; \
	fi )

ifeq ($(strip $(MULTIGRID)),no)

multigrid: ;
cleanmultigrid: ;

else ifeq ($(strip $(MULTIGRID)),dl_mg)

multigrid : objdir
	@( cd dl_mg ; \
	$(MAKE) -f $(ROOTDIR)/dl_mg/Makefile ONETEP_ROOT=$(ROOTDIR) \
	OBJDIR=$(ROOTDIR)/obj/$(ARCH) LIBDIR=$(ROOTDIR)/obj/$(ARCH) \
	ONETEP_ARCH=$(ARCH)  ONETEP_TARGET=$(TARGET) \
	PLATFORM=onetep obj_for_onetep)
cleanmultigrid :
	@( if [ -d obj/$(ARCH) ] ; then \
		cd obj/$(ARCH) ; \
		$(MAKE) -f $(ROOTDIR)/dl_mg/Makefile ONETEP_ROOT=$(ROOTDIR) MG_ROOTDIR=$(ROOTDIR)/dl_mg \
			ONETEP_ARCH=$(ARCH) PLATFORM=onetep ONETEP_TARGET=$(TARGET) clean_onetep ; \
	fi )

else 

$(error unknown value for MULTIGRID, try "no" or "dl_mg")

endif

ifeq ($(strip $(BUILD_ELPA)),no)

elpa2: ;
cleanelpa: ;

else ifeq ($(strip $(BUILD_ELPA)),yes)

elpa2 : objdir
	@ cd obj/$(ARCH); $(MAKE) -f $(ROOTDIR)/elpa/Makefile ONETEP_ROOT=$(ROOTDIR) ONETEP_ARCH=$(ARCH) VPATH=$(ROOTDIR)/elpa:$(VPATH) elpa2
cleanelpa :
	@ cd obj/$(ARCH); $(MAKE) -f $(ROOTDIR)/elpa/Makefile ONETEP_ROOT=$(ROOTDIR) ONETEP_ARCH=$(ARCH) VPATH=$(ROOTDIR)/elpa:$(VPATH) clean

else 

$(error boolean BUILD_ELPA should be either "yes" or "no".)

endif
