# List of available tasks.
.PHONY: tests serial_tests clean help

##########################
# tests and serial_tests #
##########################
# Just a helper around testcode.
# For more advanced functionality, use testcode.py directly.

tests:
ifeq ($(ARCH),)
	$(warning Please set ARCH.)
	@$(MAKE) help
else
	ARCH=$(ARCH) testcode/bin/testcode.py $(TC) --verbose
endif

serial_tests:
	$(MAKE) tests TC="-p0 $(TC)"

#########
# clean #
#########
# Remove all files apart from the input and benchmark files.
# Leave symlinks behind (used for running tests based upon output from an
# earlier test).
# This does more than testcode's tidy option, which only removes test output
# and error files.

# Note: in this section, the make macros KEEP_TESTS and CLEAN_TESTS are named
# after the concept of tests in the find utility, not in relation to
# tests in testcode.

# jme: Filenames (and globs) matching files to be preserved by "make clean"
KEEP_FILES := *.dat *.sw_ngwfs *.md *.cdft *.eda* \$$Q* \$$M* QM_*.vmatrix \
    QM_*.mpoles_for_onetep
# jme: Add "! -name" prefix to each filename/glob in KEEP_FILES to create
#      a list of (AND NOT)-separated tests for use by the find command
KEEP_TESTS := $(patsubst %,! -name %,${KEEP_FILES})

# JCW: Filenames (and globs) matching files to be removed by "make clean"
# jme: This is the list of files that, despite matching a pattern
# in KEEP_FILES, we want removed.
CLEAN_FILES := test.out.* test.err.* *inittr* spectrum.dat \
    npa_n*.dat pdos_weights.dat c2h4-group_spin_diff.cdft \
    cedft-Pt13.cdft H2O_dimer_md*.md
# JCW: Add "-o -name" prefix to each filename/glob in CLEAN_FILES to create
# JCW: a list of OR-separated tests for use by the find command
CLEAN_TESTS := $(patsubst %,-o -name %,${CLEAN_FILES})

clean:
	find test?? test_long?? \
	       -type f ${KEEP_TESTS} \
	               ! -path '*/.svn*' \
	       ${CLEAN_TESTS} \
	| xargs -i rm -v {}

########
# help #
########

help:
	@echo ""
	@echo "ONETEP QC make help"
	@echo ""
	@echo "To run all tests execute"
	@echo "> make ARCH=XXX"
	@echo "where config/conf.XXX is the config file used to compile ONETEP."
	@echo ""
	@echo "A serial version is also available:"
	@echo "> make serial_tests ARCH=XXX"
	@echo ""
	@echo "For more flexible testing, use the testcode script (testcode/bin/testcode.py) directly."
	@echo "Additional flags can be passed to testcode by setting the TC variable."
	@echo ""
	@echo "To delete all outputs in the test subdirectories execute"
	@echo "> make clean"
	@echo "(Warning: this may also delete other files you created in those subdirectories.)"
	@echo ""

