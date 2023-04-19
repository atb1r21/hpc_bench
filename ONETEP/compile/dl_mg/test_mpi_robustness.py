#!/usr/env python3

# Run test_onetep_num over a range of MPI values to test robustness of 
# implementation.

# To do:
# * Implement automatic checking of numerical results for consistency between
#   runs with different Nmpi values.
#   - Do this by parsing output files, extracting numerical results, and 
#     comparing against the correct result / results from other runs.
# * Allow options to be read in from an input file (JSON?) and passed to the
#   Fortran program for increased configurability.


import argparse
import os
import subprocess
import sys

# Create argument parser (2 positional arguments)
parser = argparse.ArgumentParser(description=\
        "Run test_onetep_num over a range of MPI values to test robustness of implementation.")
parser.add_argument('testcase',action='store',help="Name of test case (see test_onetep_num source)")
parser.add_argument('outputbasename',action='store',help="Base name of output file (a file for each Nmpi value will be created with _[Nmpi]mpi.out appended)")
parser.add_argument('--output-dir',action='store',help="Directory in which to place output files (default is current directory)", default=os.getcwd() )
parser.add_argument('--debug',action='store_true',help="Do not write to files or "+\
                    "run commands. Instead print what would have happened to stdout")

# Create a Namespace object containing the argument value
arg_obj = parser.parse_args()

# Path to test executable
test_exe_cmd = "./test_onetep_num.exe"

# Command for running with MPI
mpirun_cmd = "mpirun"
# Option for MPI command which allows number of processes to be set.
# The script will pass this option followed by the integer number of 
# MPI processes set in the parallel_info.json file for each input
# directory
mpirun_proc_opt = "-np"

# Output file extension
output_file_ext = ".out"

# List of total MPI rank values to use
nmpi_list = [ 1, 4, 9, 11, 20, 50, 99, 128 ]

# Order for defect correction
order = 8

# Check existence of output directory
assert os.path.isdir( arg_obj.output_dir ), "The output directory "+arg_obj.output_dir+\
        " does not exist. Please create it or select an alternative directory."

for nmpi in nmpi_list:
    # Run test_onetep_num for user-provided test case and all values in nmpi_list
    # Build command to run
    run_cmd_list = [ mpirun_cmd, mpirun_proc_opt, str(nmpi), test_exe_cmd, \
                     arg_obj.testcase, "1", "1", str(nmpi), str(order) ]
    output_fname = arg_obj.outputbasename+"_"+str(nmpi)+"mpi"+output_file_ext
    output_fpath = os.path.join(arg_obj.output_dir,output_fname)
    if arg_obj.debug: 
        print( "Would output to:", output_fpath )
        print( "Would run:      "," ".join(run_cmd_list) )
    else:
        # Check existence of output file
        assert not os.path.isfile(output_fpath), "Output files cannot be overwritten. "+\
                "Please remove or backup "+output_fpath+" before proceeding."
        # Run command
        with open(output_fpath, "x") as fp:
            print( " ".join(run_cmd_list)+" > "+output_fpath )
            try:
                subprocess.check_call( run_cmd_list, stdout=fp )
            except subprocess.CalledProcessError as ex:
                # Allow calls to fail and cycle to next input file, but
                # print a warning
                print( "WARNING: The command "+" ".join(run_cmd_list)+\
                       " caused an error, with return code "+str(ex.returncode) )
                print( "... proceeding to next test." )




