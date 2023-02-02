#!/bin/bash

# -------------------------------------------------
# A SLURM submission script for ONETEP on IRIDIS5.
# Supports hybrid (MPI/OMP) parallelism.
# 2021.11 Jacek Dziedzic, J.Dziedzic@soton.ac.uk
#                         University of Southampton
# -------------------------------------------------

# v1.50 (2018.01.23)
# v1.61 (2020.04.22) Fixes the location of I_MPI_PMI_LIBRARY. Needed after Iridis5 scheduler upgrade 2020.04.21.
# v2.00 (2021.03.06) Clean-up and documentation
# v2.10 (2021.11.16) Switch to Intel 2021.

# ==========================================================================================================
# Edit the following lines to your liking.
#
#SBATCH -J onetep                 # Name of the job.
#SBATCH --reservation=matlab                  # Queue. Use 'batch' for most jobs, or 'scavenger' for quick tests.
#SBATCH --qos=benchmark
#SBATCH --ntasks=10                   # Total number of MPI processes in job.
#SBATCH --nodes=1                     # Number of nodes in job.
#SBATCH --ntasks-per-node 10          # Number of MPI processes per node.
#SBATCH --cpus-per-task 4             # Number of OMP threads spawned from each MPI process.
#SBATCH --time 24:00:00                # Max time for your job (hh:mm:ss).


module load intel-compilers/2021.2.0
module load intel-mkl/2021.2.0
module load intel-mpi/2021.2.0

omp_threads_per_mpi_rank=4            # Repeat the value from 'cpus-per-task' here.

# Point this to your ONETEP executable.
onetep_exe=\
"/home/hpc/benchmarks-2022/ONETEP/source/bin/onetep.iridis5.intel21.omp.scalapack"

# Point this to your ONETEP launcher.
onetep_launcher=\
"/home/hpc/benchmarks-2022/ONETEP/source/utils/onetep_launcher"
# ==========================================================================================================

# You should not need to edit anything below this line.

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------

workdir=`pwd`
echo "--- This is the submission script, the time is `date`."
echo "--- ONETEP executable is '$onetep_exe'."
echo "--- workdir is '$workdir'."
echo "--- onetep_launcher is '$onetep_launcher'."

# Ensure exactly 1 .dat file in there.
ndats=`ls -l *dat | wc -l`

if [ "$ndats" == "0" ]; then
  echo "!!! There is no .dat file in the current directory. Aborting." >&2
  touch "%NO_DAT_FILE"
  exit 2
fi

if [ "$ndats" == "1" ]; then
  true
else
  echo "!!! More than one .dat file in the current directory, that's too many. Aborting." >&2
  touch "%MORE_THAN_ONE_DAT_FILE"
  exit 3
fi

rootname=`echo *.dat | sed -r "s/\.dat\$//"`
rootname_dat=$rootname".dat"
rootname_out=$rootname".out"
rootname_err=$rootname".err"

echo "--- The input file is $rootname_dat, the output goes to $rootname_out and errors go to $rootname_err."

# Ensure ONETEP executable is there and is indeed executable.
if [ ! -x "$onetep_exe" ]; then
  echo "!!! $onetep_exe does not exist or is not executable. Aborting!" >&2
  touch "%ONETEP_EXE_MISSING"
  exit 4
fi

# Ensure onetep_launcher is there and is indeed executable.
if [ ! -x "$onetep_launcher" ]; then
  echo "!!! $onetep_launcher does not exist or is not executable. Aborting!" >&2
  touch "%ONETEP_LAUNCHER_MISSING"
  exit 5
fi

# Dump the module list to a file.
module list 2>\$modules_loaded

ldd $onetep_exe >\$ldd

# Report details
echo "--- Number of nodes as reported by SLURM: $SLURM_JOB_NUM_NODES."
echo "--- Number of tasks as reported by SLURM: $SLURM_NTASKS."
echo "--- Using this srun executable: "`which srun`
echo "--- Executing ONETEP via $onetep_launcher."
 
# Need to set this after Iridis5 scheduler upgrade
export I_MPI_PMI_LIBRARY=/local/software/slurm/default/lib/libpmi.so

# Actually run ONETEP
########################################################################################################################################################
srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS $onetep_launcher -e $onetep_exe -t $omp_threads_per_mpi_rank $rootname_dat >$rootname_out 2>$rootname_err
########################################################################################################################################################

echo "--- srun finished at `date`."

# Check for error conditions
result=$?
if [ $result -ne 0 ]; then
  echo "!!! srun reported a non-zero exit code $result. Aborting!" >&2
  touch "%SRUN_ERROR"
  exit 6
fi

if [ -r $rootname.error_message ]; then
  echo "!!! ONETEP left an error message file. Aborting!" >&2
  touch "%ONETEP_ERROR_DETECTED"
  exit 7
fi

tail $rootname.out | grep completed >/dev/null 2>/dev/null
result=$?
if [ $result -ne 0 ]; then
  echo "!!! ONETEP calculation likely did not complete. Aborting!" >&2
  touch "%ONETEP_DID_NOT_COMPLETE"
  exit 8
fi

echo "--- Looks like everything went fine. Praise be."
touch "%DONE"

echo "--- Finished successfully at `date`."
