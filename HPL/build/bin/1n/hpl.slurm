#!/bin/sh

#SBATCH --job-name=hpl_1n
#SBATCH --output=hpl_1n.out
#SBATCH --error=hpl_1n.err
#SBATCH --qos=special
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --exclusive
#SBATCH --chdir=/home/hpc/benchmarks-2022/HPL/hpl-2.3/bin/1n
#SBATCH --time=60:00:00

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0 intel-mkl/2021.2.0

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/HPL/hpl-2.3/bin/xhpl
