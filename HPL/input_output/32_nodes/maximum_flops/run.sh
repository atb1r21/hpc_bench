#!/bin/sh

#SBATCH --job-name=hpl32
#SBATCH --partition=batch
#SBATCH --reservation=matlab
#SBATCH --nodes=32
#SBATCH --ntasks=1280
#SBATCH --time=60:00:00

module purge
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0 intel-mkl/2021.2.0

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/HPL/hpl-2.3/bin/xhpl

