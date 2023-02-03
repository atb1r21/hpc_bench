#!/bin/sh

#SBATCH --job-name=hpl_32n
#SBATCH --output=hpl_32n.out
#SBATCH --error=hpl_32n.err
#SBATCH --qos=benchmark
#SBATCH --reservation=benchmarks
# #SBATCH --partition=batch
#SBATCH --nodes=32
#SBATCH --ntasks=1280
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --chdir=/home/hpc/benchmarks-2022/HPL/hpl-2.3/bin

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0 intel-mkl/2021.2.0

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/HPL/hpl-2.3/bin/xhpl
