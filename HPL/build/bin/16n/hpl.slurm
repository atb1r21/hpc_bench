#!/bin/sh

#SBATCH --job-name=hpl_16n
#SBATCH --output=hpl_16n.out
#SBATCH --error=hpl_16n.err
#SBATCH --qos=special
#SBATCH --partition=batch
#SBATCH --nodes=16
#SBATCH --ntasks=640
#SBATCH --exclusive
#SBATCH --time=60:00:00

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0 intel-mkl/2021.2.0

cd /home/hpc/benchmarks-2022/HPL/hpl-2.3/bin/16n
mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/HPL/hpl-2.3/bin/xhpl
