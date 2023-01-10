#!/bin/sh

#SBATCH --ntasks-per-node=16     
#SBATCH --nodes=32
#SBATCH --exclusive
#SBATCH --qos=benchmark
#SBATCH --reservation=benchmarks
#SBATCH --partition=batch
#SBATCH --time=01:00:00
#SBATCH --ntasks=512
#SBATCH --output=benchio.512.out
#SBATCH --error=benchio.512.out

module purge
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0

cd /scratch/hpc

mpirun -n $SLURM_NPROCS /home/hpc/benchmarks-2022/benchio/benchio-master/shared-file/source/benchio.x
