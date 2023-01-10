#!/bin/sh

#SBATCH --ntasks-per-node=16     
#SBATCH --nodes=32
#SBATCH --exclusive
#SBATCH --reservation=benchmarks
#SBATCH --qos=benchmark
#SBATCH --partition=batch
#SBATCH --time=00:60:00
#SBATCH --ntasks=512
#SBATCH --output=benchio.512n.out
#SBATCH --error=benchio.512n.err
#SBATCH --chdir=/home/hpc/benchmarks-2022/benchio/benchio-master/file-per-process/source

module purge
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0

cd /scratch/hpc/2n

mpirun -n $SLURM_NPROCS /home/hpc/benchmarks-2022/benchio/benchio-master/file-per-process/source/benchio_fpp.x
