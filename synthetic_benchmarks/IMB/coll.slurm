#!/bin/sh

#SBATCH --job-name=coll-imb
#SBATCH --output=coll-imb.out
#SBATCH --error=coll-imb.err
#SBATCH --nodes=32
#SBATCH --ntasks=1280
#SBATCH --exclusive
#SBATCH --chdir=/home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3
#SBATCH --time=60:00:00
#SBATCH --qos=special

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3/IMB-MPI1 Allreduce Allgather Scatter Alltoall
