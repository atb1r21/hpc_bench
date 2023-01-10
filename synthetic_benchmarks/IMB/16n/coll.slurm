#!/bin/sh

#SBATCH --job-name=coll-imb
#SBATCH --ntasks-per-node=40
#SBATCH --reservation=root_143
#SBATCH --output=coll-imb-16n.out
#SBATCH --error=coll-imb-16n.err
#SBATCH --nodes=16
#SBATCH --ntasks=640
#SBATCH --exclusive
#SBATCH --chdir=/home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3/16n
#SBATCH --time=60:00:00
#SBATCH --qos=special

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0
export I_MPI_PIN_PROCESSOR_LIST="all:map=scatter"

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3/IMB-MPI1 \
       -npmin $SLURM_NTASKS -msglen length.txt -mem 4 Allreduce Allgather Scatter Alltoall
