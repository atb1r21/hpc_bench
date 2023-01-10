#!/bin/sh

#SBATCH --job-name=coll-imb
#SBATCH --ntasks-per-node=40
#SBATCH --reservation=benchmarks
#SBATCH --output=coll-imb-32n.out
#SBATCH --error=coll-imb-32n.err
#SBATCH --nodes=32
#SBATCH --ntasks=1280
#SBATCH --time=60:00:00

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0
export I_MPI_PIN_PROCESSOR_LIST="all:map=scatter"

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3/IMB-MPI1 \
       -npmin $SLURM_NTASKS -msglen length.txt -mem 5.1 Allreduce Allgather Scatter Alltoall
