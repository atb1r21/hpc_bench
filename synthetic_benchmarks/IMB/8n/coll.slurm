#!/bin/sh

#SBATCH --job-name=coll-imb
#SBATCH --output=coll-imb-8n.out
#SBATCH --error=coll-imb-8n.err
#SBATCH --reservation=root_143
#SBATCH --nodes=8
#SBATCH --ntasks=320
#SBATCH --exclusive
#SBATCH --chdir=/home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3/8n
#SBATCH --time=60:00:00
#SBATCH --qos=special

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3/IMB-MPI1 \
       -npmin $SLURM_NTASKS -msglen length.txt -mem 4 Allreduce Allgather Scatter Alltoall
