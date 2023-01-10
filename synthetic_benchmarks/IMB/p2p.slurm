#!/bin/sh

#SBATCH --job-name=p2p-imb
#SBATCH --time=00:20:00
#SBATCH --output=p2p-imb.out
#SBATCH --error=p2p-imb.err
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive

module purge 
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0

mpirun -n $SLURM_NTASKS /home/hpc/benchmarks-2022/IMB/mpi-benchmarks-IMB-v2021.3/IMB-MPI1 PingPong PingPing
