#!/bin/sh

#SBATCH --ntasks-per-node=16     
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --qos=special
#SBATCH --partition=batch
#SBATCH --time=00:60:00
#SBATCH --ntasks=16
#SBATCH --output=mdtest.1n.out
#SBATCH --error=mdtest.1n.err
#SBATCH --chdir=/home/hpc/benchmarks-2022/mdtest/ior-3.3.0/bin/1n

module purge
module load intel-mpi/2021.7.0 intel-compilers/2021.2.0
export I_MPI_PIN_PROCESSOR_LIST="all:map=scatter"

mpirun -n $SLURM_NPROCS /home/hpc/benchmarks-2022/mdtest/ior-3.3.0/bin/mdtest -d /scratch/hpc -i 10 \
                        -z 10 -I 100
 
