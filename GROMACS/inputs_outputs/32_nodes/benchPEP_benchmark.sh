#!/bin/bash
#SBATCH --nodes 32
#SBATCH -p batch
#SBATCH --tasks-per-node=40
#SBATCH --ntasks=1280

module load gromacs/2022.3/intel-2022-libs

mpirun -np 1280 gmx_mpi mdrun -ntomp 1 -deffnm benchPEP-h
