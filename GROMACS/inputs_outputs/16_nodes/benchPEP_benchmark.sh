#!/bin/bash
#SBATCH --nodes 16
#SBATCH -p batch
#SBATCH --tasks-per-node=40
#SBATCH --ntasks=640

module load gromacs/2022.3/intel-2022-libs

mpirun -np 640 gmx_mpi mdrun -ntomp 1 -deffnm benchPEP-h
