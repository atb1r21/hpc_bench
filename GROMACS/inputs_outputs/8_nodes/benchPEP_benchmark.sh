#!/bin/bash
#SBATCH --nodes 8
#SBATCH -p batch
#SBATCH --tasks-per-node=40
#SBATCH --ntasks=320

module load gromacs/2022.3/intel-2022-libs

mpirun -np 320 gmx_mpi mdrun -ntomp 1 -deffnm benchPEP-h
