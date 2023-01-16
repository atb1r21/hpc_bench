#!/bin/bash
#SBATCH --nodes 4
#SBATCH -p batch
#SBATCH --tasks-per-node=40
#SBATCH --ntasks=160
#SBATCH --res=benchmarks

module load gromacs/2022.3/intel-2022-libs

mpirun -np 160 gmx_mpi mdrun -ntomp 1 -deffnm benchPEP-h
