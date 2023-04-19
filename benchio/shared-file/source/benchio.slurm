#!/bin/sh

#SBATCH --ntasks-per-node=16     
#SBATCH --nodes=32
#SBATCH --exclusive
#SBATCH --time=00:50:00
#SBATCH --ntasks=512
#SBATCH --output=benchio.512.out
#SBATCH --error=benchio.512.out

# modify the right-hand side of --nodes= and --ntasks= (which should equal --nodes times 16)

module purge
module load openmpi

mpirun -n $SLURM_NPROCS benchio.x
