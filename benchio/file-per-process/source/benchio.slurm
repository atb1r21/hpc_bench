#!/bin/sh

#SBATCH --ntasks-per-node=16     
#SBATCH --nodes=32
#SBATCH --exclusive
#SBATCH --partition=batch
#SBATCH --time=00:50:00
#SBATCH --ntasks=512
#SBATCH --output=benchio_fpp.512.out
#SBATCH --error=benchio_fpp.512.out

# vary the value on the right-hand side of --nodes= and --ntasks= (which equals to --nodes times 16)

module purge
module load openmpi

mpirun -n $SLURM_NPROCS benchio_fpp.x
