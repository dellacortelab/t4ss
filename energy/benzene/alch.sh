#!/bin/bash --login

#SBATCH --time=1:00:00   # walltime
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=1
#SBATCH --constraint=pascal # restricts marylou usage to the m9g nodes, which are compatible with cuda > 11.4, <=12
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "alch benzene"   # job name

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

mamba activate cudaeight

python ../alch.py $1
# First tries: using martini 2 octanol
# Round two, with Martini 3 octanol:
## First: used number of beads instead of number of molecules in the .top file
## Second: "openmm.OpenMMException: All Forces must have identical exclusions"
