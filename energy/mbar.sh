#!/bin/bash --login

#SBATCH --time=0:15:00   # walltime
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=1
#SBATCH --constraint=pascal # restricts marylou usage to the m9g nodes, which are compatible with cuda > 11.4, <=12
#SBATCH --mem-per-cpu=5G   # memory per CPU core
#SBATCH -J "25k nrg mbar, plug in solvent."   # job name

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

mamba activate cudaeight

python mbar.py