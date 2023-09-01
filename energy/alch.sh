#!/bin/bash --login

#SBATCH --time=3:00:00   # walltime ## mixed precision consistently uses ~2h15m
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=1
#SBATCH --constraint=pascal # restricts marylou usage to the m9g nodes, which are compatible with cuda > 11.4, <=12
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "alch 43"   # job name

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

mamba activate cudaeight

python ../alch_43.py $1

# alch 61 has mixed precision, no mini.pdb, no defines block, 4GB, and xtc every 10k steps
# alch 25k has the above, but 10 GB and back down to 31 states (some redos went back down on memory)
# alch 43 is to try fewer coulombic lambda values and more vdw lambda values. Might as well give it a shot, I suppose.