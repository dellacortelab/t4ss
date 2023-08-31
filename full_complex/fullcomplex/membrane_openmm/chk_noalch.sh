#!/bin/bash --login

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=2
#SBATCH --constraint=pascal # restricts marylou usage to the m9g nodes, which are compatible with cuda > 11.4, <=12
#SBATCH --mem-per-cpu=70G   # memory per CPU core, the first run used about 60 GB total
#SBATCH -J "chk_noalch series, redo"   # job name
#SBATCH --mail-user=sethang@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=2

mamba activate cudaeight

python ../chk_noalch.py
