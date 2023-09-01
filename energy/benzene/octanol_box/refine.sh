#!/bin/bash --login

#SBATCH --time=1:00:00   # walltime
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=1
#SBATCH --mem-per-cpu=1G   # memory per CPU core
#SBATCH -J "Octanol minimization"   # job name
#SBATCH --mail-user=sethang@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=1

mamba activate tmpi_23


echo "Start Minimization"
gmx grompp -f minimization.mdp -c moonshine.gro -p moonshine.top -o minimization.tpr 
gmx mdrun -deffnm minimization -v -ntmpi 1

# First two tries: wrong atom line format
# Third try: needed ntmpi