#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=18  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=3
#SBATCH --mem-per-cpu=1G   # memory per CPU core
#SBATCH -J "MD reps"   # job name
#SBATCH --mail-user=sethang@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE


module purge
module load cuda/10.1 gromacs/2019.4

echo "Start Minimization"
gmx grompp -f minimization.mdp -c system.gro -p system.top -o minimization.tpr 
gmx mdrun -deffnm minimization -v -ntmpi 3

echo "Start Equilibration"
gmx grompp -f equilibration.mdp -c minimization.gro -p system.top -o equilibration.tpr -r system.gro -n index.ndx
gmx mdrun -deffnm equilibration -v -mn index.ndx -ntmpi 3

echo "Start Production"
gmx grompp -f dynamic.mdp -c equilibration.gro -p system.top -o dynamic.tpr -r system.gro -n index.ndx
gmx mdrun -deffnm dynamic -v -mn index.ndx -ntmpi 3

