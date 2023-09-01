#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=24  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gres=gpu:4
#SBATCH --mem-per-cpu=5G   # memory per CPU core
#SBATCH -J "gmx run"   # job name
#SBATCH --mail-user=starkbrenden1@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.

# Needs to be (# of cores / # of GPUs)
export OMP_NUM_THREADS=6

export GMX_MAXCONSTRWARN=10000

module purge
module load cuda/11.4 gromacs/2021.4

echo "Start Minimization"
gmx grompp -f minimization.mdp -c system.gro -p system.top -o minimization.tpr 
gmx mdrun -deffnm minimization -v -ntmpi 1

echo "Start Equilibration"
gmx grompp -f equilibration.mdp -c minimization.gro -p system.top -o equilibration.tpr -r system.gro -n index.ndx
gmx mdrun -deffnm equilibration -v -mn index.ndx -ntmpi 4

echo "Start Production"
gmx grompp -f dynamic.mdp -c equilibration.gro -p system.top -o dynamic.tpr -r system.gro -n index.ndx
gmx mdrun -deffnm dynamic -v -mn index.ndx -ntmpi 4