#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "trjconv"   # job name

module purge
module load gromacs/2021.4

echo 1 | gmx trjconv -f dynamic.xtc -s dynamic.tpr -o ca_clean.xtc -pbc nojump -skip 50
echo 1 | gmx trjconv -f minimization.gro -s dynamic.tpr -o ca_clean.pdb -pbc nojump