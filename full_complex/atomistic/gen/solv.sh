#!/bin/bash

#SBATCH --time=2:00:00   # walltime
#SBATCH --ntasks=12  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1G   # memory per CPU core

module purge
module load cuda/10.1 gromacs/2019.4

gmx solvate -cp box.gro -o solv.gro -p topol.top