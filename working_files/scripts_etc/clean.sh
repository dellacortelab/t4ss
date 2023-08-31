#!/bin/bash 


#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10M   # memory per CPU core
#SBATCH -J "gmx analysis"   # job name

module purge
module load cuda/11.4 gromacs/2021.4

echo '1 1' | gmx trjconv -f ../minimization.gro -s ../dynamic.tpr -center -pbc nojump -o clean.gro -skip 100
echo '1 1' | gmx trjconv -f ../dynamic.xtc -s ../dynamic.tpr -center -pbc nojump -o clean.xtc -skip 100

sbatch subprocess.sh dot_b 18
sbatch subprocess.sh dot_o 19
sbatch subprocess.sh o+b 20