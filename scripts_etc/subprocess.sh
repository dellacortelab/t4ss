#!/bin/bash 


#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH -J "gmx subanalysis"   # job name


module purge
module load cuda/10.1 gromacs

NAME=${1}
ID=${2}

# Create and enter directory
[ -d ${NAME} ] || mkdir ${NAME}
cd ${NAME}

# Analysis
echo ${ID} ${ID} | gmx covar -s ../clean.gro -f ../clean.xtc -n ../index.ndx

echo ${ID} ${ID} | gmx anaeig -s ../clean.gro -f ../clean.xtc -n ../index.ndx -v eigenvec.trr -eig eigenval.xvg -extr proj.pdb -first 1 -last 2 -nframes 100