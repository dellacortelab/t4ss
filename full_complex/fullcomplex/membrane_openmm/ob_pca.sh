#!/bin/bash --login


#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=1
#SBATCH --mem-per-cpu=20G   # memory per CPU core
#SBATCH -J "omm ob pca attempts"   # job name


mamba activate tmpi_23

NAME=ob
ID=12

echo $NAME
echo $ID

# Create and enter directory
[ -d ${NAME} ] || mkdir ${NAME}
cd ${NAME}

# Analysis
echo ${ID} ${ID} | gmx covar -s ../../clean.gro -f ../clean.xtc -n ../../index.ndx

echo ${ID} ${ID} | gmx anaeig -s ../../clean.gro -f ../clean.xtc -n ../../index.ndx -v eigenvec.trr -eig eigenval.xvg -extr proj.pdb -first 1 -last 2 -nframes 100

cd ..

