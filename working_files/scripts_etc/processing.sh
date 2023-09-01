#!/bin/bash 

###
# PCA preprocessing and PCA kickoff
###


#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "gmx analysis"   # job name


module purge
module load cuda/10.1 gromacs


# Clean trajectory and make index file
gmx editconf -f ../equilibration.gro -o clean.gro
gmx make_ndx -f clean.gro -o index.ndx<<EOF
a 1 - 4980 & a BB
a 4981 - 33008 & a BB
a BB
name 14 DotB-BB
name 15 DotO-BB
name 16 BB
q
EOF
echo '1 1' | gmx trjconv -f ../dynamic.xtc -s clean.gro -center -pbc nojump -o clean.xtc

sbatch subprocess.sh dot_b 14
sbatch subprocess.sh dot_o 15
sbatch subprocess.sh bb 16

