#!/bin/bash --login


#SBATCH --time=0:15:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1G   # memory per CPU core
#SBATCH -J "clean for pca"   # job name


mamba activate tmpi_23


# Clean trajectory and make index file
# echo '1 1' | gmx trjconv -f minimization.gro -s dynamic.tpr -center -pbc nojump -o clean.gro
# I used the clean.gro from the gromacs runs.
echo '1 1' | gmx trjconv -f backup2_prod.xtc -s ../system.gro -center -pbc nojump -o clean.xtc -skip 100

gmx make_ndx -f clean.gro -o index.ndx<<EOF
a 1 - 4976 & a BB
a 20791 - 47718 & a BB
name 10 B-BB
name 11 O-BB
10 | 11
name 12 OB
q
EOF

# sbatch ob.sh

