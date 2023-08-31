#!/bin/bash

# MODIFIED TO EXCLUDE MARTINIZE
###
# This version, 5/26/23, automates a few more steps in the pipeline. It takes care of:
# - rename.sh
# - insane.py
# - copy-ready .top file
# - make ndx
# - itp includes
# - .top molecule names


# What it doesn't do yet:
# - activate martinize environment
# - start simulation


#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=5120M   # memory per CPU core
#SBATCH -J "lower region martinize full"   # job name


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# ACTIVATE MARTINIZE CONDA ENV BEFORE RUNNING THIS SCRIPT

# Change dssp path
dssppath='/home/sethang/.conda/envs/second_env/bin/dssp'

echo "Start Martinize"

for chain in ../orig/*; do
    echo "Martinizing $chain"
    cp $chain .
    name=${chain%.*}
    name=${name##*/}
    echo $name
    
    
    echo "Calling insane.py"
    python2 ~/fsl_groups/fslg_dellacortelab/compute/t4ss/gitstuff/t4ss/scripts_etc/insane.py -f ${name}.pdb -o system.gro -p system.top -pbc square -d 1 -sol W -salt 0


    echo "system.top edit"
    # Get rid of martini.itp include and protein line
    sed -i '1d' system.top
    place=8
    sed -i "${place}d" system.top


    count=1
    for itp in martini*; do
        sed -i "${count}s|^|#include \"../itps/${itp}\"\n|" system.top
        ((count++))
    done
    mkdir martemp
    mv martini* martemp
    for itp in *.itp; do
        sed -i "${count}s|^|#include \"../itps/${itp}\"\n|" system.top
        ((count++))
    done
    unset count
    mv martemp/* .
    rmdir martemp

    echo "Making ndx"
    module purge
    module load gromacs/2021.4
    gmx make_ndx -f system.gro
done
    

echo "More reorganizing"
mkdir ../gmx
mv system.* ../gmx
mv index.ndx ../gmx
cp  ~/fsl_groups/fslg_dellacortelab/compute/t4ss/gitstuff/t4ss/scripts_etc/*.mdp ../gmx

echo ""
echo "Todo:"
echo "Pull in the submit.sh that you want"
echo "Submit!"
