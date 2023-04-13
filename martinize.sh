#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1G   # memory per CPU core



# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# ACTIVATE MARTINIZE CONDA ENV BEFORE RUNNING THIS SCRIPT

# Change dssp path
dssppath='~/miniconda3/envs/mart/bin/mkdssp'

echo "Start Martinize"

for chain in ../orig/*; do
    echo "Martinizing $chain"
    cp $chain .
    name=${chain%.*}
    name=${name##*/}
    echo $name
    martinize2 -f ${name}.pdb -o ${name}.top -x ${name}_cg.pdb -dssp $dssppath -p backbone -ff martini3001 -elastic -ef 700.0 -el 0.5 -eu 0.9 -ea 0 -ep 0 -scfix -cys auto -maxwarn 30
    count=$(ls chain_*.ssd | wc -l)
    if [ $count -gt 1 ]; then
        i=1
        for ssd in chain_*.ssd; do
            mv $ssd ${name}_${i}.ssd
            ((i++))
        done
    else
        mv chain_*.ssd ${name}.ssd
    fi
    count=$(ls molecule_*.itp | wc -l)
    if [ $count -gt 1 ]; then
        i=1
        for itp in molecule_*.itp; do
            mv $itp ${name}_${i}.itp
            ((i++))
        done
    else
        mv molecule_*.itp ${name}.itp
    fi
    rm ${name}.pdb
done

