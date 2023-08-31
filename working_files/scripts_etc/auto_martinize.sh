#!/bin/bash

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
    martinize2 -f ${name}.pdb -o ${name}.top -x ${name}_cg.pdb -dssp $dssppath -p backbone -ff martini3001 -elastic -ef 700.0 -el 0.5 -eu 0.9 -ea 0 -ep 0 -scfix -cys auto -maxwarn 30
    count=$(ls chain_*.ssd | wc -l)
    if [ $count -gt 1 ]; then
        for ssd in chain_*.ssd; do
            num=${ssd%.*}
            num=${num##*_}
            mv $ssd ${name}_${num}.ssd
        done
    else
        mv chain_*.ssd ${name}.ssd
    fi
    count=$(ls molecule_*.itp | wc -l)
    if [ $count -gt 1 ]; then
        for itp in molecule_*.itp; do
            num=${itp%.*}
            num=${num##*_}
            mv $itp ${name}_${num}.itp
        done
    else
        mv molecule_*.itp ${name}.itp
    fi
    rm ${name}.pdb



    # NOTE: the following works only if there is only one .pdb file
    
    echo "Starting .top fix"
    # Perform the in-place replacement using sed
    sed -i "s/molecule/${name}/g" ${name}.top
    sed -i "s/\"${name}/\"..\/itps\/${name}/g" ${name}.top
    
    echo "Calling insane.py"
    python2 ~/fsl_groups/fslg_dellacortelab/compute/t4ss/gitstuff/t4ss/scripts_etc/insane.py -f ${name}_cg.pdb -o system.gro -p system.top -pbc square -d 1 -sol W -salt 0

    echo "Copying molecules"
    echo ${name}
    # Get rid of martini.itp include and protein line
    sed -i '1d' system.top
    place=8
    sed -i "${place}d" system.top 
    sed -n "/\[ ${name}s \]/,\${/\[ ${name}s \]/d;p}" ${name}.top > molecules.txt
    while IFS= read -r line
    do
        sed -i "${place}s|^|${line}\n|" system.top
        ((place++))
    done < molecules.txt
    rm molecules.txt


    count=1
    for itp in *.itp; do
        sed -i "${count}s|^|#include \"../itps/${itp}\"\n|" system.top
        ((count++))
    done
    unset count

    echo "Making ndx"
    module purge
    module load gromacs/2021.4
    gmx make_ndx -f system.gro
done
    
echo "Calling rename.sh"
~/fsl_groups/fslg_dellacortelab/compute/t4ss/gitstuff/t4ss/scripts_etc/rename.sh

echo "Reorganizing"
mkdir ../itps
mv *.itp ../itps
cp  ~/fsl_groups/fslg_dellacortelab/compute/t4ss/gitstuff/t4ss/scripts_etc/*.itp ../itps
rm ../itps/martini_v3.0.0_phospholipids_v1.itp
mkdir ssd
mv *.ssd ssd

path="../itps/"
echo "Including"
include() {
    count=1
    for itp in ${path}${1}*; do
        sed -i "${count}s|^|#include \"${itp}\"\n|" system.top
        ((count++))
    done
}
include "martini"

echo "More reorganizing"
mkdir ../gmx
mv system.* ../gmx
mv index.ndx ../gmx
cp  ~/fsl_groups/fslg_dellacortelab/compute/t4ss/gitstuff/t4ss/scripts_etc/*.mdp ../gmx

echo ""
echo "Todo:"
echo "Pull in the submit.sh that you want"
echo "Submit!"
