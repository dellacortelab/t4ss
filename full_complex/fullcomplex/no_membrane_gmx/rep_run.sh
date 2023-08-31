#!/bin/bash

reps=${1}

for ((i=10; i<reps; i++))
do
    echo "Making directory and copying files"
    mkdir "run_${i}"
    cp index.ndx "run_${i}"
    cp *.mdp "run_${i}"
    cp system.top "run_${i}"
    cp system.gro "run_${i}"
    cd "run_${i}"
    
    echo "Submitting run $i"
    sbatch ../submit.sh
    cd ..

done
