#!/bin/bash

reps=20

for ((i=0; i<reps; i++))
do
    echo "Making directory and copying files"
    mkdir "run_${i}"
    cp index.ndx "run_${i}"
    cp *.mdp "run_${i}"
    cp system.top "run_${i}"
    cp system.gro "run_${i}"
    cd "run_${i}"
    
    echo "Submitting run $i"
    sbatch ../submit.sh $i
    cd ..

done

# Next time, execute rep_run.sh >> num_track.txt
# Or figure out a macro for the #SBATCH lines? Looks like bash doesn't do that.
