#!/bin/bash

for i in {0..30}
do
    echo "Making directory"
    mkdir "lambda_${i}"
    cd "lambda_${i}"
    
    echo "Submitting lambda $i"
    sbatch ../alch.sh $i
    cd ..

done