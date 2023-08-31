#!/bin/bash

cd first_failed
for run in run_*
do
	if [[ "$run" == "run_2" || "$run" == "run_9" || "$run" == "run_0" ]]; then
		continue
	fi
	cd $run
	sbatch ../../tenth.sh
	cd ..
done
cd ..

cd twenty_failed
for dir in run_*
do
        if [[ "$dir" == "run_3" || "$dir" == "run_5" || "$dir" == "run_8" || "$dir" == "run_16" ]]; then
		continue
	fi
        cd $dir
        sbatch ../../tenth.sh
        cd ..
done
