#!/bin/bash

# Loop through all directories beginning with "run_" (except run_6)
for dir in run_*; do
  # Check if the directory is "run_6", skip it if so
  if [[ "$dir" == "run_3" || "$dir" == "run_5" || "$dir" == "run_8" || "$dir" == "run_16" ]]; then
    continue
  fi

  # Change to the directory and execute "trim.sh"
  cp trim.sh $dir
  cd $dir
  echo "Running trim.sh in $dir..."
  sbatch trim.sh
  cd ..
done
