#!/bin/bash

# Function to delete dynamic.xtc if both files exist
delete_dynamic_xtc() {
  if [ -e "dynamic.xtc" ] && [ -e "dynamic_tenth.xtc" ]; then
    rm dynamic.xtc
    echo "Deleted dynamic.xtc in $PWD"
  fi
}

# Loop through each directory
for dir in */; do
  if [ -d "$dir" ]; then
    # Enter the directory
    cd "$dir" || continue

    # Call the delete_dynamic_xtc function for the current directory
    delete_dynamic_xtc

    # Go back to the parent directory
    cd ..
  fi
done
