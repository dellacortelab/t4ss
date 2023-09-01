#!/bin/bash

# Developed by Brenden Stark. Where the previously automated pipeline produced .itp files with generic names in them, 
# I believe this script replaces the generic names in the file with the file name itself 

for itp in ./*.itp; do
    # Rename molecule name
    name=${itp%.*}
    name=${name##*/}
    echo "renaming $name"
    sed -i "s/molecule_[0-9]\+/$name/g" $itp

    # Remove long ssd lines
    # Should remove lines longer than 2048 characters, could change to 4096
    sed '/.\{2048\}/d' $itp > ${itp}_temp
    mv ${itp}_temp $itp
done