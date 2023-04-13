#!/bin/bash

for itp in ./*.itp; do
    name=${itp%.*}
    name=${name##*/}
    echo "renaming $name"
    sed -i "s/molecule_[0-9]/$name/g" $itp
done