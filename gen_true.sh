#!/bin/bash

rm true_order.txt

while read name; do
    name=${name%_cg.pdb}
    for itp in itps/$name*; do
        echo "#include \"../$itp\"" >> true_order.txt
    done
done < order.txt

while read name; do
    name=${name%_cg.pdb}
    for itp in itps/$name*; do
        itpname=${itp%.*}
        itpname=${itpname##*/}
        echo "$itpname          1" >> true_order.txt
    done
done < order.txt