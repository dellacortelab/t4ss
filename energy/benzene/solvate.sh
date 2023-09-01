#!/bin/bash

# The octanol box:
# http://www.cgmartini.nl/index.php/example-applications2/solvent-systems

gmx solvate -cp benzene_cg.pdb -cs octanol.gro -box 4.8 4.8 4.8 -o benz_octanol.gro