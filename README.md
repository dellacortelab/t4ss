# T4SS simulations
Instructions for taking original PDB structure all the way to full course-grain simulations, using DotB and DotO complexes as examples.


# Splitting
First, the full PDB file must be split into pieces with at most 99,999 atoms, or else martinize will fail with atom indices.
(There may be another roadblock with Atom numbers, Bryce is/was investigating).
I have used both manual selection in ChimeraX for large sections, and simple bash commands for chain-by-chain splitting.
Use whatever method produces the desired sections.

For larger, multi-chain sections (which is probably most ideal), ChimeraX is easier.
First, pull up the whole structure, and then select only the chains which are part of the section (eg. select an atom, then use up arrow until entire chain is selected).
After that, simply hit file->save, making sure to check the "save selected atoms only" box.

IMPORTANT:
Depending on the initial structure used, chimera will keep the model numbers the same. This means if the atoms selected are not part of model 1 (like model 4, for example), then chimera will leave model 1 empty in the resulting pdb file.
THIS WILL BREAK MARTINIZE!!!
Martinize will fail if model 1 is empty and not the first model. To fix this, simply rename the correct model to model 1 in the pdb, and delete any extra model's and endmdl's.


# Fixing structure
For some chains, you will be ready to martinize after splitting. Others are missing atoms, and will need some cleaning to work with martinize.
Two methods are available: Manual rotamer replacement in chimera, or automatic side-chain fixing with scwrl.
If you try to martinize and it fails due to missing atoms, you will need to fix the structure.

In ChimeraX, select the problem residue(s), then go Tools->Structure Editing->Rotamers. Make sure the residue name matches, then hit apply.
After that, a box will appear with many possible rotamers. Select one (I just use the top pick), and hit "use chosen rotamer(s)", and you're done.

For scwrl, follow scwrl's included usage instructions.


# Martinize
Next, each piece must be martinized, and the resulting cg models and itps must be correctly named.
martinize.sh will do this for you. It is a messy script right now, but it should work.

For martinize.sh to work, you need to first have both martinize2 installed and dssp installed to a known path. Ideally these are installed in a new conda environment.

martinize2 will be installed with vermouth: 

    pip install vermouth
    
dssp will be installed through conda (MAKE SURE conda installs version >= 3.0.0): 
    
    conda install -c salilab dssp
    
On a mac, to get version 3.0+:

    conda install -c salilab dssp=3.0.0   

After dssp is installed, find the path to the dssp executable and replace the dssppath variable in martinize.sh.
You will also need the structure pdb in a directory named orig, and martinize.sh in a sibling directory to orig (i.e. accessible through ../orig).
Make sure your conda environment is active, and then simply run martinize.sh either as a sbatch (for large structures) or in command line.

IMPORTANT:
martinize2 looks for an executable named dssp, but conda installs an exe named mkdssp. You will need to either alias mkdssp to dssp or create a copy named dssp.

IMPORTANT:
The martinize.sh script renames the .itp files, but it doesn't change the [ moleculename ] inside the itp itself, which needs to be done separately.
I should probably just put it into the same script, but right now there is a rename.sh script to do that.


# Putting together Full Structure
Open Babel can be used to concat all the cg models together. Move all models into new directory, then run:

    obabel *_cg.pdb -O {full_name}.pdb

This will combine all the models together. To find the order they were combined in, simply grep all lines starting with COMPND:

    grep ^COMPND {full_name}.pdb > temp.txt
    cut -c 11- temp.txt > order.txt

This model order will be used in determining the order of the itps and molecule names in the .top file.


# Insane
To solvate the system and generate .gro and .top files, run insane. The command is quite simple.

    python2 insane.py -f {full_name}.pdb -o system.gro -p system.top -pbc square -d 1 -sol W -salt 0

-salt 0 neutralizes the complex with the minimum number of ions added. Change the value to desired molarity/molality (not sure which) if more ions are wanted.


# Prep .top file
Then, you need to fix the includes and molecule names in the .top file.

#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

gen_true.sh will create an .txt file with all the itp basenames in proper order.
You can copy and paste this list in both the include section, and the molecules.


# Run Simulations
Copy martini .itp files, .mdp files, and submit.sh into desired directory, and you should be good to run.
You will also need to create an index.ndx file. No new fields are required, just create the index with gmx make_ndx.

IMPORTANT:
There seems to be some disconnect between the ion names that Insane generates and the ones in the martini .itp files.
The workaround I'm using is modifying the names in the martini_v3.0.0_ions_v1.itp from NA and CL to NA+ and CL-
in both [moleculetype] and [atoms].
We can probably combine them into a single Ions group, but haven't tested that
