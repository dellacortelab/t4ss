# Contents here on Github
- Here in the readme, a summary of our work
- Contained in the directory tree:
    - Working files
        - Pipeline scripts (.sh and .py)
        - Structure files, atomistic and coarse grain (.pdb and .gro)
        - Topology files (.top and included .itp)
        - Notes for clarification (.txt)

# T4SS Summary

## Background
 ~30 years ago, some Legionella (bacterial genus) genes were found, when deleted, to make them Defective in Organelle Trafficking and/or IntraCellular Multiplication (Dot/Icm).
 This set of genes became categorized as one of two subtypes of the Type IV Secretion System, but since we only work with this subtype, we've just been saying T4SS.
 Grant Jensen (at CalTech) and collaborators applied techniques of electron microscopy to get a general spatial rendering of Legionella Pneumophila T4SS in action.
 The T4SS complex is hard to characterize, so molecular dynamics modeling is one approach that helps to check the conclusions reached and possibilities raised in the work already done by Jensen and collaborators.

 Between February and August of 2023, several people have been involved on the MD side: Brenden Stark, Cayson Hamilton, Gus Hart, Ethan Smith, Dennis Della Corte, and Bryce Hedelius. Stefano Maggi (post-doc here at BYU with Grant Jensen) has been a major liason between the Jensen model and our work with MD.
- Dennis: largely a supervisory and consulting role
- Brenden: a ton of troubleshooting to develop a consistent method that yielded models that wouldn't explode in Gromacs, steered MD with the plug
- Cayson and Gus: adding lipids to the computer model
- Ethan: applying Brenden's method to the rest of the T4SS subregions, writing some scripts to automate the method, troubleshooting alchemical free energy with the plug, transferring the system to OpenMM
- Bryce: bringing mathematical and coding skills to the OpenMM alchemical free energy problem that Ethan didn't have

## Issues dealt with
- T4SS is too big for an atomistic simulation.
	- Martini 3 coarse-grain force field
- T4SS doesn't have an existing martini model
	- martinize2 renders atomistic into coarse-grain
- T4SS is too big for martinize2
	- T4Ss was subdivided into regions that fell within the size limits
- Some subdividing methods have naming side effects that interrupt martinize2's functionality
	- Renaming is pretty easy
- Sometimes atoms are missing
	- Manually fill in rotamers with ChimeraX
- The model doesn't have solvent
	- insane.py, provided by martini, puts those in and gives you the files gromacs needs
- Insane.py gives files gromacs needs, but the files aren't quite ready for gromacs
	- Editing molecule names isn't too hard.
	- Adding the necessary "#include" lines isn't too hard either. Just keep them in the right order.
- While we can see the cg models not exploding, they're still just pieces
	- Concatenating the *_cg.pdb files is pretty easy. Just keep track of the order, MODEL number lines, and newlines added in the concatenation process
- Sometimes things explode
	- Make sure your periodic box is big enough and you've tried multiple times. Failing on one random seed doesn't mean it will fail on every random seed.
- T4SS spans two membranes. The computer model doesn't have lipid membranes
	- Insane.py can add those too. It will need to be run twice to add two membranes. Only add solvent and ions on the second run.

## Issues still pending
- Martini is pretty much made for Gromacs (They came out of the same university in the Netherlands), but our martinized T4SS runs into some issues in Gromacs. The pieces don't all stick together in longer simulations
- T4SS gets stuff from one cell into another cell, but the path through the current T4SS model is blocked by a sort of plug (IcmX).
	- Jensen et al believe the plug must come out before the T4SS can do its job
	- MD is often used to calculate the free energy of interaction between two entities. This can inform questions of association/dissociation
- Free energy
	- Steered MD, creating an artificial force to pull the plug out: computationally expensive?
	- Alchemical
- Alchemical free energy
	- Gromacs is not equipped to handle the restricted angle bending potentials assigned by martinize2
	- Changing the given potential (code 10) to a harmonic angle potential (2) makes the system explode pretty quickly. Dr. Marrink says 1 could also work.
	- OpenMM?
- Putting everything in OpenMM requires some method adaptation
	- Upside: OpenMM seemed to remove the issue of the complex dissociating from itself, and the martini_openmm adapter script gets around the restrict bending problem
	- Downside: it's hard to determine if our free energy values are accurate
	- On the side: we need to optimize the lambda schedule. The overlap matrix seems pretty useful, but it doesn't seem to be a reliable predictor of the convergence of the energy values.
- Test cases
	- We need to know what we're representing.
	- Step one for our alchemical T4SS simulation: get solvation energy for IcmX.
	- What systems have known solvation energy?
	- Benzene. But the number usually cited for that is relative, between gas and liquid. Simulating in gas is a challenge. I did find solvation energies calculated from the partition between water and octanol, and simulating in octanol might be easier. I made a martini3 octanol box, but the simulation won't run for a reason I haven't yet understood.
	- Alanine peptides have been simulated with multiple computational methods in mutual agreement. We simulated a decaalanine and found discrepancy, but the literature values depended on some parameters we didn't use the first time through.



# Brenden's notes on the pipeline to getting a working Gromacs simulation

## T4SS simulations
Instructions for taking original PDB structure all the way to full coarse-grain simulations.
I've included an example icmf directory as an example file structure.

### Steps
1. Create model file, with < 99,999 Atoms
    1a. (if necessary) Fix structure/side chains/missing atoms
2. Martinize structure(s) (using martinize.sh)
    2a. (If using multiple structures) Combine structures, keeping track of molecule order
    2b. Use rename.sh to rename molecule names within .itp files
3. Solvate using insane.py on cg file
4. Fill out .top file with .itp includes and molecule names
5. Create index.ndx file with gmx make_ndx
6. Copy all necessary files (including .mdp files) to simulation directory, if desired
    6a. Double check all paths are correct, especially paths to .itps in .top file
7. Start simulation using submit.sh (modify dssp path within submit.sh first)
8. Pray everything works :)

### Splitting
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


### Fixing structure
For some chains, you will be ready to martinize after splitting. Others are missing atoms, and will need some cleaning to work with martinize.
Two methods are available: Manual rotamer replacement in chimera, or automatic side-chain fixing with scwrl.
If you try to martinize and it fails due to missing atoms, you will need to fix the structure.

In ChimeraX, select the problem residue(s), then go Tools->Structure Editing->Rotamers. Make sure the residue name matches, then hit apply.
After that, a box will appear with many possible rotamers. Select one (I just use the top pick), and hit "use chosen rotamer(s)", and you're done.

For scwrl, follow scwrl's included usage instructions.
Also, you will have to remove a few lines at the top of the scwrl'd file, that start with SCWRL

NOTE:
It seems that scwrl'ing dotb makes martinize2 fail on HIS/HSD.
It recognizes that the HIS *should* be a HSD, but then it doesn't know what to do with HSD.
This is an interesting bug we should look into.


### Martinize
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
The `martinize.sh` script renames the `.itp` files, but it doesn't change the [ moleculename ] inside the itp itself, which needs to be done separately.
I should probably just put it into the same script, but right now there is a rename.sh script to do that.

EDIT:
It also seems that martinize doesn't like when the chain identifier is two letters. I thought I did that before, but now it's not working. 
Use obabel to quickly rewrite chain identifiers.

IMPORTANT:
martinize includes the dssp sequence in the .itp files it creates.
This can be a problem for longer chains, as the sequence can be longer than 4096 characters.
This will break some parsers (eg. insane.py), so simply delete the line from the .itp.


### Putting together Full Structure
Open Babel or cat can be used to concat all the cg models together. Move all models into new directory, then run:

    obabel *_cg.pdb -O {full_name}.pdb

or

    cat * >> {full_name}.pdb

Be careful when using cat that there is a newline after every file, otherwise gromacs will say you are missing atoms in the structure.

This will combine all the models together. If you used obabel, to find the order they were combined in, simply grep all lines starting with COMPND:

    grep ^COMPND {full_name}.pdb > temp.txt
    cut -c 11- temp.txt > order.txt

This model order will be used in determining the order of the itps and molecule names in the .top file.



### Insane
To solvate the system and generate .gro and .top files, run insane. The command is quite simple.

    python2 insane.py -f {full_name}.pdb -o system.gro -p system.top -pbc square -d 1 -sol W -salt 0

With custom box:

    python2 insane.py -f {full_name}.pdb -o system.gro -p system.top -x {x} -y {y} -z {z} -sol W -salt 0

-salt 0 neutralizes the complex with the minimum number of ions added. Change the value to desired molarity/molality (not sure which) if more ions are wanted.

NOTE:
It seems for larger systems, insane doesn't properly create the pbc, i.e. the box is too small, leading to problems.
My workaround is either to use gmx editconf, or define my own box vectors and manually put those in.
This needs more experimentation to determine what is happening.


### Prep .top file
Then, you need to fix the includes and molecule names in the .top file.

#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

DEPRECATED: gen_true.sh will create an .txt file with all the itp basenames in proper order. 
You can copy and paste this list in both the include section, and the molecules.

You can just copy the includes and molecules from the .top file created from martinizing.
Make sure the names are accurate.

NOTE:
Martinize allows identical chains to share .itp files. This can mess up the order of molecules.
Check the .top file generated by martinize2 to make sure order is ok.
The .top file will also include duplicate .itp includes, which will cause gromacs to freak out.
Remove the duplicate includes, but not the duplicate molecule types.


### Run Simulations
Copy martini `.itp` files, `.mdp` files, and submit.sh into desired directory, and you should be good to run.
You will also need to create an index.ndx file. No new fields are required, just create the index with gmx make_ndx.

IMPORTANT:
There seems to be some disconnect between the ion names that Insane generates and the ones in the martini .itp files.
The workaround I'm using is modifying the names in the martini_v3.0.0_ions_v1.itp from NA and CL to NA+ and CL-
in both [moleculetype] and [atoms].
We can probably combine them into a single Ions group, but haven't tested that


### Combining models
To create larger systems from smaller, already tested systems, simply combine the CG models using cat or obabel.
Make sure to keep track of which order you combine the models in, so that you can input the molecules in the .top file in the correct order.
(obabel will automatically label which parts came from which file, so that can be helpful for ordering things).
.itps should be able to be reused from smaller systems to combined systems, just copy them over.

After the CG models are combined, continue with the same procedure, starting with the insane.py step

NOTE:
Some models may introduce clashes by combining them together. This is where user judgement comes into play.
Either delete or move problematic chains, or otherwise move models until the run successfully makes it past the equilibration phase.

NOTE:
Sometimes when combining CG models, the system.gro and system.top will have different #'s of atoms.
This isn't very consistent, and I'm not sure where it's coming from.
My best guess would be copying the wrong cg and/or itp files.
If this happens, my best suggestion is just restarting by copying the correct files again and running from insane.
