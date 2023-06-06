import sys
import string

"""
Replaces 2-letter chain ids with one-character sequential chain ids, outputting to a new file called out.pdb.
Also creates a file called chains.txt listing which ids were replaced with which identifier.
First (and only) argument is the file to replace ids.
"""

# Get file name from first argument
file = sys.argv[1]

# Create list of valid 1-character chain ids
newids = string.printable

# Remove bad ascii characters
newids = newids.replace('/','')
newids = newids.replace('!','')
newids = newids.replace('#','')

# Read original pdb
with open(file, 'r') as f:
    lines = f.readlines()

# Start indexing file
with open('chains.txt', 'w') as f:
    f.writelines("OldChain,NewChain" + '\n')

# Loop through lines and replace ids
count = -1
oldid = ""
newlines = []
for line in lines:
    if not line.startswith('ATOM') and not line.startswith('TER'):
        newlines.append(line)
        continue

    if line[20:22] != oldid:
        count += 1
        oldid = line[20:22]
        with open('chains.txt', 'a') as f:
            f.writelines(oldid + "," + newids[count] + '\n')

    newlines.append(line[0:20] + " " + newids[count] + line[22:])

with open('out.pdb', 'w') as f:
    f.writelines(newlines)
    
