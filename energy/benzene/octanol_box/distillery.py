'''
Okey here's what we gotta do.
Col 1       Col 2   Col 3   Cols 4-6
Residue     Type    Index   XYZ
(per oco)   (of 3)  (atom)  Copy/paste, increment
--
1. Figure out the coordinates for the first one from the solvents itp
2. Check that 1183 fit in the same box size
    - Or maybe just 1000 for easy numbers?
3. Iterate in three dimensions
--
More info for the atom line
5 - 5 - 5 - 5 - 8(3)*3 - [velocities...?]

'''
import math

num_grains = 3
space_per_particle = 0.8
num_1D = 6
num_x = num_1D*2
cube_dimension = 6 * space_per_particle

def place(particle,count,i,j,k):
    if particle == 'Octanol':
        count = place('C1',count,i,j,k)
        count = place('C2',count,i,j,k)
        count = place('PC',count,i,j,k)
    else:
        resnum = int((count - 1)/ num_grains) + 1
        x,y,z = get_xyz(particle,i,j,k)
        print(format_string(resnum,particle,count,x,y,z))
        return count+1
    return count


bond_1_2 = 0.39
bond_2_3 = 0.35
bond_angle = math.radians(125)
bond_2_3_y = math.sin(bond_angle)*bond_2_3
bond_2_3_z = math.cos(bond_angle)*bond_2_3
bubble = space_per_particle
x_bubble = bubble/2
offset = 0.1
def get_xyz(particle,i,j,k):
    x=i*x_bubble + offset
    y=j*bubble + offset
    z=k*bubble + offset

    if particle == 'C1':
        pass
    elif particle == 'C2':
        z += bond_1_2
    else:
        y += bond_2_3_y
        z += bond_1_2
        z += bond_2_3_z

    x = "{:.3f}".format(x)
    y = "{:.3f}".format(y)
    z = "{:.3f}".format(z)

    return x,y,z

def format_string(resnum,particle,count,x,y,z):
    res = gro_spacer(resnum)
    oco = 'OCO  '
    par = gro_spacer(particle)
    cou = gro_spacer(count)
    coord = coord_spacer(x,y,z)

    return res + oco + par + cou + coord

def gro_spacer(input,space=5):
    i_string = str(input)
    diff = space - len(i_string)
    rest = ' ' * diff
    output = rest + i_string
    return output

def coord_spacer(x,y,z,space=8):
    x_string = gro_spacer(x,space)
    y_string = gro_spacer(y,space)
    z_string = gro_spacer(z,space)
    return x_string + y_string + z_string
    

count=1
for i in range(num_x):
    for j in range(num_1D):
        for k in range(num_1D):
            count = place('Octanol',count,i,j,k)