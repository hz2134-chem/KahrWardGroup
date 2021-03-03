#!/mnt/c/Users/zhouh/python3

#############################################################
# ./createCrystal.py *.pdb lSize hSize --coreSize=coresize --bVec=burgers_vector_size --hand=handness
############################################################

import argparse
import copy
from math import sin, cos, radians

def replicate_atoms_by_vector(atomList, num_atom_unitCell, num_molecule_unitCell, vector):
    for atom in atomList:
        atom[1] += num_atom_unitCell
        atom[4] += num_molecule_unitCell
        atom[5] += vector[0]
        atom[6] += vector[1]
        atom[7] += vector[2]
    return atomList

def printCell(atomList):
    for atom in atomList:
        oFile.write('{:<4}{:>7}{:>5}{:>4}  {:<4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6}{:>6}          {:>2}{:>14}\n'.format(*atom))

def circle(lastCell, num_atom_unitCell, num_molecule_unitCell, n, layerZ, handness=False):
    if handness:
        dis_core_shift = 0.5
        hd = -1 if handness == 'l' else 1
        if n <= coreSize:
            z = hd * bVector * c * (n / coreSize) / (6 * n)
        else:
            z = hd * bVector * c / (6 * n)
    else:
        z = 0
    turnList = [[a * cos(radians(120)), b * sin(radians(120)), z], 
                [a * cos(radians(180)), b * sin(radians(180)), z],
                [a * cos(radians(240)), b * sin(radians(240)), z],
                [a * cos(radians(300)), b * sin(radians(300)), z],
                [a * cos(radians(0)), b * sin(radians(0)), z],
                [a * cos(radians(60)), b * sin(radians(60)), z]]


    prev_unitCell = copy.deepcopy(lastCell)
    if handness:
        if n <= coreSize:
            prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, [a, 0, -hd * (2*n - 1)/(2*coreSize) * bVector * c]) 
            if n < coreSize:
                prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, [0, dis_core_shift, 0]) 
        else:
            prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, [a, 0, -hd * bVector * c])
    else:
        prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, [a, 0, 0])
    printCell(prev_unitCell)
    if handness and n < coreSize:
        prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, [0, -dis_core_shift, 0]) 
    
    for t in range(len(turnList)):
        for i in range(1, n + 1):
            prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, turnList[t]) 
            if t == len(turnList) - 1:
                if not handness:
                    if i == n:
                        continue
                else:
                    if n < coreSize:
                        if i < n:
                            prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, [0, -dis_core_shift/(n - 1),0]) 
                        else:
                            if n > 1:
                                prev_unitCell = replicate_atoms_by_vector(prev_unitCell, num_atom_unitCell, num_molecule_unitCell, [0, dis_core_shift,0]) 
                            continue
                    else:
                        if i == n:
                            continue
            printCell(prev_unitCell)
    return prev_unitCell
           
def layer(centerCell, l, handness, layerZ):
    # print the unit cell in the center first
    printCell(centerCell)

    lastCell = copy.deepcopy(centerCell)
    # replicate unit cell in counter clock circle
    for i in range(1, l):
        lastCell = circle(lastCell, num_atom_unitCell, num_molecule_unitCell, i, layerZ, handness)
    centerCell = replicate_atoms_by_vector(centerCell, num_atom_unitCell, num_molecule_unitCell, [0, 0, c])
    return centerCell 

if __name__ == "__main__":

    parser = argparse.ArgumentParser('Replicate a unit cell into a hexagonal crystal, with or without dislocation.\n./createCrystal.py unitCell.pdb 9 15 --coreSize=3 --bVec=1 --hand=r')
    parser.add_argument("file", help='The pdb file name for unit cell')
    parser.add_argument("replL", type=int, help='The repetition times in a (or b)')
    parser.add_argument("replH", type=int, help='The repetition times in c')
    parser.add_argument("--coreSize", type=int, help='The size of dislocation core')
    parser.add_argument("--bVec", type=int, help='The burgurs vector length. (in c direction)')
    parser.add_argument("--hand", choices=['l','r'], help='The handiness of the dislocation')
    args = parser.parse_args()

    iFile = open(args.file, 'r')
    fileName = args.file.split('.')[0]
    l, h = args.replL, args.replH
    coreSize, bVector = args.coreSize, args.bVec

    if coreSize:
        oFile = open(f'{fileName}_{l}{l}{h}_{coreSize}_{bVector}_{args.hand}.pdb', 'w')
    else:
        oFile = open(f'{fileName}_{l}{l}{h}.pdb', 'w')
    oFile.write(iFile.readline()) # write the first line to pdb file directly 
    
    # read unit cell parameters from input file
    tag, a, b, c, alpha, beta, gamma = iFile.readline().split()
    a, b, c = float(a), float(b), float(c)
    
    origin = [0, 0, 0]
    margin = 15
    
    # write the second line of pdb file
    oFile.write("{}{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f}\n".format(tag, 2*l*a + margin*2, 2*l*b + margin*2, h*c + margin*2, 90, 90, 90))
    
    # read atom list from input pdb file
    atomList = []
    for line in iFile.readlines():
        atomList.append(line.split())
    atomList = atomList[:-1]
    for i in range(len(atomList)):
        atomList[i][1] = int(atomList[i][1])
        atomList[i][4] = int(atomList[i][4])
        atomList[i][5] = float(atomList[i][5]) + l * a + margin
        atomList[i][6] = float(atomList[i][6]) + l * a + margin
        atomList[i][7] = float(atomList[i][7]) + margin
    num_atom_unitCell, num_molecule_unitCell = atomList[-1][1], atomList[-1][4]
    
    for z in range(h):
        atomList = layer(atomList, l, args.hand, z)
    oFile.write("END")

    iFile.close()
    oFile.close()
