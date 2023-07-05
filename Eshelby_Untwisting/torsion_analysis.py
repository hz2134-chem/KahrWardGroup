#!/usr/bin/env python3
# Need to obtain last frame of the MD simulation before using this script:
# Example:
# gpta.x -i unitCell_6615.pdb trajectory.2.dcd -top -last -o 6615_ppp_nve_com_last.pdb
#
# Usage:
# ./torsion_analysis.py last_frame.pdb side_length z_length
# ./torsion_analysis.py 6615_last_frame.pdb 6 15

import numpy as np
import argparse
from itertools import zip_longest

def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def angle_vectors(v1,v2):
    v1_u=unit_vector(v1)
    v2_u=unit_vector(v2)
    dot_prod = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return np.round(np.rad2deg(np.arccos(dot_prod)), 4)


def torsion_angle(p1, a1, a2, p2, sign=False):
    vec_1 = p1 - a1
    vec_ax1 = a2 - a1
    vec_2 = p2 - a2
    vec_ax2 = -vec_ax1

    norm_1 = np.cross(vec_1, vec_ax1)
    norm_2 = np.cross(vec_ax2, vec_2)

    if sign:
        if angle_vectors(vec_2, norm_1) > 90:
            return -angle_vectors(norm_2, norm_1)

    return angle_vectors(norm_2, norm_1)

def calc_omega_phi_psi(coords):
    '''
    Calculate 1. two Ph - CO angle: omega and phi,
              2. torsion angle between two CO bond: psi
    @param: coordinates of all atoms in Benzil
    #1 C-O: axis: C7, side: O13
    #2 C-O: axis: C14, side: O16
    #1 Ph:  axis: C1, side: C2
    #2 Ph:  axis: C15, side: C17
    '''
    omega = torsion_angle(coords[7], coords[0], coords[1], coords[2], sign=True)  # O8-C1-C2-C3
    phi = torsion_angle(coords[15], coords[13], coords[14], coords[16], sign=True) # O16-C14-C15-C17
    psi = torsion_angle(coords[15], coords[13], coords[0], coords[7])  # O16-C14-C1-O8

    return omega, phi, psi

def get_benzil_coords(f):
    coords = []
    mol_id = None
    coords_with_element = []
    for _ in range(26):
        line = f.readline()
        line = line.split()
        if mol_id == None:
            mol_id = line[4]
        else:
            assert mol_id == line[4]
        coord_str = np.array(line[5:8])
        coords.append(coord_str.astype(np.float))
        coords_with_element.append([line[10]] + line[5:8])

    return mol_id, coords, coords_with_element

def get_omega_phi_psi(f):
    benzil_coords = get_benzil_coords(f)

    return calc_omega_phi_psi(benzil_coords)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage='Get statistics on omega, phi, psi of different layers from surface\n./torsion_analysis.py unitCell.pdb 6 15')
    parser.add_argument("file", help='The pdb file that contains only last frame')
    parser.add_argument("replL", type=int, help='The repetition times in a (or b)')
    parser.add_argument("replH", type=int, help='The repetition times in c')
    parser.add_argument("--psi_threshold", type=float, help='The threshold of psi below which the molecule will be selected')

    ## variable setup
    args = parser.parse_args()
    filename = args.file
    l = args.replL
    h = args.replH
    psi_threshold = args.psi_threshold
    output_name = filename if '/' not in filename else filename.split('/')[-1]
    output_name = output_name.rsplit('.', 1)[0] + '.csv'
    skip_circle = l - 6
    skip_uc_in_layer = 1 + 6 * (skip_circle * (skip_circle + 1) // 2)
    mol_in_uc = 3

    ## begin handling files
    f = open(filename, 'r')
    out_f = open(output_name, 'w')
    out_f.write('z_layer,circle_from_surface,mol_id,omega,phi,psi\n')

    f.readline()
    f.readline()
    # Each list contains 3 sublists, which represents layer 1,2,3 from surface
    omega_list = [[] for _ in range(5)]
    phi_list   = [[] for _ in range(5)]
    psi_list   = [[] for _ in range(5)]
    extract_mol_ids = []
    coords_with_element = dict()

    # Loop through all layers
    for z_layer in range(h):
        # Skip molecules in the core
        for _ in range(26 * skip_uc_in_layer * mol_in_uc):
            f.readline()
        # Loop from inner circle to surface circle 
        # (essentially layer 3 -> 2 -> 1 from surface)
        for i in range(1,6):
            for _ in range(6 * (skip_circle + i) * mol_in_uc):
                mol_id, mol_coords, ele_coords = get_benzil_coords(f)
                omega, phi, psi = calc_omega_phi_psi(mol_coords)
                coords_with_element[mol_id] = ele_coords

                out_f.write(','.join([str(z_layer),str(6 - i), str(mol_id), str(omega), str(phi), str(psi)]) + '\n')
                omega_list[5 - i].append(omega)
                phi_list[5 - i].append(phi)
                psi_list[5 - i].append(psi)

                if psi_threshold and psi <= psi_threshold:
                    extract_mol_ids.append(mol_id)

    if f.readline().strip() != 'END':
        raise Exception('Number of unit cells not correct!')

    result = ','.join(['{}_l{}'.format(angle, layer) for angle in ('omega','phi','psi') for layer in range(1, 4)]) + '\n'

    omega_list = list(zip_longest(*omega_list))
    phi_list = list(zip_longest(*phi_list))
    psi_list = list(zip_longest(*psi_list))

    assert len(omega_list) == len(phi_list) == len(psi_list), "length of angle lists not equal"

    for i in range(len(omega_list)):
        result += ','.join([str(a) for a in omega_list[i]] + [str(a) for a in phi_list[i]] + [str(a) for a in psi_list[i]]) + '\n'

    print(result)

    if psi_threshold:
        extract_f = open(output_name.rsplit('.', 1)[0]+'.xyz','w')
        extract_f.write(str(len(extract_mol_ids) * 26) + '\n')
        extract_f.write('frame: 1\n')

        for mol_id in extract_mol_ids:
            for id_atom in range(26):
                extract_f.write(' '.join(coords_with_element[mol_id][id_atom]) + '\n')
        extract_f.close()


    f.close()
    out_f.close()
    #for i in range(3):
    #    print(len(omega_list[i]), len(phi_list[i]), len(psi_list[i]))

    #print(np.histogram(omega_list))
    #print(np.histogram(phi_list))
    #print(np.histogram(psi_list))
