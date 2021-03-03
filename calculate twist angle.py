#!/usr/bin/env python3

import argparse
from math import sqrt, acos, degrees

def read_uc(f, uc_mol):
    result = []
    for _ in range(uc_mol):
        result.append(f.readline())
    return result

def get_xyz_from_uc(uc):
    result = [0, 0, 0]
    atm_ct = len(uc)
    for i in range(atm_ct):
        coord = uc[i].split()[5:8]
        for c in range(3):
            result[c] += float(coord[c])
    for c in range(3):
        result[c] /= atm_ct
    return result

def skip_n_layers(f, n, uc_in_layer, uc_mol):
    for _ in range(n * uc_in_layer * uc_mol):
        f.readline()

def extract_coord_from_layer(f, skipCircle, uc_mol):
    result = []
    skip_uc_in_layer = 1 + 6 * ((1 + skipCircle) * skipCircle // 2)
    # Skip unit cells not considered in the layer of interest
    for _ in range(skip_uc_in_layer * uc_mol):
        f.readline()
    # Obtain the COM of six unit cell one circle below surface
    for i in range(6):
        uc = read_uc(iFile, uc_mol)
        result.append(get_xyz_from_uc(uc))
        for _ in range(skipCircle * uc_mol):
            f.readline()
    # Complete reading the surface circle
    for _ in range(6 * (skipCircle + 2) * uc_mol):
        f.readline()
    assert len(result) == 6, "Didn't extract all six corners!"
    return result

def vec_len(vec):
    return sqrt(sum([x*x for x in vec]))


def mean(angle_list):
    num_angles = len(angle_list)
    assert bool(num_angles), "Error: an empty list is provided!"

    return sum(angle_list) / num_angles
    

def stddev(angle_list):
    '''
    Given a list, calculate its corrected sample standard deviation
    @ref: https://en.wikipedia.org/wiki/Bessel%27s_correction
    '''
    num_angles = len(angle_list)
    mean_angle = mean(angle_list)

    result = 0
    for angle in angle_list:
        diff = angle - mean_angle
        result += diff * diff

    result /= (num_angles - 1)

    return sqrt(result)


def confidence_interval(mean, stddev, num_item):
    '''
    Return a tuple that contains mean and error,
    as well as low and high ends of 95% confidence_interval
    '''
    error = 1.96 * stddev / sqrt(num_item)

    return (mean, error, mean - error, mean + error)

if __name__ == "__main__":

    parser = argparse.ArgumentParser('Calculate the twisting angle of a hexagonal crystal, with or without dislocation.\n./calcTwist.py unitCell.pdb 9 15')
    parser.add_argument("file", help='The pdb file name for the molecule COMs')
    parser.add_argument("replL", type=int, help='The repetition times in a (or b)')
    parser.add_argument("replH", type=int, help='The repetition times in c')
    parser.add_argument("uc_mol_ct", type=int, help='Provide numer of molecules in a unit cell.')

    args = parser.parse_args()
    filename = args.file.split('.')[0]
    path_filename = filename.rsplit('/', 1)
    path, filename = path_filename[0], path_filename[-1]
    path = '.' if path == filename else path
    l, h = args.replL, args.replH
    uc_mol = args.uc_mol_ct
    layerL, layerH = h // 4, (h * 3) // 4 # layers to be extracted
    skipCircle = l - 3          # Circles on the layer to be extracted
    uc_in_layer = 1 + 6 * (l * (l - 1) // 2)

    iFile = open(args.file, 'r')
    #trajFile = open(f'{path}/traj_{filename}.xyz','w')
    dataFile = open(f'{path}/angle_{filename}.txt','w')

    time = -1
    angle_list = []
    while iFile.readline():
        iFile.readline()
        coord_considered = []
        skip_n_layers(iFile, layerL, uc_in_layer, uc_mol)
        coord_considered.append(extract_coord_from_layer(iFile, skipCircle, uc_mol))
        skip_n_layers(iFile, (layerH - layerL - 1), uc_in_layer, uc_mol)
        coord_considered.append(extract_coord_from_layer(iFile, skipCircle, uc_mol))
        skip_n_layers(iFile, (h - layerH - 1), uc_in_layer, uc_mol)
        if iFile.readline().strip() != 'END':
            raise Exception('Number of unit cells not correct!')

        # write coordinates into trajFile
        #trajFile.write('12\n')
        #trajFile.write(f'Timestep: {time}\n')
        #for ly in coord_considered:
        #    for unc in ly:
        #        trajFile.write('C {:>8.3f}{:>8.3f}{:>8.3f}\n'.format(*unc))

        angle_Z = 0
        for i in range(6):
            crd = coord_considered
            vec1 = [(crd[0][(i + 1) % 6][x] - crd[0][i][x]) for x in range(2)]
            vec2 = [(crd[1][(i + 1) % 6][x] - crd[1][i][x]) for x in range(2)]
            deltaZ = (crd[1][(i + 1) % 6][2] + crd[1][i][2]) / 2 - (crd[0][(i + 1) % 6][2] + crd[0][i][2]) / 2
            dot_prod = (vec1[0]*vec2[0] + vec1[1]*vec2[1])/(vec_len(vec1)*vec_len(vec2))
            dot_prod = 1 if dot_prod > 1 else dot_prod
            angle = degrees(acos(dot_prod))
            angle /= deltaZ
            angle_list.append(angle)
            angle_Z += angle
        angle_Z /= 6
        dataFile.write(f'{time} {angle_Z:6.4f}\n')
        time += 1

    sample_list = angle_list[-60:]
    num_angle = len(sample_list)
    print(f"Number of angles used in calculation: {num_angle}")
    
    mean_angle = mean(sample_list)
    stddev_angle = stddev(sample_list)
    print(f"Mean and stddev of twisting angle: {mean_angle:.5f}\t{stddev_angle:.5f}")


    iFile.close()
    #trajFile.close()
    dataFile.close()

