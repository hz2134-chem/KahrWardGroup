#!/usr/bin/env python3

import argparse

def load_mol_pe(pe_filename, pe_origin):
    pe = open(pe_filename, 'r')
    mol_ct = 0
    mol_pe = []
    for line in pe:
        if line[0] == '#':
            continue

        fields = line.strip().split()
        if len(fields) == 3:            # header of each frame
            mol_ct += int(fields[1])
        else:
            # Parse pe in pecmol.dat
            # Notice pecmol.dat contains average pe/atom within a molecule, thus need to x 26
            mol_pe.append(float(fields[5]) * 26 - pe_origin)

    pe.close()
    assert len(mol_pe) == mol_ct, f"Molecule count doesn't match!"

    return mol_pe

    
parser = argparse.ArgumentParser()
parser.add_argument('size', type=int)
parser.add_argument('layer', type=int)
parser.add_argument('--traj_filename', default='trajectory.xyz')
parser.add_argument('--pe_filename', default='pecmol.dat')
parser.add_argument('--pe_origin', type=float, default=1.25073)

args = parser.parse_args()
size = args.size
layer = args.layer
traj_filename = args.traj_filename
pe_filename = args.pe_filename
pe_origin = args.pe_origin

one_layer_mols = 3* ((size - 1) * size * 3 + 1)
skipped_mols = one_layer_mols * (layer - 1)
skipped_mols_after_select = one_layer_mols * (15 - layer)

traj = open(traj_filename,'r')
output = open('output.xyz','w')

pe_mols = load_mol_pe(pe_filename, pe_origin)

idx_atm, idx_mol = 0, 0
timestep, header_timestep = 0, False

for line in traj:
    fields = line.strip().split()
    if len(fields) != 4:
        output.write(f"{one_layer_mols*26}\n")
        header_timestep = not header_timestep
        continue
        
    elif not header_timestep:
        output.write(timestep)
        timestep += 10000 # Depends on what is your output frequency
        pre_select_lines, select_lines, post_select_lines = 0, 0, 0

        header_timestep = not header_timestep
        continue

    elif pre_select_lines < skipped_mols * 26:
        pre_select_lines += 1

    elif select_lines < one_layer_mols * 26:
        select_lines += 1
        out_line = ' '.join(fields) + f' {pe_mols[idx_mol]:.4f}\n'
        output.write(out_line)

    elif post_select_lines < skipped_mols_after_select * 26:
        post_select_lines += 1

    else:
        raise Exception("Should not reach here")

    idx_atm = (idx_atm + 1) % 26
    if idx_atm == 0:
        idx_mol += 1

assert idx_mol == len(pe_mols), f"PE file not end!\nidx_mol: {idx_mol}, len(pe_mols): {len(pe_mols)}"

traj.close()
output.close()
