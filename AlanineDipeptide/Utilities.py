import sys

import numpy as np
import pandas as pd
from tqdm import tqdm


def save_frames_as_pdb(frames, u):
    # save all frames in pdb files
    for idx in tqdm(frames.keys(), desc="Saving GRO files for each microstate"):
        u.trajectory[idx]
        u.atoms.write('../data/frames/f' + str(idx) + '.gro')


def csv_print(mat):
    np.savetxt(sys.stdout, mat, fmt='%.5f', newline="\n")


def read_dihedral_angles(phi_path, psi_path, frames=1, phi_offset=60, psi_offset=-90):
    phi_df = pd.read_csv(phi_path, header=None,
                         index_col='frame', names=['frame', 'phi'],
                         delimiter='\t').head(frames)
    psi_df = pd.read_csv(psi_path, header=None,
                         index_col='frame', names=['frame', 'psi'],
                         delimiter='\t').head(frames)

    phi_df = (((phi_df + 180) + phi_offset) % 360) - 180
    psi_df = (((psi_df + 180) + psi_offset) % 360) - 180

    both_angles = pd.concat([phi_df['phi'], psi_df['psi']], axis=1, keys=['phi', 'psi']) * -1

    return pd.DataFrame(both_angles)


def read_xvg(fname):
    """Read columns of data from file fname

    Returns numpy array of data
    """
    skip = 0
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if not line.startswith(('#', '@')):
                skip = i
                break

    return np.genfromtxt(fname, skip_header=skip, usecols=(0, 1))


def apply_offset(df, phi_offset=60, psi_offset=-90, phi_mult=1, psi_mult=1):
    phi = df[0]
    psi = df[1]

    phi = (((phi + 180) + phi_offset) % 360) - 180
    psi = (((psi + 180) + psi_offset) % 360) - 180

    phi *= phi_mult
    psi *= psi_mult

    return pd.DataFrame(pd.concat([phi, psi], axis=1))
