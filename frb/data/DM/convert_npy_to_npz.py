""" Simple script to convert .npy files to .npz files. """

import numpy as np

from frb.dm import prob_dmz

def main():

    npy_files = ['DSA_pzdm.npy',
        'parkes_mb_class_I_and_II_pzdm.npy',
        'CRAFT_class_I_and_II_pzdm.npy',
        'CRAFT_ICS_1300_pzdm.npy',
        'CRAFT_ICS_892_pzdm.npy',
        'CRAFT_ICS_1632_pzdm.npy',
        'FAST_pzdm.npy']

    # Load of z, DM_EG from CHIME
    chime_dict = prob_dmz.grab_repo_grid(prob_dmz.telescope_dict['CHIME'])

    # Parse
    z = chime_dict['z']
    DM = chime_dict['DM']

    # Loop on the files
    for npy_file in npy_files:
        # Load
        pzdm = prob_dmz.grab_repo_grid(npy_file)
        # Save
        npz_file = npy_file.replace('.npy', '.npz')
        np.savez(npz_file, pzdm=pzdm, z=z, DM=DM)
        # Report
        print(f"Wrote: {npz_file}")

# Command line execution
if __name__ == '__main__':
    main()