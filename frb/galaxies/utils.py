""" Utilities related to FRB galaxies"""

import os
import glob
from IPython import embed

try:
    from specdb.specdb import SpecDB
except ImportError:
    flg_specdb = False
else:
    flg_specdb = True


def load_specdb(specdb_file=None):
    """
    Automatically load the specDB file from $SPECDB/FRB_specDB.hdf5

    Args:
        specdb_file (str, optional):
            Over-ride the default file

    Returns:
        specdb.specdb.SpecDB:

    """
    if not flg_specdb:
        raise IOError("You must install the specdb package first!")
        return
    if specdb_file is None:
        if os.getenv('SPECDB') is None:
            raise IOError("You must set the SPECDB environmental variable")
        specdb_files = glob.glob(os.path.join(os.getenv('SPECDB'), 'FRB_specDB_*.hdf5'))
        if len(specdb_files) > 0:
            specdb_file = specdb_files[0]
            print("Loading spectra from {:s}".format(specdb_file))
        else:
            raise IOError("There are no FRB_specdb.hdf5 files in your SPECDB folder")
    # Load it up
    specDB = SpecDB(db_file=specdb_file)
    # Return
    return specDB
