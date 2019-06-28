""" Utilities related to FRB galaxies"""

import os
import warnings
from IPython import embed

try:
    from specdb.specdb import SpecDB
except ImportError:
    flg_specdb = False
else:
    flg_specdb = True

def load_specdb(specdb_file=None):
    if not flg_specdb:
        warnings.warn("You must install the specdb package first!")
        return
    if specdb_file is None:
        specdb_file = os.path.join(os.getenv('SPECDB'), 'FRB_specDB.hdf5')
    # Load it up
    specDB = SpecDB(db_file=specdb_file)
    # Return
    return specDB
