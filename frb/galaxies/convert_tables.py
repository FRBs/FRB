""" Methods to convert data into Literature tables"""
import os
from pkg_resources import resource_filename

from astropy.cosmology import Planck15

from frb import frb
from frb.galaxies import frbgalaxy

from IPython import embed

def tendulkar_nebular():
    frb121102 = frb.FRB.by_name('FRB20121102')
    # Old host file
    json_file = os.path.join(
        resource_filename('frb', 'data'),
        'Galaxies', '121102',
        'FRB121102_host.json')
    host = frbgalaxy.FRBHost.from_json(
        frb121102, json_file, cosmo=Planck15)
    embed(header='16 of convert')
    # Build  Table

if __name__ == '__main__':
    tendulkar_nebular()