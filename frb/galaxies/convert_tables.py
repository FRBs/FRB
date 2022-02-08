""" Methods to convert data into Literature tables"""
import os
from pkg_resources import resource_filename

import pandas

from frb import frb
from frb import defs
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
        frb121102, json_file, cosmo=defs.frb_cosmo) 
    # Build  Table
    df = pandas.DataFrame()
    for key in host.neb_lines.keys():
        df[key] = [host.neb_lines[key]]
    # Add
    df['ra'] = host.coord.ra.deg
    df['dec'] = host.coord.dec.deg
    df['Name'] = host.name

    # Write
    outfile = os.path.join(
        resource_filename('frb', 'data'),
        'Galaxies', 'Literature', 'tendulkar2017_nebular.csv')
    df.to_csv(outfile)

if __name__ == '__main__':
    tendulkar_nebular()