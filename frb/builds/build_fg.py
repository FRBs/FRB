""" Top-level module to build or re-build the JSON files for
FRB host galaxies"""

from pkg_resources import resource_filename
import os
import sys

from IPython import embed

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from frb.galaxies import frbgalaxy, defs
from frb.surveys import des
from frb.frb import FRB

db_path = os.getenv('FRB_GDB')
if db_path is None:
    embed(header='You need to set $FRB_GDB')


def build_fg_181112():
    # Coord from DES
    fg_coord = SkyCoord('J214923.89-525810.43', unit=(units.hourangle, units.deg))  # from DES
    frb181112 = FRB.by_name('FRB181112')

    # Instantiate
    fg_13_5 = frbgalaxy.FGGalaxy(fg_coord.ra.value, fg_coord.dec.value, '181112')
    fg_13_5.frb_coord = frb181112.coord

    # Redshift
    fg_13_5.set_z(0.36738, 'spec', err=7e-5)

    # Photometry

    # DES
    # Grab the table (requires internet)
    search_r = 2 * units.arcsec
    des_srvy = des.DES_Survey(fg_coord, search_r)
    des_tbl = des_srvy.get_catalog(print_query=True)

    fg_13_5.parse_photom(des_tbl)

    # VLT -- Lochlan 2019-05-02
    # VLT -- Lochlan 2019-06-18
    fg_13_5.photom['VLT_g'] = 21.20
    fg_13_5.photom['VLT_g_err'] = 0.05
    fg_13_5.photom['VLT_I'] = 19.20
    fg_13_5.photom['VLT_I_err'] = 0.02

    # Nebular lines
    #fg_13_5.parse_ppxf('fg_FORS2_ppxf_results.ecsv')
    fg_13_5.parse_ppxf(os.path.join(db_path, 'CRAFT', 'Prochaska2019', 'FG181112_13_5_FORS2_ppxf.ecsv'))

    # Derived quantities
    fg_13_5.calc_nebular_AV('Ha/Hb')

    # This will be an upper limit
    fg_13_5.calc_nebular_SFR('Ha')
    fg_13_5.derived['SFR_nebular_err'] = -999.

    # CIGALE
    fg_13_5.parse_cigale(os.path.join(db_path, 'CRAFT', 'Prochaska2019', 'FG181112_13_5_CIGALE.fits'))

    # Write
    path = resource_filename('frb', 'data/Galaxies/181112')
    fg_13_5.write_to_json(path=path)
    
    
def main(inflg='all'):

    if inflg == 'all':
        flg = np.sum(np.array( [2**ii for ii in range(25)]))
    else:
        flg = int(inflg)

    # 121102
    if flg & (2**0):
        build_fg_181112()


