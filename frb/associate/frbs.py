""" Sets parameters for FRB associations """
import os
from astropy import units

from IPython import embed

gdb_path = os.getenv('FRB_GDB')

base_config = dict(
    max_radius=10.,
    cut_size=None,
    deblend=False,
    cand_bright=None,
    cand_separation=None,
    skip_bayesian=False,
)

# FRB 180924
updates = dict(
    name='FRB180924',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Bannister2019', 'FRB180924_VLT_FORS2_g.fits'),
    cut_size = 30.,
    filter = 'g',
    ZP = 34.5,
    deblend=True,
    npixels=9,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.25226 * units.arcsec,
)
frb180924 = {**base_config, **updates}  # Use | in 3.9

# FRB 190523
updates = dict(
    name='FRB190523',
    image_file = os.path.join(gdb_path, 'DSA', 'Ravi2019', 'FRB190523_Keck_LRIS_R.fits'),
    cut_size = 30.,
    filter = 'R',
    ZP = 33.,
    deblend=True,
    plate_scale = 0.28 * units.arcsec,
    cand_separation=10*units.arcsec,
    npixels=9,
)
frb190523 = {**base_config, **updates}  # Use | in 3.9

