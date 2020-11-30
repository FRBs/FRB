""" Sets parameters for FRB associations

At the moment the Zero Points used include Galactic Extinction!

"""
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
    npixels=9,
)


## ##############################
## FRB 190608
#"""
#Notes:
#"""
updates = dict(
    name='FRB190608',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Unpublished', 'HG_190608_FORS2_I.fits'),
    cut_size = 34.,
    filter = 'VLT_FORS2_I',
    ZP = 27.9,  # Tied to Pan-Starrs i-band
    deblend=True,
    cand_bright=15.,
    cand_separation=10*units.arcsec,
    plate_scale=0.252202*units.arcsec,
)
frb190608 = {**base_config, **updates}  # Use | in 3.9


# ##############################
# FRB 121102
""" 
"""
updates = dict(
    name='FRB121102',
    image_file=os.path.join(gdb_path, 'Repeater', 'Tendulkar2017', 'FRB121102_GMOS_N_i.fits'),
    cut_size = 30.,
    filter = 'GMOS_N_i',
    ZP = 32.2,
    deblend=True,
    cand_bright=15.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.1459 * units.arcsec,
)
frb121102 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 180916
""" 
"""
updates = dict(
    name='FRB180916',
    image_file=os.path.join(gdb_path, 'CHIME', 'Marcote2020', 'FRB180916_GMOS_N_r.fits'),
    cut_size = 40.,
    filter = 'GMOS_N_r',
    ZP = 31.2,
    deblend=False,
    cand_bright=15.,
    cand_separation=20*units.arcsec,
    plate_scale = 0.1616 * units.arcsec,
)
frb180916 = {**base_config, **updates}  # Use | in 3.9


# ##############################
# FRB 180924
updates = dict(
    name='FRB180924',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Bannister2019', 'FRB180924_VLT_FORS2_g.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_g',
    ZP = 34.5,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.25226 * units.arcsec,
)
frb180924 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 181112
updates = dict(
    name='FRB181112',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Prochaska2019', 'FRB181112_VLT_FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 32.3,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.25226 * units.arcsec,
)
frb181112 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190102
updates = dict(
    name='FRB190102',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Macquart2020', 'FRB190102_VLT_FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 26.9,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.252 * units.arcsec,
)
frb190102 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190523
updates = dict(
    name='FRB190523',
    image_file = os.path.join(gdb_path, 'DSA', 'Ravi2019', 'FRB190523_LRIS_R.fits'),
    cut_size = 30.,
    filter = 'LRIS_R',
    ZP = 33.,
    deblend=True,
    plate_scale = 0.28 * units.arcsec,
    cand_separation=10*units.arcsec,
)
frb190523 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190611
updates = dict(
    name='FRB190611',
    image_file = os.path.join(gdb_path, 'CRAFT', 'Macquart2020', 'FRB190611_GMOS_S_i.fits'),
    cut_size = 30.,
    filter='GMOS_S_i',
    ZP=32.8,
    deblend=True,
    plate_scale = 0.160 * units.arcsec,
    cand_separation=10*units.arcsec,
)
frb190611 = {**base_config, **updates}

# ##############################
# FRB 190614
updates = dict(
    name='FRB190614',
    image_file = os.path.join(gdb_path, 'Realfast', 'Law2020', 'FRB190614_LRIS_I.fits'),
    cut_size = 30.,
    filter='LRIS_I',
    ZP=34.,
    deblend=True,
    plate_scale = 0.135 * units.arcsec,
    cand_separation=10*units.arcsec,
)
frb190614 = {**base_config, **updates}

# ##############################
# FRB 190711
"""
Notes:
   Figure out if the 2nd source is a star (HST data)
"""
updates = dict(
    name='FRB190711',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Macquart2020', 'FRB190711_GMOS_S_i.fits'),
    cut_size = 30.,
    filter = 'GMOS_S_i',
    ZP = 32.7,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.160*units.arcsec,
)
frb190711 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190714
"""
Notes:
"""
updates = dict(
    name='FRB190714',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Heintz2020', 'FRB190714_VLT_FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 27.39,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale=0.252*units.arcsec,
)
frb190714 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 191001
"""
Notes:
"""
updates = dict(
    name='FRB191001',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Heintz2020', 'FRB191001_VLT_FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 27.5,
    deblend=True,
    cand_bright=17.,
    cand_separation=10*units.arcsec,
    plate_scale=0.252*units.arcsec,
)
frb191001 = {**base_config, **updates}  # Use | in 3.9

 ##############################
# FRB 200430
"""
Notes:
"""
updates = dict(
    name='FRB200430',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Heintz2020', 'FRB200430_LRIS_I.fits'),
    cut_size = 30.,
    filter = 'LRIS_I',
    ZP = 34.2,  # Tied to Pan-Starrs i-band
    deblend=True,
    cand_bright=17.,
    cand_separation=10*units.arcsec,
    plate_scale=0.134*units.arcsec,
)
frb200430 = {**base_config, **updates}  # Use | in 3.9
#
#