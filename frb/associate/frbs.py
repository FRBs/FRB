""" Sets parameters for FRB associations

At the moment the Zero Points used include Galactic Extinction!

"""
import os
from astropy import units
import warnings

from IPython import embed

if os.getenv('FRB_GDB') is None:
    warnings.warn("FRB_GDB variable is not set.  Odds are you are doing something wrong..")
    gdb_path = ''
else:
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
    name='FRB20190608',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Unpublished', 'HG_190608_FORS2_I.fits'),
    cut_size = 34.,
    filter = 'VLT_FORS2_I',
    ZP = 27.9,  # Tied to Pan-Starrs i-band
    deblend=True,
    cand_bright=15.,
    cand_separation=10*units.arcsec,
    plate_scale=0.252202*units.arcsec,
)
frb20190608 = {**base_config, **updates}  # Use | in 3.9


# ##############################
# FRB 121102
""" 
"""
updates = dict(
    name='FRB20121102',
    image_file=os.path.join(gdb_path, 'Repeater', 'Tendulkar2017', 'FRB20121102_GMOS_N_i.fits'),
    cut_size = 30.,
    filter = 'GMOS_N_i',
    ZP = 32.2,
    deblend=True,
    cand_bright=15.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.1459 * units.arcsec,
)
frb20121102 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 180916
""" 
"""
updates = dict(
    name='FRB20180916',
    image_file=os.path.join(gdb_path, 'CHIME', 'Marcote2020', 'FRB20180916_GMOS_N_r.fits'),
    cut_size = 40.,
    filter = 'GMOS_N_r',
    ZP = 31.2,
    deblend=False,
    cand_bright=15.,
    cand_separation=20*units.arcsec,
    plate_scale = 0.1616 * units.arcsec,
)
frb20180916 = {**base_config, **updates}  # Use | in 3.9


# ##############################
# FRB 180924
updates = dict(
    name='FRB20180924',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Bannister2019', 'FRB20180924_VLT_FORS2_g.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_g',
    ZP = 34.5,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.25226 * units.arcsec,
)
frb20180924 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 181112
updates = dict(
    name='FRB20181112',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Prochaska2019', 'FRB20181112_VLT_FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 32.3,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.25226 * units.arcsec,
)
frb20181112 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190102
updates = dict(
    name='FRB20190102',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Macquart2020', 'FRB20190102_VLT_FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 26.9,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.252 * units.arcsec,
)
frb20190102 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190523
updates = dict(
    name='FRB20190523',
    image_file = os.path.join(gdb_path, 'DSA', 'Ravi2019', 'FRB20190523_LRIS_R.fits'),
    cut_size = 30.,
    filter = 'LRIS_R',
    ZP = 33.,
    deblend=True,
    plate_scale = 0.28 * units.arcsec,
    cand_separation=10*units.arcsec,
)
frb20190523 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190611
updates = dict(
    name='FRB20190611',
    image_file = os.path.join(gdb_path, 'CRAFT', 'Macquart2020', 'FRB20190611_GMOS_S_i.fits'),
    cut_size = 30.,
    filter='GMOS_S_i',
    ZP=32.8,
    deblend=True,
    plate_scale = 0.160 * units.arcsec,
    cand_separation=10*units.arcsec,
)
frb20190611 = {**base_config, **updates}

# ##############################
# FRB 190614
updates = dict(
    name='FRB20190614',
    image_file = os.path.join(gdb_path, 'Realfast', 'Law2020', 'FRB20190614_LRIS_I.fits'),
    cut_size = 30.,
    filter='LRIS_I',
    ZP=34.,
    deblend=True,
    plate_scale = 0.135 * units.arcsec,
    cand_separation=10*units.arcsec,
)
frb20190614 = {**base_config, **updates}

# ##############################
# FRB 190711
"""
Notes:
   Figure out if the 2nd source is a star (HST data)
"""
updates = dict(
    name='FRB20190711',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Macquart2020', 'FRB20190711_GMOS_S_i.fits'),
    cut_size = 30.,
    filter = 'GMOS_S_i',
    ZP = 32.7,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.160*units.arcsec,
)
frb20190711 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 190714
"""
Notes:
"""
updates = dict(
    name='FRB20190714',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Heintz2020', 'FRB20190714_VLT_FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 27.39,
    deblend=True,
    cand_bright=18.,
    cand_separation=10*units.arcsec,
    plate_scale=0.252*units.arcsec,
)
frb20190714 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 191001
"""
Notes:
"""
updates = dict(
    name='FRB20191001',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Heintz2020', 'FRB20191001_VLT-FORS2_I_BESS.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 27.5,
    deblend=True,
    cand_bright=17.,
    cand_separation=10*units.arcsec,
    plate_scale=0.252*units.arcsec,
)
frb20191001 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 191228
"""
Notes:
"""
updates = dict(
    name='FRB20191228',
    image_file=os.path.join(gdb_path, 'Realfast', 'Bhandari2021', 
                            'FRB20191228_VLT-FORS2_I.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_I',
    ZP = 27.453 - 0.039, # Reported by FORS2 QC1 Archive + IRSA Dust Tool
    deblend=True,
    cand_bright=None,
    cand_separation=10*units.arcsec,
    plate_scale=0.252*units.arcsec,
)
frb20191228 = {**base_config, **updates}  # Use | in 3.9

##############################
# FRB 200430
"""
Notes:
"""
updates = dict(
    name='FRB20200430',
    image_file=os.path.join(gdb_path, 'CRAFT', 'Heintz2020', 'FRB20200430_LRIS_I.fits'),
    cut_size = 30.,
    filter = 'LRIS_I',
    ZP = 34.2,  # Tied to Pan-Starrs i-band
    deblend=True,
    cand_bright=17.,
    cand_separation=10*units.arcsec,
    plate_scale=0.134*units.arcsec,
)
frb20200430 = {**base_config, **updates}  # Use | in 3.9

# ##############################
# FRB 200906
"""
Notes:
"""
updates = dict(
    name='FRB20200906',
    image_file=os.path.join(gdb_path, 'Realfast', 
                            'Bhandari2021', 
                            'FRB20200906_VLT-FORS2_g.fits'),
    cut_size = 30.,
    filter = 'VLT_FORS2_g',
    ZP = 27.42 - 0.016, # Tied to DES g-band
    deblend=True,
    cand_bright=None,
    cand_separation=10*units.arcsec,
    plate_scale=0.252*units.arcsec,
)
frb20200906 = {**base_config, **updates}  # Use | in 3.9

##############################
# FRB 180301
"""
Notes:
"""
updates = dict(
    name='FRB20180301',
    image_file=os.path.join(gdb_path, 'Realfast', 
                            'Bhandari2021', 'FRB20180301_GMOS_S_r.fits'),
    cut_size = 30.,
    filter = 'GMOS_S_r',
    ZP = 32.94,  # Kasper
    deblend=True,
    cand_bright=17.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.1616 * units.arcsec,
)
frb20180301 = {**base_config, **updates}  # Use | in 3.9
#
#

##############################
# FRB 201124
updates = dict(
    name='FRB20201124A',
    image_file=os.path.join(gdb_path, 'F4', 
                            'fong2021', 'FRB20201124_Pan-STARRS_r.fits'),
    cut_size = 30.,
    filter = 'Pan-STARRS_r',
    ZP = 32.29,  # Refined using catalog
    deblend=True,
    cand_bright=17.,
    cand_separation=10*units.arcsec,
    plate_scale = 0.250 * units.arcsec,
)
frb20201124 = {**base_config, **updates}  # Use | in 3.9
#
#
