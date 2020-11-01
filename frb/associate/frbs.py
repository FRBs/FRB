import os
import numpy as np

from astropy.nddata import Cutout2D
from astropy import units
from astropy.io import fits
from astropy.wcs import WCS

from frb.associate import frbassociate
from frb import frb

from IPython import embed

gdb_path = os.getenv('FRB_GDB')

base_config = dict(
    max_radius=10.,
    cut_size=None,
    npixels=5,
    deblend=False,
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
    cand_separation=10 * units.arcsec,
    npixels=9,
)
frb190523 = {**base_config, **updates}  # Use | in 3.9

def run_frb180924(image_file=None, show=False,
                      skip_bayesian=False,
                      max_radius=10.,
                      verbose=False, cut_size=30.):
    """

    Args:
        image_file:
        show:
        skip_bayesian:
        max_radius:
        verbose:
        cut_size:

    Returns:
        frbassociate.FRBAssociate

    """
    # FRB 180924
    frb180924 = frb.FRB.by_name('FRB180924')

    # Image
    hdul = fits.open(image_file)
    hdu_full = hdul[0]

    size = units.Quantity((cut_size, cut_size), units.arcsec)
    cutout = Cutout2D(hdu_full.data, frb180924.coord, size, wcs=WCS(hdu_full.header))

    # FRB Associate
    frbA_180924 = frbassociate.FRBAssociate(frb180924, max_radius=max_radius)

    # Sub image
    frbA_180924.wcs = cutout.wcs
    frbA_180924.hdu = cutout  # not really an HDU
    frbA_180924.header = hdu_full.header

    # Threshold + Segment
    frbA_180924.threshold()
    frbA_180924.segment(deblend=True, npixels=9, show=show)

    # Photometry
    frbA_180924.photometry(34.5, ifilter, show=show)
    if verbose:
        print(frbA_180924.photom[19:22]['xcentroid', 'ycentroid', 'g'])

    # Candidates
    frbA_180924.cut_candidates(plate_scale, bright_cut=18., separation=10 * units.arcsec)
    if verbose:
        print(frbA_180924.candidates[['id', ifilter, 'half_light', 'separation']])

    frbA_180924.calc_pchance()

    # Return here?
    if skip_bayesian:
        return frbA_180924

    # Finish
    return frbA_180924

