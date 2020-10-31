import os
import numpy as np

from astropy.nddata import Cutout2D
from astropy import units
from astropy.io import fits
from astropy.wcs import WCS

from frb.associate import frbassociate
from frb import frb

from IPython import embed


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
    if image_file is None:
        gdb_path = os.getenv('FRB_GDB')
        image_file = os.path.join(gdb_path, 'CRAFT', 'Bannister2019', 'FRB180924_VLT_FORS2_g.fits')
        ifilter = 'g'

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
    plate_scale = 0.25226 * units.arcsec
    frbA_180924.cut_candidates(plate_scale, bright_cut=18., separation=10 * units.arcsec)
    if verbose:
        print(frbA_180924.candidates[['id', ifilter, 'half_light', 'separation']])

    frbA_180924.calc_pchance()

    # Return here?
    if skip_bayesian:
        return frbA_180924

    # Finish
    return frbA_180924


def run_frb190523(show=False,
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
    gdb_path = os.getenv('FRB_GDB')
    image_file = os.path.join(gdb_path, 'DSA', 'Ravi2019', 'FRB190523_Keck_LRIS_R.fits')
    ifilter = 'R'

    # FRB 180924
    frb190523 = frb.FRB.by_name('FRB190523')

    # Image
    hdul = fits.open(image_file)
    hdu_full = hdul[0]

    size = units.Quantity((cut_size, cut_size), units.arcsec)
    cutout = Cutout2D(hdu_full.data, frb190523.coord, size, wcs=WCS(hdu_full.header))

    # FRB Associate
    frbA_190523 = frbassociate.FRBAssociate(frb190523, max_radius=max_radius)

    # Sub image
    frbA_190523.wcs = cutout.wcs
    frbA_190523.hdu = cutout  # not really an HDU
    frbA_190523.header = hdu_full.header

    # Threshold + Segment
    frbA_190523.threshold()
    frbA_190523.segment(deblend=True, npixels=9, show=show)

    # Photometry
    frbA_190523.photometry(33, ifilter, show=show)
    if verbose:
        print(frbA_190523.photom[19:22]['xcentroid', 'ycentroid', ifilter])

    # Candidates
    plate_scale = 0.28 * units.arcsec
    frbA_190523.cut_candidates(plate_scale, separation=10 * units.arcsec)

    # Chance probability
    frbA_190523.calc_pchance()
    if verbose:
        print(frbA_190523.candidates[['id', ifilter, 'half_light', 'separation', 'P_c']])

    # Return here?
    if skip_bayesian:
        return frbA_190523

    # Finish
    return frbA_190523
