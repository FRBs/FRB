""" Module for RM calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os

import importlib_resources

from astropy import units

import healpy as hp

def load_opperman2014():
    """
    Load the Oppermann et al. 2014 -- https://arxiv.org/abs/1404.3701
    maps Galactic Farady RM 
    and its uncertainty 

    Returns:
        healpy map, healpy map: RM and RM_err with units of rad/m^2

    """
    print("Loading RM information map from Oppermann et al. 2014")
    galactic_rm_file = importlib_resources.files('frb.data.RM')/'opp14_foreground.fits'

    # Load
    rm_sky = hp.read_map(galactic_rm_file, hdu=4)
    sig_sky = hp.read_map(galactic_rm_file, hdu=6)
    # 
    return rm_sky, sig_sky

def load_hutschenreuter2020():
    """
    Load the Hutschenreuter 2020 map -- https://wwwmpa.mpa-garching.mpg.de/~ensslin/research/data/faraday2020.html
    Galactic Farady RM and its uncertainty 

    See: https://ui.adsabs.harvard.edu/abs/2022A%26A...657A..43H/abstract
    for full details

    Returns:
        healpy map, healpy map: RM and RM_err with units of rad/m^2

    """
    print("Loading RM information map from Hutschenreuter et al. 2022, aka faraday2020v2.fits")
    galactic_rm_file = importlib_resources.files('frb.data.RM')/'faraday2020v2.fits'

    # Has it been downloaded?
    if not os.path.isfile(galactic_rm_file):
        readme_file = importlib_resources.files('frb.data.RM')/'README'
        print(f"See the README here: {readme_file}")
        raise IOError("You need to download the Hutschenreuter 2022 map to proceed, aka faraday2020v2.fits. https://wwwmpa.mpa-garching.mpg.de/ift/data/faraday2020/faraday2020v2.fits")

    # Load
    rm_sky = hp.read_map(galactic_rm_file)
    sig_sky = hp.read_map(galactic_rm_file, field=1)
    # 
    return rm_sky, sig_sky

def galactic_rm(coord, use_map=2020):
    """
    Provide an RM and error estimate for a coordinate
    on the sky. 

    Default is the new Hutschenreuter 2020 map 

    Args:
        coord (astropy.coordinates.SkyCoord): 
            Coordinate for the RM esimation
        use_map (int, optional):
            Specifies the map to use.  Options are [2014, 2020]
            Default is 2020

    Returns:
        Quantity, Quantity: RM and RM_err with units of rad/m^2

    """
    if use_map == 2014:
        rm_sky, sig_sky = load_opperman2014()
    elif use_map == 2020:
        rm_sky, sig_sky = load_hutschenreuter2020()
    else:
        raise IOError("Bad use_map input.  Allowed choices are [2014, 2020]")

    # Load
    nside = hp.get_nside(rm_sky)

    # Find the pixel
    pix = hp.ang2pix(nside, coord.galactic.l.value, coord.galactic.b.value, lonlat=True)

    # Return
    return rm_sky[pix]*units.rad/units.m**2, sig_sky[pix]*units.rad/units.m**2
