""" Methods related to nebular line analysis, e.g. dust extinction, SFR"""
import pdb

from pkg_resources import resource_filename

import numpy as np
import requests
import warnings

from xml.etree import ElementTree as ET

from scipy.interpolate import interp1d

from astropy.table import Table
from astropy import units

try:
    from linetools.lists import linelist
except ImportError:
    warnings.warn("Galaxy nebular line analysis requires linetools.  Install it if you want to use them")

# GLOBALS
Ha_Hb_intrin = 2.8  # Osterbrock
Ha_conversion = 7.9e-42 * units.Msun/units.yr   # Kennicutt 1998

def load_extinction(curve):
    """
    Load an extinction curve

    This method may move elsewhere..

    Args:
        curve (str): Name of the extinction curve
          MW = Standard Cardelli

    Returns:
        scipy.interpolate.interp1d: Interpolation function for alambda/AV

    """
    # Load the extinction curve
    if curve == 'MW':
        dust_file = resource_filename('frb', 'data/Dust/MW_dust.dat')
        MW_dust = Table.read(dust_file, format='ascii')
    else:
        raise IOError("Not ready for this extinction curve!")
    # Generate function for interpolating
    alAV = interp1d(MW_dust['wave'], MW_dust['Al_AV'])
    # Return
    return alAV


def calc_dust_extinct(neb_lines, method, curve='MW'):
    """
    Estimate the Visual extinction A_V based on input nebular emission lines

    Args:
        neb_lines (dict):  Line fluxes
        method (str): Name of the method
          Ha/Hb -- Use the Halpha/Hbeta ratio and standard intrinsic flux
        curve (str): Extinction curve to use

    Returns:
        float: A_V in magnitudes

    """


    # Dust extinction curve
    alAV = load_extinction(curve)

    # Which ratio?
    if method == 'Ha/Hb':
        wave1 = 6564.6  # redder
        wave2 = 4862.7
        #
        F1_obs = neb_lines['Halpha']
        F2_obs = neb_lines['Hbeta']
        #
        pair = True
    else:
        print("Not prepared for this method of analysis: {}".format(method))
        raise IOError("See docs for available options")

    if not pair:
        raise IOError("Not ready for this mode")

    # Extinction
    a1AV = alAV(wave1)
    a2AV = alAV(wave2)

    # Observed ratio
    fratio_obs = F1_obs/F2_obs

    # Calculate using intrinsic ratio
    AV = 2.5 * np.log10(Ha_Hb_intrin/fratio_obs) / (a1AV - a2AV)

    # Return
    return AV


def calc_lum(neb_lines, line, z, cosmo, AV=None, curve='MW'):
    """
    Calculate the line luminosity (and error) from input nebular line emission

    Error is -999.*erg/s if input line flux has negative error

    Args:
        neb_lines (dict):  Observed line fluxes and errors
        line (str): Line to analyze
        z (float):  Emission redshift -- for Luminosity distance
        cosmo (astropy.cosmology.FLRW): Cosmology
        AV (float, optional):  Visual extinction, if supplied will apply
        curve (str):  Name of the extinction curve.  Only used if A_V is supplied

    Returns:
        Quantity, Quantity:  Luminosity, sig(Luminosity)
    """
    # Grab rest wavelength (vacuum)
    llist = linelist.LineList('Galaxy')
    wave = llist[line]['wrest']

    # Dust correct?
    if AV is not None:
        alAV = load_extinction(curve)
        al = AV * alAV(wave.to('Angstrom').value)
    else:
        al = 0.

    # Cosmology
    DL = cosmo.luminosity_distance(z)

    # Luminosity
    flux = neb_lines[line]
    Lum = flux * units.erg/units.s/units.cm**2 * 10**(al/2.5) * (4*np.pi * DL**2)

    # Error
    if neb_lines[line+'_err'] > 0.:
        flux_err = neb_lines[line+'_err']
        Lum_err = flux_err * units.erg/units.s/units.cm**2 * 10**(al/2.5) * (4*np.pi * DL**2)
    else:
        Lum_err = -999 * units.erg/units.s
    # Return
    return Lum.to('erg/s'), Lum_err.to('erg/s')


def calc_SFR(neb_lines, method, z, cosmo, AV=None, curve='MW'):
    """
    Calculate the SFR from input nebular line emission

    Args:
        neb_lines (dict):  Observed line fluxes
        method (str): Method for deriving the SFR
          Ha -- Use the Halpha line flux
        z (float):  Emission redshift -- for Luminosity distance
        cosmo (astropy.cosmology.FLRW): Cosmology
        AV (float, optional):  Visual extinction, if supplied will apply
        curve (str):  Name of the extinction curve.  Only used if A_V is supplied

    Returns:
        Quantity:  SFR with units of Msun/yr

    """
    if method == 'Ha':
        line = 'Halpha'
        conversion = Ha_conversion
    elif method == 'Hb':
        line = 'Hbeta'
        conversion = Ha_conversion * Ha_Hb_intrin
    else:
        raise IOError("Not prepared for method: {}".format(method))

    # Luminosity
    Lum, Lum_err = calc_lum(neb_lines, line, z, cosmo, AV=AV, curve=curve)

    # SFR
    SFR = Lum.to('erg/s').value * conversion

    return SFR

def get_ebv(coords,definition="SFD",region=5*units.deg):
    """
    Get the E(B-V) value and statistic from the Milky way dust extinction
    within the query region around the input coordinate
    
    Args:
        coords: Astropy SkyCoord
            Input celestial coordinates
        definition: str, optional
            Can be either "SFD" or "SandF". They stand for the 
            definitions of E(B-V) according to either Schlegel et al. 1998 (ApJ 500, 525)
            or Schlafly and Finkbeiner 2011 (ApJ 737, 103) respectively
        region: Astropy Angle (Quantity), optional
            Angular radius around the input coordinate where
            the query is run to obtain statistics. Must be between
            2 deg and 37.5 deg. Default value: 5 deg.
        get_stats: bool, optional
            If true, also returns a dict with the statistics of E(B-V)
            within the query region.
    Returns:
        ebvdict: dict
            Dict with E(B-V) at pixel, mean, std, min and max values in
            the query region. All values are in mags.
    """
    import pdb
    assert definition in ['SFD','SandF'], "definition can only be one of 'SFD' and 'SandF'"
    assert (region>2*units.deg) & (region<37.5*units.deg), "Search radius must be between 3 and 37.5 degrees"

    ra,dec = str(coords.ra.value),str(coords.dec.value)
    radius = str(region.to(units.deg).value)
    query_url = \
        "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={:s}+{:s}+equ+J2000&regSize={:s}".format(ra,dec,radius)
    result = requests.get(query_url)
    #pdb.set_trace()
    tree = ET.ElementTree(ET.fromstring(result.content.decode("ascii")))
    root = tree.getroot()

    #By default it looks at the first entry with the name "results"
    #This corresponds to the reddening data.
    #TODO: make this smarter. It should look for the reddening results
    #by name.
    statchild = root.find('result/statistics')

    ebvdict = {}
    for elem in statchild.findall('*'):
        if definition in elem.tag:
            ebvdict[elem.tag.replace(definition,'')] = elem.text.split()[0]
    return ebvdict
