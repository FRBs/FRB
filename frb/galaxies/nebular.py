""" Methods related to nebular line analysis, e.g. dust extinction, SFR"""

import numpy as np
import requests
import warnings

from xml.etree import ElementTree as ET

from astropy.table import Table
from astropy import units

try:
    import dust_extinction
except ImportError:
    warnings.warn("Galaxy nebular line analysis requires dust_extionction.  Install it if you want to use them")


try:
    from linetools.lists import linelist
except ImportError:
    warnings.warn("Galaxy nebular line analysis requires linetools.  Install it if you want to use them")

# GLOBALS
Ha_Hb_intrin = 2.87  # Osterbrock 2006 Book
Hb_Hg_intrin = 1./0.466  # Osterbrock 2006 Book
salpeter_to_chabrier_SFR = 0.63
Oii_conversion = salpeter_to_chabrier_SFR * 1.4e-41 * units.Msun/units.yr # Kennicutt 1998
Oii_conversion_err = salpeter_to_chabrier_SFR * 0.4e-41 * units.Msun/units.yr
Ha_conversion = salpeter_to_chabrier_SFR * 7.9e-42 * units.Msun/units.yr   # Kennicutt 1998 + Chabrier,
# e.g. https://ned.ipac.caltech.edu/level5/March14/Madau/Madau3.html


def calc_dust_extinct(neb_lines, method):
    """
    Estimate the Visual extinction A_V based on input nebular emission lines

    Uses the Gordon 2024


    Args:
        neb_lines (dict):  Line fluxes and errors (keys: line name, line name + '_err')
        method (str): Name of the method
          Ha/Hb -- Use the Halpha/Hbeta ratio and standard intrinsic flux
        curve (str): Extinction curve to use


    Returns:
        float: A_V in magnitudes
        float or None: 1-sigma uncertainty on A_V (None if line errors unavailable)

    """
    if method == 'Ha/Hb':
        wave1 = 6564.6  # redder
        wave2 = 4862.7
        #
        F1_obs = neb_lines['Halpha']
        F2_obs = neb_lines['Hbeta']
        F1_err = neb_lines.get('Halpha_err', -1.)
        F2_err = neb_lines.get('Hbeta_err', -1.)
        #
        pair = True
        intrinsic = Ha_Hb_intrin
    elif method == 'Hb/Hg':
        wave1 = 4862.7 #redder
        wave2 = 4341.7
        #
        F1_obs = neb_lines['Hbeta']
        F2_obs = neb_lines['Hgamma']
        F1_err = neb_lines.get('Hbeta_err', -1.)
        F2_err = neb_lines.get('Hgamma_err', -1.)
        #
        pair = True
        intrinsic = Hb_Hg_intrin
    else:
        print("Not prepared for this method of analysis: {}".format(method))
        raise IOError("See docs for available options")

    if not pair:
        raise IOError("Not ready for this mode")

    # Extinction
    #a1AV = extinction.fm07(np.atleast_1d(wave1), 1.0)[0]
    #a2AV = extinction.fm07(np.atleast_1d(wave2), 1.0)[0]
    extmod = dust_extinction.parameter_averages.G23(Rv=3.1)
    a1AV = extmod(np.atleast_1d(wave1*units.AA))[0]  # *units.AA)
    a2AV = extmod(np.atleast_1d(wave2*units.AA))[0]  # *units.AA)

    # Observed ratio
    fratio_obs = F1_obs/F2_obs

    # Calculate using intrinsic ratio
    AV = 2.5 * np.log10(intrinsic/fratio_obs) / (a1AV - a2AV)

    # Propagate line flux errors through Balmer decrement
    # AV = 2.5*log10(intrinsic*F2/F1) / (a1AV-a2AV)
    # dAV/dF1 = -2.5 / (ln10 * F1 * (a1AV-a2AV))
    # dAV/dF2 = +2.5 / (ln10 * F2 * (a1AV-a2AV))
    if F1_err > 0. and F2_err > 0.:
        denom = np.log(10) * (a1AV - a2AV)
        dAV_dF1 = -2.5 / (denom * F1_obs)
        dAV_dF2 = +2.5 / (denom * F2_obs)
        AV_err = np.sqrt((dAV_dF1 * F1_err)**2 + (dAV_dF2 * F2_err)**2)
    else:
        AV_err = None

    # Return
    return AV, AV_err


def calc_logOH(neb_lines, method):
    """ 
    Estimate the oxygen abundance based on the input nebular emission line fluxes
    For now based on the O3N2 calibration from https://ui.adsabs.harvard.edu/abs/2018AJ....155...82H/abstract


    Args:
        neb_lines (dict):  Line fluxes
        method (str): Name of the method
          O3N2 -- Use the O3N2 calibration from Hirschauer+18


    Returns:
        tuple: 12+log(O/H), sigma+, sigma-
    """

    if method == 'O3N2':
        # Check for all lines
        req_lines = ['[NII] 6584','Halpha','[OIII] 5007','Hbeta']
        for iline in req_lines:
            if iline not in neb_lines.keys():
                print("One or more lines missing for logOH calculation.  Returning None's")
                return None, None, None
        # Proceed.  pPXF reports the total doublet flux for [NII] and [OIII];
        # the 6584/5007 components are 75% of the total (3:1 ratio within each doublet).
        x0 = neb_lines['[NII] 6584']  * 0.75 / neb_lines['Halpha']
        y0 = neb_lines['[OIII] 5007'] * 0.75 / neb_lines['Hbeta']
        x0_err = x0 * np.sqrt((neb_lines['[NII] 6584_err'] / neb_lines['[NII] 6584'])**2 + (neb_lines['Halpha_err'] / neb_lines['Halpha'])**2)
        y0_err = y0 * np.sqrt((neb_lines['[OIII] 5007_err'] / neb_lines['[OIII] 5007'])**2 + (neb_lines['Hbeta_err'] / neb_lines['Hbeta'])**2)
    else:
        raise IOError("Not ready for this method for logOH")
        
    # Calculate O3N2 (linear)
    o3n2 = y0 / x0
    o3n2_err = o3n2 * np.sqrt((x0_err/x0)**2 + (y0_err/y0)**2)

    # Log it
    log_o3n2 = np.log10(o3n2)
    log_03n2_errp = np.log10(o3n2 + o3n2_err) - log_o3n2
    log_03n2_errm = log_o3n2 - np.log10(o3n2 - o3n2_err)

    # Hirschauer+18 O3N2 calibration
    logOH = 8.987 - 0.297*log_o3n2 - 0.0592*(log_o3n2)**2 - 0.009*(log_o3n2)**3
    logOH_errp = log_03n2_errp
    logOH_errm = log_03n2_errm      
 
    # Return 
    return logOH, logOH_errp, logOH_errm


def calc_lum(neb_lines, line, z, cosmo, AV=None, AV_err=None):
    """
    Calculate the line luminosity (and error) from input nebular line emission

    Error is -999.*erg/s if input line flux has negative error


    Args:
        neb_lines (dict):  Observed line fluxes and errors
        line (str): Line to analyze
        z (float):  Emission redshift -- for Luminosity distance
        cosmo (astropy.cosmology.FLRW): Cosmology
        AV (float, optional):  Visual extinction, if supplied will apply
        AV_err (float, optional):  1-sigma uncertainty on AV; propagated into Lum_err


    Returns:
        Quantity, Quantity:  Luminosity, sig(Luminosity)
    """
    # Grab rest wavelength (vacuum)
    llist = linelist.LineList('Galaxy')
    wave = llist[line]['wrest']

    # Dust correct?
    AlAV = 0.
    if AV is not None:
        #al = extinction.fm07(np.atleast_1d(wave.to('Angstrom').value), AV)[0]
        extmod = dust_extinction.parameter_averages.G23(Rv=3.1)
        AlAV = float(extmod(np.atleast_1d(wave))[0])
        al = AlAV * AV
    else:
        al = 0.

    # Cosmology
    DL = cosmo.luminosity_distance(z)

    # Luminosity
    flux = neb_lines[line]
    Lum = flux * units.erg/units.s/units.cm**2 * 10**(al/2.5) * (4*np.pi * DL**2)

    # Error: flux contribution
    if neb_lines[line+'_err'] > 0.:
        flux_err = neb_lines[line+'_err']
        Lum_err = flux_err * units.erg/units.s/units.cm**2 * 10**(al/2.5) * (4*np.pi * DL**2)
        # AV uncertainty contribution: dLum/dAV = Lum * ln(10)/2.5 * AlAV
        if AV is not None and AV_err is not None:
            Lum_err_AV = (abs(Lum.to('erg/s').value) * np.log(10) / 2.5 * AlAV * AV_err
                          ) * units.erg/units.s
            Lum_err = np.sqrt(Lum_err.to('erg/s')**2 + Lum_err_AV**2)
    else:
        Lum_err = -999 * units.erg/units.s
    # Return
    return Lum.to('erg/s'), Lum_err.to('erg/s')


def calc_SFR(neb_lines, method, z, cosmo, AV=None, AV_err=None, curve='MW'):
    """
    Calculate the SFR from input nebular line emission


    Args:
        neb_lines (dict):  Observed line fluxes
        method (str): Method for deriving the SFR
          Ha -- Use the Halpha line flux
        z (float):  Emission redshift -- for Luminosity distance
        cosmo (astropy.cosmology.FLRW): Cosmology
        AV (float, optional):  Visual extinction, if supplied will apply
        AV_err (float, optional):  1-sigma uncertainty on AV; propagated into SFR_err
        curve (str):  Name of the extinction curve.  Only used if A_V is supplied


    Returns:
        Quantity, Quantity:  SFR, SFR_err with units of Msun/yr

    """
    if method == 'Ha':
        line = 'Halpha'
        conversion = Ha_conversion
    elif method == 'Hb':
        line = 'Hbeta'
        conversion = Ha_conversion * Ha_Hb_intrin
    elif method == 'Oii':
        line = '[OII] 3727'
        conversion = Oii_conversion
    else:
        raise IOError("Not prepared for method: {}".format(method))

    # Luminosity
    Lum, Lum_err = calc_lum(neb_lines, line, z, cosmo, AV=AV, AV_err=AV_err)

    # SFR
    SFR = Lum.to('erg/s').value * conversion
    SFR_err = Lum_err.to('erg/s').value * conversion
    # [OII] calibration uncertainty: σ_conv contribution is L * σ_k (in quadrature)
    if method == 'Oii' and Lum_err.to('erg/s').value > 0:
        sfr_conv_err = Lum.to('erg/s').value * Oii_conversion_err
        SFR_err = np.sqrt(SFR_err**2 + sfr_conv_err**2)

    return SFR, SFR_err

def get_ebv(coords,definition="SandF",
            region=5*units.deg,get_ext_table=False):
    """
    Get the E(B-V) value and statistic from the Milky way dust extinction
    within the query region around the input coordinate

    
    Args:
        coords (Astropy SkyCoord):
            Input celestial coordinates
        definition (str, optional):
            Can be either "SFD" or "SandF". They stand for the 
            definitions of E(B-V) according to either Schlegel et al. 1998 (ApJ 500, 525)
            or Schlafly and Finkbeiner 2011 (ApJ 737, 103) respectively
        region (Astropy Angle (Quantity), optional):
            Angular radius around the input coordinate where
            the query is run to obtain statistics. Must be between
            2 deg and 37.5 deg. Default value: 5 deg.
        get_stats: bool, optional
            If true, also returns a dict with the statistics of E(B-V)
            within the query region.
        get_ext_table: bool, optional
            If true, also returns the table with A/E(B-V) ratios
            for multiple filters.

    Returns:
        dict:
            Dict with E(B-V) at refPixelValue, meanValue, std, minValue and maxValue in
            the query region. All values are in mags.
    """
    assert definition in ['SFD','SandF'], "definition can only be one of 'SFD' and 'SandF'"
    assert (region>2*units.deg) & (region<37.5*units.deg), "Search radius must be between 3 and 37.5 degrees"

    # Coords
    ra,dec = str(coords.ra.value),str(coords.dec.value)
    radius = str(region.to(units.deg).value)
    query_url = \
        "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={:s}+{:s}+equ+J2000&regSize={:s}".format(ra,dec,radius)
    result = requests.get(query_url)
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
            ebvdict[elem.tag.replace(definition,'')] = float(elem.text.split()[0])
    
    if get_ext_table:
        table_url = root.find('result/data/table').text.split()[0]
        ext_table = Table.read(table_url,format="ascii")
        return ebvdict,ext_table
    
    return ebvdict
