"""
Module to generate CIGALE .ini file 
"""

import numpy as _np
from astropy.table import Table as _Table
from frb.galaxies import defs

    

def gen_cigale_ini(frbgal, ID=None,filename='data.fits'):
    """
    Generates the input data file for CIGALE
    given the photometric points and redshift
    of a galaxy
    Parameters
    ----------
    frbgal: FRBGalaxy or FRBHost object
        A Galaxy object having redshift estimates
        and photometric data available.
    ID: str, optional
        An ID for the galaxy. If none, "GalaxyA" is assigned.
    filename: str, optional
        Name of fits file (with path if needed) to store data in.
        Default value is 'data.fits'
    """
    assert (frbgal.photom != {}),"No photometry found. CIGALE cannot be run."
    assert (frbgal.redshift != {}),"No redshift found. CIGALE cannot be run"

    photom = frbgal.photom
    if ID is None:
        ID = "GalaxyA"
    photom['ID'] = ID
    photom['redshift'] = frbgal.redshift
    
    #Convert DES fluxes to mJy
    for band in defs.DES_bands:
        colname = "DES_"+band
        photom[colname] = 3630780.5*10**(photom[colname]/-2.5)
        photom[colname+"_err"] = photom[colname+"_err"]/1.087*photom[colname]
    
    #Convert WISE fluxes to mJy
    wise_fnu0 = [309.54,171.787,31.674,8.363] #http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2flux
    for band,zpt in zip(defs.WISE_bands,wise_fnu0):
        photom[band] = zpt*10**(-photom[band]/2.5)
        errname = band+"_err"
        if photom[errname]!=-999.0:
            photom[errname] =-99.0
        else:
            photom[errname] = photom[errname]/1.087*photom[band]
    
    #Write to file
    photom = _Table([photom])
    photom.write(filename,format="fits")