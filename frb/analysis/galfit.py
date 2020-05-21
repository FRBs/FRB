"""
Automate the GALFIT analysis procedure. GALFIT
is a free, closed-source software that can be found
here:
https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html

Here is the user manual:
https://users.obs.carnegiescience.edu/peng/work/galfit/README.pdf

This module is mainly for a python wrapper to automate
using GALFIT. It also contains additional functions
to parse the output automatically.
"""

import numpy as np
import pandas as pd
import sys, os, subprocess
import warnings
from astropy.io import fits
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS

def genconf(imgfile:str, psffile:str,
            configfile:str=None, outfile:str=None,
            finesample:int = 1, badpix:str = "none",
            constraints:str = "none",
            region:tuple = None, convobox:tuple = (100,100),
            zeropoint:float = 25.0, platescale:float = 0.125,
            position:tuple = None, int_mag:float = None,
            r_e:float = None, n:float = 1.0, b/a:float = 0.5,
            pa:float = 0, skip_sky:bool = False):
    """
    Creates a configuration file for GALFIT. Conventionally,
    GALFIT is run using the command: `galfit <config-file>`.
    Args
    ----
    imgfile (str): path to the image fits file.
    psffile (str): path to the PSF model fits file.
    configfile (str, optional): path to configuration file to
        be created via this function.
    outfile (str, optional): path to GALFIT's output fits file.
        Defaults to `out.fits`
    finesample (int, optional): The PSF fine-sampling factor.
        Assumes no fine-sampling (i.e. a factor of 1) by default.
    badpix (str, optional): File containing a list of bad pixels.
        Assumes all pixels a re fine by default.
    constraints (str, optional): File containing fit parameter
        bounds. Check out this example constraints file
        https://users.obs.carnegiescience.edu/peng/work/galfit/EXAMPLE.CONSTRAINTS
        to learn how to use one.
    region (tuple, optional): Pixel coordinate bounds for the fitting.
        Required format: (xmin, xmax, ymin, ymax). All of the
        elements must be ints!
    convobox (tuple, optional): Size of the convolution box
        used to assess the model fit chi-squared value.
        Required format: (x_box:int , y_box:int)
    zeropoint (float, optional): Zeropoint of the image
        to compute the model magnitude.
    platescale (float, optional): Size of pixel in arcsec.
        Assumes 0.125''/pixel by default.
    position (tuple, optional): Guess for the centroid
        of the model fit. Required format (x_cen:float, y_cen:float).
        Assumes the center of the image is the initial guess by
        default.
    int_mag (float, optional): Initial guess for the 
        integrated magnitude. Taken as the magnitude
        corresponding to the sum of all counts in the
        region to be fit.
    r_e (float, optional): Initial guess for
        the half light radius in pixels. Assumes
        half the size of the fitting region by default.
    n (float, optional): initial guess for the Sersic
        index.
    b/a (float, optional): Initial guess for
        the ratio of minor to major axis of the fit model.
    PA (float, optional): Initial guess for the 
        angle of the major axis counter clockwise 
        relative to the vertical.
    skip_sky (bool, optional): Do you also want
        to fit a constant sky background? Set to 
        false if your sky background is 0.
    Returns
    -------
    configfile (str): Path to the configuration
        file.
    """
    if configfile is None:
        warnings.warn("Creating a configuration file here")
        configfile = "galfit.feedme"
    with open(configfile,"w+") as fstream:
        #Image parameters.
        fstream.write("""===============================================================================\n
        # IMAGE and GALFIT CONTROL PARAMETERS\n
        """)
        assert os.path.isfile(imgfile), "Invalid image file path {:s}".format(imgfile)
        fstream.write("A) {:s}  # Input data image (FITS file)\n".format(imgfile))
        if outfile is None:
            warnings.warn("Creating output file here")
            outfile = "out.fits"
        fstream.write("B) {:s}  # Output data image block\n".format(outfile))
        fstream.write("C) none  # Sigma image name (made from data if blank or 'none')\n")
        assert os.path.isfile(psffile), "Invalid psf file path {:s}".format(psffile)
        fstream.write("D) {:s}  # Input PSF file\n".format(psffile))
        fstream.write("E) {:d}  # PSF fine sampling factor\n".format(finesample))
        fstream.write("F) {:s}  #Bad pixel mask\n".format(badpix))
        fstream.write("G) {:s}  # File with parameter constraints (ASCII file)\n".format(constraints))
        if region is None:
            input_str = input("Input image region to fit (xmin xmax ymin ymax): ").split()
            while len(input_str)!=4:
                input_str = input("Invalid region input. Insert again (xmin xmax ymin ymax): ").split()
            region = [int(bound) for bound in input_str]
        fstream.write("H) {:d} {:d} {:d} {:d}   # Image region to fit (xmin xmax ymin ymax)\n".format(region[0],region[1],region[2], region[3]))
        assert len(convobox)==2, "Only two integers specify convolution box dimensions"
        fstream.write("I) {:d} {:d} # Size of convolution box (x y)\n".format(convobox[0],convobox[1]))
        fstream.write("J) {:f}  # Photometric zeropoint (mag)\n".format(zeropoint))
        fstream.write("K) {:f} {:f} # Plate scale (dx dy) [arcsec/pixel]\n".format(platescale,platescale))
        fstream.write("O) regular   # Display type (regular, curses, both)\n")
        fstream.write("P) 0 # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")
        fstream.write("\n# INITIAL FITTING PARAMETERS\n")
        fstream.write("# Component number: 1\n")
        fstream.write("0) sersic # Component type\n")
        img = fits.getdata(imgfile)
        if position is None:
            pos_x, pos_y = img.shape/2 
        else:
            pos_x, pos_y = position
        fstream.write("1) {:d} {:d} 1 1 # position x y\n".format(pos_x,pos_y))
        if int_mag is None:
            int_mag = zeropoint-2.5*np.log10(np.sum(img[region[2]:region[3],region[0]:region[1]]))
        fstream.write("3) {:f} 1 # Integrated magnitude\n".format(int_mag))
        if r_e is None:
            r_e = (region[3]-region[2]+region[1]-region[0])/2/3 #A third of the average region dimension
        fstream.write("4) {:f} 1 # effective radius (pix)\n".format(r_e))
        fstream.write("5) {:f} 1 # sersic index\n".format(n))
        for num in [6,7,8]:
            fstream.write("{:d}) 0.0000 0 # ----\n".format(num))
        fstream.write("9) {:f} 1 # Axis ratio (b/a)\n".format(b/a))
        fstream.write("10) {:f} 1 # Position angle (PA) [deg: Up=0, left =90]\n".format(pa))
        fstream.write("Z) 0 # Skip this model? (yes=1,no=0)\n\n")
        if not skip_sky:
            fstream.write("# Component number: 2\n")
            fstream.write("0) sky # component type\n")
            _, median, _ = sigma_clipped_stats(img)
            fstream.write("1) {:f} 1 # Sky background\n".format(median))
            fstream.write("2) 0 0 # dsky/dx\n")
            fstream.write("3) 0 0 # dsky/dy\n")
            fstream.write("Z) 0 # Skip this model\n")
        fstream.write("================================================================================\n")
    return configfile

def run(imgfile, psffile,**kwargs):
    """
    Run galfit. The arguments are
    the same as genconf.
    """
    configfile = genconf(imgfile,psffile,**kwargs)
    os.system("galfit {:s}".format(configfile))
    return 0

def pix2coord(sersic_tab, wcs, platescale):
    """
    Takes the output table from galfit's
    fit.log file and converts all pixel
    measurements to physical measurements.
    Args
    ----
    sersic_tab (Table): Raw sersic
        fit parameter table
    wcs (WCS): transformation to
        go from pixels to angular units
    platescale (Quantity): angular size
        per pixel.
    Returns
    -------
    sky_tab (Table): Sersic table
        translated to angular units
        on the sky.
    """

    sky_tab = Table()
    platescale = platescale.to(u.arcsec)

    # First read the invariants
    sky_tab['b/a'] = sersic_tab['b/a']
    sky_tab['b/a_err'] = sersic_tab['b/a_err']
    sky_tab['n'] = sersic_tab['n']
    sky_tab['n_err'] = sersic_tab['n_err']
    sky_tab['PA'] = sersic_tab['PA']
    sky_tab['PA_err'] = sersic_tab['PA_err']
    sky_tab['mag'] = sersic_tab['mag']
    sky_tab['mag_err'] = sersic_tab['mag_err']

    # Read in centroids next
    xpix, ypix = sersic_tab['x']-1, sersic_tab['y']-1
    centr_coords = wcs.pixel_to_world(xpix, ypix)
    sky_tab['ra'] = centr_coords.ra.value
    sky_tab['dec'] = centr_coords.dec.value

    sky_tab['ra_err'] = sersic_tab['x_err']*platescale
    sky_tab['dec_err'] = sersic_tab['y_err']*platescale

    # Read in r_e next.
    sky_tab['reff_ang'] = sersic_tab['reff_ang']*platescale
    sky_tab['reff_ang_err'] = sersic_tab['reff_ang_err']*platescale

    # Reorder columns
    sky_tab = sky_tab['ra','ra_err', 'dec','dec_err', 'mag','mag_err','reff_ang', 'reff_ang_err', 'n', 'n_err','b/a','b/a_err','PA','PA_err']

    return sky_tab

def surf_brightness(coord, sky_tab):
    """
    Estimates the surface brightness from
    the sersic model fit at the input
    coordinate location.
    Args
    ----
    coord (SkyCoord): Target coordinate
    sky_tab (Table): Table of best fit
        sersic model parameters in
        angular coordinates.
    Returns
    ------
    surf_brightness (m)
    """
