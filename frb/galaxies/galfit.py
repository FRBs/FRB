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
import re

def write_cutout(cutout, filename = "cutout.fits", overwrite=False):
    """
    Takes an astropy Cutout2D
    object and writes it to disk.
    Args:
        cutout (Cutout2D): Prefaerably
            has WCS info in it.
        filename (str, optional): FITS
            file into which data is saved.
            Defaults to "cutout.fits" in
            the current working directory.
        overwrite (bool, optional): Do you
            want to overwrite filename if
            it already exists?
    Returns:
    """
    hdr = cutout.wcs.to_header()
    imghdu = fits.PrimaryHDU(cutout.data, hdr)
    hdulist = fits.HDUList([imghdu])
    hdulist.writeto(filename, overwrite=overwrite)
    return

def _genconf(imgfile:str, psffile:str,
            configfile:str=None, outdir:str=None, outfile:str=None,
            finesample:int = 1, badpix:str = "none",
            constraints:str = "none",
            region:tuple = None, convobox:tuple = (100,100),
            zeropoint:float = 25.0, platescale:float = 0.125,
            position:tuple = None, int_mag:float = None,
            r_e:float = None, n:float = 1.0, axis_ratio:float = 0.5,
            pa:float = 0, skip_sky:bool = False):
    """
    Creates a configuration file for GALFIT. Conventionally,
    GALFIT is run using the command: `galfit <config-file>`.
    Args:
        imgfile (str): path to the image fits file.
        psffile (str): path to the PSF model fits file.
        outdir (str): Name of output directory. Default
            value is 'galfit_out` in the current directory.
        configfile (str, optional): path to configuration file to
            be created via this function. Defaults to `galfit.feedme`
            in <outdir>.
        outfile (str, optional): name of GALFIT's output fits file.
            Defaults to `out.fits` in <outdir>.
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
        axis_ratio (float, optional): Initial guess for
            the ratio of minor to major axis of the fit model.
        PA (float, optional): Initial guess for the 
            angle of the major axis counter clockwise 
            relative to the vertical.
        skip_sky (bool, optional): Do you also want
            to fit a constant sky background? Set to 
            false if your sky background is 0.
    Returns:
        configfile (str): Path to the configuration
            file.
    """
    if outdir is None:
        outdir = "galfit_out"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if configfile is None:
        warnings.warn("Creating a configuration file here")
        configfile = "galfit.feedme"
    
    with open(os.path.join(outdir,configfile),"w+") as fstream:
        #Image parameters.
        fstream.write("""===============================================================================\n
        # IMAGE and GALFIT CONTROL PARAMETERS\n
        """)
        img = fits.getdata(imgfile)
        fstream.write("A) {:s}  # Input data image (FITS file)\n".format(imgfile))
        if outfile is None:
            warnings.warn("Creating output file here")
            outfile = "out.fits"
        fstream.write("B) {:s}  # Output data image block\n".format(outfile))
        fstream.write("C) none  # Sigma image name (made from data if blank or 'none')\n")
        fstream.write("D) {:s}  # Input PSF file\n".format(psffile))
        fstream.write("E) {:d}  # PSF fine sampling factor\n".format(finesample))
        fstream.write("F) {:s}  #Bad pixel mask\n".format(badpix))
        fstream.write("G) {:s}  # File with parameter constraints (ASCII file)\n".format(constraints))
        if region is None:
            input_str = input("Input image region to fit (xmin xmax ymin ymax) or use 'all' ").split()
            #import pdb; pdb.set_trace()
            if input_str == ['all']:
                region = (0, img.shape[1]-1, 0, img.shape[0]-1)
            elif len(input_str)==4:
                try:
                    region = [int(bound) for bound in input_str]
                except:
                    raise RuntimeError("Invalid region input")
            else:
                raise RuntimeError("Invalid region input")
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
        if position is None:
            pos_x, pos_y = int(img.shape[0]/2), int(img.shape[1]/2)  
        else:
            pos_x, pos_y = position
        fstream.write("1) {:d} {:d} 1 1 # position x y\n".format(pos_x,pos_y))

        # What is the integrated magnitude?
        if int_mag is None:
            warnings.warn("No guess given for integrated magnitude. Proceeding with sum within region.")
            int_mag = zeropoint-2.5*np.log10(np.sum(img[region[2]:region[3],region[0]:region[1]]))
        fstream.write("3) {:f} 1 # Integrated magnitude\n".format(int_mag))
        if r_e is None:
            warnings.warn("Guess for r_e not given. This might not converge.")
            r_e = (region[3]-region[2]+region[1]-region[0])/2/3 #A third of the average region dimension
        fstream.write("4) {:f} 1 # effective radius (pix)\n".format(r_e))
        fstream.write("5) {:f} 1 # sersic index\n".format(n))
        for num in [6,7,8]:
            fstream.write("{:d}) 0.0000 0 # ----\n".format(num))
        fstream.write("9) {:f} 1 # Axis ratio (b/a)\n".format(axis_ratio))
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

def read_fitlog(outfile, initfile):
    """
    Reads the output fit.log
    file and returns the sersic fit parameters
    """
    # TODO: create a regex expression to look for blocks
    # in the fit.log file. Then create an expression
    # to lock for model sub-blocks within each block. 
    lines = [line.rstrip('\n') for line in open(outfile)]
    instance = []  # going to put all instances of use of input file here
    for kk, line in enumerate(lines):  # for all lines in fit.log
        if 'Init. par. file : ' + str(initfile) in line:  # if the it is from the input file used
            for pp, item in enumerate(lines[kk:kk + 10]):
                if 'sersic' in item:  # look for the sersic fit model
                    # This assumes each fitting had a single
                    # sersic model fit.
                    # TODO: accommodate more complex fits.
                    instance.append((item, pp, lines[kk:kk + 10][pp + 1]))  # and keep the model and error

    fit_out, line, err_out = instance[-1]  # only keep the latest one

    # Use regex to identify all floats in fit_out and err_out
    float_regex = r"\d+\.\d+" # Any string of the form <int>.<int> (which is just a float)

    # Find all instances of float in the strings and convert from string to float
    fitparams = np.array(re.findall(float_regex, fit_out)).astype('float')
    fiterrs = np.array(re.findall(float_regex, err_out)).astype('float')

    #Store values
    pix_dict = {'x':fitparams[0], 'y':fitparams[1],
              'mag':fitparams[2], 'reff':fitparams[3],
              'n':fitparams[4], 'b/a':fitparams[5], 'PA':fitparams[6],
              'x_err':fiterrs[0], 'y_err':fiterrs[1],
              'mag_err':fiterrs[2], 'reff_err':fiterrs[3],
              'n_err':fiterrs[4], 'b/a_err':fiterrs[5], 'PA_err':fiterrs[6],
              }

    return pix_dict

def pix2coord(pix_dict, wcs, platescale, table=False):
    """
    Takes the output table from galfit's
    fit.log file and converts all pixel
    measurements to physical measurements.
    Args:
        pix_dict (dict): Raw sersic
            fit parameter dict from fit.log
        wcs (WCS): Input image WCS.
        platescale (float): angular size
            per pixel in arcsec.
        table (bool, optional): Return as
            a table instead?
    Returns:
        sky_dict (dict/Table): pix_dict
            translated to angular units
            on the sky.
    """

    sky_dict = {}

    # Centroids
    xpix, ypix = pix_dict['x']-1, pix_dict['y']-1
    centr_coords = wcs.pixel_to_world(xpix, ypix)
    sky_dict['ra'] = centr_coords.ra.value
    sky_dict['dec'] = centr_coords.dec.value
    sky_dict['ra_err'] = pix_dict['x_err']*platescale
    sky_dict['dec_err'] = pix_dict['y_err']*platescale
    # Magnitude
    sky_dict['mag'] = pix_dict['mag']
    sky_dict['mag_err'] = pix_dict['mag_err']
    # Half-light radius
    sky_dict['reff_ang'] = pix_dict['reff']*platescale
    sky_dict['reff_ang_err'] = pix_dict['reff_err']*platescale
    # Sersic index
    sky_dict['n'] = pix_dict['n']
    sky_dict['n_err'] = pix_dict['n_err']
    # Axis ratio
    sky_dict['b/a'] = pix_dict['b/a']
    sky_dict['b/a_err'] = pix_dict['b/a_err']
    # Sky position angle
    sky_dict['PA'] = pix_dict['PA']
    sky_dict['PA_err'] = pix_dict['PA_err']

    if table:
        sky_dict = Table([sky_dict])
        # Reorder because dict-> Table messes it up
        sky_dict = sky_dict['ra','ra_err','dec','dec_err','mag',
                            'mag_err','reff_ang','reff_ang_err',
                            'n','n_err','b/a','b/a_err','PA','PA_err']
    return sky_dict

def run(imgfile, psffile, platescale=0.125, **kwargs):
    """
    Run galfit. 
    
    Args:
        imgfile (str): path to the image fits file.
        psffile (str): path to the PSF model fits file.
        platescale (float, optional): Size of pixel in arcsec.
            Assumes 0.125''/pixel by default.

    Valid kwargs:
        outdir (str): Name of output directory. Default
            value is 'galfit_out` in the current directory.
        configfile (str, optional): path to configuration file to
            be created via this function. Defaults to `<outdir>/galfit.feedme`.
        outfile (str, optional): path to GALFIT's output fits file.
            Defaults to `<outdir>/out.fits`
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
        axis_ratio (float, optional): Initial guess for
            the ratio of minor to major axis of the fit model.
        PA (float, optional): Initial guess for the 
            angle of the major axis counter clockwise 
            relative to the vertical.
        skip_sky (bool, optional): Do you also want
            to fit a constant sky background? Set to 
            false if your sky background is 0.
    Returns:
        configfile (str): Path to the configuration
            file.
    """

    # Check input paths first
    assert os.path.isfile(imgfile), "Invalid image file path {:s}".format(imgfile)
    assert os.path.isfile(psffile), "Invalid psf file path {:s}".format(psffile)

    # Is a badpix map given?
    if 'badpix' in kwargs.keys():
        assert os.path.isfile(kwargs['badpix']), "Invalid bad pixel map path {:s}".format(kwargs['badpix'])
        # Use absolute path henceforth
        kwargs['badpix'] = os.path.abspath(kwargs['badpix'])
    
    # Is a constraints file given?
    if 'constraints' in kwargs.keys():
        assert os.path.isfile(kwargs['constraints']), "Invalid constraints file path {:s}".format(kwargs['constraints'])
        # Use absolute path henceforth
        kwargs['constraints'] = os.path.abspath(kwargs['constraints'])

    # Abspaths for image and psf files too.
    imgfile = os.path.abspath(imgfile)
    psffile = os.path.abspath(psffile)
    
    # Generate Galfit config file
    configfile = _genconf(imgfile, psffile, platescale=platescale,**kwargs)

    # Go to the output directory
    if 'outdir' not in kwargs:
        kwargs['outdir'] = "galfit_out"
    # Keep track of the original directory
    curdir = os.path.abspath(os.path.curdir)
    # Move to outdir
    os.chdir(kwargs['outdir'])
    # Run galfit
    return_value = os.system("galfit {:s}".format(configfile))

    if return_value!=0:
        # Something broken?
        warnings.warn("Something went wrong with the fit. Check terminal output.")
        os.chdir(curdir)
        return return_value
    # Read fit.log and get the fit results
    pix_dict = read_fitlog("fit.log", configfile)
    # Convert to sky angular units and stuff it into the output
    # fits file
    hdr = fits.getheader(imgfile)
    wcs = WCS(hdr)
    sky_tab = pix2coord(pix_dict,wcs, platescale, table=True)
    if 'outfile' not in kwargs:
        kwargs['outfile'] = 'out.fits'
    fitloghdu = fits.BinTableHDU(sky_tab,name="FITPARAMS")
    hdulist = fits.open(kwargs['outfile'])
    # This needs to be done for some blackbox
    # in astropy to not barf.
    for idx in [2,3]:
        hdulist[idx].header.insert('OBJECT',('PCOUNT',0))
        hdulist[idx].header.insert('OBJECT',('GCOUNT',1))

    # Dump the table in the fits file
    hdulist.append(fitloghdu)
    # Overwrite
    hdulist.writeto(kwargs['outfile'], overwrite=True)
    os.chdir(curdir)
    return return_value

def surf_brightness(coord, sky_dict):
    """
    Estimates the surface brightness from
    the sersic model fit at the input
    coordinate location.
    Args
    ----
    coord (SkyCoord): Target coordinate
    sky_dict (Table): Table of best fit
        sersic model parameters in
        angular coordinates.
    Returns
    ------
    surf_brightness (mag/arcsec**2):
    """
    
