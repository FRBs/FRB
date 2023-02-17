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
import sys, os, subprocess
import warnings

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import Cutout2D

import re

def get_platescale(wcs:WCS)->float:
    """
    Extract the plate-scale from
    a WCS object
    Args:
        wcs (WCS): A celestial WCS object
            extracted from a FITS header.
    Returns:
        platescale (float): arcsec/pixel
    """
    assert wcs.is_celestial or wcs.has_celestial, "Can't get a plate scale from a non-celestial WCS"

    if not wcs.is_celestial:
        wcs = wcs.celestial
    
    platescale = np.mean(np.sum(wcs.pixel_scale_matrix**2, axis=0)**0.5)*3600 #arcsec
    return platescale

def write_cutout(cutout:Cutout2D, filename:str = "cutout.fits", overwrite:bool=False, exptime:float=None):
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
        exptime (float, optional): exposure
            time in seconds. If given, this
            is inserted as a header keyword
            EXPTIME. 
    Returns:
    """
    hdr = cutout.wcs.to_header()
    if exptime is not None:
        assert type(exptime)==float, "exposure time should be a float"
        hdr['EXPTIME'] = exptime
        # Note to user: Make sure the image is not normalized when
        # inserting an exposure time. i.e. if the reduced image is normalised,
        # multiply it by the EXPTIME and then pass on the EXPTIME to
        # write_cutout.
    imghdu = fits.PrimaryHDU(cutout.data, hdr)
    hdulist = fits.HDUList([imghdu])
    hdulist.writeto(filename, overwrite=overwrite)
    return


def _genconf(imgfile:str, psffile:str=None, mode=0,
            configfile:str=None, cdkfile:str=None, outdir:str=None, outfile:str=None,
            noisefile:str=None,
            finesample:int = 1, badpix:str = "none",
            constraints:str = "none",
            region:tuple = None, convobox:tuple = (100,100),
            zeropoint:float = 25.0,
            position:tuple = None, int_mag:float = None,
            r_e:float = None, n:float = 1.0, axis_ratio:float = 0.5,
            pa:float = 0, skip_sky:bool = False)->str:
    """
    Creates a configuration file for GALFIT. Conventionally,
    GALFIT is run using the command: `galfit <config-file>`.
    Args:
        imgfile (str): path to the image fits file.
        psffile (str, optional): path to the PSF model fits file.
            If nothing is given, the fit is performed without
            a PSF model.
        mode (int, optional): 0=optimize, 1=model, 2=imgblock, 3=subcomps.
            Which mode would you like galfit to run in?
        outdir (str): Name of output directory. Default
            value is 'galfit_out` in the current directory.
        configfile (str, optional): path to configuration file to
            be created via this function. Defaults to `galfit.feedme`
            in <outdir>.
        cdkfile (str, optional): Path to the charge diffusion kernel
            file. Useful is you have oversampled HST PSFs from 
            Tiny Tim.
        outfile (str, optional): name of GALFIT's output fits file.
            Defaults to `out.fits` in <outdir>.
        noisefile (str, optional): If you'd like to use an
            image noise estimate of your own instead of what Galfit
            generates. 
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
        pa (float, optional): Initial guess for the 
            angle of the major axis counter clockwise 
            relative to the vertical.
        skip_sky (bool, optional): Do you also want
            to fit a constant sky background? Set to 
            false if your sky background is 0.
    Returns:
        configfile (str): Path to the configuration
            file.
    """
    
    # Run checks for file paths and use default paths when
    # none are given.
    if outdir is None:
        outdir = "galfit_out"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if configfile is None:
        warnings.warn("Creating a configuration file here")
        configfile = "galfit.feedme"
    
    # Copy PSF file to outdir
    if isinstance(psffile, str):
        os.system("cp {:s} {:s}".format(psffile, os.path.join(outdir,"psffile.fits")))
        psffile = "psffile.fits"
    # If CDK file exists, do the same
    if cdkfile is not None:
        os.system("cp {:s} {:s}".format(cdkfile, os.path.join(outdir,"cdkfile.txt")))
        cdkfile = "cdkfile.txt"
    #Also for noise file
    if noisefile is not None:
        os.system("cp {:s} {:s}".format(noisefile, os.path.join(outdir,"noisefile.fits")))
        noisefile = "noisefile.fits"
    
    
    # Begin writing config file
    with open(os.path.join(outdir,configfile),"w+") as fstream:
        #Image parameters.
        fstream.write("""===============================================================================\n
        # IMAGE and GALFIT CONTROL PARAMETERS\n
        """)
        img, hdr = fits.getdata(imgfile, header=True)
        # Image file
        fstream.write("A) {:s}  # Input data image (FITS file)\n".format(imgfile))
        if outfile is None:
            warnings.warn("Creating output file here")
            outfile = "out.fits"
        # Output file
        fstream.write("B) {:s}  # Output data image block\n".format(outfile))
        # Include sigma image if provided
        fstream.write(f"C) {noisefile}  # Sigma image name (made from data if blank or 'none')\n")
        # PSF file
        if cdkfile is None:
            fstream.write("D) {:s}  # Input PSF file\n".format(str(psffile)))
        else:
            fstream.write("D) {:s}  # Input PSF file\n".format(str(psffile)+" "+cdkfile))
        # PSF fine-sampling
        fstream.write("E) {:d}  # PSF fine sampling factor\n".format(finesample))
        # Bad pixel map
        fstream.write("F) {:s}  #Bad pixel mask\n".format(badpix))
        # Parameter constraints file
        fstream.write("G) {:s}  # File with parameter constraints (ASCII file)\n".format(constraints))
        # Fit region in pixels
        if region is None:
            # Fit full region by default
            warnings.warn("Using full region to fit. May fail if other objects in the region aren't masked.")
            region = (0, img.shape[1]-1, 0, img.shape[0]-1)

        fstream.write("H) {:d} {:d} {:d} {:d}   # Image region to fit (xmin xmax ymin ymax)\n".format(region[0],region[1],region[2], region[3]))

        # Convolution box
        assert len(convobox)==2, "Only two integers specify convolution box dimensions"
        fstream.write("I) {:d} {:d} # Size of convolution box (x y)\n".format(convobox[0],convobox[1]))
        # Photometric zero point for the image
        fstream.write("J) {:f}  # Photometric zeropoint (mag)\n".format(zeropoint))
        # Platescale of the image
        # Assuming square pixels
        platescale = get_platescale(WCS(hdr))
        fstream.write("K) {:f} {:f} # Plate scale (dx dy) [arcsec/pixel]\n".format(platescale,platescale))
        # Display type? Not sure what this is.
        fstream.write("O) regular   # Display type (regular, curses, both)\n")
        # Only fit for now.
        fstream.write(f"P) {mode} # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")
        ################################################################################3
        # Fit parameters
        fstream.write("\n# INITIAL FITTING PARAMETERS\n")
        # Only 1 sersic component for now.
        fstream.write("# Component number: 1\n")
        fstream.write("0) sersic # Component type\n")

        # Centroid
        if position is None:
            pos_x, pos_y = int(img.shape[0]/2), int(img.shape[1]/2)  
        else:
            pos_x, pos_y = position
        fstream.write("1) {:d} {:d} 1 1 # position x y\n".format(pos_x,pos_y))

        # What is the integrated magnitude?
        if int_mag is None:
            warnings.warn("No guess given for integrated magnitude. Proceeding with sum within region.")
            cropped_img = img[region[2]:region[3],region[0]:region[1]]
            # Better guess if other sources are masked out in badpix
            if badpix !='none':
                bad_pix = fits.getdata(badpix)
                # Convert to boolean
                bad_pix = np.where(bad_pix==0, True, False)
                # Crop to region
                bad_pix = bad_pix[region[2]:region[3],region[0]:region[1]]
            else:
                bad_pix = True
            int_mag = zeropoint-2.5*np.log10(np.sum(cropped_img[bad_pix]))
        fstream.write("3) {:f} 1 # Integrated magnitude\n".format(int_mag))
        # Effective radius
        if r_e is None:
            warnings.warn("Guess for r_e not given. This might not converge.")
            r_e = (region[3]-region[2]+region[1]-region[0])/2/3 #A third of the average region dimension
        fstream.write("4) {:f} 1 # effective radius (pix)\n".format(r_e))
        # Sersic index
        fstream.write("5) {:f} 1 # sersic index\n".format(n))
        # Blanks
        for num in [6,7,8]:
            fstream.write("{:d}) 0.0000 0 # ----\n".format(num))
        # b/a
        fstream.write("9) {:f} 1 # Axis ratio (b/a)\n".format(axis_ratio))
        # PA of the galaxy
        fstream.write("10) {:f} 1 # Position angle (PA) [deg: Up=0, left =90]\n".format(pa))
        # Don't skip sersic fitting at any cost.
        fstream.write("Z) 0 # Skip this model? (yes=1,no=0)\n\n")
        # Skip sky fitting though?
        if not skip_sky:
            # Just a flat sky. No fancy gradients.
            fstream.write("# Component number: 2\n")
            fstream.write("0) sky # component type\n")
            _, median, _ = sigma_clipped_stats(img)
            fstream.write("1) {:f} 1 # Sky background\n".format(median))
            fstream.write("2) 0 0 # dsky/dx\n")
            fstream.write("3) 0 0 # dsky/dy\n")
            fstream.write("Z) 0 # Skip this model\n")
        fstream.write("================================================================================\n")
    # Done.
    return configfile

def read_fitlog(outfile:str, initfile:str, twocomponent:bool=False)->dict:
    """
    Reads the output fit.log
    file and returns the sersic fit parameters
    Args:
        outfile (str): Path and name of the fit log
            file.
        initfile (str): Path to the config file
            used to produce the log entry. This
            is used to distinguish multiple
            entried in the logfile.
    Returns:
        pix_dict (dict): A dict containing
            the GALFIT best fit parameters
            and their uncertainties.
    """
    # TODO: create a regex string to look for blocks
    # in the fit.log file. Then create an expression
    # to lock for model sub-blocks within each block. 
    lines = [line.rstrip('\n') for line in open(outfile)]
    instance = []  # going to put all instances of use of input file here
    for kk, line in enumerate(lines):  # for all lines in fit.log
        if 'Init. par. file : ' + str(initfile.split("/")[-1]) in line:  # if the it is from the input file used
            for pp, item in enumerate(lines[kk:kk + 10]):
                if 'sersic' in item:  # look for the sersic fit model
                    # This assumes each fitting had a single
                    # sersic model fit.
                    # TODO: accommodate more complex fits.
                    instance.append((item, pp, lines[kk:kk + 10][pp + 1]))  # and keep the model and error
    # Use regex to identify all floats in fit_out and err_out
    float_regex = r"[+-]?\d+\.\d+" # Any string of the form <int>.<int> (which is just a float)

    if not twocomponent: # Assumes a single component fit
        # only keep the latest one in case there are multiple runs
        fit_out, _, err_out = instance[-1]

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
    else: # If there are two components
        fit_out_1, _, err_out_1 = instance[-2]
        fit_out_2, _, err_out_2 = instance[-1]

        # Find all instances of float in the strings and convert from string to float
        fitparams_1 = np.array(re.findall(float_regex, fit_out_1)).astype('float')
        fiterrs_1 = np.array(re.findall(float_regex, err_out_1)).astype('float')
        fitparams_2 = np.array(re.findall(float_regex, fit_out_2)).astype('float')
        fiterrs_2 = np.array(re.findall(float_regex, err_out_2)).astype('float')

        fitparams = np.vstack([fitparams_1, fitparams_2])
        fiterrs = np.vstack([fiterrs_1, fiterrs_2])
        #Store values
        pix_dict = {'x':fitparams[:,0], 'y':fitparams[:,1],
                'mag':fitparams[:,2], 'reff':fitparams[:,3],
                'n':fitparams[:,4], 'b/a':fitparams[:,5], 'PA':fitparams[:,6],
                'x_err':fiterrs[:,0], 'y_err':fiterrs[:,1],
                'mag_err':fiterrs[:,2], 'reff_err':fiterrs[:,3],
                'n_err':fiterrs[:,4], 'b/a_err':fiterrs[:,5], 'PA_err':fiterrs[:,6],
                }
    return pix_dict

def pix2coord(pix_dict:dict, wcs:WCS, table:bool=False,multicomponent:bool=False)->dict:
    """
    Takes the output table from galfit's
    fit.log file and converts all pixel
    measurements to physical measurements.
    Args:
        pix_dict (dict): Raw sersic
            fit parameter dict from fit.log
        wcs (WCS): Input image WCS.
        table (bool, optional): Return as
            a table instead?
        multicomponent (bool, optional): If
            true, this snippet accounts out for
            multiple sersic profiles in the same
            model.
    Returns:
        sky_dict (dict/Table): pix_dict
            translated to angular units
            on the sky.
    """

    sky_dict = {}
    platescale = get_platescale(wcs)

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
    pix_matrix = wcs.pixel_scale_matrix
    cdelt1 = np.sqrt(np.sum(pix_matrix[:,0]**2))
    sin_theta = pix_matrix[1,0]/cdelt1
    cos_theta = pix_matrix[0,0]/cdelt1
    north_angle = 180+np.arctan2(sin_theta,cos_theta)*180/np.pi

    sky_dict['PA'] = (pix_dict['PA']-north_angle)%360
    sky_dict['PA_err'] = pix_dict['PA_err']

    if table:
        if multicomponent:
            sky_dict = Table(sky_dict)
        else:
            sky_dict = Table([sky_dict])
        # Reorder because dict-> Table messes it up
        sky_dict = sky_dict['ra','ra_err','dec','dec_err','mag',
                            'mag_err','reff_ang','reff_ang_err',
                            'n','n_err','b/a','b/a_err','PA','PA_err']
    return sky_dict

def run(imgfile:str, psffile:str=None, **kwargs)->int:
    """
    Run galfit. 
    
    Args:
        imgfile (str): path to the image fits file.
        psffile (str): path to the PSF model fits file.

    Valid kwargs:
        mode (int, optional): 0=optimize, 1=model, 2=imgblock, 3=subcomps.
            Which mode would you like galfit to run in?
        outdir (str): Name of output directory. Default
            value is 'galfit_out` in the current directory.
        configfile (str, optional): path to configuration file to
            be created via this function. Defaults to `<outdir>/galfit.feedme`.
        cdkfile (str, optional): Path to the charge diffusion kernel
            file. Useful is you have oversampled HST PSFs from 
            Tiny Tim.
        outfile (str, optional): path to GALFIT's output fits file.
            Defaults to `<outdir>/out.fits`
        noisefile (str, optional): If you'd like to use an
            image noise estimate of your own instead of what Galfit
            generates. 
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
        pa (float, optional): Initial guess for the 
            angle of the major axis counter clockwise 
            relative to the vertical.
        skip_sky (bool, optional): Do you also want
            to fit a constant sky background? Set to 
            false if your sky background is 0.
    Returns:
        fit_outcome (int): An int encoding the success/failure
            of the fitting procedure. 0 corresponds
            to a successful fit. See GALFIT documentation
            to learn what other values stand for.
    """

    # Check input paths first
    assert os.path.isfile(imgfile), "Invalid image file path {:s}".format(imgfile)

    if isinstance(psffile, str):
        assert os.path.isfile(psffile), "Invalid psf file path {:s}".format(psffile)
        psffile = os.path.abspath(psffile)

    # Is a CDK file given?
    if 'cdkfile' in kwargs.keys():
        assert os.path.isfile(kwargs['cdkfile']), "Invalid CDK file path {:s}".format(kwargs['cdkfile'])
        # Use absolute path henceforth
        kwargs['cdkfile'] = os.path.abspath(kwargs['cdkfile'])

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
    
    
    # Generate Galfit config file
    configfile = _genconf(imgfile, psffile,**kwargs)

    # Go to the output directory
    if 'outdir' not in kwargs:
        kwargs['outdir'] = "galfit_out"
    # Keep track of the original directory
    curdir = os.path.abspath(os.path.curdir)
    # Move to outdir
    os.chdir(kwargs['outdir'])
    # Run galfit
    fit_outcome = os.system("galfit {:s}".format(configfile))
    if fit_outcome!=0:
        # TODO: This doesn't work right now because
        # Galfit is still returning 0 on crashing.
        # Something broken?
        warnings.warn("Something went wrong with the fit. Check terminal output.")
        os.chdir(curdir)
        return fit_outcome
    # Read fit.log and get the fit results
    # Temporary fix for the crash:
    if kwargs['mode']==0:
        try:
            pix_dict = read_fitlog("fit.log", configfile)
        except FileNotFoundError:
            print("""
                Doh!  GALFIT crashed because at least one of the model parameters 
                is bad.  The most common causes are: effective radius too small/big,
                component is too far outside of fitting region (also check fitting
                region), model mag too faint, axis ratio too small, Sersic index
                too small/big, Nuker powerlaw too small/big.  If frustrated or 
                problem should persist, email for help or report problem to: 
                                    Chien.Y.Peng@gmail.com 


                GALFIT Version 3.0.5 -- Apr. 23, 2013
                """)
            os.chdir(curdir)
            return 1
        # Convert to sky angular units and stuff it into the output
        # fits file
        hdr = fits.getheader(imgfile)
        wcs = WCS(hdr)
        sky_tab = pix2coord(pix_dict, wcs, table=True)
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
    else:
        os.chdir(curdir)
    return fit_outcome