"""
Utilities to handle KCWI datacubes
"""

import numpy as np
import warnings

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, utils as wcsutils
from astropy import units as u
from scipy.interpolate import interp1d

try:
    import sep
except ImportError:
    raise ImportError("Requirement unmet: sep. Run `pip install sep`")

try:
    from spectral_cube import SpectralCube
except ImportError:
    raise ImportError("Requirement unmet: SpectralCube. Run `pip install spectral-cube`.")
try:
    import pyregion as pyreg
except ImportError:
    raise ImportError("Requirement unmet: pyregion. Run `pip install pyregion`.")

import glob, os, sys

def _air_to_vac(wave):
    """
    Helper function to convert wavelengths
    in air to vacuum.
    """
    sigma2 = (1e4/wave)**2
    fact = 1.+5.792105e-2/(238.0185-sigma2)+1.67917e-3/(57.362-sigma2)
    return wave*fact 

def _spectral_tile(array, cube):
    """
    Helper function that tiles 1D array of
    size cube.spectral_axis.shape and tiles it
    to the same shape as the cube.
    """
    return np.tile(array, (cube.shape[2],cube.shape[1],1)).T

def _clean_wave(cube):
    """
    If there are "good" wavelengths defined
    in the header, return a subcube filtering
    out the bad wavelengths. 
    """
    if 'WAVEGOOD0' in list(cube.header.keys()):
        wlow = cube.header['WAVEGOOD0']
        whigh = cube.header['WAVEGOOD1']
        clean_cube = cube.spectral_slab(wlow, whigh)
        return clean_cube
    else:
        return cube

def _interp_trans(transfile, kind= "cubic", fill_value=0, **readkw):
    """
    Interpolate a transmission curve from
    a file.
    Args:
        transfile (str): Path to transmission curve file.
            The file should contain two columns in this order:
            1) wavelength in angstroms
            2) intensity transmission in fractions.
        kind (str or int, optional): input to scipy.interpolate.interp1d.
        fill_value (array-like, optional): input to scipy.interpolate.interp1d
            Check out this link for details
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
        **readkw: Keyword arguments for reading the
            input file using astropy Table.
    Returns:
        transfunc (function): interpolation
    """
    wave, transval = Table.read(transfile, **readkw)
    trans = interp1d(wave, transval, kind=kind, fill_value=fill_value)
    return trans

def silence_warnings(warncategory):
    """
    To silence spectral cube warnings.
    Check out spectral_cube.utils for
    warning categories.
    """
    warnings.filterwarnings('ignore', category=warncategory, append=True)

def wave_mask(cube, mask_1d):
    """
    Mask out wavelengths using
    a 1D grid. Values corresponding
    to "False" are masked out.
    """
    assert len(mask_1d) == len(cube.spectral_axis), "Mask length ({:d}) doesn't match cube's spectral axis length ({:d}).".format(len(mask_1d), len(cube.spectral_axis))

    assert mask_1d.dtype == bool, "Mask must be a boolean array."

    mm = _spectral_tile(mask_1d, cube)
    return cube.with_mask(mm)


def get_img(cubefile, wlow = None, whigh = None, trans_curve = None, save = None, overwrite = False):
    """
    Flatten cube along wavelength and produce a 2D
    image.
    Args:
        cubefile (str): Path to the datacube
        wlow, whigh (Quantity, optional): wavelength 
            limits (with astropy units) to flatten between.
            If nothing is given, the cube is checked for the
            WAVGOOD keywords and flattened between
            them. If they don't exist, it's flattened fully.
        filter (function, optional): transmission
            curve as a function of wavelength. Should be able
            to take vector inputs and produce vector outputs.
            We recommend passing a function produced by scipy
            interpolation method. Wavelength is assumed to
            be in angstroms.
        save (str, optional): Path to file to be
            saved to.
        overwrite (bool, optional): Overwrite existing
            file?
    Returns:
        img (Spectral Cube Projection): Flattened 2D image
    """
    # Read in datacube
    cube = SpectralCube.read(cubefile)

    # Create a truncated cube based on wlow and whigh
    if not wlow:
        try:
            wlow = cube.header['WAVGOOD0']*cube.spectral_axis.unit
        except KeyError:
            wlow = cube.spectral_extrema[0]
    
    if not whigh:
        try:
            whigh = cube.header['WAVGOOD1']*cube.spectral_axis.unit
        except KeyError:
            whigh = cube.spectral_extrema[1]
    
    goodcube = cube.spectral_slab(wlow, whigh)

    # Do you want to use a filter?
    if trans_curve:
        # Compute transmission curve for cube wavelengths
        trans = trans_curve(goodcube.spectral_axis.value)

        # Create a 3D array of trans stacked in the same
        # shape as the cube spatial dimensions.
        # TODO: make this more elegant.
        tt = _spectral_tile(trans, cube)
        goodcube = goodcube*tt
    
    # Make image
    img = goodcube.sum(axis = 0, how = "ray")
    if save:
        img.write(save, overwrite = overwrite)
    
    return img

def spec_from_mask(cube, mask_arr, varcube=None, kind="mean", how="cube"):
    """
    Extract a spectrum from a cube
    within a mask.
    Args:
        cube (Spectral Cube): A datacube object
        mask_arr (numpy array): A 2D boolean array
        varcube (Spectral Cube, optional): Variance cube
        kind (str, optional): median or mean
        how (str, optional): "cube" or "slice". Load 
            the entire masked cube to memory when 
            computing spectra or compute it looping
            over slices.
    Returns:
        spec (OneDSpectrum): The extracted spectrum.
        var (OneDSpectrum): Variance in spectrum.
            Only returned if varcube is supplied.
    """
    assert mask_arr.dtype == bool, "Masks must be boolean. int type masks make computation slow."

    assert how in ["cube", "slice"], "You can either take the full cube or compute the spectrum slice-by-slice."

    masked_cube = cube.subcube_from_mask(mask_arr)
    masked_cube = _clean_wave(masked_cube)

    #TODO: add more methods of obtaining the central estimate
    if kind is "mean":
        spec = masked_cube.mean(axis=(1,2), how = how)
    # if kind is max:
    # if kind is something_else:

    if varcube:
        masked_var = varcube.subcube_from_mask(mask_arr)
        var = _clean_wave(masked_cube)
        var = masked_var.mean(axis = (1,2), how = how)
        return spec, var
    
    return spec

def spec_from_ellipse(cube, varcube = None,
                      x0 = 0., y0 = 0., a = 1.,
                      b = 1., theta = 0., r = 1.):
    
    """
    Get the spectrum within an elliptical region
    Args:
        cube (Spectral Cube): A datacube object
        varcube (Spectral Cube, optional): Variance cube
        x0, y0 (float, optional): Centroid of ellipse
        a, b (float, optional): semi-major and semi-minor axes
        theta (float, optional): rotation angle of the semi-major
            axis from the positive x axis.
        r (float, optional): Scaling factor for a and b.
            If not 1, a = r*a and b = r*b. 
    Returns:
        spec (OneDSpectrum): The extracted spectrum.
        var (OneDSpectrum): Variance in spectrum.
            Only returned if varcube is supplied.
    """
    #TODO: Use photutils aperture object for this.
    # Create 2D mask first
    mask = np.zeros(cube.shape[1:], dtype=np.bool)
    sep.mask_ellipse(mask, x0, y0, a, b, theta,r)
    
    return spec_from_mask(cube, mask, varcube)

def find_sources(imgfile, nsig = 1.5, minarea = 10., clean=True, deblend_cont = 0.0001, regfile = None, write = None, bkgsub = True):
    """
    Find sources in the whitelight image
    using SExtractor.
    Args:
        imgfile (str): An image fits file
        n_sig (float, optional): Detection threshold in units
            of sky background rms.
        minarea (float, optional): minimum area in pixels
            to be considered a valid detection.
        clean (bool, optional): Perform cleaning?
        deblend_cont (float, optional): Minimum contrast ratio
            used for object deblending. Default is 0.005.
            To entirely disable deblending, set to 1.0.
        regfile (str, optional): A ds9 region file of
            areas to be masked out.
        write (str, optional): write extracted object table
            to this path.
        bkgsub (bool, optional): perform background subtraction?
            Default is set to true.
    Returns:
        objects (Table): Summary table of detected objects.
        segmap (ndarray): Segmentation map.
    """
    # Get whitelight image
    hdulist = fits.open(imgfile)
    white = hdulist[0]

    data = white.data
    data = data.byteswap().newbyteorder() # sep requires this
    
    # Make a mask if available
    if regfile:
        reg = pyreg.open(regfile).as_imagecoord(white.header)
        mask = reg.get_filter().mask(data)
    else:
        mask = None

    # Characterize sky
    bkg = sep.Background(data, mask = mask)

    # Subtract background? 
    if bkgsub:
        bkg.subfrom(data)

    # Compute source detection threshold
    thresh = nsig*bkg.globalrms
    # Extract sources
    objects, segmap = sep.extract(data, thresh = thresh, mask = mask,
                                  deblend_cont = deblend_cont,
                                  minarea = minarea, clean = clean,
                                  segmentation_map=True)
    if write:
        Table(objects).write(write, overwrite = True)
    
    return Table(objects), segmap

def _make_marz(cube, speclist, varspeclist, objects,marzfile="marzfile.fits", tovac = True):
    """
    Helper function to create a MARZ input file
    """
    # TODO: Actually compute sky background
    nobjs = len(objects)
    wcsinfo = cube.wcs.celestial

    sky = np.zeros_like(speclist)

    # Convert wavelengths from air to vacuum
    wave = cube.spectral_axis.value
    if tovac:
        wave = _air_to_vac(wave)
    wavelist = np.tile(wave, (nobjs,1))

    # Set infs to nan
    speclist[np.isinf(speclist)] = np.nan
    varspeclist[np.isinf(varspeclist)] = np.nan

    # Create HDUs
    extnames = ['INTENSITY', 'VARIANCE', 'SKY', 'WAVELENGTH']
    datalists = [speclist, varspeclist, sky, wavelist]

    marz_hdu = fits.HDUList()
    for ext, data in zip(extnames, datalists):
        hdu = fits.ImageHDU(data)
        hdu.header.set('extname', ext)
        marz_hdu.append(hdu)

    # Create object table

    ids = np.arange(nobjs)+1
    x = objects['x'].data
    y = objects['y'].data
    coords = wcsutils.pixel_to_skycoord(x,y,wcsinfo)
    ra = coords.ra.value
    dec = coords.dec.value
    types = ('P_'*nobjs).split('_')[:-1]

    cols = []
    colnames = ['source_id', 'RA', 'DEC', 'X', 'Y', 'TYPE']
    formats = ['80A', 'D', 'D', 'J', 'J', '1A']
    coldata = [ids, ra, dec, x, y, types]
    for cname, form, cdat in zip(colnames, formats, coldata):
        cols.append(fits.Column(cname, form, array = cdat))
    coldefs = fits.ColDefs(cols)
    tabhdu = fits.BinTableHDU.from_columns(coldefs)
    tabhdu.header.set('extname', 'FIBRES')
    marz_hdu.append(tabhdu)

    marz_hdu.writeto(marzfile, overwrite = True)
    return

def get_source_spectra(cubefile, varfile, objects, outdir = "spectra/", marzfile = None, tovac = True):
    """
    Extract spectra of sources found using SExtractor
    from datacube.
    Args:
        cubefile (str): A datacube fits file
        varfile (str): Variance datacube fits file
        objects (Table): Table of extracted objects produced
            by sep.extract
        outdir (str, optional): directory to store spectra
        marzfile (str, optional): name of MARZ file to dump
            all spectra into. File creation is skipped if
            a name is not supplied.
        tovac (bool, optional): Covert wavelengths to vacuum.
    Returns:
        speclist (ndarray): A 2D array of with an extracted
            spectrum in each row.
        varspeclist (ndarray): Similarly designed array with
            variance information.
        wave (1D Quantity array): Wavelength array.
    """

    # Preliminaries
    nobjs = len(objects)
    cube = SpectralCube.read(cubefile)
    varcube = SpectralCube.read(varfile)

    wave = cube.spectral_axis.value

    # Convert to vacuum wavelengths?
    if tovac:
        wave = _air_to_vac(wave)
    # Prepare an HDU in advance
    wavehdu = fits.ImageHDU(wave)
    wavehdu.header.set('extname', 'WAVELENGTH')

    # Initialize output lists 
    speclist = np.zeros([nobjs, len(wave)])
    varspeclist = np.zeros_like(speclist)

    # Create output folder?
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for idx, obj in enumerate(objects):
        # TODO: Make a specdb file instead of saving
        # individual spectra
        spec, varspec = spec_from_ellipse(cube, varcube,
                                          obj['x'], obj['y'],
                                          obj['a'], obj['b'],
                                          obj['theta'], r = 2)

        # Produce spectrum fits file
        spechdu = fits.PrimaryHDU(spec.data, header=spec.header)
        spechdu.header.set('extname', 'SPEC')
        varhdu = fits.ImageHDU(varspec.data, header=varspec.header)
        varhdu.header.set('extname', 'VAR')
        hdulist = fits.HDUList([spechdu, varhdu, wavehdu])
        specfile_name = outdir+str(idx)+"_spec1d.fits"
        hdulist.writeto(specfile_name, overwrite=True)

        # Append spectrum to list
        speclist[idx] = spec.data
        varspeclist[idx] = varspec.data

    if marzfile:
        _make_marz(cube, speclist, varspeclist, objects, outdir+marzfile, tovac=tovac)
    return speclist, varspeclist, wave
