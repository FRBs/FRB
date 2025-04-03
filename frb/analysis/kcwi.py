"""
Utilities to handle KCWI datacubes
"""

from textwrap import fill
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
    from spectral_cube import BooleanArrayMask
except ImportError:
    raise ImportError("Requirement unmet: SpectralCube. Run `pip install spectral-cube`.")
try:
    import pyregion as pyreg
except ImportError:
    raise ImportError("Requirement unmet: pyregion. Run `pip install pyregion`.")

import os

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
    transtable = Table.read(transfile, **readkw)
    wave, transval = transtable['col1'].data, transtable['col2'].data
    trans = interp1d(wave, transval, kind=kind, fill_value=fill_value)
    return trans

def silence_warnings(warncategory):
    """
    To silence spectral cube warnings.
    Check out spectral_cube.utils for
    warning categories.
    Args:
        warncategory (Warnings category): category of Warnings you want
            to silence.
    """
    warnings.filterwarnings('ignore', category=warncategory, append=True)
    return

def _air_to_vac(wave):
        """
        Helper function to convert wavelengths
        in air to vacuum.
        Args:
            wave (float or ndarray): wavelength in air in angstroms
        Returns:
            wave_vac (float or ndarray): wavelength in vacuum in angstroms
        """
        sigma2 = (1e4/wave)**2
        fact = 1.+5.792105e-2/(238.0185-sigma2)+1.67917e-3/(57.362-sigma2)
        wave_vac = wave*fact
        return wave_vac

def _read_in_data(datafile, reduction_pipeline='pypeit'):
        """
        Read in data from a KCWI datacube.            flux_hdr.set('CRVAL3', flux_hdr['CRVAL3']*1e10, "[Angstrom] Reference value for wavelength")
            flux_hdr.set('CDELT3', flux_hdr['CDELT3']*1e10, "[Angstrom] Coordinate increment at reference point")        Args:
            datafile (str): Path to datacube fits file.
            reduction_pipeline (str, optional): Reduction pipeline
                used to reduce the datacube. Default is 'pypeit'.
                Also allowed 'kcwidrp'. Other pipelines are not
                supported yet.
        Returns:
            flux (array): 3D array of fluxes.
            std (array): 3D Std. dev. array.
            wcs (WCS): WCS object.
            mask (array): 3D boolean mask.
        """
        hdulist = fits.open(datafile)
        if reduction_pipeline == 'pypeit':
            flux, flux_hdr = hdulist['FLUX'].data, hdulist['FLUX'].header.copy()
            std = hdulist['SIG'].data
            flux_hdr.set('CRVAL3', flux_hdr['CRVAL3']*1e10, "[AA] Reference value for wavelength")
            flux_hdr.set('CDELT3', flux_hdr['CDELT3']*1e10, "[AA] Coordinate increment at reference point")
            flux_hdr.set('CUNIT3', 'Angstrom', "Units of coordinate increment and value")
            wcs = WCS(flux_hdr)
            mask = BooleanArrayMask(np.logical_not(hdulist['BPM'].data.astype('bool')), wcs=wcs)
            fluxunit = 1e-17*u.erg*u.s**-1*u.cm**-2*u.AA**-1
        elif reduction_pipeline == 'kcwidrp':
            flux, flux_hdr = hdulist['PRIMARY'].data, hdulist['PRIMARY'].header
            std = hdulist['UNCERT'].data
            wcs = WCS(hdulist['PRIMARY'].header)
            mask = BooleanArrayMask(np.logical_not(hdulist['MASK'].data.astype('bool')), wcs=wcs)
            fluxunit = 1e-16*u.erg*u.s**-1*u.cm**-2*u.AA**-1
        else:
            raise ValueError("Reduction pipeline not supported")
        return flux, flux_hdr, std, wcs, mask, fluxunit

def find_sources(imgfile, nsig = 1.5,
                    minarea = 10., clean=True, deblend_cont = 0.0001, regfile = None, write = None, bkgsub = True):
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
                used for object deblending. Default is 0.0001.
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
            # Compute background again
            bkg = sep.Background(data, mask = mask)


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

class KCWIDatacube():
    """
    A class to handle KCWI datacubes from PypeIt reductions.
    """

    def __init__(self, datafile,
                 reduction_pipeline='pypeit',
                 allow_huge_operations=True,**kwargs):

        # Read in data        
        flux, flux_hdr, std, wcs, mask, fluxunit = _read_in_data(datafile, reduction_pipeline=reduction_pipeline)
        self.cube = SpectralCube(data = flux, wcs = wcs,
                                mask = mask, header = flux_hdr,
                                meta = {'reduction_pipeline': reduction_pipeline,
                                        'datafile': os.path.abspath(datafile),
                                        'fluxunit': fluxunit},
                                allow_huge_operations=allow_huge_operations,
                                fill_value=np.nan,
                                **kwargs)
        self.varcube = SpectralCube(data = std**2, wcs = wcs, mask=mask, header = flux_hdr, fill_value=np.nan,
                                    allow_huge_operations = allow_huge_operations)

    def _spectral_tile(self, array):
        """
        Helper function that tiles 1D array of
        size cube.spectral_axis.shape and tiles it
        to the same shape as the cube.
        Args:
            array (ndarray): the 1D array that needs to be tiled.
        Returns:
            tiled_array (ndarray): an array with the 1D array tiled
                in the spatial dimension. This has the same shape
                as the cube.
        """
        tiled_array = np.tile(array, (self.cube.shape[2],self.cube.shape[1],1)).T
        return tiled_array

    def wave_mask(self, mask_1d):
        """
        Mask out wavelengths using
        a 1D grid. Values corresponding
        to "False" are masked out.
        Args:
            mask_1D (bool ndarray): 1D boolean array
                of same length as cube.spectral_axis
        Returns:
            masked_cube (Spectral cube): masked datacube.
            masked_varcube (Spectral cube): masked variance cube.
        """
        assert len(mask_1d) == len(self.cube.spectral_axis), "Mask length ({:d}) doesn't match cube's spectral axis length ({:d}).".format(len(mask_1d), len(self.cube.spectral_axis))

        assert mask_1d.dtype == bool, "Mask must be a boolean array."

        mm = self._spectral_tile(mask_1d)
        return self.cube.with_mask(mm), self.varcube.with_mask(mm)


    def get_img(self, wlow = None, whigh = None, 
                transfile = None, how = "cube",
                bkgsub = False, save = None,
                overwrite = False, bkgsubkw= {}, trans_readkw = {}):    
        """
        Flatten cube along wavelength and produce a 2D
        image.
        Args:
            wlow, whigh (Quantity, optional): wavelength 
                limits (with astropy units) to flatten between.
                If nothing is given, the cube is checked for the
                WAVGOOD keywords and flattened between
                them. If they don't exist, it's flattened fully.
            transfile (str, optional): file containing transmission
                curve as a function of wavelength. Wavelength should
                be in angstroms and transmission in fractions.
            how (str, optional): "cube", "slice" or "ray". How do
                you want to load the cube to memory?
                "cube" loads the whole thing for summing. "slice"
                and "ray" do it slicewise or spectral-ray-wise.
            bkgsub (bool, optional): Subtract background continuum?
            bkgsubkw (dict, optional): Keyword args to be passed to sep.Background
                for background estimation.
            trans_readkw (dict, optional): Keyword arguments for reading
                the transmission curve file using astropy Table.
            save (str, optional): Path to file to be
                saved to.
            overwrite (bool, optional): Overwrite existing
                file?
        Returns:
            img (Spectral Cube Projection): Flattened 2D image
        """
        assert how in ["cube", "slice", "ray"], "Invalid summing method. Choose one of 'cube', 'slice' and 'ray'."

        # Create a truncated cube based on wlow and whigh
        if wlow == None:
            wlow = self.cube.spectral_extrema[0]
        
        if whigh == None:
            whigh = self.cube.spectral_extrema[1]
        
        goodcube = self.cube.spectral_slab(wlow, whigh)

        # Do you want to use a filter?
        if transfile:
            # Compute transmission curve for cube wavelengths
            transfunc = _interp_trans(transfile, **trans_readkw)
            trans = transfunc(goodcube.spectral_axis.value)

            # Create a 3D array of trans stacked in the same
            # shape as the cube spatial dimensions.
            tt = self._spectral_tile(trans)
            goodcube = goodcube*tt
        
        # Make image
        img = goodcube.sum(axis = 0, how = how)
        if bkgsub:
            bkg = sep.Background(img.value, **bkgsubkw)
            img = img - bkg*img.unit
        if save:
            img.write(save, overwrite = overwrite)
        # TODO: Also return the std image, mask?
        return img

    def spec_from_spatial_mask(self, mask_arr, how="cube"):
        """
        Extract a spectrum from a cube
        within a spatial mask.
        Args:
            cube (Spectral Cube): A datacube object
            mask_arr (numpy array): A 2D boolean array. A
                spectrum is extracted from the regions
                corresponding to True.
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

        masked_cube = self.cube.subcube_from_mask(mask_arr)

        #TODO: add more methods of obtaining the central estimate
        spec = masked_cube.mean(axis=(1,2), how = how)
        # if kind is max:
        # if kind is something_else:
        masked_var = self.varcube.subcube_from_mask(mask_arr)
        # TODO: Find out how to estimate the variance of a sample median.
        varspec = masked_var.mean(axis = (1,2), how = how)
        return spec, varspec

    def spec_from_ellipse(self,
                        x0 = 0., y0 = 0., a = 1.,
                        b = 1., theta = 0., r = 1.):
        
        """
        Get the spectrum within an elliptical region
        Args:
            cube (Spectral Cube): A datacube object
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
        mask = np.zeros(self.shape[1:], dtype=np.bool)
        sep.mask_ellipse(mask, x0, y0, a, b, theta,r)
        
        return self.spec_from_spatial_mask(mask)


    def _make_marz(self, speclist, varspeclist, objects,marzfile="marzfile.fits", tovac = True):
        """
        Helper function to create a MARZ input file
        Args:
            cube (Spectral cube): Datacube
            speclist (list): list of spectra. i.e.
                Spectral cube Projections of same
                shape as cube.spectral_axis.
            varspeclist (list): list of spectral
                variances (same object class and
                shape as speclist elements)
            objects (astropy Table): Table of
                objects detected using sep.extract
            marzfile (str, optional): Name of
                output MARZ fits file.
            tovac (bool, optional): Convert wavelengths
                to vacuum?
        """
        # TODO: Actually compute sky background
        nobjs = len(objects)
        wcsinfo = self.wcs.celestial

        sky = np.zeros_like(speclist)

        # Convert wavelengths from air to vacuum
        wave = self.spectral_axis.value
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
        return marz_hdu
    def get_source_spectra(self, objects, outdir = "spectra/", marzfile = None, tovac = True):
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

        wave = self.spectral_axis.value

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
            spec, varspec = self.spec_from_ellipse(obj['x'], obj['y'],
                                            obj['a'], obj['b'],
                                            obj['theta'], r = 2)

            # Produce spectrum fits file
            spechdu = fits.PrimaryHDU(spec.data, header=spec.header)
            spechdu.header.set('extname', 'SPEC')
            varhdu = fits.ImageHDU(varspec.data, header=varspec.header)
            varhdu.header.set('extname', 'VAR')
            hdulist = fits.HDUList([spechdu, varhdu, wavehdu])
            specfile_name = os.path.join(outdir,str(idx)+"_spec1d.fits")
            hdulist.writeto(specfile_name, overwrite=True)

            # Append spectrum to list
            speclist[idx] = spec.data
            varspeclist[idx] = varspec.data

        if marzfile:
            marzfile = os.path.join(outdir, marzfile)
            self._make_marz(speclist=speclist,
                            varspeclist=varspeclist,
                            objects=objects,
                            marzfile=marzfile, tovac=tovac)
        return speclist, varspeclist, wave
