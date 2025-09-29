"""
Utilities to handle KCWI datacubes
"""

import numpy as np
import warnings

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.wcs import WCS, utils as wcsutils
from astropy import units as u
from scipy.interpolate import interp1d

try:
    from photutils.segmentation import SourceCatalog, SourceFinder
    from photutils.aperture import EllipticalAperture
    from photutils.background import Background2D, MedianBackground
except ImportError:
    raise ImportError("Requirement unmet: photutils. Run `pip install photutils`")

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
                    minarea = 10.,
                    deblend_cont = 0.0001, cont_levels = 32,
                    connectivity = 8,
                    regfile = None, write = None):
        """
        Find sources in the whitelight image
        using Photutils segmentation.
        Args:
            imgfile (str): An image fits file
            n_sig (float, optional): Detection threshold in units
                of sky background rms.
            minarea (float, optional): minimum area in pixels
                to be considered a valid detection.
            deblend_cont (float, optional): Minimum contrast ratio
                used for object deblending. Default is 0.0001.
                To entirely disable deblending, set to 1.0.
            cont_levels (int, optional): Number of levels
                for deblending.
            connectivity (int, optional): Connectivity for
                segmentation map. Default is 8.
            regfile (str, optional): A ds9 region file of
                areas to be masked out.
            write (str, optional): write extracted object table
                to this path.
        Returns:
            objects (Table): Summary table of detected objects.
            segmap (ndarray): Segmentation map.
        """
        # Get whitelight image
        hdulist = fits.open(imgfile)
        white = hdulist[0]

        data = white.data
        wcs = WCS(white.header)
        #data = data.astype(data.dtype.newbyteorder("=")) # sep requires this
        # Keep this here just in case you want to migrate back to sep in the
        # future 
        finder = SourceFinder(npixels=minarea, contrast=deblend_cont,
                              progress_bar=False,nlevels=cont_levels,
                              connectivity=connectivity)
        
        _, med, std = sigma_clipped_stats(data)
        threshold = med+nsig*std
        
        # Make a mask if available
        if regfile:
            reg = pyreg.open(regfile).as_imagecoord(white.header)
            mask = reg.get_filter().mask(data)
        else:
            mask = None

        segmap = finder(data, threshold=threshold, mask=mask)
        # Extract sources
        objects = SourceCatalog(data, segmap, wcs = wcs).to_table()
        if write:
            Table(objects).write(write, overwrite = True)
        
        return objects, segmap

def splice_spectra(blue_spec, red_spec, blue_wave, red_wave):
    """
    Splice together spectra of the same objects from the blue
    and red channels.
    """

    # Clean up the spectra so that the region of overlap is
    # set to NaN.
    delta_lambda_blue = np.nanmedian(np.diff(blue_wave))
    min_red_wav = red_wave[0]
    max_blue_wav = blue_wave[-1]
    # Create a wave grid for the overlap region
    overlap_wave = np.arange(max_blue_wav, min_red_wav, delta_lambda_blue)

    # Concatentate the wavelengths and fluxes
    spliced_wave = np.concatenate((blue_wave, overlap_wave, red_wave))
    splice_spec = np.concatenate((blue_spec, np.full(overlap_wave.shape, np.nan), red_spec))
    return spliced_wave, splice_spec

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
                as the cube (nwave, nx, ny).
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
            bkgsubkw (dict, optional): Keyword args to be passed to Background2D
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
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            bkg = Background2D(img.value, sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator, **bkgsubkw)
            img = img - bkg.background*img.unit
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
        spec = masked_cube.sum(axis=(1,2), how = how)
        # if kind is max:
        # if kind is something_else:
        masked_var = self.varcube.subcube_from_mask(mask_arr)
        varspec = masked_var.sum(axis = (1,2), how = how) # Var of the sum = Sum of var assuming each pix is uncorrelated (not true, really).
        return spec, varspec

    def spec_from_ellipse(self,
                        x0 = 0., y0 = 0., a = 1.,
                        b = 1., theta = 0.):
        
        """
        Get the spectrum within an elliptical region
        Args:
            cube (Spectral Cube): A datacube object
            x0, y0 (float, optional): Centroid of ellipse in pixels.
            a, b (float, optional): semi-major and semi-minor axes in pixels.
            theta (float, optional): rotation angle of the semi-major
                axis from the positive x axis.
        Returns:
            spec (OneDSpectrum): The extracted spectrum.
            var (OneDSpectrum): Variance in spectrum.
                Only returned if varcube is supplied.
        """
        #TODO: Use photutils aperture object for this.
        # Create 2D mask first
        ellipse = EllipticalAperture((x0, y0), a, b, theta)
        # Create a mask from the ellipse
        aper_mask = ellipse.to_mask(method='center').to_image(self.cube.shape[1:]).astype(bool)
        
        return self.spec_from_spatial_mask(aper_mask)


    def _make_marz(self, speclist, varspeclist, objects,  wave=None, marzfile="marzfile.fits", vac_wave = True):
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
                objects detected using SourceCatalog.
            wave (ndarray, optional): Wavelength array.
                If None, the cube.spectral_axis
                is used.
            marzfile (str, optional): Name of
                output MARZ fits file.
            vac_wave (bool, optional): Is the wavelength
                in vaccum? True by default for KCWI.
        """
        # TODO: Actually compute sky background
        nobjs = len(objects)
        wcsinfo = self.cube.wcs.celestial

        sky = np.zeros_like(speclist)

        # Convert wavelengths from air to vacuum
        if wave is None:
            wave = self.cube.spectral_axis.value
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
        
        if vac_wave:
            marz_hdu['INTENSITY'].header.set('VACUUM', 'True', "Wavelength is in vacuum?")
        else:
            marz_hdu['INTENSITY'].header.set('VACUUM', 'False', "Wavelength is in vacuum?")

        # Create object table

        ids = np.arange(nobjs)+1
        x = objects['xcentroid'].data
        y = objects['ycentroid'].data
        coords = wcsutils.pixel_to_skycoord(x,y,wcsinfo)
        ra = coords.ra.to('rad').value
        dec = coords.dec.to('rad').value
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
    def get_source_spectra(self, objects, outdir = "spectra/", marzfile = None, vac_wave=True, wvlims= None):
        """
        Extract spectra of sources found using SExtractor
        from datacube.
        Args:
            cubefile (str): A datacube fits file
            varfile (str): Variance datacube fits file
            objects (Table): Table of extracted objects produced
                by SourceCatalog.
            outdir (str, optional): directory to store spectra
            marzfile (str, optional): name of MARZ file to dump
                all spectra into. File creation is skipped if
                a name is not supplied.
            vac_wave (bool, optional): Wavelengths are in vacuum.
            wvlims (tuple, optional): Wavelength limits to
                extract spectra between. If None, the full
                wavelength range is used. Expects floats in angstroms.
        Returns:
            speclist (ndarray): A 2D array of with an extracted
                spectrum in each row.
            varspeclist (ndarray): Similarly designed array with
                variance information.
            wave (1D Quantity array): Wavelength array.
        """

        # Preliminaries
        nobjs = len(objects)

        wave = self.cube.spectral_axis.to('AA').value

        # Mask out wavelengths outside the limits
        if wvlims!=None:
            assert len(wvlims) == 2, "wavelength limits should be a tuple of length 2"
            # Check if the spectrum is within the limits
            wvmask = (wave >= wvlims[0]) & (wave <= wvlims[1])
            if np.all(~wvmask):
                raise ValueError("No pixels in the spectrum are within the wavelength limits.")
            else:
                wave = wave[wvmask]
                
        else:
            wvmask = np.ones_like(wave, dtype=bool)
        wavehdu = fits.ImageHDU(wave)
        wavehdu.header.set('extname', 'WAVELENGTH')


        # Initialize output lists
        speclist = np.zeros([nobjs, len(wave)])
        varspeclist = np.zeros_like(speclist)

        # Create output folder?
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        for obj in objects:
            # Instead of relying on the x,y coordinates,
            # Use the sky centroid coordinates to define the ellipse.
            # This way, forced extraction works even if the object
            # detection is from a different KCWI image.
            x, y = wcsutils.skycoord_to_pixel(obj['sky_centroid'], self.cube.wcs.celestial)
            spec, varspec = self.spec_from_ellipse(x, y,
                                            2*obj['semimajor_sigma'].value, 2*obj['semiminor_sigma'].value,
                                            obj['orientation'])
                    

            # Produce spectrum fits file
            spechdu = fits.PrimaryHDU(spec.value[wvmask], header=spec.header)
            spechdu.header.set('extname', 'SPEC')
            varhdu = fits.ImageHDU(spec.value[wvmask], header=varspec.header)
            varhdu.header.set('extname', 'VAR')
            hdulist = fits.HDUList([spechdu, varhdu, wavehdu])
            specfile_name = os.path.join(outdir,str(obj['label'])+"_spec1d.fits")
            hdulist.writeto(specfile_name, overwrite=True)

            # Append spectrum to list
            speclist[obj['label']-1] = spec.value[wvmask]
            varspeclist[obj['label']-1] = varspec.value[wvmask]

        if marzfile:
            marzfile = os.path.join(outdir, marzfile)
            self._make_marz(speclist=speclist,
                            varspeclist=varspeclist,
                            objects=objects,
                            wave=wave,
                            marzfile=marzfile, vac_wave=vac_wave)
        return speclist, varspeclist, wave
