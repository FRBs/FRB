import warnings
import numpy as np
from matplotlib import pyplot as plt

from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clipped_stats
from astropy import units
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from astropy.visualization import LogStretch
from astropy import wcs as astropy_wcs
from astropy.table import Table
from astropy.coordinates import SkyCoord
import importlib_resources

import pandas 

from frb import frb
import astropy_healpix
from astropath import chance
from astropath import path
from frb.galaxies import photom
from frb.galaxies import nebular
from frb.galaxies.galfit import write_cutout

import photutils

from IPython import embed

class FRBAssociate(path.PATH):
    """
    Class that guides the PATH analysis for an FRB
    It is instantiated with an FRB object 

    This is mainly used to do the candidate slurping
    from images.

    Use PATH methods for the Bayesian analysis
    

    Args:
        frb (frb.frb.FRB): FRB object
        image_file (str, optional): Name of image file
        max_radius (float, optional): Maximum radius for analysis (arcsec)

    Attributes:
        hdu (fits.HDU: FITS header-data unit
        photom (pandas.DataFrame):  Photometry table
        candidate (pandas.DataFrame):  Candidates table
            Note, while this is derived from photom, it is a *separate* copy
        Pchance (np.ndarray): Chance probability
    """

    def __init__(self, frb, image_file=None, max_radius=1e9):
        """

        """
        # Init PATH
        super().__init__()

        self.frb = frb
        self.image_file = image_file
        self.max_radius = max_radius

        # Attributes
        self.hdu = None
        self.wcs = None
        self.theta_max = None
        self.Pchance = None
        self.Sigma_m = None
        self.theta_prior = None

        self.photom = None
        self.candidates = None

        # Internals
        self.exc_per = 10.  # exclude_percentile for 2D Background

    @property
    def sigR(self):
        return np.sqrt(self.frb.sig_a * self.frb.sig_b) * units.arcsec

    @property
    def frb_eellipse(self):
        return dict(a=self.frb.sig_a,
                    b=self.frb.sig_b,
                    theta=self.frb.eellipse['theta'])

    def load_image(self):
        """
        Load the image from self.image_file

        Returns:

        """

        if self.image_file is None:
            raise IOError("Set image_file before calling this method")
        self.hdu = fits.open(self.image_file)[0]
        self.wcs = astropy_wcs.WCS(self.hdu.header)
        self.header = self.hdu.header

    def make_host_cutout(self, imgdata, wcs, size=5. * units.arcsec)->Cutout2D:
        """
        Make a cutout of the image around the FRB
        and write to the data directory under "Galaxies".

        Args:
            imgdata (np.ndarray): Image data
            wcs (astropy.wcs.WCS): WCS of the image
            size (Quantity, optional): Size of the cutout
        
        Returns:
            Cutout2D: Cutout of the image
        """

        # Cutout
        host_gal = self.frb.grab_host()
        cutout = Cutout2D(imgdata, host_gal.coord, size, wcs=wcs)
        _, med, std = sigma_clipped_stats(cutout.data)

        # Make the figure
        fig,ax = plt.subplots(1, 1, figsize=(6, 6), subplot_kw={'projection': cutout.wcs})
        norm = ImageNormalize(stretch=LogStretch(), vmin=med+std, vmax=med+20*std)

        ax.imshow(cutout.data, origin='lower', cmap='hot', norm=norm)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec', labelpad=-2)
        ax.set_title(self.frb.frb_name+" Host")
        ax.set_xlim(0, cutout.data.shape[1])
        ax.set_ylim(0, cutout.data.shape[0])

        # Write
        output_file = importlib_resources.files(f'frb.data.Galaxies.{self.frb.frb_name[3:]}')/f'{self.frb.frb_name}_cutout.png'
        plt.subplots_adjust(left=0.2, bottom=0.15, right=0.95, top=0.95)
        fig.savefig(output_file, dpi=300)

        return cutout

    def calc_pchance(self, ndens_eval='driver', extinction_correct=False):
        """
        Calculate the Pchance values for the candidates

        self.Pchance filled in place
        Added as P_c to candidates Table

        Args:
            ndens_eval (str, optional): Source of number density evaluation.
                See frb.associate.chance.pchance for options
            extinction_correct (bool, optional):
                If True, apply an extinction correction to the photometry based
                on the coordinates.

        """
        # Correct for extinction
        if extinction_correct and self.filter+'_orig' not in self.candidates.keys():
            ebv = nebular.get_ebv(self.frb.coord, definition='SandF')['meanValue']
            linear_ext = photom.extinction_correction(self.filter, ebv)
            self.candidates[self.filter+'_orig'] = self.candidates[self.filter].values.copy()
            self.candidates[self.filter] += -2.5*np.log10(linear_ext)
        # Do it
        self.Pchance, self.Sigma_m = chance.pchance(self.candidates[self.filter],
                                        self.candidates['separation'],
                                        self.candidates['ang_size'],
                                        self.sigR.to('arcsec').value, ndens_eval=ndens_eval)

        # Add to table
        self.candidates['P_c'] = self.Pchance
        self.candidates['Sigma_m'] = self.Sigma_m


    def cut_candidates(self, plate_scale, 
                       bright_cut:float=None, 
                       separation:float=None):
        """
        Cut down to candidates

        self.candidates is made in place

        Args:
            plate_scale (float or str):
                If str -- grab the value from the Header with this key
                If float -- use this value (arcsec)
            bright_cut (float, optional):
                Cut objects on this magnitude
            separation (float, optional):
                Cut objects on this angular separation (arcsec)

        """

        # Plate scale
        if isinstance(plate_scale, str):
            plate_scale = self.header[plate_scale]

        if self.photom is None:
            raise ValueError("photom table not built!")

        cands = np.ones(len(self.photom), dtype=bool)

        # Cut on brightness?
        if bright_cut is not None:
            good_bright = self.photom[self.filter] > bright_cut
            cands &= good_bright

        # Candidate table
        self.candidates = self.photom[cands].copy()

        # Add coords
        coords = astropy_wcs.utils.pixel_to_skycoord(
            self.candidates['xcentroid'],
            self.candidates['ycentroid'],
            self.wcs)
        # Insist on ICRS
        coords = coords.transform_to('icrs')

        self.candidates['ra'] = coords.ra.value
        self.candidates['dec'] = coords.dec.value
        self.candidates['coords'] = coords

        # Separation
        seps = self.frb.coord.separation(coords)
        self.candidates['separation'] = seps.to('arcsec').value

        # Cut on separation?
        if separation is not None:
            cut_seps = seps < separation
            self.candidates = self.candidates[cut_seps]

        # Half light
        self.candidates['ang_size'] = self.candidates['semimajor_sigma'] * plate_scale

    def photometry(self, ZP, ifilter, radius=3., show=False, outfile=None):
        """
        Perform photometry

        Fills self.photom in place
        
        Half-light radii:
            https://iopscience.iop.org/article/10.1086/444475/pdf

        Args:
            ZP (float):
                Zero point magnitude
            ifilter (str): Filter name to be used in the anaysis
            radius (float, optional):
                Scaling for semimajor/minor axes for Elliptical apertures
            show:
            outfile:
        """

        # Init
        if self.segm is None:
            raise ValueError("segm not set!")
        if self.hdu is None:
            self.load_image()

        # Zero point
        if isinstance(ZP, str):
            ZP = self.header[ZP]

        self.cat = photutils.segmentation.SourceCatalog(
            self.hdu.data - self.bkg.background,
            self.segm,
            background=self.bkg.background)

        # Apertures
        apertures = []
        for obj in self.cat:
            position = np.transpose((obj.xcentroid, obj.ycentroid))
            a = obj.semimajor_sigma.value * radius
            b = obj.semiminor_sigma.value * radius
            theta = obj.orientation.to(units.rad).value
            apertures.append(photutils.aperture.EllipticalAperture(position, a, b, theta=theta))
        self.apertures = apertures

        # Magnitudes
        self.filter = ifilter
        self.photom = Table(self.cat.to_table()).to_pandas()
        self.photom[ifilter] = -2.5 * np.log10(self.photom['segment_flux']) + ZP

        # Add in ones lost in the pandas conversion Kron
        for key in ['kron_radius']:
            self.photom[key] = getattr(self.cat, key).value # pixel

        # Plot?
        if show or outfile is not None:
            norm = ImageNormalize(stretch=SqrtStretch())
            fig = plt.figure(figsize=(6, 6))

            # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
            plt.clf()
            ax1 = plt.gca()
            ax1.imshow(self.hdu.data, origin='lower', cmap='Greys_r', norm=norm)
            ax1.set_title('Data')
            #
            for aperture in apertures:
                aperture.plot(axes=ax1, color='white', lw=1.5)
            if outfile is not None:  # This must come first
                plt.savefig(outfile, dpi=300)
            if show:
                plt.show()


    def segment(self, nsig=3., xy_kernel=(3,3), npixels=3, show=False, outfile=None,
                deblend=False):
        """
        Generate the segment image

        Args:
            nsig (float):
                Kernel parameter
            xy_kernel:
                Kernel parameter
            npixels:
            show:
            outfile:
            deblend (bool, optional):
                Run deblend algorithm too

        Returns:

        """

        if self.thresh_img is None:
            raise ValueError("Threshold image not set!")
        if self.hdu is None:
            self.load_image()

        # Kernel
        sigma = nsig * gaussian_fwhm_to_sigma
        self.kernel = Gaussian2DKernel(sigma, x_size=xy_kernel[0], y_size=xy_kernel[1])  # Might need a somewhat larger one
        self.kernel.normalize()

        # Segment
        try:
            self.segm = photutils.segmentation.detect_sources(self.hdu.data, self.thresh_img,
                                             npixels=npixels, 
                                             kernel=self.kernel)
        except TypeError:
            warnings.warn("Support for older photutils versions will be deprecated. Upgrade to >1.8", category=DeprecationWarning)
            convolved_data = convolve(self.hdu.data, self.kernel)
            self.segm = photutils.segmentation.detect_sources(convolved_data, self.thresh_img,
                                             npixels=npixels)

        # Debelnd?
        if deblend:
            try:
                segm_deblend = photutils.segmentation.deblend_sources(self.hdu.data, self.segm,
                                                     npixels=npixels,
                                                     kernel=self.kernel,
                                                     nlevels=32,
                                                     contrast=0.001)
            except TypeError:
                warnings.warn("Support for older photutils versions will be deprecated. Upgrade to >1.8", category=DeprecationWarning)
                segm_deblend = photutils.segmentation.deblend_sources(convolved_data, self.segm,
                                                     npixels=npixels,
                                                     nlevels=32,
                                                     contrast=0.001)
            self.orig_segm = self.segm.copy()
            self.segm = segm_deblend


        # Show?
        if show or outfile is not None:
            fig = plt.figure(figsize=(6, 6))

            ax = plt.gca()
            cmap = self.segm.make_cmap()
            ax.imshow(self.segm, origin='lower', cmap=cmap, interpolation='nearest')
            ax.set_title('Segmentation Image')

            if outfile is not None:  # This must come first
                plt.savefig(outfile, dpi=300)
            if show:
                plt.show()

    def threshold(self, nsig=1.5, box_size=(50,50), filter_size=(3,3)):
        """
        Generate threshold image

        self.thresh_img is set in place

        Args:
            nsig (float, optional):
                Primary threshold parameter
            box_size (tuple):
                Primary Background2D parameter
            filter_size (tuple):
                Primary Background2D parameter
        Returns:

        """

        if self.hdu is None:
            self.load_image()

        # Background
        bkg_estimator = photutils.MedianBackground()
        #embed(header='threshold 402 frbassociate')
        self.bkg = photutils.Background2D(self.hdu.data, box_size,
                                          filter_size=filter_size,
                                          bkg_estimator=bkg_estimator,
                                          exclude_percentile=self.exc_per)

        # Threshold
        self.thresh_img = self.bkg.background + (nsig * self.bkg.background_rms)

    def view_candidates(self):
        """
        Convenience method to show candidate table
        """
        items = ['id', self.filter, 'ang_size', 'separation', 'P_c']
        for add_on in ['P_O', 'P_Ox']:
            if add_on in self.candidates.keys():
                items += [add_on]
        print(self.candidates[items])

    def __repr__(self):
        txt = '<{:s}: {}'.format(self.__class__.__name__, self.frb.frb_name)
        # Finish
        txt = txt + '>'
        return (txt)


def run_individual(config, prior:dict=None, show=False, 
                   skip_bayesian=False, 
                   verbose=False,
                   loc:dict=None,
                   posterior_method:str='fixed',
                   extinction_correct=False,
                   generate_png:bool=False,
                   FRB:frb.FRB=None,
                   internals:dict=None,
                   debug:bool=False):
    """
    Run through the steps leading up to Bayes

    Args:
        config (dict):  Runs the PATH analysis
            keys:
                name (str): Name of the FRB, e.g. FRB20121102
                max_radius (float): Maximum radius in arcsec for the analysis
                cut_size (float): Size to trim the image down to
                deblend (bool): If True, apply the deblending algorithm
                npixels (int): number of pixels in image segmentation
                filter (str): filter name
                ZP (float): Zero point value (magnitudes)
                plate_scale (float): Plate scale in arcsec
                cand_bright (float): Sources brighter than this are assumed stars and ignored
                host_cut_size (float, optional): Size of the host cutout; required if generate_png=True
        prior(dict, optional):
            Contains information on the priors
        posterior_method(str, optional):
            Method of calculation
        loc (dict, optional):
            Contains the localization
        show (bool, optional):
            Show a few things on the screen
        skip_bayesian (bool, optional):
            Skip the Bayesian part, i.e. only do the setup
        extinction_correct (bool, optional):
            If True, correct for Galactic extinction
        FRB (frb.FRB, optional):
            FRB object
        internals(dict, optional):
            Attributes to set in the FRBA object
        generate_png (bool, optional):
            Generate PNGs of the cutouts
        verbose (bool, optional):
    """
    if not skip_bayesian and prior == None:
        raise IOError("Must specify the priors if you are running the Bayesian analysis")
    # FRB
    if FRB is None:
        FRB = frb.FRB.by_name(config['name'])

    # FRB Associate
    frbA= FRBAssociate(FRB, max_radius=config['max_radius'])

    # Internals
    if internals is not None:
        for key in internals.keys():
            setattr(frbA, key, internals[key])
    
    if 'cand_file' in config.keys():
        # Read
        frbA.candidates = pandas.read_csv(config['cand_file'])
        # Coords
        frbA.coords = SkyCoord(ra=frbA.candidates.ra.values,
                               dec=frbA.candidates.dec.values,
                               unit='deg')
        # Add separation?
        if 'separation' not in frbA.candidates.keys():
            seps = frbA.frb.coord.separation(frbA.coords)
            frbA.candidates['separation'] = seps.to('arcsec').value
        
    elif 'image_file' in config.keys():
        # Image
        print("Using image {}".format(config['image_file']))


        # Load image
        hdul = fits.open(config['image_file'])
        # A hack for some image packing
        if 'SCI' in [ihdu.name for ihdu in hdul]:
            hdu_full = hdul['SCI']
        else:
            hdu_full = hdul[0]

        if config['cut_size'] is not None:
            size = units.Quantity((config['cut_size'], config['cut_size']), units.arcsec)
            cutout = Cutout2D(hdu_full.data, FRB.coord, size, wcs=WCS(hdu_full.header))
            frbA.wcs = cutout.wcs
            frbA.hdu = cutout  # not really an HDU
        else:
            frbA.hdu = hdu_full  # not really an HDU
            frbA.wcs = WCS(hdu_full.header)

        frbA.header = hdu_full.header

        # Make a cutout of the host
        if generate_png:
            frbA.make_host_cutout(frbA.hdu.data, wcs = frbA.wcs, size=config['host_cut_size']*units.arcsec)

        # Threshold + Segment
        frbA.threshold()
        frbA.segment(deblend=config['deblend'], npixels=config['npixels'], show=show,
                     xy_kernel=config['xy_kernel'])

        # Photometry
        ZP = config['ZP'] if 'ZP' in config.keys() else frbA.header['ZP']
        frbA.photometry(ZP, config['filter'], show=show)
        if verbose:
            print(frbA.photom[['xcentroid', 'ycentroid', config['filter']]])

        #embed(header='threshold 546 frbassociate')
        # Candidates
        frbA.cut_candidates(config['plate_scale'], bright_cut=config['cand_bright'],
                        separation=config['cand_separation'])

        # Chance probability
        frbA.calc_pchance(ndens_eval='driver', extinction_correct=extinction_correct)

        frbA.candidates['mag'] = frbA.candidates[frbA.filter]

    #if verbose:
    #    print(frbA.candidates[['id', config['filter'], 'ang_size', 'separation', 'P_c']])

    # BAYESIAN 
    if skip_bayesian:
        return frbA

    # Init
    frbA.init_cand_coords()

    # Set priors
    frbA.init_cand_prior('inverse', P_U=prior['U'])
    frbA.init_theta_prior(prior['theta']['method'], 
                            prior['theta']['max'],
                            prior['theta']['scale'])

    # Localization
    if loc is None:
        if 'hpix_file' in config.keys():
            # Load healpix
            hpix = Table.read(config['hpix_file'])
            header = fits.open(config['hpix_file'])[1].header

            nside = 2**header['MOCORDER']

            # Normalize
            healpix = astropy_healpix.HEALPix(nside=nside)
            norm =  np.sum(hpix['PROBDENSITY']) *  healpix.pixel_area.to('arcsec**2').value
            hpix['PROBDENSITY'] /= norm

            # Set
            localiz = dict(type='healpix',
                        healpix_data=hpix, 
                        healpix_nside=nside,
                        healpix_ordering='NUNIQ',
                        healpix_coord='C')            
            frbA.init_localization('healpix', **localiz)
        else:
            frbA.init_localization('eellipse', 
                            center_coord=frbA.frb.coord,
                            eellipse=frbA.frb_eellipse)
    else:                    
        frbA.init_localization(loc['type'], **loc)
    
    # Calculate priors
    frbA.calc_priors()                            

    # Calculate p(O_i|x)
    debug = True
    frbA.calc_posteriors(posterior_method, 
                         box_hwidth=frbA.max_radius,
                         max_radius=frbA.max_radius, # For unseen prior
                         debug=debug)

    # Reverse Sort
    frbA.candidates = frbA.candidates.sort_values('P_Ox', ascending=False)

    # Finish
    return frbA