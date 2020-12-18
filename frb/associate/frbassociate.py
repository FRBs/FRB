import warnings
import numpy as np
from matplotlib import pyplot as plt

from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy import units
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from astropy import wcs as astropy_wcs

from frb import frb
from astropath import bayesian
from astropath import chance
from frb.galaxies import photom
from frb.galaxies import nebular

import photutils

from IPython import embed

class FRBAssociate():
    """

    Attributes:
        hdu:
        photom (pandas.DataFrame):  Photometry table
        candidate (pandas.DataFrame):  Candidates table
            Note, while this is derived from photom, it is a *separate* copy
        Pchance (np.ndarray): Chance probability
        Sigma_m (np.ndarray): Surface density of sources on the sky
    """

    def __init__(self, frb, image_file=None, max_radius=1e9):
        """

        Args:
            frb:
            image_file:
            max_radius (float, optional):
                Maximum radius for analysis (arcsec)
        """
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

    def calc_pchance(self, ndens_eval='driver', extinction_correct=False):
        """
        Calculate the Pchance values for the candidates

        self.Pchance filled in place
        Addes as P_c to candidates

        Args:
            ndens_eval:

        Returns:

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
                                        self.candidates['half_light'],
                                        self.sigR.to('arcsec').value, ndens_eval=ndens_eval)

        # Add to table
        self.candidates['P_c'] = self.Pchance
        self.candidates['Sigma_m'] = self.Sigma_m

    def calc_priors(self, prior_U, method='linear'):
        """
        Calulate the priors based on Pchance
        and the input method

        prior_Mi and prior_S are set in place
        The candidates table is updated with P_O

        Args:
            prior_S (float):
                If the input value is <0, use the P_c product else use the input value
            method (str, optional):
                'linear':  P(O) = 1 - Pchance
                'inverse':  P(O) = 1 / Pchance
        """
        if self.Pchance is None:
            raise IOError("Set Pchance before calling this method")

        # TODO -- Move this into Bayesian
        if prior_U < 0.:
            self.prior_U = np.product(self.candidates['P_c'])
        else:
            self.prior_U = prior_U
        # Raw priors
        self.raw_prior_Oi = bayesian.raw_prior_Oi(self.Pchance, self.Sigma_m, method)

        # Normalize
        self.prior_Oi = bayesian.renorm_priors(self.raw_prior_Oi, self.prior_U)

        # Add to table
        self.candidates['P_O'] = self.prior_Oi

    def calc_POx(self, **kwargs):
        """
        Calculate p(O|x) by running through
        the series of:
            self.calc_pxO()
            self.calc_pxS()
            self.calc_px()

        Values are stored in self.P_Oix
        and the candidates table as P_Ox
        """

        # Intermediate steps
        self.calc_pxO(**kwargs)
        self.calc_pxU()
        self.calc_px()

        # Finish
        self.P_Oix = self.prior_Oi * self.p_xOi / self.p_x
        self.candidates['P_Ox'] = self.P_Oix

        # P(S|x)
        self.P_Ux = self.prior_U * self.p_xU / self.p_x

    def calc_pxO(self, **kwargs):
        """
        Calculate p(x|O) and assign to p_xOi

        self.p_xOi is set in place
        """

        # This hack sets the minimum localization to 0.2''
        # TODO -- Do this better
        eellipse = self.frb.eellipse.copy()
        eellipse['a'] = max(self.frb.sig_a, 0.2)
        warnings.warn("Need to improve the hack above")

        # Do it
        self.p_xOi = bayesian.px_Oi(self.max_radius,
                                    self.frb.coord,
                                    eellipse,
                                    self.candidates['coords'].values,
                                    self.theta_prior, **kwargs)

    def calc_pxU(self):
        """
        Calculate p(x|U) and assign to p_xU
        """
        self.p_xU = bayesian.px_Oi(self.max_radius,
                                   self.frb.coord,
                                   self.frb_eellipse,
                                   SkyCoord([self.frb.coord]),
                                   self.theta_prior)[0]

    def calc_px(self):
        """
        Calculate p(x)

        Returns:
            np.ndarray:

        """
        self.p_x = self.prior_U * self.p_xU + np.sum(self.prior_Oi * self.p_xOi)

    def cut_candidates(self, plate_scale, bright_cut=None, separation=None):
        """
        Cut down to candidates

        self.candidates is made in place

        Args:
            plate_scale:
            bright_cut:
            separation:

        Returns:

        """

        # Zero point
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

        self.candidates['ra'] = coords.ra
        self.candidates['dec'] = coords.dec
        self.candidates['coords'] = coords

        # Separation
        seps = self.frb.coord.separation(coords)
        self.candidates['separation'] = seps.to('arcsec')

        # Cut on separation?
        if separation is not None:
            cut_seps = seps < separation
            self.candidates = self.candidates[cut_seps]

        # Half light
        self.candidates['half_light'] = self.candidates['semimajor_axis_sigma'] * plate_scale

    def photometry(self, ZP, ifilter, radius=3., show=False, outfile=None):
        """
        Perform photometry

        Fills self.photom in place
        
        Half-light radii:
            https://iopscience.iop.org/article/10.1086/444475/pdf

        Args:
            ZP (float):
                Zero point magnitude
            ifilter (str):
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

        self.cat = photutils.source_properties(self.hdu.data - self.bkg.background,
                                               self.segm,
                                               kron_params=('mask', 2.5, 1.0, 'exact', 5),
                                               background=self.bkg.background,
                                               filter_kernel=self.kernel)

        # Apertures
        apertures = []
        for obj in self.cat:
            position = np.transpose((obj.xcentroid.value, obj.ycentroid.value))
            a = obj.semimajor_axis_sigma.value * radius
            b = obj.semiminor_axis_sigma.value * radius
            theta = obj.orientation.to(units.rad).value
            apertures.append(photutils.EllipticalAperture(position, a, b, theta=theta))
        self.apertures = apertures

        # Magnitudes
        self.filter = ifilter
        self.photom = self.cat.to_table().to_pandas()
        self.photom[ifilter] = -2.5 * np.log10(self.photom['source_sum']) + ZP

        # Kron
        for key in ['kron_radius']:
            self.photom[key] = getattr(self.cat, key)

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

    def run_bayes(self, prior, fill_half_light=True, verbose=True):
        """
        Run the main steps of the Bayesian analysis

        self.POx is filled in place
        and added to self.candidates

        Args:
            prior (dict):
            fill_half_light (bool, optional):
                Add half_light radii to the prior dict
            verbose (bool, optional):

        Returns:

        """
        # Checks
        if self.Pchance is None:
            raise IOError("Set Pchance before calling this method")
        # Priors
        self.calc_priors(prior['S'], method=prior['M'])
        # Fill?
        if fill_half_light:
            prior['theta']['r_half'] = self.candidates['half_light'].values
        self.set_theta_prior(prior['theta'])
        # P(O|x)
        self.calc_POx()
        #
        if verbose:
            print("All done with Bayes")

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
        self.segm = photutils.detect_sources(self.hdu.data, self.thresh_img,
                                             npixels=npixels, filter_kernel=self.kernel)

        # Debelnd?
        if deblend:
            segm_deblend = photutils.deblend_sources(self.hdu.data, self.segm,
                                                     npixels=npixels,
                                                     filter_kernel=self.kernel,
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

    def set_theta_prior(self, theta_dict):
        """

        Args:
            theta_dict (dict):

        """
        self.theta_max = theta_dict['max']
        self.theta_prior = theta_dict

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
        self.bkg = photutils.Background2D(self.hdu.data, box_size,
                                          filter_size=filter_size,
                                          bkg_estimator=bkg_estimator)

        # Threshold
        self.thresh_img = self.bkg.background + (nsig * self.bkg.background_rms)

    def view_candidates(self):
        """
        Convenience method to show candidate table

        Args:
            nsig:
            box_size:
            filter_size:

        Returns:

        """
        items = ['id', self.filter, 'half_light', 'separation', 'P_c']
        for add_on in ['P_O', 'P_Ox']:
            if add_on in self.candidates.keys():
                items += [add_on]
        print(self.candidates[items])

    def __repr__(self):
        txt = '<{:s}: {}'.format(self.__class__.__name__, self.frb.frb_name)
        # Finish
        txt = txt + '>'
        return (txt)


def run_individual(config, show=False, skip_bayesian=False, verbose=False):
    """
    Run through the steps leading up to Bayes

    Args:
        config:
        show:
        skip_bayesian:
        verbose:

    Returns:

    """
    # FRB
    FRB = frb.FRB.by_name(config['name'])

    # FRB Associate
    frbA= FRBAssociate(FRB, max_radius=config['max_radius'])

    # Image
    print("Using image {}".format(config['image_file']))
    hdul = fits.open(config['image_file'])
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

    # Threshold + Segment
    frbA.threshold()
    frbA.segment(deblend=config['deblend'], npixels=config['npixels'], show=show)

    # Photometry
    frbA.photometry(config['ZP'], config['filter'], show=show)
    if verbose:
        print(frbA.photom[['xcentroid', 'ycentroid', config['filter']]])

    # Candidates
    frbA.cut_candidates(config['plate_scale'], bright_cut=config['cand_bright'],
                        separation=config['cand_separation'])

    # Chance probability
    frbA.calc_pchance(ndens_eval='driver')

    if verbose:
        print(frbA.candidates[['id', config['filter'], 'half_light', 'separation', 'P_c']])

    # Return here?
    if skip_bayesian:
        return frbA

    # BAYESIAN GOES HERE....

    # Finish
    return frbA
