import os
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
from frb.associate import bayesian

import photutils

from IPython import embed

class FRBAssociate():

    def __init__(self, frb, image_file=None, max_radius=1e9):
        self.frb = frb
        self.image_file = image_file
        self.max_radius = max_radius

        # Attributes
        self.hdu = None
        self.wcs = None
        self.theta_max = None
        self.Pchance = None

    @property
    def sigR(self):
        return np.sqrt(self.frb.sig_a * self.frb.sig_b) * units.arcsec

    @property
    def frb_eellipse(self):
        return dict(a=self.frb.sig_a,
                    b=self.frb.sig_b,
                    theta=self.frb.eellipse['theta'])

    def load_image(self):

        if self.image_file is None:
            raise IOError("Set image_file before calling this method")
        self.hdu = fits.open(self.image_file)[0]
        self.wcs = astropy_wcs.WCS(self.hdu.header)
        self.header = self.hdu.header

    def calc_pchance(self, ndens_eval='bloom'):
        self.Pchance = bayesian.pchance(self.candidates[self.filter].data,
                                        self.candidates['separation'].to('arcsec').value,
                                        self.candidates['half_light'].value,
                                        self.sigR.to('arcsec').value, ndens_eval=ndens_eval)
        # Add to table
        self.candidates['P_c'] = self.Pchance

    def calc_priors(self, prior_S, method='linear'):

        if self.Pchance is None:
            raise IOError("Set Pchance before calling this method")

        if prior_S < 0.:
            self.prior_S = np.product(self.candidates['P_c'])
        else:
            self.prior_S = prior_S
        # Raw priors
        self.raw_prior_Mi = bayesian.raw_prior_Mi(self.Pchance, method)

        # Normalize
        self.prior_Mi = bayesian.renorm_priors(self.raw_prior_Mi, self.prior_S)

        # Add to table
        self.candidates['P_M'] = self.prior_Mi

    def calc_PMx(self):
        """
        Calculate p(M|x) by running through
        the series of:
            self.calc_pxM()
            self.calc_pxS()
            self.calc_px()

        Values are stored in self.P_Mix
        and the candidates table as P_Mx
        """

        # Intermediate steps
        self.calc_pxM()
        self.calc_pxS()
        self.calc_px()

        # Finish
        self.P_Mix = self.prior_Mi * self.p_xMi / self.p_x
        self.candidates['P_Mx'] = self.P_Mix

        # P(S|x)
        self.P_Sx = self.prior_S * self.p_xS / self.p_x

    def calc_pxM(self):
        """
        Calculate p(x|M) and assign to p_xMi
        """
        self.p_xMi = bayesian.px_Mi(self.max_radius,
                              self.frb.coord, self.frb_eellipse,
                              self.candidates['coords'],
                              self.theta_prior)

    def calc_pxS(self):
        """
        Calculate p(x|S) and assign to p_xS
        """
        self.p_xS = bayesian.px_Mi(self.max_radius,
                                   self.frb.coord,
                                   self.frb_eellipse,
                                   SkyCoord([self.frb.coord]),
                                   self.theta_prior)[0]

    def calc_px(self):
        self.p_x = self.prior_S * self.p_xS + np.sum(self.prior_Mi * self.p_xMi)

    def cut_candidates(self, plate_scale, bright_cut=None, separation=None):

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
        self.candidates = self.photom[cands]

        # Add coords
        coords = astropy_wcs.utils.pixel_to_skycoord(
            self.candidates['xcentroid'].value,
            self.candidates['ycentroid'].value,
            self.wcs)
        # Insist on ICRS
        coords = coords.transform_to('icrs')

        self.candidates['ra'] = coords.ra.value
        self.candidates['dec'] = coords.dec.value
        self.candidates['coords'] = coords

        # Separation
        seps = self.frb.coord.separation(self.candidates['coords'])
        self.candidates['separation'] = seps.to('arcsec')

        # Cut on separation?
        if separation is not None:
            cut_seps = seps < separation
            self.candidates = self.candidates[cut_seps]

        # Half light
        self.candidates['half_light'] = self.candidates['semimajor_axis_sigma'].value * plate_scale

    def photometry(self, ZP, filter, radius=3., show=False, outfile=None):

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

        # Magnitudes
        self.filter = filter
        self.photom = self.cat.to_table()
        self.photom[filter] = -2.5 * np.log10(self.photom['source_sum']) + ZP

        # Plot?
        if show or outfile is not None:
            norm = ImageNormalize(stretch=SqrtStretch())
            fig = plt.figure(figsize=(9, 9))

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
            xy_kernel:
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
            fig = plt.figure(figsize=(7, 7))

            ax = plt.gca()
            cmap = self.segm.make_cmap()
            ax.imshow(self.segm, origin='lower', cmap=cmap, interpolation='nearest')
            ax.set_title('Segmentation Image')

            if outfile is not None:  # This must come first
                plt.savefig(outfile, dpi=300)
            if show:
                plt.show()

    def set_theta_prior(self, theta_dict):
        self.theta_max = theta_dict['max']
        self.theta_prior = theta_dict

    def threshold(self, nsig=1.5):

        if self.hdu is None:
            self.load_image()

        # Background
        bkg_estimator = photutils.MedianBackground()
        self.bkg = photutils.Background2D(self.hdu.data, (50, 50), filter_size=(3, 3),
                                          bkg_estimator=bkg_estimator)

        # Threshold
        self.thresh_img = self.bkg.background + (nsig * self.bkg.background_rms)

    def __repr__(self):
        txt = '<{:s}: {}'.format(self.__class__.__name__, self.frb.frb_name)
        # Finish
        txt = txt + '>'
        return (txt)


def run_individual(config, show=False, skip_bayesian=False, verbose=False):
    """

    Returns:
        FRBAssociate:

    """
    # FRB
    FRB = frb.FRB.by_name(config['name'])

    # FRB Associate
    frbA= FRBAssociate(FRB, max_radius=config['max_radius'])

    # Image
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
    frbA.cut_candidates(config['plate_scale'], separation=config['cand_separation'])

    # Chance probability
    frbA.calc_pchance()
    if verbose:
        print(frbA.candidates[['id', config['filter'], 'half_light', 'separation', 'P_c']])

    # Return here?
    if skip_bayesian:
        return frbA

    # Finish
    return frbA
