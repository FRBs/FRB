import numpy as np
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy import units
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from astropy import wcs as astropy_wcs

import photutils

from IPython import embed

class FRBAssociate():

    def __init__(self, frb, image_file=None):
        self.frb = frb
        self.image_file = image_file

        # Attributes
        self.hdu = None
        self.wcs = None

    @property
    def sigR(self):

    def load_image(self):

        if self.image_file is None:
            raise IOError("Set image_file before calling this method")
        self.hdu = fits.open(self.image_file)[0]
        self.wcs = astropy_wcs.WCS(self.hdu.header)

    def cut_candidates(self, plate_scale, bright_cut=None):

        # Zero point
        if isinstance(plate_scale, str):
            plate_scale = self.hdu.header[plate_scale]

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

        self.candidates['ra'] = coords.ra.value
        self.candidates['dec'] = coords.dec.value
        self.candidates['coords'] = coords

        # Separation
        seps = self.frb.coord.separation(self.candidates['coords'])
        self.candidates['separation'] = seps.to('arcsec')

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
            ZP = self.hdu.header[ZP]

        self.cat = photutils.source_properties(self.hdu.data, self.segm, filter_kernel=self.kernel)

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

    def segment(self, nsig=3., xy_kernel=(3,3), npixels=3, show=False, outfile=None):

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

    def threshold(self, nsig=1.5):

        if self.hdu is None:
            self.load_image()

        # Background
        bkg_estimator = photutils.MedianBackground()
        self.bkg = photutils.Background2D(self.hdu.data, (50, 50), filter_size=(3, 3),
                                          bkg_estimator=bkg_estimator)

        # Threshold
        self.thresh_img = self.bkg.background + (nsig * self.bkg.background_rms)

if __name__ == '__main__':
    from frb import frb

    # Testing

    # #########################
    # Prep
    frb180924 = frb.FRB.by_name('FRB180924')

    # Instantiate
    frbA_180924 = FRBAssociate(frb180924, image_file='dev/FRB180924_DESr.fits')

    # Threshold
    frbA_180924.threshold()

    # Segment
    frbA_180924.segment()#show=True, outfile='dev/FRB180924_segm.png')

    # Photometry
    frbA_180924.photometry('MAGZERO', 'r')#, show=True, outfile='dev/FRB180924_aper.png')

    # Candidates
    plate_scale = frbA_180924.hdu.header['CD2_2'] * 3600. * units.arcsec  # arcsec
    frbA_180924.cut_candidates(plate_scale, bright_cut=18.1)


    # #########################
    # Bayesian time

    embed(header='180 of frbassociate')