""" Module for generating a finder chart """
import os
import numpy as np

from IPython import embed

from matplotlib import pyplot as plt
from matplotlib import font_manager
import matplotlib.cm as cm

from PIL import Image

from astropy import units
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.stats import sigma_clipped_stats
from astropy.visualization import LogStretch, mpl_normalize as mplnorm
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import ICRS
from astropy.time import Time

from linetools import utils as ltu

from frb.figures import utils
from frb.surveys import survey_utils
from frb.surveys import images

try:
     from photutils import SkyRectangularAperture
except ImportError:
    flag_photu = False
    print('Install the photutils package to be able to add a slit to an image')
else:
    flag_photu = True


def from_hdu(hdu, title, **kwargs):
    """
    Convenience function to handle an HDU to generate a finder chart

    see generate() for more details

    Args:
        hdul (astropy.io.fits.PrimaryHDU):
        title (str):

    Returns:
        matplotlib.pyplot.figure, matplotlib.pyplot.Axis: figure generated

    """
    image = hdu.data
    wcs = WCS(hdu.header)

    return generate(image, wcs, title, **kwargs)


def generate(image, wcs, title, flip_ra=False, flip_dec=False,
             log_stretch=False,
             cutout=None,
             primary_coord=None, secondary_coord=None,
             third_coord=None, slit=None,
             vmnx=None, extra_text=None, outfile=None, figsize=None):
    """
    Basic method to generate a Finder chart figure

    Args:
        image (np.ndarray):
          Image for the finder
        wcs (astropy.wcs.WCS):
          WCS solution
        title (str):
          Title; typically the name of the primary source
        flip_ra (bool, default False):
            Flip the RA (x-axis). Useful for southern hemisphere finders.
        flip_dec (bool, default False):
            Flip the Dec (y-axis). Useful for southern hemisphere finders.
        log_stretch (bool, optional):
            Use a log stretch for the image display
        cutout (tuple, optional):
            SkyCoord (center coordinate) and Quantity (image angular size)
            for a cutout from the input image.
        primary_coord (astropy.coordinates.SkyCoord, optional):
          If provided, place a mark in red at this coordinate
        secondary_coord (astropy.coordinates.SkyCoord, optional):
          If provided, place a mark in cyan at this coordinate
          Assume it is an offset star (i.e. calculate offsets)
        third_coord (astropy.coordinates.SkyCoord, optional):
          If provided, place a mark in yellow at this coordinate
        slit (tuple, optional):
          If provided, places a rectangular slit with specified
          coordinates, width, length, and position angle on image (from North to East)
          [SkyCoords('21h44m25.255s',-40d54m00.1s', frame='icrs'), 1*u.arcsec, 10*u.arcsec, 20*u.deg]
        vmnx (tuple, optional):
          Used for scaling the image.  Otherwise, the image
          is analyzed for these values.
        extra_text : str
          Extra text to be added at the bottom of the Figure.
          e.g. `DSS r-filter`
        outfile (str, optional):
          Filename for the figure.  File type will be according
          to the extension
        figsize (tuple, optional):
          tuple to define the figure size. It goes within the fig=plt.figure(figsize=figsize)

    Returns:
        matplotlib.pyplot.figure, matplotlib.pyplot.Axis

    """

    utils.set_mplrc()

    plt.clf()
    
    if figsize is not None:
        fig = plt.figure(figsize=figsize)
    else: 
        fig = plt.figure(figsize=(8.5,10.5))
        # fig.set_size_inches(7.5,10.5)

    # Cutout?
    if cutout is not None:
        cutout_img = Cutout2D(image, cutout[0], cutout[1], wcs=wcs)
        # Overwrite
        wcs = cutout_img.wcs
        image = cutout_img.data

    # Axis
    ax = fig.add_axes([0.10, 0.20, 0.75, 0.5], projection=wcs)

    # Show
    if log_stretch:
        norm = mplnorm.ImageNormalize(stretch=LogStretch())
    else:
        norm = None
    cimg = ax.imshow(image, cmap='Greys', norm=norm)

    # Flip so RA increases to the left
    if flip_ra is True:
        ax.invert_xaxis()
    if flip_dec is True:
        ax.invert_yaxis()

    # N/E
    overlay = ax.get_coords_overlay('icrs')

    overlay['ra'].set_ticks(color='white')
    overlay['dec'].set_ticks(color='white')

    overlay['ra'].set_axislabel('Right Ascension')
    overlay['dec'].set_axislabel('Declination')

    overlay.grid(color='green', linestyle='solid', alpha=0.5)

    # Contrast
    if vmnx is None:
        mean, median, stddev = sigma_clipped_stats(image)  # Also set clipping level and number of iterations here if necessary
        #
        vmnx = (median-stddev, median + 2*stddev)  # sky level - 1 sigma and +2 sigma above sky level
        print("Using vmnx = {} based on the image stats".format(vmnx))
    cimg.set_clim(vmnx[0], vmnx[1])

    # Add Primary
    if primary_coord is not None:
        c = SphericalCircle((primary_coord.ra, primary_coord.dec),
                        2*units.arcsec, transform=ax.get_transform('icrs'),
                        edgecolor='red', facecolor='none')
        ax.add_patch(c)
        # Text
        jname = ltu.name_from_coord(primary_coord)
        ax.text(0.5,1.34, jname, fontsize=28,
             horizontalalignment='center',transform=ax.transAxes)

    # Secondary
    if secondary_coord is not None:
        c_S1 = SphericalCircle((secondary_coord.ra, secondary_coord.dec),
                           2*units.arcsec, transform=ax.get_transform('icrs'),
                           edgecolor='cyan', facecolor='none')
        ax.add_patch(c_S1)
        # Text
        jname = ltu.name_from_coord(secondary_coord)
        ax.text(0.5,1.24, jname, fontsize=22, color='blue',
                horizontalalignment='center',transform=ax.transAxes)
        # Print offsets
        if primary_coord is not None:
            sep = primary_coord.separation(secondary_coord).to('arcsec')
            PA = primary_coord.position_angle(secondary_coord)
            # RA/DEC
            dec_off = np.cos(PA) * sep # arcsec
            ra_off = np.sin(PA) * sep # arcsec (East is *higher* RA)
            ax.text(0.5, 1.22, 'Offset from Ref. Star (cyan) to Target (red):\nRA(to targ) = {:.2f}  DEC(to targ) = {:.2f}'.format(
                -1*ra_off.to('arcsec'), -1*dec_off.to('arcsec')),
                     fontsize=15, horizontalalignment='center',transform=ax.transAxes, color='blue', va='top')
    # Add tertiary
    if third_coord is not None:
        c = SphericalCircle((third_coord.ra, third_coord.dec),
                            2*units.arcsec, transform=ax.get_transform('icrs'),
                            edgecolor='yellow', facecolor='none')
        ax.add_patch(c)
    
    # Slit
    if ((slit is not None) and (flag_photu is True)):
        # List of values - [coordinates, width, length, PA],
        # e.g. [SkyCoords('21h44m25.255s',-40d54m00.1s', frame='icrs'), 1*u.arcsec, 10*u.arcsec, 20*u.deg]
        slit_coords, width, length, pa = slit
        
        if flip_ra:  # NT: not sure this is needed...
            pa = -1*pa
        
        pa_deg = pa.to('deg').value
        # according to SkyRectangularAperture the theta is for the "width", thus we should conform that geometry here.
        # for PA = 0 we want the slit to be oriented North-South, thus we need to rotate by +90 degrees and apply - sign 
        aper = SkyRectangularAperture(positions=slit_coords, w=width, h=length, theta=-1*(pa+90*units.deg))  # PA should increase from North to East 
        
        apermap = aper.to_pixel(wcs)
        
        apermap.plot(color='purple', lw=1)
        
        plt.text(0.5, -0.15, 'Slit PA={} deg'.format(pa_deg), color='purple',
                 fontsize=15, ha='center', va='top', transform=ax.transAxes)
    
    if ((slit is not None) and (flag_photu is False)):
        raise IOError('Slit cannot be placed without photutils package')
    
    # Title
    ax.text(0.5, 1.44, title, fontsize=32, horizontalalignment='center', transform=ax.transAxes)

    # Extra text
    if extra_text is not None:
        ax.text(-0.1, -0.3, extra_text, fontsize=20, horizontalalignment='left', transform=ax.transAxes)
    # Sources

    # Labels
    #ax.set_xlabel(r'\textbf{DEC (EAST direction)}')
    #ax.set_ylabel(r'\textbf{RA (SOUTH direction)}')

    if outfile is not None:
        plt.savefig(outfile, dpi=400)
        plt.close()
    else:
        plt.show()

    # Return
    return fig, ax


#### ###############################
def sdss_dss(coord, title, show_circ=True, EPOCH=None, imsize=5.*units.arcmin, outfile=None,
         show_slit=None, show_another=None, cradius=None, in_img=None, dss_only=False,
         vmnx=(None,None)):
    """
    Generate a SDSS or DSS finder

    Pulled from xastropy

    Args:
        coord (SkyCoord):
        title (str):
        show_circ (bool, optional):
            Show a yellow circle on the target
        EPOCH (float, optional):
            Defaults to 2000.
        imsize (Quantity or Angle):
            Size of the finder
        show_slit (list, optional):
            Show a slit with [width, length, PA]
        show_another (bool, optional):
            Not yet functional
        cradius (Quantity, optional):
            Circle radius, only shown if show_circ is True.
            Default is imsize/50.
        in_img (PIL image, optional):
        dss_only (bool, optional):
            Only pull from DSS
        vmnx (tuple, optional):
            vmin, vmax
    """
    # Init
    vimsize = imsize.to('arcmin').value
    if cradius is None:
        cradius = vimsize / 50.
    else:
        cradius = cradius.to('arcmin').value

    # Precess (as necessary)
    if EPOCH is not None:
        # Precess to 2000.
        tEPOCH = Time(EPOCH, format='jyear', scale='utc')
        # Load into astropy
        icrs = ICRS(ra=coord.ra.value, dec=coord.dec.value, equinox=tEPOCH, unit='deg')
        # Precess
        newEPOCH = Time(2000., format='jyear', scale='utc')
        coord = icrs.precess_to(newEPOCH)

    # Grab the Image
    if in_img is None:
        img, oBW = getjpg(coord, imsize, dss_only=dss_only)
    else:
        img = in_img
        oBW = True

    # Generate the plot
    plt.clf()
    fig = plt.figure(dpi=700)
    fig.set_size_inches(8.0, 10.5)

    # Font
    plt.rcParams['font.family'] = 'times new roman'
    ticks_font = font_manager.FontProperties(family='times new roman', style='normal',
                                             size=16, weight='normal', stretch='normal')
    ax = plt.gca()
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)

    # Image
    if oBW == 1:
        cmm = cm.Greys
    else:
        cmm = None
    plt.imshow(img, cmap=cmm, aspect='equal', extent=(-vimsize / 2., vimsize / 2,
                                                      -vimsize / 2., vimsize / 2),
               vmin=vmnx[0], vmax=vmnx[1])

    # Axes
    plt.xlim(-vimsize / 2., vimsize / 2.)
    plt.ylim(-vimsize / 2., vimsize / 2.)

    # Label
    plt.xlabel('Relative ArcMin', fontsize='large')
    xpos = 0.12 * vimsize
    ypos = 0.02 * vimsize
    plt.text(-vimsize / 2. - xpos, 0., 'EAST', rotation=90., fontsize=20)
    plt.text(0., vimsize / 2. + ypos, 'NORTH', fontsize=20, horizontalalignment='center')

    # import pdb; pdb.set_trace()

    # Circle
    if show_circ:
        circle = plt.Circle((0, 0), cradius, color='y', fill=False)
        plt.gca().add_artist(circle)

    # Second Circle
    if show_another is not None:
        raise NotImplementedError
        # Coordinates
        canother = coord.to_coord(show_another)
        embed(header='338')
        # Offsets
        xanother = -1 * off[0].to('arcmin').value
        yanother = off[1].to('arcmin').value
        square = matplotlib.patches.Rectangle((xanother - cradius,
                                               yanother - cradius),
                                              cradius * 2, cradius * 2, color='cyan', fill=False)
        plt.gca().add_artist(square)
        plt.text(0.5, 1.24, str(nm), fontsize=32,
                 horizontalalignment='center', transform=ax.transAxes)
        plt.text(0.5, 1.18, 'RA (J2000) = ' + str(obj['RAS']) +
                 '  DEC (J2000) = ' + str(obj['DECS']), fontsize=22,
                 horizontalalignment='center', transform=ax.transAxes)
        plt.text(0.5, 1.12, 'RA(other) = {:s}  DEC(other) = {:s}'.format(
            canother.ra.to_string(unit=astrou.hour, pad=True, sep=':', precision=2),
            canother.dec.to_string(pad=True, alwayssign=True, sep=':', precision=1)),
                 fontsize=22, horizontalalignment='center', transform=ax.transAxes,
                 color='blue')
        plt.text(0.5, 1.06, 'RA(to targ) = {:.2f}  DEC(to targ) = {:.2f} PA={:g}'.format(
            -1 * off[0].to('arcsec'), -1 * off[1].to('arcsec'), PA),
                 fontsize=18, horizontalalignment='center', transform=ax.transAxes)
    else:
        # Title
        plt.text(0.5, 1.24, str(title), fontsize=32,
                 horizontalalignment='center', transform=ax.transAxes)
        ras = coord.ra.to_string(unit=units.hour,sep=':',pad=True,precision=2)
        decs = coord.dec.to_string(sep=':',pad=True, alwayssign=True,precision=1)
        plt.text(0.5, 1.16, 'RA (J2000) = ' + ras, fontsize=28,
                 horizontalalignment='center', transform=ax.transAxes)
        plt.text(0.5, 1.10, 'DEC (J2000) = ' + decs, fontsize=28,
                 horizontalalignment='center', transform=ax.transAxes)

    # Show slit??
    if show_slit is not None:
        # List of values - [width, length, PA],
        # e.g. [1*u.arcsec, 10*u.arcsec, 20*u.deg]
        w, l, pa = show_slit
        w_arcmin = w.to('arcmin').value
        l_arcmin = l.to('arcmin').value
        pa_deg = pa.to('deg').value
        pa_rad = pa_deg * np.pi / 180.
        # get the new position of the lower-left corner of rectangle given the PA
        y_new = 0. - 0.5 * (l_arcmin * np.sin(np.pi / 2. - pa_rad) + w_arcmin * np.sin(pa_rad))
        x_new = 0. + 0.5 * (l_arcmin * np.cos(np.pi / 2. - pa_rad) - w_arcmin * np.cos(pa_rad))
        xy = (x_new, y_new)  # xy of lower-left corner (after rotation)
        box = plt.Rectangle(xy, w_arcmin, l_arcmin, color='k', angle=pa_deg, fill=False, lw=0.5)
        plt.gca().add_artist(box)
        plt.text(0.5, 0.05, 'Slit PA={}deg'.format(pa_deg),
                 fontsize=15, ha='center', va='top', transform=ax.transAxes)

    if outfile is not None:
        print("Writing: {}".format(outfile))
        plt.savefig(outfile)
        plt.close()
    else:
        plt.show()


# ##########################################
def getjpg(coord, imsize, dss_only=False):
    """
    Grab an SDSS or DSS image

    Args:
        coord (SkyCoord):
        imsize (Angle or Quantity):
            image size
        dss_only (bool, optional):
            Only pull from DSS

    Returns:
        PIL, bool:  Image, flag indicating if the image is B&W
    """
    sdss = survey_utils.load_survey_by_name('SDSS', coord, 0.02*units.deg)
    # Dummy call to see if SDSS covers it
    if dss_only:
        cat = None
    else:
        cat = sdss.get_catalog()
    # Pull from DSS?
    if cat is None:
        print("No SDSS Image;  Querying DSS")
        BW = True
        url = dsshttp(coord, imsize)  # DSS
        img = images.grab_from_url(url)
    else:
        BW = False
        img, _ = sdss.get_cutout(imsize)
    # Return
    return img, BW


def dsshttp(coord, imsize):
    """
    Generate URL for DSS

    Args:
        coord (SkyCoord):
        imsize (Angle or Quantity):
            image size

    Returns:
        str:  URL

    """
    # https://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r=00:42:44.35&d=+41:16:08.6&e=J2000&h=15.0&w=15.0&f=gif&c=none&fov=NONE&v3=

    Equinox = 'J2000'
    dss = 'poss2ukstu_red'
    url = "http://archive.stsci.edu/cgi-bin/dss_search?"
    url += "v=" + dss + '&r=' + str(coord.ra.value) + '&d=' + str(coord.dec.value)
    url += "&e=" + Equinox
    url += '&h=' + str(imsize) + "&w=" + str(imsize)
    url += "&f=gif"
    url += "&c=none"
    url += "&fov=NONE"
    url += "&v3="

    return url
