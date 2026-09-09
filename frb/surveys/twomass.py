"""
Slurp data from 2MASS catalog.

"""

import gzip
import io
import warnings
import numpy as np

import requests
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from ..galaxies.defs import MASS_bands
from astroquery.ipac.irsa import Irsa

from frb.surveys import surveycoord, catalog_utils

try:
    from photutils.aperture import (
        CircularAperture, CircularAnnulus, aperture_photometry
    )
    _HAS_PHOTUTILS = True
    try:
        from photutils.aperture import ApertureStats
        _HAS_APERTURE_STATS = True
    except ImportError:
        _HAS_APERTURE_STATS = False
except ImportError:
    _HAS_PHOTUTILS = False
    _HAS_APERTURE_STATS = False

try:
    from reproject import reproject_interp
    from reproject.mosaicking import reproject_and_coadd
    _HAS_REPROJECT = True
except ImportError:
    _HAS_REPROJECT = False


# IRSA IBE root URL for 2MASS All-Sky image retrieval
_IRSA_IBE_ROOT = "https://irsa.ipac.caltech.edu/ibe"

# Vega zero-point fluxes in Jy; used for Vega-to-AB conversion
_TWOMASS_VEGAFLUX = {'j': 1594.0, 'h': 1024.0, 'k': 666.7}

# Nominal Vega magnitude zero points used when MAGZP is absent from the header
_TWOMASS_DEFAULT_MAGZP = {'j': 20.9, 'h': 20.4, 'k': 19.9}

# Metadata columns needed to build the data URL, rank tiles by how well they
# cover the target, and harmonize zero points before mosaicking.
_IBE_COLUMNS = ("ordate,hemisphere,scanno,fname,ra,dec,"
                "ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,magzp")

# 2MASS All-Sky Atlas images are 1"/pixel
_TWOMASS_PIXSCALE = 1.0


def _tile_margin(row, ra_deg, dec_deg):
    """
    Half-width of the largest box centered on the target that fits in a tile.

    2MASS All-Sky Atlas images are 512x1024 pixels (8.5' x 17') and overlap, so
    most positions are covered by several tiles. IRSA returns them in an order
    that has nothing to do with where the target falls, and a tile whose seam
    runs through the target yields a cutout clipped at the array boundary. This
    ranks candidate tiles so the best-centered one can be chosen.

    Args:
        row (astropy.table.Row): IBE metadata row with ra1..ra4/dec1..dec4.
        ra_deg (float): Target RA in degrees.
        dec_deg (float): Target Dec in degrees.

    Returns:
        float: Distance in arcsec from the target to the nearest tile edge, or
        -1. if the target lies outside the tile footprint.
    """
    cosd = np.cos(np.radians(dec_deg))
    dx = np.array([float(row['ra{:d}'.format(i)]) - ra_deg for i in (1, 2, 3, 4)])
    dx = ((dx + 180.) % 360. - 180.) * cosd * 3600.
    dy = np.array([float(row['dec{:d}'.format(i)]) - dec_deg
                   for i in (1, 2, 3, 4)]) * 3600.

    # crota2 is ~0.002-0.009 deg for these tiles, so treating the footprint as
    # an axis-aligned box is accurate to well under a pixel.
    if not (dx.min() <= 0 <= dx.max() and dy.min() <= 0 <= dy.max()):
        return -1.
    return float(min(-dx.min(), dx.max(), -dy.min(), dy.max()))


def _select_tiles(tab, band, ra_deg, dec_deg):
    """
    Rank the atlas tiles of one band, best-centered first.

    Tiles that do not contain the target (margin < 0) are kept and sorted last:
    they cannot serve as the single-tile cutout, but they are exactly what fills
    out the corners of a mosaic, so dropping them here would leave NaN gaps.
    Callers wanting only tiles that contain the target should filter on
    ``margin >= 0``.

    Args:
        tab (astropy.table.Table): IBE metadata rows.
        band (str): 2MASS band -- 'j', 'h' or 'k'.
        ra_deg (float): Target RA in degrees.
        dec_deg (float): Target Dec in degrees.

    Returns:
        list: (row, margin) pairs for the requested band, sorted by descending
        margin.
    """
    scored = []
    for row in tab:
        if not str(row['fname']).strip().lower().startswith(band):
            continue
        scored.append((row, _tile_margin(row, ra_deg, dec_deg)))
    scored.sort(key=lambda rm: rm[1], reverse=True)
    return scored


def _tile_url(row):
    """Build the IBE data URL for one 2MASS atlas tile."""
    return (
        "{root}/data/twomass/allsky/allsky/"
        "{ordate:06d}{hemi}/s{scanno:03d}/image/{fname}"
    ).format(
        root=_IRSA_IBE_ROOT,
        ordate=int(row['ordate']),
        hemi=str(row['hemisphere']).strip(),
        scanno=int(row['scanno']),
        fname=str(row['fname']).strip(),
    )


def _fetch_fits(url, timeout=180):
    """GET a (possibly gzipped) FITS file and return its PrimaryHDU, or None."""
    r = requests.get(url, timeout=timeout)
    if not r.ok:
        warnings.warn("2MASS download failed (HTTP {:d}): {:s}".format(
            r.status_code, url))
        return None
    raw = r.content
    # Decompress if the server returned gzip despite a gzip=false request
    if len(raw) >= 2 and raw[0] == 0x1F and raw[1] == 0x8B:
        raw = gzip.decompress(raw)
    return fits.open(io.BytesIO(raw))[0]


def _mosaic_tiles(scored, ra_deg, dec_deg, size_arcsec):
    """
    Mosaic overlapping 2MASS atlas tiles onto a grid centered on the target.

    Used when no single tile can supply the requested cutout centered on the
    target. Tiles are rescaled to a common zero point before coadding so the
    returned MAGZP stays meaningful; the measured spread between neighbouring
    tiles is <0.001 mag, so this is a sub-0.1% correction.

    Note that full tiles are downloaded rather than server-side cutouts: IBE's
    cutout-on-URL service returns HTTP 500 when the requested center falls
    outside the tile, which is common for tiles that merely overlap the box.

    Args:
        scored (list): (row, margin) pairs from _select_tiles.
        ra_deg (float): Target RA in degrees.
        dec_deg (float): Target Dec in degrees.
        size_arcsec (float): Requested cutout size in arcsec.

    Returns:
        fits.PrimaryHDU or None: Mosaicked cutout centered on the target.
    """
    from astropy.wcs import WCS

    ref_row = scored[0][0]
    zp_ref = float(ref_row['magzp'])

    inputs = []
    names = []
    ref_hdr = None
    for row, _ in scored:
        hdu = _fetch_fits(_tile_url(row))
        if hdu is None or hdu.data is None:
            continue
        # Put every tile on the reference tile's zero point
        scale = 10 ** (-0.4 * (zp_ref - float(row['magzp'])))
        inputs.append((np.array(hdu.data, dtype=float) * scale, WCS(hdu.header)))
        names.append(str(row['fname']).strip())
        if str(row['fname']).strip() == str(ref_row['fname']).strip():
            ref_hdr = hdu.header

    if not inputs:
        warnings.warn("No 2MASS tiles could be downloaded for the mosaic.")
        return None

    npix = int(round(size_arcsec / _TWOMASS_PIXSCALE))
    out_wcs = WCS(naxis=2)
    out_wcs.wcs.crpix = [npix / 2. + 0.5, npix / 2. + 0.5]
    out_wcs.wcs.crval = [ra_deg, dec_deg]
    out_wcs.wcs.cdelt = [-_TWOMASS_PIXSCALE / 3600., _TWOMASS_PIXSCALE / 3600.]
    out_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    # match_background removes the per-tile sky offsets; the resulting absolute
    # sky level is arbitrary, which is fine since photometry below subtracts a
    # local annulus background.
    mosaic, _ = reproject_and_coadd(
        inputs, out_wcs, shape_out=(npix, npix),
        reproject_function=reproject_interp,
        match_background=True, combine_function='mean',
    )

    hdr = out_wcs.to_header()
    hdr['MAGZP'] = zp_ref
    if ref_hdr is not None:
        for key in ('FILTER', 'ORIGIN', 'SKYSIG'):
            if key in ref_hdr:
                hdr[key] = ref_hdr[key]
    hdr['HISTORY'] = 'Mosaic of 2MASS tiles: ' + ', '.join(names)
    return fits.PrimaryHDU(data=mosaic, header=hdr)


# Define the data model for 2MASS data
# XSC magnitudes use the k20fe aperture column names (e.g. j_m_k20fe, j_msig_k20fe).
# PSC magnitudes use the simpler j_m / j_msigcom names; those are overridden below in
# get_catalog when the PSC fallback is triggered.
photom = {}
photom['2MASS'] = {}
for band in MASS_bands:
    photom["2MASS"]["2MASS"+'_{:s}'.format(band)] = '{:s}_m_k20fe'.format(band)
    photom["2MASS"]["2MASS"+'_{:s}_err'.format(band)] = '{:s}_msig_k20fe'.format(band)
    photom["2MASS"]["2MASS_ID"] = 'designation'
photom["2MASS"]['ra'] = 'ra'
photom["2MASS"]['dec'] = 'dec'

# Define the default set of query fields
# http://tdc-www.harvard.edu/catalogs/tmpsc.format.html For PSC
# https://www.ipac.caltech.edu/2mass/releases/second/doc/ancillary/xscformat.html For XSC
_DEFAULT_query_fields = ['designation','survey','ra','dec']
_DEFAULT_query_fields +=['{:s}_m_k20fe'.format(band) for band in MASS_bands]
_DEFAULT_query_fields +=['{:s}_msig_k20fe'.format(band) for band in MASS_bands]

class TwoMASS_Survey(surveycoord.SurveyCoord):
    """
    A class to access all the catalogs hosted on the
    IRSA database. Inherits from SurveyCoord. This
    is a super class not meant for use by itself and
    instead meant to instantiate specific children
    classes like TwoMASS_Survey
    """
    def __init__(self,coord,radius,**kwargs):
        surveycoord.SurveyCoord.__init__(self,coord,radius,**kwargs)

        self.Survey = "2MASS"

    def get_catalog(self, query_fields=None, use_image_photom=True,
                    imsize=1*u.arcmin, aperture_radius=8*u.arcsec,
                    sky_inner=12*u.arcsec, sky_outer=14*u.arcsec):
        """
        Query a catalog in the IRSA 2MASS survey for
        photometry.


        Args:
            query_fields: list, optional
                A list of query fields to
                get in addition to the
                default fields.
            use_image_photom: bool, optional
                If True and the source is only found in the PSC (not XSC),
                download 2MASS images via IRSA IBE and perform aperture
                photometry to replace the unreliable PSC catalog magnitudes.
                Requires photutils. Default True.
            imsize: Quantity, optional
                Angular image size for cutout download when
                use_image_photom is True. Default 1 arcmin.
            aperture_radius: Quantity, optional
                Source aperture radius for image photometry.
                Default 5 arcsec.
            sky_inner: Quantity, optional
                Inner radius of the sky background annulus.
                Default 8 arcsec.
            sky_outer: Quantity, optional
                Outer radius of the sky background annulus.
                Default 12 arcsec.


        Returns:
            catalog: astropy.table.Table
                Contains all query results
        """

        if query_fields is None:
            query_fields = _DEFAULT_query_fields
        else:
            query_fields = _DEFAULT_query_fields+query_fields

        # First query the extended source catalog
        # Fields described here: http://tdc-www.harvard.edu/catalogs/tmx.format.html
        ret = Irsa.query_region(self.coord, radius=self.radius, spatial='Cone',
                                catalog="fp_xsc")
        psc_only = len(ret) == 0

        if psc_only:
            print("\tNo sources found in 2MASS Extended Source Catalog. " \
            "Querying Point Source Catalog instead.")
            # If fp_xsc is empty, query the psc catalog
            ret = Irsa.query_region(self.coord, radius=self.radius, spatial='Cone',
                                    catalog="fp_psc")
            for band in MASS_bands: # Rename columns for mags for PSC
                photom["2MASS"]["2MASS"+'_{:s}'.format(band)] = '{:s}_m'.format(band.lower())
                photom["2MASS"]["2MASS"+'_{:s}_err'.format(band)] = '{:s}_msigcom'.format(band.lower())

        pdict = photom['2MASS'].copy()

        photom_catalog = catalog_utils.clean_cat(ret,pdict) # rename columns

        photom_catalog.keep_columns(list(pdict.keys())) # Keep only the columns we care about

        # Remove duplicate entries.
        photom_catalog = catalog_utils.remove_duplicates(photom_catalog, "2MASS_ID")

        self.catalog = catalog_utils.sort_by_separation(photom_catalog, self.coord,
                                                        radec=('ra','dec'), add_sep=True)

        self.convert_to_AB()

        # PSC magnitudes are aperture-based and unreliable for extended/galaxy sources.
        # Replace the nearest match's mags with aperture photometry on 2MASS image atlas.
        if psc_only and use_image_photom and len(self.catalog) > 0:
            if not _HAS_PHOTUTILS:
                warnings.warn(
                    "photutils is required for 2MASS image photometry. "
                    "Falling back to unreliable PSC catalog magnitudes."
                )
            else:
                print("\tPSC magnitudes are unreliable for galaxies; "
                      "performing aperture photometry on 2MASS images.")
                img_mags = self.get_image_photom(
                    imsize=imsize, aperture_radius=aperture_radius,
                    sky_inner=sky_inner, sky_outer=sky_outer,
                )
                # Overwrite the nearest source (row 0) with image-derived AB mags.
                # Other rows retain (unreliable) PSC mags.
                for band in MASS_bands:
                    mag_col = '2MASS_{:s}'.format(band)
                    err_col = '2MASS_{:s}_err'.format(band)
                    if mag_col in self.catalog.colnames:
                        self.catalog[mag_col][0] = img_mags.get(mag_col, np.nan)
                    if err_col in self.catalog.colnames:
                        self.catalog[err_col][0] = img_mags.get(err_col, np.nan)
                # Pin the position to self.coord so cross-matching with other surveys
                # succeeds in search_all_surveys. The image photometry is measured at
                # self.coord, not at the (potentially offset) PSC catalog position.
                self.catalog['ra'][0]  = self.coord.ra.deg
                self.catalog['dec'][0] = self.coord.dec.deg

        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        #Validate
        self.validate_catalog()

        #Return
        return self.catalog.copy()

    def get_cutout(self, imsize, band='j'):
        """
        Download a single-band 2MASS FITS image cutout via the IRSA IBE API.

        Queries the 2MASS All-Sky metadata table for every atlas tile covering
        self.coord and picks the one on which the target sits furthest from any
        edge, then downloads a centered FITS cutout at the requested angular
        size. The returned HDU preserves the original FITS header including the
        MAGZP calibration keyword.

        Atlas tiles are only 8.5' x 17', so a large enough request cannot be
        satisfied by any single tile centered on the target. In that case the
        overlapping tiles are mosaicked instead (requires the reproject
        package); without reproject the best single tile is returned and the
        cutout is trimmed on the short side.

        Args:
            imsize (Quantity): Angular size of the cutout.
            band (str): 2MASS band — 'j', 'h', or 'k'.

        Returns:
            fits.PrimaryHDU or None: FITS HDU with image data and header,
            or None if no image was found at this position.
        """
        band = band.lower()
        if band not in ('j', 'h', 'k'):
            raise ValueError("band must be one of 'j', 'h', 'k'")

        ra_deg  = self.coord.ra.deg
        dec_deg = self.coord.dec.deg
        size_arcsec = imsize.to(u.arcsec).value

        # Query IRSA IBE metadata table: get atlas tile info for this position.
        # SIZE is the requested box rather than 0.0 so that tiles which overlap
        # the cutout without containing its center are available to the mosaic.
        search_url = "{:s}/search/twomass/allsky/allsky".format(_IRSA_IBE_ROOT)
        params = {
            "POS":     "{:.6f},{:.6f}".format(ra_deg, dec_deg),
            "SIZE":    "{:.6f}".format(size_arcsec / 3600.),
            "ct":      "csv",
            "columns": _IBE_COLUMNS,
        }
        r = requests.get(search_url, params=params, timeout=120)
        if not r.ok:
            warnings.warn(
                "IRSA IBE metadata query failed for 2MASS {:s} band "
                "(HTTP {:d}).".format(band, r.status_code)
            )
            return None

        tab = Table.read(io.BytesIO(r.content), format="ascii.csv")
        if len(tab) == 0:
            warnings.warn(
                "No 2MASS atlas tile found at this position for {:s} band.".format(band)
            )
            return None

        # Rank this band's tiles by how far the target sits from their edges.
        # Tiles that merely overlap the box (margin < 0) are retained for the
        # mosaic; only those containing the target can serve a single cutout.
        scored = _select_tiles(tab, band, ra_deg, dec_deg)
        covering = [(row, margin) for row, margin in scored if margin >= 0]
        if not covering:
            warnings.warn(
                "No 2MASS {:s}-band image found at this position.".format(band)
            )
            return None

        chosen, margin = covering[0]

        if margin < size_arcsec / 2.:
            # No single tile can supply the request centered on the target
            if _HAS_REPROJECT:
                hdu = _mosaic_tiles(scored, ra_deg, dec_deg, size_arcsec)
                if hdu is not None:
                    self.cutout      = hdu
                    self.cutout_size = imsize
                    return hdu
            else:
                warnings.warn(
                    "reproject is not installed; the 2MASS {:s}-band cutout "
                    "will be trimmed to the best single atlas tile.".format(band)
                )

        # Construct the IBE data URL with on-the-fly cutout parameters
        data_url = _tile_url(chosen) + (
            "?center={ra:.6f},{dec:.6f}&size={size:.2f}arcsec&gzip=false"
        ).format(ra=ra_deg, dec=dec_deg, size=size_arcsec)

        hdu = _fetch_fits(data_url, timeout=120)
        if hdu is None:
            return None

        self.cutout      = hdu
        self.cutout_size = imsize
        return hdu

    def get_image_photom(self, imsize=1*u.arcmin, aperture_radius=8*u.arcsec,
                         sky_inner=12*u.arcsec, sky_outer=14*u.arcsec,
                         bands=None):
        """
        Download 2MASS images and perform circular aperture photometry
        at self.coord in J, H, and K bands.

        Sky background is estimated from the pixel median in a circular
        annulus and subtracted before the source flux is summed.
        Calibration uses the MAGZP keyword from the FITS header (Vega
        system); a nominal fallback zero point is used when absent.
        Returned magnitudes are in the AB system.

        Args:
            imsize (Quantity): Angular size of the image to retrieve.
                Default 1 arcmin.
            aperture_radius (Quantity): Radius of the source aperture.
                Default 5 arcsec.
            sky_inner (Quantity): Inner radius of the background annulus.
                Default 8 arcsec.
            sky_outer (Quantity): Outer radius of the background annulus.
                Default 12 arcsec.
            bands (list, optional): Subset of bands to process.
                Defaults to all three: ['j', 'h', 'k'].

        Returns:
            dict: Keys are '2MASS_j', '2MASS_h', '2MASS_k' and their
                '_err' counterparts. Values are AB magnitudes (float);
                999. where photometry could not be performed.
        """
        if not _HAS_PHOTUTILS:
            raise ImportError(
                "photutils is required for 2MASS image photometry."
            )

        from astropy.wcs import WCS
        from astropy.wcs.utils import proj_plane_pixel_scales

        if bands is None:
            bands = MASS_bands

        results = {}

        for band in bands:
            hdu = self.get_cutout(imsize, band)

            if hdu is None or hdu.data is None:
                results['2MASS_{:s}'.format(band)]     = 999.
                results['2MASS_{:s}_err'.format(band)] = 999.
                continue

            wcs  = WCS(hdu.header)
            data = np.array(hdu.data, dtype=float)

            # Pixel coordinates of the target
            x, y = wcs.world_to_pixel(self.coord)
            x, y = float(x), float(y)

            # Verify the target falls within the image footprint
            ny, nx = data.shape
            if not (0 <= x < nx and 0 <= y < ny):
                warnings.warn(
                    "Target falls outside the 2MASS {:s}-band image.".format(band)
                )
                results['2MASS_{:s}'.format(band)]     = 999.
                results['2MASS_{:s}_err'.format(band)] = 999.
                continue

            # Pixel scale (arcsec/pixel), averaged over both axes
            pscale = np.mean(proj_plane_pixel_scales(wcs)) * 3600.0

            ap_r   = aperture_radius.to(u.arcsec).value / pscale
            sky_ri = sky_inner.to(u.arcsec).value / pscale
            sky_ro = sky_outer.to(u.arcsec).value / pscale

            aperture = CircularAperture((x, y), r=ap_r)
            annulus  = CircularAnnulus((x, y), r_in=sky_ri, r_out=sky_ro)

            # Sky background: median and rms from the annulus pixels
            if _HAS_APERTURE_STATS:
                sky_stats  = ApertureStats(data, annulus)
                bkg_median = float(sky_stats.median)
                bkg_rms    = float(sky_stats.std)
            else:
                ann_mask   = annulus.to_mask(method='center')
                sky_pix    = ann_mask.multiply(data)
                sky_pix    = sky_pix[sky_pix != 0]
                bkg_median = float(np.nanmedian(sky_pix))
                bkg_rms    = float(np.nanstd(sky_pix))

            # Background-subtracted aperture sum
            data_sub = data - bkg_median
            phot_tab = aperture_photometry(data_sub, aperture)
            total_dn = float(phot_tab['aperture_sum'][0])

            # Noise: sky fluctuations propagated over the aperture area
            noise_dn = bkg_rms * np.sqrt(aperture.area)

            # Magnitude calibration via header MAGZP (Vega); fall back to nominal
            magzp   = hdu.header.get('MAGZP', None)
            exptime = float(hdu.header.get('EXPTIME', 1.0))
            if magzp is None:
                warnings.warn(
                    "MAGZP not found in 2MASS {:s}-band image header; "
                    "using nominal zero point.".format(band)
                )
                magzp = _TWOMASS_DEFAULT_MAGZP[band]
            magzp = float(magzp)

            if total_dn > 0 and noise_dn > 0:
                # Vega magnitude from calibrated pixel sum
                m_vega = magzp - 2.5 * np.log10(total_dn / exptime)
                m_err  = (2.5 / np.log(10)) * (noise_dn / total_dn)
                # Convert Vega -> AB
                m_ab = m_vega - 2.5 * np.log10(_TWOMASS_VEGAFLUX[band] / 3630.7805)
            else:
                m_ab  = 999.
                m_err = 999.

            results['2MASS_{:s}'.format(band)]     = m_ab
            results['2MASS_{:s}_err'.format(band)] = m_err

        return results

    def convert_to_AB(self):
        """Convert from 2MASS internal to AB magnitudes in the catalog."""
        # Convert to AB mag
        fnu0 = {'2MASS_j':1594,
                '2MASS_h':1024,
                '2MASS_k':666.7}

        for band in MASS_bands:
            filt = '2MASS_{:s}'.format(band.lower())
            if filt in self.catalog.columns:
                self.catalog[filt] -= 2.5*np.log10(fnu0[filt]/3630.7805)
            else:
                raise ValueError(f"Column {filt} not found in catalog.")

        return self.catalog
