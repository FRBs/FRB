"""
Slurp data from 2MASS catalog.

"""

import numpy as np
import warnings

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io.fits import PrimaryHDU
from astropy.table import Table
from ..galaxies.defs import MASS_bands
from astroquery.ipac.irsa import Irsa

from frb.surveys import surveycoord,catalog_utils
from frb.surveys.skyview import SkyView_Survey


# Define the data model for 2MASS data
photom = {}
photom['2MASS'] = {}
for band in MASS_bands:
    photom["2MASS"]["2MASS"+'_{:s}'.format(band)] = '{:s}_m'.format(band.lower()) # Many options for apertures, tbd
    photom["2MASS"]["2MASS"+'_{:s}_err'.format(band)] = '{:s}_msig'.format(band.lower())
    photom["2MASS"]["2MASS_ID"] = 'designation'
photom["2MASS"]['ra'] = 'ra'
photom["2MASS"]['dec'] = 'dec'

# Define the default set of query fields
# http://tdc-www.harvard.edu/catalogs/tmpsc.format.html For PSC
# https://www.ipac.caltech.edu/2mass/releases/second/doc/ancillary/xscformat.html For XSC
_DEFAULT_query_fields = ['designation','survey','ra','dec']
_DEFAULT_query_fields +=['{:s}_m'.format(band) for band in MASS_bands]
_DEFAULT_query_fields +=['{:s}_msig'.format(band) for band in MASS_bands]

class TwoMASS_Survey(SkyView_Survey):
    """
    A class to access the 2MASS catalogs hosted on the
    IRSA database and 2MASS images via SkyView.
    Inherits from SkyView_Survey.

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """
    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        SkyView_Survey.__init__(self, coord, radius, '2mass', **kwargs)

        self.Survey = "2MASS"
        self.survey = '2MASS'
    
    def get_catalog(self, query_fields: list[str] | None = None) -> Table:
        """
        Query a catalog in the IRSA 2MASS survey for
        photometry.

        The extended source catalog is searched first; the point source
        catalog is used only if the former has no sources in the cone.
        Magnitudes are converted to AB.

        Args:
            query_fields (list of str, optional): A list of query fields to
                get in addition to the default fields.

        Returns:
            astropy.table.Table: Contains all query results
        """

        if query_fields is None:
            query_fields = _DEFAULT_query_fields
        else:
            query_fields = _DEFAULT_query_fields+query_fields
        
        data = {}
        data['ra'] = self.coord.ra.value
        data['dec'] = self.coord.dec.value
        data['radius'] = self.radius.to(u.deg).value
        data['columns'] = query_fields
        data['format'] = 'csv'

        # First query the extended source catalog
        # Fields described here: http://tdc-www.harvard.edu/catalogs/tmx.format.html
        ret = Irsa.query_region(self.coord, radius=self.radius, spatial='Cone',
                                catalog="fp_xsc")
        isempty = len(ret) == 0

        if isempty:
            # If fp_xsc is empty, query the psc catalog
            ret = Irsa.query_region(self.coord, radius=self.radius, spatial='Cone',
                                    catalog="fp_psc")
            for band in MASS_bands: # Rename columns for mags for PSC
                photom["2MASS"]["2MASS"+'_{:s}'.format(band)] = '{:s}_m'.format(band.lower())
                photom["2MASS"]["2MASS"+'_{:s}_err'.format(band)] = '{:s}_msigcom'.format(band.lower())
        else: # if XSC is not empty, rename columns for mags for XSC
            # Instead of _m and _msig, it's _m_fe and _msig_fe for fiducial elliptical Kron
            for band in MASS_bands:
                photom["2MASS"]["2MASS"+'_{:s}'.format(band)] = '{:s}_m_fe'.format(band.lower())
                photom["2MASS"]["2MASS"+'_{:s}_err'.format(band)] = '{:s}_msig_fe'.format(band.lower())

        pdict = photom['2MASS'].copy()
        
        photom_catalog = catalog_utils.clean_cat(ret, pdict, mask_photometry=True) # rename columns

        photom_catalog.keep_columns(list(pdict.keys())) # Keep only the columns we care about

        # Remove duplicate entries.
        photom_catalog = catalog_utils.remove_duplicates(photom_catalog, "2MASS_ID")

        self.catalog = catalog_utils.sort_by_separation(photom_catalog, self.coord,
                                                        radec=('ra','dec'), add_sep=True)
        
        self.convert_to_AB()

        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        #Validate
        self.validate_catalog()

        #Return
        return self.catalog.copy()
        
    def convert_to_AB(self) -> Table:
        """
        Convert from 2MASS internal to AB magnitudes in the catalog.

        ``self.catalog`` is modified in place.

        Returns:
            astropy.table.Table: ``self.catalog`` with AB magnitudes.

        Raises:
            ValueError: If a 2MASS magnitude column is missing from the catalog.

        """
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

    def get_image(self, imsize: u.Quantity, band: str = 'J') -> PrimaryHDU | None:
        """
        Retrieve a SkyView FITS image for 2MASS.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            band (str, optional): 'J', 'H' or 'K'.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        """
        return SkyView_Survey.get_image(self, imsize=imsize, band=band)

    def get_cutout(self, imsize: u.Quantity, band: str = 'J') -> PrimaryHDU | None:
        """
        Deprecated alias for FITS image retrieval.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            band (str, optional): 'J', 'H' or 'K'.

        Returns:
            astropy.io.fits.PrimaryHDU or None: See :meth:`get_image`.

        """
        warnings.warn(
            "get_cutout() returns FITS products for this survey and is deprecated; use get_image() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_image(imsize=imsize, band=band)