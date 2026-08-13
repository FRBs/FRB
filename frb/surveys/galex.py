"""
Slurp data from GALEX catalog using the MAST API.

"""

import warnings

import numpy as np

from ..galaxies.defs import GALEX_bands
from astroquery.mast import Catalogs
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.table import vstack

from frb.surveys import surveycoord,catalog_utils
from frb.surveys.skyview import SkyView_Survey

import os

# Define the data model for GALEX data
photom = {}
photom['GALEX'] = {}
for band in GALEX_bands:
    photom["GALEX"]["GALEX"+'_{:s}'.format(band)] = '{:s}_mag'.format(band.lower())
    photom["GALEX"]["GALEX"+'_{:s}_err'.format(band)] = '{:s}_magerr'.format(band.lower())
    photom["GALEX"]["GALEX_ID"] = 'objID'
photom["GALEX"]['ra'] = 'ra'
photom["GALEX"]['dec'] = 'dec'

# Define the default set of query fields
# See: http://www.galex.caltech.edu/researcher/files/mcat_columns_long.txt
# for additional Fields
_DEFAULT_query_fields = ['distance_arcmin','objID','survey','ra','dec','e_bv']
_DEFAULT_query_fields +=['{:s}_mag'.format(band) for band in GALEX_bands]
_DEFAULT_query_fields +=['{:s}_magerr'.format(band) for band in GALEX_bands]


def _remove_duplicates_preserve_order(table, id_column):
    """Remove duplicate rows while keeping the first occurrence in order."""
    if len(table) == 0:
        return table

    _, unique_indices = np.unique(np.asarray(table[id_column]), return_index=True)
    unique_indices.sort()
    return table[unique_indices]


def _build_galex_tile_centers(coord, radius, tile_radius):
    """Build a tangent-plane grid of tile centers that fully covers the cone.

    The centers are laid out on a square lattice with spacing no larger than
    ``sqrt(2) * tile_radius``. That guarantees the tiles overlap enough to avoid
    gaps in the search geometry, while the extra outer ring ensures the edge of
    the requested cone is also covered.
    """
    radius = u.Quantity(radius, u.deg)
    tile_radius = u.Quantity(tile_radius, u.deg)
    if tile_radius <= 0 * u.deg:
        raise ValueError("tile_radius must be positive")

    step = tile_radius * np.sqrt(2.0)
    span = radius + tile_radius
    offsets = np.arange(-span.to_value(u.deg), span.to_value(u.deg) + 0.5 * step.to_value(u.deg),
                        step.to_value(u.deg))

    offset_frame = SkyOffsetFrame(origin=coord)
    tile_centers = []
    for lon_offset in offsets:
        for lat_offset in offsets:
            tile_centers.append(
                SkyCoord(lon=lon_offset * u.deg, lat=lat_offset * u.deg, frame=offset_frame).transform_to(coord.frame)
            )

    return tile_centers

class GALEX_Survey(SkyView_Survey):
    """
    A class to access all the catalogs hosted on the
    MAST database. Inherits from SurveyCoord. This
    is a super class not meant for use by itself and
    instead meant to instantiate specific children
    classes like GALEX_Survey
    """
    def __init__(self,coord,radius,**kwargs):
        SkyView_Survey.__init__(self, coord, radius, 'galex', **kwargs)

        self.Survey = "GALEX"
        self.survey = 'GALEX'
    
    def get_catalog(self, query_fields=None, print_query=False, tile=False, tile_radius=0.25*u.deg):
        """
        Query a catalog in the MAST GALEX database for
        photometry.


        Args:
            query_fields: list, optional
                A list of query fields to
                get in addition to the
                default fields.
            tile: bool, optional
                If True, split the search into smaller tiled cone searches and
                merge the results.
            tile_radius: Quantity, optional
                Radius of each tile cone search when tiling is enabled.

        
        Returns:
            catalog: astropy.table.Table
                Contains all query results
        """
        if query_fields is None:
            query_fields = _DEFAULT_query_fields
        else:
            query_fields = _DEFAULT_query_fields+query_fields

        pdict = photom['GALEX'].copy()

        if tile:
            tile_catalogs = []
            for tile_center in _build_galex_tile_centers(self.coord, self.radius, tile_radius):
                ret = Catalogs.query_region(tile_center, radius=tile_radius, catalog="GALEX")
                photom_catalog = catalog_utils.clean_cat(ret, pdict.copy(), mask_photometry=True)
                photom_catalog.keep_columns(list(pdict.keys()))
                tile_catalogs.append(photom_catalog)

            if len(tile_catalogs) == 0:
                photom_catalog = catalog_utils.ensure_empty_schema(Table(), list(pdict.keys()))
            elif len(tile_catalogs) == 1:
                photom_catalog = tile_catalogs[0]
            else:
                photom_catalog = vstack(tile_catalogs, metadata_conflicts='silent')
        else:
            ret = Catalogs.query_region(self.coord, radius=self.radius, catalog="GALEX")
            photom_catalog = catalog_utils.clean_cat(ret, pdict, mask_photometry=True)
            photom_catalog.keep_columns(list(pdict.keys()))

        photom_catalog = catalog_utils.sort_by_separation(photom_catalog, self.coord,
                                                          radec=('ra','dec'), add_sep=True)

        # Remove duplicates after sorting so the closest copy of each object wins.
        photom_catalog = _remove_duplicates_preserve_order(photom_catalog, "GALEX_ID")

        self.catalog = photom_catalog
        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        #Validate
        self.validate_catalog()

        #Return
        return self.catalog.copy()

    def get_image(self, imsize, band='NUV'):
        """Retrieve a SkyView FITS image for GALEX."""
        return SkyView_Survey.get_image(self, imsize=imsize, band=band)

    def get_cutout(self, imsize, band='NUV'):
        """Deprecated alias for FITS image retrieval."""
        warnings.warn(
            "get_cutout() returns FITS products for this survey and is deprecated; use get_image() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_image(imsize=imsize, band=band)