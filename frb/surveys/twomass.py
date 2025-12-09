"""
Slurp data from 2MASS catalog.

"""

import numpy as np

from astropy import units as u,utils as astroutils
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, join
from ..galaxies.defs import MASS_bands
from astroquery.ipac.irsa import Irsa

from frb.surveys import surveycoord,catalog_utils,images

import os

# Define the data model for 2MASS data
photom = {}
photom['2MASS'] = {}
for band in MASS_bands:
    photom["2MASS"]["2MASS"+'_{:s}'.format(band)] = '{:s}_m'.format(band) # Many options for apertures, tbd
    photom["2MASS"]["2MASS"+'_{:s}_err'.format(band)] = '{:s}_msig'.format(band)
    photom["2MASS"]["2MASS_ID"] = 'designation'
photom["2MASS"]['ra'] = 'ra'
photom["2MASS"]['dec'] = 'dec'

# Define the default set of query fields
# http://tdc-www.harvard.edu/catalogs/tmpsc.format.html For PSC
# https://www.ipac.caltech.edu/2mass/releases/second/doc/ancillary/xscformat.html For XSC
_DEFAULT_query_fields = ['designation','survey','ra','dec']
_DEFAULT_query_fields +=['{:s}_m'.format(band) for band in MASS_bands]
_DEFAULT_query_fields +=['{:s}_msig'.format(band) for band in MASS_bands]

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
    
    def get_catalog(self,query_fields=None):
        """
        Query a catalog in the IRSA 2MASS survey for
        photometry.

        Args:
            query_fields: list, optional
                A list of query fields to
                get in addition to the
                default fields.
        
        Returns:
            catalog: astropy.table.Table
                Contains all query results
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
        
        photom_catalog = catalog_utils.clean_cat(ret,pdict) # rename columns

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