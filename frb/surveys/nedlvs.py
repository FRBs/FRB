#!/bin/env python3
import os
from . import surveycoord
from astropy.table import Table
from frb.defs import frb_cosmo
from astropy.coordinates import SkyCoord
from astropy.cosmology import Cosmology
from astropy import units as u

import numpy as np


class NEDLVS(surveycoord.SurveyCoord):
    """
    A class to handle local NED Local Volume Sample queries.
    This requires the LVS table to be downloaded
    from https://ned.ipac.caltech.edu/NED::LVS/
    and linked via the environment variable NEDLVS

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.units.Quantity, optional): Search radius around the
            coordinate. Must not exceed 90 deg.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used to convert
            redshifts to distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: u.Quantity = 90.*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        assert 'NEDLVS' in os.environ, "NEDLVS environment variable not set. Please download the LVS table from https://ned.ipac.caltech.edu/NED::LVS/ and set the environment variable NEDLVS to the path of the downloaded file."        
        self.survey = 'NEDLVS'
        self.datapath = os.environ['NEDLVS']
        if cosmo is None:
            self.cosmo = frb_cosmo
        else:
            self.cosmo = cosmo

        # Read in the data and store in memory
        self.datatab = Table.read(self.datapath)
        self.datatab['coord'] = SkyCoord(self.datatab['ra'], self.datatab['dec'], unit="deg")
        self.datatab['ang_sep'] = self.coord.separation(self.datatab['coord']).to('arcmin')

        # Set redshift distances using the cosmology of choice
        redshift_dist_sources = self.datatab['DistMpc_method']=='Redshift'
        self.datatab['DistMpc'][redshift_dist_sources] = self.cosmo.luminosity_distance(self.datatab['z'][redshift_dist_sources])
        self.datatab['phys_sep'] = self.datatab['DistMpc']*u.Mpc*np.sin(self.datatab['ang_sep'].to('rad').value)
    
    def get_column_names(self) -> list[str]:
        """
        Get the names of the columns of the NEDLVS table.

        Returns:
            list of str: Column names that can be passed as ``query_fields``
            to :meth:`get_catalog`.

        """
        return self.datatab.colnames

    def get_catalog(self, z_lim: float = np.inf,
                    impact_par_lim: u.Quantity = np.inf*u.Mpc,
                    query_fields: list[str] | None = None,
                    print_query: bool = False) -> Table:
        """
        Get the catalog of objects within the given limits of redshift, impact parameter, and angular separation.
        The angular separation limit is the search radius of the survey (``self.radius``).

        Args:
            z_lim (float, optional): The maximum redshift of the objects to include in the catalog.
            impact_par_lim (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.
            query_fields (list of str, optional): The fields to include in the catalog. If None, the default fields are used.
            print_query (bool, optional): Print the query limits and fields.

        Returns:
            astropy.table.Table: A table of objects within the given limits.

        Raises:
            AssertionError: If a requested field is not in the NEDLVS table
                or ``self.radius`` is greater than 90 deg.

        """
        if query_fields is None:
            query_fields = ['objname', 'ra', 'dec', 'ebv', 'z', 'z_unc', 'z_tech', 'DistMpc', 'DistMpc_unc', 'DistMpc_method', 'Mstar', 'Mstar_unc', 'ang_sep', 'phys_sep'] 
        else:
            assert np.isin(query_fields, self.datatab.colnames).all(), "One or more of the requested fields is not in the NEDLVS table. Check the column names with get_column_names()."
        # ...
        if print_query:
            print(f"Querying NEDLVS for objects within {self.radius} of {self.coord} with z < {z_lim} and impact parameter < {impact_par_lim}.")
            print(f"Query fields: {query_fields}")
        distance_cut = self.datatab['DistMpc']<self.cosmo.luminosity_distance(z_lim).to('Mpc').value #Only need foreground objects
        valid_distances = self.datatab['DistMpc']>0 # Exclude weird sources with negative distances
        phys_sep_cut = self.datatab['phys_sep']<impact_par_lim # Impact param within limit

         # Make sure the earth is not between the FRB and the galaxy
        assert self.radius<=90*u.deg, "The radius of the search cone is too large. Please set it to a value less than 90 degrees."

        ang_sep_cut = self.datatab['ang_sep']<self.radius
        is_nearby_fg = valid_distances&distance_cut & phys_sep_cut & ang_sep_cut
        
        close_by = self.datatab[is_nearby_fg][query_fields]

        self.catalog = close_by
        
        # Normalize and validate (preserves ang_sep/phys_sep, adds canonical separation if needed)
        self.validate_catalog()
        
        return self.catalog