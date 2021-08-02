""" Methods to load up additional data """

from pkg_resources import resource_filename
import os

import pandas

from astropy import units

def load_mannings2021():
    """ Load a table of measurements from the 
    Mannings+2021 paper

    Returns:
        tuple: measurements (pandas.DataFrame), units and description (dict)
    """
    mannings2021_file = os.path.join(resource_filename('frb','data'),
                                          'Galaxies','Additional','Mannings2021', 
                                          'galaxy_measurements.csv')
    mannings2021 = pandas.read_csv(mannings2021_file)                                        

    # Units and comments
    tbl_units = {}
    tbl_units['UV_SB'] = dict(unit=units.Unit('mJy/arcsec**2'),
                              comment='UV SB in UV_filter at FRB position')
    tbl_units['UV_SB_err'] = dict(unit=units.Unit('mJy/arcsec**2'),
                              comment='Error in UV_SB')
    tbl_units['SSFR'] = dict(unit=units.Unit('M_sun/yr/kpc**2'),
                              comment='Surface density of Star formation at FRB position')
    tbl_units['SSFR_err'] = dict(unit=units.Unit('M_sun/yr/kpc**2'),
                              comment='Error in SSFR')
    tbl_units['SMStar'] = dict(unit=units.Unit('M_sun/kpc**2'),
                              comment='Surface density of Stars at FRB position')
    tbl_units['SMStar_err'] = dict(unit=units.Unit('M_sun/kpc**2'),
                              comment='Error in SMStar')
    tbl_units['reff_iso'] = dict(unit=units.Unit('arcsec'),
                              comment='Effective Isophotal radius (angular)')
    tbl_units['reff_iso_err'] = dict(unit=units.Unit('arcsec'),
                              comment='Error in reff_iso')
    tbl_units['UVff'] = dict(unit=None, comment='Fractional flux at FRB location in UV_filter')
    tbl_units['IRff'] = dict(unit=None, comment='Fractional flux at FRB location in IR_filter')
    tbl_units['IRfe'] = dict(unit=None, comment='Enclosed flux at FRB location in IR_filter')
    tbl_units['IRlim'] = dict(unit=None, comment='Limiting IR flux (mag) at FRB location after subtracting off galfit model.')

    # Return
    return mannings2021, tbl_units
