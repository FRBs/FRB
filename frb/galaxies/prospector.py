""" WARNING -- WIP -- Methods related to running and reading Prospector outputs """

from frb.galaxies import nebular
from frb.galaxies import utils as galaxy_utils
from frb.galaxies import defs

# TODO -- THIS NEEDS MORE FILTERS
prospector_filter_dict = {}
for key in defs.SDSS_bands:
    prospector_filter_dict['SDSS_{:s}'.format(key)] = 'sdss_{:s}'.format(key)

def load_spec_from_specdb(host_obj, INSTR:str, deredden=True,
                          ebv_method = 'SandF'):
    # Grab it
    meta, spectrum = host_obj.get_metaspec(instr=INSTR)

    # De-redden?
    if deredden:
        EBV = nebular.get_ebv(host_obj.coord, definition=ebv_method)['meanValue'] 
        spectrum = galaxy_utils.deredden_spec(spectrum, EBV)
    
    # Return
    return spectrum.wavelength.values, spectrum.flux.values, spectrum.sig.values


def grab_magnitudes(host_obj, FRB_filters:list):
    # Load up magnitudes from a host object

    prospector_filters, m_AB, m_AB_err = [], [], []

    for ifilt in FRB_filters:
        # Vette
        assert ifilt in defs.valid_filters
        # Translate
        prospector_filters.append(prospector_filter_dict[ifilt])
        # Grab data
        m_AB.append(host_obj.photom[ifilt])
        m_AB_err.append(host_obj.photom[ifilt+'_err'])
    
    # Return
    return prospector_filters, m_AB, m_AB_err
    
    


