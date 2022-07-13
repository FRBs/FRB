# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table

from linetools.spectra import xspectrum1d

from frb.galaxies import frbgalaxy, defs, utils
from frb.frb import FRB

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_frbhost():
    repeater_coord = SkyCoord('05h31m58.698s +33d8m52.59s', frame='icrs')  # Use as host here
    # Instantiate
    frb121102 = FRB.by_name('FRB20121102A')
    host121102 = frbgalaxy.FRBHost(repeater_coord.ra.value, repeater_coord.dec.value, frb121102)
    # Redshift
    host121102.set_z(0.19273, 'spec', err=0.00008)
    # Add a few nebular lines  (Tendulkar+17)
    neb_lines = {}
    neb_lines['Halpha'] = 0.652e-16
    neb_lines['Halpha_err'] = 0.009e-16
    neb_lines['Halpha_Al'] = 0.622
    #
    neb_lines['Hbeta'] = 0.118e-16
    neb_lines['Hbeta_err'] = 0.011e-16
    neb_lines['Hbeta_Al'] = 0.941
    AV = 2.42
    # Deal with Galactic extinction
    for key in neb_lines.keys():
        if '_err' in key:
            continue
        if 'Al' in key:
            continue
        # Ingest
        host121102.neb_lines[key] = neb_lines[key] * 10 ** (neb_lines[key + '_Al'] * AV / 2.5)
        host121102.neb_lines[key + '_err'] = neb_lines[key + '_err'] * 10 ** (neb_lines[key + '_Al'] * AV / 2.5)
    # Vette
    for key in host121102.neb_lines.keys():
        if '_err' in key:
            continue
        assert key in defs.valid_neb_lines
    # Morphology
    host121102.morphology['reff_ang'] = 0.41
    host121102.morphology['reff_ang_err'] = 0.06
    #
    host121102.morphology['n'] = 2.2
    host121102.morphology['n_err'] = 1.5
    #
    host121102.morphology['b/a'] = 0.25
    host121102.morphology['b/a_err'] = 0.13
    # Vet
    assert host121102.vet_one('morphology')
    # Vet all
    assert host121102.vet_all()
    # Write
    outfile = data_path('test_frbhost.json')
    host121102.write_to_json(outfile=outfile)


def test_read_frbhost():
    # This test will fail if the previous failed
    outfile = data_path('test_frbhost.json')
    frb121102 = FRB.by_name('FRB20121102A')
    host121102 = frbgalaxy.FRBHost.from_json(frb121102, 
                                             outfile)
    # Test
    assert host121102.frb.frb_name == 'FRB20121102A'
    assert np.isclose(host121102.morphology['b/a'], 0.25)
    assert host121102.vet_all()


def test_luminosity():
    # This test will fail if the previous failed
    outfile = data_path('test_frbhost.json')
    frb121102 = FRB.by_name('FRB20121102A')
    host121102 = frbgalaxy.FRBHost.from_json(frb121102, outfile)
    Lum_Ha, Lum_Ha_err = host121102.calc_nebular_lum('Halpha')
    # Test
    assert Lum_Ha.unit == units.erg/units.s
    assert np.isclose(Lum_Ha.value, 2.9447660789951146e+40)

    # Remove
    os.remove(outfile)


def test_get_spectra():
    """
    Check grabbing spectra from the specDB

    Note:  This requires that specdb be installed..
    """
    try:
        import specdb
    except ImportError:
        assert True
        return
    # Check for specDB file
    if utils.load_specdb() is None:
        assert True
        return
    # Do it!
    frb180924 = FRB.by_name('FRB20180924B')
    host180924 = frb180924.grab_host()
    meta, xspec = host180924.get_metaspec()
    # Test
    assert isinstance(meta, Table)
    assert isinstance(xspec, xspectrum1d.XSpectrum1D)
    assert xspec.nspec == 1
    #
    meta, xspec = host180924.get_metaspec(instr='MUSE')
    # Test
    assert isinstance(xspec, xspectrum1d.XSpectrum1D)
    assert xspec.nspec == 1
    #
    meta, xspecs = host180924.get_metaspec(return_all=True)
    assert xspecs.nspec == 2


def test_table():
    frbs, hosts = utils.list_of_hosts()

    host_tbl, host_units = utils.build_table_of_hosts()

    # Tests
    assert len(host_tbl) > 15
    assert np.isclose(np.max(host_tbl.P_Ox), 1.)
