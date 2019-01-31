""" Top-level module to build or re-build the JSON files for
FRB host galaxies"""

from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
from astropy import units

from frb.galaxies import frbgalaxy, defs

def frb_121102():
    FRB_coord = SkyCoord('05h31m58.698s +33d8m52.59s', frame='icrs')
    # Eyeball Tengulkar+17 PA
    gal_coord = FRB_coord.directional_offset_by(-45 * units.deg, 286e-3 * units.arcsec)

    # Instantiate
    host121102 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, '121102')

    # Redshift
    host121102.set_z(0.19273, 'spec', err=0.00008)

    # Nebular lines
    neb_lines = {}
    neb_lines['Ha'] = 0.652e-16
    neb_lines['Ha_err'] = 0.009e-16
    neb_lines['Ha_Al'] = 0.622
    #
    neb_lines['Hb'] = 0.118e-16
    neb_lines['Hb_err'] = 0.011e-16
    neb_lines['Hb_Al'] = 0.941
    #
    neb_lines['[OIII] 5007'] = 0.575e-16
    neb_lines['[OIII] 5007_err'] = 0.011e-16
    neb_lines['[OIII] 5007_Al'] = 0.911

    AV = 2.42

    # Extinction correct
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
    # Vette
    for key in host121102.morphology.keys():
        if '_err' in key:
            continue
        assert key in defs.valid_morphology

    # Derived quantities
    host121102.derived['M_r'] = -17.0  # AB; Tengulkar+17
    host121102.derived['SFR_nebular'] = 0.23  # MSun/yr; Tengulkar+17
    host121102.derived['Mstar'] = 5.5e7  # Msun; Tengulkar+17
    host121102.derived['Mstar_err'] = 1.5e7  # Msun; Tengulkar+17
    # Vette
    for key in host121102.derived.keys():
        if '_err' in key:
            continue
        assert key in defs.valid_derived

    # Write
    path = resource_filename('frb', 'data/Galaxies/121102')
    host121102.write_to_json(path=path, overwrite=True)


# Command line execution
if __name__ == '__main__':
    # FRB 121102
    frb_121102()




