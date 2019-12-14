""" Module to generate individual FRB files """

from astropy import units

from frb import frb

def frb_121102():
    frb121102 = frb.FRB('FRB121102', 'J053158.7+330852.5',
                        558.1*units.pc/units.cm**3,
                        z_frb=0.19273)
    # NE2001
    frb121102.set_dmISM()
    # Error ellipse
    frb121102.set_ee(0.1, 0.1, 0., 95.)
    # Write
    frb121102.write_to_json()
    # Test
    frb121102.from_json('FRB121102.json')


# Command line execution
if __name__ == '__main__':
    # FRB 121102
    frb_121102()



