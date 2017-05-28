""" Module for calculations related to FRB experiments
"""

from __future__ import print_function, absolute_import, division, unicode_literals


def loadjson(filename):
    """
    Parameters
    ----------
    filename : str

    Returns
    -------
    obj : dict

    """
    import json, gzip
    #
    if filename.endswith('.gz'):
        with gzip.open(filename, "rb") as f:
            obj = json.loads(f.read().decode("ascii"))
    else:
        with open(filename, 'rt') as fh:
            obj = json.load(fh)

    return obj


def loadyaml(filename):
    from astropy.io.misc import yaml as ayaml
    # Read yaml
    with open(filename, 'r') as infile:
        data = ayaml.load(infile)
    # Return
    return data



def radec_to_coord(radec):
    """ Converts one of many of Celestial Coordinates
    `radec` formats to an astropy SkyCoord object. Assumes
    J2000 equinox.

    Parameters
    ----------
    radec : str or tuple or SkyCoord or list
        Examples:
        'J124511+144523',
        '124511+144523',
        'J12:45:11+14:45:23',
        ('12:45:11','+14:45:23')
        ('12 45 11', +14 45 23)
        ('12:45:11','14:45:23')  -- Assumes positive DEC
        (123.123, 12.1224) -- Assumed deg
        [(123.123, 12.1224), (125.123, 32.1224)]

    Returns
    -------
    coord : SkyCoord
      Converts to astropy.coordinate.SkyCoord (as needed)
      Returns a SkyCoord array if input is a list


    """
    from astropy.coordinates import SkyCoord

    # RA/DEC
    if isinstance(radec, (tuple)):
        if isinstance(radec[0], basestring):
            if radec[1][0] not in ['+', '-']:  #
                DEC = '+'+radec[1]
                warnings.warn("Assuming your DEC is +")
            else:
                DEC = radec[1]
            #
            coord = SkyCoord(radec[0]+DEC, frame='fk5',
                                  unit=(u.hourangle, u.deg))
        else:
            coord = SkyCoord(ra=radec[0], dec=radec[1], unit='deg')
    elif isinstance(radec,SkyCoord):
        coord = radec
    elif isinstance(radec,basestring):
        # Find first instance of a number (i.e. strip J, SDSS, etc.)
        for ii in range(len(radec)):
            if radec[ii].isdigit():
                break
        radec = radec[ii:]
        #
        if ':' in radec:
            coord = SkyCoord(radec, frame='fk5', unit=(u.hourangle, u.deg))
        else:  # Add in :
            if ('+' in radec) or ('-' in radec):
                sign = max(radec.find('+'), radec.find('-'))
            else:
                raise ValueError("radec must include + or - for DEC")
            newradec = (radec[0:2]+':'+radec[2:4]+':'+radec[4:sign+3] +':'+radec[sign+3:sign+5]+':'+radec[sign+5:])
            coord = SkyCoord(newradec, frame='fk5', unit=(u.hourangle, u.deg))
    elif isinstance(radec,list):
        clist = []
        for item in radec:
            clist.append(radec_to_coord(item))
        # Convert to SkyCoord array
        ras = [ii.fk5.ra.value for ii in clist]
        decs = [ii.fk5.dec.value for ii in clist]
        return SkyCoord(ra=ras, dec=decs, unit='deg')
    else:
        raise IOError("Bad input type for radec")
    # Return
    return coord
