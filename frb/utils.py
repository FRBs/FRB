""" Module for calculations related to FRB experiments
"""
import os
import numpy as np

from astropy import units
import warnings

import json, gzip

from IPython import embed

# Simple method to help with value/units in
def assign_value(tfrb, key, ilist, tbl_units):
    """
    Assign a value into a list, dealing with astropy.units.Quantity objects

    The input list and dict of units may be modified in place

    Args:
        tfrb (frb.frb.FRB):
        key (str):
        ilist (list):
            Input list
        tbl_units (dict):
            Input dict of units

    """
    val, unit = get_valunit(getattr(tfrb, key))
    ilist.append(val)
    # Deal with units
    if key not in tbl_units.keys():
        if unit is not None:
            tbl_units[key] = unit.to_string()
        else:
            tbl_units[key] = None
    else:
        if unit is not None:
            assert tbl_units[key] == unit.to_string()


def get_valunit(item):
    """
    Grab the value and unit of an input item allowing
    for an astropy.units.Quantity

    Args:
        item (units.Quantity or any Python object):

    Returns:
        value, unit

    """
    if isinstance(item, units.Quantity):
        return item.value, item.unit
    else:
        return item, None


def jsonify(obj, debug=False):
    """ Recursively process an object so it can be serialised in json
    format.

    WARNING - the input object may be modified if it's a dictionary or
    list!

    Parameters
    ----------
    obj : any object
    debug : bool, optional

    Returns
    -------
    obj - the same obj is json_friendly format (arrays turned to
    lists, np.int64 converted to int, np.float64 to float, and so on).

    """
    if isinstance(obj, np.float64):
        obj = float(obj)
    elif isinstance(obj, np.float32):
        obj = float(obj)
    elif isinstance(obj, np.int32):
        obj = int(obj)
    elif isinstance(obj, np.int64):
        obj = int(obj)
    elif isinstance(obj, np.int16):
        obj = int(obj)
    elif isinstance(obj, np.bool_):
        obj = bool(obj)
    elif isinstance(obj, np.bytes_):
        obj = str(obj)
    elif isinstance(obj, units.Quantity):
        if obj.size == 1:
            obj = dict(value=obj.value, unit=obj.unit.to_string())
        else:
            obj = dict(value=obj.value.tolist(), unit=obj.unit.to_string())
    elif isinstance(obj, np.ndarray):  # Must come after Quantity
        obj = obj.tolist()
    elif isinstance(obj, dict):
        for key, value in obj.items():
            obj[key] = jsonify(value, debug=debug)
    elif isinstance(obj, list):
        for i,item in enumerate(obj):
            obj[i] = jsonify(item, debug=debug)
    elif isinstance(obj, tuple):
        obj = list(obj)
        for i,item in enumerate(obj):
            obj[i] = jsonify(item, debug=debug)
        obj = tuple(obj)
    elif isinstance(obj, units.Unit):
        obj = obj.name
    elif obj is units.dimensionless_unscaled:
        obj = 'dimensionless_unit'

    if debug:
        print(type(obj))
    return obj


def loadjson(filename):
    """
    Parameters
    ----------
    filename : str

    Returns
    -------
    obj : dict

    """
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


def savejson(filename, obj, overwrite=False, indent=None, easy_to_read=False,
             **kwargs):
    """ Save a python object to filename using the JSON encoder.

    Parameters
    ----------
    filename : str
    obj : object
      Frequently a dict
    overwrite : bool, optional
    indent : int, optional
      Input to json.dump
    easy_to_read : bool, optional
      Another approach and obj must be a dict
    kwargs : optional
      Passed to json.dump

    Returns
    -------

    """
    import io

    if os.path.lexists(filename) and not overwrite:
        raise IOError('%s exists' % filename)
    if easy_to_read:
        if not isinstance(obj, dict):
            raise IOError("This approach requires obj to be a dict")
        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write(json.dumps(obj, sort_keys=True, indent=4,
                               separators=(',', ': '), **kwargs))
    else:
        if filename.endswith('.gz'):
            with gzip.open(filename, 'wt') as fh:
                json.dump(obj, fh, indent=indent, **kwargs)
        else:
            with open(filename, 'wt') as fh:
                json.dump(obj, fh, indent=indent, **kwargs)

def name_from_coord(coord, precision=(2,1)):
    """ Generate a standard JXXXXXX.XX+XXXXXX.X name from a SkyCoord object

    Parameters
    ----------
    coord : SkyCoord
    precision : tuple, optional
      Number of decimal places to include in name

    Returns
    -------
    name : str
      In JXX format
    """
    name = 'J{:s}{:s}'.format(coord.ra.to_string(unit=units.hour,sep='',pad=True,precision=precision[0]),
            coord.dec.to_string(sep='',pad=True,alwayssign=True,precision=precision[1]))
    # Return
    return name

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
        if isinstance(radec[0], str):
            if radec[1][0] not in ['+', '-']:  #
                DEC = '+'+radec[1]
                warnings.warn("Assuming your DEC is +")
            else:
                DEC = radec[1]
            #
            coord = SkyCoord(radec[0]+DEC, frame='icrs',
                                  unit=(units.hourangle, units.deg))
        else:
            coord = SkyCoord(ra=radec[0], dec=radec[1], unit='deg')
    elif isinstance(radec,SkyCoord):
        coord = radec
    elif isinstance(radec,str):
        # Find first instance of a number (i.e. strip J, SDSS, etc.)
        for ii in range(len(radec)):
            if radec[ii].isdigit():
                break
        radec = radec[ii:]
        #
        if ':' in radec:
            coord = SkyCoord(radec, frame='icrs', unit=(units.hourangle, units.deg))
        else:  # Add in :
            if ('+' in radec) or ('-' in radec):
                sign = max(radec.find('+'), radec.find('-'))
            else:
                raise ValueError("radec must include + or - for DEC")
            newradec = (radec[0:2]+':'+radec[2:4]+':'+radec[4:sign+3] +':'+radec[sign+3:sign+5]+':'+radec[sign+5:])
            coord = SkyCoord(newradec, frame='icrs', unit=(units.hourangle, units.deg))
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


def Tsky(nu):
    """ Sky temperature
     Tsky for all other surveys has been evaluated assuming an average sky
     temperature of 34 K at 408 MHz and a spectral index of -2.6
    Follows Haslam et al. 1982

    Parameters
    ----------
    nu : Quantity

    Returns
    -------
    Tsky : Quantity


    """
    # TODO  -- Need some guidance here
    return 34*units.K * (nu/(408*units.MHz))**(-2.6)

def parse_frb_name(name:str, prefix='FRB'):
    """Parse the incoming name to generate a 'proper' FRB name

    Args:
        name (str): [description]
        prefix (str, optional): [description]. Defaults to 'FRB'.

    Raises:
        IOError: [description]

    Returns:
        str: The proper FRB name
    """
    
    if name[0:3] == 'FRB':
        # Strip FRB
        return parse_frb_name(name[3:], prefix=prefix)
    elif len(name) in [6,7]: # Original YYMMDD or YYMMDDL  (with L a 'letter')
        # Append 20
        return parse_frb_name('20'+name, prefix=prefix)
    elif len(name) in [8,9]: # YYYYMMDD or YYYYMMDDL
        # All set
        return prefix+name
    else:
        raise IOError(f"Not prepared for this type of format: {name}")



def log10_to_linear_errors(value, log_upper_error, log_lower_error):
    """
    Convert log10 errors to linear errors.
    
    Parameters:
    -----------
    value : float
        The actual value of the quantity (not log10)
    log_upper_error : float
        Upper error in log10 space
    log_lower_error : float
        Lower error in log10 space
    
    Returns:
    --------
    linear_upper_error : float
        Upper error in linear space
    linear_lower_error : float
        Lower error in linear space
    """
    # Calculate the log10 of the value
    log_value = np.log10(value)
    
    # Calculate upper and lower bounds in log10 space
    log_upper_bound = log_value + log_upper_error
    log_lower_bound = log_value - log_lower_error
    
    # Convert back to linear space
    upper_bound = 10**log_upper_bound
    lower_bound = 10**log_lower_bound
    
    # Calculate linear errors
    linear_upper_error = upper_bound - value
    linear_lower_error = value - lower_bound
    
    return linear_upper_error, linear_lower_error
