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
