""" Module for I/O activities in the FRB repo
"""
import importlib_resources
import gzip
import json
import os



def load_dla_fits():
    dla_fit_file = str(importlib_resources.files('frb.data.IGM')/'dla_fits.json')
    dla_fits = loadjson(dla_fit_file)
    # Return
    return dla_fits


def loadjson(filename):
    """ Load a python object saved in JSON
    Parameters
    ----------
    filename : str

    Returns
    -------
    obj : unknown
      Often a dict
    """
    if filename.endswith('.gz'):
        with gzip.open(filename, "rb") as f:
            obj = json.loads(f.read().decode("ascii"))
    else:
        with open(filename, 'rt') as fh:
            obj = json.load(fh)

    return obj

