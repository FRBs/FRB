""" Module for an FRB event
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from . import utils

class FRB(object):
    """
    """
    def __init__(self, S, nu_c, DM, coord=None):
        """
        Parameters
        ----------
        S : Quantity
          Source density of the burst
        width : Quantity
          Width
        nu_c : Quantity
          Centre frequency
        DM : Quantity
        coord : multi-format, optional
          RA/DEC in one of many formats (see utils.radec_to_coord)
        """
        self.S = S
        self.nu_c = nu_c
        self.DM = DM
        # Coord
        if coord is not None:
            self.coord = utils.radec_to_coord(coord)
        #

    def set_width(self, wtype, value, overwrite=False):
        """ Set a Width value
        Parameters
        ----------
        wtype : str
          Indicates the type of width
        value : Quantity
        overwrite : bool, optional
        """
        # Restrict types
        assert wtype in ['Wi', 'Wb', 'Wobs']  # Intrinsic, broadened, observed
        # Check
        if hasattr(self, wtype) and (not overwrite):
            raise IOError("Width type {:s} is already set!".format(wtype))
        # Vette
        try:
            value.to('s')
        except:
            raise IOError("Bad Quantity for value")
        # Set
        setattr(self, wtype, value)


