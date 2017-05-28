""" Module for an FRB event
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from . import utils

class FRB(object):
    """
    """
    def __init__(self, S, width, nu_c, DM, coord=None):
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
        self.width = width
        self.nu_c = nu_c
        self.DM = DM
        # Coord
        if coord is not None:
            self.coord = utils.radec_to_coord(coord)
        #

