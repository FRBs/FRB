""" Module for an FRB event
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from pkg_resources import resource_filename

from frb import utils

class FRB(object):
    """
    """
    def __init__(self, S, Wi, nu, coord=None):
        """
        Parameters
        ----------
        S : Quantity
          Source density of the burst
        Wi : Quantity
        nu : Quantity
          Characteristic frequency of the event
        coord
        """
        self.S = S  # Source density (e.g. Jy)
        #
        self.setup()

    def setup(self):
        """ Load the characteristics of the experiment
        """
        if self.name == "CHIME":
            self.data_file=resource_filename('frbdm', 'data/experiments/chime.yaml')
            self.data = utils.loadyaml(self.data_file)
        else:
            raise IOError("Not ready for this Experiment: {:s}".format(self.name))

    def signal_to_noise(self):
        """
        Returns
        -------
        s2n : float

        """