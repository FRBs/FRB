""" Module for calculations related to FRB experiments
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from pkg_resources import resource_filename

from . import utils

class Experiment(object):
    """
    """
    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
        """
        self.name = name
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