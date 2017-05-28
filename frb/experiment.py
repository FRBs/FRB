""" Module for calculations related to FRB experiments
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

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
            self.data_file=resource_filename('frb', 'data/experiments/chime.yaml')
            self.data = utils.loadyaml(self.data_file)
        else:
            raise IOError("Not ready for this Experiment: {:s}".format(self.name))

    def signal_to_noise(self, frb, beta=1., T_Sky=None):
        """
        Follows Cordes & McLaughlin 2003

        Parameters
        ----------
        frb : FRB
        beta : float, optional
           Factor for digitization losses

        Returns
        -------
        s2n : float

        """
        # T_Sky
        if T_Sky is None:
            T_Sky = utils.Tsky(frb.nu_c)
        #
        sqrt_term = np.sqrt(frb.Wb/(self.data['np']*self.data['Dnu']))
        # Here we go
        s2n = frb.S * self.data['G'] / self.data['beta'] / (
            self.data['Trec'] + T_Sky) / sqrt_term
        # Return
        return s2n