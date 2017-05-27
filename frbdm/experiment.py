""" Module for calculations related to FRB experiments
"""

from __future__ import print_function, absolute_import, division, unicode_literals

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
        setup()

    def setup(self):
        """ Load the characteristics of the experiment
        """
        if self.name == "CHIME":

