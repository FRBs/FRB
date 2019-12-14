""" Module for calculations related to FRB experiments
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from pkg_resources import resource_filename

from astropy import units as u

from frb import utils

class Experiment(object):
    """
    """
    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
          See YAML files in data/experiment
        """
        self.name = name
        #
        self.setup()

    def setup(self):
        """ Load the characteristics of the experiment
        """
        self.data_file=resource_filename('frb', 'data/experiments/{:s}.yaml'.format(
            self.name.lower()))
        self.data = utils.loadyaml(self.data_file)

    def signal_to_noise(self, frb, beta=1., T_Sky=None, t_scatt=None):
        """
        Follows Cordes & McLaughlin 2003

        Parameters
        ----------
        frb : FRB
        beta : float, optional
           Factor for digitization losses
        t_scatt : Quantity, optional
          Scattering time

        Returns
        -------
        s2n : float

        """
        # TODO -- Add t_samp to experiment data
        t_samp = 0 * u.s
        # t_scatt
        if t_scatt is None:
            try:
                t_scatt = frb.t_scatt
            except AttributeError:
                t_scatt = 0.*u.s
        # t_chan  (Lorimer & Kramer 2005)
        t_chan = 8.3e-6*u.s * (self.data['Dnu'].to('MHz').value/self.data['Channels']) * (
            frb.nu_c.to('GHz').value)**(-3) * frb.dm.to('pc/cm**3').value
        # Wb
        Wb = np.sqrt(frb.Wi**2 + t_chan**2 + t_samp**2 + t_scatt**2)
        # T_Sky
        if T_Sky is None:
            T_Sky = utils.Tsky(frb.nu_c)
        #
        sqrt_term = np.sqrt(Wb/(self.data['np']*self.data['Dnu']))
        # Here we go
        s2n = frb.S * self.data['G'] * frb.Wi / beta / (
            self.data['Trec'] + T_Sky) / sqrt_term
        # Return
        return s2n.decompose()

    def __repr__(self):
        txt = '<{:s}: name={:s} data={}'.format(
            self.__class__.__name__, self.name, self.data)
        # Finish
        txt = txt + '>'
        return (txt)
