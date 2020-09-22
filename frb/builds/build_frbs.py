""" Module to generate individual FRB files """

from pkg_resources import resource_filename

import numpy as np

from astropy import units
from astropy.coordinates import SkyCoord

from frb import frb


def frb_121102():
    """
    FRB 121102
        All of the data currently comes from Tendulkar et al. 2017
        https://ui.adsabs.harvard.edu/abs/2017ApJ...834L...7T/abstract
    """
    frb121102 = frb.FRB('FRB121102', 'J053158.7+330852.5',
                        558.1*units.pc/units.cm**3,
                        z_frb=0.19273)
    # NE2001
    frb121102.set_DMISM()
    # Error ellipse
    frb121102.set_ee(0.1, 0.1, 0., 95.)
    # References
    frb121102.refs = ['Tendulkar2017']
    # Write
    frb121102.write_to_json()
    # Test
    frb121102.from_json('FRB121102.json')


def frb_180916():
    """
     FRB 180916.J0158+65
        All of the data currently comes from Marcote et al. 2020
        https://ui.adsabs.harvard.edu/abs/2020Natur.577..190M/abstract
    """
    coord = SkyCoord("01h58m00.75017s 65d43m00.3152s", frame='icrs')
    name = 'FRB180916'
    frb180916 = frb.FRB(name, coord,
                        348.76*units.pc / units.cm**3,
                        z_frb=0.0337)
    # Error ellipse
    frb180916.set_ee(0.0023, 0.0023, 0., 68.)
    # Error in DM
    frb180916.DM_err = 0.10 * units.pc / units.cm**3
    # NE2001
    frb180916.set_DMISM()
    #
    #frb180924.fluence = 16 * units.Jy * units.ms
    #frb180924.fluence_err = 1 * units.Jy * units.ms
    #
    # References
    frb180916.refs = ['Marcote2020']
    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb180916.write_to_json(path=path)

def frb_180924():
    """
    FRB 180924
        All of the data currently comes from Bannister et al. 2019
        https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B/abstract
        Now including updated values on FRB position, DM, pulse width 
        and scattering time from Day et al. 2020
    """
    frb180924 = frb.FRB('FRB180924', 'J214425.255-405400.10',
                        362.16*units.pc / units.cm**3,
                        z_frb=0.3212)
    # Error in DM (Day 2019)
    frb180924.DM_err = 0.01 * units.pc / units.cm**3

    # NE2001
    frb180924.set_DMISM()

    # FRB properties
    frb180924.fluence = 16 * units.Jy * units.ms    # From Bhandari+20
    frb180924.fluence_err = 1 * units.Jy * units.ms # -- //-- 
    frb180924.RM = 22 * units.rad / units.m**2 
    frb180924.RM_err = 2 * units.rad / units.m**2
    frb180924.tau = 0.68 * units.ms
    frb180924.tau_err = 0.03 * units.ms
    frb180924.lpol = 80.  # %
    frb180924.lpol_err = 10.
    # Error ellipse
    frb180924.set_ee(a=100./1e3, b=100./1e3, theta=0., cl=68.) #statistical
    frb180924.set_ee(a=0.09, b=0.09, theta=0, cl=68., stat=False) #systematic

    # References
    frb180924.refs = ['Bannister2019', 'Day2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb180924.write_to_json(path=path)

def frb_181112():
    """
    Generate the JSON file for FRB 181112
        All of the data comes from Prochaska+2019, Science, in press
    Returns:

    """
    frb181112 = frb.FRB('FRB181112', 'J214923.63-525815.33',
                        589.27 * units.pc / units.cm**3,
                        z_frb=0.4755)
    # Error in DM
    frb181112.DM_err = 0.03 * units.pc / units.cm**3
    # Error ellipse
    frb181112.set_ee(a=555.30/1e3, b=152.93/1e3, theta=120.15, cl=68.)

    # RM
    frb181112.RM = 10.9 * units.rad / units.m**2
    frb181112.RM_err = 0.9 * units.rad / units.m**2
    frb181112.fluence = 26. * units.Jy * units.ms
    frb181112.fluence_err = 3. * units.Jy * units.ms    

    # NE2001
    frb181112.set_DMISM()

    # References
    frb181112.refs = ['Prochaska2019']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb181112.write_to_json(path=path)

def frb_190102():
    """Bhandari+20, ApJL --
    Now including updated values on FRB position, DM, pulse width 
    and scattering time from Day et al. 2020
    """
    fname = 'FRB190102'
    #wv_oiii = 6466.48
    #z_OIII = wv_oiii / 5008.239 - 1
    frb190102 = frb.FRB(fname, 'J212939.76-792832.5', # Day+20
                        364.545 * units.pc / units.cm**3,
                        z_frb=0.2912) # Updated redshift
    # Error ellipse [REQUIRED]
    frb190102.set_ee(0.1, 0.1, 0., 95.) #statistical
    frb190102.set_ee(0.4, 0.5, 0., 95., stat=False) #systematic

    # Error in DM
    frb190102.DM_err = 0.004 * units.pc / units.cm**3

    # NE2001
    frb190102.set_DMISM()
    
     # FRB properties
    frb190102.fluence = 14 * units.Jy * units.ms    # From Bhandari+20
    frb190102.fluence_err = 1 * units.Jy * units.ms # -- //-- 
    frb190102.RM = -105 * units.rad / units.m**2
    frb190102.RM_err = 1 * units.rad / units.m**2
    frb190102.tau = 0.041 * units.ms
    frb190102.tau_err = 0.003 * units.ms
    
    # References
    frb190102.refs = ['Bhandari2020', 'Day2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190102.write_to_json(path=path)


def frb_190523():
    """Ravi+19, Nature --
        https://ui.adsabs.harvard.edu/abs/2019Natur.572..352R/abstract
    """
    fname = 'FRB190523'
    frb190523 = frb.FRB(fname, 'J134815.6+722811',
                        760.8 * units.pc / units.cm**3,
                        z_frb=0.660)
    # Error ellipse [REQUIRED]
    frb190523.set_ee(4, 1.5, 0., 340) # Halfing the 3x8'' at 95% reported by Ravi et al.
    # Error in DM
    frb190523.DM_err = 0.6 * units.pc / units.cm**3

    # Fluence
    frb190523.fluence = 280 * units.Jy * units.ms
    #frb180924.fluence_err = 1 * units.Jy * units.ms

    # NE2001
    frb190523.set_DMISM()

    # RM
    #frb190102.RM = 10 * units.rad / units.m**2
    #frb190102.RM_err = 1 * units.rad / units.m**2

    # References
    frb190523.refs = ['Ravi2019']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190523.write_to_json(path=path)

def frb_190608():
    """Bhandari+20, ApJL,
    Day+20  https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract
    Macquart al. Nature
    """
    fname = 'FRB190608'
    frb190608 = frb.FRB(fname, "J221604.77-075353.7",  # Pulled from Slack on 2020 Mar 18
                        340.05 * units.pc / units.cm**3,
                        z_frb=0.1177805)  # Taken from the SDSS table
    # Error ellipse [REQUIRED]
    frb190608.set_ee(0.19315, 0.18, 0., 68.) # Statistsical
    frb190608.set_ee(0.178292, 0.18, 0., 68., stat=False)  # Systematic
    # Error in DM
    frb190608.DM_err = 0.6 * units.pc / units.cm**3

    # NE2001
    frb190608.set_DMISM()

    # FRB properties
    frb190608.fluence = 26 * units.Jy * units.ms    # From Bhandari+20
    frb190608.fluence_err = 4 * units.Jy * units.ms # -- //-- 
    frb190608.RM = 353 * units.rad / units.m**2
    frb190608.RM_err = 2 * units.rad / units.m**2
    frb190608.tau = 3.3 * units.ms
    frb190608.tau_err = 0.2 * units.ms

    # References
    frb190608.refs = ['Bhandari2020','Day2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190608.write_to_json(path=path)


def frb_190611():
    """
    Macquart al. Nature
    Day et al. 2020  https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract
    """
    FRB_190611_coord = SkyCoord('J212258.91-792351.3',  # Day+2020
                                unit=(units.hourangle, units.deg))
    frb190611 = frb.FRB('FRB190611', FRB_190611_coord,
                        332.63 * units.pc / units.cm**3,
                        z_frb=0.3778)  # Bright
    # Error ellipse [REQUIRED]
    frb190611.set_ee(0.7, 0.7, 0., 68.) #statistical
    frb190611.set_ee(0.4, 0.3, 0., 68, stat=False ) #systematic
    # Error in DM
    frb190611.DM_err = 0.04 * units.pc / units.cm**3

    # NE2001
    frb190611.set_DMISM()
    
    # FRB properties
    frb190611.RM = 20 * units.rad / units.m**2
    frb190611.RM_err = 4 * units.rad / units.m**2
    frb190611.tau = 0.18 * units.ms
    frb190611.tau_err = 0.02 * units.ms
    frb190611.fluence = 10 * units.Jy * units.ms
    frb190611.fluence_err = 2 * units.Jy * units.ms    

    # References
    frb190611.refs = ['MacQuart2020', 'Day2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190611.write_to_json(path=path)

def frb_190614():
    """
    Law+2020 https://ui.adsabs.harvard.edu/abs/2020ApJ...899..161L/abstract

    """
    fname = 'FRB190614'
    FRB_190614_coord = SkyCoord(ra=65.07552, dec=73.70674, unit='deg')
    frb190614 = frb.FRB(fname, FRB_190614_coord,
                        959 * units.pc / units.cm ** 3)
    # Error ellipse [REQUIRED]
    frb190614.set_ee(a=0.8, b=0.4, theta=67., cl=68.)
    # Error in DM
    frb190614.DM_err = 1 * units.pc / units.cm ** 3

    # NE2001
    frb190614.set_DMISM()
    # RM
    # frb190102.RM = 10 * units.rad / units.m**2
    # frb190102.RM_err = 1 * units.rad / units.m**2

    # References
    frb190614.refs = ['Law2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190614.write_to_json(path=path)


def frb_190711():
    """MacQuart+20, Day+2020
      https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract
    """
    fname = 'FRB190711'
    frb190711 = frb.FRB(fname, 'J215740.68-802128.8',  # MacQuarter+2020, Day+2020
                        587.9 * units.pc / units.cm ** 3,    # Day+2020
                        z_frb=0.52172)
    # Error ellipse
    frb190711.set_ee(0.3, 0.3, 0., 68.)  # (Statistical)
    frb190711.set_ee(0.38, 0.3, 0., 68., stat=False)  # Systematic

    # Error in DM
    frb190711.DM_err = 1 * units.pc / units.cm ** 3

    # NE2001
    frb190711.set_DMISM()
    # RM -- Day+2020
    frb190711.RM = 9 * units.rad / units.m**2
    frb190711.RM_err = 2 * units.rad / units.m**2
    frb190711.fluence = 34 * units.Jy * units.ms
    frb190711.fluence_err = 3 * units.Jy * units.ms

    # References
    frb190711.refs = ['MacQuart2020', 'Day2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190711.write_to_json(path=path)


def frb_190714():
    """
    Day+2020 https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract
    Heintz+2020


    """
    # Taken from Slack; Cherie posted on 19 June 2020
    FRB_190714_coord = SkyCoord('J121555.12-130115.7',
                                unit=(units.hourangle, units.deg))
    frb190714 = frb.FRB('FRB190714', FRB_190714_coord,
                        504.13 * units.pc / units.cm ** 3,
                        z_frb=0.2365)
    # Error ellipse
    frb190714.set_ee(0.011, 0.1, 0., 68.)
    frb190714.set_ee(0.022, 0.2, 0., 68., stat=False)  # Systematic

    # Error in DM
    frb190714.DM_err = 0.1 * units.pc / units.cm ** 3

    # NE2001
    frb190714.set_DMISM()

    # RM

    # References
    frb190714.refs = ['Heintz2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190714.write_to_json(path=path)


def frb_191001():
    """
    Bhandari+2020b  https://ui.adsabs.harvard.edu/abs/2020arXiv200812488B/abstract

    """
    # Taken from Bhandari, Table 1
    FRB_191001_coord = SkyCoord("21h33m24.373s -54d44m51.86s", frame='icrs')
    frb191001 = frb.FRB('FRB191001', FRB_191001_coord,
                        507.90 * units.pc / units.cm ** 3,
                        z_frb=0.2340)
    # Error ellipse [REQUIRED]
    frb191001.set_ee(0.006, 0.13, 90., 68.)  # This is statistical + systematic in quadrature
    # Error in DM
    frb191001.DM_err = 0.07 * units.pc / units.cm ** 3

    # NE2001
    frb191001.set_DMISM()
    # RM
    frb191001.RM = 54.86 * units.rad / units.m**2
    frb191001.RM_err = 0.50 * units.rad / units.m**2
    # Fluence
    frb191001.fluence = 143 * units.Jy * units.ms
    frb191001.fluence_err = 15 * units.Jy * units.ms
    # Spectral energy density

    # Scattering time at 1 GHz
    frb191001.tau = 2.19 * units.ms
    frb191001.tau_err = 0.23 *units.ms

    # References
    frb191001.refs = ['Bhandari2020b']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb191001.write_to_json(path=path)



def main(inflg='all'):

    if inflg == 'all':
        flg = np.sum(np.array( [2**ii for ii in range(25)]))
    else:
        flg = int(inflg)

    # 121102
    if flg & (2**0):
        frb_121102()

    # 180924
    if flg & (2**1):  # 2
        frb_180924()

    # 181112
    if flg & (2**2): #4
        frb_181112()

    # 195023
    if flg & (2**3): #8
        frb_190523()

    # 190608
    if flg & (2**4): #16
        frb_190608()

    # 190102
    if flg & (2**5): #32
        frb_190102()

    # 190711
    if flg & (2**6): # 64
        frb_190711()

    # 180916
    if flg & (2**7): # 128
        frb_180916()

    # 190611
    if flg & (2**8):  # 256
        frb_190611()

    # FRB 190614      # 512
    if flg & (2**9):
        frb_190614()

    # FRB 190714      # 1024
    if flg & (2**10):
        frb_190714()

    # FRB 191001      # 2048
    if flg & (2**11):
        frb_191001()


# Command line execution
#  Only for testing
#  Use the Build script to build
if __name__ == '__main__':
    # FRB 121102
    frb_121102()



