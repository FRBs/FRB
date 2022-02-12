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
        or Chatterjee et al. 2017
    Update for Marcote et al. 2017
        05h31m58.7013s +/- 4 mas
        +33d08m52.5536s +/- 2.3 mas
    """
    frb121102 = frb.FRB('FRB20121102', 'J053158.7013+330852.5536',
                        558.1*units.pc/units.cm**3,
                        z_frb=0.19273, repeater=True)
    # NE2001
    frb121102.set_DMISM()

    # RM
    frb121102.RM = 1.e5 * units.rad / units.m**2

    # Pulse properties
    frb121102.set_pulse(1*units.GHz,
        Wi=3.0*units.ms,
        Wi_err=0.5*units.ms,
        tscatt=0.024*units.ms) # no errors given

    # Error ellipse
    frb121102.set_ee(0.004, 0.002, theta=90., cl=68.)
    frb121102.set_ee(a=0.0015, b=0.0015, theta=0., cl=68.,stat=False) # Marcote priv. corr

    # References
    frb121102.refs = ['Tendulkar2017', 'Marcote2017']
    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb121102.write_to_json(path=path)
    # Test
    #frb121102.from_json('FRB121102.json')


def frb_180301():
    """
    Bhandari+2021
    """
    frbname = 'FRB20180301'

    FRB_180301_coord = SkyCoord("06h12m54.44s +04d40m15.8s", frame='icrs') # Bhandari+2021
    frb180301 = frb.FRB(frbname, FRB_180301_coord,
                        536 * units.pc / units.cm ** 3,
                        z_frb=0.33044, repeater=True)  # Slack posting
    
    # Error ellipse (Statistical) # 
    frb180301.set_ee(0.01142,0.00825, 0., 68.)
    # Error ellipse (Systematic)
    frb180301.set_ee(0.619, 0.603, 0., 68., stat=False)


    # Error in DM
    frb180301.DM_err = 13 * units.pc / units.cm ** 3

    # NE2001
    frb180301.set_DMISM()
    # RM
    # frb190102.RM = 10 * units.rad / units.m**2

    # References
    frb180301.refs = ['Bhandari2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb180301.write_to_json(path=path)

def frb_180916():
    """
     FRB 180916.J0158+65
        All of the data currently comes from Marcote et al. 2020
        https://ui.adsabs.harvard.edu/abs/2020Natur.577..190M/abstract
    """
    coord = SkyCoord("01h58m00.75017s 65d43m00.3152s", frame='icrs')
    name = 'FRB20180916'
    frb180916 = frb.FRB(name, coord,
                        348.76*units.pc / units.cm**3,
                        z_frb=0.0337, repeater=True)
    # Error ellipse
    frb180916.set_ee(0.0011, 0.0011, theta=0., cl=68.)
    frb180916.set_ee(a=0.002, b=0.002, theta=0., cl=68.,stat=False)  # Marcote priv. corr
    # Error in DM
    frb180916.DM_err = 0.10 * units.pc / units.cm**3
    # NE2001
    frb180916.set_DMISM()
   
    # RM and fluence
    frb180916.fluence = 2.53*units.Jy*units.ms # Brightest
    frb180916.RM = -114.6 * units.rad / units.m**2 # From CHIME/FRB 2019
    frb180916.RM_err = 0.6 * units.rad / units.m**2
    
    # Pulse properties
    frb180916.set_pulse(1*units.GHz,
        Wi=1.66*units.ms,
        Wi_err=0.05*units.ms, 
        tscatt=0.0027*units.ms) # no errors given
    

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

        Uncertinties in FRB localization updated by Day et al. 2021
    """
    frb180924 = frb.FRB('FRB20180924', 'J214425.255-405400.10',
                        362.16*units.pc / units.cm**3,
                        z_frb=0.3212, repeater=False)

    # Error ellipse 
    frb180924.set_ee(a=0.07, b=0.07, theta=0., cl=68.) # Statistical (Day+2020)
    frb180924.set_ee(a=0.1611, b=0.1611, theta=0., cl=68.,stat=False) # Systematic (Day+2021)

    # Error in DM (Day 2019)
    frb180924.DM_err = 0.01 * units.pc / units.cm**3

    # NE2001
    frb180924.set_DMISM()

    # RM, fluence and polarization
    frb180924.fluence = 16 * units.Jy * units.ms    # From Bhandari+20
    frb180924.fluence_err = 1 * units.Jy * units.ms # -- //-- 
    frb180924.RM = 22 * units.rad / units.m**2
    frb180924.RM_err = 2 * units.rad / units.m**2
    frb180924.lpol = 80.  # %
    frb180924.lpol_err = 10.

    # Pulse properties
    frb180924.set_pulse(1.2725*units.GHz,
        Wi=0.09*units.ms,
        Wi_err=0.04*units.ms,
        tscatt=0.68*units.ms,
        tscatt_err=0.03*units.ms)
 

    # References
    frb180924.refs = ['Bannister2019', 'Day2020', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb180924.write_to_json(path=path)

def frb_181112():
    """
    Generate the JSON file for FRB 181112
        All of the data comes from Prochaska+2019, Science, in press
    Returns:

    FRB position and localization updated by Day+2021

    """
    frb181112 = frb.FRB('FRB20181112', 
                        'J214923.63-525815.4',
                        589.27 * units.pc / units.cm**3,
                        z_frb=0.4755, repeater=False)
    # Error in DM
    frb181112.DM_err = 0.03 * units.pc / units.cm**3
    # Error ellipse -- Updated by Day+2021
    frb181112.set_ee(a=555.30/1e3, b=152.93/1e3, theta=120.15, cl=68.) # Statistical (Prochaska+2019)
    frb181112.set_ee(a=3.2*1.79, b=0.8*1.79, theta=120.15, cl=68.,stat=False) # Systematic

    # RM and fluence
    frb181112.RM = 10.5 * units.rad / units.m**2
    frb181112.RM_err = 0.4 * units.rad / units.m**2
    frb181112.fluence = 20.2 * units.Jy * units.ms
    frb181112.fluence_err = 0.1 * units.Jy * units.ms

    # Pulse properties
    frb181112.set_pulse(1.2725*units.GHz,
        Wi=0.016*units.ms,
        Wi_err=0.001*units.ms,
        tscatt=0.021*units.ms,
        tscatt_err=0.001*units.ms)

    # NE2001
    frb181112.set_DMISM()

    # References
    frb181112.refs = ['Prochaska2019', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb181112.write_to_json(path=path)

def frb_190102():
    """Bhandari+20, ApJL --
    Now including updated values on FRB position, DM, pulse width 
    and scattering time from Day et al. 2020
    Update systematic uncertainty in FRB localization from Day+2021
    """
    fname = 'FRB20190102'
    #wv_oiii = 6466.48
    #z_OIII = wv_oiii / 5008.239 - 1
    frb190102 = frb.FRB(fname, 'J212939.76-792832.5', # Day+20
                        364.545 * units.pc / units.cm**3,
                        z_frb=0.2912, repeater=False) # Updated redshift
    # Error ellipse [REQUIRED]
    frb190102.set_ee(0.21, 0.17, theta=0., cl=68.) # Statistical (Day+2020)
    frb190102.set_ee(0.936, 0.7876, theta=0., cl=68., stat=False) # Systematic (Day+2021)

    # Error in DM
    frb190102.DM_err = 0.004 * units.pc / units.cm**3

    # NE2001
    frb190102.set_DMISM()
    
     # RM and fluence
    frb190102.fluence = 14 * units.Jy * units.ms    # From Bhandari+20
    frb190102.fluence_err = 1 * units.Jy * units.ms # -- //-- 
    frb190102.RM = -105 * units.rad / units.m**2
    frb190102.RM_err = 1 * units.rad / units.m**2

    # Pulse properties
    frb190102.set_pulse(1.2725*units.GHz,
        Wi=0.053*units.ms,
        Wi_err=0.002*units.ms,
        tscatt=0.041*units.ms,
        tscatt_err=0.003*units.ms)   
 
    # References
    frb190102.refs = ['Bhandari2020', 'Day2020', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190102.write_to_json(path=path)


def frb_190523():
    """Ravi+19, Nature --
        https://ui.adsabs.harvard.edu/abs/2019Natur.572..352R/abstract
    """
    fname = 'FRB20190523'
    frb190523 = frb.FRB(fname, 'J134815.6+722811',
                        760.8 * units.pc / units.cm**3,
                        z_frb=0.660, repeater=False)
    # Error ellipse [REQUIRED]
    frb190523.set_ee(4, 1.5, cl=68., theta=340.) # Halfing the 3x8'' at 95% reported by Ravi et al.
    # Error in DM
    frb190523.DM_err = 0.6 * units.pc / units.cm**3

    # FRB properties
    frb190523.fluence = 280 * units.Jy * units.ms
    #frb180924.fluence_err = 1 * units.Jy * units.ms
    frb190523.tau = 1.4 * units.ms
    frb190523.tau_err = 0.2 * units.ms

    # NE2001
    frb190523.set_DMISM()

    # Pulse properties
    frb190523.set_pulse(1.*units.GHz,
        Wi=0.42*units.ms,
        Wi_err=0.05*units.ms,
        tscatt=1.4*units.ms,
        tscatt_err=0.2*units.ms)

    # References
    frb190523.refs = ['Ravi2019']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190523.write_to_json(path=path)


def frb_190608():
    """Bhandari+20, ApJL,
    Day+20  https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract
    Macquart al. Nature

    Update systematic uncertainty in FRB localization from Day+2021
    """
    fname = 'FRB20190608'
    frb190608 = frb.FRB(fname, "J221604.77-075353.7",  # Pulled from Slack on 2020 Mar 18
                        340.05 * units.pc / units.cm**3,
                        z_frb=0.1177805, repeater=False)  # Taken from the SDSS table
    # Error ellipse [REQUIRED]
    frb190608.set_ee(0.19, 0.18, theta=90., cl=68.) # Statistsical (Day+2020)
    frb190608.set_ee(0.33, 0.30, theta=90., cl=68., stat=False) # Systematic
    
    # Error in DM
    frb190608.DM_err = 0.6 * units.pc / units.cm**3

    # NE2001
    frb190608.set_DMISM()

    # RM and fluence 
    frb190608.fluence = 26 * units.Jy * units.ms    # From Macquart+20
    frb190608.fluence_err = 4 * units.Jy * units.ms # -- //-- 
    frb190608.RM = 353 * units.rad / units.m**2 # From Day+20
    frb190608.RM_err = 2 * units.rad / units.m**2

    # Pulse properties
    frb190608.set_pulse(1.2725*units.GHz,
        Wi=1.1*units.ms,
        Wi_err=0.2*units.ms,
        tscatt=3.3*units.ms,
        tscatt_err=0.2*units.ms)

    # References
    frb190608.refs = ['Bhandari2020','Day2020', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190608.write_to_json(path=path)


def frb_190611():
    """
    Macquart al. Nature
    Day et al. 2020  https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract

    Update systematic uncertainty in FRB localization from Day+2021
    """
    FRB_190611_coord = SkyCoord('J212258.94-792351.3',  # Day+2021
                                unit=(units.hourangle, units.deg))
    frb190611 = frb.FRB('FRB20190611', FRB_190611_coord,
                        332.63 * units.pc / units.cm**3,
                        z_frb=0.3778, repeater=False)  # Bright
    # Error ellipse 
    frb190611.set_ee(0.34, 0.32, theta=0., cl=68.) # Statistical (Day+2020)
    frb190611.set_ee(1.12, 1.07, theta=0., cl=68., stat=False) # Systematic (Day+2021)

    # Error in DM
    frb190611.DM_err = 0.04 * units.pc / units.cm**3

    # NE2001
    frb190611.set_DMISM()
    
    # RM and fluence
    frb190611.RM = 20 * units.rad / units.m**2 # From Day+20
    frb190611.RM_err = 4 * units.rad / units.m**2
    frb190611.fluence = 10 * units.Jy * units.ms # From Macquart+20
    frb190611.fluence_err = 2 * units.Jy * units.ms    

    # Pulse properties
    frb190611.set_pulse(1.2725*units.GHz,
        Wi=0.09*units.ms,
        Wi_err=0.02*units.ms,
        tscatt=0.18*units.ms,
        tscatt_err=0.02*units.ms)

    # References
    frb190611.refs = ['MacQuart2020', 'Day2020', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190611.write_to_json(path=path)

def frb_190614():
    """
    Law+2020 https://ui.adsabs.harvard.edu/abs/2020ApJ...899..161L/abstract

    """
    fname = 'FRB20190614'
    FRB_190614_coord = SkyCoord(ra=65.07552, dec=73.70674, unit='deg')
    frb190614 = frb.FRB(fname, FRB_190614_coord,
                        959 * units.pc / units.cm ** 3, repeater=False)
    # Error ellipse [REQUIRED]
    frb190614.set_ee(a=0.8, b=0.4, theta=67., cl=68.)
    # Error in DM
    frb190614.DM_err = 1 * units.pc / units.cm ** 3

    # NE2001
    frb190614.set_DMISM()
    
    # FRB properties
    frb190614.fluence = 0.62 * units.Jy * units.ms
    frb190614.fluence_err = 0.07 * units.Jy * units.ms

    # References
    frb190614.refs = ['Law2020']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190614.write_to_json(path=path)


def frb_190711():
    """MacQuart+20, Day+2020
      https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract

    Updated FRB localization + error from Day+2021
    """
    fname = 'FRB20190711'
    frb190711 = frb.FRB(fname,'J215740.62-802128.8',  # MacQuarter+2020, Day+2020
                        587.9 * units.pc / units.cm ** 3,    # Day+2020
                        z_frb=0.52172, repeater=True)
    # Error ellipse
    frb190711.set_ee(0.12, 0.075, theta=90., cl=68.)  # Statistical (Day+2021)
    frb190711.set_ee(0.646, 0.563, theta=90., cl=68., stat=False) # Systematic (Day+2021)

    # Error in DM
    frb190711.DM_err = 1 * units.pc / units.cm ** 3

    # NE2001
    frb190711.set_DMISM()
    # RM and fluence -- Day+2020
    frb190711.RM = 9 * units.rad / units.m**2 # Day+20
    frb190711.RM_err = 2 * units.rad / units.m**2
    frb190711.fluence = 34 * units.Jy * units.ms # Macquart+20
    frb190711.fluence_err = 3 * units.Jy * units.ms

    # References
    frb190711.refs = ['MacQuart2020', 'Day2020', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190711.write_to_json(path=path)


def frb_190714():
    """
    Day+2020 https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3335D/abstract
    Heintz+2020

    Updated FRB localization + error from Day+2021
    """
    # Taken from Slack; Cherie posted on 19 June 2020
    # Updated as per Day 2021
    FRB_190714_coord = SkyCoord('J121555.13-130115.6',
                                unit=(units.hourangle, units.deg))
    frb190714 = frb.FRB('FRB20190714', FRB_190714_coord,
                        504.13 * units.pc / units.cm ** 3,
                        z_frb=0.2365, repeater=False)

    # Error ellipse (Day+2021)
    frb190714.set_ee(0.17, 0.10, theta=90., cl=68.) # Statistical
    frb190714.set_ee(0.5191, 0.376, theta=90., cl=68., stat=False)  # Systematic

    # Error in DM
    frb190714.DM_err = 0.1 * units.pc / units.cm ** 3

    # NE2001
    frb190714.set_DMISM()

    # Fluence
    frb190714.fluence = 12 * units.Jy * units.ms # From Cherie on Slack (2/10 - 2020)
    frb190714.fluence_err = 2 * units.Jy * units.ms

    # References
    frb190714.refs = ['Heintz2020', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb190714.write_to_json(path=path)


def frb_191001():
    """
    Bhandari+2020b  https://ui.adsabs.harvard.edu/abs/2020arXiv200812488B/abstract

    Updated FRB localization + error from Day+2021
    """
    # Taken from Bhandari, Table 1
    # Updated Day 2021
    FRB_191001_coord = SkyCoord("21h33m24.313s -54d44m51.86s", frame='icrs')
    frb191001 = frb.FRB('FRB20191001', FRB_191001_coord,
                        507.90 * units.pc / units.cm ** 3,
                        z_frb=0.2340, repeater=False)
    # Error ellipse [REQUIRED]
    frb191001.set_ee(0.13, 0.08, theta=90., cl=68.) # Statistical  -- Day+2021
    frb191001.set_ee(0.1737, 0.160, theta=90., cl=68., stat=False)  # Systematic
    # Error in DM
    frb191001.DM_err = 0.07 * units.pc / units.cm ** 3

    # NE2001
    frb191001.set_DMISM()
    # RM and fluence (Bhandari+20b)
    frb191001.RM = 55.5 * units.rad / units.m**2
    frb191001.RM_err = 0.9 * units.rad / units.m**2
    frb191001.fluence = 143 * units.Jy * units.ms
    frb191001.fluence_err = 15 * units.Jy * units.ms

    # Pulse properties
    frb191001.set_pulse(0.920*units.GHz,
        Wi=0.22*units.ms,
        Wi_err=0.03*units.ms,
        tscatt=3.3*units.ms,
        tscatt_err=0.2*units.ms)

    # References
    frb191001.refs = ['Bhandari2020b', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb191001.write_to_json(path=path)


def frb_191228():
    """
    Bhandari+2021

    Updated FRB localization + error from Day+2021
    """
    frbname = 'FRB20191228'
    FRB_191228_coord = SkyCoord('22h57m43.30s -29d35m38.7s',  frame='icrs') # Taken from Slack on 11 Oct 2020
    frb191228 = frb.FRB(frbname, FRB_191228_coord,
                        298 * units.pc / units.cm ** 3,
                        repeater=False)  # First Slack posting
    # Error ellipse : Day+2021
    frb191228.set_ee(0.34, 0.34, 0., 68.)
    frb191228.set_ee(0.830, 0.823, 0., 68., stat=False)
    # Error in DM
    frb191228.DM_err = 0.05 * units.pc / units.cm ** 3

    # NE2001
    frb191228.set_DMISM()
    # RM
    # frb190102.RM = 10 * units.rad / units.m**2
    # frb190102.RM_err = 1 * units.rad / units.m**2

    # References
    frb191228.refs = ['Bhandari2021', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb191228.write_to_json(path=path)

def frb_20200120E():
    """M81 + globular cluster
    """
    frbname = 'FRB20200120E'

    FRB_20200120E_coord = SkyCoord("09h57m54.68s  +68d49m08.0s", frame='icrs')    # From EVN localisation
    frb20200120E = frb.FRB(frbname, FRB_20200120E_coord,
                        87.818 * units.pc / units.cm ** 3,z_frb=0.0008,repeater=True)  
    # Error ellipse (Statistical)
    frb20200120E.set_ee(0.4,0.4, 0., 68.)
    # Error ellipse (Systematic)
    frb20200120E.set_ee(0., 0., 0., 68., stat=False)

    # Error in DM
    frb20200120E.DM_err = 0.007 * units.pc / units.cm ** 3

    # RM
    frb20200120E.RM = -29.8 * units.rad / units.m**2
    frb20200120E.RM_err = 0.5 * units.rad / units.m**2

    # NE2001
    frb20200120E.set_DMISM()
    # RM
    # frb190102.RM = 10 * units.rad / units.m**2
    # frb190102.RM_err = 1 * units.rad / units.m**2

    # References
    frb20200120E.refs = ['Bhardwaj2021','Kirsten2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb20200120E.write_to_json(path=path)

def frb_171020():
    frbname = 'FRB20171020'

    FRB_171020_coord = SkyCoord("22h15m18.55s  -19d40m11.23s", frame='icrs')    # From Shannon+18
    frb171020 = frb.FRB(frbname, FRB_171020_coord,
                        114.1 * units.pc / units.cm ** 3,z_frb=0.00867,repeater=False)  
    # Error ellipse (Statistical)
    frb171020.set_ee(600,600, 0., 68.)
    # Error ellipse (Systematic)
    frb171020.set_ee(0., 0., 0., 68., stat=False)

    # Error in DM
    frb171020.DM_err = 0.2 * units.pc / units.cm ** 3

    # NE2001
    frb171020.set_DMISM()
    # RM
    # frb190102.RM = 10 * units.rad / units.m**2
    # frb190102.RM_err = 1 * units.rad / units.m**2

    # References
    frb171020.refs = ['Shannon2018','Mahony2018']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb171020.write_to_json(path=path)


def frb_200430():
    """
    Heintz+2020
    Bhandari+2021

    Updated FRB localization + error from Day+2021
    """
    frbname = 'FRB20200430'

    FRB_200430_coord = SkyCoord('15h18m49.54s +12d22m36.3s', frame='icrs')
    frb200430 = frb.FRB(frbname, FRB_200430_coord,
                        380 * units.pc / units.cm ** 3,  # First Slack posting
                        z_frb = 0.161,
                        repeater=False)
    # Error ellipse (Statistical)
    frb200430.set_ee(0.24, 0.17, 0., 68.)
    # Error ellipse (Systematic)
    frb200430.set_ee(0.98, 0.251, 0., 68., stat=False)

    # Error in DM
    #frb191001.DM_err = 1 * units.pc / units.cm ** 3

    # NE2001
    frb200430.set_DMISM()
    # RM
    # frb190102.RM = 10 * units.rad / units.m**2
    # frb190102.RM_err = 1 * units.rad / units.m**2

    # References
    frb200430.refs = ['Bhandari2021', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb200430.write_to_json(path=path)


def frb_200906():
    """
    Bhandari+2021

    Updated FRB localization + error from Day+2021
    """
    frbname = 'FRB20200906'

    FRB_200906_coord = SkyCoord("03h33m59.08s -14d04m59.46s", frame='icrs')  # Day+2021
    frb200906 = frb.FRB(frbname, FRB_200906_coord,
                        577.84 * units.pc / units.cm**3,
                        z_frb=0.36879,
                        repeater=False)  # Slack posting
    # Error ellipse (Statistical)
    frb200906.set_ee(0.11, 0.102, 0., 68.)
    # Error ellipse (Systematic)
    frb200906.set_ee(0.55, 0.340, 0., 68., stat=False)

    # Error in DM
    frb200906.DM_err = 0.02 * units.pc / units.cm ** 3

    # NE2001
    frb200906.set_DMISM()
    # RM
    # frb190102.RM = 10 * units.rad / units.m**2
    # frb190102.RM_err = 1 * units.rad / units.m**2

    # References
    frb200906.refs = ['Bhandari2021', 'Day2021']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb200906.write_to_json(path=path)

def frb_201124():
    """
    ATELs only so far

    """
    # ATEL 14603 (VLBI)
    FRB_201124_coord = SkyCoord("05h08m03.5077s 26d03m38.504s", frame='icrs')
    frb201124 = frb.FRB('FRB20201124', FRB_201124_coord,
                        411. * units.pc / units.cm ** 3,
                        z_frb=0.0982, repeater=True)

    # Error ellipse [REQUIRED]
    frb201124.set_ee(0.004, 0.004, theta=0., cl=68.) # ATEL

    # Error in DM
    #frb191001.DM_err = 0.07 * units.pc / units.cm ** 3

    # NE2001
    frb201124.set_DMISM()

    # RM (Kumar+21)
    frb201124.RM = -613 * units.rad / units.m**2
    frb201124.RM_err = 2 * units.rad / units.m**2

    # RM and fluence (Bhandari+20b)
    #frb201124.RM = 55.5 * units.rad / units.m**2
    #frb201124.RM_err = 0.9 * units.rad / units.m**2
    #frb201124.fluence = 143 * units.Jy * units.ms
    #frb201124.fluence_err = 15 * units.Jy * units.ms

    # Pulse properties
    #frb191001.set_pulse(0.920*units.GHz,
    #    Wi=0.22*units.ms,
    #    Wi_err=0.03*units.ms,
    #    tscatt=3.3*units.ms,
    #    tscatt_err=0.2*units.ms)

    # References
    #frb191001.refs = ['Bhandari2020b']

    # Write
    path = resource_filename('frb', 'data/FRBs')
    frb201124.write_to_json(path=path)
    


def main(inflg='all'):

    if inflg == 'all':
        flg = np.sum(np.array( [2**ii for ii in range(25)]))
    else:
        flg = int(inflg)

    # 121102
    if flg & (2**0):  # 1
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

    # FRB 201124a      # 4096
    if flg & (2**12):
        frb_201124()

    # FRB 180301      # 8192
    if flg & (2**13):
        frb_180301()

    # FRB 191228      # 16384
    if flg & (2**14):
        frb_191228()

    # FRB 200430      # 32768
    if flg & (2**15):
        frb_200430()

    # FRB 200906      # 65536
    if flg & (2**16):
        frb_200906()

    # FRB 200906      
    if flg & (2**17):  # 131072
        frb_20200120E()

    # FRB 20171020
    if flg & (2**18):
        frb_171020()


# Command line execution
#  Only for testing
#  Use the Build script to build
if __name__ == '__main__':
    # FRB 121102
    #frb_190611()
    #frb_190711()
    #frb_190714()
    #frb_191001()
    #frb_180301()
    #frb_200906()
    #frb_191228()
    #frb_171020()
    frb_20200120E()
    #frb_201124()
