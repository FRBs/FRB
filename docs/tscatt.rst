.. highlight:: rest

********************
Turbulent Scattering
********************

.. _Macquart&Koay13: http://adsabs.harvard.edu/abs/2013ApJ...776..125M

The formalism developed by `Macquart&Koay13`_ has been ingested
into the *Turbulence* class for calculations of temporal
smearing and angular broadening.

.. _turbulence:

Turbulence
==========

This class describes the turbulence from a single
object with properties listed in the following Table:

========== =============== =================================== ====================
Parameter  Unit            Description                         Required?
========== =============== =================================== ====================
ne         density; cm**-3 Fiducial value of electron density  Yes
l0         length; AU      Inner length scale                  Yes
L0         length; pc      Outer length scale                  Yes
zL         unitless        Redshift of scattering medium       Yes
beta       unitless        Exponent of turbulence              Yes
DL         length; kpc     Thickness of the object             Required for SM
lobs       length; cm      Wavelength of interest              Input for tau, theta
========== =============== =================================== ====================


Here is a simple instantiation::

    def_l0 = 1*u.AU
    def_L0 = 0.001 * u.pc
    def_ne = 1e-2 / u.cm**3
    def_DL = 1 * u.kpc
    def_zL = 1.
    turb = Turbulence(def_ne, def_l0, def_L0)

The resultant object can then be used for scattering
calculations.  Note that in this instantiation,
beta has a default value of 11/3 (Kolmogorov) and
zL=0.

Temporal Smearing
-----------------

With a :ref:`turbulence` object defined as above,
one can calculate the temporal smearing of a source
at a given redshift and for a given observed wavelength
(or frequency).  Here is an example::

    from astropy import constants as const
    turb = Turbulence(def_ne, def_l0, def_L0, def_zL)
    inu = 1 * u.GHz
    lobs = (const.c/inu).to('cm')
    # SM, rdiff
    turb.set_SM_obj(def_DL)
    turb.set_rdiff(lobs)
    # tau
    zsource = 2.
    tau = turb.temporal_smearing(lobs, zsource)
    print(tau)  # Given in ms

Angular Broadening
------------------

With a :ref:`turbulence` object defined as above,
one can calculate the angular broadening of a source
at a given redshift and for a given observed wavelength
(or frequency).  Here is an example::

    theta = turb.angular_broadening(lobs, zsource) # Returned in arcsec


Intervening Galaxies
--------------------

One may generate a set of random draws for temporal smearing
due to the ISM of intervening galaxies towards a source with given
redshift using the *monte_tau* method::

    from frb.dlas import monte_tau
    zeval = np.array([0.,1.,2.])
    taus = monte_tau(np.array(zeval))  # Returned unitless but in ms

