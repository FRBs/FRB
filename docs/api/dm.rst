frb.dm - Dispersion Measure Module
===================================

The ``frb.dm`` package provides tools for estimating and modeling DM
contributions from the IGM, host galaxy, and related components.

.. automodule:: frb.dm
   :members:
   :undoc-members:
   :show-inheritance:

Submodules
----------

.. toctree::
   :maxdepth: 1

   dm.cosmic
   dm.dm_ism_healpix_map
   dm.host
   dm.igm
   dm.mcmc
   dm.prob_dmz

.. note::
   Optional dependencies: ``pymc3`` is required for MCMC workflows in
   ``frb.dm.mcmc``.
