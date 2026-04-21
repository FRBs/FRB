Halo Mass Function
==================

.. currentmodule:: frb.halos.hmf

.. automodule:: frb.halos.hmf
   :members:
   :undoc-members:
   :show-inheritance:

.. warning::
   This module is **deprecated**. Much of this module's functionality has been moved to :mod:`frb.halos.models`.
   
   **Migration Guide:**
   
   - ``init_hmf()`` → :func:`frb.halos.models.init_hmf`
   - ``frac_in_halos()`` → :func:`frb.halos.models.frac_in_halos`
   - ``halo_incidence()`` → :func:`frb.halos.models.halo_incidence`
   - ``build_grid()`` → :func:`frb.halos.models.build_grid`
   
   Please update your imports to use :mod:`frb.halos.models` instead.

.. note::
   This module requires the optional ``hmf_emulator`` package.
   Install it with: ``pip install hmf_emulator``
