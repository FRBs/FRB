**************
Installing frb
**************

This document describes how to install the `frb`
repository.

Quick Start
===========

The recommended installation path is now through ``pip`` using the
repository's ``pyproject.toml`` metadata.

We recommend Python 3.11 or later in a fresh virtual environment or
conda environment.

For example, with ``venv``::

	git clone https://github.com/FRBs/FRB.git
	cd FRB
	python -m venv .venv
	source .venv/bin/activate
	pip install --upgrade pip
	pip install -e .

or with conda::

	conda create -n frb python=3.11 -y
	conda activate frb
	git clone https://github.com/FRBs/FRB.git
	cd FRB
	pip install --upgrade pip
	pip install -e .

This installs the base package plus its core dependencies, including
the git-based requirements such as ``linetools``, ``ne2001``, and
``astropath``.

.. note::

	``setup.py`` is kept for backward compatibility, but the supported
	installation workflow is ``pip install -e .``.
	Likewise, ``frb/requirements.txt`` and ``frb/optional_requirements.txt``
	are now reference files rather than the primary installation source.

Base Dependencies
-----------------

The base install is defined in ``pyproject.toml`` and includes the main
runtime requirements for the package.  At the time of writing, these include
packages such as:

* `python <http://www.python.org/>`_ 3.11 or later
* `numpy <http://www.numpy.org/>`_ 2.2 or later
* `scipy <http://www.scipy.org/>`_ 1.17 or later
* `astropy <http://www.astropy.org/>`_ 7.1 or later
* `pandas <https://pandas.pydata.org/>`_ 2.2 or later
* `matplotlib <https://matplotlib.org/>`_ 3.7 or later
* `healpy <https://healpy.readthedocs.io/en/latest/index.html>`_ 1.19 or later
* `requests <https://requests.readthedocs.io/>`_ 2.18 or later
* `dust_extinction <https://dust-extinction.readthedocs.io/en/latest/>`_
* `photutils <https://photutils.readthedocs.io/en/stable/>`_
* `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ 0.4.11 or later
* `astro-datalab <https://github.com/noaodatalab/datalab/>`_
* `pyvo <https://pyvo.readthedocs.io/en/latest/>`_ 1.5.3 or later
* `ligo.skymap <https://lscsoft.docs.ligo.org/ligo.skymap/>`_ 2.3.0 or later
* `numba <https://numba.readthedocs.io/>`_ 0.50 or later
* `tqdm <https://tqdm.github.io/>`_
* `linetools <https://github.com/linetools/linetools>`_
* `ne2001 <https://github.com/FRBs/ne2001.git>`_
* `astropath <https://github.com/FRBs/astropath>`_

Extras
------

The current ``optional`` packages covers:

* `ppxf <https://github.com/SunilSimha/frb_ppxf>`_
* `pymc3 <https://pypi.org/project/pymc3/>`_
* `pcigale <https://cigale.lam.fr/>`_
* `pathos <https://pypi.org/project/pathos/>`_
* `hmf_emulator <https://github.com/AemulusProject/hmf_emulator>`_
* `specdb <https://github.com/specdb/specdb>`_
* `pyregion <https://pypi.org/project/pyregion/>`_
* `spectral-cube <https://pypi.org/project/spectral-cube/>`_

We do not include these as dependencies as they may require additional system setup.
These are however necessary to run specific modules of the code-base,
and we recommend installing them if you intend to use those modules.

The following feature-specific dependencies may still require separate,
manual installation:

* `FRB-pulsars <https://github.com/FRBs/pulsars>`_ for ``frb.surveys.psrcat``
* `asymmetric_kde <https://github.com/tillahoffmann/asymmetric_kde>`_ for ``frb.dm_kde``
* `scikit-image <https://scikit-image.org/>`_ for some image-analysis workflows

The ``hmf_emulator`` package is especially important for halo-related modules,
but it may require additional local system setup depending on your platform.

Environment Variables
---------------------

Most users do not need any extra environment variables for a basic install.
Some workflows do:

* ``FRB_GDB`` is needed by several build scripts and tests that expect access
	to the FRB database or generated data products.
* ``EAZYDIR`` is needed for the EAZY wrappers described below.

Our CIGALE wrappers use custom filter files not
provided by their current release (e.g DES, Pan-STARRS).
See the instructions for adding those as needed.

Installing frb
==============

If you already cloned the repository, install it from the checkout root with::

	pip install -e .

This installs the package in editable mode and exposes the scripts in your
active Python environment.

Adding additional CIGALE filter files
=====================================

We have had to add new filters to CIGALE (e.g. from
DES and Pan-STARRS).
You can find these filter files in
`frb.data.analysis.CIGALE`.
If you use any of those surveys,
**our wrappers will not work without them.**

Here are the steps to update CIGALE:

* cd frb/data/analysis/CIGALE
* Add any desired filter into your CIGALE code base with:  `pcigale-filters add file.dat`

Note that DECaLs uses the BASS and MzLS data files.

EAZY setup
==========

In order to perform photo-z estimation
with EAZY using our wrappers, the following
changes need to be made.

* Add an environment variable `EAZYDIR` that points to your EAZY installation. Add this to your `bashrc`::

	export EAZYDIR="/path/to/eazy-photoz/"

* Locate the `templates` folder in `$EAZYDIR` and edit the paths present in `*.spectra.param`. Replace all SED file paths with the absolute paths. For instance, in `$EAZYDIR/templates/eazy_v1.3.spectra.param`, replace::

	templates/EAZY_v1.1_lines/eazy_v1.1_sed1.dat

with::

	/path/to/eazy-photoz/templates/EAZY_v1.1_lines/eazy_v1.1_sed1.dat


.. _download-public:


pPXF
====

Our pPXF wrapper currently uses an older version of the code 
(v 6.7.17) and a few custom files. These are made available in a separate repository,
and you can install it with::

	pip install git+https://github.com/SunilSimha/frb_ppxf.git

If instead you want to install this manually and apply the relevevant patches yourself,
see  the InstallNotes in this
`Google Drive <https://drive.google.com/drive/folders/1_nu8IiBm0-dnkpoKBcoXyQuqbsrYHNXh?usp=sharing>`_.
Contact JXP if you have any questions about this process.