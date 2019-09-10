.. highlight:: rest

**************
Installing frb
**************

This document describes how to install the `frb`
repository.  We also describe
:ref:`download-public`.

Installing Dependencies
=======================
We have and will continue to keep the number of dependencies low.
There are a few standard packages that must be installed
and one package `linetools` under review for
`astropy` affiliated status.

In general, we recommend that you use Anaconda for the majority of
these installations.

Detailed installation instructions are presented below:

Python Dependencies
-------------------

frb depends on the following list of Python packages.

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_
to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 3.6 or later
* `numpy <http://www.numpy.org/>`_ version 1.14 or later
* `astropy <http://www.astropy.org/>`_ version 3.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.19 or later
* `healpy <https://healpy.readthedocs.io/en/latest/index.html>`_ version 1.12 or later

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python|numpy|astropy|scipy"

If the packages have been installed, this command should print
out all the packages and their version numbers.

The following packages are required to access surveys (e.g. SDSS, DES)
for data that may be associated to an FRB:

* `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ v0.3.8 (v0.3.9 currently fails)
* `datalab-client <https://github.com/noaodatalab/datalab/>`_ v1.1 or later
* `pyvo <https://pyvo.readthedocs.io/en/latest/>`_  version 0.9.2 or later
* `PIL <https://pillow.readthedocs.io/en/5.3.x/>`_  version 5.3 or later (only for SDSS cutouts)
* `requests <https://pillow.readthedocs.io/en/5.3.x/>`_  version 2.18 or later

The following are required to run spectral line analysis (e.g. frb.galaxies.nebular):

* `linetools <https://github.com/linetools/linetools>`_  version 0.3 or later

The following are required to run some of the halo codes:

* `hmf_emulator <https://github.com/profxj/hmf_emulator.git>`_  WARNING: This is JXP's fork.
* george :: Use pip
* `class <https://github.com/lesgourg/class_public>`_ version 2.7 or greater

The following are required to build host galaxy objects:

* `pPXF <https://pypi.org/project/ppxf/>`_ version 6.7 or greater

For pPXF, you will also likely need to modify the standard install
to use the Chabrier libraries.  See the InstallNotes in this
`Google Drive <https://drive.google.com/drive/folders/1_nu8IiBm0-dnkpoKBcoXyQuqbsrYHNXh?usp=sharing>`_.


Installing frb
==============

Presently, you must download the code from github::

	#go to the directory where you would like to install specdb.
	git clone https://github.com/FRBs/FRB.git

From there, you can build and install with::

	cd FRB
	python setup.py install  # or use develop


This should install the package and scripts.
Make sure that your PATH includes the standard
location for Python scripts (e.g. ~/anaconda/bin)


.. _download-public:


