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

* `python <http://www.python.org/>`_ versions 2.7, or 3.6 or later
* `numpy <http://www.numpy.org/>`_ version 1.11 or later
* `astropy <http://www.astropy.org/>`_ version 2.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.19 or later

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python|numpy|astropy|scipy"

If the packages have been installed, this command should print
out all the packages and their version numbers.

Installing frb
==============

Presently, you must download the code from github::

	#go to the directory where you would like to install specdb.
	git clone https://github.com/FRB/FRB.git

From there, you can build and install with

	cd FRB
	python setup.py install  # or use develop


This should install the package and scripts.
Make sure that your PATH includes the standard
location for Python scripts (e.g. ~/anaconda/bin)


.. _download-public:


