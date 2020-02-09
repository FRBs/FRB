********
Database
********

These docs describe several utility methods to access the
:doc:`data` held in this Repository.  This includes both FRBs
and their host galaxies.

The data and measurements for each are stored in the
Repository in individual JSON files.  These are convenient
for machines but not so useful for humans.  You can
see examples of accessing them in Python in the
`FRB <https://github.com/FRBs/FRB/blob/master/docs/nb/FRB_Event.ipynb>`_
and
`Hosts <https://github.com/FRBs/FRB/blob/master/docs/nb/FRB_Galaxies.ipynb>`_
Notebooks.

There is also a Database Notebook showing examples
described here.

It is intended to pivot the tables on the "FRB" tag.

====
FRBs
====

Use the `build_table_of_frbs()` method to generate a
`pandas <https://pandas.pydata.org/>`_ table of
key quantities::

    from frb import frb
    frb_tbl, tbl_units = frb.build_table_of_frbs()

The `frb_tbl` is a Pandas table containing items like the name,
RA, DEC, DM, etc.  `tbl_units` contains the units of
each column.

=====
Hosts
=====

Use the `build_table_of_hosts()` method to generate a
`pandas <https://pandas.pydata.org/>`_ table of
key quantities::

    from frb.galaxies import frbgalaxy
    host_tbl, tbl_units = frbgalaxy.build_table_of_hosts()

This includes items like the photometry, nebular emission
line fluxes, and derived quantities (e.g. stellar mass).


====
Misc
====

In pandas, here is how you would merge the FRBs and
Hosts tables::

    import pandas as pd
    joined_tbl = pd.merge(frb_tbl, host_tbl, on='FRB', how='outer')

That will give you a large table with NaN's for masked
values.
