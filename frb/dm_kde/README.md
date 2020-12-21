### For KDE analysis on FRBs and pulsars used in Platts, Prochaska and Law (2020).

Pulsar data is available in the FRBs/pulsars/pulsars/data/atnf_cat repo and FRB data is avialable in the FRBs/FRB/frb/data/FRBs repo. 

First sort data for processing (removes sources near Milky Way disk and the Magnellanic Clounds) by running `sort_transient_data.py`

`transient_kdes.py` generates KDEs for pulsars and FRBs from observations and finds their maxima/minima.

`dm_frb_sim.py` simulates a PDF of FRB DMs, applies KDE and finds the minima.
