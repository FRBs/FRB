### For KDE analysis on FRBs and pulsars used in Platts, Prochaska and Law (2020).

Pulsar data is available in the [FRBs/pulsars/pulsars/data/atnf_cat](https://github.com/FRBs/pulsars/tree/master/pulsars/data/atnf_cat) repo and FRB data is avialable in the [FRBs/FRB/frb/data/FRBs](https://github.com/FRBs/FRB/tree/master/frb/data/FRBs) repo.

Requires the [`asymmetric_kde`](https://github.com/tillahoffmann/asymmetric_kde) package.

First sort data for processing (removes sources near Milky Way disk and the Magnellanic Clounds) by running `sort_transient_data.py`

`transient_kdes.py` generates KDEs for pulsars and FRBs from observations and finds their maxima/minima.

`dm_frb_sim.py` simulates a PDF of FRB DMs, applies KDE and finds the minima.
