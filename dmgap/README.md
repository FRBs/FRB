# Dispersion Measure Gap

This code uses pulsar and FRB data to place constraints on the disperison measure (DM) of the Milky Way galactic halo. Density estimation using field theory ([DEFT](https://arxiv.org/pdf/1804.01932.pdf)) is invoked with the [Python](https://www.python.org/) package [SUFTWARE](https://suftware.readthedocs.io/en/latest/).


## Density Estimation Using Field Theory
To find constraints on the DM gap, follow the steps below:

1. Clean input FRB and pulsar data from transient_data/. This may take a little while.
    
        $ python data_sort.py
       
2. Calculate and plot DEFT PDFs

        $ python deft_obs_calcs_and_plots.py

3. Find DM values for different confidence intervals

        $ python dm_obs_ci.py


## FRB Simulations Using DEFT
To evaluate the expected performance of DEFT as more data becomes available, we simulate a PDF of FRB DMs. By taking draws of  different sizes and applying DEFT to approximate the PDF, one can see how the DEFT approximatiom of the simulated PDF improves as sample size increases.

1. Run simulation to generate FRB PDF

        $ python frb_sim.py

2. Calculate DEFT approximations using different sample sizes 

        $ python deft_sim_calcs.py

3. Plot the PDFs

        $ python deft_sim_plots.py


## Comparing DEFT to KDE with FRB Data
We note that kernel density estimation (KDE)&mdash;a traditional density estimation technique&mdash;may present itself as a more obvious choice. Here we provide simulations to motivate the use of DEFT in this context.
