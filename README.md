# S_cerevisiae_DNA_rep_model

We created a stochastic model for S. cerevisiae whole-genome replication in which origins compete to associate with limtted firing factors, needed for activation, which then recycle for reuse.

This repository contains the Beacon Calculus (bcs) script used for the model. Beacon Calculus is an open access software available at https://github.com/MBoemo/bcs.git (Boemo et al., 2020). 

This repository also contains example Python sripts used to extract various measures of DNA replication dynamics from the simulation output.

## Descriptions of files:

### Beacon calculus script

200FF0p05d_fitted.bc : Model with 200 firing factors and recycling rates of 0.05. This is the main model variation used for most analysis. The origin firing rates in this script have already been fitted to experimental replication timing data (Muller et al., 2014).

### Python scripts

fitting_060624.py : Fits the model iteratively to experimental replication timing data.

running_080624.py : Runs a specified number of simulations of the model.

repTime_080624.py : Calculates the simulated mean, std, median, upper quartile, and lower quartile of replication timing at each kb of the genome.

IOD_Sphase_060624.py : Calculates the simulated inter-origin distances (IOD) and durations of the simulations (which represents the simulated length of S-phase).

efficiencies_060624.py : Calculates the simulated efficiency of each origin.

RFD_080624.py : Calculates the simulated replication fork directionality (RFD) at each kb of the genome.

activeForks_080624.py : Calculates the simulated mean number of active replication forks over time.

firingTimes_080624.py : Records the times at which each origin fires in the different simulations (from which their firing time distributions can later be plotted).

replicons_080624.py : Calculates the simulated mean replicon length of each origin.

freeFF_080624.py : Calculates the mean number of availabe firing factors at different times over the cause of the simulations.


The experimentally determined replication times were sourced from (Muller et al., 2014). Linear interpolation was used to estimate the experimental replication timings at each kb.

## References:

Boemo, M.A., Cardelli, L. and Nieduszynski, C.A., 2020. The Beacon Calculus: A formal method for the flexible and concise modelling of biological systems. PLoS computational biology, 16(3), p.e1007651.

Müller, C.A., Hawkins, M., Retkute, R., Malla, S., Wilson, R., Blythe, M.J., Nakato, R., Komata, M., Shirahige, K., de Moura, A.P. and Nieduszynski, C.A., 2014. The dynamics of genome replication using deep sequencing. Nucleic acids research, 42(1), pp.e3-e3.
