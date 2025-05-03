# S_cerevisiae_DNA_rep_model

We created a stochastic model for Saccharomyces cerevisiae whole-genome replication in which origins compete to associate with limited firing factors, needed for activation, which then recycle for reuse.

This repository contains the Beacon Calculus (bcs) script used for the model. Beacon Calculus is an open access software available at https://github.com/MBoemo/bcs.git (Boemo et al., 2020). 

(bcs must be cloned into the same folder as the Python scripts in this repository, therefore users are advised to clone this repository first, then clone and make bcs inside it).

This repository also contains example Python scripts used to extract various measures of DNA replication dynamics from the simulation output.

## Descriptions of files:

### Experimental data files (data subfolder)

trep.wig : Replication timing data from (Müller et al., 2014).

oriDB_S_cerevisiae.txt : origins from the OriDB (Siow et al. 2012).

completeExpRepTime.csv : Processed replication timing data used in subsequent analysis (interpolated to estimate the timings for every kb). 

origin_positions.csv : Processed origin data used in subsequent analysis (only including "confirmed" and "likely origins").



### Beacon Calculus script for full model (bcs_scripts subfolder)

200FF0p05d_fitted.bc : Model with 200 firing factors and recycling rates of 0.05. This is the main model variation used for most analysis. The origin firing rates in this script have already been fitted to experimental replication timing data (Müller et al., 2014).

200FF0p05d_mapOri_fitted.bc : Same as above, except origin and replication fork processes have an additional "ori" parameter which allows the origin forks originate from to be tracked for computing replicon lengths and the number of active replication forks. 

### Example model outputs (bcs_output subfolder)

200FF0p05d_fitted_s10.simulation.bcs : Output from 10 simulations of the fitted model

200FF0p05d_mapOris_fitted_s10.simulation.bcs : Output from 10 simulations of a fitted version of the model with the additional "ori" parameter.

### Python scripts for full model

processing_expRepTimes.py : Processes the experimentally determined replication timings from trep.wig (sourced from Müller et al., 2014). Uses linear interpolation to estimate the replication timings at each kb.

prossessing_originPositions.py : Extracts the positions of 'Confirmed' and 'Likely' origins from the OriDB, oriDB_S_cerevisiae.txt (Siow et al. 2012).

writing_bcs_scripts.py : Writes Beacon Calculus scripts with custom parameters.

writing_bcs_scripts_mapOris.py : Same as above but for writing Beacon Calculus scripts with the additional "ori" parameter (for computing replicon lengths and active replication forks).

fitting.py : Fits the model iteratively to experimental replication timing data.

fitting_mapOris.py : Same as above but for fitting Beacon Calculus scripts with the addition "ori" parameter.

running.py : Runs a specified number of simulations of the model.

repTime.py : Calculates the simulated mean, std, median, upper quartile, and lower quartile of replication timing at each kb of the genome.

IOD_Sphase.py : Calculates the simulated inter-origin distances (IOD) and durations of the simulations (which represents the simulated length of S-phase).

efficiencies.py : Calculates the simulated efficiency of each origin.

RFD.py : Calculates the simulated replication fork directionality (RFD) at each kb of the genome.

activeForks.py : Calculates the simulated mean number of active replication forks over time.

firingTimes.py : Records the times at which each origin fires in the different simulations (from which their firing time distributions can later be plotted).

replicons.py : Calculates the simulated mean replicon length of each origin.

freeFF.py : Calculates the mean number of available firing factors at different times over the cause of the simulations.

## Instructions for running code

1: Download the experimental replication timings from (Müller et al., 2014) and the positions of origins from the OriDB (Siow et al. 2012), (or use trep.wig and oriDB_S_cerevisiae.txt provided).

2: Process the data into a compatible format by running: processing_expRepTimes.py and processing_originPositions.py (Optional if using pre-processed data: completeExpRepTime.csv and origin_positions.csv).

3: Generate a Beacon Calculus script with custom parameters by running: writing_bcs_scripts.py and writing_bcs_scripts_mapOris.py (optional).

4: Fit Beacon Calculus scripts to experimental replication timings by running: fitting.py and fitting_mapOris.py (not necessary if using the pre-fitted scripts 200FF0p05d_fitted.bc and 200FF0p05d_mapOris_fitted.bc).

5: Run Beacon Calculus simulations by running: running.py. (500 simulations are recommended).   

6: Analyse the model output using the remaining Python scripts:

(repTime.py, IOD_Sphase.py, efficiencies.py, RFD.py, activeForks.py, firingTimes.py, replicons.py,freeFF.py)

Results from analyses are saved to the "output" subfolder

## Chromosome II example (chr2_example subfolder)

Code for a smaller version of the model, just for chromosome II, is provided as an example designed to be run locally. The following scripts are versions of the above, but modified for chromosome II. 

### Beacon calculus script for chromosome II model

15FF0p05d_unfitted_chrII.bc : Unfitted version of the model for chromosome II only. 15 firing factors are used based on the proportion of origins on chromosome II (needs fitting by running fitting_chrII.py).

15FF0p05d_unfitted_mapOri_chrII.bc : Unfitted version of the chromosome II model with the additional "ori" parameter for computing replicon lengths and numbers of active replication forks (needs fitting by running fitting_mapOris_chrII.py).

### Python scripts for chromosome II model

fitting_chrII.py 

fitting_mapOris_chrII.py 

repTime_chrII.py 

IOD_Sphase_chrII.py 

efficiencies_chrII.py 

RFD_chrII.py 

activeForks_chrII.py 

firingTimes_chrII.py 

replicons_chrII.py 

freeFF_chrII.py 

(Running these scripts relies on the processed datasets). 

The chromosome II Beacon calculus scripts can by run using running.py (same as for the full model).

## References:

Boemo, M.A., Cardelli, L. and Nieduszynski, C.A., 2020. The Beacon Calculus: A formal method for the flexible and concise modelling of biological systems. PLoS computational biology, 16(3), p.e1007651.

Müller, C.A., Hawkins, M., Retkute, R., Malla, S., Wilson, R., Blythe, M.J., Nakato, R., Komata, M., Shirahige, K., de Moura, A.P. and Nieduszynski, C.A., 2014. The dynamics of genome replication using deep sequencing. Nucleic acids research, 42(1), pp.e3-e3.

Siow, C.C., Nieduszynska, S.R., Müller, C.A. and Nieduszynski, C.A., 2012. OriDB, the DNA replication origin database updated and extended. Nucleic acids research, 40(D1), pp.D682-D686.

