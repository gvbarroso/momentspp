# Predicting diversity along chromosomes with moments++ and bgshr
Scripts for generating figures of the manuscript introducing moments++ (Barroso & Ragsdale 2025).

These are presented here as examples. For reference, the directory 'fwdpy11_simulations' contains scripts that we used to simulate data with fwdpy11 software (https://github.com/molpopgen/fwdpy11).

Executing the notebooks in the 'figures' directory requires the following:

* Pre-computed lookup tables (in the 'lookup_tables' directory).
* The Omini recombination map (https://www.sciencedirect.com/science/article/pii/S0002929720301634?via%3Dihub) and a bed file with positions of exonic sites (human chr 2), both of which are found in the 'data' directory.
* The results from forward-in-time simulations with fwdpy11 (presented in the 'simulation_outputs' directory).

Using moments++ output (stored in the lookup tables) to predict B-value maps requires the `bgshr` python package (https://github.com/apragsdale/bgshr/tree/main/bgshr).

After installation of `bgshr` (and activation of the associated environment), reproducing the figures in the original publication amounts to navigating to the 'figures' directory, opening each \*.ipynb or *.R file individually (e.g.., `jupyter-lab cBGS_dynamics.ipynb`) and executing all code chunks in order. 
Besides `bgshr`, other python or R dependencies to execute these scripts are listed at the top of each file.

*B-map\_diff\_time.ipynb* generates Figures S4 and S5. \
*cBGS\_dynamics.ipynb* generates Figure S3. \
*dfe-fitting\_2\_epochs.ipynb* generates Figures 4 and S11. \
*human-chr2-comparison.ipynb* generates Figures 3 and S6. \
*dfe-fitting\_4\_epochs.ipynb* generates Figures 5 and S12. \
*steady-state-validation.ipynb* generates Figures 2. 

The following R scripts generate the remaining supplemental figures: \
*B\_iters.R* generates Figure S1. \
*two\_locus\_plots.R* generates Figures S8, S9 and S10. 

*bench\_1-pop.Rmd* generates Figure S7. \
**NOTE**: execution of *bench\_1-pop.Rmd* requires `moments.TwoLocus` -- the python sofware we benchmark against. This is part of the `moments package`: https://github.com/apragsdale/moments

*two\_locus\_sims\_plots.R* generates Figures A1 and A2. \
**NOTE**: execution of two\_locus\_sims\_plots.R is based on a massive set of simulations performed with `twoLocusSim`, a companion software to `moments++` and part of this repository. Hosting these simulated data online (46 Gb compressed) seems impractical, unfortunately.
