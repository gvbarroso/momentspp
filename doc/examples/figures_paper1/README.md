# Predicting diversity along chromosomes with moments++ and bgshr
Scripts for generating figures of the manuscript introducing moments++ (Barroso & Ragsdale 2025).

These are presented here as examples. For reference, the directory 'fwdpy11_simulations' contains scripts that we used to simulate data with fwdpy11 software (https://github.com/molpopgen/fwdpy11)

Executing the notebooks in the 'figures' directory requires the following data:

Pre-computed lookup tables (in the 'lookup_tables' directory).

The OMNI recombination map (https://www.sciencedirect.com/science/article/pii/S0002929720301634?via%3Dihub) and a bed file with positions of exonic sites (human chr 2), both of which are found in the 'data' directory.

The results from forward-in-time simulations performed with fwdpy11 (presented in the 'simulation_outputs' directory).

Using moments++ output (stored in the lookup tables) to predict B-value maps requires the bgshr python package (https://github.com/apragsdale/bgshr/tree/main/bgshr).

After installation of bgshr (and potentially activation of the associanted conda environment), reproducing the figures in the original momentspp publication amounts to navigating to the 'figures' directory, opening each \*.ipynb file individually (e.g.., `jupyter-lab cBGS_dynamics.ipynb`) and executing all cells in order. 
