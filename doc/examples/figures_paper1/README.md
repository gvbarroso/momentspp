# Predicting diversity along chromosomes with moments++ and bgshr
Scripts for generating figures of the manuscript introducing moments++ (Barroso & Ragsdale 2025).

These are presented here as examples. The directory 'fwdpy11_simulations' contains scripts to simulate data with fwdpy11 software (https://github.com/molpopgen/fwdpy11)

Executing the notebooks in the 'figures' directory requires the following data:
'lookup_tables' contains pre-computed lookup tables
'data' contains the OMNI recombination map (https://www.sciencedirect.com/science/article/pii/S0002929720301634?via%3Dihub) and a bed file with positions of exonic sites (human chr 2).
'simulation_outputs' were generated with fwdpy11 and are presented as a convenience.

Using moments++ output (stored in the lookup tables) to predict B-value maps requires the bgshr python package (https://github.com/apragsdale/bgshr/tree/main/bgshr).
