# mpp_figures
Scripts for generating figures of the manuscript introducing moments++ (Barroso & Ragsdale 2025).

These are presented here as examples.

Running the scripts in the 'figures' directory requires the following data:
'lookup_tables' contains pre-computed lookup tables
'data' contains the OMINI recombination map (Yoruba population) and a bed file with positions of exonic sites (chr 2).
'simulation_outputs' were generated with fwdpy11.

For using the two-locus moments++ predictions (contained in the lookup tables) to compute B-value maps, we will need the bgshr python package (https://github.com/apragsdale/bgshr/tree/main/bgshr).

Installation:

