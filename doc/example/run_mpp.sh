#!/bin/bash

# Date created: March 13 2025
# Created by: Gustavo Valadares Barroso

# This script runs moments++ with parameters specified in test_model.yaml
# test_model.yaml is a demes file augmented by listing the rates of recombination, mutation and the selection coefficient of a pure two-locus model.
# This short example shows that the options file (opt.bpp) can be quite flexible:

momentspp params=opt.bpp F=test_model.yaml O=25 T=1000
