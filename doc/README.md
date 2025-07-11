# A starter guide for creating the input files for moments++

`moments++` basic functionally is built around the bio++ libraries (https://github.com/BioPP)
As such, input parameters are specified in an 'options file'. Let's call this options file opt.bpp. To run moments++ from the command-line:

```
momentspp param=opt.bpp
```

The options file for `moments++ v0.1` is quite simple. Although under the hood there are a number of options that enable flexibility for handling more complicated situations, these have default parameters that we don't have to worry about here. The options that we want to specify in opt.bpp are the following:

## number of computing threads

The number of threads that `moments++` can use (INTEGER). This allows the Eigen3 libraries to paralellize the linear algebra computations. By default, `moments++` will use all threads available in the machine.

```
number_threads = 4
```

## time steps to sample temporal heterozigosity values

In a non-equilibrum model, this specifies the periodicity (in units of generations) with which `moments++` outputs values of $p(1-p)$ and $q(1-q)$ (after the first population size change, forwards-in-time). The default value is zero, meaning these temporal dynamics of diversity are NOT recorded and only the statistics at the most recent time are written to \*\_expectations.txt. Specify this with e.g. (INTEGER):

```
time_steps = 1000
```

## order of 1-2p factors

The count of $1-2p$ factors to be included in the basis **v**. As discussed in the original publication of `moments++` (Barroso & Ragsdale, 2025), this is mainly a function of the strength of selection, but it also depends on the demography ($N[t]$) and the mutation rate $\mu$. It is therefore tricky to suggest a rule of thumb, but a good starting value is something between 2 and 4 times $N_e \times s$. (The example given in 'two_locus_time.Rmd' suggests how to adjust the order as such). Specify this option in opt.bpp with e.g. (INTEGER):

```
factor_order = 50
```

## the yaml (demes) file with the prescribed evolutionary model

This is the file that contains all the parameters describing the model for which `moments++` will output expected two-locus statistics. Think about it as an augmented demes file (https://popsim-consortium.github.io/demes-spec-docs/main/tutorial.html) where the mutation rate $\mu$, recombination rate $r$ and selection coefficient $s$ are specified as metadata. What makes this an *augmented* demes file (as opposed to a regular demes file) is that the demes standard prohibits the specification of metadata as dictionaries. It turns out that it is *really* convenient for `moments++` to have $\mu$, $r$ and $s$ be specified as dictionaries, so we go ahead an do it, at the cost of this explanation. This makes it straightforwad to write models with e.g. epoch-specific selection coefficients. Noting that `moments++` is currently restricted to single-population models, specify this option with e.g. (STRING):

```
demes_file = model.yaml
```

The 'test_run' directory contains an example file ('test_model.yaml'), whereas 'two_locus_time.Rmd' shows how to automatize the writing of augmented demes files.

Finally, the options file supports passing bash-style variables as input. For example, if opt.bpp looks like the following

```
number_threads = 8
time_steps = $(T)
factor_order = $(O)
demes_file = $(F)
```

we can execute, e.g.:

```
momentspp param=opt.bpp T=2000 O=75 F=model_pop1.yaml
```

Once again, an example of exploiting this functionality for automatizing the assembly of a lookup table is given in 'two_locus_time.Rmd'
