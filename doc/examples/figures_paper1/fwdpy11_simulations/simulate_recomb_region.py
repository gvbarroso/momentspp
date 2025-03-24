import numpy as np
import fwdpy11
import time
from dataclasses import dataclass
import sys
from typing import List
from collections import defaultdict
import pickle
import moments


@dataclass
class Recorder:
    burnin: int
    spacing: int

    def __call__(self, pop, sampler):
        if pop.generation % 1000 == 0:
            print(f"{time.ctime()}, at generation {pop.generation}")
        if pop.generation >= self.burnin and pop.generation % self.spacing == 0:
            sampler.assign(range(pop.N))


def runsim(seed, Ne=1000, s=-0.01, u=1e-8, r=1e-8, L=100000, num_reps=1):
    """
    Single region, constant recombination, mutation rates,

    Ne: The effective population size
    R: The total recombination rate to the right of the non-recombining locus
    U: The total mutation rate (u*L) within the recombining locus
    s: The selection coefficients of the selected alleles

    Returns a tree sequence
    """
    # Set the rng with the given seed
    rng = fwdpy11.GSLrng(seed)
    # The parameters for the simulation

    U = u * L
    R = r * L

    burnin = 20 * Ne
    spacing = Ne // 10
    simlen = burnin + (num_reps - 1) * spacing
    print(f"{time.ctime()}, total simulation length:", simlen)

    pdict = {
        # Multiplicative selection model
        "gvalue": fwdpy11.Multiplicative(2.0),
        # Rates: (neutral, selection, recombination)
        "rates": (0.0, U, None),
        "nregions": [],
        # Selection within the non-recombining locus
        "sregions": [fwdpy11.ConstantS(0, L, 1, s, 1.0)],
        # Recombination to the right of this locus
        "recregions": [fwdpy11.PoissonInterval(0, L, R)],
        # Evolve a single deme of size N for 20*N generations
        "demography": fwdpy11.ForwardDemesGraph.tubes([Ne], simlen),
        "simlen": simlen,
    }
    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation(Ne, L)
    recorder = Recorder(burnin=burnin, spacing=spacing)
    print(f"{time.ctime()}, starting simulation")
    fwdpy11.evolvets(
        rng,
        pop,
        params,
        100,
        recorder=recorder,
        suppress_table_indexing=True,
    )
    ts = pop.dump_tables_to_tskit()
    return ts


def get_diversity(ts, L=100000, winsize=1):
    """
    Linked neutral sites are between 1 and 2.
    """
    assert ts.sequence_length == L

    print(f"{time.ctime()}, getting diversity at each timepoint")
    # get sample times
    sample_sets = defaultdict(list)
    for samp in ts.samples():
        t = ts.node(samp).time
        sample_sets[t].append(samp)

    # set up windows
    windows = [0, L / 2 - winsize / 2, L / 2 + winsize / 2, L]

    print(f"{time.ctime()}, created windows")
    # get divergences after simplifying to samples at each sample time
    divs = np.zeros(len(sample_sets.keys()))
    sfs = np.zeros((len(sample_sets.keys()), len(sample_sets[0]) + 1))
    for i, t in enumerate(sorted(sample_sets.keys())):
        ts_simp = ts.simplify(samples=sample_sets[t])
        branch_div = ts_simp.diversity(windows=windows, mode="branch")
        branch_sfs = ts_simp.allele_frequency_spectrum(
            windows=windows, mode="branch", polarised=True
        )
        divs[i] = branch_div[1]
        sfs[i] = branch_sfs[1]
    print(f"{time.ctime()}, computed diversity")
    return divs, sfs


if __name__ == "__main__":
    Ne = 1000
    r = 1e-8
    L = 100000

    n = 20

    seed = int(sys.argv[1]) + 1
    s = float(sys.argv[2])
    u = float(sys.argv[3])

    num_reps = 1000

    ts = runsim(seed, Ne=Ne, s=s, u=u, r=r, L=L, num_reps=num_reps)
    divs, sfs = get_diversity(ts, L=L)

    ave_div = np.mean(divs)
    ave_sfs = np.mean(sfs, axis=0)

    sfs_proj = moments.Spectrum(ave_sfs).project([n])

    assert np.isclose(sfs_proj.pi(), ave_div)

    data = {"div": ave_div, "sfs": sfs_proj.data}
    fname = f"div_sfs.Ne_{Ne}.r_{r}.L_{L}.u_{u}.s_{s}.seed_{seed}.pkl"

    with open(fname, "wb+") as fout:
        pickle.dump(data, fout)
    print(f"{time.ctime()}, saved data. Done!")
