"""
Run simulation with recombination map, exon regions, DFE, and (possibly) size change.
"""

import numpy as np
import fwdpy11
import time
from dataclasses import dataclass
import sys
from collections import defaultdict
import demes
import os
import gzip


@dataclass
class Recorder:
    burnin: int
    spacing: int

    def __call__(self, pop, sampler):
        if pop.generation % 1000 == 0:
            print(f"{time.ctime()}, at generation {pop.generation}")
        if pop.generation >= self.burnin and pop.generation % self.spacing == 0:
            sampler.assign(range(pop.N))


def build_demographic_model(Nsim, sample_gens):
    b = demes.Builder(time_units="generations")
    if sample_gens == 0:
        epochs = [dict(start_size=Nsim)]
    else:
        epochs = [
            dict(start_size=Nsim, end_time=sample_gens),
            dict(start_size=Nsim),
        ]
    b.add_deme(
        name="deme0",
        epochs=epochs,
    )
    return b.resolve()


def build_rec_regions(rec_map, L, Q=1):
    rec_pos = [0]
    rates = [0]
    col = None
    with gzip.open(rec_map, "rb") as fin:
        for line in fin:
            l = line.decode().split()
            if l[0] == "Position(bp)":
                col = l.index("Rate(cM/Mb)")
            else:
                assert col is not None
                rate = float(l[col]) / 1e6 / 100
                if rate == rates[-1]:
                    continue
                else:
                    rec_pos.append(int(l[0]))
                    rates.append(float(l[col]) / 1e6 / 100)
    if L is None:
        raise ValueError("Need to provide chromosome length")
    rec_pos.append(L)
    assert len(rec_pos) == len(rates) + 1
    rec_regions = []
    for left, right, rate in zip(rec_pos[:-1], rec_pos[1:], rates):
        rec_regions.append(
            #  Scale by Q
            fwdpy11.PoissonInterval(beg=left, end=right, mean=(right - left) * rate * Q)
        )

    return rec_regions


def build_sel_regions(exons, mean, shape, Q=1):
    exon_left = []
    exon_right = []
    with open(exons, "r") as fin:
        for line in fin:
            c, l, r = line.split()
            l = int(l)
            r = int(r)
            if len(exon_left) == 0:
                exon_left.append(l)
                exon_right.append(r)
            elif l <= exon_right[-1]:
                assert l >= exon_left[-1]
                exon_right[-1] = r
            else:
                exon_left.append(l)
                exon_right.append(r)
    Lexons = np.sum(np.array(exon_right) - np.array(exon_left))
    sel_regions = []
    for l, r in zip(exon_left, exon_right):
        sel_regions.append(
            fwdpy11.GammaS(
                beg=l,
                end=r,
                weight=1,
                mean=-mean * Q,  # scale by Q
                shape_parameter=shape,
                h=1,
            )
        )
    return Lexons, sel_regions


def runsim(seed, Ne=10000, chrom=2, u=1.5e-8, L=None, Q=1, n_samples=1):
    """
    If Q is not 1, Ne is divided by Q, u and r are increased by Q, as are s values
    Returns a tree sequence
    """
    # Set the rng with the given seed
    rng = fwdpy11.GSLrng(seed)
    # The parameters for the simulation

    N = int(Ne / Q)

    # DFE mean, shape from Kim et al 2017 (LuCamp)
    shape = 0.215
    scale = 562.1
    mean = scale * shape / (2 * Ne)

    exon_path = "data_files"
    exon_file = f"chr{chrom}_1-315215_exons.bed"
    exons = os.path.join(exon_path, exon_file)
    exon_regions = []
    Lexons, sel_regions = build_sel_regions(exons, mean, shape, Q=Q)

    rec_map_path = "data_files"
    rec_map_file = f"YRI-{chrom}-1-315215-final.txt.gz"
    rec_map = os.path.join(rec_map_path, rec_map_file)
    rec_regions = build_rec_regions(rec_map, L, Q=Q)

    # Rates
    nonsyn_fac = 2.31 / (2.31 + 1)
    U = u * nonsyn_fac * Lexons * Q

    burnin = 20
    sample_spacing = N // 10
    sample_gens = (n_samples - 1) * sample_spacing
    simlen = burnin * N + sample_gens
    print(f"{time.ctime()}, total simulation length:", simlen)

    g = build_demographic_model(N, sample_gens=sample_gens)
    demog = fwdpy11.discrete_demography.from_demes(g, burnin=burnin)

    pdict = {
        # Multiplicative selection model
        "gvalue": fwdpy11.Multiplicative(2.0),
        # Rates: (neutral, selection, recombination)
        "rates": (0.0, U, None),
        "nregions": [],
        "sregions": sel_regions,
        "recregions": rec_regions,
        "demography": demog,
        "simlen": simlen,
    }
    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation(N, L)
    recorder = Recorder(burnin=burnin * N, spacing=sample_spacing)
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


def get_diversity(ts, winsize=1000):
    """ """
    windows = np.arange(0, ts.sequence_length, 1000)
    if windows[-1] != ts.sequence_length:
        windows = np.concatenate((windows, [ts.sequence_length]))

    print(f"{time.ctime()}, getting diversity at each timepoint")
    # get sample times
    sample_sets = defaultdict(list)
    for samp in ts.samples():
        t = ts.node(samp).time
        sample_sets[t].append(samp)

    div = np.zeros((len(sample_sets.keys()), len(windows) - 1))

    # get diversity after simplifying to samples at each sample time
    for i, t in enumerate(sorted(sample_sets.keys())):
        ts_simp = ts.simplify(samples=sample_sets[t])
        branch_div = ts_simp.diversity(windows=windows, mode="branch")
        div[i] = branch_div
    print(f"{time.ctime()}, computed diversity")
    return div


def main(seed, Q):
    chrom = 2
    # set the length for the size you want to simulate
    chrom_lengths = {2: 315215}

    ts = runsim(
        seed, Ne=10000, chrom=chrom, u=1.5e-8, L=chrom_lengths[chrom], Q=Q, n_samples=100
    )

    div = get_diversity(ts)
    ave_div = div.mean(axis=0)
    return ave_div


if __name__ == "__main__":
    seed = int(sys.argv[1]) + 1
    Q = float(sys.argv[2])
    ave_div = main(seed, Q)

    output = f"ave_div.1kb.chr2.1-315215.seed_{seed}.Q_{Q}.np.gz"
    np.savetxt(output, ave_div)
