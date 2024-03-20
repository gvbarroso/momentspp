import numpy as np
import gzip
import sys

# run script as "python process_roulette.py {chrom}"

chrom = int(sys.argv[1])

# got these from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/
chrom_lengths = {2: 242193529, 8: 145138636}
L = chrom_lengths[chrom]

# https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate
# guess that rescaling factor has an extra factor of 2 (to make align with 1.2e-8)
rescale_fac = 1.015e-7 / 2

fname = f"{chrom}_rate_v5.2_TFBS_correction_all.vcf.gz"

rates = np.zeros(L, dtype=np.float32)
# for printing progress
last = -1

with gzip.open(fname, "rb") as fin:
    for line in fin:
        l = line.decode()
        if l.startswith("#"):
            continue
        spl = l.split()
        pos = int(spl[1])
        data = spl[7].split("MR=")[1]
        if ";" in data:
            data = data.split(";")[0]
        rates[pos - 1] += float(data)
        if pos % 1000000 == 0 and pos != last:
            last = pos
            print("At pos", pos, "of", L)

rates *= rescale_fac
print("Average rate:", np.mean(rates[rates > 0]))


def ave_rates(rates, spacing, tol=0.95):
    i = 0
    rates_ave = []
    end = False
    while not end:
        l = i * spacing
        r = (i + 1) * spacing
        nsites = np.sum(rates[l:r] > 0)
        sumrates = np.sum(rates[l:r])
        if nsites / spacing < tol:
            rates_ave.append(np.nan)
        else:
            rates_ave.append(sumrates / nsites)
        i += 1
        if r >= len(rates):
            end = True
    return np.array(rates_ave)


rates_10kb = ave_rates(rates, 10000)
rates_100kb = ave_rates(rates, 100000)
rates_1Mb = ave_rates(rates, 1000000)
