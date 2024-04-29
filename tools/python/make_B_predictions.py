import pandas
from scipy import interpolate
import numpy as np
import matplotlib.pylab as plt

df = pandas.read_csv("hr_time.csv.gz")

# just the expansion scenario
df = df[df["N1"] == 100000]

# parameters
u = 1e-8
r = 1e-8
L = 100000
Ne = 10000

gen = 50000

pi0 = 2 * Ne * u

s = -0.001

# set up windows and "exons" -- assertions are just to evenly grid the sequence length
spacing = 1000
assert L % spacing == 0

exon_L = 100
assert L % exon_L == 0

exon_mids = np.arange(exon_L / 2, L - exon_L / 2 + 1, exon_L)
assert len(exon_mids) == L // exon_L

xx = np.arange(0, L + 1, spacing)

# build cubic splines
s_vals = np.array(sorted(list(set(df["lookup_s"]))))
s_splines = {}
for s_val in s_vals:
    df_s = df[df["lookup_s"] == s_val]
    rs = np.array(df_s["lookup_r"])
    Hs = np.array(df_s[str(gen)])
    s_splines[s_val] = interpolate.CubicSpline(rs, Hs, bc_type="natural")


#### compute B values

# initial B value
B_vals = np.zeros(len(xx))
for i, x in enumerate(xx):
    B_vals[i] = np.prod((s_splines[s](np.abs(x - exon_mids) * r) / pi0) ** exon_L)


def B_interference_iteration(xx, B_func, u, r, s):
    """
    Note that this is really rough and could be made better in a lot of ways.....

    xx: positions to compute B value for
    B_func: cubic spline function interpolating B values from those computed at xx
    u: initial mutation rate
    r: initial recombination rate
    s: selection coefficient
    """
    B_mids = B_func(exon_mids)
    # adjusted cumulative r
    r_exon = B_mids * r
    r_cum = [0]
    for r_e in r_exon:
        r_cum.append(r_cum[-1] + r_e * exon_L)
    r_cum = np.array(r_cum)

    r_mids = (r_cum[1:] + r_cum[:-1]) / 2
    r_xx = np.concatenate(
        ([0], r_cum[1:].reshape(L // spacing, spacing // exon_L)[:, -1])
    )

    # adjusted mut within each exon
    u_exon = B_mids * u
    # adjusted s within each exon
    s_exon = B_mids * s
    ## note in future would pick weights between two closest?
    # s_closest = [s_vals[np.argmin(np.abs(s_vals - s_e))] for s_e in s_exon]
    s_low = np.array([s_vals[np.where(s_e >= s_vals)[0][-1]] for s_e in s_exon])
    s_high = np.array([s_vals[np.where(s_e < s_vals)[0][0]] for s_e in s_exon])
    p_high = (s_exon - s_low) / (s_high - s_low)
    p_low = 1 - p_high

    B_vals_interference = np.zeros(len(xx))
    for i, x in enumerate(xx):
        B_i = 1
        r_x = r_xx[i]
        r_dists = np.abs(r_mids - r_x)
        for j, (r_e, u_e) in enumerate(zip(r_dists, u_exon)):
            # B_i *= (s_splines[s_e](r_e) / pi0) ** (exon_L * u_e / u)
            p_l = p_low[j]
            p_h = p_high[j]
            s_l = s_low[j]
            s_h = s_high[j]
            B_l = (s_splines[s_l](r_e) / pi0) ** (exon_L * u_e / u)
            B_h = (s_splines[s_h](r_e) / pi0) ** (exon_L * u_e / u)
            B_i *= p_l * B_l + p_h * B_h
        B_vals_interference[i] = B_i
    return B_vals_interference


# run for some iterations
B_iter = {0: B_vals}
for i in range(1, 11):
    print("running iteration", i)
    B_func = interpolate.CubicSpline(xx, B_iter[i - 1], bc_type="natural")
    B_iter[i] = B_interference_iteration(xx, B_func, u, r, s)


# plot stuff
data = (
    np.loadtxt(
        f"../steady_state/Bvalues.steady_state.N_1000.s_{s*10}.u_1e-07.r_1e-07.txt"
    )
    / 0.0004
)
plt.plot(xx, data, "o--", label="Data")

for i in range(11):
    plt.plot(xx, B_iter[i], label="Iteration " + str(i))

plt.legend()
plt.ylabel("B")
plt.xlabel("Position")
plt.show()
