import numpy as np
import numpy.random as npr
import moments.LD

rng = npr.default_rng()


# get statistics
def getD(X):
    return X[:, 0] * X[:, 3] - X[:, 1] * X[:, 2]


def getp(X):
    return X[:, 0] + X[:, 1]


def getq(X):
    return X[:, 0] + X[:, 2]


def getHl(Xl):
    return Xl * (1 - Xl)


def getHr(Xr):
    return Xr * (1 - Xr)


def getDsqr(X):
    D = getD(X)
    return D ** 2


def getDz(X):
    D = getD(X)
    p = getp(X)
    q = getq(X)
    return D * (1 - 2 * p) * (1 - 2 * q)


def getpi2(X):
    p = getp(X)
    q = getq(X)
    return p * (1 - p) * q * (1 - q)


# evolve system
def mutate(Xl, Xr, X, Ne, u, L=1):
    # new single mutations at left locus
    num_left_mut = rng.poisson(2 * Ne * u * L)
    Xl = np.concatenate((Xl, np.ones(num_left_mut) / 2 / Ne))

    # new single mutations at right locus
    num_right_mut = rng.poisson(2 * Ne * u * L)
    Xr = np.concatenate((Xr, np.ones(num_right_mut) / 2 / Ne))

    # mutate against segregating left loci in Xl
    for _ in range(L):
        new_muts = rng.random(len(Xl)) < 2 * Ne * u
        if np.any(new_muts):
            for j in np.where(new_muts)[0]:
                if rng.random() < Xl[j]:
                    # falls on bA background
                    X = np.concatenate(
                        (X, [[1 / 2 / Ne, Xl[j] - 1 / 2 / Ne, 0, 1 - Xl[j]]])
                    )
                else:
                    # falls on ab background
                    X = np.concatenate(
                        (X, [[0, Xl[j], 1 / 2 / Ne, 1 - Xl[j] - 1 / 2 / Ne]])
                    )

    # mutate against segregating right loci in Xr
    for _ in range(L):
        new_muts = rng.random(len(Xr)) < 2 * Ne * u
        if np.any(new_muts):
            for j in np.where(new_muts)[0]:
                if rng.random() < Xr[j]:
                    # falls on aB background
                    X = np.concatenate(
                        (X, [[1 / 2 / Ne, 0, Xr[j] - 1 / 2 / Ne, 1 - Xr[j]]])
                    )
                else:
                    # falls on ab background
                    X = np.concatenate(
                        (X, [[0, 1 / 2 / Ne, Xr[j], 1 - Xr[j] - 1 / 2 / Ne]])
                    )

    return Xl, Xr, X


def recombine(X, r):
    D = getD(X)
    fac = np.outer(D, [-1, 1, 1, -1])
    Xrec = X + r * fac
    return Xrec


def select(Xl, s):
    # s < 0 implies negative selection
    Xl = ((1 + s) * Xl) / (Xl * (1 + s) + 1 - Xl)
    return Xl


def drift(Xl, Xr, X, Ne):
    assert np.allclose(X.sum(axis=1), 1)
    Xl = rng.binomial(2 * Ne, Xl) / 2 / Ne
    Xr = rng.binomial(2 * Ne, Xr) / 2 / Ne
    # sometime floating point stuff gives machine precision slightly outside of [0,1]
    X[X < 0] = 0
    X[X > 1] = 1
    X = rng.multinomial(2 * Ne, X) / 2 / Ne
    return Xl, Xr, X


def evolve(Xl, Xr, X, Ne, u, r, L, s=0):
    X = recombine(X, r)
    Xl = select(Xl, s)
    Xl, Xr, X = drift(Xl, Xr, X, Ne)
    Xl, Xr, X = mutate(Xl, Xr, X, Ne, u, L=L)
    return Xl, Xr, X


def cleanup(Xl, Xr, X):
    # remove fixed or lost in Xl and Xr
    fixed = np.logical_or(Xl == 0, Xl == 1)
    Xl = Xl.compress(1 - fixed)
    fixed = np.logical_or(Xr == 0, Xr == 1)
    Xr = Xr.compress(1 - fixed)
    # remove any with p or q at 0 or 1
    p = getp(X)
    q = getq(X)
    fixed_p = np.logical_or(p == 0, p == 1)
    fixed_q = np.logical_or(q == 0, q == 1)
    fixed = np.logical_or(fixed_p, fixed_q)
    X = X.compress(1 - fixed, axis=0)
    return Xl, Xr, X


# parameters
Ne = 2000
u = 1e-6
r = 1e-4
s = -1e-4
rho = 4 * Ne * r
theta = 4 * Ne * u
print("theta   = ", theta)
print("theta^2 = ", theta ** 2)

y = moments.LD.Demographics1D.snm(theta=theta, rho=rho)
H = y.H()[0]
LD = y.LD()[0]
print(f"E[Hl]   = {H/2:0.4f}")
print(f"E[Hr]   = {H/2:0.4f}")
print(f"E[D^2]  = {LD[0]:0.8f}")
print(f"E[Dz]   = {LD[1]:0.8f}")
print(f"E[pi2]  = {LD[2]:0.8f}")

# X stores haplotype frequencies
# M has shape nx2, indicates which sites have mutated, n=len(X)
# F keeps track of fixed entries in X, which can be replaced
X = np.empty((0, 4))
Xl = np.empty(0)
Xr = np.empty(0)

L = 50
print("    L   = ", L)

# burn in for 40 * Ne generations
for gen in range(40 * Ne):
    Xl, Xr, X = evolve(Xl, Xr, X, Ne, u, r, L, s)
    Xl, Xr, X = cleanup(Xl, Xr, X)

# start collecting statistics
# run for sum number of generations
Hl = 0
Hr = 0
D2 = 0
Dz = 0
pi2 = 0

gens = 100_000_000
c = 0
for gen in range(gens):
    Xl, Xr, X = evolve(Xl, Xr, X, Ne, u, r, L, s)
    Xl, Xr, X = cleanup(Xl, Xr, X)
    Hl += getHl(Xl).sum()
    Hr += getHr(Xr).sum()
    D2 += getDsqr(X).sum()
    Dz += getDz(X).sum()
    pi2 += getpi2(X).sum()
    c += 1
    if c % (100 * Ne) == 0:
        print(
            f"gen = {c}",
            f"Hl = {Hl/c/L:0.4f}",
            f"Hr = {Hr/c/L:0.4f}",
            f"D^2 = {D2/c/L**2:0.8f}",
            f"Dz = {Dz/c/L**2:0.8f}",
            f"pi2 = {pi2/c/L**2:0.8f}",
        )
