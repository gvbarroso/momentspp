"""
From the README:

Output from running moments++ in VERBOSE mode.

Individual transision matrices:

model_1_e_0_op_0.csv --> selection
model_1_e_0_op_1.csv --> recombination
model_1_e_0_op_2.csv --> mutation
model_1_e_0_op_3.csv --> drift

Full transition matrix M:

model_1_e_0_O_25_transitions.csv --> M = D * U * R * S

~~~~~~~~~~~~

Install the requirements: `pip install -r requirements.txt`
"""

import numpy as np
from scipy.sparse.linalg import factorized
from scipy.sparse import csc_matrix
from scipy.sparse import identity


def load_matrices():
    ## TODO: arguments for which model to load
    sel = np.loadtxt("model_1_e_0_op_0.csv", delimiter=",")
    rec = np.loadtxt("model_1_e_0_op_1.csv", delimiter=",")
    mut = np.loadtxt("model_1_e_0_op_2.csv", delimiter=",")
    drift = np.loadtxt("model_1_e_0_op_3.csv", delimiter=",")
    return sel, rec, mut, drift


def get_parameters():
    ## TODO: parse the demes file
    ## TODO: fix moments++ handling of metadata to be demes-compliant
    s = -0.0001
    r = 1e-8
    u = 1e-6
    Ne = 1e4
    return Ne, s, r, u


def build_transition_matrix(nu=1):
    ## nu is the relative size change
    sel, rec, mut, drift = load_matrices()
    Ne, s, r, u = get_parameters()

    # Rescale parameters by 2Ne, and drift by nu, the relative size.
    sel *= 2 * Ne
    rec *= 2 * Ne
    mut *= 2 * Ne
    drift *= 2 * Ne / nu

    # Identify and remove the identity row/column
    I_idx = np.where(np.sum(np.abs(sel + drift + rec + mut), axis=1) == 0)[0]
    mut_1loc = mut[:, I_idx]

    P = sel + drift + rec + mut

    P = np.delete(P, I_idx, axis=1)
    P = np.delete(P, I_idx, axis=0)
    mut_1loc = np.delete(mut_1loc, I_idx)

    # We'll make use of the scipy sparse matrix representations
    P_sp = csc_matrix(P)
    return P_sp, mut_1loc


def steady_state():
    P_sp, mut_1loc = build_transition_matrix(nu=1)
    # This comes from y' = (I+P).y + mut, and then setting y'=y
    y = factorized(P_sp)(-mut_1loc)
    return y


def integrate(y, nu=1, T=1, dt=0.001):
    P_sp, mut_1loc = build_transition_matrix(nu=nu)
    t_elapsed = 0
    while t_elapsed + dt < T:
        dt = min(dt, T - t_elapsed)
        # set up the forward and backward transitions for the Crank-Nicolson algorithm
        A_fd = identity(P_sp.shape[0], format="csc") + dt / 2 * P_sp
        A_bd = factorized(identity(P_sp.shape[0], format="csc") - dt / 2 * P_sp)
        # forward step
        y = A_fd.dot(y) + dt * mut_1loc
        # backward step
        y = A_bd(y)
        # increment time
        t_elapsed += dt
    return y


def write_output(outfile, y):
    fname = outfile.split(".txt")[0] + ".ode.txt"
    i = 0
    with open(outfile, "r") as fin, open(fname, "w+") as fout:
        for line in fin:
            l = line.strip()
            if l.startswith("I"):
                val = 1
            else:
                val = y[i]
                i += 1
            l += ", " + str(val) + "\n"
            fout.write(l)


if __name__ == "__main__":
    # steady state
    y = steady_state()

    # size change with relative size nu, for time T (measured in 2Nanc generations)
    nu = 10
    T = 0.05
    dt = 0.002  # a smaller dt will be more accurate, but will involve more time-steps
    y2 = steady_state()
    y2 = integrate(y2, nu=nu, T=T, dt=0.005)

    # print(y2 - y)
    # write_output("model_1_O_25_expectations.txt", y)
