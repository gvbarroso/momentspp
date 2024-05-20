"""
Run as: python bmap_prediction.py {s} {n_iter}

n_iter is optional

Functions for predicting B-maps from a data table. At a minumum, these should have
columns:

Na uL uR r s t pi0 piL piR

Na: ancestral Ne
uL: mutation rate at left locus
uR: mutation rate at right locus (often a single value, 1e-8)
r: recombination rate
s: selection coefficient
t: time since steady state
pi0: neutral expected diversity (without selection or BGS)
piL: predicted diversity at left locus
piR: predicted diversity at right locus
"""

import pandas
from scipy import interpolate
import numpy as np
import sys
import matplotlib.pylab as plt


def generate_cubic_splines(df_sub):
    """
    df_sub is the dataframe subsetted to a single Na, uR and t

    Cubic spline functions are created over r, for each combination of
    uL and s, returning the fractional reduction based on piR and pi0.

    This function returns the (sorted) arrays of uL and s and the
    dictionary of cubic spline functions with keys (uL, s).
    """
    # Check that only
    assert len(np.unique(df_sub["Na"])) == 1
    assert len(np.unique(df_sub["uR"])) == 1
    assert len(np.unique(df_sub["t"])) == 1

    # Get arrays of selection and mutation values
    s_vals = np.array(sorted(list(set(df_sub["s"]))))
    u_vals = np.array(sorted(list(set(df_sub["uL"]))))

    # Store cubic splines of fractional reduction for each pair of u and s, over r
    splines = {}
    for u in u_vals:
        for s in s_vals:
            key = (u, s)
            # Subset to given s and u values
            df_s_u = df_sub[(df_sub["s"] == s) & (df_sub["uR"] == u)]
            rs = np.array(df_s_u["r"])
            piR = np.array(df_s_u["Hr"])
            pi0 = np.array(df_s_u["pi0"])
            fs = piR / pi0
            splines[key] = interpolate.CubicSpline(rs, fs, bc_type="natural")

    return u_vals, s_vals, splines


def get_B_uniform_rmap(
    x, x_elem, L_elem, u_vals, s_vals, splines, u=1e-8, s=-0.001, r=1e-8, max_r=1e-2
):
    """
    For a given position x within the region, compute the reduction in diversity
    due to linked elements. The midpoints of these elements are given by the array
    x_elem, and L_elem is the length of each of these. We assume each element has
    the same length.

    x is the position in the chromosome
    x_elem is the midpoint of each element
    L_elem is the length of the element
    u_vals, s_vals, and splines is the output of generate_cubic_splines()

    with max_r, we don't account for any effects from regions farther away than
    this distance, as the cubic splines are not stable beyond the endpoints

    returns a B map that as a cubic spline interpolated function
    """
    # remove distances greater than max r value
    dists = np.abs(x - x_elem) * r
    dists = dists[dists <= max_r]
    uL = 1e-8
    assert len(u_vals) == 1 and u_vals[0] == uL
    B = u / uL

    if s in s_vals:
        return np.prod(splines[(uL, s)](dists) ** (B * L_elem))
    else:
        if s < s_vals[0] or s > s_vals[-1]:
            raise ValueError(f"s={s} is outside of computed s-value range")
        s0 = s_vals[np.where(s > s_vals)[0][-1]]
        s1 = s_vals[np.where(s > s_vals)[0][-1] + 1]
        p1 = (s - s0) / (s1 - s0)
        p0 = 1 - p1
        return p0 * np.prod(splines[(uL, s0)](dists) ** (B * L_elem)) + p1 * np.prod(
            splines[(uL, s1)](dists) ** (B * L_elem)
        )


def get_Bmap_uniform(
    xx, x_elem, L_elem, u_vals, s_vals, splines, s=-0.001, u=1e-8, r=1e-8
):
    B_vals = np.zeros(len(xx))
    for i, x in enumerate(xx):
        B_vals[i] = get_B_uniform_rmap(
            x, x_elem, L_elem, u_vals, s_vals, splines, s=s, u=u, r=r
        )
    Bmap = interpolate.CubicSpline(xx, B_vals, bc_type="natural")
    return Bmap


def get_Bmap_adjusted(
    xx,
    x_elem,
    L_elem,
    u_vals,
    s_vals,
    splines,
    Bmap,
    u=1e-8,
    s=-0.001,
    r=1e-8,
    max_r=1e-2,
):
    """
    Same arguments as get_B_uniform_rmap(), but pass an Bmap to adjust parameters.

    returns a B map that as a cubic spline interpolated function
    """
    uL = 1e-8
    # adjust parameters based on input Bmap
    r_xx = update_rmap(Bmap, xx, r)
    r_elem = update_rmap(Bmap, x_elem, r)
    B_elem = Bmap(x_elem)
    s_elem = B_elem * s

    # build the adjusted Bmap
    B_vals_adjust = np.ones(len(xx))
    for i, r_x in enumerate(r_xx):
        r_dists = np.abs(r_elem - r_xx[i])
        for j, r_e in enumerate(r_dists):
            if r_e > max_r:
                # don't compute for distances greater than max_r
                continue
            B_e = B_elem[j]
            s_e = s_elem[j]
            L_e = L_elem * B_e * u / uL
            if s_e in s_vals:
                B_vals_adjust[i] *= splines[(uL, s_e)](r_e) ** L_e
            else:
                if s_e < s_vals[0] or s_e > s_vals[-1]:
                    raise ValueError(f"s={s_e} is outside of computed s-value range")
                s0 = s_vals[np.where(s_e > s_vals)[0][-1]]
                s1 = s_vals[np.where(s_e > s_vals)[0][-1] + 1]
                p1 = (s_e - s0) / (s1 - s0)
                p0 = 1 - p1
                B0 = splines[(uL, s0)](r_e) ** L_e
                B1 = splines[(uL, s1)](r_e) ** L_e
                B_vals_adjust[i] *= p0 * B0 + p1 * B1
    Bmap_adjust = interpolate.CubicSpline(xx, B_vals_adjust, bc_type="natural")
    return Bmap_adjust


def update_rmap(Bmap, xs, r):
    """
    Bmap: interpolated B value function
    xs: an array of positions
    r: initial (constant) recombination rate

    returns the cumulative r values at each x in xs
    """
    rs = np.zeros(len(xs))
    for i, x in enumerate(xs):
        rs[i] = Bmap.integrate(0, x) * r
    return rs


def update_umap(Bmap, x_elem, u):
    B_elem = Bmap(x_elem)
    u_elem = B_elem * u
    return u_elem
