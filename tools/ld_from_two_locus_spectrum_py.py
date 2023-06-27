import moments.TwoLocus
import numpy as np

Ne = 1e4
s = -1e-4
u = 1 / 4 / Ne
r = 1e-5

gammas = [2 * Ne * s, 2 * Ne * s, 0] # selection on [AB, Ab, aB] haplotypes, resp
rho = 4 * Ne * r
theta = 4 * Ne * u

if np.abs(2 * Ne * s) > 20:
    print("Strong selection, approximations may break down")

print("2Ns   =", 2 * Ne * s)
print("rho   =", rho)
print("theta =", theta)


n_samples = 40

F = moments.TwoLocus.Demographics.equilibrium(
    n_samples, rho=rho, theta=theta, sel_params=gammas
)


# compute LD stats from F
print(f"D^2   = {F.D2():}")
print(f"Dz    = {F.Dz():}")
print(f"pi2   = {F.pi2():}")
print(f"D     = {F.D():}")

# single locus
fs = moments.Spectrum(
    moments.LinearSystem_1D.steady_state_1D(n_samples, gamma=2 * Ne * s) * theta
)

fs_left = moments.Spectrum(F[0, :, 0])
fs_right = moments.Spectrum(F[0, 0, :])

assert np.allclose(fs_left, fs)
print(f"pi(A) = {fs_left.pi()}")
print(f"pi(B) = {fs_right.pi()}")
