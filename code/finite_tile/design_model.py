"""
Design model for a resonant inductive-grid (wire-medium, epsilon-near-zero) booster
for axion-photon conversion at ~0.5-5 meV.

1D model along x (emission direction). Static field B0 || z, wires || z, E_z polarisation.
Units: SI for electromagnetics; axion source normalised so that a bare mirror emits
|E_out| = E0 (boost beta = |E_out|/E0, beta^2 = P/P_dish).

Field in a vacuum gap:   E(x) = -E0 + A exp(i k x) + B exp(-i k x),  time dependence exp(-i w t)
Mirror at x=0 (perfect conductor):  E(0) = 0
Grid at x_j with sheet impedance Z_g:  E continuous;  (i/(w mu0)) [E'(x+) - E'(x-)] = E(x_j)/Z_g
Output region x > x_N: E = -E0 + T exp(i k x)  (outgoing only).
Grid impedance (E || thin wires, radius r, pitch a; Marcuvitz):
    Z_g = -i w L_s + R_g,  L_s = mu0 a/(2 pi) ln(a/(2 pi r)),  R_g = R_s a/(2 pi r)
Homogenised limit (spacing d << lambda): eps = 1 - wp^2/(w(w + i Gamma)),
    wp^2 = 1/(eps0 d L_s),  Gamma = R_g/L_s,  Q = wp/Gamma = w mu0 r ln(a/2 pi r)/R_s.
"""
import numpy as np

c = 2.99792458e8
mu0 = 4e-7 * np.pi
eps0 = 1.0 / (mu0 * c**2)
Z0 = mu0 * c
h = 6.62607015e-34
kB = 1.380649e-23
eV = 1.602176634e-19

def nu_of_m(m_eV):
    return m_eV * eV / h

def Rs_copper(nu, scale=1.0):
    """Cryogenic high-RRR Cu surface resistance in the anomalous-skin regime,
    anchored to R_s = 2.9 mOhm at 3.9 GHz with nu^(2/3) scaling (Posen et al.,
    arXiv:2201.10733, as used in the phase-2 vetting). 'scale' multiplies for roughness/field."""
    return 2.9e-3 * (nu / 3.9e9) ** (2.0 / 3.0) * scale

def grid_impedance(nu, a, r, Rs):
    w = 2 * np.pi * nu
    Ls = mu0 * a / (2 * np.pi) * np.log(a / (2 * np.pi * r))
    Rg = Rs * a / (2 * np.pi * r)
    return -1j * w * Ls + Rg, Ls, Rg

def plasma_freq(a, r, d):
    Ls = mu0 * a / (2 * np.pi) * np.log(a / (2 * np.pi * r))
    return 1.0 / np.sqrt(eps0 * d * Ls) / (2 * np.pi)

def boost_grid_stack(nu, positions, Zg):
    """beta (complex E_out/E0) for grids at positions x_1<...<x_N (m) above a mirror at x=0.
    Zg: scalar or array (per grid) sheet impedance at frequency nu. Solves the linear system
    for gap amplitudes directly (exact within the 1D sheet model)."""
    w = 2 * np.pi * nu
    k = w / c
    x = np.asarray(positions, float)
    N = len(x)
    Zg = np.broadcast_to(np.asarray(Zg, complex), (N,))
    # unknowns: (A_0,B_0) for gap [0,x1], ..., (A_{N-1},B_{N-1}) for gap [x_{N-1},x_N], T for output
    nunk = 2 * N + 1
    M = np.zeros((nunk, nunk), complex)
    rhs = np.zeros(nunk, complex)
    E0 = 1.0
    row = 0
    # mirror: -E0 + A0 + B0 = 0
    M[row, 0] = 1; M[row, 1] = 1; rhs[row] = E0; row += 1
    for j in range(N):
        xj = x[j]
        e_p, e_m = np.exp(1j * k * xj), np.exp(-1j * k * xj)
        iL = 2 * j  # left gap index
        # field just left: -E0 + A e_p + B e_m ; derivative: ik(A e_p - B e_m)
        if j < N - 1:
            iR = 2 * (j + 1)
            # continuity: left = right (the -E0 cancels)
            M[row, iL] = e_p; M[row, iL + 1] = e_m; M[row, iR] = -e_p; M[row, iR + 1] = -e_m
            row += 1
            # jump: (i/(w mu0)) [ik(AR e_p - BR e_m) - ik(AL e_p - BL e_m)] = (-E0 + AL e_p + BL e_m)/Z
            fac = (1j / (w * mu0)) * 1j * k
            M[row, iR] = fac * e_p; M[row, iR + 1] = -fac * e_m
            M[row, iL] = -fac * e_p - e_p / Zg[j]; M[row, iL + 1] = fac * e_m - e_m / Zg[j]
            rhs[row] = -E0 / Zg[j]
            row += 1
        else:
            iT = 2 * N
            # continuity with output: AL e_p + BL e_m = T e_p
            M[row, iL] = e_p; M[row, iL + 1] = e_m; M[row, iT] = -e_p; row += 1
            fac = (1j / (w * mu0)) * 1j * k
            M[row, iT] = fac * e_p
            M[row, iL] = -fac * e_p - e_p / Zg[j]; M[row, iL + 1] = fac * e_m - e_m / Zg[j]
            rhs[row] = -E0 / Zg[j]
            row += 1
    sol = np.linalg.solve(M, rhs)
    return sol[2 * N] / E0

def boost_homogeneous_slab(nu, L, nu_p, Q, eps_b=1.0):
    """Exact 1D slab of Drude medium eps = eps_b - wp^2/(w(w+i Gamma)) on a mirror, thickness L.
    Field in slab: -E0/eps + A e^{i n k x} + B e^{-i n k x}; vacuum outside."""
    w = 2 * np.pi * nu; wp = 2 * np.pi * nu_p; G = wp / Q
    eps = eps_b - wp**2 / (w * (w + 1j * G))
    n = np.sqrt(eps + 0j)
    if n.imag < 0: n = -n
    k = w / c
    # unknowns A,B,T
    ep, em = np.exp(1j * n * k * L), np.exp(-1j * n * k * L)
    M = np.array([[1, 1, 0],
                  [ep, em, -np.exp(1j * k * L)],
                  [n * ep, -n * em, -np.exp(1j * k * L)]], complex)
    rhs = np.array([1 / eps, (1 / eps - 1), 0.0], complex)
    A, B, T = np.linalg.solve(M, rhs)
    return T

def nbar(nu, T):
    x = h * nu / (kB * T)
    return 1.0 / np.expm1(x) if x < 700 else np.exp(-x)
