"""Small numbers quoted in the Introduction and End Matter B that are not produced by the main solvers.

  1. Introduction: unloaded Q of a square copper-wire lattice (ALPHA lattice expression, ln(a/2 pi r) + 0.5275),
     Q_u ~ 510 (m_a/meV)^(-2/3), ~14 % above Q_m of Eq. (2), and the field-amplitude 1/e length
     l_amp = lambda sqrt(2 Q_u) / 2 pi  (6.3 mm at 1 meV, 0.4 mm at 8 meV).
  2. End Matter B: frequency-integrated boost per element, int |beta|^2 dnu/nu over 0.5-1.5 nu_0,
     for sapphire disks (eps = 9.4, half-wave configuration) and for the 1-2 meV grid cartridge (lossless and copper).
  3. End Matter B: scan integral I ~ pi |beta_pk|^4 / (4 Q_r) of the Table I cartridges and the resulting
     dish-vs-grid improvement in chi, [(I / 5 ln 2) (1/2)^2 (eta_ap xi)^2]^(1/4).
Usage (from code/validation): python endmatter_checks.py
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from reproduce_final import chi_their   # halo factor xi of SM S2

c, Z0 = 299792458.0, 376.73
def Rs_paper(f): return 83.6e-3 * (f / 250e9) ** (2 / 3)


def beta_layers(k, eps, z, Y):
    """1D mirror-backed multilayer driven by a uniform axion source (E0 = 1) confined to [0, z[-1]].
    eps[j]: permittivity between z[j] and z[j+1]; Y: {boundary index: normalized sheet admittance}.
    Returns the outgoing amplitude beta (bare mirror: |beta| = 1)."""
    nL = len(eps); n = np.sqrt(np.array(eps, complex))
    M = np.zeros((2 * nL + 1, 2 * nL + 1), complex); b = np.zeros(2 * nL + 1, complex)
    M[0, 0] = M[0, 1] = 1; b[0] = 1 / eps[0]; row = 1
    for j in range(nL):
        s = z[j + 1]; e = np.exp(1j * n[j] * k * s); y = Y.get(j + 1, 0); pj = -1 / eps[j]
        E = [(2 * j, e), (2 * j + 1, 1 / e)]; dE = [(2 * j, 1j * n[j] * k * e), (2 * j + 1, -1j * n[j] * k / e)]
        if j < nL - 1:
            q = n[j + 1]; e2 = np.exp(1j * q * k * s)
            for i, v in E: M[row, i] += v
            M[row, 2 * j + 2] -= e2; M[row, 2 * j + 3] -= 1 / e2; b[row] = -1 / eps[j + 1] - pj; row += 1
            M[row, 2 * j + 2] += 1j * q * k * e2; M[row, 2 * j + 3] -= 1j * q * k / e2
        else:
            for i, v in E: M[row, i] += v
            M[row, 2 * nL] = -1; b[row] = -pj; row += 1
            M[row, 2 * nL] = 1j * k
        for i, v in dE: M[row, i] -= v
        for i, v in E: M[row, i] += 1j * k * y * v
        b[row] = -1j * k * y * pj; row += 1
    return np.linalg.solve(M, b)[2 * nL]


def per_element(fn, f0, nel):
    f = np.linspace(0.5, 1.5, 6001) * f0
    return np.trapezoid(np.array([abs(fn(2 * np.pi * x / c)) ** 2 for x in f]) / f, f) / nel


if __name__ == "__main__":
    # 1. Q_u and l_amp (resonance at the plasma frequency, r/a = 0.03, d = a)
    for m in (1.0, 8.0):
        f = m * 241.8e9; k = 2 * np.pi * f / c; ell = np.log(1 / (2 * np.pi * 0.03))
        a = np.sqrt(2 * np.pi / ell) / k; Qm = k * 0.03 * a * ell / (Rs_paper(f) / Z0)
        Qu = Qm * np.sqrt((ell + 0.5275) / ell); lam = c / f
        print(f"m = {m:.0f} meV: Q_m = {Qm:.0f}, Q_u = {Qu:.0f} (Q_u/Q_m = {Qu/Qm:.3f}), "
              f"l_amp = {lam*np.sqrt(2*Qu)/(2*np.pi)*1e3:.2f} mm")
    # 2. per-element frequency-integrated boost
    a = 237.69e-6; r = 0.03 * a; ell = np.log(a / (2 * np.pi * r)); N = 11; f0 = 341.6e9
    def grid(scale):
        def fn(k):
            f = k * c / (2 * np.pi); zg = scale * Rs_paper(f) / Z0 * a / (2 * np.pi * r) - 1j * k * a * ell / (2 * np.pi)
            return beta_layers(k, [1.0] * N, a * np.arange(N + 1), {i: 1 / zg for i in range(1, N + 1)})
        return fn
    print(f"1-2 meV grids: int|beta|^2 dnu/nu per grid = {per_element(grid(0.0), f0, N):.2f} (lossless), "
          f"{per_element(grid(1.0), f0, N):.2f} (copper R_s)")
    for nd in (5, 10, 20):
        ed = 9.4; lam = c / 20e9; eps, z = [], [0.0]
        for _ in range(nd):
            eps += [1.0, ed]; z += [z[-1] + lam / 2, z[-1] + lam / 2 + lam / (2 * np.sqrt(ed))]
        print(f"sapphire, {nd:2d} half-wave disks: per disk = "
              f"{per_element(lambda k: beta_layers(k, eps, np.array(z), {}), 20e9, nd):.2f}")
    # 3. scan integral of Table I and dish comparison (xi of SM S2 at the geometric cartridge centres)
    tab = [("0.25-0.5", 0.25, 585, 5530, 0.10), ("0.5-1", 0.5, 379, 3040, 0.10), ("1-2", 1.0, 239, 1650, 0.10),
           ("2-4", 2.0, 146, 883, 0.10), ("4-8", 4.0, 123, 641, 0.05)]
    for band, lo, Qr, b2, W in tab:
        I = np.pi * b2 ** 2 / (4 * Qr); xi = chi_their(lo * np.sqrt(2), w=W / 2, L=W)
        print(f"{band:9s} I = {I:8.0f}   xi = {xi:.2f}   chi improvement over a dish = "
              f"{(I / (5*np.log(2)) * 0.25 * (0.65*xi) ** 2) ** 0.25:.2f}")
