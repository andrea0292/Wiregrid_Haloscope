"""Boost factor as an overlap integral with the receiving mode (SM section and figure fig_overlap).

Receiving mode E_rec(s): the steady-state field when a unit plane wave e^{-ik(s-L)}, polarized along the wires,
is incident from the detector side on the passive stack (PEC mirror at s=0, impedance sheets at s=jd).
Reciprocity / 1D Green's function:  beta = (ik/2) * int_0^L E_rec(s) ds   (source confined to the stack,
beta referenced to the emitting face as in End Matter A).  Checked here against the direct solver of
reproduce_final.py (beta_conf) over the spectrum.

Usage (from code/validation):  python overlap_receiving_mode.py  ->  ../../results/fig_overlap.(pdf|png), printed numbers
"""
import os, sys, importlib.util, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..', '..')
def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m
sys.argv = sys.argv[:1]
rf = _load('rf', os.path.join(ROOT, 'code', 'reproduce_final.py'))
fr = _load('fr', os.path.join(ROOT, 'code', 'figures_reach.py'))
c0, mu0, hP, eV = rf.c0, rf.mu0, 6.62607015e-34, 1.602176634e-19

# 1-2 meV benchmark cartridge at its geometric centre (Table I)
a = 237.7e-6; r = 0.03*a; d = a; N = 11; L = N*d
Rs = lambda nu: 83.6e-3*(nu/250e9)**(2/3)
nu_c = 1.41421e-3*eV/hP

def rec_field(nu, s):
    """E_rec(s) for an incident wave e^{-ik(s-L)}; returns (E_rec, k, r) with r the reflection amplitude at s=L."""
    k = 2*np.pi*nu/c0; w = 2*np.pi*nu
    Zg = rf.grid_impedance(nu, a, r, Rs(nu))[0]; g = -1j*w*mu0/Zg
    xs = d*np.arange(1, N+1)
    A, B = 1/(2j), -1/(2j); coef = [(A, B)]            # u = sin(ks) between mirror and first grid
    for xj in xs:                                      # E continuous, dE/ds jumps by g E at each grid
        e = np.exp(1j*k*xj); u = A*e + B/e; du = 1j*k*(A*e - B/e) + g*u
        A = (1j*k*u + du)/(2j*k)/e; B = (1j*k*u - du)/(2j*k)*e; coef.append((A, B))
    A_out, B_out = coef[-1]
    norm = B_out*np.exp(-1j*k*L)                       # coefficient of e^{-ik(s-L)}
    idx = np.searchsorted(xs, s, side='left')
    Aa = np.array([coef[i][0] for i in idx]); Bb = np.array([coef[i][1] for i in idx])
    rr = A_out*np.exp(1j*k*L)/norm
    return (Aa*np.exp(1j*k*s) + Bb*np.exp(-1j*k*s))/norm, k, rr

def beta_overlap(nu, n=400):
    s = np.linspace(0, L, n*N + 1); E, k, _ = rec_field(nu, s)
    return 1j*k/2*np.trapezoid(E, s)

def beta_direct(nu):  # direct solver, phase referenced to the face s = L
    return rf.beta_conf(nu, a, r, d, N, Rs)*np.exp(2j*np.pi*nu/c0*L)

def find_peak(f0, h=2e-3):
    from scipy.optimize import minimize_scalar
    return minimize_scalar(lambda v: -abs(beta_overlap(v))**2, bracket=(f0*(1-h), f0, f0*(1+h))).x

if __name__ == '__main__':
    f = np.linspace(0.95, 1.05, 801)*nu_c
    dev = max(abs(beta_overlap(v)/beta_direct(v) - 1) for v in f[::8])
    nu1 = find_peak(nu_c); nu2 = find_peak(1.0928*nu_c); lam1 = c0/nu1
    for lab, nu in (('nu1', nu1), ('nu2', nu2)):
        s = np.linspace(0, L, 200001); E, k, rr = rec_field(nu, s)
        print(f'{lab}: {nu/1e9:.2f} GHz = {nu*hP/eV*1e3:.4f} meV, |beta|^2 = {abs(beta_overlap(nu))**2:.1f}, '
              f'|r|^2 = {abs(rr)**2:.4f}, max|E_rec| = {abs(E).max():.2f}, |<E_rec>| = {abs(E.mean()):.2f}')
    print(f'max |beta_overlap/beta_direct - 1| over 0.95-1.05 nu_c: {dev:.1e}; lambda1 = {lam1*1e3:.4f} mm, L/lambda1 = {L/lam1:.3f}, '
          f'kd = {2*np.pi*d/lam1:.3f}, Rs(nu1) = {Rs(nu1)*1e3:.1f} mOhm')

    fr.style(); C = fr.SET1
    s = np.linspace(0, L + 0.65*lam1, 30001); S = s/lam1
    E1, k1, _ = rec_field(nu1, s); E2, k2, _ = rec_field(nu2, s)
    Em = np.exp(-1j*k1*(s - L)) - np.exp(1j*k1*(s + L))     # bare mirror, same incident normalization
    fig, axs = plt.subplots(2, 1, figsize=(7, 9.6), sharex=True)
    for ax in axs:
        ax.axvspan(0, L/lam1, color='0.93', zorder=0)
        for xg in d*np.arange(1, N+1)/lam1: ax.axvline(xg, color='0.6', lw=0.9, zorder=1)
        ax.axvline(0, color='k', lw=4, zorder=2)
    ax = axs[0]
    ax.plot(S, abs(E1), color=C[0], label=r'fundamental mode, $\nu_1$')
    ax.plot(S, abs(E2), color=C[1], label=r'next resonance, $\nu_2$')
    ax.plot(S, abs(Em), color='0.35', ls='--', lw=1.8, label=r'mirror alone, $\nu_1$')
    ax.set_ylim(0, 10.2); ax.set_ylabel(r'$|E_{\rm rec}(s)|$')
    ax.legend(fontsize=16, loc='upper right', handlelength=1.6, borderaxespad=0.3); ax.text(0.02, 0.93, '(a)', transform=ax.transAxes, fontsize=19)
    ax = axs[1]
    for E, k, col, ls in ((E1, k1, C[0], '-'), (E2, k2, C[1], '-'), (Em, k1, '0.35', '--')):
        src = np.where(s <= L, E, 0)                                    # source confined to the stack
        cum = 1j*k/2*np.concatenate([[0], np.cumsum(0.5*(src[1:] + src[:-1])*np.diff(s))])
        ax.plot(S, abs(cum)**2, color=col, ls=ls, lw=2.5 if ls == '-' else 1.8)
    ax.set_yscale('log'); ax.set_ylim(1e-2, 5e3)
    ax.set_ylabel(r'$\left|\frac{ik}{2}\int_0^{s}E_{\rm rec}\,\mathrm{d}s^\prime\right|^2$')
    ax.set_xlabel(r'$s/\lambda_1$'); ax.set_xlim(-0.04, S[-1])
    ax.text(0.02, 0.93, '(b)', transform=ax.transAxes, fontsize=19)
    ax.text(S[-1] - 0.04, abs(beta_overlap(nu1))**2*1.35, r'$|\beta_1|^2=1651$', color=C[0], ha='right', fontsize=17)
    ax.text(S[-1] - 0.04, abs(beta_overlap(nu2))**2*1.35, r'$|\beta_2|^2=108$', color=C[1], ha='right', fontsize=17)
    fig.tight_layout(h_pad=0.6)
    out = os.path.join(ROOT, 'results'); os.makedirs(out, exist_ok=True)
    for ext in ('pdf', 'png'): fig.savefig(os.path.join(out, 'fig_overlap.' + ext), dpi=200)
