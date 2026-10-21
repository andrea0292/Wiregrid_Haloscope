"""Fig. S1 (version with stack modes m = 1, 2, 3): real part of the receiving mode and the running overlap integral.

Uses the solver in overlap_receiving_mode.py (1-2 meV cartridge of Table I, impedance-sheet model).
  m = 1 : fundamental resonance nu_1 (one internal half-wavelength across the stack)
  m = 2 : next stack resonance, antisymmetric ("dark": its overlap with the uniform axion source cancels)
  m = 3 : next resonance with nonzero boost
Usage (from code/validation): python overlap_modes_figure.py -> ../../results/fig_overlap_modes.(pdf|png)
"""
import os, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import overlap_receiving_mode as OV

nu_c = OV.nu_c; L = OV.L; c0 = OV.c0
nu1 = OV.find_peak(nu_c); nu3 = OV.find_peak(1.0928*nu_c)
nu2 = minimize_scalar(lambda v: abs(OV.rec_field(v, np.array([L]))[2])**2,      # reflection dip of the dark mode
                      bracket=(1.034*nu_c, 1.0356*nu_c, 1.037*nu_c)).x
lam1 = c0/nu1
s = np.linspace(0, L + 0.65*lam1, 30001); S = s/lam1; ins = s <= L

def referenced(nu):
    """E_rec times a constant phase: the overlap phase (m = 1, 3) or the phase at the maximum (dark mode),
    with the sign chosen so that the lobe next to the mirror is positive."""
    E, k, rr = OV.rec_field(nu, s)
    ov = np.trapezoid(E[ins], s[ins])
    big = abs(ov) > 0.05*np.trapezoid(abs(E[ins]), s[ins])
    Ep = E*np.exp(-1j*(np.angle(ov) if big else np.angle(E[np.argmax(abs(E)*ins)])))
    if Ep.real[np.argmax(abs(Ep[s < 0.3*L]))] < 0: Ep = -Ep
    return Ep, k, rr

if __name__ == '__main__':
    OV.fr.style(); C = OV.fr.SET1
    Em = np.exp(-1j*2*np.pi*nu1/c0*(s - L)) - np.exp(1j*2*np.pi*nu1/c0*(s + L))
    Em = Em*np.exp(-1j*np.angle(Em[np.argmax(abs(Em[s < 0.5*lam1]))]))
    rows = [(r'$m=1$, $\nu_1$', nu1, C[0], '-'), (r'$m=2$, $\nu_2=1.036\,\nu_1$', nu2, C[3], '-'),
            (r'$m=3$, $\nu_3=1.093\,\nu_1$', nu3, C[1], '-')]
    fig, axs = plt.subplots(2, 1, figsize=(7, 9.6), sharex=True)
    for ax in axs:
        ax.axvspan(0, L/lam1, color='0.93', zorder=0)
        for xg in OV.d*np.arange(1, OV.N+1)/lam1: ax.axvline(xg, color='0.6', lw=0.9, zorder=1)
        ax.axvline(0, color='k', lw=4, zorder=2)
    axs[0].axhline(0, color='k', lw=0.8)
    for lab, nu, col, ls in rows:
        Ep, k, rr = referenced(nu)
        axs[0].plot(S, Ep.real, color=col, ls=ls, label=lab)
        cum = 1j*k/2*np.concatenate([[0], np.cumsum(0.5*(np.where(ins, Ep, 0)[1:] + np.where(ins, Ep, 0)[:-1])*np.diff(s))])
        axs[1].plot(S, abs(cum)**2, color=col, ls=ls)
        print(lab, f'nu/nu1={nu/nu1:.4f} |beta|^2={abs(cum[-1])**2:.3g} 1-|r|^2={1-abs(rr)**2:.3f} max Re={Ep.real.max():.2f}')
    axs[0].plot(S, Em.real, color='0.35', ls='--', lw=1.8, label=r'mirror alone, $\nu_1$')
    k1 = 2*np.pi*nu1/c0
    cm = 1j*k1/2*np.concatenate([[0], np.cumsum(0.5*(np.where(ins, Em, 0)[1:] + np.where(ins, Em, 0)[:-1])*np.diff(s))])
    axs[1].plot(S, abs(cm)**2, color='0.35', ls='--', lw=1.8)
    axs[0].set_ylim(-7.5, 11.5); axs[0].set_ylabel(r'Re$\,[E_{\rm rec}(s)\,e^{-i\phi_0}]$')
    axs[0].legend(fontsize=14.5, loc='upper right', ncol=2, handlelength=1.4, columnspacing=1.0, borderaxespad=0.3)
    axs[0].text(0.02, 0.06, '(a)', transform=axs[0].transAxes, fontsize=19)
    axs[1].set_yscale('log'); axs[1].set_ylim(1e-2, 5e3); axs[1].set_xlim(-0.04, S[-1])
    axs[1].set_ylabel(r'$\left|\frac{ik}{2}\int_0^{s}E_{\rm rec}\,\mathrm{d}s^\prime\right|^2$'); axs[1].set_xlabel(r'$s/\lambda_1$')
    axs[1].text(0.02, 0.93, '(b)', transform=axs[1].transAxes, fontsize=19)
    fig.tight_layout(h_pad=0.6)
    out = os.path.join(HERE, '..', '..', 'results'); os.makedirs(out, exist_ok=True)
    for ext in ('pdf', 'png'): fig.savefig(os.path.join(out, 'fig_overlap_modes.' + ext), dpi=200)
