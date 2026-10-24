"""Fig. S1 of the Letter (paper/fig_overlap_modes.png): receiving modes m=1,2,3 with real (solid) and imaginary (dotted) parts.
Run from code/validation/."""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import overlap_modes_figure as F
OV = F.OV; s = F.s; S = F.S; ins = F.ins; L = F.L; lam1 = F.lam1
OV.fr.style(); C = OV.fr.SET1
Em = np.exp(-1j*2*np.pi*F.nu1/F.c0*(s - L)) - np.exp(1j*2*np.pi*F.nu1/F.c0*(s + L))
Em = Em*np.exp(-1j*np.angle(Em[np.argmax(abs(Em[s < 0.5*lam1]))]))
rows = [(r'$m=1$, $\nu_1$', F.nu1, C[0]), (r'$m=2$, $\nu_2=1.036\,\nu_1$', F.nu2, C[3]),
        (r'$m=3$, $\nu_3=1.093\,\nu_1$', F.nu3, C[1])]
fig, axs = plt.subplots(2, 1, figsize=(7, 9.6), sharex=True)
for ax in axs:
    ax.axvspan(0, L/lam1, color='0.93', zorder=0)
    for xg in OV.d*np.arange(1, OV.N+1)/lam1: ax.axvline(xg, color='0.6', lw=0.9, zorder=1)
    ax.axvline(0, color='k', lw=4, zorder=2)
axs[0].axhline(0, color='k', lw=0.8)
for lab, nu, col in rows:
    Ep, k, rr = F.referenced(nu)
    axs[0].plot(S, Ep.real, color=col, ls='-', label=lab)
    axs[0].plot(S, Ep.imag, color=col, ls=':', lw=2.0)
    cum = 1j*k/2*np.concatenate([[0], np.cumsum(0.5*(np.where(ins, Ep, 0)[1:] + np.where(ins, Ep, 0)[:-1])*np.diff(s))])
    axs[1].plot(S, abs(cum)**2, color=col)
    E = Ep[ins]
    print(lab, f'max|Im|/max|Re|={abs(E.imag).max()/abs(E.real).max():.3f}  |intE|/int|E|={abs(np.trapezoid(E, s[ins]))/np.trapezoid(abs(E), s[ins]):.3f}')
axs[0].plot(S, Em.real, color='0.35', ls='--', lw=1.8, label=r'mirror alone, $\nu_1$')
axs[0].plot([], [], color='0.35', ls=':', lw=2.0, label='imaginary parts')
k1 = 2*np.pi*F.nu1/F.c0
cm = 1j*k1/2*np.concatenate([[0], np.cumsum(0.5*(np.where(ins, Em, 0)[1:] + np.where(ins, Em, 0)[:-1])*np.diff(s))])
axs[1].plot(S, abs(cm)**2, color='0.35', ls='--', lw=1.8)
axs[0].set_ylim(-7.5, 13.0); axs[0].set_ylabel(r'$E_{\rm rec}(s)\,e^{-i\phi_0}$')
axs[0].legend(fontsize=13.5, loc='upper right', ncol=2, handlelength=1.4, columnspacing=1.0, borderaxespad=0.3)
axs[0].text(0.02, 0.06, '(a)', transform=axs[0].transAxes, fontsize=19)
axs[1].set_yscale('log'); axs[1].set_ylim(1e-2, 5e3); axs[1].set_xlim(-0.04, S[-1])
axs[1].set_ylabel(r'$\left|\frac{ik}{2}\int_0^{s}E_{\rm rec}\,\mathrm{d}s^\prime\right|^2$'); axs[1].set_xlabel(r'$s/\lambda_1$')
axs[1].text(0.02, 0.93, '(b)', transform=axs[1].transAxes, fontsize=19)
fig.tight_layout(h_pad=0.6)
for out in (os.path.join(HERE, '..', '..', 'paper'), os.path.join(HERE, '..', '..', 'results')):
    fig.savefig(os.path.join(out, 'fig_overlap_modes.png'), dpi=200); fig.savefig(os.path.join(out, 'fig_overlap_modes.pdf'))
