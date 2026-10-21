"""Fig. S-FDTD: finite wire length. Reads ../../../data/fdtd3d/*.npz, writes ../../../results/fig_fdtd3d_finite_wires.(pdf|png)."""
import os, sys, importlib.util, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..', '..', '..')
spec = importlib.util.spec_from_file_location('fr', os.path.join(ROOT, 'code', 'figures_reach.py'))
fr = importlib.util.module_from_spec(spec); spec.loader.exec_module(fr); fr.style(); SET1 = fr.SET1
D = os.path.join(ROOT, 'data', 'fdtd3d'); OUT = os.path.join(ROOT, 'results'); os.makedirs(OUT, exist_ok=True)
cols = {'bars': SET1[1], 'open': SET1[2], 'walls': SET1[0]}
labs = {'1-2': r'1–2 meV cartridge ($W=114\,\lambda$)', '0.25-0.5': r'0.25–0.5 meV cartridge ($W=28.5\,\lambda$)'}
fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))
for r, band in enumerate(['1-2', '0.25-0.5']):
    P = np.load(os.path.join(D, f'{band}_per.npz')); nuc = float(P['nu_c']); ref = P['b2_mid'].max()
    ax, ax2 = axs[r, 0], axs[r, 1]
    ax.plot(P['f']/nuc, P['b2_mid']/ref, 'k', lw=1.6, label='infinite wires')
    for fr_ in ['bars', 'open', 'walls']:
        M = np.load(os.path.join(D, f'{band}_{fr_}_W100mm.npz'))
        ax.plot(M['f']/nuc, M['b2_G']/ref, color=cols[fr_], label=f'{fr_}, Gaussian-collected')
        ax.plot(M['f']/nuc, M['b2_tot']/P['b2_tot'].max(), color=cols[fr_], ls='--', lw=1.6)
        ax2.plot(M['y']*1e3, M['prof']/P['prof'][0], color=cols[fr_], lw=1.8, label=fr_)
    ax.axhline(0.89, color='0.5', ls=':', lw=1.5)
    ax.set_xlim(0.995, 1.02); ax.set_ylim(0, 1.05)
    ax.set_xlabel(r'$\nu/\nu_c$'); ax.set_ylabel(r'$\beta^2/\beta^2_{\infty,\rm peak}$'); ax.set_title(labs[band], fontsize=16)
    ax2.axhline(1, color='k', lw=1); ax2.set_xlim(-55, 55)
    ax2.set_xlabel(r'$y$ along wires [mm]'); ax2.set_ylabel(r'$|E_{\rm out}(y)|^2/|E_{\infty}|^2$ at peak')
    ax2.set_title('field across the aperture', fontsize=16)
axs[0, 0].legend(fontsize=12, loc='upper left'); axs[0, 1].legend(fontsize=13, loc='lower center')
axs[0, 0].text(1.0125, 0.83, 'ideal uniform\naperture 0.89', fontsize=11, color='0.4')
fig.tight_layout()
for ext in ('pdf', 'png'): fig.savefig(os.path.join(OUT, 'fig_fdtd3d_finite_wires.' + ext), dpi=200)
