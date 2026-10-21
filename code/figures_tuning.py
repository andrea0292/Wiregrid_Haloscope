"""Fig. 4 of the Letter (v14): (a) optimal number of grids, (b) boost along the tuning paths.

Usage (from package/code):  python figures_tuning.py  -> ../paper/fig2.pdf (+ .png)
Inputs: ../data/nsweep.csv (N-sweep per cartridge, from nsweep_fig2_fig4.py)
        ../data/fig2b_tuning_paths.csv (sheet model and circular-wire tuning paths)
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figures_reach import style, SET1

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "..", "paper")
BANDS = ["0.25-0.5", "0.5-1", "1-2", "2-4", "4-8"]
COL = dict(zip(BANDS, [SET1[1], SET1[2], SET1[3], SET1[4], SET1[0]]))
NDES = dict(zip(BANDS, [15, 13, 11, 9, 8]))


def th_b(x): return x ** (2 / 3) / (1 + x) ** 2
def th_I(x): return x ** (4 / 3) / (1 + x) ** 3


def main():
    style()
    ns = pd.read_csv(os.path.join(D, "nsweep.csv"))
    tp = pd.read_csv(os.path.join(D, "fig2b_tuning_paths.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.9), gridspec_kw=dict(wspace=0.28))
    ax = axs[0]
    xx = np.geomspace(0.03, 60, 500)
    ax.axvline(0.5, ymin=0.5, color='0.65', lw=1.1, ls='--', zorder=1)
    ax.axvline(0.8, ymin=0.5, color='0.65', lw=1.1, zorder=1)
    ax.plot(xx, th_I(xx) / th_I(0.8), 'k-', lw=2.2, zorder=3, label=r'scan response $\mathcal{I}$')
    ax.plot(xx, th_b(xx) / th_b(0.5), 'k--', lw=2.2, zorder=3, label=r'peak boost $\beta^2$')
    for b in BANDS:
        g = ns[ns.band == b].sort_values("N"); c = COL[b]
        ax.plot(g.x, g.I / g.I.max(), 'o', color=c, ms=6.5, zorder=4)
        ax.plot(g.x, g.beta2 / g.beta2.max(), 's', mfc='w', mec=c, mew=1.5, ms=6, zorder=4)
        s = g[g.N == NDES[b]]
        ax.plot(s.x, s.I / g.I.max(), '*', color=c, ms=18, mec='k', mew=0.8, zorder=6)
    ax.set_xscale('log'); ax.set_xlim(0.05, 40); ax.set_ylim(0, 1.08)
    ax.set_xticks([0.1, 1, 10]); ax.set_xticklabels(['0.1', '1', '10'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r'$\gamma_{\rm rad}/\gamma_{\rm ohm}$')
    ax.set_ylabel('relative to maximum')
    # fractions above the frame, direction labels below the axis: nothing overlaps data or curves
    tr = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(0.5, 1.015, r'$\frac{1}{2}$', transform=tr, ha='center', va='bottom', fontsize=17, color='0.35')
    ax.text(0.8, 1.015, r'$\frac{4}{5}$', transform=tr, ha='center', va='bottom', fontsize=17, color='0.35')
    ax.text(0.0, -0.205, r'$\leftarrow$ thicker: ohmic', transform=ax.transAxes, ha='left', va='top', fontsize=13, color='0.35')
    ax.text(1.0, -0.205, r'thinner: radiative $\rightarrow$', transform=ax.transAxes, ha='right', va='top', fontsize=13, color='0.35')
    ax.legend(loc='lower center', fontsize=13, handlelength=1.8, bbox_to_anchor=(0.45, 0.0), borderaxespad=0.3)
    ax.text(0.04, 0.93, '(a)', transform=ax.transAxes, fontsize=17)

    ax = axs[1]
    for b in BANDS:
        c = COL[b]
        s = tp[(tp.band == b) & (tp.kind == 'sheet')].sort_values('m_meV')
        w = tp[(tp.band == b) & (tp.kind == 'circular_wire')].sort_values('m_meV')
        ax.plot(s.m_meV, s.beta2, color=c, lw=2.5)
        ax.plot(w.m_meV, w.beta2, 'o', mfc='w', mec=c, mew=2, ms=8)
        lo = float(b.split('-')[0])
        lab = rf'$N={NDES[b]}$' if b == BANDS[0] else rf'${NDES[b]}$'
        ax.text(lo * (1.25 if b != BANDS[0] else 1.32), s.beta2.max() * 1.18, lab, color=c, fontsize=16)
    ax.plot([], [], 'k-', lw=2.5, label='impedance sheets')
    ax.plot([], [], 'o', mfc='w', mec='k', mew=2, ms=8, label='circular wires')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.24, 8.3); ax.set_ylim(450, 1.1e4)
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(['0.25', '0.5', '1', '2', '4', '8'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r'$m_a$ [meV]'); ax.set_ylabel(r'$\beta^2$')
    ax.legend(loc='lower left', fontsize=14)
    ax.text(0.9, 0.9, '(b)', transform=ax.transAxes, fontsize=17)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, 'fig2.' + ext), dpi=150)
    print('written', OUT)


if __name__ == '__main__':
    main()
